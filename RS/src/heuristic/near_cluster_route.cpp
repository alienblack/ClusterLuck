#include "heuristic/near_cluster_route.h"

#include <algorithm>
#include <climits>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "forest/validator.h"
#include "heuristic/greedy_merge.h"
#include "io/newick.h"
#include "tree/induced.h"

namespace {

struct NearCandidate {
    std::vector<Leaf> leaves;
    int spill = 0;
    int trimmed = 0;
    int score = 0;
};

struct PseudoCandidate {
    std::vector<Leaf> leaves;
    std::vector<std::vector<Leaf>> components;
    int spill = 0;
    int score = 0;
};

constexpr CompId kMixedComponent = -2;

struct NearCladeRegion {
    std::vector<Leaf> seed_leaves;
    int common_size = 0;
    int spill = 0;
    int score = 0;
};

struct DepthWindowRegion {
    std::vector<Leaf> seed_leaves;
    int common_size = 0;
    int spill = 0;
    int min_depth = 0;
    int max_depth = 0;
    int score = 0;
};

struct RegionSolveStats {
    int bad_region = 0;
    int parse_fail = 0;
    int local_rebuild_fail = 0;
    int local_invalid = 0;
    int restart_gate_reject = 0;
    int slack_reject = 0;
    int duplicate_leaf = 0;
    int coverage_fail = 0;
    int global_rebuild_fail = 0;
    int global_invalid = 0;
    int accepted = 0;
};

struct RegionPatch {
    std::vector<CompId> region;
    std::vector<std::vector<Leaf>> replacement;
    int region_count = 0;
    int replacement_count = 0;
    int score = 0;
};

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

void maybe_dump_forest(const char* prefix_name, const Forest& forest,
                       const Tree& t1, int serial) {
    const char* prefix = std::getenv(prefix_name);
    if (prefix == nullptr || *prefix == '\0') return;
    std::string path = std::string(prefix) + "_" +
                       std::to_string(forest.active_count) + "_" +
                       std::to_string(serial) + ".out";
    std::ofstream out(path);
    if (!out) return;
    out << forest_to_output(forest, t1);
}

std::string leaf_key(const std::vector<Leaf>& leaves) {
    std::string key;
    for (Leaf leaf : leaves) {
        key += std::to_string(leaf);
        key.push_back(',');
    }
    return key;
}

bool same_induced(const Tree& t1, const Tree& t2,
                  const std::vector<Leaf>& leaves) {
    if (leaves.empty()) return false;
    try {
        InducedInfo i1 = compute_induced(t1, leaves);
        InducedInfo i2 = compute_induced(t2, leaves);
        return i1.signature == i2.signature;
    } catch (...) {
        return false;
    }
}

CompId subtree_component(const Tree& tree, const Forest& forest, NodeId node,
                         std::vector<std::pair<CompId, CompId>>& pairs) {
    if (tree.is_leaf(node)) {
        Leaf leaf = tree.label[node];
        if (leaf <= 0 || leaf >= static_cast<int>(forest.leaf_comp.size())) {
            return kMixedComponent;
        }
        return forest.leaf_comp[leaf];
    }
    CompId left = subtree_component(tree, forest, tree.left[node], pairs);
    CompId right = subtree_component(tree, forest, tree.right[node], pairs);
    if (left >= 0 && right >= 0 && left != right &&
        forest.comps[left].active && forest.comps[right].active) {
        pairs.push_back({left, right});
    }
    if (left == right) return left;
    return kMixedComponent;
}

std::vector<std::pair<CompId, CompId>> active_sibling_pairs(
    const Tree& tree, const Forest& forest) {
    std::vector<std::pair<CompId, CompId>> pairs;
    if (tree.root >= 0) subtree_component(tree, forest, tree.root, pairs);
    return pairs;
}

uint64_t pair_key(CompId a, CompId b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<unsigned>(a)) << 32) |
           static_cast<unsigned>(b);
}

void add_region_component(std::vector<CompId>& region,
                          std::vector<unsigned char>& seen,
                          const Forest& forest, CompId id, int limit) {
    if (id < 0 || id >= static_cast<CompId>(forest.comps.size()) ||
        !forest.comps[id].active || seen[id]) {
        return;
    }
    if (static_cast<int>(region.size()) >= limit) return;
    seen[id] = 1;
    region.push_back(id);
}

void add_subtree_components(std::vector<CompId>& region,
                            std::vector<unsigned char>& seen,
                            const Forest& forest, const Tree& tree,
                            NodeId root, int limit) {
    if (root < 0) return;
    std::vector<NodeId> stack{root};
    while (!stack.empty() && static_cast<int>(region.size()) < limit) {
        NodeId node = stack.back();
        stack.pop_back();
        if (tree.is_leaf(node)) {
            CompId id = forest.leaf_comp[tree.label[node]];
            add_region_component(region, seen, forest, id, limit);
            continue;
        }
        stack.push_back(tree.right[node]);
        stack.push_back(tree.left[node]);
    }
}

void add_local_components(std::vector<CompId>& region,
                          std::vector<unsigned char>& seen,
                          const Forest& forest, const Tree& tree,
                          NodeId root, int limit) {
    if (root < 0) return;
    NodeId node = root;
    int rounds = 0;
    while (node >= 0 && rounds < 3 && static_cast<int>(region.size()) < limit) {
        add_subtree_components(region, seen, forest, tree, node, limit);
        node = tree.parent[node];
        ++rounds;
    }
}

void add_path_side_components(std::vector<CompId>& region,
                              std::vector<unsigned char>& seen,
                              const Forest& forest, const Tree& tree,
                              NodeId from, NodeId stop, int limit) {
    NodeId node = from;
    while (node >= 0 && node != stop && static_cast<int>(region.size()) < limit) {
        NodeId parent = tree.parent[node];
        if (parent < 0) break;
        NodeId side = tree.left[parent] == node ? tree.right[parent] : tree.left[parent];
        add_subtree_components(region, seen, forest, tree, side, limit);
        node = parent;
    }
}

std::vector<CompId> edge_blockers(const std::vector<NodeId>& edges,
                                  const std::vector<CompId>& owner,
                                  CompId allowed_a, CompId allowed_b) {
    std::vector<CompId> blockers;
    for (NodeId edge : edges) {
        if (edge < 0 || edge >= static_cast<NodeId>(owner.size())) continue;
        CompId id = owner[edge];
        if (id >= 0 && id != allowed_a && id != allowed_b) {
            blockers.push_back(id);
        }
    }
    std::sort(blockers.begin(), blockers.end());
    blockers.erase(std::unique(blockers.begin(), blockers.end()), blockers.end());
    return blockers;
}

std::vector<CompId> build_conflict_region(const Forest& forest,
                                          const Tree& source,
                                          const Tree& other,
                                          bool source_is_t1,
                                          CompId a, CompId b,
                                          int max_components) {
    std::vector<CompId> region;
    std::vector<unsigned char> seen(forest.comps.size(), 0);
    add_region_component(region, seen, forest, a, max_components);
    add_region_component(region, seen, forest, b, max_components);

    std::vector<Leaf> merged = forest.comps[a].leaves;
    merged.insert(merged.end(), forest.comps[b].leaves.begin(),
                  forest.comps[b].leaves.end());
    std::sort(merged.begin(), merged.end());
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

    NodeId source_lca = source.lca_leaves(merged);
    NodeId other_lca = other.lca_leaves(merged);
    add_local_components(region, seen, forest, source, source_lca, max_components);
    add_local_components(region, seen, forest, other, other_lca, max_components);

    NodeId other_a = source_is_t1 ? forest.comps[a].root2 : forest.comps[a].root1;
    NodeId other_b = source_is_t1 ? forest.comps[b].root2 : forest.comps[b].root1;
    add_path_side_components(region, seen, forest, other, other_a, other_lca,
                             max_components);
    add_path_side_components(region, seen, forest, other, other_b, other_lca,
                             max_components);

    try {
        InducedInfo i1 = compute_induced(source_is_t1 ? source : other, merged);
        InducedInfo i2 = compute_induced(source_is_t1 ? other : source, merged);
        if (i1.signature == i2.signature) {
            const std::vector<CompId>& owner1 =
                source_is_t1 ? forest.owner1 : forest.owner2;
            const std::vector<CompId>& owner2 =
                source_is_t1 ? forest.owner2 : forest.owner1;
            std::vector<CompId> b1 = edge_blockers(i1.used_edges, owner1, a, b);
            std::vector<CompId> b2 = edge_blockers(i2.used_edges, owner2, a, b);
            for (CompId id : b1) {
                add_region_component(region, seen, forest, id, max_components);
            }
            for (CompId id : b2) {
                add_region_component(region, seen, forest, id, max_components);
            }
        }
    } catch (...) {
    }
    return region;
}

void collect_subtree_leaves(const Tree& tree, NodeId node,
                            std::vector<std::vector<Leaf>>& out) {
    if (tree.is_leaf(node)) {
        out[node].push_back(tree.label[node]);
        return;
    }
    collect_subtree_leaves(tree, tree.left[node], out);
    collect_subtree_leaves(tree, tree.right[node], out);
    out[node].reserve(out[tree.left[node]].size() + out[tree.right[node]].size());
    out[node].insert(out[node].end(), out[tree.left[node]].begin(),
                     out[tree.left[node]].end());
    out[node].insert(out[node].end(), out[tree.right[node]].begin(),
                     out[tree.right[node]].end());
    std::sort(out[node].begin(), out[node].end());
}

std::vector<std::vector<Leaf>> all_clade_leaves(const Tree& tree) {
    std::vector<std::vector<Leaf>> out(tree.node_count());
    if (tree.root >= 0) collect_subtree_leaves(tree, tree.root, out);
    return out;
}

std::vector<Leaf> intersection_sorted(const std::vector<Leaf>& a,
                                      const std::vector<Leaf>& b) {
    std::vector<Leaf> out;
    out.reserve(std::min(a.size(), b.size()));
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(out));
    return out;
}

std::vector<Leaf> union_sorted(const std::vector<Leaf>& a,
                               const std::vector<Leaf>& b) {
    std::vector<Leaf> out;
    out.reserve(a.size() + b.size());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                   std::back_inserter(out));
    return out;
}

std::string restricted_newick_rec(const Tree& tree, NodeId node,
                                  const std::vector<int>& original_to_local) {
    if (tree.is_leaf(node)) {
        Leaf leaf = tree.label[node];
        if (leaf > 0 && leaf < static_cast<int>(original_to_local.size()) &&
            original_to_local[leaf] > 0) {
            return std::to_string(original_to_local[leaf]);
        }
        return "";
    }
    std::string left =
        restricted_newick_rec(tree, tree.left[node], original_to_local);
    std::string right =
        restricted_newick_rec(tree, tree.right[node], original_to_local);
    if (left.empty()) return right;
    if (right.empty()) return left;
    return "(" + left + "," + right + ")";
}

std::string restricted_newick(const Tree& tree, const std::vector<Leaf>& leaves,
                              const std::vector<int>& original_to_local) {
    NodeId root = tree.lca_leaves(leaves);
    return restricted_newick_rec(tree, root, original_to_local) + ";";
}

bool try_trim_to_valid(const Tree& t1, const Tree& t2,
                       const std::vector<Leaf>& base, int min_size,
                       int max_trim, std::vector<Leaf>& out,
                       int& trimmed) {
    if (static_cast<int>(base.size()) < min_size) return false;
    if (same_induced(t1, t2, base)) {
        out = base;
        trimmed = 0;
        return true;
    }
    if (max_trim <= 0) return false;

    const int n = static_cast<int>(base.size());
    for (int i = 0; i < n; ++i) {
        if (n - 1 < min_size) break;
        std::vector<Leaf> trial;
        trial.reserve(n - 1);
        for (int a = 0; a < n; ++a) {
            if (a != i) trial.push_back(base[a]);
        }
        if (same_induced(t1, t2, trial)) {
            out = std::move(trial);
            trimmed = 1;
            return true;
        }
    }
    if (max_trim <= 1) return false;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (n - 2 < min_size) return false;
            std::vector<Leaf> trial;
            trial.reserve(n - 2);
            for (int a = 0; a < n; ++a) {
                if (a != i && a != j) trial.push_back(base[a]);
            }
            if (same_induced(t1, t2, trial)) {
                out = std::move(trial);
                trimmed = 2;
                return true;
            }
        }
    }
    if (max_trim <= 2 || n > env_int("RS_NEAR_CLUSTER_TRIM3_MAX_SIZE", 45)) {
        return false;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                if (n - 3 < min_size) return false;
                std::vector<Leaf> trial;
                trial.reserve(n - 3);
                for (int a = 0; a < n; ++a) {
                    if (a != i && a != j && a != k) trial.push_back(base[a]);
                }
                if (same_induced(t1, t2, trial)) {
                    out = std::move(trial);
                    trimmed = 3;
                    return true;
                }
            }
        }
    }
    return false;
}

std::vector<NearCandidate> generate_near_candidates(
    const Tree& t1, const Tree& t2, int min_island_size, int max_spill,
    int max_trim, int max_candidates) {
    std::vector<std::vector<Leaf>> c1 = all_clade_leaves(t1);
    std::vector<std::vector<Leaf>> c2 = all_clade_leaves(t2);
    std::vector<NearCandidate> candidates;
    std::unordered_set<std::string> seen;

    for (NodeId u = 0; u < t1.node_count(); ++u) {
        if (t1.is_leaf(u)) continue;
        const std::vector<Leaf>& a = c1[u];
        if (static_cast<int>(a.size()) < min_island_size) continue;
        for (NodeId v = 0; v < t2.node_count(); ++v) {
            if (t2.is_leaf(v)) continue;
            const std::vector<Leaf>& b = c2[v];
            if (static_cast<int>(b.size()) < min_island_size) continue;
            std::vector<Leaf> common = intersection_sorted(a, b);
            if (static_cast<int>(common.size()) < min_island_size) continue;
            int spill = static_cast<int>(a.size() + b.size() - 2 * common.size());
            if (spill > max_spill) continue;

            std::vector<Leaf> island;
            int trimmed = 0;
            if (!try_trim_to_valid(t1, t2, common, min_island_size, max_trim,
                                   island, trimmed)) {
                continue;
            }
            std::string key = leaf_key(island);
            if (!seen.insert(key).second) continue;
            NearCandidate candidate;
            candidate.leaves = std::move(island);
            candidate.spill = spill;
            candidate.trimmed = trimmed;
            candidate.score =
                static_cast<int>(candidate.leaves.size()) * 10 -
                18 * candidate.spill - 25 * candidate.trimmed;
            candidates.push_back(std::move(candidate));
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const NearCandidate& a, const NearCandidate& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.leaves.size() != b.leaves.size()) {
                      return a.leaves.size() > b.leaves.size();
                  }
                  if (a.spill != b.spill) return a.spill < b.spill;
                  if (a.trimmed != b.trimmed) return a.trimmed < b.trimmed;
                  return a.leaves < b.leaves;
              });
    if (static_cast<int>(candidates.size()) > max_candidates) {
        candidates.resize(max_candidates);
    }
    return candidates;
}

std::vector<NearCandidate> generate_raw_intersections(
    const Tree& t1, const Tree& t2, int min_island_size, int max_spill,
    int max_candidates) {
    std::vector<std::vector<Leaf>> c1 = all_clade_leaves(t1);
    std::vector<std::vector<Leaf>> c2 = all_clade_leaves(t2);
    std::vector<NearCandidate> candidates;
    std::unordered_set<std::string> seen;

    for (NodeId u = 0; u < t1.node_count(); ++u) {
        if (t1.is_leaf(u)) continue;
        const std::vector<Leaf>& a = c1[u];
        if (static_cast<int>(a.size()) < min_island_size) continue;
        for (NodeId v = 0; v < t2.node_count(); ++v) {
            if (t2.is_leaf(v)) continue;
            const std::vector<Leaf>& b = c2[v];
            if (static_cast<int>(b.size()) < min_island_size) continue;
            std::vector<Leaf> common = intersection_sorted(a, b);
            if (static_cast<int>(common.size()) < min_island_size) continue;
            int spill = static_cast<int>(a.size() + b.size() - 2 * common.size());
            if (spill > max_spill) continue;
            std::string key = leaf_key(common);
            if (!seen.insert(key).second) continue;
            NearCandidate candidate;
            candidate.leaves = std::move(common);
            candidate.spill = spill;
            candidate.score = static_cast<int>(candidate.leaves.size()) * 10 -
                              12 * candidate.spill;
            candidates.push_back(std::move(candidate));
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const NearCandidate& a, const NearCandidate& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.leaves.size() != b.leaves.size()) {
                      return a.leaves.size() > b.leaves.size();
                  }
                  if (a.spill != b.spill) return a.spill < b.spill;
                  return a.leaves < b.leaves;
              });
    if (static_cast<int>(candidates.size()) > max_candidates) {
        candidates.resize(max_candidates);
    }
    return candidates;
}

bool solve_pseudo_island(const Tree& t1, const Tree& t2,
                         const std::vector<Leaf>& leaves, int restart_budget,
                         std::vector<std::vector<Leaf>>& components) {
    if (leaves.size() <= 1) return false;
    std::vector<int> original_to_local(t1.n + 1, -1);
    std::vector<Leaf> local_to_original(leaves.size() + 1, -1);
    for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
        original_to_local[leaves[i]] = i + 1;
        local_to_original[i + 1] = leaves[i];
    }
    Tree local_t1;
    Tree local_t2;
    try {
        local_t1 = parse_newick_tree_h45_order(
            restricted_newick(t1, leaves, original_to_local),
            static_cast<int>(leaves.size()));
        local_t2 = parse_newick_tree_h45_order(
            restricted_newick(t2, leaves, original_to_local),
            static_cast<int>(leaves.size()));
    } catch (...) {
        return false;
    }

    Forest local = make_singleton_forest(static_cast<int>(leaves.size()));
    if (!rebuild_forest_data(local, local_t1, local_t2)) return false;
    local = run_h45_greedy_merge(std::move(local), local_t1, local_t2, 2, true, true);
    local = exact_output_repair(std::move(local), local_t1, local_t2);
    if (restart_budget > 0) {
        Forest restarted =
            run_h45_split_restarts(local, local_t1, local_t2, 2, restart_budget,
                                   static_cast<int>(leaves.size()) * 19, 0);
        restarted = exact_output_repair(std::move(restarted), local_t1, local_t2);
        if (validate_forest(restarted, local_t1, local_t2).ok &&
            restarted.active_count < local.active_count) {
            local = std::move(restarted);
        }
    }
    if (!validate_forest(local, local_t1, local_t2).ok) return false;

    components.clear();
    for (CompId id : active_components(local)) {
        std::vector<Leaf> component;
        for (Leaf local_leaf : local.comps[id].leaves) {
            component.push_back(local_to_original[local_leaf]);
        }
        std::sort(component.begin(), component.end());
        components.push_back(std::move(component));
    }
    std::sort(components.begin(), components.end());
    return !components.empty() &&
           static_cast<int>(components.size()) < static_cast<int>(leaves.size());
}

std::vector<PseudoCandidate> generate_pseudo_candidates(
    const Tree& t1, const Tree& t2, int min_island_size, int max_spill,
    int max_raw_candidates, int local_restart_budget) {
    std::vector<NearCandidate> raw =
        generate_raw_intersections(t1, t2, min_island_size, max_spill,
                                   max_raw_candidates);
    std::vector<PseudoCandidate> out;
    for (const NearCandidate& candidate : raw) {
        std::vector<std::vector<Leaf>> components;
        if (!solve_pseudo_island(t1, t2, candidate.leaves, local_restart_budget,
                                 components)) {
            continue;
        }
        int gain = static_cast<int>(candidate.leaves.size()) -
                   static_cast<int>(components.size());
        PseudoCandidate pseudo;
        pseudo.leaves = candidate.leaves;
        pseudo.components = std::move(components);
        pseudo.spill = candidate.spill;
        pseudo.score = 100 * gain + static_cast<int>(candidate.leaves.size()) -
                       10 * candidate.spill;
        out.push_back(std::move(pseudo));
    }
    std::sort(out.begin(), out.end(),
              [](const PseudoCandidate& a, const PseudoCandidate& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.leaves.size() != b.leaves.size()) {
                      return a.leaves.size() > b.leaves.size();
                  }
                  if (a.components.size() != b.components.size()) {
                      return a.components.size() < b.components.size();
                  }
                  return a.leaves < b.leaves;
              });
    return out;
}

bool solve_region_replacement(const Forest& base,
                              const std::vector<CompId>& region,
                              const Tree& t1, const Tree& t2,
                              int local_restart_budget,
                              int local_slack,
                              Forest& out,
                              RegionSolveStats* stats = nullptr) {
    if (region.empty()) {
        if (stats != nullptr) ++stats->bad_region;
        return false;
    }
    std::vector<unsigned char> in_region(base.comps.size(), 0);
    std::vector<Leaf> leaves;
    for (CompId id : region) {
        if (id < 0 || id >= static_cast<CompId>(base.comps.size()) ||
            !base.comps[id].active) {
            if (stats != nullptr) ++stats->bad_region;
            return false;
        }
        in_region[id] = 1;
        leaves.insert(leaves.end(), base.comps[id].leaves.begin(),
                      base.comps[id].leaves.end());
    }
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    if (leaves.size() <= 1) {
        if (stats != nullptr) ++stats->bad_region;
        return false;
    }

    std::vector<int> original_to_local(t1.n + 1, -1);
    std::vector<Leaf> local_to_original(leaves.size() + 1, -1);
    for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
        original_to_local[leaves[i]] = i + 1;
        local_to_original[i + 1] = leaves[i];
    }

    Tree local_t1;
    Tree local_t2;
    try {
        local_t1 = parse_newick_tree_h45_order(
            restricted_newick(t1, leaves, original_to_local),
            static_cast<int>(leaves.size()));
        local_t2 = parse_newick_tree_h45_order(
            restricted_newick(t2, leaves, original_to_local),
            static_cast<int>(leaves.size()));
    } catch (...) {
        if (stats != nullptr) ++stats->parse_fail;
        return false;
    }

    Forest local = make_singleton_forest(static_cast<int>(leaves.size()));
    if (!rebuild_forest_data(local, local_t1, local_t2)) {
        if (stats != nullptr) ++stats->local_rebuild_fail;
        return false;
    }
    local = run_h45_greedy_merge(std::move(local), local_t1, local_t2, 2, true, true);
    local = exact_output_repair(std::move(local), local_t1, local_t2);
    if (local_restart_budget > 0) {
        Forest restarted = run_h45_split_restarts(
            local, local_t1, local_t2, 2, local_restart_budget,
            static_cast<int>(leaves.size()) * 41, 0);
        restarted = exact_output_repair(std::move(restarted), local_t1, local_t2);
        if (validate_forest(restarted, local_t1, local_t2).ok &&
            restarted.active_count < local.active_count) {
            local = std::move(restarted);
        }
    }
    if (!validate_forest(local, local_t1, local_t2).ok) {
        if (stats != nullptr) ++stats->local_invalid;
        return false;
    }
    if (local.active_count >
        static_cast<int>(region.size()) + std::max(0, local_slack)) {
        if (stats != nullptr) ++stats->slack_reject;
        return false;
    }

    Forest trial;
    trial.leaf_comp.assign(t1.n + 1, -1);
    for (CompId id = 0; id < static_cast<CompId>(base.comps.size()); ++id) {
        if (!base.comps[id].active || in_region[id]) continue;
        Component component;
        component.leaves = base.comps[id].leaves;
        component.repr = base.comps[id].repr;
        for (Leaf leaf : component.leaves) {
            if (trial.leaf_comp[leaf] != -1) {
                if (stats != nullptr) ++stats->duplicate_leaf;
                return false;
            }
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (CompId id : active_components(local)) {
        Component component;
        for (Leaf local_leaf : local.comps[id].leaves) {
            component.leaves.push_back(local_to_original[local_leaf]);
        }
        std::sort(component.leaves.begin(), component.leaves.end());
        component.repr = component.leaves.front();
        for (Leaf leaf : component.leaves) {
            if (trial.leaf_comp[leaf] != -1) {
                if (stats != nullptr) ++stats->duplicate_leaf;
                return false;
            }
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (Leaf leaf = 1; leaf <= t1.n; ++leaf) {
        if (trial.leaf_comp[leaf] == -1) {
            if (stats != nullptr) ++stats->coverage_fail;
            return false;
        }
    }
    if (!rebuild_forest_data(trial, t1, t2)) {
        if (stats != nullptr) ++stats->global_rebuild_fail;
        return false;
    }
    if (!validate_forest(trial, t1, t2).ok) {
        if (stats != nullptr) ++stats->global_invalid;
        return false;
    }
    out = std::move(trial);
    if (stats != nullptr) ++stats->accepted;
    return true;
}

bool solve_region_replacement_portfolio(const Forest& base,
                                        const std::vector<CompId>& region,
                                        const Tree& t1, const Tree& t2,
                                        int local_restart_budget,
                                        int local_slack,
                                        Forest& out,
                                        RegionSolveStats* stats = nullptr) {
    if (region.empty()) {
        if (stats != nullptr) ++stats->bad_region;
        return false;
    }
    std::vector<unsigned char> in_region(base.comps.size(), 0);
    std::vector<Leaf> leaves;
    for (CompId id : region) {
        if (id < 0 || id >= static_cast<CompId>(base.comps.size()) ||
            !base.comps[id].active) {
            if (stats != nullptr) ++stats->bad_region;
            return false;
        }
        in_region[id] = 1;
        leaves.insert(leaves.end(), base.comps[id].leaves.begin(),
                      base.comps[id].leaves.end());
    }
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    if (leaves.size() <= 1) {
        if (stats != nullptr) ++stats->bad_region;
        return false;
    }

    std::vector<int> original_to_local(t1.n + 1, -1);
    std::vector<Leaf> local_to_original(leaves.size() + 1, -1);
    for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
        original_to_local[leaves[i]] = i + 1;
        local_to_original[i + 1] = leaves[i];
    }

    Tree local_t1;
    Tree local_t2;
    try {
        local_t1 = parse_newick_tree_h45_order(
            restricted_newick(t1, leaves, original_to_local),
            static_cast<int>(leaves.size()));
        local_t2 = parse_newick_tree_h45_order(
            restricted_newick(t2, leaves, original_to_local),
            static_cast<int>(leaves.size()));
    } catch (...) {
        if (stats != nullptr) ++stats->parse_fail;
        return false;
    }

    Forest local_best;
    int best_count = INT_MAX;
    std::vector<Forest> local_by_mode;
    std::vector<int> modes;
    if (env_int("RS_SHALLOW_REPLACE_CORE_MODES", 0) != 0) {
        modes = {14, 2, 20, 1};
    } else {
        modes = {2, 14, 20, 1, 5, 6, 0};
    }
    int variant = 0;
    for (int mode : modes) {
        Forest local = make_singleton_forest(static_cast<int>(leaves.size()));
        if (!rebuild_forest_data(local, local_t1, local_t2)) {
            if (stats != nullptr) ++stats->local_rebuild_fail;
            continue;
        }
        local = run_h45_greedy_merge(std::move(local), local_t1, local_t2,
                                     mode, true, true);
        local = exact_output_repair(std::move(local), local_t1, local_t2);
        if (validate_forest(local, local_t1, local_t2).ok &&
            local.active_count < best_count) {
            local_best = local;
            best_count = local.active_count;
        }
        local_by_mode.push_back(std::move(local));
        ++variant;
    }

    const int restart_gate_margin =
        env_int("RS_SHALLOW_REPLACE_RESTART_GATE_MARGIN", -1);
    if (local_restart_budget > 0 && restart_gate_margin >= 0 &&
        best_count > static_cast<int>(region.size()) +
                         std::max(0, local_slack) + restart_gate_margin) {
        if (stats != nullptr) ++stats->restart_gate_reject;
        return false;
    }

    if (local_restart_budget > 0) {
        variant = 0;
        const bool restart_from_best =
            env_int("RS_SHALLOW_REPLACE_RESTART_FROM_BEST", 0) != 0;
        for (int mode_index = 0; mode_index < static_cast<int>(modes.size());
             ++mode_index) {
            int mode = modes[mode_index];
            const int per_mode_budget =
                std::max(1, local_restart_budget /
                                static_cast<int>(modes.size()));
            Forest restart_seed =
                restart_from_best ? local_best : local_by_mode[mode_index];
            if (restart_seed.active_count == 0) {
                restart_seed = make_singleton_forest(static_cast<int>(leaves.size()));
                if (!rebuild_forest_data(restart_seed, local_t1, local_t2)) {
                    if (stats != nullptr) ++stats->local_rebuild_fail;
                    ++variant;
                    continue;
                }
            }
            Forest restarted = run_h45_split_restarts(
                restart_seed, local_t1, local_t2, mode, per_mode_budget,
                static_cast<int>(leaves.size()) * 53 + 997 * variant,
                env_int("RS_SHALLOW_REPLACE_LOCAL_EXTRA_SPLITS", 1),
                true);
            restarted = exact_output_repair(std::move(restarted), local_t1, local_t2);
            if (validate_forest(restarted, local_t1, local_t2).ok &&
                restarted.active_count < best_count) {
                local_best = std::move(restarted);
                best_count = local_best.active_count;
            }
            ++variant;
        }
    }
    if (best_count == INT_MAX ||
        !validate_forest(local_best, local_t1, local_t2).ok) {
        if (stats != nullptr) ++stats->local_invalid;
        return false;
    }
    if (local_best.active_count >
        static_cast<int>(region.size()) + std::max(0, local_slack)) {
        if (stats != nullptr) ++stats->slack_reject;
        return false;
    }

    Forest trial;
    trial.leaf_comp.assign(t1.n + 1, -1);
    for (CompId id = 0; id < static_cast<CompId>(base.comps.size()); ++id) {
        if (!base.comps[id].active || in_region[id]) continue;
        Component component;
        component.leaves = base.comps[id].leaves;
        component.repr = base.comps[id].repr;
        for (Leaf leaf : component.leaves) {
            if (trial.leaf_comp[leaf] != -1) {
                if (stats != nullptr) ++stats->duplicate_leaf;
                return false;
            }
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (CompId id : active_components(local_best)) {
        Component component;
        for (Leaf local_leaf : local_best.comps[id].leaves) {
            component.leaves.push_back(local_to_original[local_leaf]);
        }
        std::sort(component.leaves.begin(), component.leaves.end());
        component.repr = component.leaves.front();
        for (Leaf leaf : component.leaves) {
            if (trial.leaf_comp[leaf] != -1) {
                if (stats != nullptr) ++stats->duplicate_leaf;
                return false;
            }
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (Leaf leaf = 1; leaf <= t1.n; ++leaf) {
        if (trial.leaf_comp[leaf] == -1) {
            if (stats != nullptr) ++stats->coverage_fail;
            return false;
        }
    }
    if (!rebuild_forest_data(trial, t1, t2)) {
        if (stats != nullptr) ++stats->global_rebuild_fail;
        return false;
    }
    if (!validate_forest(trial, t1, t2).ok) {
        if (stats != nullptr) ++stats->global_invalid;
        return false;
    }
    out = std::move(trial);
    if (stats != nullptr) ++stats->accepted;
    return true;
}

std::vector<CompId> components_intersecting_leaves(const Forest& forest,
                                                   const std::vector<Leaf>& leaves,
                                                   int max_components) {
    std::vector<CompId> region;
    std::vector<unsigned char> seen(forest.comps.size(), 0);
    for (Leaf leaf : leaves) {
        if (leaf <= 0 || leaf >= static_cast<int>(forest.leaf_comp.size())) {
            continue;
        }
        add_region_component(region, seen, forest, forest.leaf_comp[leaf],
                             max_components + 1);
        if (static_cast<int>(region.size()) > max_components) {
            return {};
        }
    }
    return region;
}

bool edge_close_region(std::vector<CompId>& region, const Forest& forest,
                       const Tree& t1, const Tree& t2, int max_components) {
    std::vector<unsigned char> seen(forest.comps.size(), 0);
    for (CompId id : region) {
        if (id < 0 || id >= static_cast<CompId>(forest.comps.size()) ||
            !forest.comps[id].active) {
            return false;
        }
        seen[id] = 1;
    }
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<Leaf> leaves;
        for (CompId id : region) {
            leaves.insert(leaves.end(), forest.comps[id].leaves.begin(),
                          forest.comps[id].leaves.end());
        }
        std::sort(leaves.begin(), leaves.end());
        leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
        if (leaves.size() <= 1) return false;

        InducedInfo i1;
        InducedInfo i2;
        try {
            i1 = compute_induced(t1, leaves);
            i2 = compute_induced(t2, leaves);
        } catch (...) {
            return false;
        }
        for (NodeId edge : i1.used_edges) {
            if (edge < 0 || edge >= static_cast<NodeId>(forest.owner1.size())) {
                continue;
            }
            CompId id = forest.owner1[edge];
            if (id >= 0 && !seen[id]) {
                if (static_cast<int>(region.size()) >= max_components) return false;
                seen[id] = 1;
                region.push_back(id);
                changed = true;
            }
        }
        for (NodeId edge : i2.used_edges) {
            if (edge < 0 || edge >= static_cast<NodeId>(forest.owner2.size())) {
                continue;
            }
            CompId id = forest.owner2[edge];
            if (id >= 0 && !seen[id]) {
                if (static_cast<int>(region.size()) >= max_components) return false;
                seen[id] = 1;
                region.push_back(id);
                changed = true;
            }
        }
    }
    return true;
}

bool extract_region_patch(const Forest& base, const Forest& trial,
                          const std::vector<CompId>& region, int n,
                          RegionPatch& patch) {
    std::vector<unsigned char> region_leaf(n + 1, 0);
    std::vector<CompId> sorted_region = region;
    std::sort(sorted_region.begin(), sorted_region.end());
    sorted_region.erase(std::unique(sorted_region.begin(), sorted_region.end()),
                        sorted_region.end());
    for (CompId id : sorted_region) {
        if (id < 0 || id >= static_cast<CompId>(base.comps.size()) ||
            !base.comps[id].active) {
            return false;
        }
        for (Leaf leaf : base.comps[id].leaves) {
            if (leaf < 1 || leaf > n) return false;
            region_leaf[leaf] = 1;
        }
    }

    std::vector<std::vector<Leaf>> replacement;
    for (CompId id : active_components(trial)) {
        bool touches = false;
        bool outside = false;
        for (Leaf leaf : trial.comps[id].leaves) {
            if (leaf < 1 || leaf > n) return false;
            if (region_leaf[leaf]) {
                touches = true;
            } else {
                outside = true;
            }
        }
        if (!touches) continue;
        if (outside) return false;
        std::vector<Leaf> component = trial.comps[id].leaves;
        std::sort(component.begin(), component.end());
        replacement.push_back(std::move(component));
    }
    if (replacement.empty()) return false;
    std::sort(replacement.begin(), replacement.end());

    patch.region = std::move(sorted_region);
    patch.replacement = std::move(replacement);
    patch.region_count = static_cast<int>(patch.region.size());
    patch.replacement_count = static_cast<int>(patch.replacement.size());
    patch.score = 20 * (patch.region_count - patch.replacement_count) +
                  patch.region_count;
    return true;
}

bool patches_disjoint(const RegionPatch& a, const RegionPatch& b) {
    return std::find_first_of(a.region.begin(), a.region.end(),
                              b.region.begin(), b.region.end()) ==
           a.region.end();
}

bool apply_region_patches(const Forest& base,
                          const std::vector<const RegionPatch*>& patches,
                          const Tree& t1, const Tree& t2, Forest& out) {
    std::vector<unsigned char> removed(base.comps.size(), 0);
    for (const RegionPatch* patch : patches) {
        for (CompId id : patch->region) {
            if (id < 0 || id >= static_cast<CompId>(base.comps.size()) ||
                !base.comps[id].active || removed[id]) {
                return false;
            }
            removed[id] = 1;
        }
    }

    Forest trial;
    trial.leaf_comp.assign(t1.n + 1, -1);
    for (CompId id = 0; id < static_cast<CompId>(base.comps.size()); ++id) {
        if (!base.comps[id].active || removed[id]) continue;
        Component component;
        component.leaves = base.comps[id].leaves;
        component.repr = base.comps[id].repr;
        for (Leaf leaf : component.leaves) {
            if (leaf < 1 || leaf > t1.n || trial.leaf_comp[leaf] != -1) {
                return false;
            }
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (const RegionPatch* patch : patches) {
        for (const std::vector<Leaf>& leaves : patch->replacement) {
            Component component;
            component.leaves = leaves;
            component.repr = component.leaves.front();
            for (Leaf leaf : component.leaves) {
                if (leaf < 1 || leaf > t1.n || trial.leaf_comp[leaf] != -1) {
                    return false;
                }
                trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
            }
            trial.comps.push_back(std::move(component));
            ++trial.active_count;
        }
    }
    for (Leaf leaf = 1; leaf <= t1.n; ++leaf) {
        if (trial.leaf_comp[leaf] == -1) return false;
    }
    if (!rebuild_forest_data(trial, t1, t2)) return false;
    if (!validate_forest(trial, t1, t2).ok) return false;
    out = std::move(trial);
    return true;
}

std::vector<NearCladeRegion> generate_near_clade_regions(
    const Tree& t1, const Tree& t2, int min_common, int max_spill,
    int candidate_limit) {
    std::vector<std::vector<Leaf>> c1 = all_clade_leaves(t1);
    std::vector<std::vector<Leaf>> c2 = all_clade_leaves(t2);
    std::vector<NearCladeRegion> candidates;
    std::unordered_set<std::string> seen;

    for (NodeId u = 0; u < t1.node_count(); ++u) {
        if (t1.is_leaf(u)) continue;
        const std::vector<Leaf>& a = c1[u];
        if (static_cast<int>(a.size()) < min_common) continue;
        for (NodeId v = 0; v < t2.node_count(); ++v) {
            if (t2.is_leaf(v)) continue;
            const std::vector<Leaf>& b = c2[v];
            if (static_cast<int>(b.size()) < min_common) continue;
            std::vector<Leaf> common = intersection_sorted(a, b);
            if (static_cast<int>(common.size()) < min_common) continue;
            int spill = static_cast<int>(a.size() + b.size() -
                                         2 * common.size());
            if (spill > max_spill) continue;
            std::vector<Leaf> seed = union_sorted(a, b);
            std::string key = leaf_key(seed);
            if (!seen.insert(key).second) continue;

            NearCladeRegion region;
            region.seed_leaves = std::move(seed);
            region.common_size = static_cast<int>(common.size());
            region.spill = spill;
            region.score = 25 * region.common_size - 9 * spill -
                           static_cast<int>(region.seed_leaves.size());
            candidates.push_back(std::move(region));
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const NearCladeRegion& a, const NearCladeRegion& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.common_size != b.common_size) {
                      return a.common_size > b.common_size;
                  }
                  if (a.spill != b.spill) return a.spill < b.spill;
                  return a.seed_leaves < b.seed_leaves;
              });
    if (candidate_limit >= 0 &&
        static_cast<int>(candidates.size()) > candidate_limit) {
        candidates.resize(candidate_limit);
    }
    return candidates;
}

Forest run_near_clade_polluter_route(Forest initial, const Tree& t1,
                                     const Tree& t2, int max_leaves,
                                     int min_common, int max_spill,
                                     int candidate_limit,
                                     int max_region_components,
                                     int max_region_leaves,
                                     int local_restart_budget,
                                     int local_slack,
                                     int final_restart_budget,
                                     int target_count) {
    if (t1.n > max_leaves || min_common < 2 || max_spill < 0 ||
        candidate_limit == 0) {
        return initial;
    }
    const bool debug = std::getenv("RS_NEAR_CLUSTER_DEBUG") != nullptr;
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = std::move(initial);

    std::vector<NearCladeRegion> candidates = generate_near_clade_regions(
        t1, t2, min_common, max_spill, candidate_limit);
    int considered = 0;
    int solved = 0;
    int improved = 0;
    int dump_serial = 0;
    int skipped_region_size = 0;
    int skipped_region_leaves = 0;
    int skipped_edge_closure = 0;
    int combo_considered = 0;
    int combo_valid = 0;
    int combo_improved = 0;
    RegionSolveStats stats;
    std::vector<RegionPatch> patches;
    const int patch_pool_limit =
        env_int("RS_NEAR_CLADE_PATCH_POOL_LIMIT", 48);
    for (const NearCladeRegion& candidate : candidates) {
        ++considered;
        std::vector<CompId> region = components_intersecting_leaves(
            best, candidate.seed_leaves, max_region_components);
        if (region.size() < 2) {
            ++skipped_region_size;
            continue;
        }
        if (!edge_close_region(region, best, t1, t2, max_region_components)) {
            ++skipped_edge_closure;
            continue;
        }
        std::vector<Leaf> region_leaves;
        for (CompId id : region) {
            region_leaves.insert(region_leaves.end(), best.comps[id].leaves.begin(),
                                 best.comps[id].leaves.end());
        }
        std::sort(region_leaves.begin(), region_leaves.end());
        region_leaves.erase(std::unique(region_leaves.begin(), region_leaves.end()),
                            region_leaves.end());
        if (static_cast<int>(region_leaves.size()) > max_region_leaves) {
            ++skipped_region_leaves;
            continue;
        }

        Forest trial;
        if (!solve_region_replacement(best, region, t1, t2, local_restart_budget,
                                      local_slack, trial, &stats)) {
            continue;
        }
        ++solved;
        if (static_cast<int>(patches.size()) < patch_pool_limit) {
            RegionPatch patch;
            if (extract_region_patch(best, trial, region, t1.n, patch)) {
                patches.push_back(std::move(patch));
            }
        }
        for (int mode : {14, 2, 1}) {
            Forest merged = run_h45_greedy_merge(trial, t1, t2, mode, true, true);
            merged = exact_output_repair(std::move(merged), t1, t2);
            if (!validate_forest(merged, t1, t2).ok) continue;
            if (final_restart_budget > 0) {
                Forest restarted = run_h45_split_restarts(
                    merged, t1, t2, mode, final_restart_budget,
                    31000 + mode * 211 + considered, 1, false, target_count, 0.0);
                restarted = exact_output_repair(std::move(restarted), t1, t2);
                if (validate_forest(restarted, t1, t2).ok &&
                    restarted.active_count < merged.active_count) {
                    merged = std::move(restarted);
                }
            }
            if (merged.active_count < best.active_count) {
                best = std::move(merged);
                ++improved;
                maybe_dump_forest("RS_SHALLOW_REPLACE_DUMP_PREFIX", best, t1,
                                  ++dump_serial);
                if (debug) {
                    std::cerr << "rs_near_clade_polluter_improve components="
                              << best.active_count
                              << " common=" << candidate.common_size
                              << " spill=" << candidate.spill
                              << " region_components=" << region.size()
                              << " region_leaves=" << region_leaves.size()
                              << " mode=" << mode << '\n';
                }
                if (target_count > 0 && best.active_count <= target_count) {
                    return best;
                }
                break;
            }
        }
    }
    if (improved == 0 &&
        env_int("RS_NEAR_CLADE_PATCH_COMBO_ROUTE", 1) != 0 &&
        patches.size() >= 2) {
        std::sort(patches.begin(), patches.end(),
                  [](const RegionPatch& a, const RegionPatch& b) {
                      if (a.score != b.score) return a.score > b.score;
                      if (a.replacement_count != b.replacement_count) {
                          return a.replacement_count < b.replacement_count;
                      }
                      return a.region < b.region;
                  });
        const int combo_limit =
            env_int("RS_NEAR_CLADE_PATCH_COMBO_LIMIT", 600);
        const int combo_restart_budget =
            env_int("RS_NEAR_CLADE_PATCH_COMBO_RESTART_BUDGET",
                    final_restart_budget);
        for (int i = 0; i < static_cast<int>(patches.size()) &&
                        combo_considered < combo_limit; ++i) {
            for (int j = i + 1; j < static_cast<int>(patches.size()) &&
                            combo_considered < combo_limit; ++j) {
                if (!patches_disjoint(patches[i], patches[j])) continue;
                ++combo_considered;
                Forest trial;
                std::vector<const RegionPatch*> combo = {&patches[i], &patches[j]};
                if (!apply_region_patches(best, combo, t1, t2, trial)) continue;
                ++combo_valid;
                for (int mode : {14, 2}) {
                    Forest merged =
                        run_h45_greedy_merge(trial, t1, t2, mode, true, true);
                    merged = exact_output_repair(std::move(merged), t1, t2);
                    if (!validate_forest(merged, t1, t2).ok) continue;
                    if (combo_restart_budget > 0) {
                        Forest restarted = run_h45_split_restarts(
                            merged, t1, t2, mode, combo_restart_budget,
                            42000 + combo_considered * 17 + mode, 1, false,
                            target_count, 0.0);
                        restarted =
                            exact_output_repair(std::move(restarted), t1, t2);
                        if (validate_forest(restarted, t1, t2).ok &&
                            restarted.active_count < merged.active_count) {
                            merged = std::move(restarted);
                        }
                    }
                    if (merged.active_count < best.active_count) {
                        best = std::move(merged);
                        ++combo_improved;
                        if (debug) {
                            std::cerr
                                << "rs_near_clade_patch_combo_improve components="
                                << best.active_count
                                << " combo_considered=" << combo_considered
                                << " mode=" << mode << '\n';
                        }
                        return best;
                    }
                }
            }
        }
    }
    if (debug) {
        std::cerr << "rs_near_clade_polluter_done components="
                  << best.active_count
                  << " candidates=" << candidates.size()
                  << " considered=" << considered
                  << " solved=" << solved
                  << " improved=" << improved
                  << " patches=" << patches.size()
                  << " combo_considered=" << combo_considered
                  << " combo_valid=" << combo_valid
                  << " combo_improved=" << combo_improved
                  << " skip_region_size=" << skipped_region_size
                  << " skip_region_leaves=" << skipped_region_leaves
                  << " skip_edge_closure=" << skipped_edge_closure
                  << " bad_region=" << stats.bad_region
                  << " parse_fail=" << stats.parse_fail
                  << " local_rebuild_fail=" << stats.local_rebuild_fail
                  << " local_invalid=" << stats.local_invalid
                  << " restart_gate_reject=" << stats.restart_gate_reject
                  << " slack_reject=" << stats.slack_reject
                  << " duplicate_leaf=" << stats.duplicate_leaf
                  << " coverage_fail=" << stats.coverage_fail
                  << " global_rebuild_fail=" << stats.global_rebuild_fail
                  << " global_invalid=" << stats.global_invalid
                  << " accepted=" << stats.accepted << '\n';
    }
    return best;
}

Forest run_conflict_polluter_route(Forest initial, const Tree& t1,
                                   const Tree& t2, int max_leaves,
                                   int pair_limit, int max_region_components,
                                   int max_region_leaves,
                                   int local_restart_budget,
                                   int local_slack, int final_restart_budget,
                                   int target_count) {
    if (t1.n > max_leaves || pair_limit <= 0 || max_region_components <= 1 ||
        max_region_leaves <= 1) {
        return initial;
    }
    const bool debug = std::getenv("RS_NEAR_CLUSTER_DEBUG") != nullptr;
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = std::move(initial);

    int considered = 0;
    int solved = 0;
    int improved = 0;
    int skipped_region_size = 0;
    int skipped_region_leaves = 0;
    int skipped_edge_closure = 0;
    RegionSolveStats stats;
    for (int pass = 0; pass < env_int("RS_CONFLICT_POLLUTER_PASSES", 2); ++pass) {
        std::vector<std::pair<CompId, CompId>> pairs =
            active_sibling_pairs(t1, best);
        std::vector<std::pair<CompId, CompId>> p2 =
            active_sibling_pairs(t2, best);
        pairs.insert(pairs.end(), p2.begin(), p2.end());
        std::vector<std::pair<CompId, CompId>> unique_pairs;
        std::unordered_set<uint64_t> seen;
        for (auto pair : pairs) {
            if (pair.first == pair.second) continue;
            uint64_t key = pair_key(pair.first, pair.second);
            if (seen.insert(key).second) unique_pairs.push_back(pair);
        }
        std::sort(unique_pairs.begin(), unique_pairs.end(),
                  [&](const auto& a, const auto& b) {
                      int sa = static_cast<int>(best.comps[a.first].leaves.size() +
                                                best.comps[a.second].leaves.size());
                      int sb = static_cast<int>(best.comps[b.first].leaves.size() +
                                                best.comps[b.second].leaves.size());
                      if (sa != sb) return sa > sb;
                      return pair_key(a.first, a.second) < pair_key(b.first, b.second);
                  });
        if (static_cast<int>(unique_pairs.size()) > pair_limit) {
            unique_pairs.resize(pair_limit);
        }

        bool pass_improved = false;
        for (auto [a, b] : unique_pairs) {
            if (!best.comps[a].active || !best.comps[b].active) continue;
            for (int side = 0; side < 2; ++side) {
                std::vector<CompId> region = build_conflict_region(
                    best, side == 0 ? t1 : t2, side == 0 ? t2 : t1,
                    side == 0, a, b, max_region_components);
                ++considered;
                if (static_cast<int>(region.size()) < 2) {
                    ++skipped_region_size;
                    continue;
                }
                if (!edge_close_region(region, best, t1, t2,
                                       max_region_components)) {
                    ++skipped_edge_closure;
                    continue;
                }
                std::vector<Leaf> leaves;
                for (CompId id : region) {
                    leaves.insert(leaves.end(), best.comps[id].leaves.begin(),
                                  best.comps[id].leaves.end());
                }
                std::sort(leaves.begin(), leaves.end());
                leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
                if (static_cast<int>(leaves.size()) > max_region_leaves) {
                    ++skipped_region_leaves;
                    continue;
                }

                Forest trial;
                if (!solve_region_replacement(best, region, t1, t2,
                                              local_restart_budget, local_slack,
                                              trial, &stats)) {
                    continue;
                }
                ++solved;
                for (int mode : {14, 2, 1}) {
                    Forest merged =
                        run_h45_greedy_merge(trial, t1, t2, mode, true, true);
                    merged = exact_output_repair(std::move(merged), t1, t2);
                    if (!validate_forest(merged, t1, t2).ok) continue;
                    if (final_restart_budget > 0) {
                        Forest restarted = run_h45_split_restarts(
                            merged, t1, t2, mode, final_restart_budget,
                            23000 + pass * 997 + mode * 37, 1, false,
                            target_count, 0.0);
                        restarted = exact_output_repair(std::move(restarted), t1, t2);
                        if (validate_forest(restarted, t1, t2).ok &&
                            restarted.active_count < merged.active_count) {
                            merged = std::move(restarted);
                        }
                    }
                    if (merged.active_count < best.active_count) {
                        best = std::move(merged);
                        ++improved;
                        pass_improved = true;
                        if (debug) {
                            std::cerr
                                << "rs_conflict_polluter_improve components="
                                << best.active_count
                                << " region_components=" << region.size()
                                << " region_leaves=" << leaves.size()
                                << " mode=" << mode << '\n';
                        }
                        if (target_count > 0 &&
                            best.active_count <= target_count) {
                            return best;
                        }
                        break;
                    }
                }
                if (pass_improved) break;
            }
            if (pass_improved) break;
        }
        if (!pass_improved) break;
    }
    if (debug) {
        std::cerr << "rs_conflict_polluter_done components=" << best.active_count
                  << " considered=" << considered
                  << " solved=" << solved
                  << " improved=" << improved
                  << " skip_region_size=" << skipped_region_size
                  << " skip_region_leaves=" << skipped_region_leaves
                  << " skip_edge_closure=" << skipped_edge_closure
                  << " bad_region=" << stats.bad_region
                  << " parse_fail=" << stats.parse_fail
                  << " local_rebuild_fail=" << stats.local_rebuild_fail
                  << " local_invalid=" << stats.local_invalid
                  << " restart_gate_reject=" << stats.restart_gate_reject
                  << " slack_reject=" << stats.slack_reject
                  << " duplicate_leaf=" << stats.duplicate_leaf
                  << " coverage_fail=" << stats.coverage_fail
                  << " global_rebuild_fail=" << stats.global_rebuild_fail
                  << " global_invalid=" << stats.global_invalid
                  << " accepted=" << stats.accepted << '\n';
    }
    return best;
}

std::vector<DepthWindowRegion> generate_depth_window_regions(
    const Tree& t1, const Tree& t2, int min_depth, int max_depth,
    int min_seed_leaves, int max_seed_leaves, int candidate_limit) {
    std::vector<std::vector<Leaf>> c1 = all_clade_leaves(t1);
    std::vector<std::vector<Leaf>> c2 = all_clade_leaves(t2);
    std::vector<NodeId> nodes1;
    std::vector<NodeId> nodes2;
    for (NodeId u = 0; u < t1.node_count(); ++u) {
        if (t1.is_leaf(u)) continue;
        int size = static_cast<int>(c1[u].size());
        if (t1.depth[u] >= min_depth && t1.depth[u] <= max_depth &&
            size >= min_seed_leaves && size <= max_seed_leaves) {
            nodes1.push_back(u);
        }
    }
    for (NodeId v = 0; v < t2.node_count(); ++v) {
        if (t2.is_leaf(v)) continue;
        int size = static_cast<int>(c2[v].size());
        if (t2.depth[v] >= min_depth && t2.depth[v] <= max_depth &&
            size >= min_seed_leaves && size <= max_seed_leaves) {
            nodes2.push_back(v);
        }
    }

    std::vector<DepthWindowRegion> out;
    std::unordered_set<std::string> seen;
    auto add_candidate = [&](std::vector<Leaf> seed, int common, int spill,
                             int dmin, int dmax) {
        if (static_cast<int>(seed.size()) < min_seed_leaves ||
            static_cast<int>(seed.size()) > max_seed_leaves) {
            return;
        }
        std::sort(seed.begin(), seed.end());
        seed.erase(std::unique(seed.begin(), seed.end()), seed.end());
        std::string key = leaf_key(seed);
        if (!seen.insert(key).second) return;
        DepthWindowRegion region;
        region.seed_leaves = std::move(seed);
        region.common_size = common;
        region.spill = spill;
        region.min_depth = dmin;
        region.max_depth = dmax;
        region.score = 40 * common - 5 * spill +
                       static_cast<int>(region.seed_leaves.size());
        out.push_back(std::move(region));
    };

    if (env_int("RS_SHALLOW_REPLACE_SINGLE_CLADE", 1) != 0) {
        for (NodeId u : nodes1) {
            add_candidate(c1[u], static_cast<int>(c1[u].size()), 0,
                          t1.depth[u], t1.depth[u]);
        }
        for (NodeId v : nodes2) {
            add_candidate(c2[v], static_cast<int>(c2[v].size()), 0,
                          t2.depth[v], t2.depth[v]);
        }
    }
    for (NodeId u : nodes1) {
        std::vector<DepthWindowRegion> best_for_u;
        for (NodeId v : nodes2) {
            std::vector<Leaf> common = intersection_sorted(c1[u], c2[v]);
            int common_size = static_cast<int>(common.size());
            if (common_size < min_seed_leaves / 3) continue;
            int spill = static_cast<int>(c1[u].size() + c2[v].size() -
                                         2 * common.size());
            if (spill > max_seed_leaves) continue;
            std::vector<Leaf> seed = union_sorted(c1[u], c2[v]);
            if (static_cast<int>(seed.size()) > max_seed_leaves) continue;
            DepthWindowRegion region;
            region.seed_leaves = std::move(seed);
            region.common_size = common_size;
            region.spill = spill;
            region.min_depth = std::min(t1.depth[u], t2.depth[v]);
            region.max_depth = std::max(t1.depth[u], t2.depth[v]);
            region.score = 55 * common_size - 8 * spill +
                           static_cast<int>(region.seed_leaves.size()) -
                           3 * std::abs(t1.depth[u] - t2.depth[v]);
            best_for_u.push_back(std::move(region));
        }
        std::sort(best_for_u.begin(), best_for_u.end(),
                  [](const DepthWindowRegion& a, const DepthWindowRegion& b) {
                      if (a.score != b.score) return a.score > b.score;
                      if (a.common_size != b.common_size) {
                          return a.common_size > b.common_size;
                      }
                      if (a.spill != b.spill) return a.spill < b.spill;
                      return a.seed_leaves < b.seed_leaves;
                  });
        int keep_per_node = env_int("RS_SHALLOW_REPLACE_KEEP_PER_NODE", 4);
        for (int i = 0; i < static_cast<int>(best_for_u.size()) &&
                        i < keep_per_node; ++i) {
            add_candidate(std::move(best_for_u[i].seed_leaves),
                          best_for_u[i].common_size, best_for_u[i].spill,
                          best_for_u[i].min_depth, best_for_u[i].max_depth);
        }
    }

    std::sort(out.begin(), out.end(),
              [](const DepthWindowRegion& a, const DepthWindowRegion& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.common_size != b.common_size) {
                      return a.common_size > b.common_size;
                  }
                  if (a.spill != b.spill) return a.spill < b.spill;
                  if (a.seed_leaves.size() != b.seed_leaves.size()) {
                      return a.seed_leaves.size() > b.seed_leaves.size();
                  }
                  return a.seed_leaves < b.seed_leaves;
              });
    if (candidate_limit >= 0 &&
        static_cast<int>(out.size()) > candidate_limit) {
        out.resize(candidate_limit);
    }
    return out;
}

Forest forest_from_islands(const std::vector<std::vector<Leaf>>& islands, int n,
                           const Tree& t1, const Tree& t2) {
    Forest forest;
    forest.leaf_comp.assign(n + 1, -1);
    std::vector<unsigned char> covered(n + 1, 0);
    for (const std::vector<Leaf>& island : islands) {
        Component component;
        component.leaves = island;
        component.repr = component.leaves.front();
        for (Leaf leaf : component.leaves) {
            if (leaf < 1 || leaf > n || covered[leaf]) return Forest{};
            covered[leaf] = 1;
            forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        }
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    for (Leaf leaf = 1; leaf <= n; ++leaf) {
        if (covered[leaf]) continue;
        Component component;
        component.leaves.push_back(leaf);
        component.repr = leaf;
        forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    if (!rebuild_forest_data(forest, t1, t2)) return Forest{};
    if (!validate_forest(forest, t1, t2).ok) return Forest{};
    return forest;
}

std::vector<std::vector<Leaf>> select_islands(
    const std::vector<NearCandidate>& candidates, int n, const Tree& t1,
    const Tree& t2, int max_islands) {
    std::vector<std::vector<Leaf>> selected;
    std::vector<unsigned char> covered(n + 1, 0);
    for (const NearCandidate& candidate : candidates) {
        if (static_cast<int>(selected.size()) >= max_islands) break;
        bool overlaps = false;
        for (Leaf leaf : candidate.leaves) {
            if (covered[leaf]) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;

        std::vector<std::vector<Leaf>> trial_selected = selected;
        trial_selected.push_back(candidate.leaves);
        Forest trial = forest_from_islands(trial_selected, n, t1, t2);
        if (trial.active_count == 0) continue;

        selected.push_back(candidate.leaves);
        for (Leaf leaf : candidate.leaves) covered[leaf] = 1;
    }
    return selected;
}

std::vector<std::vector<Leaf>> select_pseudo_components(
    const std::vector<PseudoCandidate>& candidates, int n, const Tree& t1,
    const Tree& t2, int max_islands) {
    std::vector<std::vector<Leaf>> selected_components;
    std::vector<unsigned char> covered(n + 1, 0);
    int selected_islands = 0;
    for (const PseudoCandidate& candidate : candidates) {
        if (selected_islands >= max_islands) break;
        bool overlaps = false;
        for (Leaf leaf : candidate.leaves) {
            if (covered[leaf]) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;

        std::vector<std::vector<Leaf>> trial = selected_components;
        trial.insert(trial.end(), candidate.components.begin(),
                     candidate.components.end());
        Forest trial_forest = forest_from_islands(trial, n, t1, t2);
        if (trial_forest.active_count == 0) continue;

        selected_components = std::move(trial);
        for (Leaf leaf : candidate.leaves) covered[leaf] = 1;
        ++selected_islands;
    }
    return selected_components;
}

}  // namespace

Forest run_shallow_depth_replace_route(Forest initial, const Tree& t1,
                                       const Tree& t2, int max_leaves,
                                       int min_depth, int max_depth,
                                       int min_seed_leaves,
                                       int max_seed_leaves,
                                       int candidate_limit,
                                       int max_region_components,
                                       int max_region_leaves,
                                       int local_restart_budget,
                                       int local_slack,
                                       int final_restart_budget,
                                       int target_count,
                                       double max_seconds) {
    if (t1.n > max_leaves || min_depth < 0 || max_depth < min_depth ||
        min_seed_leaves < 2 || max_seed_leaves < min_seed_leaves ||
        candidate_limit == 0 || max_region_components <= 1 ||
        max_region_leaves <= 1) {
        return initial;
    }
    const auto start_time = std::chrono::steady_clock::now();
    auto time_left = [&]() {
        if (max_seconds <= 0.0) return true;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        return elapsed < max_seconds;
    };
    auto remaining_seconds = [&]() {
        if (max_seconds <= 0.0) return 0.0;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        return std::max(0.1, max_seconds - elapsed);
    };

    const bool debug = std::getenv("RS_SHALLOW_REPLACE_DEBUG") != nullptr;
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = exact_output_repair(std::move(initial), t1, t2);
    if (!validate_forest(best, t1, t2).ok) return best;

    std::vector<DepthWindowRegion> candidates = generate_depth_window_regions(
        t1, t2, min_depth, max_depth, min_seed_leaves, max_seed_leaves,
        candidate_limit);
    const int filter_common_min =
        env_int("RS_SHALLOW_REPLACE_FILTER_COMMON_MIN", -1);
    const int filter_common_max =
        env_int("RS_SHALLOW_REPLACE_FILTER_COMMON_MAX", INT_MAX);
    const int filter_spill_min =
        env_int("RS_SHALLOW_REPLACE_FILTER_SPILL_MIN", -1);
    const int filter_spill_max =
        env_int("RS_SHALLOW_REPLACE_FILTER_SPILL_MAX", INT_MAX);
    const int filter_region_min =
        env_int("RS_SHALLOW_REPLACE_FILTER_REGION_MIN", -1);
    const int filter_region_max =
        env_int("RS_SHALLOW_REPLACE_FILTER_REGION_MAX", INT_MAX);
    RegionSolveStats stats;
    std::vector<RegionPatch> patches;
    const int patch_pool_limit = env_int("RS_SHALLOW_REPLACE_PATCH_POOL_LIMIT", 96);
    int considered = 0;
    int solved = 0;
    int improved = 0;
    int dump_serial = 0;
    int skipped_region_size = 0;
    int skipped_region_leaves = 0;
    int skipped_edge_closure = 0;

    for (const DepthWindowRegion& candidate : candidates) {
        if (!time_left()) break;
        ++considered;
        std::vector<CompId> region = components_intersecting_leaves(
            best, candidate.seed_leaves, max_region_components);
        if (static_cast<int>(region.size()) < 2) {
            ++skipped_region_size;
            continue;
        }
        if (candidate.common_size < filter_common_min ||
            candidate.common_size > filter_common_max ||
            candidate.spill < filter_spill_min ||
            candidate.spill > filter_spill_max ||
            static_cast<int>(region.size()) < filter_region_min ||
            static_cast<int>(region.size()) > filter_region_max) {
            ++skipped_region_size;
            continue;
        }
        if (!edge_close_region(region, best, t1, t2, max_region_components)) {
            ++skipped_edge_closure;
            continue;
        }
        std::vector<Leaf> region_leaves;
        for (CompId id : region) {
            region_leaves.insert(region_leaves.end(), best.comps[id].leaves.begin(),
                                 best.comps[id].leaves.end());
        }
        std::sort(region_leaves.begin(), region_leaves.end());
        region_leaves.erase(std::unique(region_leaves.begin(), region_leaves.end()),
                            region_leaves.end());
        if (static_cast<int>(region_leaves.size()) > max_region_leaves) {
            ++skipped_region_leaves;
            continue;
        }

        Forest trial;
        if (!solve_region_replacement_portfolio(
                best, region, t1, t2, local_restart_budget, local_slack,
                trial, &stats)) {
            continue;
        }
        ++solved;
        RegionPatch patch;
        if (static_cast<int>(patches.size()) < patch_pool_limit &&
            extract_region_patch(best, trial, region, t1.n, patch)) {
            patches.push_back(std::move(patch));
        }

        for (int mode : {14, 2, 20, 1}) {
            if (!time_left()) break;
            Forest merged = run_h45_greedy_merge(trial, t1, t2, mode, true, true);
            merged = exact_output_repair(std::move(merged), t1, t2);
            if (!validate_forest(merged, t1, t2).ok) continue;
            if (final_restart_budget > 0 && time_left()) {
                Forest restarted = run_h45_split_restarts(
                    merged, t1, t2, mode, final_restart_budget,
                    61000 + considered * 53 + mode * 997, 1, true,
                    target_count, remaining_seconds());
                restarted = exact_output_repair(std::move(restarted), t1, t2);
                if (validate_forest(restarted, t1, t2).ok &&
                    restarted.active_count < merged.active_count) {
                    merged = std::move(restarted);
                }
            }
            if (merged.active_count < best.active_count) {
                best = std::move(merged);
                ++improved;
                maybe_dump_forest("RS_SHALLOW_REPLACE_DUMP_PREFIX", best, t1,
                                  ++dump_serial);
                if (debug) {
                    std::cerr << "rs_shallow_replace_improve components="
                              << best.active_count
                              << " region_components=" << region.size()
                              << " region_leaves=" << region_leaves.size()
                              << " seed_leaves=" << candidate.seed_leaves.size()
                              << " common=" << candidate.common_size
                              << " spill=" << candidate.spill
                              << " depth=" << candidate.min_depth << "-"
                              << candidate.max_depth
                              << " mode=" << mode << '\n';
                }
                if (target_count > 0 && best.active_count <= target_count) {
                    return best;
                }
                break;
            }
        }
    }

    if (env_int("RS_SHALLOW_REPLACE_PATCH_COMBO_ROUTE", 1) != 0 &&
        patches.size() >= 2 && time_left()) {
        std::sort(patches.begin(), patches.end(),
                  [](const RegionPatch& a, const RegionPatch& b) {
                      if (a.score != b.score) return a.score > b.score;
                      if (a.replacement_count != b.replacement_count) {
                          return a.replacement_count < b.replacement_count;
                      }
                      return a.region < b.region;
                  });
        int combo_limit = env_int("RS_SHALLOW_REPLACE_PATCH_COMBO_LIMIT", 1200);
        int combo_considered = 0;
        for (int i = 0; i < static_cast<int>(patches.size()) &&
                        combo_considered < combo_limit && time_left(); ++i) {
            for (int j = i + 1; j < static_cast<int>(patches.size()) &&
                            combo_considered < combo_limit && time_left(); ++j) {
                if (!patches_disjoint(patches[i], patches[j])) continue;
                ++combo_considered;
                Forest combo;
                std::vector<const RegionPatch*> selected = {&patches[i], &patches[j]};
                if (!apply_region_patches(best, selected, t1, t2, combo)) continue;
                for (int mode : {14, 2, 20}) {
                    if (!time_left()) break;
                    Forest merged =
                        run_h45_greedy_merge(combo, t1, t2, mode, true, true);
                    merged = exact_output_repair(std::move(merged), t1, t2);
                    if (!validate_forest(merged, t1, t2).ok) continue;
                    if (final_restart_budget > 0 && time_left()) {
                        Forest restarted = run_h45_split_restarts(
                            merged, t1, t2, mode, final_restart_budget,
                            71000 + combo_considered * 131 + mode, 1, true,
                            target_count, remaining_seconds());
                        restarted = exact_output_repair(std::move(restarted), t1, t2);
                        if (validate_forest(restarted, t1, t2).ok &&
                            restarted.active_count < merged.active_count) {
                            merged = std::move(restarted);
                        }
                    }
                    if (merged.active_count < best.active_count) {
                        best = std::move(merged);
                        ++improved;
                        maybe_dump_forest("RS_SHALLOW_REPLACE_DUMP_PREFIX", best,
                                          t1, ++dump_serial);
                        if (debug) {
                            std::cerr
                                << "rs_shallow_replace_combo_improve components="
                                << best.active_count
                                << " combo_considered=" << combo_considered
                                << " mode=" << mode << '\n';
                        }
                        if (target_count > 0 && best.active_count <= target_count) {
                            return best;
                        }
                        break;
                    }
                }
            }
        }
    }

    if (debug) {
        std::cerr << "rs_shallow_replace_done components=" << best.active_count
                  << " candidates=" << candidates.size()
                  << " considered=" << considered
                  << " solved=" << solved
                  << " improved=" << improved
                  << " patches=" << patches.size()
                  << " skip_region_size=" << skipped_region_size
                  << " skip_region_leaves=" << skipped_region_leaves
                  << " skip_edge_closure=" << skipped_edge_closure
                  << " bad_region=" << stats.bad_region
                  << " parse_fail=" << stats.parse_fail
                  << " local_rebuild_fail=" << stats.local_rebuild_fail
                  << " local_invalid=" << stats.local_invalid
                  << " slack_reject=" << stats.slack_reject
                  << " duplicate_leaf=" << stats.duplicate_leaf
                  << " coverage_fail=" << stats.coverage_fail
                  << " global_rebuild_fail=" << stats.global_rebuild_fail
                  << " global_invalid=" << stats.global_invalid
                  << " accepted=" << stats.accepted << '\n';
    }
    return best;
}

Forest run_near_cluster_route(Forest initial, const Tree& t1, const Tree& t2,
                              int max_leaves, int min_island_size,
                              int max_spill, int max_trim,
                              int max_candidates, int max_islands,
                              int restart_budget, int target_count) {
    if (t1.n > max_leaves || min_island_size < 2 || max_spill < 0 ||
        max_candidates < 0 || max_islands < 0) {
        return initial;
    }
    const bool debug = std::getenv("RS_NEAR_CLUSTER_DEBUG") != nullptr;
    if (!validate_forest(initial, t1, t2).ok) return initial;
    Forest best = std::move(initial);

    if (env_int("RS_NEAR_CLADE_POLLUTER_ROUTE", 1) != 0) {
        best = run_near_clade_polluter_route(
            std::move(best), t1, t2,
            env_int("RS_NEAR_CLADE_POLLUTER_MAX_LEAVES", max_leaves),
            env_int("RS_NEAR_CLADE_POLLUTER_MIN_COMMON", min_island_size),
            env_int("RS_NEAR_CLADE_POLLUTER_MAX_SPILL", max_spill),
            env_int("RS_NEAR_CLADE_POLLUTER_CANDIDATES", 600),
            env_int("RS_NEAR_CLADE_POLLUTER_MAX_REGION_COMPONENTS", 18),
            env_int("RS_NEAR_CLADE_POLLUTER_MAX_REGION_LEAVES", 75),
            env_int("RS_NEAR_CLADE_POLLUTER_LOCAL_RESTART_BUDGET", 32),
            env_int("RS_NEAR_CLADE_POLLUTER_LOCAL_SLACK", 3),
            env_int("RS_NEAR_CLADE_POLLUTER_FINAL_RESTART_BUDGET", 64),
            target_count);
        if (target_count > 0 && best.active_count <= target_count) {
            return best;
        }
    }

    if (env_int("RS_CONFLICT_POLLUTER_ROUTE", 1) != 0) {
        best = run_conflict_polluter_route(
            std::move(best), t1, t2,
            env_int("RS_CONFLICT_POLLUTER_MAX_LEAVES", max_leaves),
            env_int("RS_CONFLICT_POLLUTER_PAIR_LIMIT", 80),
            env_int("RS_CONFLICT_POLLUTER_MAX_REGION_COMPONENTS", 14),
            env_int("RS_CONFLICT_POLLUTER_MAX_REGION_LEAVES", 60),
            env_int("RS_CONFLICT_POLLUTER_LOCAL_RESTART_BUDGET", 32),
            env_int("RS_CONFLICT_POLLUTER_LOCAL_SLACK", 2),
            env_int("RS_CONFLICT_POLLUTER_FINAL_RESTART_BUDGET", 64),
            target_count);
        if (target_count > 0 && best.active_count <= target_count) {
            return best;
        }
    }

    if (env_int("RS_NEAR_CLUSTER_SEED_ROUTE", 1) == 0 ||
        max_candidates == 0 || max_islands == 0) {
        return best;
    }

    std::vector<PseudoCandidate> pseudo_candidates = generate_pseudo_candidates(
        t1, t2, min_island_size, max_spill,
        env_int("RS_NEAR_CLUSTER_RAW_CANDIDATES", 80),
        env_int("RS_NEAR_CLUSTER_LOCAL_RESTART_BUDGET", 64));
    std::vector<std::vector<Leaf>> islands =
        select_pseudo_components(pseudo_candidates, t1.n, t1, t2, max_islands);

    std::vector<NearCandidate> candidates =
        generate_near_candidates(t1, t2, min_island_size, max_spill, max_trim,
                                 max_candidates);
    if (islands.empty()) {
        islands = select_islands(candidates, t1.n, t1, t2, max_islands);
    }
    if (debug) {
        int island_leaves = 0;
        for (const auto& island : islands) {
            island_leaves += static_cast<int>(island.size());
        }
        std::cerr << "rs_near_cluster_candidates count=" << candidates.size()
                  << " pseudo=" << pseudo_candidates.size()
                  << " selected=" << islands.size()
                  << " selected_leaves=" << island_leaves << '\n';
    }
    if (islands.empty()) return best;

    Forest seed = forest_from_islands(islands, t1.n, t1, t2);
    if (seed.active_count == 0) return best;

    const int modes[] = {2, 14, 1, 5};
    for (int mode : modes) {
        Forest trial_seed = seed;
        Forest trial =
            run_h45_greedy_merge(std::move(trial_seed), t1, t2, mode, true, true);
        trial = exact_output_repair(std::move(trial), t1, t2);
        if (validate_forest(trial, t1, t2).ok &&
            trial.active_count < best.active_count) {
            best = std::move(trial);
            if (debug) {
                std::cerr << "rs_near_cluster_greedy_improve mode=" << mode
                          << " components=" << best.active_count << '\n';
            }
            if (target_count > 0 && best.active_count <= target_count) return best;
        }
        if (restart_budget > 0) {
            Forest restarted = run_h45_split_restarts(
                best, t1, t2, mode, restart_budget, 17000 + mode * 313, 1,
                false, target_count, 0.0);
            restarted = exact_output_repair(std::move(restarted), t1, t2);
            if (validate_forest(restarted, t1, t2).ok &&
                restarted.active_count < best.active_count) {
                best = std::move(restarted);
                if (debug) {
                    std::cerr << "rs_near_cluster_restart_improve mode=" << mode
                              << " components=" << best.active_count << '\n';
                }
                if (target_count > 0 && best.active_count <= target_count) {
                    return best;
                }
            }
        }
    }
    if (debug) {
        std::cerr << "rs_near_cluster_done components=" << best.active_count << '\n';
    }
    return best;
}
