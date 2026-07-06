#include "heuristic/greedy_merge.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <unordered_set>

#include "forest/validator.h"
#include "heuristic/candidates.h"
#include "tree/induced.h"
#include "util/signal_handler.h"

namespace {

int h45_env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

struct MergeTest {
    bool ok = false;
    std::vector<Leaf> leaves;
    InducedInfo i1;
    InducedInfo i2;
};

std::vector<Leaf> merged_leaves_of(const Component& a, const Component& b) {
    std::vector<Leaf> merged;
    merged.reserve(a.leaves.size() + b.leaves.size());
    std::merge(a.leaves.begin(), a.leaves.end(), b.leaves.begin(), b.leaves.end(),
               std::back_inserter(merged));
    return merged;
}

bool edges_available(const std::vector<NodeId>& edges, const std::vector<CompId>& owner,
                     CompId a, CompId b) {
    for (NodeId edge : edges) {
        CompId current = owner[edge];
        if (current != -1 && current != a && current != b) return false;
    }
    return true;
}

MergeTest test_merge(const Forest& forest, const Tree& t1, const Tree& t2,
                     CompId a, CompId b) {
    MergeTest test;
    test.leaves = merged_leaves_of(forest.comps[a], forest.comps[b]);
    test.i1 = compute_induced(t1, test.leaves);
    test.i2 = compute_induced(t2, test.leaves);
    if (test.i1.signature != test.i2.signature) return test;
    if (!edges_available(test.i1.used_edges, forest.owner1, a, b)) return test;
    if (!edges_available(test.i2.used_edges, forest.owner2, a, b)) return test;
    test.ok = true;
    return test;
}

void clear_edges(const std::vector<NodeId>& edges, std::vector<CompId>& owner) {
    for (NodeId edge : edges) owner[edge] = -1;
}

void claim_edges(const std::vector<NodeId>& edges, std::vector<CompId>& owner, CompId id) {
    for (NodeId edge : edges) owner[edge] = id;
}

void apply_merge(Forest& forest, CompId a, CompId b, MergeTest test) {
    clear_edges(forest.comps[a].edges1, forest.owner1);
    clear_edges(forest.comps[a].edges2, forest.owner2);
    clear_edges(forest.comps[b].edges1, forest.owner1);
    clear_edges(forest.comps[b].edges2, forest.owner2);

    Component& kept = forest.comps[a];
    Leaf repr = kept.repr;
    kept.leaves = std::move(test.leaves);
    kept.repr = repr > 0 ? repr : kept.leaves.front();
    kept.root1 = test.i1.root;
    kept.root2 = test.i2.root;
    kept.signature1 = std::move(test.i1.signature);
    kept.signature2 = std::move(test.i2.signature);
    kept.edges1 = std::move(test.i1.used_edges);
    kept.edges2 = std::move(test.i2.used_edges);

    forest.comps[b].active = false;
    forest.comps[b].edges1.clear();
    forest.comps[b].edges2.clear();
    for (Leaf x : kept.leaves) forest.leaf_comp[x] = a;
    claim_edges(kept.edges1, forest.owner1, a);
    claim_edges(kept.edges2, forest.owner2, a);
    --forest.active_count;
}

CompId append_exact_component(Forest& forest, std::vector<Leaf> leaves, Leaf repr,
                              const Tree& t1, const Tree& t2) {
    if (leaves.empty()) return -1;

    InducedInfo i1 = compute_induced(t1, leaves);
    InducedInfo i2 = compute_induced(t2, leaves);
    if (i1.signature != i2.signature) return -1;
    if (!edges_available(i1.used_edges, forest.owner1, -2, -2)) return -1;
    if (!edges_available(i2.used_edges, forest.owner2, -2, -2)) return -1;

    CompId id = static_cast<CompId>(forest.comps.size());
    Component component;
    component.leaves = std::move(leaves);
    component.repr = repr > 0 ? repr : component.leaves.front();
    component.root1 = i1.root;
    component.root2 = i2.root;
    component.signature1 = std::move(i1.signature);
    component.signature2 = std::move(i2.signature);
    component.edges1 = std::move(i1.used_edges);
    component.edges2 = std::move(i2.used_edges);
    for (Leaf x : component.leaves) forest.leaf_comp[x] = id;
    claim_edges(component.edges1, forest.owner1, id);
    claim_edges(component.edges2, forest.owner2, id);
    forest.comps.push_back(std::move(component));
    ++forest.active_count;
    return id;
}

using H45Hash = uint64_t;

H45Hash h45_combine_hashes(std::vector<H45Hash>& values) {
    std::sort(values.begin(), values.end());
    H45Hash hash = 1469598103934665603ULL;
    for (H45Hash value : values) {
        hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return hash;
}

H45Hash h45_combine_two_hashes(H45Hash a, H45Hash b) {
    if (b < a) std::swap(a, b);
    H45Hash hash = 1469598103934665603ULL;
    hash ^= a + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    hash ^= b + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    return hash;
}

H45Hash h45_canonical_restriction_hash_recursive(const Tree& tree, NodeId node,
                                                 const std::vector<int>& mark,
                                                 int token) {
    if (tree.is_leaf(node)) {
        Leaf leaf = tree.label[node];
        if (leaf > 0 && leaf < static_cast<int>(mark.size()) &&
            mark[leaf] == token) {
            return static_cast<H45Hash>(leaf) * 1315423911ULL;
        }
        return 0;
    }

    H45Hash left = tree.left[node] >= 0
                       ? h45_canonical_restriction_hash_recursive(
                             tree, tree.left[node], mark, token)
                       : 0;
    H45Hash right = tree.right[node] >= 0
                        ? h45_canonical_restriction_hash_recursive(
                              tree, tree.right[node], mark, token)
                        : 0;
    if (left == 0) return right;
    if (right == 0) return left;
    return h45_combine_two_hashes(left, right);
}

H45Hash h45_canonical_restriction_hash(const Tree& tree, NodeId root,
                                       const std::vector<int>& mark, int token) {
    if (root < 0 || root >= tree.node_count()) return 0;
    if (tree.node_count() < 10000) {
        return h45_canonical_restriction_hash_recursive(tree, root, mark, token);
    }
    struct Frame {
        NodeId node;
        bool exiting;
    };
    std::vector<Frame> stack;
    std::vector<H45Hash> value(tree.node_count(), 0);
    stack.push_back({root, false});
    while (!stack.empty()) {
        Frame frame = stack.back();
        stack.pop_back();
        NodeId node = frame.node;
        if (frame.exiting) {
            if (tree.is_leaf(node)) {
                Leaf leaf = tree.label[node];
                if (leaf > 0 && leaf < static_cast<int>(mark.size()) &&
                    mark[leaf] == token) {
                    value[node] = static_cast<H45Hash>(leaf) * 1315423911ULL;
                } else {
                    value[node] = 0;
                }
                continue;
            }
            std::vector<H45Hash> children;
            if (tree.left[node] >= 0 && value[tree.left[node]] != 0) {
                children.push_back(value[tree.left[node]]);
            }
            if (tree.right[node] >= 0 && value[tree.right[node]] != 0) {
                children.push_back(value[tree.right[node]]);
            }
            if (children.empty()) value[node] = 0;
            else if (children.size() == 1) value[node] = children.front();
            else value[node] = h45_combine_hashes(children);
            continue;
        }
        stack.push_back({node, true});
        if (!tree.is_leaf(node)) {
            if (tree.right[node] >= 0) stack.push_back({tree.right[node], false});
            if (tree.left[node] >= 0) stack.push_back({tree.left[node], false});
        }
    }
    return value[root];
}

bool h45_component_topology_valid(const std::vector<Leaf>& leaves,
                                  const Tree& t1, const Tree& t2,
                                  NodeId root1, NodeId root2) {
    if (leaves.empty() || root1 < 0 || root2 < 0) return false;
    std::vector<int> mark(t1.n + 1, 0);
    for (Leaf leaf : leaves) {
        if (leaf <= 0 || leaf > t1.n) return false;
        mark[leaf] = 1;
    }
    H45Hash h1 = h45_canonical_restriction_hash(t1, root1, mark, 1);
    H45Hash h2 = h45_canonical_restriction_hash(t2, root2, mark, 1);
    return h1 != 0 && h1 == h2;
}

bool h45_paths_allowed(const Tree& tree, const std::vector<CompId>& owner,
                       const std::vector<Leaf>& leaves, NodeId root,
                       CompId allowed_a, CompId allowed_b) {
    for (Leaf leaf : leaves) {
        NodeId node = tree.leaf_node[leaf];
        while (node >= 0 && node != root) {
            CompId current = owner[node];
            if (current != -1 && current != allowed_a && current != allowed_b) {
                return false;
            }
            node = tree.parent[node];
        }
    }
    return true;
}

std::vector<NodeId> h45_mark_paths(const Tree& tree, std::vector<CompId>& owner,
                                   const std::vector<Leaf>& leaves, NodeId root,
                                   CompId id) {
    std::vector<NodeId> edges;
    for (Leaf leaf : leaves) {
        NodeId node = tree.leaf_node[leaf];
        while (node >= 0 && node != root) {
            owner[node] = id;
            edges.push_back(node);
            node = tree.parent[node];
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

CompId h45_try_merge(Forest& forest, const Tree& t1, const Tree& t2,
                     CompId a, CompId b) {
    if (a < 0 || b < 0 || a == b) return -1;
    if (a >= static_cast<CompId>(forest.comps.size()) ||
        b >= static_cast<CompId>(forest.comps.size())) {
        return -1;
    }
    if (!forest.comps[a].active || !forest.comps[b].active) return -1;

    std::vector<Leaf> merged = forest.comps[a].leaves;
    merged.insert(merged.end(), forest.comps[b].leaves.begin(),
                  forest.comps[b].leaves.end());
    NodeId root1 = -1;
    NodeId root2 = -1;
    std::vector<NodeId> edges1;
    std::vector<NodeId> edges2;
    if (std::getenv("RS_H45_EXACT_MERGE") != nullptr) {
        InducedInfo i1 = compute_induced(t1, merged);
        InducedInfo i2 = compute_induced(t2, merged);
        if (i1.signature != i2.signature) return -1;
        if (!edges_available(i1.used_edges, forest.owner1, a, b)) return -1;
        if (!edges_available(i2.used_edges, forest.owner2, a, b)) return -1;
        root1 = i1.root;
        root2 = i2.root;
        edges1 = std::move(i1.used_edges);
        edges2 = std::move(i2.used_edges);
    } else {
        root1 = t1.lca_leaves(merged);
        root2 = t2.lca_leaves(merged);
        if (!h45_component_topology_valid(merged, t1, t2, root1, root2)) return -1;
        if (!h45_paths_allowed(t1, forest.owner1, merged, root1, a, b)) return -1;
        if (!h45_paths_allowed(t2, forest.owner2, merged, root2, a, b)) return -1;
    }

    Leaf repr = forest.comps[a].repr > 0 ? forest.comps[a].repr
                                         : forest.comps[a].leaves.front();
    Component component;
    component.leaves = std::move(merged);
    component.repr = repr;
    component.root1 = root1;
    component.root2 = root2;

    forest.comps[a].active = false;
    forest.comps[b].active = false;
    forest.comps[a].edges1.clear();
    forest.comps[a].edges2.clear();
    forest.comps[b].edges1.clear();
    forest.comps[b].edges2.clear();

    CompId id = static_cast<CompId>(forest.comps.size());
    for (Leaf x : component.leaves) forest.leaf_comp[x] = id;
    if (std::getenv("RS_H45_EXACT_MERGE") != nullptr) {
        component.edges1 = std::move(edges1);
        component.edges2 = std::move(edges2);
        claim_edges(component.edges1, forest.owner1, id);
        claim_edges(component.edges2, forest.owner2, id);
    } else {
        component.edges1 = h45_mark_paths(t1, forest.owner1, component.leaves, root1, id);
        component.edges2 = h45_mark_paths(t2, forest.owner2, component.leaves, root2, id);
    }
    forest.comps.push_back(std::move(component));
    if (std::getenv("RS_TRACE_N6_MERGES") != nullptr && t1.n == 6) {
        const Component& merged_component = forest.comps[id];
        std::cerr << "rs_trace_merge n=6 A=" << a
                  << " B=" << b
                  << " C=" << id
                  << " reprA=" << merged_component.repr
                  << " leaves=";
        for (Leaf leaf : merged_component.leaves) std::cerr << leaf << ',';
        std::cerr << '\n';
    }
    --forest.active_count;
    return id;
}

Forest split_component_to_singletons(Forest forest, CompId id,
                                     const Tree& t1, const Tree& t2) {
    Component& old = forest.comps[id];
    std::vector<Leaf> leaves = old.leaves;
    if (!old.active || leaves.size() < 2) return forest;

    old.active = false;
    old.edges1.clear();
    old.edges2.clear();
    --forest.active_count;

    for (Leaf x : leaves) {
        Component c;
        c.leaves.push_back(x);
        c.repr = x;
        forest.leaf_comp[x] = static_cast<CompId>(forest.comps.size());
        forest.comps.push_back(std::move(c));
        ++forest.active_count;
    }
    rebuild_forest_data(forest, t1, t2);
    return forest;
}

std::vector<CompId> largest_components(const Forest& forest, int limit) {
    std::vector<CompId> ids = active_components(forest);
    std::sort(ids.begin(), ids.end(), [&](CompId a, CompId b) {
        size_t sa = forest.comps[a].leaves.size();
        size_t sb = forest.comps[b].leaves.size();
        if (sa != sb) return sa > sb;
        return a < b;
    });
    ids.erase(std::remove_if(ids.begin(), ids.end(), [&](CompId id) {
        return forest.comps[id].leaves.size() < 2;
    }), ids.end());
    if (static_cast<int>(ids.size()) > limit) ids.resize(limit);
    return ids;
}

std::vector<std::pair<CompId, CompId>> largest_component_pairs(const Forest& forest,
                                                               int top_components,
                                                               int limit) {
    std::vector<CompId> ids = largest_components(forest, top_components);
    std::vector<std::pair<CompId, CompId>> pairs;
    for (size_t i = 0; i < ids.size(); ++i) {
        for (size_t j = i + 1; j < ids.size(); ++j) {
            pairs.push_back({ids[i], ids[j]});
        }
    }
    std::sort(pairs.begin(), pairs.end(), [&](const auto& x, const auto& y) {
        size_t sx = forest.comps[x.first].leaves.size() + forest.comps[x.second].leaves.size();
        size_t sy = forest.comps[y.first].leaves.size() + forest.comps[y.second].leaves.size();
        if (sx != sy) return sx > sy;
        return x < y;
    });
    if (static_cast<int>(pairs.size()) > limit) pairs.resize(limit);
    return pairs;
}

std::vector<CompId> unique_sorted_components(std::vector<CompId> ids) {
    ids.erase(std::remove_if(ids.begin(), ids.end(), [](CompId id) {
        return id < 0;
    }), ids.end());
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

Forest split_components_to_singletons(Forest forest, const std::vector<CompId>& ids,
                                      const Tree& t1, const Tree& t2) {
    std::vector<CompId> split_ids = unique_sorted_components(ids);
    if (split_ids.empty()) return forest;

    for (CompId id : split_ids) {
        if (id < 0 || id >= static_cast<CompId>(forest.comps.size())) continue;
        Component& old = forest.comps[id];
        if (!old.active || old.leaves.size() < 2) continue;
        std::vector<Leaf> leaves = old.leaves;
        old.active = false;
        old.edges1.clear();
        old.edges2.clear();
        --forest.active_count;

        for (Leaf leaf : leaves) {
            Component component;
            component.leaves.push_back(leaf);
            component.repr = leaf;
            forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
            forest.comps.push_back(std::move(component));
            ++forest.active_count;
        }
    }
    rebuild_forest_data(forest, t1, t2);
    return forest;
}

std::vector<CompId> edge_blockers(const std::vector<NodeId>& edges,
                                  const std::vector<CompId>& owner,
                                  CompId a, CompId b) {
    std::vector<CompId> blockers;
    for (NodeId edge : edges) {
        CompId id = owner[edge];
        if (id != -1 && id != a && id != b) blockers.push_back(id);
    }
    return unique_sorted_components(std::move(blockers));
}

struct ConflictTarget {
    std::vector<CompId> split;
    int score = 0;
};

void add_conflict_target(std::vector<ConflictTarget>& targets,
                         std::vector<CompId> split,
                         const Forest& forest,
                         int base_score) {
    split = unique_sorted_components(std::move(split));
    if (split.empty()) return;
    int leaf_count = 0;
    for (CompId id : split) {
        if (id >= 0 && id < static_cast<CompId>(forest.comps.size()) &&
            forest.comps[id].active) {
            leaf_count += static_cast<int>(forest.comps[id].leaves.size());
        }
    }
    if (leaf_count == 0) return;
    targets.push_back({std::move(split), base_score + leaf_count});
}


CompId subtree_component_plain(const Tree& tree, const Forest& forest, NodeId v,
                               std::vector<std::pair<CompId, CompId>>& pairs) {
    if (tree.is_leaf(v)) return forest.leaf_comp[tree.label[v]];

    CompId a = subtree_component_plain(tree, forest, tree.left[v], pairs);
    CompId b = subtree_component_plain(tree, forest, tree.right[v], pairs);
    if (a == b) return a;
    if (a >= 0 && b >= 0 && forest.comps[a].active && forest.comps[b].active) {
        if (b < a) std::swap(a, b);
        pairs.push_back({a, b});
    }
    return -2;
}

std::vector<std::pair<CompId, CompId>> sibling_component_pairs(const Forest& forest,
                                                               const Tree& tree) {
    std::vector<std::pair<CompId, CompId>> pairs;
    subtree_component_plain(tree, forest, tree.root, pairs);
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

bool contract_common_cherries(Forest& forest, const Tree& t1, const Tree& t2) {
    bool any = false;
    while (true) {
        std::vector<std::pair<CompId, CompId>> p1 = sibling_component_pairs(forest, t1);
        std::vector<std::pair<CompId, CompId>> p2 = sibling_component_pairs(forest, t2);
        std::sort(p2.begin(), p2.end());
        bool changed = false;
        for (const auto& pair : p1) {
            if (!std::binary_search(p2.begin(), p2.end(), pair)) continue;
            if (!forest.comps[pair.first].active || !forest.comps[pair.second].active) {
                continue;
            }
            MergeTest test = test_merge(forest, t1, t2, pair.first, pair.second);
            if (!test.ok) continue;
            apply_merge(forest, pair.first, pair.second, std::move(test));
            any = true;
            changed = true;
            break;
        }
        if (!changed) break;
    }
    return any;
}

void collect_conflict_targets_from_pairs(const Forest& forest,
                                         const Tree& source_tree,
                                         const Tree& other_tree,
                                         const std::vector<CompId>& other_owner,
                                         bool other_is_t1,
                                         std::vector<ConflictTarget>& targets) {
    for (const auto& pair : sibling_component_pairs(forest, source_tree)) {
        CompId a = pair.first;
        CompId b = pair.second;
        std::vector<Leaf> merged = merged_leaves_of(forest.comps[a], forest.comps[b]);
        InducedInfo source = compute_induced(source_tree, merged);
        InducedInfo other = compute_induced(other_tree, merged);
        if (source.signature == other.signature &&
            edges_available(other.used_edges, other_owner, a, b)) {
            continue;
        }

        std::vector<CompId> blockers = edge_blockers(other.used_edges, other_owner, a, b);
        std::sort(blockers.begin(), blockers.end(), [&](CompId x, CompId y) {
            size_t sx = forest.comps[x].leaves.size();
            size_t sy = forest.comps[y].leaves.size();
            if (sx != sy) return sx < sy;
            return x < y;
        });

        int base = source.signature == other.signature ? 1800 : 1200;
        std::vector<CompId> split;
        for (CompId id : blockers) {
            if (static_cast<int>(split.size()) >= 3) break;
            split.push_back(id);
        }
        if (!split.empty()) add_conflict_target(targets, split, forest, base);

        split.push_back(a);
        add_conflict_target(targets, split, forest, base - 150);
        split.back() = b;
        add_conflict_target(targets, split, forest, base - 150);

        if (blockers.empty()) {
            add_conflict_target(targets, {a}, forest, 700);
            add_conflict_target(targets, {b}, forest, 700);
        }

        (void)other_is_t1;
    }
}

std::vector<ConflictTarget> conflict_targets(const Forest& forest, const Tree& t1,
                                             const Tree& t2, int limit,
                                             int max_split_components) {
    std::vector<ConflictTarget> targets;
    collect_conflict_targets_from_pairs(forest, t1, t2, forest.owner2, false, targets);
    collect_conflict_targets_from_pairs(forest, t2, t1, forest.owner1, true, targets);

    std::map<std::vector<CompId>, int> best_scores;
    for (ConflictTarget& target : targets) {
        target.split = unique_sorted_components(std::move(target.split));
        if (target.split.empty() ||
            static_cast<int>(target.split.size()) > max_split_components) {
            continue;
        }
        auto it = best_scores.find(target.split);
        if (it == best_scores.end() || target.score > it->second) {
            best_scores[target.split] = target.score;
        }
    }

    std::vector<ConflictTarget> deduped;
    deduped.reserve(best_scores.size());
    for (const auto& item : best_scores) deduped.push_back({item.first, item.second});
    std::sort(deduped.begin(), deduped.end(), [&](const auto& x, const auto& y) {
        if (x.score != y.score) return x.score > y.score;
        if (x.split.size() != y.split.size()) return x.split.size() < y.split.size();
        return x.split < y.split;
    });
    if (static_cast<int>(deduped.size()) > limit) deduped.resize(limit);
    return deduped;
}

Forest run_one_route(Forest forest, const Tree& t1, const Tree& t2,
                     const GreedyConfig& config, CandidateMode mode,
                     int jitter_seed, int jitter_span) {
    if (config.enable_common_cherry_contraction) {
        contract_common_cherries(forest, t1, t2);
    }
    while (true) {
        bool changed = false;

        for (const MergeCandidate& candidate :
             generate_candidates(forest, t1, t2, config.all_pair_cap, mode,
                                 jitter_seed, jitter_span)) {
            MergeTest test = test_merge(forest, t1, t2, candidate.a, candidate.b);
            if (!test.ok) continue;
            apply_merge(forest, candidate.a, candidate.b, std::move(test));
            if (config.enable_common_cherry_contraction) {
                contract_common_cherries(forest, t1, t2);
            }
            changed = true;
            break;
        }

        if (!changed) break;
    }
    return forest;
}

Forest repair_by_single_splits(Forest best, const Tree& t1, const Tree& t2,
                               const GreedyConfig& config, bool debug) {
    struct Route {
        CandidateMode mode;
        int seed;
        int jitter;
    };
    const Route repair_routes[] = {
        {CandidateMode::WindowWide, 701, 350},
        {CandidateMode::WindowWide, 809, 350},
        {CandidateMode::WindowScored, 911, 350},
    };

    for (int pass = 0; pass < config.repair_passes; ++pass) {
        bool improved = false;
        std::vector<CompId> targets = largest_components(best, config.repair_top_components);
        if (debug) {
            std::cerr << "repair_start pass=" << pass
                      << " targets=" << targets.size()
                      << " best=" << best.active_count << "\n";
        }
        for (CompId id : targets) {
            Forest split = split_component_to_singletons(best, id, t1, t2);
            for (const Route& route : repair_routes) {
                Forest candidate = run_one_route(split, t1, t2, config, route.mode,
                                                 route.seed + pass * 1009, route.jitter);
                bool ok = validate_forest(candidate, t1, t2).ok;
                if (debug) {
                    std::cerr << "repair pass=" << pass
                              << " split=" << id
                              << " route=" << candidate_mode_name(route.mode)
                              << " seed=" << route.seed + pass * 1009
                              << " components=" << candidate.active_count
                              << " valid=" << ok << "\n";
                }
                if (ok && candidate.active_count < best.active_count) {
                    best = std::move(candidate);
                    if (config.publish_improvements) {
                        publish_best_output(forest_to_output(best, t1));
                    }
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
        if (!improved) break;
    }
    return best;
}

Forest repair_by_pair_splits(Forest best, const Tree& t1, const Tree& t2,
                             const GreedyConfig& config, bool debug) {
    struct Route {
        CandidateMode mode;
        int seed;
        int jitter;
    };
    const Route repair_routes[] = {
        {CandidateMode::WindowWide, 1201, 450},
        {CandidateMode::WindowScored, 1303, 450},
    };

    for (int pass = 0; pass < config.repair_passes; ++pass) {
        bool improved = false;
        std::vector<std::pair<CompId, CompId>> pairs =
            largest_component_pairs(best, config.repair_pair_top_components,
                                    config.repair_pair_limit);
        if (debug) {
            std::cerr << "pair_repair_start pass=" << pass
                      << " pairs=" << pairs.size()
                      << " best=" << best.active_count << "\n";
        }
        for (const auto& pair : pairs) {
            Forest split = split_component_to_singletons(best, pair.first, t1, t2);
            split = split_component_to_singletons(std::move(split), pair.second, t1, t2);
            for (const Route& route : repair_routes) {
                Forest candidate = run_one_route(split, t1, t2, config, route.mode,
                                                 route.seed + pass * 2003, route.jitter);
                bool ok = validate_forest(candidate, t1, t2).ok;
                if (debug) {
                    std::cerr << "pair_repair pass=" << pass
                              << " split=" << pair.first << "," << pair.second
                              << " route=" << candidate_mode_name(route.mode)
                              << " seed=" << route.seed + pass * 2003
                              << " components=" << candidate.active_count
                              << " valid=" << ok << "\n";
                }
                if (ok && candidate.active_count < best.active_count) {
                    best = std::move(candidate);
                    if (config.publish_improvements) {
                        publish_best_output(forest_to_output(best, t1));
                    }
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
        if (!improved) break;
    }
    return best;
}

Forest repair_by_conflict_targets(Forest best, const Tree& t1, const Tree& t2,
                                  const GreedyConfig& config, bool debug) {
    struct Route {
        CandidateMode mode;
        int seed;
        int jitter;
    };
    const Route repair_routes[] = {
        {CandidateMode::WindowWide, 2201, 500},
        {CandidateMode::WindowScored, 2303, 500},
        {CandidateMode::WindowTight, 2407, 350},
    };

    for (int pass = 0; pass < config.repair_passes; ++pass) {
        bool improved = false;
        std::vector<ConflictTarget> targets =
            conflict_targets(best, t1, t2, config.conflict_target_limit,
                             config.conflict_max_split_components);
        if (debug) {
            std::cerr << "conflict_repair_start pass=" << pass
                      << " targets=" << targets.size()
                      << " best=" << best.active_count << "\n";
        }
        for (const ConflictTarget& target : targets) {
            Forest split = split_components_to_singletons(best, target.split, t1, t2);
            for (const Route& route : repair_routes) {
                Forest candidate = run_one_route(split, t1, t2, config, route.mode,
                                                 route.seed + pass * 3001, route.jitter);
                bool ok = validate_forest(candidate, t1, t2).ok;
                if (debug) {
                    std::cerr << "conflict_repair pass=" << pass
                              << " split_size=" << target.split.size()
                              << " score=" << target.score
                              << " route=" << candidate_mode_name(route.mode)
                              << " seed=" << route.seed + pass * 3001
                              << " components=" << candidate.active_count
                              << " valid=" << ok << "\n";
                }
                if (ok && candidate.active_count < best.active_count) {
                    best = std::move(candidate);
                    if (config.publish_improvements) {
                        publish_best_output(forest_to_output(best, t1));
                    }
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
        if (!improved) break;
    }
    return best;
}

constexpr int H45_LOCAL_COMPONENT_WINDOW = 24;
constexpr int H45_NEIGHBOR_WALK_LIMIT = 20;
constexpr int H45_RANDOM_CANDIDATE_MULTIPLIER = 2;
constexpr int H45_MIN_MERGE_CHECKS = 10000;
constexpr int H45_MERGE_CHECK_MULTIPLIER = 80;
constexpr int H45_JITTER_CYCLE_LENGTH = 24;
constexpr int H45_CONFLICT_RESTART_INTERVAL = 2;
constexpr int H45_LOOKAHEAD_RESTART_INTERVAL = 16;
constexpr int H45_LOOKAHEAD_TOP_CANDIDATES = 8;
constexpr int H45_LOOKAHEAD_NEIGHBOR_LIMIT = 16;
constexpr int H45_MAX_LOOKAHEAD_CHOICES = 4;
constexpr int H45_CONFLICT_COMPONENT_WINDOW = 32;
constexpr int H45_MAX_CONFLICT_PAIRS_PER_TREE = 512;
constexpr int H45_MAX_CONFLICT_POOL_SIZE = 256;

struct H45SearchState {
    int jitter = 0;
    int merge_score_mode = 2;
    int restart_seed_offset = 0;
    int extra_split_count = 0;
    std::mt19937 rng{1};
};

struct H45Candidate {
    int score = 0;
    Leaf a = -1;
    Leaf b = -1;
};

struct H45CandidateGreater {
    bool operator()(const H45Candidate& x, const H45Candidate& y) const {
        return x.score > y.score;
    }
};

class H45MergeQueue {
public:
    bool empty() const { return heap_.empty(); }
    const H45Candidate& top() const { return heap_.front(); }

    void push(const H45Candidate& value) {
        heap_.push_back(value);
        size_t child = heap_.size() - 1;
        while (child > 0) {
            size_t parent = (child - 1) / 2;
            if (!comp_(heap_[parent], heap_[child])) break;
            std::swap(heap_[parent], heap_[child]);
            child = parent;
        }
    }

    void pop() {
        if (heap_.empty()) return;
        H45Candidate value = heap_.back();
        heap_.pop_back();
        if (heap_.empty()) return;

        size_t hole = 0;
        while (true) {
            size_t left = 2 * hole + 1;
            if (left >= heap_.size()) break;
            size_t right = left + 1;
            size_t child = right < heap_.size() && comp_(heap_[left], heap_[right])
                               ? right
                               : left;
            heap_[hole] = heap_[child];
            hole = child;
        }
        while (hole > 0) {
            size_t parent = (hole - 1) / 2;
            if (!comp_(heap_[parent], value)) break;
            heap_[hole] = heap_[parent];
            hole = parent;
        }
        heap_[hole] = value;
    }

private:
    std::vector<H45Candidate> heap_;
    H45CandidateGreater comp_;
};

uint64_t repr_pair_key(Leaf a, Leaf b) {
    Leaf lo = std::min(a, b);
    Leaf hi = std::max(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(lo)) << 32) |
           static_cast<uint32_t>(hi);
}

uint64_t comp_pair_key(CompId a, CompId b) {
    CompId lo = std::min(a, b);
    CompId hi = std::max(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(lo)) << 32) |
           static_cast<uint32_t>(hi);
}

int leaf_dist(const Tree& tree, Leaf a, Leaf b) {
    if (a <= 0 || b <= 0 || a >= static_cast<int>(tree.leaf_node.size()) ||
        b >= static_cast<int>(tree.leaf_node.size())) {
        return INT_MAX / 4;
    }
    NodeId va = tree.leaf_node[a];
    NodeId vb = tree.leaf_node[b];
    NodeId w = tree.lca(va, vb);
    if (w < 0) return INT_MAX / 4;
    return tree.depth[va] + tree.depth[vb] - 2 * tree.depth[w];
}

int h45_base_merge_score(Leaf a, Leaf b, const Tree& t1, const Tree& t2) {
    return std::max(leaf_dist(t1, a, b), leaf_dist(t2, a, b));
}

bool roots_are_siblings(const Tree& tree, NodeId a, NodeId b) {
    return a >= 0 && b >= 0 && a != b && tree.parent[a] >= 0 &&
           tree.parent[a] == tree.parent[b];
}

bool roots_are_cut_one_b_like(const Tree& tree, NodeId a, NodeId b) {
    if (a < 0 || b < 0 || a == b) return false;
    NodeId pa = tree.parent[a];
    NodeId pb = tree.parent[b];
    if (pa < 0 || pb < 0) return false;
    NodeId ppa = tree.parent[pa];
    NodeId ppb = tree.parent[pb];
    return (ppa >= 0 && ppa == pb) || (ppb >= 0 && ppb == pa);
}

int h45_component_merge_score(const Forest& forest, CompId a, CompId b,
                              const Tree& t1, const Tree& t2,
                              H45SearchState& search) {
    const Component& ca = forest.comps[a];
    const Component& cb = forest.comps[b];
    int score = 0;
    if (ca.root1 >= 0 && cb.root1 >= 0 && ca.root2 >= 0 && cb.root2 >= 0) {
        NodeId l1 = t1.lca(ca.root1, cb.root1);
        NodeId l2 = t2.lca(ca.root2, cb.root2);
        int d1 = t1.depth[ca.root1] + t1.depth[cb.root1] - 2 * t1.depth[l1];
        int d2 = t2.depth[ca.root2] + t2.depth[cb.root2] - 2 * t2.depth[l2];
        int max_distance = std::max(d1, d2);
        int total_distance = d1 + d2;
        int combined_size =
            std::min(255, static_cast<int>(ca.leaves.size() + cb.leaves.size()));
        if (search.merge_score_mode == 5) {
            score = h45_base_merge_score(ca.repr, cb.repr, t1, t2);
        } else if (search.merge_score_mode == 1) {
            score = total_distance;
        } else if (search.merge_score_mode == 2) {
            score = 32 * max_distance + total_distance;
        } else if (search.merge_score_mode == 3) {
            score = 128 * max_distance + total_distance - combined_size;
        } else if (search.merge_score_mode == 4) {
            score = 128 * total_distance - combined_size;
        } else if (search.merge_score_mode == 6) {
            int common_sibling = (roots_are_siblings(t1, ca.root1, cb.root1) ? 1 : 0) +
                                 (roots_are_siblings(t2, ca.root2, cb.root2) ? 1 : 0);
            int imbalance = std::abs(d1 - d2);
            int local_depth = std::min(t1.depth[l1], t2.depth[l2]);
            score = 512 * max_distance + 32 * imbalance + total_distance -
                    12 * local_depth - 3 * combined_size;
            if (common_sibling == 2) score -= 1000000;
            else if (common_sibling == 1) score -= 128;
        } else if (search.merge_score_mode == 14) {
            int sibling1 = roots_are_siblings(t1, ca.root1, cb.root1) ? 1 : 0;
            int sibling2 = roots_are_siblings(t2, ca.root2, cb.root2) ? 1 : 0;
            int cut_one_like =
                (roots_are_cut_one_b_like(t1, ca.root1, cb.root1) ? 1 : 0) +
                (roots_are_cut_one_b_like(t2, ca.root2, cb.root2) ? 1 : 0);
            int imbalance = std::abs(d1 - d2);
            int local_depth = std::min(t1.depth[l1], t2.depth[l2]);
            score = 384 * max_distance + 24 * total_distance + 48 * imbalance -
                    10 * local_depth - 4 * combined_size;
            if (sibling1 && sibling2) score -= 1000000;
            else if (sibling1 || sibling2) score += 2048;
            if (cut_one_like == 1) score += 512;
            else if (cut_one_like == 2) score -= 64;
        } else {
            score = max_distance;
        }
    } else {
        score = h45_base_merge_score(ca.repr, cb.repr, t1, t2);
    }
    if (search.jitter > 0) {
        score += static_cast<int>(search.rng() % static_cast<unsigned>(search.jitter));
    }
    return score;
}

std::vector<CompId> h45_local_components(const Tree& tree, const Forest& forest,
                                         NodeId root, int limit) {
    std::vector<CompId> result;
    if (root < 0 || root >= tree.node_count() || limit <= 0) return result;
    std::vector<unsigned char> seen(forest.comps.size(), 0);
    std::vector<NodeId> stack{root};
    while (!stack.empty() && static_cast<int>(result.size()) < limit) {
        NodeId v = stack.back();
        stack.pop_back();
        if (tree.is_leaf(v)) {
            CompId id = forest.leaf_comp[tree.label[v]];
            if (id >= 0 && id < static_cast<CompId>(forest.comps.size()) &&
                forest.comps[id].active && !seen[id]) {
                seen[id] = 1;
                result.push_back(id);
            }
            continue;
        }
        stack.push_back(tree.right[v]);
        stack.push_back(tree.left[v]);
    }
    return result;
}

void h45_push_candidate(H45MergeQueue& pq, std::unordered_set<uint64_t>& seen,
                        const Forest& forest, const Tree& t1, const Tree& t2,
                        H45SearchState& search, CompId a, CompId b) {
    if (a < 0 || b < 0 || a == b) return;
    if (!forest.comps[a].active || !forest.comps[b].active) return;
    Leaf ra = forest.comps[a].repr;
    Leaf rb = forest.comps[b].repr;
    uint64_t key = repr_pair_key(ra, rb);
    if (!seen.insert(key).second) return;
    pq.push({h45_component_merge_score(forest, a, b, t1, t2, search), ra, rb});
}

void h45_seed_candidates(H45MergeQueue& pq, std::unordered_set<uint64_t>& seen,
                         const Forest& forest, const Tree& t1, const Tree& t2,
                         H45SearchState& search) {
    for (const Tree* tree : {&t1, &t2}) {
        for (NodeId v = 0; v < tree->node_count(); ++v) {
            if (tree->is_leaf(v)) continue;
            std::vector<CompId> ids =
                h45_local_components(*tree, forest, v,
                                     H45_LOCAL_COMPONENT_WINDOW);
            for (size_t i = 0; i < ids.size(); ++i) {
                for (size_t j = i + 1; j < ids.size(); ++j) {
                    h45_push_candidate(pq, seen, forest, t1, t2, search,
                                       ids[i], ids[j]);
                }
            }
        }
    }
    if (t1.n < 2) return;
    unsigned range = static_cast<unsigned>(t1.n);
    unsigned mask = 1;
    while (mask < range - 1) mask = (mask << 1) | 1u;
    auto random_leaf = [&]() {
        unsigned value = 0;
        do {
            value = search.rng() & mask;
        } while (value >= range);
        return static_cast<Leaf>(value + 1);
    };
    int trials = t1.n * H45_RANDOM_CANDIDATE_MULTIPLIER;
    for (int t = 0; t < trials; ++t) {
        Leaf la = random_leaf();
        Leaf lb = random_leaf();
        if (la == lb) continue;
        CompId a = forest.leaf_comp[la];
        CompId b = forest.leaf_comp[lb];
        h45_push_candidate(pq, seen, forest, t1, t2, search, a, b);
    }
}

bool h45_active_candidate_components(const H45Candidate& candidate,
                                     const Forest& forest,
                                     CompId& a, CompId& b) {
    if (candidate.a <= 0 || candidate.b <= 0 ||
        candidate.a >= static_cast<int>(forest.leaf_comp.size()) ||
        candidate.b >= static_cast<int>(forest.leaf_comp.size())) {
        return false;
    }
    a = forest.leaf_comp[candidate.a];
    b = forest.leaf_comp[candidate.b];
    return a >= 0 && b >= 0 && a != b && forest.comps[a].active &&
           forest.comps[b].active;
}

bool h45_pop_active_candidate(H45MergeQueue& pq, const Forest& forest,
                              H45Candidate& out) {
    while (!pq.empty()) {
        out = pq.top();
        pq.pop();
        CompId a = -1;
        CompId b = -1;
        if (h45_active_candidate_components(out, forest, a, b)) return true;
    }
    return false;
}

struct H45SmallSet {
    CompId ids[3] = {-1, -1, -1};
    int count = 0;
    bool overflow = false;
};

void h45_small_add(H45SmallSet& set, CompId id) {
    if (id < 0 || set.overflow) return;
    for (int i = 0; i < set.count; ++i) {
        if (set.ids[i] == id) return;
    }
    if (set.count == 3) {
        set.overflow = true;
        return;
    }
    set.ids[set.count++] = id;
}

void h45_small_merge(H45SmallSet& a, const H45SmallSet& b) {
    if (a.overflow) return;
    if (b.overflow) {
        a.overflow = true;
        return;
    }
    for (int i = 0; i < b.count; ++i) h45_small_add(a, b.ids[i]);
}

H45SmallSet h45_collect_cherries(const Tree& tree, const Forest& forest,
                                 NodeId v,
                                 std::vector<std::pair<CompId, CompId>>& pairs,
                                 std::unordered_set<uint64_t>& seen) {
    H45SmallSet here;
    if (tree.is_leaf(v)) {
        CompId id = forest.leaf_comp[tree.label[v]];
        if (id >= 0 && forest.comps[id].active) h45_small_add(here, id);
        return here;
    }
    H45SmallSet left = h45_collect_cherries(tree, forest, tree.left[v], pairs, seen);
    H45SmallSet right = h45_collect_cherries(tree, forest, tree.right[v], pairs, seen);
    h45_small_merge(here, left);
    h45_small_merge(here, right);
    if (!here.overflow && here.count == 2) {
        CompId a = here.ids[0];
        CompId b = here.ids[1];
        uint64_t key = comp_pair_key(a, b);
        if (seen.insert(key).second) pairs.push_back({a, b});
    }
    return here;
}

std::vector<std::pair<CompId, CompId>> h45_active_cherry_pairs(const Tree& tree,
                                                               const Forest& forest) {
    std::vector<std::pair<CompId, CompId>> pairs;
    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<size_t>(tree.n) * 2 + 1);
    if (tree.root >= 0) h45_collect_cherries(tree, forest, tree.root, pairs, seen);
    return pairs;
}

bool h45_contract_common_cherries(Forest& forest, const Tree& t1, const Tree& t2) {
    bool any = false;
    while (true) {
        std::vector<std::pair<CompId, CompId>> p1 = h45_active_cherry_pairs(t1, forest);
        std::vector<std::pair<CompId, CompId>> p2 = h45_active_cherry_pairs(t2, forest);
        std::unordered_set<uint64_t> in_t2;
        in_t2.reserve(p2.size() * 2 + 1);
        for (const auto& p : p2) in_t2.insert(comp_pair_key(p.first, p.second));

        bool changed = false;
        for (const auto& p : p1) {
            if (!in_t2.count(comp_pair_key(p.first, p.second))) continue;
            if (h45_try_merge(forest, t1, t2, p.first, p.second) >= 0) {
                changed = true;
                any = true;
            }
        }
        if (!changed) break;
    }
    return any;
}

void h45_push_neighbors(H45MergeQueue& pq, std::unordered_set<uint64_t>& seen,
                        std::unordered_set<uint64_t>& failed,
                        const Forest& forest, const Tree& t1, const Tree& t2,
                        H45SearchState& search, CompId focus) {
    if (focus < 0 || !forest.comps[focus].active) return;
    std::unordered_set<CompId> refreshed;
    refreshed.reserve(H45_LOCAL_COMPONENT_WINDOW * 4 + 1);
    for (const Tree* tree : {&t1, &t2}) {
        std::unordered_set<NodeId> touched;
        for (Leaf leaf : forest.comps[focus].leaves) {
            NodeId v = tree->leaf_node[leaf];
            for (int walk = 0; v >= 0 && walk < H45_NEIGHBOR_WALK_LIMIT; ++walk) {
                if (touched.insert(v).second) {
                    std::vector<CompId> ids =
                        h45_local_components(*tree, forest, v,
                                             H45_LOCAL_COMPONENT_WINDOW);
                    for (CompId other : ids) {
                        if (other == focus || !forest.comps[other].active) continue;
                        if (!refreshed.insert(other).second) continue;
                        uint64_t key = repr_pair_key(forest.comps[focus].repr,
                                                     forest.comps[other].repr);
                        if (failed.erase(key)) seen.erase(key);
                        h45_push_candidate(pq, seen, forest, t1, t2, search,
                                           focus, other);
                    }
                }
                v = tree->parent[v];
            }
        }
    }
}

void h45_strip_cached_data(Forest& forest) {
    for (Component& component : forest.comps) {
        component.signature1.clear();
        component.signature2.clear();
        component.edges1.clear();
        component.edges2.clear();
    }
}

std::vector<CompId> h45_nearby_components(const Forest& forest,
                                           const Tree& t1, const Tree& t2,
                                           CompId focus) {
    std::vector<CompId> result;
    if (focus < 0 || focus >= static_cast<CompId>(forest.comps.size()) ||
        !forest.comps[focus].active) {
        return result;
    }
    std::unordered_set<CompId> seen;
    for (const Tree* tree : {&t1, &t2}) {
        std::unordered_set<NodeId> touched;
        for (Leaf leaf : forest.comps[focus].leaves) {
            NodeId v = tree->leaf_node[leaf];
            while (v >= 0) {
                if (touched.insert(v).second) {
                    std::vector<CompId> ids =
                        h45_local_components(*tree, forest, v,
                                             H45_LOCAL_COMPONENT_WINDOW);
                    for (CompId other : ids) {
                        if (other == focus || !forest.comps[other].active ||
                            !seen.insert(other).second) {
                            continue;
                        }
                        result.push_back(other);
                        if (static_cast<int>(result.size()) >=
                            H45_LOOKAHEAD_NEIGHBOR_LIMIT) {
                            return result;
                        }
                    }
                }
                v = tree->parent[v];
            }
        }
    }
    return result;
}

CompId h45_best_followup_merge(Forest& forest, const Tree& t1, const Tree& t2,
                               CompId focus) {
    std::vector<CompId> neighbors = h45_nearby_components(forest, t1, t2, focus);
    CompId best_neighbor = -1;
    int best_score = INT_MAX;
    int best_size = -1;
    for (CompId other : neighbors) {
        Forest trial = forest;
        CompId merged = h45_try_merge(trial, t1, t2, focus, other);
        if (merged < 0) continue;
        int score = h45_base_merge_score(forest.comps[focus].repr,
                                         forest.comps[other].repr, t1, t2);
        int size = static_cast<int>(trial.comps[merged].leaves.size());
        if (score < best_score ||
            (score == best_score && size > best_size) ||
            (score == best_score && size == best_size &&
             (best_neighbor < 0 ||
              forest.comps[other].repr < forest.comps[best_neighbor].repr))) {
            best_neighbor = other;
            best_score = score;
            best_size = size;
        }
    }
    if (best_neighbor < 0) return -1;
    return h45_try_merge(forest, t1, t2, focus, best_neighbor);
}

long long h45_lookahead_score(const H45Candidate& candidate,
                              const Forest& forest,
                              const Tree& t1, const Tree& t2) {
    constexpr long long invalid = (1LL << 60);
    CompId a = -1;
    CompId b = -1;
    if (!h45_active_candidate_components(candidate, forest, a, b)) return invalid;
    Forest trial = forest;
    CompId focus = h45_try_merge(trial, t1, t2, a, b);
    if (focus < 0) return invalid;
    CompId followup = h45_best_followup_merge(trial, t1, t2, focus);
    int followups = followup >= 0 ? 1 : 0;
    if (followup >= 0) focus = followup;
    int final_size = static_cast<int>(trial.comps[focus].leaves.size());
    return candidate.score - 1000000LL * followups - 10000LL * final_size;
}

bool h45_choose_lookahead_candidate(H45MergeQueue& pq, const Forest& forest,
                                    const Tree& t1, const Tree& t2,
                                    std::unordered_set<uint64_t>& failed,
                                    H45Candidate& chosen) {
    constexpr long long invalid = (1LL << 60);
    while (!pq.empty()) {
        std::vector<H45Candidate> batch;
        H45Candidate candidate;
        while (static_cast<int>(batch.size()) < H45_LOOKAHEAD_TOP_CANDIDATES &&
               h45_pop_active_candidate(pq, forest, candidate)) {
            batch.push_back(candidate);
        }
        if (batch.empty()) return false;

        int best_index = -1;
        int best_base_score = batch.front().score;
        long long best_score = invalid;
        std::vector<char> keep(batch.size(), 1);
        for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
            if (batch[i].score != best_base_score) continue;
            long long score = h45_lookahead_score(batch[i], forest, t1, t2);
            if (score >= invalid) {
                keep[i] = 0;
                failed.insert(repr_pair_key(batch[i].a, batch[i].b));
                continue;
            }
            if (best_index < 0 || score < best_score ||
                (score == best_score && batch[i].score < batch[best_index].score) ||
                (score == best_score && batch[i].score == batch[best_index].score &&
                 batch[i].a < batch[best_index].a) ||
                (score == best_score && batch[i].score == batch[best_index].score &&
                 batch[i].a == batch[best_index].a &&
                 batch[i].b < batch[best_index].b)) {
                best_index = i;
                best_score = score;
            }
        }
        if (best_index < 0) {
            for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
                if (batch[i].score != best_base_score) pq.push(batch[i]);
            }
            continue;
        }
        for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
            if (i != best_index && keep[i]) pq.push(batch[i]);
        }
        chosen = batch[best_index];
        return true;
    }
    return false;
}

std::vector<CompId> h45_conflict_split_candidates(const Forest& forest,
                                                  const Tree& t1,
                                                  const Tree& t2,
                                                  bool use_fpt_conflict) {
    std::vector<std::pair<CompId, CompId>> p1 = h45_active_cherry_pairs(t1, forest);
    std::vector<std::pair<CompId, CompId>> p2 = h45_active_cherry_pairs(t2, forest);
    std::unordered_set<uint64_t> in_t1;
    std::unordered_set<uint64_t> in_t2;
    for (const auto& p : p1) in_t1.insert(comp_pair_key(p.first, p.second));
    for (const auto& p : p2) in_t2.insert(comp_pair_key(p.first, p.second));

    std::vector<int> score(forest.comps.size(), 0);
    auto add_pair = [&](const Tree& other, bool other_is_t1, CompId a, CompId b) {
        if (a < 0 || b < 0 || !forest.comps[a].active || !forest.comps[b].active) {
            return;
        }
        NodeId ra = other_is_t1 ? forest.comps[a].root1 : forest.comps[a].root2;
        NodeId rb = other_is_t1 ? forest.comps[b].root1 : forest.comps[b].root2;
        NodeId l = other.lca(ra, rb);
        bool scored = false;
        if (l >= 0) {
            std::vector<CompId> ids =
                h45_local_components(other, forest, l, H45_CONFLICT_COMPONENT_WINDOW);
            for (CompId id : ids) {
                if (id < 0 || id >= static_cast<CompId>(forest.comps.size()) ||
                    !forest.comps[id].active || forest.comps[id].leaves.size() <= 1) {
                    continue;
                }
                int bonus = (id == a || id == b) ? 4 : 16;
                score[id] += bonus +
                             std::min(64, static_cast<int>(forest.comps[id].leaves.size()));
                scored = true;
            }
        }
        if (!scored) {
            for (CompId id : {a, b}) {
                if (forest.comps[id].leaves.size() > 1) {
                    score[id] += 4 +
                                 std::min(64, static_cast<int>(forest.comps[id].leaves.size()));
                }
            }
        }
    };
    auto add_fpt_score = [&](CompId id, int bonus, std::vector<int>& scores) {
        if (id < 0 || id >= static_cast<CompId>(forest.comps.size()) ||
            !forest.comps[id].active || forest.comps[id].leaves.size() <= 1) {
            return;
        }
        scores[id] += bonus +
                      std::min(96, static_cast<int>(forest.comps[id].leaves.size()));
    };
    auto score_path_side = [&](const Tree& other, NodeId child, NodeId ancestor,
                               CompId a, CompId b, int& side_hits) {
        NodeId v = child;
        int guard = 0;
        int max_depth = 0;
        for (int depth : other.depth) max_depth = std::max(max_depth, depth);
        while (v >= 0 && v != ancestor && guard++ <= max_depth + 2) {
            NodeId p = other.parent[v];
            if (p < 0) break;
            for (NodeId side : {other.left[p], other.right[p]}) {
                if (side == v || side < 0) continue;
                std::vector<CompId> ids =
                    h45_local_components(other, forest, side,
                                         H45_CONFLICT_COMPONENT_WINDOW);
                for (CompId id : ids) {
                    if (id == a || id == b) continue;
                    add_fpt_score(id, 96, score);
                    ++side_hits;
                }
            }
            v = p;
        }
    };
    auto add_fpt_pair = [&](const Tree& other, bool other_is_t1, CompId a, CompId b) {
        NodeId ra = other_is_t1 ? forest.comps[a].root1 : forest.comps[a].root2;
        NodeId rb = other_is_t1 ? forest.comps[b].root1 : forest.comps[b].root2;
        NodeId l = other.lca(ra, rb);
        if (l < 0 || ra < 0 || rb < 0) return;

        int d = other.depth[ra] + other.depth[rb] - 2 * other.depth[l];
        int side_hits = 0;
        score_path_side(other, ra, l, a, b, side_hits);
        score_path_side(other, rb, l, a, b, side_hits);

        std::vector<CompId> overlap =
            h45_local_components(other, forest, l, H45_CONFLICT_COMPONENT_WINDOW);
        for (CompId id : overlap) {
            if (id == a || id == b) continue;
            add_fpt_score(id, 32 + std::min(96, 4 * d), score);
        }

        if (side_hits == 0) {
            add_fpt_score(a, 16 + 2 * d, score);
            add_fpt_score(b, 16 + 2 * d, score);
            return;
        }
        bool cut_one_like = roots_are_cut_one_b_like(other, ra, rb);
        if (cut_one_like) {
            add_fpt_score(a, 12, score);
            add_fpt_score(b, 12, score);
        } else if (d > 4) {
            add_fpt_score(a, 24, score);
            add_fpt_score(b, 24, score);
        }
    };
    int t1_only = 0;
    int t2_only = 0;
    int examined = 0;
    for (const auto& p : p1) {
        if (examined >= H45_MAX_CONFLICT_PAIRS_PER_TREE) break;
        if (!in_t2.count(comp_pair_key(p.first, p.second))) {
            add_pair(t2, false, p.first, p.second);
            if (use_fpt_conflict) add_fpt_pair(t2, false, p.first, p.second);
            ++t1_only;
            ++examined;
        }
    }
    for (const auto& p : p2) {
        if (examined >= H45_MAX_CONFLICT_PAIRS_PER_TREE * 2) break;
        if (!in_t1.count(comp_pair_key(p.first, p.second))) {
            add_pair(t1, true, p.first, p.second);
            if (use_fpt_conflict) add_fpt_pair(t1, true, p.first, p.second);
            ++t2_only;
            ++examined;
        }
    }
    std::vector<CompId> ids;
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        if (score[id] > 0 && forest.comps[id].active &&
            forest.comps[id].leaves.size() > 1) {
            ids.push_back(id);
        }
    }
    std::sort(ids.begin(), ids.end(), [&](CompId a, CompId b) {
        if (score[a] != score[b]) return score[a] > score[b];
        if (forest.comps[a].leaves.size() != forest.comps[b].leaves.size()) {
            return forest.comps[a].leaves.size() > forest.comps[b].leaves.size();
        }
        return a < b;
    });
    if (ids.size() > H45_MAX_CONFLICT_POOL_SIZE) {
        ids.resize(H45_MAX_CONFLICT_POOL_SIZE);
    }
    if (use_fpt_conflict && std::getenv("H45_FPT_CONFLICT_TRACE") != nullptr) {
        static int trace_count = 0;
        if (trace_count++ < 32) {
            std::cerr << "fpt_conflict_pool t1_pairs=" << p1.size()
                      << " t2_pairs=" << p2.size()
                      << " t1_only=" << t1_only
                      << " t2_only=" << t2_only
                      << " selected=" << ids.size();
            for (int i = 0; i < std::min<int>(8, ids.size()); ++i) {
                CompId id = ids[i];
                std::cerr << " cid" << i << '=' << id
                          << ":score=" << score[id]
                          << ":leaves=" << forest.comps[id].leaves.size();
            }
            std::cerr << '\n';
        }
    }
    return ids;
}

Forest h45_compact_forest(const Forest& forest, const Tree& t1, const Tree& t2) {
    Forest compact;
    compact.leaf_comp.assign(t1.n + 1, -1);
    compact.owner1.assign(t1.node_count(), -1);
    compact.owner2.assign(t2.node_count(), -1);
    compact.comps.reserve(forest.active_count);
    for (const Component& component : forest.comps) {
        if (!component.active) continue;
        append_exact_component(compact, component.leaves, component.repr, t1, t2);
    }
    return compact.active_count == forest.active_count ? compact : forest;
}

Forest h45_greedy_merge_impl(Forest forest, const Tree& t1, const Tree& t2,
                             H45SearchState& search, bool use_lookahead,
                             bool reconsider_failed) {
    if (!rebuild_forest_data(forest, t1, t2)) return forest;
    h45_strip_cached_data(forest);
    h45_contract_common_cherries(forest, t1, t2);

    H45MergeQueue pq;
    std::unordered_set<uint64_t> seen;
    std::unordered_set<uint64_t> failed;
    seen.reserve(static_cast<size_t>(t1.n) * 8 + 1);
    failed.reserve(static_cast<size_t>(t1.n) * 2 + 1);
    h45_seed_candidates(pq, seen, forest, t1, t2, search);

    int max_checks = std::max(H45_MIN_MERGE_CHECKS, t1.n * H45_MERGE_CHECK_MULTIPLIER);
    int checks = 0;
    int lookahead_choices = 0;
    while (!pq.empty() && checks < max_checks) {
        H45Candidate candidate;
        if (use_lookahead && lookahead_choices < H45_MAX_LOOKAHEAD_CHOICES) {
            if (!h45_choose_lookahead_candidate(pq, forest, t1, t2,
                                                failed, candidate)) {
                break;
            }
        } else if (!h45_pop_active_candidate(pq, forest, candidate)) {
            break;
        }
        CompId a = -1;
        CompId b = -1;
        if (!h45_active_candidate_components(candidate, forest, a, b)) continue;
        ++checks;
        CompId merged = h45_try_merge(forest, t1, t2, a, b);
        if (merged >= 0) {
            failed.erase(repr_pair_key(candidate.a, candidate.b));
            if (use_lookahead && lookahead_choices < H45_MAX_LOOKAHEAD_CHOICES) {
                ++lookahead_choices;
            }
            if (h45_contract_common_cherries(forest, t1, t2)) {
                pq = H45MergeQueue();
                seen.clear();
                failed.clear();
                h45_seed_candidates(pq, seen, forest, t1, t2, search);
            } else {
                h45_push_neighbors(pq, seen, failed, forest, t1, t2, search, merged);
            }
        } else if (reconsider_failed) {
            failed.insert(repr_pair_key(candidate.a, candidate.b));
        }
    }
    return forest;
}

struct H45BlockerRepairCandidate {
    CompId a = -1;
    CompId b = -1;
    std::vector<CompId> blockers;
    int score = 0;
};

std::vector<CompId> h45_merge_blockers(const Forest& forest, const Tree& t1,
                                       const Tree& t2, CompId a, CompId b) {
    if (a < 0 || b < 0 || a == b ||
        a >= static_cast<CompId>(forest.comps.size()) ||
        b >= static_cast<CompId>(forest.comps.size()) ||
        !forest.comps[a].active || !forest.comps[b].active) {
        return {};
    }
    std::vector<Leaf> leaves = merged_leaves_of(forest.comps[a], forest.comps[b]);
    InducedInfo i1 = compute_induced(t1, leaves);
    InducedInfo i2 = compute_induced(t2, leaves);
    if (i1.signature != i2.signature) return {};

    std::vector<CompId> blockers;
    auto collect = [&](const std::vector<NodeId>& edges,
                       const std::vector<CompId>& owner) {
        for (NodeId edge : edges) {
            CompId current = owner[edge];
            if (current >= 0 && current != a && current != b &&
                current < static_cast<CompId>(forest.comps.size()) &&
                forest.comps[current].active &&
                forest.comps[current].leaves.size() > 1) {
                blockers.push_back(current);
            }
        }
    };
    collect(i1.used_edges, forest.owner1);
    collect(i2.used_edges, forest.owner2);
    std::sort(blockers.begin(), blockers.end());
    blockers.erase(std::unique(blockers.begin(), blockers.end()), blockers.end());
    return blockers;
}

Forest h45_blocker_repair_impl(Forest initial, const Tree& t1, const Tree& t2,
                               int max_components, int max_blockers,
                               int candidate_limit, int restart_budget,
                               int target_count, double max_seconds) {
    if (max_components <= 0 || max_blockers <= 0 || candidate_limit <= 0) {
        return initial;
    }
    const auto start_time = std::chrono::steady_clock::now();
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = initial;
    if (best.active_count > max_components) return best;

    std::vector<CompId> active = active_components(best);
    std::vector<H45BlockerRepairCandidate> candidates;
    H45SearchState score_search;
    score_search.merge_score_mode = 2;
    bool timed_out = false;

    for (size_t i = 0; i < active.size(); ++i) {
        if (max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) break;
        }
        for (size_t j = i + 1; j < active.size(); ++j) {
            if (max_seconds > 0.0 && ((j & 31u) == 0u)) {
                double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start_time).count();
                if (elapsed >= max_seconds) {
                    timed_out = true;
                    break;
                }
            }
            CompId a = active[i];
            CompId b = active[j];
            std::vector<CompId> blockers = h45_merge_blockers(best, t1, t2, a, b);
            if (blockers.empty() ||
                static_cast<int>(blockers.size()) > max_blockers) {
                continue;
            }
            int blocker_leaves = 0;
            bool usable = true;
            for (CompId id : blockers) {
                if (id < 0 || id >= static_cast<CompId>(best.comps.size()) ||
                    !best.comps[id].active || best.comps[id].leaves.size() <= 1) {
                    usable = false;
                    break;
                }
                blocker_leaves += static_cast<int>(best.comps[id].leaves.size());
            }
            if (!usable) continue;

            int pair_score = h45_component_merge_score(best, a, b, t1, t2,
                                                       score_search);
            int merged_size =
                static_cast<int>(best.comps[a].leaves.size() +
                                 best.comps[b].leaves.size());
            H45BlockerRepairCandidate candidate;
            candidate.a = a;
            candidate.b = b;
            candidate.blockers = std::move(blockers);
            candidate.score = 1000000 * static_cast<int>(candidate.blockers.size()) +
                              1000 * blocker_leaves + pair_score -
                              8 * merged_size;
            candidates.push_back(std::move(candidate));
        }
        if (timed_out) break;
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const H45BlockerRepairCandidate& x,
                 const H45BlockerRepairCandidate& y) {
                  if (x.score != y.score) return x.score < y.score;
                  if (x.blockers.size() != y.blockers.size()) {
                      return x.blockers.size() < y.blockers.size();
                  }
                  if (x.a != y.a) return x.a < y.a;
                  return x.b < y.b;
              });
    if (static_cast<int>(candidates.size()) > candidate_limit) {
        candidates.resize(candidate_limit);
    }

    const int modes[] = {14, 1, 2};
    int variant = 0;
    for (const H45BlockerRepairCandidate& candidate : candidates) {
        if (max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) break;
        }
        Forest seed = split_components_to_singletons(best, candidate.blockers, t1, t2);
        for (int mode : modes) {
            Forest trial = run_h45_greedy_merge(seed, t1, t2, mode, true, true);
            trial = exact_output_repair(std::move(trial), t1, t2);
            if (validate_forest(trial, t1, t2).ok &&
                trial.active_count < best.active_count) {
                best = std::move(trial);
                if (target_count > 0 && best.active_count <= target_count) return best;
            }
        }
        if (restart_budget > 0) {
            Forest restarted = run_h45_split_restarts(
                seed, t1, t2, 14, restart_budget,
                9100 + 97 * variant, 1, true, target_count,
                max_seconds > 0.0
                    ? std::max(0.1, max_seconds -
                                      std::chrono::duration<double>(
                                          std::chrono::steady_clock::now() -
                                          start_time).count())
                    : 0.0);
            restarted = exact_output_repair(std::move(restarted), t1, t2);
            if (validate_forest(restarted, t1, t2).ok &&
                restarted.active_count < best.active_count) {
                best = std::move(restarted);
                if (target_count > 0 && best.active_count <= target_count) return best;
            }
        }
        ++variant;
    }
    return best;
}

std::vector<CompId> h45_window_components(const Forest& forest,
                                          const Tree& t1, const Tree& t2,
                                          CompId focus, int window_limit) {
    std::vector<CompId> window;
    if (focus < 0 || focus >= static_cast<CompId>(forest.comps.size()) ||
        !forest.comps[focus].active || window_limit <= 0) {
        return window;
    }
    std::vector<unsigned char> seen(forest.comps.size(), 0);
    auto add = [&](CompId id) {
        if (id < 0 || id >= static_cast<CompId>(forest.comps.size()) ||
            !forest.comps[id].active || forest.comps[id].leaves.size() <= 1 ||
            seen[id]) {
            return;
        }
        seen[id] = 1;
        window.push_back(id);
    };

    add(focus);
    std::vector<CompId> frontier{focus};
    for (size_t cursor = 0;
         cursor < frontier.size() && static_cast<int>(window.size()) < window_limit;
         ++cursor) {
        for (CompId id : h45_nearby_components(forest, t1, t2, frontier[cursor])) {
            if (static_cast<int>(window.size()) >= window_limit) break;
            size_t before = window.size();
            add(id);
            if (window.size() != before) frontier.push_back(id);
        }
    }
    std::sort(window.begin(), window.end(), [&](CompId a, CompId b) {
        size_t sa = forest.comps[a].leaves.size();
        size_t sb = forest.comps[b].leaves.size();
        if (sa != sb) return sa > sb;
        return a < b;
    });
    if (static_cast<int>(window.size()) > window_limit) {
        window.resize(window_limit);
    }
    return window;
}

std::vector<std::vector<CompId>> h45_tree_component_regions(const Forest& forest,
                                                            const Tree& tree,
                                                            int max_region_components) {
    std::vector<std::vector<CompId>> regions;
    if (max_region_components <= 1) return regions;
    std::function<std::vector<CompId>(NodeId)> visit = [&](NodeId node) {
        if (tree.is_leaf(node)) {
            Leaf leaf = tree.label[node];
            if (leaf > 0 && leaf < static_cast<int>(forest.leaf_comp.size())) {
                CompId id = forest.leaf_comp[leaf];
                if (id >= 0 && id < static_cast<CompId>(forest.comps.size()) &&
                    forest.comps[id].active) {
                    return std::vector<CompId>{id};
                }
            }
            return std::vector<CompId>{};
        }
        std::vector<CompId> merged;
        for (NodeId child : {tree.left[node], tree.right[node]}) {
            if (child < 0) continue;
            std::vector<CompId> part = visit(child);
            merged.insert(merged.end(), part.begin(), part.end());
            if (static_cast<int>(merged.size()) > max_region_components * 2) {
                break;
            }
        }
        std::sort(merged.begin(), merged.end());
        merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
        if (static_cast<int>(merged.size()) >= 2 &&
            static_cast<int>(merged.size()) <= max_region_components) {
            regions.push_back(merged);
        }
        return merged;
    };
    visit(tree.root);
    std::sort(regions.begin(), regions.end(),
              [](const std::vector<CompId>& a, const std::vector<CompId>& b) {
                  if (a.size() != b.size()) return a.size() < b.size();
                  return a < b;
              });
    regions.erase(std::unique(regions.begin(), regions.end()), regions.end());
    return regions;
}

Forest h45_window_repair_impl(Forest initial, const Tree& t1, const Tree& t2,
                              int max_components, int candidate_limit,
                              int window_limit, int restart_budget,
                              int target_count, double max_seconds) {
    if (max_components <= 0 || candidate_limit <= 0 || window_limit <= 0) {
        return initial;
    }
    const auto start_time = std::chrono::steady_clock::now();
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = initial;
    if (best.active_count > max_components) return best;

    std::vector<CompId> focuses = h45_conflict_split_candidates(best, t1, t2, true);
    if (focuses.empty()) {
        focuses = largest_components(best, candidate_limit);
    }
    if (static_cast<int>(focuses.size()) > candidate_limit) {
        focuses.resize(candidate_limit);
    }

    int variant = 0;
    const int modes[] = {14, 1, 2, 6};
    std::unordered_set<uint64_t> seen_windows;
    for (CompId focus : focuses) {
        if (max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) break;
        }
        std::vector<CompId> window =
            h45_window_components(best, t1, t2, focus, window_limit);
        if (window.empty()) continue;

        std::vector<H45Hash> ids;
        ids.reserve(window.size());
        for (CompId id : window) ids.push_back(static_cast<H45Hash>(id + 1));
        uint64_t key = h45_combine_hashes(ids);
        if (!seen_windows.insert(key).second) continue;

        Forest seed = split_components_to_singletons(best, window, t1, t2);
        for (int mode : modes) {
            Forest trial = run_h45_greedy_merge(seed, t1, t2, mode, true, true);
            trial = exact_output_repair(std::move(trial), t1, t2);
            if (validate_forest(trial, t1, t2).ok &&
                trial.active_count < best.active_count) {
                best = std::move(trial);
                if (target_count > 0 && best.active_count <= target_count) return best;
            }
        }
        if (restart_budget > 0) {
            double remaining = 0.0;
            if (max_seconds > 0.0) {
                double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start_time).count();
                remaining = std::max(0.1, max_seconds - elapsed);
            }
            Forest restarted = run_h45_split_restarts(
                seed, t1, t2, 14, restart_budget, 12000 + 131 * variant,
                1, true, target_count, remaining);
            restarted = exact_output_repair(std::move(restarted), t1, t2);
            if (validate_forest(restarted, t1, t2).ok &&
                restarted.active_count < best.active_count) {
                best = std::move(restarted);
                if (target_count > 0 && best.active_count <= target_count) return best;
            }
        }
        ++variant;
    }
    return best;
}

Forest h45_cross_recombine_repair_impl(Forest initial, const Tree& t1,
                                       const Tree& t2, int max_components,
                                       int candidate_limit, int window_limit,
                                       int restart_budget, int target_count,
                                       double max_seconds) {
    if (max_components <= 0 || candidate_limit <= 0 || window_limit <= 0) {
        return initial;
    }
    const auto start_time = std::chrono::steady_clock::now();
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = initial;
    if (best.active_count > max_components) return best;

    std::vector<CompId> focuses = h45_conflict_split_candidates(best, t1, t2, true);
    std::vector<CompId> large = largest_components(best, candidate_limit * 2);
    focuses.insert(focuses.end(), large.begin(), large.end());
    std::vector<unsigned char> seen(best.comps.size(), 0);
    std::vector<CompId> unique_focus;
    for (CompId id : focuses) {
        if (id < 0 || id >= static_cast<CompId>(best.comps.size()) ||
            !best.comps[id].active || best.comps[id].leaves.size() <= 1 ||
            seen[id]) {
            continue;
        }
        seen[id] = 1;
        unique_focus.push_back(id);
        if (static_cast<int>(unique_focus.size()) >= candidate_limit) break;
    }

    const int modes[] = {14, 0, 1, 2, 5, 6};
    int variant = 0;
    for (CompId focus : unique_focus) {
        if (max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) break;
        }
        std::vector<CompId> region =
            h45_window_components(best, t1, t2, focus, window_limit);
        if (region.empty()) region.push_back(focus);
        if (std::find(region.begin(), region.end(), focus) == region.end()) {
            region.push_back(focus);
        }

        std::vector<std::vector<CompId>> split_sets;
        split_sets.push_back({focus});
        split_sets.push_back(region);

        for (const std::vector<CompId>& split_set : split_sets) {
            if (max_seconds > 0.0) {
                double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start_time).count();
                if (elapsed >= max_seconds) break;
            }
            Forest seed =
                split_components_to_singletons(best, split_set, t1, t2);
            for (int mode : modes) {
                Forest trial = run_h45_greedy_merge(seed, t1, t2, mode, true, true);
                trial = exact_output_repair(std::move(trial), t1, t2);
                if (validate_forest(trial, t1, t2).ok &&
                    trial.active_count < best.active_count) {
                    best = std::move(trial);
                    if (target_count > 0 && best.active_count <= target_count) {
                        return best;
                    }
                }
            }
            if (restart_budget > 0) {
                double remaining = 0.0;
                if (max_seconds > 0.0) {
                    double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start_time).count();
                    remaining = std::max(0.1, max_seconds - elapsed);
                }
                Forest restarted = run_h45_split_restarts(
                    seed, t1, t2, 14, restart_budget,
                    17000 + 211 * variant,
                    std::max(1, std::min(4, static_cast<int>(split_set.size()))),
                    true, target_count, remaining);
                restarted = exact_output_repair(std::move(restarted), t1, t2);
                if (validate_forest(restarted, t1, t2).ok &&
                    restarted.active_count < best.active_count) {
                    best = std::move(restarted);
                    if (target_count > 0 && best.active_count <= target_count) {
                        return best;
                    }
                }
            }
            ++variant;
        }
    }
    return best;
}

struct H45LocalColumn {
    std::vector<Leaf> leaves;
    uint64_t mask = 0;
    std::vector<NodeId> edges1;
    std::vector<NodeId> edges2;
};

bool h45_edges_available_against_outside(const std::vector<NodeId>& edges,
                                         const std::vector<CompId>& owner,
                                         const std::vector<unsigned char>& in_region) {
    for (NodeId edge : edges) {
        CompId id = owner[edge];
        if (id >= 0 &&
            (id >= static_cast<CompId>(in_region.size()) || !in_region[id])) {
            return false;
        }
    }
    return true;
}

uint64_t h45_leaf_mask(const std::vector<Leaf>& leaves,
                       const std::vector<int>& leaf_pos) {
    uint64_t mask = 0;
    for (Leaf leaf : leaves) {
        if (leaf <= 0 || leaf >= static_cast<int>(leaf_pos.size())) return 0;
        int pos = leaf_pos[leaf];
        if (pos < 0 || pos >= 64) return 0;
        mask |= (uint64_t{1} << pos);
    }
    return mask;
}

void h45_add_local_column(std::vector<H45LocalColumn>& columns,
                          std::unordered_set<uint64_t>& seen_masks,
                          const std::vector<Leaf>& raw_leaves,
                          const std::vector<int>& leaf_pos,
                          const std::vector<unsigned char>& in_region,
                          const Forest& forest, const Tree& t1, const Tree& t2) {
    if (raw_leaves.empty()) return;
    std::vector<Leaf> leaves = raw_leaves;
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    uint64_t mask = h45_leaf_mask(leaves, leaf_pos);
    if (mask == 0 || seen_masks.count(mask)) return;

    InducedInfo i1 = compute_induced(t1, leaves);
    InducedInfo i2 = compute_induced(t2, leaves);
    if (i1.signature != i2.signature) return;
    if (!h45_edges_available_against_outside(i1.used_edges, forest.owner1, in_region)) {
        return;
    }
    if (!h45_edges_available_against_outside(i2.used_edges, forest.owner2, in_region)) {
        return;
    }

    H45LocalColumn column;
    column.leaves = std::move(leaves);
    column.mask = mask;
    column.edges1 = std::move(i1.used_edges);
    column.edges2 = std::move(i2.used_edges);
    seen_masks.insert(mask);
    columns.push_back(std::move(column));
}

std::vector<uint64_t> h45_region_clade_masks(const Tree& tree,
                                             const std::vector<Leaf>& region_leaves,
                                             const std::vector<int>& leaf_pos) {
    std::vector<uint64_t> masks;
    masks.reserve(tree.node_count());
    for (NodeId node = 0; node < tree.node_count(); ++node) {
        uint64_t mask = 0;
        for (Leaf leaf : region_leaves) {
            NodeId leaf_node = tree.leaf_node[leaf];
            if (tree.is_ancestor(node, leaf_node)) {
                mask |= (uint64_t{1} << leaf_pos[leaf]);
            }
        }
        if ((mask & (mask - 1)) != 0) masks.push_back(mask);
    }
    std::sort(masks.begin(), masks.end());
    masks.erase(std::unique(masks.begin(), masks.end()), masks.end());
    return masks;
}

std::vector<Leaf> h45_leaves_from_mask(uint64_t mask,
                                       const std::vector<Leaf>& region_leaves) {
    std::vector<Leaf> leaves;
    while (mask != 0) {
        int bit = __builtin_ctzll(mask);
        leaves.push_back(region_leaves[bit]);
        mask &= mask - 1;
    }
    return leaves;
}

bool h45_column_edges_conflict(const H45LocalColumn& column,
                               const std::vector<unsigned char>& used1,
                               const std::vector<unsigned char>& used2) {
    for (NodeId edge : column.edges1) {
        if (edge >= 0 && edge < static_cast<NodeId>(used1.size()) && used1[edge]) {
            return true;
        }
    }
    for (NodeId edge : column.edges2) {
        if (edge >= 0 && edge < static_cast<NodeId>(used2.size()) && used2[edge]) {
            return true;
        }
    }
    return false;
}

void h45_mark_column_edges(const H45LocalColumn& column,
                           std::vector<unsigned char>& used1,
                           std::vector<unsigned char>& used2,
                           unsigned char value) {
    for (NodeId edge : column.edges1) {
        if (edge >= 0 && edge < static_cast<NodeId>(used1.size())) used1[edge] = value;
    }
    for (NodeId edge : column.edges2) {
        if (edge >= 0 && edge < static_cast<NodeId>(used2.size())) used2[edge] = value;
    }
}

bool h45_solve_local_exact_cover(
    const std::vector<H45LocalColumn>& columns,
    const std::vector<std::vector<int>>& by_leaf,
    int leaf_count,
    int incumbent_count,
    const Tree& t1,
    const Tree& t2,
    std::vector<int>& best_solution,
    double max_seconds,
    const std::chrono::steady_clock::time_point& start_time) {
    const uint64_t full_mask =
        leaf_count == 64 ? ~uint64_t{0} : ((uint64_t{1} << leaf_count) - 1);
    std::vector<int> current;
    std::vector<unsigned char> used1(t1.node_count(), 0);
    std::vector<unsigned char> used2(t2.node_count(), 0);
    int best_count = incumbent_count;
    bool improved = false;
    int nodes = 0;
    const int node_limit =
        std::max(10000, h45_env_int("RS_LOCAL_EXCHANGE_NODE_LIMIT", 250000));

    std::function<void(uint64_t)> dfs = [&](uint64_t covered) {
        if (++nodes > node_limit) return;
        if (max_seconds > 0.0 && (nodes & 1023) == 0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) return;
        }
        if (static_cast<int>(current.size()) >= best_count) return;
        uint64_t uncovered = full_mask & ~covered;
        if (uncovered == 0) {
            best_count = static_cast<int>(current.size());
            best_solution = current;
            improved = true;
            return;
        }

        int best_leaf = -1;
        int best_options = INT_MAX;
        uint64_t scan = uncovered;
        while (scan != 0) {
            int leaf = __builtin_ctzll(scan);
            int options = 0;
            for (int col : by_leaf[leaf]) {
                const H45LocalColumn& column = columns[col];
                if ((column.mask & covered) == 0 &&
                    !h45_column_edges_conflict(column, used1, used2)) {
                    ++options;
                }
            }
            if (options < best_options) {
                best_options = options;
                best_leaf = leaf;
                if (options <= 1) break;
            }
            scan &= scan - 1;
        }
        if (best_leaf < 0 || best_options == 0) return;

        for (int col : by_leaf[best_leaf]) {
            const H45LocalColumn& column = columns[col];
            if ((column.mask & covered) != 0) continue;
            if (h45_column_edges_conflict(column, used1, used2)) continue;
            current.push_back(col);
            h45_mark_column_edges(column, used1, used2, 1);
            dfs(covered | column.mask);
            h45_mark_column_edges(column, used1, used2, 0);
            current.pop_back();
        }
    };

    dfs(0);
    return improved;
}

Forest h45_apply_local_columns(const Forest& forest,
                               const std::vector<CompId>& region,
                               const std::vector<H45LocalColumn>& columns,
                               const std::vector<int>& selected,
                               const Tree& t1, const Tree& t2) {
    std::vector<unsigned char> in_region(forest.comps.size(), 0);
    for (CompId id : region) {
        if (id >= 0 && id < static_cast<CompId>(in_region.size())) in_region[id] = 1;
    }

    Forest out;
    out.leaf_comp.assign(t1.n + 1, -1);
    out.owner1.assign(t1.node_count(), -1);
    out.owner2.assign(t2.node_count(), -1);
    out.comps.reserve(forest.active_count - static_cast<int>(region.size()) +
                      static_cast<int>(selected.size()));
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        if (!forest.comps[id].active || in_region[id]) continue;
        if (append_exact_component(out, forest.comps[id].leaves, forest.comps[id].repr,
                                   t1, t2) < 0) {
            return forest;
        }
    }
    for (int col : selected) {
        if (append_exact_component(out, columns[col].leaves, columns[col].leaves.front(),
                                   t1, t2) < 0) {
            return forest;
        }
    }
    return out;
}

Forest h45_local_exchange_repair_impl(Forest initial, const Tree& t1, const Tree& t2,
                                      int max_components, int candidate_limit,
                                      int window_limit, int, int, int target_count,
                                      double max_seconds) {
    if (max_components <= 0 || candidate_limit <= 0 || window_limit <= 0) {
        return initial;
    }
    const auto start_time = std::chrono::steady_clock::now();
    auto timed_out = [&]() {
        if (max_seconds <= 0.0) return false;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        return elapsed >= max_seconds;
    };
    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = std::move(initial);
    if (best.active_count > max_components) return best;
    const bool debug = std::getenv("RS_LOCAL_EXCHANGE_DEBUG") != nullptr;
    const int max_region_leaves = h45_env_int("RS_LOCAL_EXCHANGE_MAX_LEAVES", 46);
    const int max_columns = h45_env_int("RS_LOCAL_EXCHANGE_MAX_COLUMNS", 6000);
    const int max_component_union = h45_env_int("RS_LOCAL_EXCHANGE_COMPONENT_UNION", 3);
    const int max_tree_regions = h45_env_int("RS_LOCAL_EXCHANGE_TREE_REGIONS", 0);

    std::vector<CompId> focuses = h45_conflict_split_candidates(best, t1, t2, true);
    std::vector<CompId> large = largest_components(best, candidate_limit);
    focuses.insert(focuses.end(), large.begin(), large.end());
    std::vector<unsigned char> seen_focus(best.comps.size(), 0);
    std::vector<CompId> unique_focus;
    for (CompId id : focuses) {
        if (id < 0 || id >= static_cast<CompId>(best.comps.size()) ||
            !best.comps[id].active || seen_focus[id]) {
            continue;
        }
        seen_focus[id] = 1;
        unique_focus.push_back(id);
        if (static_cast<int>(unique_focus.size()) >= candidate_limit) break;
    }

    std::vector<std::vector<CompId>> candidate_regions;
    candidate_regions.reserve(unique_focus.size() + std::max(0, max_tree_regions) * 2);
    for (CompId focus : unique_focus) {
        std::vector<CompId> region =
            h45_window_components(best, t1, t2, focus, window_limit);
        if (region.empty()) region.push_back(focus);
        if (std::find(region.begin(), region.end(), focus) == region.end()) {
            region.push_back(focus);
        }
        candidate_regions.push_back(std::move(region));
    }
    if (max_tree_regions > 0) {
        std::vector<std::vector<CompId>> tree_regions =
            h45_tree_component_regions(best, t1, window_limit);
        std::vector<std::vector<CompId>> tree_regions2 =
            h45_tree_component_regions(best, t2, window_limit);
        tree_regions.insert(tree_regions.end(), tree_regions2.begin(), tree_regions2.end());
        std::sort(tree_regions.begin(), tree_regions.end(),
                  [](const std::vector<CompId>& a, const std::vector<CompId>& b) {
                      if (a.size() != b.size()) return a.size() < b.size();
                      return a < b;
                  });
        tree_regions.erase(std::unique(tree_regions.begin(), tree_regions.end()),
                           tree_regions.end());
        if (static_cast<int>(tree_regions.size()) > max_tree_regions) {
            tree_regions.resize(max_tree_regions);
        }
        candidate_regions.insert(candidate_regions.end(), tree_regions.begin(),
                                 tree_regions.end());
    }

    std::unordered_set<uint64_t> seen_regions;
    int region_index = 0;
    for (std::vector<CompId> region : candidate_regions) {
        if (timed_out()) break;
        std::sort(region.begin(), region.end());
        region.erase(std::unique(region.begin(), region.end()), region.end());

        std::vector<Leaf> region_leaves;
        for (CompId id : region) {
            region_leaves.insert(region_leaves.end(), best.comps[id].leaves.begin(),
                                 best.comps[id].leaves.end());
        }
        std::sort(region_leaves.begin(), region_leaves.end());
        region_leaves.erase(std::unique(region_leaves.begin(), region_leaves.end()),
                            region_leaves.end());
        if (region_leaves.empty() ||
            static_cast<int>(region_leaves.size()) > std::min(63, max_region_leaves)) {
            continue;
        }

        std::vector<H45Hash> region_ids;
        for (CompId id : region) region_ids.push_back(static_cast<H45Hash>(id + 1));
        uint64_t region_key = h45_combine_hashes(region_ids);
        if (!seen_regions.insert(region_key).second) continue;

        std::vector<int> leaf_pos(t1.n + 1, -1);
        for (int i = 0; i < static_cast<int>(region_leaves.size()); ++i) {
            leaf_pos[region_leaves[i]] = i;
        }
        std::vector<unsigned char> in_region(best.comps.size(), 0);
        for (CompId id : region) in_region[id] = 1;

        std::vector<H45LocalColumn> columns;
        std::unordered_set<uint64_t> seen_masks;
        seen_masks.reserve(8192);
        for (Leaf leaf : region_leaves) {
            h45_add_local_column(columns, seen_masks, {leaf}, leaf_pos, in_region,
                                 best, t1, t2);
        }
        for (CompId id : region) {
            h45_add_local_column(columns, seen_masks, best.comps[id].leaves,
                                 leaf_pos, in_region, best, t1, t2);
        }

        for (int i = 0; i < static_cast<int>(region_leaves.size()); ++i) {
            for (int j = i + 1; j < static_cast<int>(region_leaves.size()); ++j) {
                h45_add_local_column(columns, seen_masks,
                                     {region_leaves[i], region_leaves[j]},
                                     leaf_pos, in_region, best, t1, t2);
            }
        }
        if (region_leaves.size() <= 36) {
            for (int i = 0; i < static_cast<int>(region_leaves.size()); ++i) {
                for (int j = i + 1; j < static_cast<int>(region_leaves.size()); ++j) {
                    for (int k = j + 1; k < static_cast<int>(region_leaves.size()); ++k) {
                        h45_add_local_column(columns, seen_masks,
                                             {region_leaves[i], region_leaves[j],
                                              region_leaves[k]},
                                             leaf_pos, in_region, best, t1, t2);
                    }
                }
            }
        }

        std::vector<uint64_t> clades1 =
            h45_region_clade_masks(t1, region_leaves, leaf_pos);
        std::vector<uint64_t> clades2 =
            h45_region_clade_masks(t2, region_leaves, leaf_pos);
        for (uint64_t mask : clades1) {
            h45_add_local_column(columns, seen_masks,
                                 h45_leaves_from_mask(mask, region_leaves),
                                 leaf_pos, in_region, best, t1, t2);
        }
        for (uint64_t mask : clades2) {
            h45_add_local_column(columns, seen_masks,
                                 h45_leaves_from_mask(mask, region_leaves),
                                 leaf_pos, in_region, best, t1, t2);
        }
        for (uint64_t a : clades1) {
            for (uint64_t b : clades2) {
                uint64_t mask = a & b;
                if ((mask & (mask - 1)) != 0) {
                    h45_add_local_column(columns, seen_masks,
                                         h45_leaves_from_mask(mask, region_leaves),
                                         leaf_pos, in_region, best, t1, t2);
                }
                uint64_t left_only = a & ~b;
                if ((left_only & (left_only - 1)) != 0) {
                    h45_add_local_column(columns, seen_masks,
                                         h45_leaves_from_mask(left_only, region_leaves),
                                         leaf_pos, in_region, best, t1, t2);
                }
                uint64_t right_only = b & ~a;
                if ((right_only & (right_only - 1)) != 0) {
                    h45_add_local_column(columns, seen_masks,
                                         h45_leaves_from_mask(right_only, region_leaves),
                                         leaf_pos, in_region, best, t1, t2);
                }
                uint64_t both = a | b;
                if ((both & (both - 1)) != 0) {
                    h45_add_local_column(columns, seen_masks,
                                         h45_leaves_from_mask(both, region_leaves),
                                         leaf_pos, in_region, best, t1, t2);
                }
            }
        }

        if (max_component_union >= 2) {
            for (int i = 0; i < static_cast<int>(region.size()); ++i) {
                for (int j = i + 1; j < static_cast<int>(region.size()); ++j) {
                    std::vector<Leaf> leaves = best.comps[region[i]].leaves;
                    leaves.insert(leaves.end(), best.comps[region[j]].leaves.begin(),
                                  best.comps[region[j]].leaves.end());
                    h45_add_local_column(columns, seen_masks, leaves, leaf_pos,
                                         in_region, best, t1, t2);
                    if (max_component_union >= 3) {
                        for (int k = j + 1; k < static_cast<int>(region.size()); ++k) {
                            std::vector<Leaf> leaves3 = leaves;
                            leaves3.insert(leaves3.end(),
                                           best.comps[region[k]].leaves.begin(),
                                           best.comps[region[k]].leaves.end());
                            h45_add_local_column(columns, seen_masks, leaves3,
                                                 leaf_pos, in_region, best, t1, t2);
                        }
                    }
                }
            }
        }
        if (static_cast<int>(columns.size()) > max_columns) {
            std::sort(columns.begin(), columns.end(),
                      [](const H45LocalColumn& a, const H45LocalColumn& b) {
                          int pa = __builtin_popcountll(a.mask);
                          int pb = __builtin_popcountll(b.mask);
                          if (pa != pb) return pa > pb;
                          return a.mask < b.mask;
                      });
            columns.resize(max_columns);
        }
        std::sort(columns.begin(), columns.end(),
                  [](const H45LocalColumn& a, const H45LocalColumn& b) {
                      int pa = __builtin_popcountll(a.mask);
                      int pb = __builtin_popcountll(b.mask);
                      if (pa != pb) return pa > pb;
                      return a.mask < b.mask;
                  });

        std::vector<std::vector<int>> by_leaf(region_leaves.size());
        for (int i = 0; i < static_cast<int>(columns.size()); ++i) {
            uint64_t mask = columns[i].mask;
            while (mask != 0) {
                int bit = __builtin_ctzll(mask);
                by_leaf[bit].push_back(i);
                mask &= mask - 1;
            }
        }

        std::vector<int> selected;
        bool improved = h45_solve_local_exact_cover(
            columns, by_leaf, static_cast<int>(region_leaves.size()),
            static_cast<int>(region.size()), t1, t2, selected,
            max_seconds, start_time);
        if (debug) {
            std::cerr << "rs_local_exchange_region index=" << region_index
                      << " comps=" << region.size()
                      << " leaves=" << region_leaves.size()
                      << " columns=" << columns.size()
                      << " improved=" << (improved ? 1 : 0)
                      << " selected=" << selected.size() << '\n';
        }
        ++region_index;
        if (!improved || selected.empty() ||
            static_cast<int>(selected.size()) >= static_cast<int>(region.size())) {
            continue;
        }

        Forest trial = h45_apply_local_columns(best, region, columns, selected, t1, t2);
        trial = exact_output_repair(std::move(trial), t1, t2);
        ValidationResult check = validate_forest(trial, t1, t2);
        if (check.ok && trial.active_count < best.active_count) {
            best = std::move(trial);
            if (debug) {
                std::cerr << "rs_local_exchange_improve components="
                          << best.active_count << '\n';
            }
            if (target_count > 0 && best.active_count <= target_count) return best;
            return best;
        }
    }
    return best;
}

}  // namespace

Forest exact_output_repair(Forest forest, const Tree& t1, const Tree& t2) {
    Forest repaired;
    repaired.leaf_comp.assign(t1.n + 1, -1);
    repaired.owner1.assign(t1.node_count(), -1);
    repaired.owner2.assign(t2.node_count(), -1);
    repaired.comps.reserve(forest.comps.size());
    for (const Component& component : forest.comps) {
        if (!component.active) continue;
        if (append_exact_component(repaired, component.leaves, component.repr, t1, t2) >= 0) {
            continue;
        }
        for (Leaf leaf : component.leaves) {
            append_exact_component(repaired, std::vector<Leaf>{leaf}, leaf, t1, t2);
        }
    }
    return repaired;
}

Forest run_h45_greedy_merge(Forest forest, const Tree& t1, const Tree& t2,
                            int score_mode, bool use_lookahead,
                            bool reconsider_failed) {
    H45SearchState search;
    search.merge_score_mode = score_mode;
    return h45_greedy_merge_impl(std::move(forest), t1, t2, search,
                                 use_lookahead, reconsider_failed);
}

Forest run_h45_split_restarts(Forest initial, const Tree& t1, const Tree& t2,
                              int score_mode, int attempts,
                              int seed_offset, int extra_splits,
                              bool use_fpt_conflict,
                              int target_count,
                              double max_seconds) {
    if (attempts <= 0) return initial;
    const auto start_time = std::chrono::steady_clock::now();
    Forest best = initial;
    Forest current = std::move(initial);
    int best_count = best.active_count;
    int current_count = current.active_count;
    std::mt19937 pick_rng(987654321u + 104729u * static_cast<unsigned>(seed_offset));
    std::vector<CompId> conflict_pool;
    int conflict_pool_current_count = -1;

    for (int attempt = 0; attempt < attempts; ++attempt) {
        if (max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed >= max_seconds) break;
        }
        std::vector<CompId> splittable;
        for (CompId id = 0; id < static_cast<CompId>(current.comps.size()); ++id) {
            if (current.comps[id].active && current.comps[id].leaves.size() > 1) {
                splittable.push_back(id);
            }
        }
        if (splittable.empty()) break;

        std::vector<CompId> split_ids;
        if ((attempt % H45_CONFLICT_RESTART_INTERVAL) == 0) {
            if (conflict_pool_current_count != current_count) {
                conflict_pool =
                    h45_conflict_split_candidates(current, t1, t2, use_fpt_conflict);
                conflict_pool_current_count = current_count;
            }
            for (size_t s = 0; s < conflict_pool.size(); ++s) {
                CompId id = conflict_pool[((attempt / H45_CONFLICT_RESTART_INTERVAL) + s +
                                           static_cast<size_t>(seed_offset)) %
                                          conflict_pool.size()];
                if (id >= 0 && id < static_cast<CompId>(current.comps.size()) &&
                    current.comps[id].active && current.comps[id].leaves.size() > 1) {
                    split_ids.push_back(id);
                    break;
                }
            }
        }
        if (split_ids.empty()) {
            std::sort(splittable.begin(), splittable.end(), [&](CompId a, CompId b) {
                return current.comps[a].leaves.size() >
                       current.comps[b].leaves.size();
            });
            split_ids.push_back(splittable[(attempt + seed_offset) % splittable.size()]);
        }
        if (splittable.size() > 1 &&
            (attempt % H45_CONFLICT_RESTART_INTERVAL) != 0) {
            CompId id = splittable[pick_rng() % splittable.size()];
            if (id != split_ids[0]) split_ids.push_back(id);
        }
        for (int extra = 0; extra < extra_splits &&
                            split_ids.size() < splittable.size(); ++extra) {
            CompId id = splittable[pick_rng() % splittable.size()];
            if (std::find(split_ids.begin(), split_ids.end(), id) == split_ids.end()) {
                split_ids.push_back(id);
            }
        }

        Forest trial = split_components_to_singletons(current, split_ids, t1, t2);
        H45SearchState search;
        search.merge_score_mode = score_mode;
        search.restart_seed_offset = seed_offset;
        int mode = attempt % H45_JITTER_CYCLE_LENGTH;
        if (mode < 8) search.jitter = 3 + mode;
        else if (mode < 16) search.jitter = 12 + 3 * (mode - 8);
        else search.jitter = 40 + 8 * (mode - 16);
        search.rng.seed(7654321u + 1000003u * static_cast<unsigned>(seed_offset) +
                        65537u * static_cast<unsigned>(attempt));
        bool lookahead =
            (attempt % H45_LOOKAHEAD_RESTART_INTERVAL) ==
            H45_LOOKAHEAD_RESTART_INTERVAL - 1;
        trial = h45_greedy_merge_impl(std::move(trial), t1, t2, search,
                                      lookahead, lookahead);

        int trial_count = trial.active_count;
        if (trial_count <= current_count) {
            current = trial;
            current_count = trial_count;
            if (current.comps.size() > static_cast<size_t>(current_count) * 4) {
                current = h45_compact_forest(current, t1, t2);
            }
            conflict_pool_current_count = -1;
            if (current_count < best_count) {
                best = current;
                best_count = current_count;
                if (target_count > 0 && best_count <= target_count) break;
            }
        }
    }
    return best;
}

Forest run_h45_blocker_repair(Forest initial, const Tree& t1, const Tree& t2,
                              int max_components, int max_blockers,
                              int candidate_limit, int restart_budget,
                              int target_count, double max_seconds) {
    return h45_blocker_repair_impl(std::move(initial), t1, t2, max_components,
                                   max_blockers, candidate_limit, restart_budget,
                                   target_count, max_seconds);
}

Forest run_h45_window_repair(Forest initial, const Tree& t1, const Tree& t2,
                             int max_components, int candidate_limit,
                             int window_limit, int restart_budget,
                             int target_count, double max_seconds) {
    return h45_window_repair_impl(std::move(initial), t1, t2, max_components,
                                  candidate_limit, window_limit, restart_budget,
                                  target_count, max_seconds);
}

Forest run_h45_local_exchange_repair(Forest initial, const Tree& t1,
                                     const Tree& t2, int max_components,
                                     int candidate_limit, int window_limit,
                                     int max_blockers, int restart_budget,
                                     int target_count, double max_seconds) {
    return h45_local_exchange_repair_impl(
        std::move(initial), t1, t2, max_components, candidate_limit,
        window_limit, max_blockers, restart_budget, target_count, max_seconds);
}

Forest run_h45_cross_recombine_repair(Forest initial, const Tree& t1,
                                      const Tree& t2, int max_components,
                                      int candidate_limit, int window_limit,
                                      int restart_budget, int target_count,
                                      double max_seconds) {
    return h45_cross_recombine_repair_impl(
        std::move(initial), t1, t2, max_components, candidate_limit,
        window_limit, restart_budget, target_count, max_seconds);
}

Forest run_h45_targeted_split_repack(Forest initial, const Tree& t1,
                                     const Tree& t2, int max_components,
                                     int candidate_limit, int group_width,
                                     int restart_budget, int passes,
                                     int target_count, double max_seconds) {
    if (max_components <= 0 || candidate_limit <= 0 || group_width <= 0 ||
        passes <= 0) {
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

    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = std::move(initial);
    if (best.active_count > max_components) return best;
    const bool debug = std::getenv("RS_SPLIT_REPACK_DEBUG") != nullptr;
    const int modes[] = {14, 0, 2, 1, 5};
    int variant = 0;

    for (int pass = 0; pass < passes && time_left(); ++pass) {
        std::vector<CompId> candidates;
        for (CompId id : active_components(best)) {
            if (best.comps[id].leaves.size() > 1) candidates.push_back(id);
        }
        std::sort(candidates.begin(), candidates.end(), [&](CompId a, CompId b) {
            size_t sa = best.comps[a].leaves.size();
            size_t sb = best.comps[b].leaves.size();
            if (sa != sb) return sa > sb;
            return a < b;
        });
        if (static_cast<int>(candidates.size()) > candidate_limit) {
            candidates.resize(candidate_limit);
        }
        if (debug) {
            std::cerr << "rs_split_repack_start pass=" << pass
                      << " components=" << best.active_count
                      << " candidates=" << candidates.size() << '\n';
        }

        bool improved = false;
        for (size_t i = 0; i < candidates.size() && time_left(); ++i) {
            std::vector<CompId> split_ids;
            split_ids.push_back(candidates[i]);
            for (int extra = 1; extra < group_width &&
                                i + static_cast<size_t>(extra) < candidates.size();
                 ++extra) {
                split_ids.push_back(candidates[i + static_cast<size_t>(extra)]);
            }

            Forest seed = split_components_to_singletons(best, split_ids, t1, t2);
            if (debug) {
                int leaf_count = 0;
                for (CompId id : split_ids) {
                    leaf_count += static_cast<int>(best.comps[id].leaves.size());
                }
                std::cerr << "rs_split_repack_candidate variant=" << variant
                          << " split_components=" << split_ids.size()
                          << " leaves=" << leaf_count
                          << " seed_components=" << seed.active_count << '\n';
            }

            for (int mode : modes) {
                if (!time_left()) break;
                Forest trial = run_h45_greedy_merge(seed, t1, t2, mode, true, true);
                trial = exact_output_repair(std::move(trial), t1, t2);
                ValidationResult check = validate_forest(trial, t1, t2);
                if (debug) {
                    std::cerr << "rs_split_repack_trial variant=" << variant
                              << " mode=" << mode
                              << " components=" << trial.active_count
                              << " valid=" << (check.ok ? 1 : 0) << '\n';
                }
                if (check.ok && trial.active_count < best.active_count) {
                    best = std::move(trial);
                    improved = true;
                    if (debug) {
                        std::cerr << "rs_split_repack_improve components="
                                  << best.active_count << '\n';
                    }
                    if (target_count > 0 && best.active_count <= target_count) {
                        return best;
                    }
                    break;
                }
            }

            if (!improved && restart_budget > 0 && time_left()) {
                const int restart_modes[] = {14, 0, 2};
                for (int mode : restart_modes) {
                    if (!time_left()) break;
                    Forest trial = run_h45_split_restarts(
                        seed, t1, t2, mode, restart_budget,
                        23000 + 997 * variant + 37 * mode, 1, true,
                        target_count, remaining_seconds());
                    trial = exact_output_repair(std::move(trial), t1, t2);
                    ValidationResult check = validate_forest(trial, t1, t2);
                    if (debug) {
                        std::cerr << "rs_split_repack_restart variant="
                                  << variant
                                  << " mode=" << mode
                                  << " components=" << trial.active_count
                                  << " valid=" << (check.ok ? 1 : 0) << '\n';
                    }
                    if (check.ok && trial.active_count < best.active_count) {
                        best = std::move(trial);
                        improved = true;
                        if (debug) {
                            std::cerr << "rs_split_repack_improve components="
                                      << best.active_count << '\n';
                        }
                        if (target_count > 0 &&
                            best.active_count <= target_count) {
                            return best;
                        }
                        break;
                    }
                }
            }
            ++variant;
            if (improved) break;
        }
        if (!improved) break;
    }
    return best;
}

Forest run_greedy_merge(Forest forest, const Tree& t1, const Tree& t2,
                        const GreedyConfig& config) {
    const bool debug = std::getenv("RS_DEBUG_ROUTES") != nullptr;
    struct Route {
        CandidateMode mode;
        int seed;
        int jitter;
    };
    std::vector<Route> routes = {
        {CandidateMode::SiblingAllPairs, 0, 0},
        {CandidateMode::WindowScored, 0, 0},
        {CandidateMode::WindowTight, 0, 0},
        {CandidateMode::WindowWide, 0, 0},
        {CandidateMode::WindowScored, 101, 250},
        {CandidateMode::WindowScored, 211, 250},
        {CandidateMode::WindowScored, 307, 250},
        {CandidateMode::WindowWide, 401, 250},
        {CandidateMode::WindowWide, 503, 250},
        {CandidateMode::WindowWide, 607, 250},
    };
    if (config.enable_distance_routes) {
        routes.push_back({CandidateMode::LocalDistance, 0, 0});
        routes.push_back({CandidateMode::LocalDistance, 503, 250});
        routes.push_back({CandidateMode::DistanceScored, 0, 0});
        routes.push_back({CandidateMode::DistanceScored, 701, 250});
    }

    Forest best;
    bool have_best = false;
    for (const Route& route : routes) {
        Forest candidate = run_one_route(forest, t1, t2, config, route.mode,
                                         route.seed, route.jitter);
        bool ok = validate_forest(candidate, t1, t2).ok;
        if (debug) {
            std::cerr << "route " << candidate_mode_name(route.mode)
                      << " seed=" << route.seed
                      << " jitter=" << route.jitter
                      << " components=" << candidate.active_count
                      << " valid=" << ok << "\n";
        }
        if (ok && (!have_best || candidate.active_count < best.active_count)) {
            best = std::move(candidate);
            if (config.publish_improvements) {
                publish_best_output(forest_to_output(best, t1));
            }
            have_best = true;
        }
    }
    if (have_best) {
        best = repair_by_conflict_targets(std::move(best), t1, t2, config, debug);
        best = repair_by_single_splits(std::move(best), t1, t2, config, debug);
        best = repair_by_conflict_targets(std::move(best), t1, t2, config, debug);
        return repair_by_pair_splits(std::move(best), t1, t2, config, debug);
    }
    return forest;
}

Forest run_split_restarts(Forest initial, const Tree& t1, const Tree& t2,
                          const GreedyConfig& config, int attempts,
                          int seed_offset, int extra_splits) {
    if (attempts <= 0) return initial;
    struct Route {
        CandidateMode mode;
        int seed;
        int jitter;
    };
    const Route routes[] = {
        {CandidateMode::WindowWide, 3001, 450},
        {CandidateMode::WindowScored, 4001, 350},
        {CandidateMode::WindowTight, 5003, 250},
        {CandidateMode::WindowWide, 6007, 700},
    };

    Forest best = std::move(initial);
    Forest current = best;
    int current_count = current.active_count;

    for (int attempt = 0; attempt < attempts; ++attempt) {
        std::vector<CompId> split_ids;
        if (attempt % 3 == 0) {
            std::vector<ConflictTarget> targets =
                conflict_targets(current, t1, t2, config.conflict_target_limit,
                                 config.conflict_max_split_components);
            if (!targets.empty()) {
                split_ids = targets[(attempt / 3 + seed_offset) % targets.size()].split;
            }
        }
        if (split_ids.empty()) {
            std::vector<CompId> ids = largest_components(current, config.repair_top_components);
            if (ids.empty()) break;
            split_ids.push_back(ids[(attempt + seed_offset) % ids.size()]);
            for (int extra = 0; extra < extra_splits && static_cast<int>(split_ids.size()) < static_cast<int>(ids.size()); ++extra) {
                CompId id = ids[(attempt * 17 + seed_offset + extra * 5) % ids.size()];
                if (std::find(split_ids.begin(), split_ids.end(), id) == split_ids.end()) {
                    split_ids.push_back(id);
                }
            }
        }

        Forest split = split_components_to_singletons(current, split_ids, t1, t2);
        const Route& route = routes[(attempt + seed_offset) % 4];
        Forest candidate = run_one_route(split, t1, t2, config, route.mode,
                                         route.seed + seed_offset * 1009 + attempt * 37,
                                         route.jitter);
        if ((attempt % 4) == 3) {
            Forest portfolio = run_greedy_merge(split, t1, t2, config);
            if (validate_forest(portfolio, t1, t2).ok &&
                portfolio.active_count < candidate.active_count) {
                candidate = std::move(portfolio);
            }
        }
        ValidationResult check = validate_forest(candidate, t1, t2);
        if (!check.ok) continue;

        int candidate_count = candidate.active_count;
        if (candidate_count <= current_count ||
            candidate_count <= best.active_count + 1) {
            current = candidate;
            current_count = candidate_count;
        }
        if (candidate_count < best.active_count) {
            best = std::move(candidate);
            if (config.publish_improvements) {
                publish_best_output(forest_to_output(best, t1));
            }
        }
    }
    return best;
}
