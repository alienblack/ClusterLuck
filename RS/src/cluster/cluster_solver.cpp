#include "cluster/cluster_solver.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "cluster/cluster_points.h"
#include "forest/validator.h"
#include "heuristic/greedy_merge.h"
#include "kernel/kernel_reduce.h"
#include "io/newick.h"
#include "io/pace_io.h"
#include "tree/induced.h"
#include "util/signal_handler.h"
#include "whidden/whidden_bridge.h"

namespace {

thread_local int g_cluster_solver_depth = 0;
thread_local int g_nested_local_calls = 0;

struct LocalSolution {
    std::vector<std::vector<std::string>> components;
    int score = 0;
};

struct BoundaryAttachStats {
    int candidates = 0;
    int topology_ok = 0;
    int edge_free = 0;
    int valid = 0;
    int group_candidates = 0;
    int group_valid = 0;
    int topology_fail = 0;
    int edge_blocked_first = 0;
    int edge_blocked_second = 0;
    int full_validation_fail = 0;
    int best_valid_size = -1;
    CompId best_valid_component = -1;
    int best_group_size = -1;
};

struct MutableLocalTree {
    struct Node {
        int parent = -1;
        int left = -1;
        int right = -1;
        int label = 0;
        bool alive = true;
    };

    int n = 0;
    int root = -1;
    std::vector<Node> nodes;
    std::vector<int> leaf_node;

    static MutableLocalTree from_tree(const Tree& tree) {
        MutableLocalTree out;
        out.n = tree.n;
        out.root = tree.root;
        out.nodes.resize(tree.node_count());
        out.leaf_node.assign(tree.n + 1, -1);
        for (int v = 0; v < tree.node_count(); ++v) {
            out.nodes[v].parent = tree.parent[v];
            out.nodes[v].left = tree.left[v];
            out.nodes[v].right = tree.right[v];
            out.nodes[v].label = tree.label[v];
            if (tree.label[v] > 0 && tree.label[v] <= tree.n) {
                out.leaf_node[tree.label[v]] = v;
            }
        }
        return out;
    }

    bool is_leaf(int node) const {
        return node >= 0 && nodes[node].alive &&
               nodes[node].left == -1 && nodes[node].right == -1;
    }

    int other_child(int parent_node, int child) const {
        if (nodes[parent_node].left == child) return nodes[parent_node].right;
        if (nodes[parent_node].right == child) return nodes[parent_node].left;
        return -1;
    }

    int leaf_node_of(int label) const {
        if (label <= 0 || label >= static_cast<int>(leaf_node.size())) return -1;
        int node = leaf_node[label];
        if (node < 0 || !nodes[node].alive) return -1;
        return node;
    }

    int cherry_sibling(int label) const {
        int node = leaf_node_of(label);
        if (node < 0) return -1;
        int parent_node = nodes[node].parent;
        if (parent_node < 0 || !nodes[parent_node].alive) return -1;
        int sibling = other_child(parent_node, node);
        if (!is_leaf(sibling)) return -1;
        return nodes[sibling].label;
    }

    void delete_leaf(int label) {
        int node = leaf_node_of(label);
        if (node < 0) return;
        int parent_node = nodes[node].parent;
        if (parent_node < 0) {
            nodes[node].alive = false;
            leaf_node[label] = -1;
            root = -1;
            return;
        }
        int sibling = other_child(parent_node, node);
        int grandparent = nodes[parent_node].parent;
        if (grandparent < 0) {
            root = sibling;
            nodes[sibling].parent = -1;
        } else {
            if (nodes[grandparent].left == parent_node) {
                nodes[grandparent].left = sibling;
            } else if (nodes[grandparent].right == parent_node) {
                nodes[grandparent].right = sibling;
            } else {
                throw std::runtime_error("bad local kernel parent pointer");
            }
            nodes[sibling].parent = grandparent;
        }
        nodes[node].alive = false;
        nodes[parent_node].alive = false;
        leaf_node[label] = -1;
    }

    std::string newick_rec(int node, const std::vector<int>& compact) const {
        if (node < 0 || !nodes[node].alive) return "";
        if (is_leaf(node)) return std::to_string(compact[nodes[node].label]);
        return "(" + newick_rec(nodes[node].left, compact) + "," +
               newick_rec(nodes[node].right, compact) + ")";
    }
};

struct LocalCherryKernel {
    bool changed = false;
    int reduced_n = 0;
    int reductions = 0;
    std::string first_newick;
    std::string second_newick;
    std::vector<std::vector<int>> payloads;
};

double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
}

LocalCherryKernel build_local_cherry_kernel(const Tree& local_t1,
                                            const Tree& local_t2,
                                            int solve_n,
                                            int protected_leaf) {
    MutableLocalTree t1 = MutableLocalTree::from_tree(local_t1);
    MutableLocalTree t2 = MutableLocalTree::from_tree(local_t2);
    std::vector<char> active(solve_n + 1, 1);
    std::vector<std::vector<int>> payloads(solve_n + 1);
    for (int leaf = 1; leaf <= solve_n; ++leaf) payloads[leaf] = {leaf};
    int active_count = solve_n;
    int reductions = 0;

    bool changed = true;
    while (changed && active_count > 1) {
        changed = false;
        for (int x = 1; x <= solve_n; ++x) {
            if (!active[x] || x == protected_leaf) continue;
            int sibling1 = t1.cherry_sibling(x);
            int sibling2 = t2.cherry_sibling(x);
            if (sibling1 < 0 || sibling1 != sibling2) continue;
            int y = sibling1;
            if (y == protected_leaf || y <= 0 || y > solve_n || !active[y]) continue;

            active[y] = 0;
            --active_count;
            ++reductions;
            payloads[x].insert(payloads[x].end(), payloads[y].begin(),
                               payloads[y].end());
            std::sort(payloads[x].begin(), payloads[x].end());
            payloads[x].erase(std::unique(payloads[x].begin(), payloads[x].end()),
                              payloads[x].end());
            t1.delete_leaf(y);
            t2.delete_leaf(y);
            changed = true;
            break;
        }
    }

    LocalCherryKernel kernel;
    kernel.changed = reductions > 0;
    kernel.reductions = reductions;
    if (!kernel.changed) return kernel;

    std::vector<int> compact(solve_n + 1, 0);
    kernel.payloads.push_back({});
    for (int leaf = 1; leaf <= solve_n; ++leaf) {
        if (!active[leaf]) continue;
        compact[leaf] = static_cast<int>(kernel.payloads.size());
        kernel.payloads.push_back(payloads[leaf]);
    }
    kernel.reduced_n = static_cast<int>(kernel.payloads.size()) - 1;
    kernel.first_newick = t1.newick_rec(t1.root, compact) + ";";
    kernel.second_newick = t2.newick_rec(t2.root, compact) + ";";
    return kernel;
}

Forest expand_local_kernel_forest(const Forest& reduced,
                                  const LocalCherryKernel& kernel,
                                  int original_n,
                                  const Tree& local_t1,
                                  const Tree& local_t2) {
    Forest expanded;
    expanded.leaf_comp.assign(original_n + 1, -1);
    std::vector<unsigned char> seen(original_n + 1, 0);
    for (CompId id : active_components(reduced)) {
        std::vector<int> leaves;
        for (int reduced_leaf : reduced.comps[id].leaves) {
            if (reduced_leaf <= 0 ||
                reduced_leaf >= static_cast<int>(kernel.payloads.size())) {
                throw std::runtime_error("bad reduced local kernel leaf");
            }
            leaves.insert(leaves.end(), kernel.payloads[reduced_leaf].begin(),
                          kernel.payloads[reduced_leaf].end());
        }
        std::sort(leaves.begin(), leaves.end());
        leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
        Component component;
        component.leaves = std::move(leaves);
        component.repr = component.leaves.front();
        for (int leaf : component.leaves) {
            if (leaf < 1 || leaf > original_n || seen[leaf]) {
                throw std::runtime_error("bad expanded local kernel cover");
            }
            seen[leaf] = 1;
            expanded.leaf_comp[leaf] = static_cast<CompId>(expanded.comps.size());
        }
        expanded.comps.push_back(std::move(component));
        ++expanded.active_count;
    }
    for (int leaf = 1; leaf <= original_n; ++leaf) {
        if (!seen[leaf]) throw std::runtime_error("missing expanded local kernel leaf");
    }
    if (!rebuild_forest_data(expanded, local_t1, local_t2)) {
        throw std::runtime_error("expanded local kernel rebuild failed");
    }
    return expanded;
}

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

bool env_int_list_contains(const char* name, int target) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return false;
    const char* cursor = raw;
    while (*cursor != '\0') {
        while (*cursor == ',' || std::isspace(static_cast<unsigned char>(*cursor))) {
            ++cursor;
        }
        if (*cursor == '\0') break;
        char* end = nullptr;
        long value = std::strtol(cursor, &end, 10);
        if (end == cursor) {
            while (*cursor != '\0' && *cursor != ',') ++cursor;
            continue;
        }
        if (value == target) return true;
        cursor = end;
    }
    return false;
}

bool contains_leaf(const Component& component, Leaf leaf) {
    return std::find(component.leaves.begin(), component.leaves.end(), leaf) !=
           component.leaves.end();
}

std::string tree_newick_from_node(const Tree& tree, NodeId node) {
    if (tree.is_leaf(node)) return std::to_string(tree.label[node]);
    return "(" + tree_newick_from_node(tree, tree.left[node]) + "," +
           tree_newick_from_node(tree, tree.right[node]) + ")";
}

std::string tree_to_newick(const Tree& tree) {
    return tree_newick_from_node(tree, tree.root) + ";";
}

bool build_forest_from_local_components(const Forest& source, int solve_n,
                                        const Tree& local_t1, const Tree& local_t2,
                                        CompId boundary_component,
                                        CompId attach_to) {
    Forest trial;
    trial.leaf_comp.assign(solve_n + 1, -1);
    for (CompId id : active_components(source)) {
        if (id == boundary_component) continue;
        Component component;
        component.leaves = source.comps[id].leaves;
        if (id == attach_to) component.leaves.push_back(solve_n);
        std::sort(component.leaves.begin(), component.leaves.end());
        component.leaves.erase(std::unique(component.leaves.begin(),
                                           component.leaves.end()),
                               component.leaves.end());
        for (Leaf leaf : component.leaves) {
            if (leaf < 1 || leaf > solve_n) return false;
            if (trial.leaf_comp[leaf] != -1) return false;
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }
    for (Leaf leaf = 1; leaf <= solve_n; ++leaf) {
        if (trial.leaf_comp[leaf] == -1) return false;
    }
    return rebuild_forest_data(trial, local_t1, local_t2) &&
           validate_forest(trial, local_t1, local_t2).ok;
}

bool build_forest_from_boundary_group(const Forest& source, int solve_n,
                                      const Tree& local_t1,
                                      const Tree& local_t2,
                                      CompId boundary_component,
                                      const std::vector<CompId>& group) {
    std::vector<unsigned char> in_group(source.comps.size(), 0);
    std::vector<Leaf> merged_leaves;
    merged_leaves.push_back(solve_n);
    for (CompId id : group) {
        if (id < 0 || id >= static_cast<CompId>(source.comps.size()) ||
            !source.comps[id].active || id == boundary_component) {
            return false;
        }
        in_group[id] = 1;
        merged_leaves.insert(merged_leaves.end(), source.comps[id].leaves.begin(),
                             source.comps[id].leaves.end());
    }
    std::sort(merged_leaves.begin(), merged_leaves.end());
    merged_leaves.erase(std::unique(merged_leaves.begin(), merged_leaves.end()),
                        merged_leaves.end());

    Forest trial;
    trial.leaf_comp.assign(solve_n + 1, -1);
    for (CompId id : active_components(source)) {
        if (id == boundary_component) continue;
        if (id >= 0 && id < static_cast<CompId>(in_group.size()) && in_group[id]) {
            continue;
        }
        Component component;
        component.leaves = source.comps[id].leaves;
        for (Leaf leaf : component.leaves) {
            if (leaf < 1 || leaf > solve_n) return false;
            if (trial.leaf_comp[leaf] != -1) return false;
            trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
        }
        trial.comps.push_back(std::move(component));
        ++trial.active_count;
    }

    Component boundary_group;
    boundary_group.leaves = std::move(merged_leaves);
    for (Leaf leaf : boundary_group.leaves) {
        if (leaf < 1 || leaf > solve_n) return false;
        if (trial.leaf_comp[leaf] != -1) return false;
        trial.leaf_comp[leaf] = static_cast<CompId>(trial.comps.size());
    }
    trial.comps.push_back(std::move(boundary_group));
    ++trial.active_count;

    for (Leaf leaf = 1; leaf <= solve_n; ++leaf) {
        if (trial.leaf_comp[leaf] == -1) return false;
    }
    return rebuild_forest_data(trial, local_t1, local_t2) &&
           validate_forest(trial, local_t1, local_t2).ok;
}

std::vector<CompId> external_edge_blockers(const std::vector<NodeId>& edges,
                                           const std::vector<CompId>& owner,
                                           CompId allowed_a,
                                           CompId allowed_b) {
    std::vector<CompId> blockers;
    for (NodeId edge : edges) {
        if (edge < 0 || edge >= static_cast<NodeId>(owner.size())) continue;
        CompId current = owner[edge];
        if (current != -1 && current != allowed_a && current != allowed_b) {
            blockers.push_back(current);
        }
    }
    std::sort(blockers.begin(), blockers.end());
    blockers.erase(std::unique(blockers.begin(), blockers.end()), blockers.end());
    return blockers;
}

int count_external_edge_blockers(const std::vector<NodeId>& edges,
                                 const std::vector<CompId>& owner,
                                 CompId allowed_a, CompId allowed_b) {
    return static_cast<int>(
        external_edge_blockers(edges, owner, allowed_a, allowed_b).size());
}

bool probe_boundary_attach(const Forest& source, int solve_n,
                           const Tree& local_t1, const Tree& local_t2,
                           CompId boundary_component, CompId attach_to,
                           BoundaryAttachStats& stats) {
    ++stats.candidates;
    std::vector<Leaf> leaves = source.comps[attach_to].leaves;
    leaves.push_back(solve_n);
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());

    InducedInfo i1 = compute_induced(local_t1, leaves);
    InducedInfo i2 = compute_induced(local_t2, leaves);
    if (i1.signature != i2.signature) {
        ++stats.topology_fail;
        return false;
    }
    ++stats.topology_ok;

    int blocked_first = count_external_edge_blockers(
        i1.used_edges, source.owner1, boundary_component, attach_to);
    int blocked_second = count_external_edge_blockers(
        i2.used_edges, source.owner2, boundary_component, attach_to);
    if (blocked_first > 0) ++stats.edge_blocked_first;
    if (blocked_second > 0) ++stats.edge_blocked_second;
    if (blocked_first > 0 || blocked_second > 0) return false;
    ++stats.edge_free;

    bool valid = build_forest_from_local_components(
        source, solve_n, local_t1, local_t2, boundary_component, attach_to);
    if (!valid) {
        ++stats.full_validation_fail;
        return false;
    }
    ++stats.valid;
    int size = static_cast<int>(source.comps[attach_to].leaves.size());
    if (size > stats.best_valid_size) {
        stats.best_valid_size = size;
        stats.best_valid_component = attach_to;
    }
    return true;
}

CompId local_root_component(const Forest& forest, const Tree& local_t1,
                            const Tree& local_t2) {
    for (CompId id : active_components(forest)) {
        const Component& component = forest.comps[id];
        if (component.root1 == local_t1.root && component.root2 == local_t2.root) {
            return id;
        }
    }
    for (CompId id : active_components(forest)) {
        if (forest.comps[id].root1 == local_t1.root) return id;
    }
    std::vector<CompId> ids = active_components(forest);
    return ids.empty() ? -1 : ids.front();
}

int count_cluster_marker_payloads(const LocalProblem& problem,
                                  const KernelMap* kernel_map = nullptr) {
    if (kernel_map == nullptr) {
        int count = 0;
        for (size_t i = 1; i < problem.atom_names.size(); ++i) {
            const std::string& atom = problem.atom_names[i];
            if (!atom.empty() && atom[0] == 'C') ++count;
        }
        return count;
    }

    int count = 0;
    for (const auto& entry : kernel_map->active_payload) {
        bool has_marker = false;
        for (int local_leaf : entry.second) {
            if (local_leaf <= 0 ||
                local_leaf >= static_cast<int>(problem.atom_names.size())) {
                continue;
            }
            const std::string& atom = problem.atom_names[local_leaf];
            if (!atom.empty() && atom[0] == 'C') {
                has_marker = true;
                break;
            }
        }
        if (has_marker) ++count;
    }
    return count;
}

Forest solve_numeric_problem(int solve_n, int cluster_id, bool residual,
                             const Tree& local_t1, const Tree& local_t2,
                             int cluster_marker_count = 0) {
    const bool debug = std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr;
    const auto solve_start = std::chrono::steady_clock::now();
    Forest forest = make_singleton_forest(solve_n);
    if (!rebuild_forest_data(forest, local_t1, local_t2)) {
        throw std::runtime_error("local singleton forest invalid");
    }

    const int score_mode = env_int("RS_LOCAL_SCORE_MODE", 2);
    Forest singleton_raw =
        run_h45_greedy_merge(std::move(forest), local_t1, local_t2,
                             score_mode, false, false);
    Forest candidate = exact_output_repair(singleton_raw, local_t1, local_t2);
    if (debug) {
        std::cerr << "rs_local_stage id=" << cluster_id
                  << " stage=singleton components=" << candidate.active_count
                  << " seconds=" << seconds_since(solve_start) << '\n';
    }

    const int default_restart1_budget = residual
        ? env_int("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", 512)
        : env_int("RS_LOCAL_CLUSTER_RESTART1_BUDGET", 512);
    const int default_restart2_budget = residual
        ? env_int("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", 512)
        : env_int("RS_LOCAL_CLUSTER_RESTART2_BUDGET", 0);
    const int restart1_budget =
        std::getenv("RS_LOCAL_RESTART1_BUDGET") != nullptr
            ? env_int("RS_LOCAL_RESTART1_BUDGET", default_restart1_budget)
            : default_restart1_budget;
    const int restart2_budget =
        std::getenv("RS_LOCAL_RESTART2_BUDGET") != nullptr
            ? env_int("RS_LOCAL_RESTART2_BUDGET", default_restart2_budget)
            : default_restart2_budget;
    const int restart2_min_leaves =
        std::getenv("RS_LOCAL_RESTART2_MIN_LEAVES") != nullptr
            ? env_int("RS_LOCAL_RESTART2_MIN_LEAVES", 0)
            : (residual
                   ? env_int("RS_LOCAL_RESIDUAL_RESTART2_MIN_LEAVES", 0)
                   : env_int("RS_LOCAL_CLUSTER_RESTART2_MIN_LEAVES", 0));
    const int default_portfolio_budget = residual
        ? env_int("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", 384)
        : env_int("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", 192);
    const int portfolio_budget =
        std::getenv("RS_LOCAL_PORTFOLIO_BUDGET") != nullptr
            ? env_int("RS_LOCAL_PORTFOLIO_BUDGET", default_portfolio_budget)
            : default_portfolio_budget;
    const int default_portfolio_min_leaves = residual ? 150 : 50;
    const int portfolio_min_leaves =
        std::getenv("RS_LOCAL_PORTFOLIO_MIN_LEAVES") != nullptr
            ? env_int("RS_LOCAL_PORTFOLIO_MIN_LEAVES",
                      default_portfolio_min_leaves)
            : (residual
                   ? env_int("RS_LOCAL_RESIDUAL_PORTFOLIO_MIN_LEAVES",
                             default_portfolio_min_leaves)
                   : env_int("RS_LOCAL_CLUSTER_PORTFOLIO_MIN_LEAVES",
                             default_portfolio_min_leaves));
    int effective_restart1_budget = restart1_budget;
    int effective_restart2_budget = restart2_budget;
    int effective_portfolio_budget = portfolio_budget;
    if (env_int("RS_LOCAL_COMPACT_BUDGETS", 0) != 0) {
        const int threshold = env_int("RS_LOCAL_COMPACT_BUDGET_THRESHOLD", 200);
        if (solve_n > threshold) {
            effective_restart1_budget =
                env_int("RS_LOCAL_COMPACT_BIG_RESTART1", 128);
            effective_restart2_budget =
                env_int("RS_LOCAL_COMPACT_BIG_RESTART2", 0);
        } else {
            effective_restart1_budget =
                env_int("RS_LOCAL_COMPACT_SMALL_RESTART1", 512);
            effective_restart2_budget =
                env_int("RS_LOCAL_COMPACT_SMALL_RESTART2", 512);
        }
        effective_portfolio_budget =
            env_int("RS_LOCAL_COMPACT_PORTFOLIO_BUDGET", 0);
    }
    bool child_marker_light = false;
    if (!residual && env_int("RS_CHILD_MARKER_LIGHT_RESTART", 0) != 0) {
        const int marker_threshold =
            env_int("RS_CHILD_MARKER_LIGHT_RESTART_THRESHOLD", 25);
        const int solve_n_threshold =
            env_int("RS_CHILD_MARKER_LIGHT_RESTART_MIN_SOLVE_N", 1000);
        if (solve_n > solve_n_threshold &&
            cluster_marker_count > marker_threshold) {
            child_marker_light = true;
            effective_restart1_budget =
                env_int("RS_CHILD_MARKER_LIGHT_RESTART_BUDGET", 256);
            effective_restart2_budget = 0;
            effective_portfolio_budget = 0;
        }
    }

    auto stage_start = std::chrono::steady_clock::now();
    const int restart1_mode = env_int("RS_LOCAL_RESTART1_MODE", score_mode);
    const int restart2_mode = env_int("RS_LOCAL_RESTART2_MODE", score_mode);
    Forest restarted = run_h45_split_restarts(singleton_raw, local_t1, local_t2,
                                              restart1_mode, effective_restart1_budget,
                                              cluster_id * 17, 0);
    if (effective_restart1_budget > 0) {
        restarted = exact_output_repair(std::move(restarted), local_t1, local_t2);
        if (validate_forest(restarted, local_t1, local_t2).ok &&
            restarted.active_count < candidate.active_count) {
            candidate = std::move(restarted);
        }
    }
    if (debug) {
        std::cerr << "rs_local_stage id=" << cluster_id
                  << " stage=restart1 budget=" << effective_restart1_budget
                  << " base_budget=" << restart1_budget
                  << " child_markers=" << cluster_marker_count
                  << " child_marker_light=" << (child_marker_light ? 1 : 0)
                  << " mode=" << restart1_mode
                  << " best=" << candidate.active_count
                  << " seconds=" << seconds_since(stage_start) << '\n';
    }

    stage_start = std::chrono::steady_clock::now();
    const bool run_restart2 =
        effective_restart2_budget > 0 && solve_n >= restart2_min_leaves;
    if (run_restart2) {
        restarted = run_h45_split_restarts(candidate, local_t1, local_t2,
                                          restart2_mode, effective_restart2_budget,
                                          193 + cluster_id * 31, 0);
        restarted = exact_output_repair(std::move(restarted), local_t1, local_t2);
        if (validate_forest(restarted, local_t1, local_t2).ok &&
            restarted.active_count < candidate.active_count) {
            candidate = std::move(restarted);
        }
    }
    if (debug) {
        std::cerr << "rs_local_stage id=" << cluster_id
                  << " stage=restart2 budget=" << effective_restart2_budget
                  << " base_budget=" << restart2_budget
                  << " min_leaves=" << restart2_min_leaves
                  << " ran=" << (run_restart2 ? 1 : 0)
                  << " mode=" << restart2_mode
                  << " best=" << candidate.active_count
                  << " seconds=" << seconds_since(stage_start) << '\n';
    }

    if (solve_n >= portfolio_min_leaves && effective_portfolio_budget > 0) {
        stage_start = std::chrono::steady_clock::now();
        struct Variant {
            int seed_offset;
            int mode;
            bool start_from_best;
        };
        const Variant variants[] = {
            {701 + cluster_id * 37, score_mode, true},
            {991 + cluster_id * 53, 2, false},
            {1301 + cluster_id * 71, 5, true},
        };

        for (const Variant& variant : variants) {
            Forest seed = variant.start_from_best ? candidate : singleton_raw;
            Forest trial = run_h45_split_restarts(std::move(seed), local_t1, local_t2,
                                                  variant.mode, effective_portfolio_budget,
                                                  variant.seed_offset, 0);
            trial = exact_output_repair(std::move(trial), local_t1, local_t2);
            if (validate_forest(trial, local_t1, local_t2).ok &&
                trial.active_count < candidate.active_count) {
                candidate = std::move(trial);
            }
        }
        if (debug) {
            std::cerr << "rs_local_stage id=" << cluster_id
                      << " stage=portfolio budget=" << effective_portfolio_budget
                      << " base_budget=" << portfolio_budget
                      << " best=" << candidate.active_count
                      << " seconds=" << seconds_since(stage_start) << '\n';
        }
    }
    if (env_int("RS_NESTED_LOCAL", 0) != 0 && g_cluster_solver_depth == 1) {
        const int nested_min_n = env_int("RS_NESTED_LOCAL_MIN_N", 70);
        const int nested_max_n = env_int("RS_NESTED_LOCAL_MAX_N", 250);
        const int nested_max_calls = env_int("RS_NESTED_LOCAL_MAX_CALLS", 3);
        const bool include_residual =
            env_int("RS_NESTED_LOCAL_INCLUDE_RESIDUAL", 0) != 0;
        const bool id_selected =
            env_int_list_contains("RS_NESTED_LOCAL_IDS", cluster_id) ||
            std::getenv("RS_NESTED_LOCAL_IDS") == nullptr;
        if ((!residual || include_residual) &&
            solve_n >= nested_min_n &&
            (nested_max_n <= 0 || solve_n <= nested_max_n) &&
            g_nested_local_calls < nested_max_calls &&
            id_selected) {
            ++g_nested_local_calls;
            stage_start = std::chrono::steady_clock::now();
            Forest nested;
            bool nested_ok = run_cluster_solver(local_t1, local_t2, nested);
            if (nested_ok) {
                nested = exact_output_repair(std::move(nested), local_t1, local_t2);
                ValidationResult nested_check =
                    validate_forest(nested, local_t1, local_t2);
                if (nested_check.ok && nested.active_count < candidate.active_count) {
                    if (debug) {
                        std::cerr << "rs_local_stage id=" << cluster_id
                                  << " stage=nested_rs improved=1 previous="
                                  << candidate.active_count
                                  << " best=" << nested.active_count
                                  << " seconds=" << seconds_since(stage_start)
                                  << '\n';
                    }
                    candidate = std::move(nested);
                } else if (debug) {
                    std::cerr << "rs_local_stage id=" << cluster_id
                              << " stage=nested_rs improved=0 candidate="
                              << (nested_check.ok ? nested.active_count : -1)
                              << " best=" << candidate.active_count
                              << " valid=" << (nested_check.ok ? 1 : 0)
                              << " seconds=" << seconds_since(stage_start)
                              << '\n';
                }
            } else if (debug) {
                std::cerr << "rs_local_stage id=" << cluster_id
                          << " stage=nested_rs skipped=1 seconds="
                          << seconds_since(stage_start) << '\n';
            }
        }
    }
    ValidationResult check = validate_forest(candidate, local_t1, local_t2);
    if (!check.ok) {
        throw std::runtime_error("local solve invalid: " + check.message);
    }
    if (debug) {
        std::cerr << "rs_local_stage id=" << cluster_id
                  << " stage=total components=" << candidate.active_count
                  << " seconds=" << seconds_since(solve_start) << '\n';
    }
    return candidate;
}

std::vector<std::string> component_atoms(const Component& component,
                                         const LocalProblem& problem,
                                         int boundary_leaf) {
    std::vector<std::string> atoms;
    for (Leaf leaf : component.leaves) {
        if (leaf == boundary_leaf) continue;
        if (leaf <= 0 || leaf >= static_cast<int>(problem.atom_names.size())) {
            throw std::runtime_error("local leaf outside atom map");
        }
        std::string atom = problem.atom_names[leaf];
        if (!atom.empty() && atom[0] != 'C' && atom != "__BOUNDARY__") {
            char* end = nullptr;
            long internal = std::strtol(atom.c_str(), &end, 10);
            if (end != atom.c_str() && *end == '\0' && internal >= 0 &&
                internal < static_cast<long>(problem.internal_to_original.size())) {
                atom = std::to_string(problem.internal_to_original[internal]);
            }
        }
        atoms.push_back(std::move(atom));
    }
    std::sort(atoms.begin(), atoms.end());
    return atoms;
}

Forest split_leaf_to_singleton(Forest source, Leaf leaf,
                               const Tree& local_t1,
                               const Tree& local_t2) {
    Forest out;
    out.leaf_comp.assign(local_t1.n + 1, -1);
    for (CompId id : active_components(source)) {
        const Component& original = source.comps[id];
        if (!contains_leaf(original, leaf)) {
            Component component;
            component.leaves = original.leaves;
            out.comps.push_back(std::move(component));
            ++out.active_count;
            continue;
        }

        std::vector<Leaf> rest;
        rest.reserve(original.leaves.size());
        for (Leaf x : original.leaves) {
            if (x != leaf) rest.push_back(x);
        }
        if (!rest.empty()) {
            Component component;
            component.leaves = std::move(rest);
            out.comps.push_back(std::move(component));
            ++out.active_count;
        }
        Component singleton;
        singleton.leaves.push_back(leaf);
        out.comps.push_back(std::move(singleton));
        ++out.active_count;
    }
    rebuild_forest_data(out, local_t1, local_t2);
    return out;
}

LocalSolution solve_prebuilt_local_route(const LocalProblem& problem,
                                         int cluster_id,
                                         bool attach_boundary_marker,
                                         bool residual,
                                         int route_seed_id) {
    Tree local_t1 =
        parse_newick_tree_h45_order(problem.first_newick, problem.solve_n);
    Tree local_t2 =
        parse_newick_tree_h45_order(problem.second_newick, problem.solve_n);
    Forest forest;
    const int boundary_leaf_for_kernel =
        attach_boundary_marker ? problem.solve_n : -1;
    const bool use_local_full_kernel =
        env_int("RS_LOCAL_FULL_KERNEL", 0) != 0 &&
        problem.solve_n >= env_int("RS_LOCAL_FULL_KERNEL_MIN_N", 2);
    const bool use_local_cherry_kernel =
        env_int("RS_LOCAL_CHERRY_KERNEL", 0) != 0 &&
        problem.solve_n >= env_int("RS_LOCAL_CHERRY_KERNEL_MIN_N", 2);
    if (use_local_full_kernel) {
        Instance local_instance;
        local_instance.n = problem.solve_n;
        local_instance.t1 = local_t1;
        local_instance.t2 = local_t2;
        auto kernel = kernel_reduce_instance(
            local_instance,
            std::getenv("RS_LOCAL_FULL_KERNEL_VERBOSE") != nullptr,
            boundary_leaf_for_kernel);
        if (kernel.has_value()) {
            int cluster_marker_count =
                count_cluster_marker_payloads(problem, &kernel->map);
            Forest reduced_forest = solve_numeric_problem(
                kernel->reduced.n, route_seed_id, residual, kernel->reduced.t1,
                kernel->reduced.t2, cluster_marker_count);
            forest = expand_kernel_reduced_forest(reduced_forest, kernel->map,
                                                  local_instance);
            ValidationResult kernel_check =
                validate_forest(forest, local_t1, local_t2);
            if (!kernel_check.ok) {
                throw std::runtime_error("local full kernel invalid: " +
                                         kernel_check.message);
            }
            if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                std::cerr << "rs_local_full_kernel id=" << cluster_id
                          << " route_id=" << route_seed_id
                          << " original_n=" << problem.solve_n
                          << " reduced_n=" << kernel->reduced.n
                          << " subtree=" << kernel->subtree_deleted
                          << " chain=" << kernel->chain_deleted
                          << " three_two=" << kernel->param_reduction
                          << " marker_payloads=" << cluster_marker_count
                          << " components=" << forest.active_count << '\n';
            }
        } else {
            int cluster_marker_count = count_cluster_marker_payloads(problem);
            forest = solve_numeric_problem(problem.solve_n, route_seed_id,
                                           residual, local_t1, local_t2,
                                           cluster_marker_count);
        }
    } else if (use_local_cherry_kernel) {
        LocalCherryKernel kernel = build_local_cherry_kernel(
            local_t1, local_t2, problem.solve_n, boundary_leaf_for_kernel);
        if (kernel.changed) {
            Tree reduced_t1 =
                parse_newick_tree_h45_order(kernel.first_newick, kernel.reduced_n);
            Tree reduced_t2 =
                parse_newick_tree_h45_order(kernel.second_newick, kernel.reduced_n);
            int cluster_marker_count = count_cluster_marker_payloads(problem);
            Forest reduced_forest = solve_numeric_problem(
                kernel.reduced_n, route_seed_id, residual, reduced_t1, reduced_t2,
                cluster_marker_count);
            forest = expand_local_kernel_forest(
                reduced_forest, kernel, problem.solve_n, local_t1, local_t2);
            ValidationResult kernel_check =
                validate_forest(forest, local_t1, local_t2);
            if (!kernel_check.ok) {
                throw std::runtime_error("local cherry kernel invalid: " +
                                         kernel_check.message);
            }
            if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                std::cerr << "rs_local_cherry_kernel id=" << cluster_id
                          << " route_id=" << route_seed_id
                          << " original_n=" << problem.solve_n
                          << " reduced_n=" << kernel.reduced_n
                          << " reductions=" << kernel.reductions
                          << " components=" << forest.active_count << '\n';
            }
        } else {
            int cluster_marker_count = count_cluster_marker_payloads(problem);
            forest = solve_numeric_problem(problem.solve_n, route_seed_id,
                                           residual, local_t1, local_t2,
                                           cluster_marker_count);
        }
    } else {
        int cluster_marker_count = count_cluster_marker_payloads(problem);
        forest = solve_numeric_problem(problem.solve_n, route_seed_id,
                                       residual, local_t1, local_t2,
                                       cluster_marker_count);
    }

    LocalSolution solution;
    const int debug_id = attach_boundary_marker ? cluster_id : route_seed_id;
    const int boundary_leaf = attach_boundary_marker ? problem.solve_n : -1;
    const int force_singleton_cluster =
        env_int("RS_FORCE_BOUNDARY_SINGLETON_CLUSTER", -1);
    const bool force_boundary_singleton =
        attach_boundary_marker &&
        (cluster_id == force_singleton_cluster ||
         env_int_list_contains("RS_FORCE_BOUNDARY_SINGLETON_CLUSTERS",
                               cluster_id));
    if (force_boundary_singleton) {
        int before_size = -1;
        for (CompId id : active_components(forest)) {
            if (contains_leaf(forest.comps[id], boundary_leaf)) {
                before_size = static_cast<int>(forest.comps[id].leaves.size());
                break;
            }
        }
        forest = split_leaf_to_singleton(std::move(forest), boundary_leaf,
                                         local_t1, local_t2);
        ValidationResult forced_check =
            validate_forest(forest, local_t1, local_t2);
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
            std::cerr << "rs_forced_boundary_singleton id=" << cluster_id
                      << " before_size=" << before_size
                      << " components=" << forest.active_count
                      << " valid=" << (forced_check.ok ? 1 : 0)
                      << " reason=" << forced_check.message << '\n';
        }
    }
    CompId boundary_component = -1;
    CompId root_component = -1;
    std::vector<CompId> boundary_group;
    if (attach_boundary_marker) {
        for (CompId id : active_components(forest)) {
            if (contains_leaf(forest.comps[id], boundary_leaf)) {
                boundary_component = id;
                break;
            }
        }
        if (boundary_component < 0) {
            throw std::runtime_error("boundary component missing");
        }

        root_component = boundary_component;
        if (forest.comps[boundary_component].leaves.size() == 1 &&
            !force_boundary_singleton) {
            const bool debug_boundary_attach =
                std::getenv("RS_DEBUG_BOUNDARY_ATTACH") != nullptr;
            const bool run_group_attach =
                env_int("RS_BOUNDARY_GROUP_ATTACH", 0) != 0;
            BoundaryAttachStats attach_stats;
            int best_size = -1;
            for (CompId id : active_components(forest)) {
                if (id == boundary_component) continue;
                bool valid_attach = false;
                if (debug_boundary_attach || run_group_attach) {
                    valid_attach = probe_boundary_attach(
                        forest, problem.solve_n, local_t1, local_t2,
                        boundary_component, id, attach_stats);
                } else {
                    valid_attach = build_forest_from_local_components(
                        forest, problem.solve_n, local_t1, local_t2,
                        boundary_component, id);
                }
                if (valid_attach) {
                    int size = static_cast<int>(forest.comps[id].leaves.size());
                    if (size > best_size) {
                        best_size = size;
                        root_component = id;
                    }
                }
            }
            if (root_component == boundary_component && run_group_attach) {
                const int max_group =
                    env_int("RS_BOUNDARY_GROUP_ATTACH_MAX_COMPONENTS", 8);
                int best_group_count = -1;
                int best_group_leaf_count = -1;
                for (CompId id : active_components(forest)) {
                    if (id == boundary_component) continue;
                    std::vector<Leaf> leaves = forest.comps[id].leaves;
                    leaves.push_back(problem.solve_n);
                    std::sort(leaves.begin(), leaves.end());
                    leaves.erase(std::unique(leaves.begin(), leaves.end()),
                                 leaves.end());
                    InducedInfo i1 = compute_induced(local_t1, leaves);
                    InducedInfo i2 = compute_induced(local_t2, leaves);
                    if (i1.signature != i2.signature) continue;

                    std::vector<CompId> group;
                    group.push_back(id);
                    std::vector<CompId> b1 = external_edge_blockers(
                        i1.used_edges, forest.owner1, boundary_component, id);
                    std::vector<CompId> b2 = external_edge_blockers(
                        i2.used_edges, forest.owner2, boundary_component, id);
                    group.insert(group.end(), b1.begin(), b1.end());
                    group.insert(group.end(), b2.begin(), b2.end());
                    std::sort(group.begin(), group.end());
                    group.erase(std::unique(group.begin(), group.end()),
                                group.end());
                    if (static_cast<int>(group.size()) <= 1 ||
                        static_cast<int>(group.size()) > max_group) {
                        continue;
                    }
                    ++attach_stats.group_candidates;
                    if (!build_forest_from_boundary_group(
                            forest, problem.solve_n, local_t1, local_t2,
                            boundary_component, group)) {
                        continue;
                    }
                    ++attach_stats.group_valid;
                    int leaf_count = 0;
                    for (CompId member : group) {
                        leaf_count +=
                            static_cast<int>(forest.comps[member].leaves.size());
                    }
                    if (static_cast<int>(group.size()) > best_group_count ||
                        (static_cast<int>(group.size()) == best_group_count &&
                         leaf_count > best_group_leaf_count)) {
                        best_group_count = static_cast<int>(group.size());
                        best_group_leaf_count = leaf_count;
                        boundary_group = std::move(group);
                    }
                }
                attach_stats.best_group_size = best_group_count;
                if (!boundary_group.empty()) {
                    root_component = -1;
                }
            }
            if (debug_boundary_attach) {
                std::cerr << "rs_boundary_attach id=" << cluster_id
                          << " local_components=" << forest.active_count
                          << " candidates=" << attach_stats.candidates
                          << " topology_ok=" << attach_stats.topology_ok
                          << " edge_free=" << attach_stats.edge_free
                          << " valid=" << attach_stats.valid
                          << " group_candidates="
                          << attach_stats.group_candidates
                          << " group_valid=" << attach_stats.group_valid
                          << " topology_fail=" << attach_stats.topology_fail
                          << " edge_blocked_first="
                          << attach_stats.edge_blocked_first
                          << " edge_blocked_second="
                          << attach_stats.edge_blocked_second
                          << " full_validation_fail="
                          << attach_stats.full_validation_fail
                          << " selected_size=" << best_size
                          << " selected_group_size="
                          << attach_stats.best_group_size << '\n';
            }
        }

        std::vector<std::string> boundary_atoms;
        if (!boundary_group.empty()) {
            Component merged_boundary;
            for (CompId id : boundary_group) {
                merged_boundary.leaves.insert(
                    merged_boundary.leaves.end(), forest.comps[id].leaves.begin(),
                    forest.comps[id].leaves.end());
            }
            std::sort(merged_boundary.leaves.begin(),
                      merged_boundary.leaves.end());
            merged_boundary.leaves.erase(
                std::unique(merged_boundary.leaves.begin(),
                            merged_boundary.leaves.end()),
                merged_boundary.leaves.end());
            boundary_atoms =
                component_atoms(merged_boundary, problem, boundary_leaf);
        } else {
            boundary_atoms =
                component_atoms(forest.comps[root_component], problem, boundary_leaf);
        }
        solution.components.push_back(std::move(boundary_atoms));
    } else {
        root_component = local_root_component(forest, local_t1, local_t2);
        if (root_component < 0) {
            throw std::runtime_error("local root component missing");
        }
        std::vector<std::string> root_atoms =
            component_atoms(forest.comps[root_component], problem, boundary_leaf);
        if (!root_atoms.empty()) solution.components.push_back(std::move(root_atoms));
    }

    std::vector<unsigned char> skip_group(forest.comps.size(), 0);
    for (CompId id : boundary_group) {
        if (id >= 0 && id < static_cast<CompId>(skip_group.size())) {
            skip_group[id] = 1;
        }
    }
    for (CompId id : active_components(forest)) {
        if (id == boundary_component || id == root_component) continue;
        if (id >= 0 && id < static_cast<CompId>(skip_group.size()) &&
            skip_group[id]) {
            continue;
        }
        if (attach_boundary_marker && contains_leaf(forest.comps[id], boundary_leaf)) {
            continue;
        }
        std::vector<std::string> atoms =
            component_atoms(forest.comps[id], problem, boundary_leaf);
        if (!atoms.empty()) solution.components.push_back(std::move(atoms));
    }
    solution.score = static_cast<int>(solution.components.size());
    if (std::getenv("RS_DUMP_LOCAL_SOLUTIONS") != nullptr &&
        (debug_id == 2 || debug_id == 3 || debug_id == 13 ||
         debug_id == 15)) {
        std::cerr << "rs_local_solution_begin id=" << debug_id
                  << " components=" << solution.components.size()
                  << " boundary_atoms="
                  << (solution.components.empty() ? 0 : solution.components[0].size())
                  << '\n';
        for (size_t i = 0; i < solution.components.size(); ++i) {
            std::cerr << "rs_local_component id=" << debug_id
                      << " component=" << i
                      << " atom_count=" << solution.components[i].size()
                      << " atoms=";
            for (const std::string& atom : solution.components[i]) {
                std::cerr << atom << ',';
            }
            std::cerr << '\n';
        }
        std::cerr << "rs_local_solution_end id=" << debug_id << '\n';
    }
    return solution;
}

LocalSolution solve_local_route(const Tree& t1, const Tree& t2,
                                const ClusterPlan& plan, int cluster_id,
                                bool attach_boundary_marker,
                                int route_seed_id) {
    LocalProblem problem =
        build_local_problem(t1, t2, plan, cluster_id, attach_boundary_marker);
    return solve_prebuilt_local_route(problem, cluster_id, attach_boundary_marker,
                                      !attach_boundary_marker, route_seed_id);
}

LocalSolution solve_local_route_with_suppressed(
    const Tree& t1, const Tree& t2, const ClusterPlan& plan, int cluster_id,
    bool attach_boundary_marker, int route_seed_id,
    const std::vector<unsigned char>& suppressed_clusters) {
    LocalProblem problem = build_local_problem_with_open_children_and_suppressed(
        t1, t2, plan, cluster_id, attach_boundary_marker, {},
        suppressed_clusters);
    if (std::getenv("RS_DEBUG_BOTTOMUP_MARKERS") != nullptr) {
        int cluster_atoms = 0;
        int ordinary_atoms = 0;
        for (size_t i = 1; i < problem.atom_names.size(); ++i) {
            const std::string& atom = problem.atom_names[i];
            if (atom == "__BOUNDARY__") continue;
            if (!atom.empty() && atom[0] == 'C') {
                ++cluster_atoms;
            } else {
                ++ordinary_atoms;
            }
        }
        int suppressed_count = 0;
        for (unsigned char value : suppressed_clusters) {
            if (value) ++suppressed_count;
        }
        std::cerr << "rs_bottomup_problem id=" << cluster_id
                  << " route_id=" << route_seed_id
                  << " attach_boundary=" << (attach_boundary_marker ? 1 : 0)
                  << " solve_n=" << problem.solve_n
                  << " cluster_markers=" << cluster_atoms
                  << " ordinary_atoms=" << ordinary_atoms
                  << " suppressed_count=" << suppressed_count << '\n';
    }
    return solve_prebuilt_local_route(problem, cluster_id, attach_boundary_marker,
                                      !attach_boundary_marker, route_seed_id);
}

bool consider_probe_candidate(Forest candidate, const Tree& t1, const Tree& t2,
                              Forest& best, const char* label, bool debug,
                              std::chrono::steady_clock::time_point start) {
    candidate = exact_output_repair(std::move(candidate), t1, t2);
    ValidationResult check = validate_forest(candidate, t1, t2);
    if (debug) {
        std::cerr << "rs_whole_probe_stage label=" << label
                  << " components=" << candidate.active_count
                  << " valid=" << (check.ok ? 1 : 0)
                  << " seconds=" << seconds_since(start) << '\n';
    }
    if (check.ok && candidate.active_count < best.active_count) {
        best = std::move(candidate);
        return true;
    }
    return false;
}

bool forest_from_leaf_components(const std::vector<std::vector<int>>& components,
                                 int n, const Tree& t1, const Tree& t2,
                                 Forest& forest) {
    forest = Forest{};
    forest.leaf_comp.assign(n + 1, -1);
    for (const std::vector<int>& leaves : components) {
        if (leaves.empty()) continue;
        Component component;
        component.leaves = leaves;
        std::sort(component.leaves.begin(), component.leaves.end());
        component.leaves.erase(std::unique(component.leaves.begin(),
                                           component.leaves.end()),
                               component.leaves.end());
        component.repr = component.leaves.front();
        for (Leaf leaf : component.leaves) {
            if (leaf < 1 || leaf > n || forest.leaf_comp[leaf] != -1) {
                return false;
            }
            forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        }
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    for (Leaf leaf = 1; leaf <= n; ++leaf) {
        if (forest.leaf_comp[leaf] == -1) return false;
    }
    return rebuild_forest_data(forest, t1, t2) &&
           validate_forest(forest, t1, t2).ok;
}

std::vector<CompId> largest_splittable_components(const Forest& forest, int limit) {
    std::vector<CompId> ids;
    for (CompId id : active_components(forest)) {
        if (forest.comps[id].leaves.size() > 1) ids.push_back(id);
    }
    std::sort(ids.begin(), ids.end(), [&](CompId a, CompId b) {
        size_t sa = forest.comps[a].leaves.size();
        size_t sb = forest.comps[b].leaves.size();
        if (sa != sb) return sa > sb;
        return a < b;
    });
    if (static_cast<int>(ids.size()) > limit) ids.resize(limit);
    return ids;
}

Forest split_selected_to_singletons(Forest forest, std::vector<CompId> ids,
                                    const Tree& t1, const Tree& t2) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    for (CompId id : ids) {
        if (id < 0 || id >= static_cast<CompId>(forest.comps.size())) continue;
        Component& old = forest.comps[id];
        if (!old.active || old.leaves.size() <= 1) continue;
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

bool run_standalone_whole_probe(const Tree& t1, const Tree& t2,
                                const Forest& cluster_forest,
                                Forest& out, bool debug) {
    if (env_int("RS_WHOLE_PROBE", 1) == 0) return false;
    const int min_leaves = env_int("RS_WHOLE_PROBE_MIN_LEAVES", 0);
    const int max_leaves = env_int("RS_WHOLE_PROBE_MAX_LEAVES", 0);
    const int max_cluster_components =
        env_int("RS_WHOLE_PROBE_MAX_CLUSTER_COMPONENTS", 350);
    if (min_leaves > 0 && t1.n < min_leaves) return false;
    if (max_leaves > 0 && t1.n > max_leaves) return false;
    if (max_cluster_components > 0 &&
        cluster_forest.active_count > max_cluster_components) return false;

    const auto start = std::chrono::steady_clock::now();
    Forest best = cluster_forest;
    if (debug) {
        std::cerr << "rs_whole_probe_start leaves=" << t1.n
                  << " cluster_components=" << cluster_forest.active_count
                  << '\n';
    }

    Forest singleton = make_singleton_forest(t1.n);
    if (!rebuild_forest_data(singleton, t1, t2)) return false;

    const int greedy_mode = env_int("RS_WHOLE_GREEDY_MODE", 1);
    Forest whole_seed =
        run_h45_greedy_merge(std::move(singleton), t1, t2,
                             greedy_mode, false, false);
    const int greedy_upper_bound = std::max(0, whole_seed.active_count - 1);
    consider_probe_candidate(whole_seed, t1, t2, best,
                             "singleton", debug, start);
    Forest restart_base = best;
    bool have_restart_base = false;

    const int singleton_gap = whole_seed.active_count - cluster_forest.active_count;
    const int min_singleton_gap = env_int(
        "RS_WHOLE_MIN_SINGLETON_GAP",
        std::max(8, cluster_forest.active_count / 20));
    const bool try_whidden = singleton_gap >= min_singleton_gap;
    if (debug && !try_whidden) {
        std::cerr << "rs_whole_probe_stage label=whidden_skipped"
                  << " singleton_gap=" << singleton_gap
                  << " min_gap=" << min_singleton_gap
                  << " seconds=" << seconds_since(start) << '\n';
    }

    std::vector<std::vector<int>> whidden_components;
    if (try_whidden && whidden_seed_components(
            tree_to_newick(t1), tree_to_newick(t2), t1.n,
            greedy_upper_bound,
            env_int("RS_WHOLE_WHIDDEN_TUNE", 30),
            env_int("RS_WHOLE_WHIDDEN_SPLIT", 25),
            env_int("RS_WHOLE_WHIDDEN_SEED", 7),
            whidden_components)) {
        Forest whidden_forest;
        if (forest_from_leaf_components(whidden_components, t1.n, t1, t2,
                                        whidden_forest)) {
            consider_probe_candidate(whidden_forest, t1, t2, best,
                                     "whidden", debug, start);
            Forest whidden_repaired =
                run_h45_greedy_merge(whidden_forest, t1, t2,
                                     env_int("RS_WHOLE_REPAIR_MODE", 1),
                                     true, true);
            consider_probe_candidate(whidden_repaired, t1, t2, best,
                                     "whidden_repair", debug, start);
            restart_base = std::move(whidden_repaired);
            have_restart_base = true;
        } else if (debug) {
            std::cerr << "rs_whole_probe_stage label=whidden_convert_failed"
                      << " seconds=" << seconds_since(start) << '\n';
        }
    } else if (debug && try_whidden) {
        std::cerr << "rs_whole_probe_stage label=whidden_failed"
                  << " seconds=" << seconds_since(start) << '\n';
    }

    const int repair_mode = env_int("RS_WHOLE_REPAIR_MODE", 1);
    Forest repaired =
        run_h45_greedy_merge(best, t1, t2, repair_mode, true, true);
    consider_probe_candidate(repaired, t1, t2, best,
                             "repair_best", debug, start);

    const int restart_budget = env_int("RS_WHOLE_RESTART_BUDGET", 6144);
    if (restart_budget > 0 && have_restart_base) {
        const int restart_mode = env_int("RS_WHOLE_RESTART_MODE", 1);
        const int restart_seed = env_int("RS_WHOLE_RESTART_SEED_OFFSET", 0);
        const int extra_splits = env_int("RS_WHOLE_EXTRA_SPLITS", 0);
        const int target = env_int(
            "RS_WHOLE_TARGET", std::max(1, cluster_forest.active_count - 2));
        const int seconds = env_int("RS_WHOLE_RESTART_SECONDS", 20);
        Forest restarted = run_h45_split_restarts(
            have_restart_base ? restart_base : best, t1, t2,
            restart_mode, restart_budget, restart_seed,
            extra_splits, false, target, static_cast<double>(seconds));
        consider_probe_candidate(restarted, t1, t2, best,
                                 "restart_best", debug, start);
    }

    const int perturb_variants = env_int("RS_WHOLE_PERTURB_VARIANTS", 0);
    const int perturb_top = env_int("RS_WHOLE_PERTURB_TOP", 16);
    const int perturb_margin = env_int("RS_WHOLE_PERTURB_MARGIN", 8);
    const int perturb_restart_budget =
        env_int("RS_WHOLE_PERTURB_RESTART_BUDGET", 2048);
    const int perturb_seconds = env_int("RS_WHOLE_PERTURB_SECONDS", 4);
    std::vector<CompId> large =
        largest_splittable_components(cluster_forest, perturb_top);
    for (int variant = 0;
         variant < perturb_variants && !large.empty() &&
         best.active_count > env_int(
             "RS_WHOLE_TARGET", std::max(1, cluster_forest.active_count - 2));
         ++variant) {
        std::vector<CompId> split_ids;
        split_ids.push_back(large[(variant * 5) % large.size()]);
        if (large.size() > 1 && (variant & 1)) {
            split_ids.push_back(large[(variant * 7 + 1) % large.size()]);
        }
        Forest seed = split_selected_to_singletons(cluster_forest, split_ids, t1, t2);
        const int mode = (variant % 3 == 0) ? 1 : ((variant % 3 == 1) ? 2 : 14);
        seed = run_h45_greedy_merge(std::move(seed), t1, t2, mode, true, true);
        seed = exact_output_repair(std::move(seed), t1, t2);
        ValidationResult seed_check = validate_forest(seed, t1, t2);
        if (debug) {
            std::cerr << "rs_whole_probe_perturb variant=" << variant
                      << " split_count=" << split_ids.size()
                      << " mode=" << mode
                      << " components=" << seed.active_count
                      << " valid=" << (seed_check.ok ? 1 : 0)
                      << " seconds=" << seconds_since(start) << '\n';
        }
        if (!seed_check.ok ||
            seed.active_count > cluster_forest.active_count + perturb_margin) {
            continue;
        }
        Forest restarted = run_h45_split_restarts(
            seed, t1, t2, env_int("RS_WHOLE_RESTART_MODE", 1),
            perturb_restart_budget, variant * 97,
            env_int("RS_WHOLE_EXTRA_SPLITS", 0), false,
            env_int("RS_WHOLE_TARGET",
                    std::max(1, cluster_forest.active_count - 2)),
            static_cast<double>(perturb_seconds));
        consider_probe_candidate(restarted, t1, t2, best,
                                 "perturb_restart", debug, start);
    }

    if (best.active_count < cluster_forest.active_count) {
        out = std::move(best);
        if (debug) {
            std::cerr << "rs_whole_probe_selected components="
                      << out.active_count
                      << " seconds=" << seconds_since(start) << '\n';
        }
        return true;
    }
    if (debug) {
        std::cerr << "rs_whole_probe_no_improve components="
                  << best.active_count
                  << " seconds=" << seconds_since(start) << '\n';
    }
    return false;
}

bool parse_cluster_atom(const std::string& atom, int& cluster_id) {
    if (atom.size() < 2 || atom[0] != 'C') return false;
    int value = 0;
    for (size_t i = 1; i < atom.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(atom[i]))) return false;
        value = value * 10 + (atom[i] - '0');
    }
    cluster_id = value;
    return value > 0;
}

void expand_atom_variant(const std::string& atom,
                         const std::vector<LocalSolution>& bridge_solved,
                         const std::vector<unsigned char>& nonbridge,
                         std::vector<unsigned char>& visiting,
                         std::vector<std::vector<int>>& boundary_cache,
                         std::vector<int>& leaves) {
    int cluster_id = 0;
    if (!parse_cluster_atom(atom, cluster_id)) {
        leaves.push_back(std::stoi(atom));
        return;
    }
    if (cluster_id <= 0 ||
        cluster_id >= static_cast<int>(bridge_solved.size())) {
        throw std::runtime_error("cluster atom out of range");
    }
    if (cluster_id < static_cast<int>(nonbridge.size()) &&
        nonbridge[cluster_id]) {
        return;
    }
    if (boundary_cache[cluster_id].empty() &&
        !bridge_solved[cluster_id].components.empty() &&
        !bridge_solved[cluster_id].components[0].empty()) {
        if (visiting[cluster_id]) throw std::runtime_error("cluster expansion cycle");
        visiting[cluster_id] = 1;
        std::vector<int> expanded;
        for (const std::string& nested : bridge_solved[cluster_id].components[0]) {
            expand_atom_variant(nested, bridge_solved, nonbridge, visiting,
                                boundary_cache, expanded);
        }
        std::sort(expanded.begin(), expanded.end());
        expanded.erase(std::unique(expanded.begin(), expanded.end()), expanded.end());
        boundary_cache[cluster_id] = std::move(expanded);
        visiting[cluster_id] = 0;
    }
    leaves.insert(leaves.end(), boundary_cache[cluster_id].begin(),
                  boundary_cache[cluster_id].end());
}

std::vector<int> expand_component_variant(
    const std::vector<std::string>& atoms,
    const std::vector<LocalSolution>& bridge_solved,
    const std::vector<unsigned char>& nonbridge,
    std::vector<unsigned char>& visiting,
    std::vector<std::vector<int>>& boundary_cache) {
    std::vector<int> leaves;
    for (const std::string& atom : atoms) {
        expand_atom_variant(atom, bridge_solved, nonbridge, visiting,
                            boundary_cache, leaves);
    }
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    return leaves;
}

Forest forest_from_components(const std::vector<std::vector<int>>& components,
                              int n, const Tree& t1, const Tree& t2) {
    Forest forest;
    forest.leaf_comp.assign(n + 1, -1);
    for (const std::vector<int>& leaves : components) {
        if (leaves.empty()) continue;
        Component component;
        component.leaves = leaves;
        for (Leaf leaf : leaves) {
            if (leaf < 1 || leaf > n) throw std::runtime_error("global leaf out of range");
            if (forest.leaf_comp[leaf] != -1) {
                throw std::runtime_error("duplicate global leaf");
            }
            forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        }
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    if (!rebuild_forest_data(forest, t1, t2)) {
        throw std::runtime_error("expanded forest data invalid");
    }
    return forest;
}

Forest add_missing_singletons(Forest forest, int n, const Tree& t1,
                              const Tree& t2) {
    std::vector<unsigned char> seen(n + 1, 0);
    for (const Component& component : forest.comps) {
        if (!component.active) continue;
        for (Leaf leaf : component.leaves) {
            if (leaf >= 1 && leaf <= n) seen[leaf] = 1;
        }
    }
    if (forest.leaf_comp.size() < static_cast<size_t>(n + 1)) {
        forest.leaf_comp.resize(n + 1, -1);
    }
    for (Leaf leaf = 1; leaf <= n; ++leaf) {
        if (seen[leaf]) continue;
        Component component;
        component.leaves.push_back(leaf);
        component.repr = leaf;
        forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    rebuild_forest_data(forest, t1, t2);
    return forest;
}

}  // namespace

bool run_cluster_solver(const Tree& t1, const Tree& t2, Forest& out) {
    try {
        ++g_cluster_solver_depth;
        struct DepthGuard {
            ~DepthGuard() { --g_cluster_solver_depth; }
        } depth_guard;

        std::vector<ClusterCandidate> candidates = find_cluster_candidates(t1, t2);
        ClusterPlan plan = build_cluster_plan(t1, candidates);
        if (plan.clusters.empty()) return false;

        const bool debug = std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr;
        if (debug) {
            std::cerr << "rs_cluster_solver_start proper_clusters="
                      << plan.clusters.size()
                      << " residual_leaves=" << plan.residual_leaves << '\n';
        }

        auto publish_valid = [&](const Forest& forest, const char* label) {
            ValidationResult check = validate_forest(forest, t1, t2);
            if (!check.ok) return false;
            if (g_cluster_solver_depth == 1) {
                publish_best_output(forest_to_output(forest, t1));
            }
            if (debug) {
                std::cerr << "rs_cluster_publish label=" << label
                          << " components=" << forest.active_count << '\n';
            }
            return true;
        };

        std::vector<LocalSolution> solved(plan.residual_id + 1);
        const bool bottomup_cleanup =
            env_int("RS_BOTTOMUP_CLEANUP", 0) != 0;
        std::vector<unsigned char> bottomup_sealed(plan.residual_id + 1, 0);
        if (bottomup_cleanup) {
            std::vector<std::vector<int>> children(plan.residual_id + 1);
            std::vector<int> top_level;
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.parent > 0 &&
                    cluster.parent < static_cast<int>(children.size())) {
                    children[cluster.parent].push_back(cluster.id);
                } else if (cluster.parent == 0) {
                    top_level.push_back(cluster.id);
                }
            }
            for (std::vector<int>& group : children) {
                std::sort(group.begin(), group.end(), [&](int a, int b) {
                    int sa = plan.clusters[a - 1].original_leaves;
                    int sb = plan.clusters[b - 1].original_leaves;
                    if (sa != sb) return sa < sb;
                    return a < b;
                });
            }
            std::sort(top_level.begin(), top_level.end(), [&](int a, int b) {
                int sa = plan.clusters[a - 1].original_leaves;
                int sb = plan.clusters[b - 1].original_leaves;
                if (sa != sb) return sa < sb;
                return a < b;
            });
            std::vector<int> bottomup_order;
            std::vector<unsigned char> visited(plan.residual_id + 1, 0);
            auto dfs = [&](auto&& self, int cluster_id) -> void {
                if (cluster_id <= 0 ||
                    cluster_id >= static_cast<int>(visited.size()) ||
                    visited[cluster_id]) {
                    return;
                }
                visited[cluster_id] = 1;
                for (int child : children[cluster_id]) self(self, child);
                bottomup_order.push_back(cluster_id);
            };
            for (int cluster_id : top_level) dfs(dfs, cluster_id);
            for (const PlannedCluster& cluster : plan.clusters) {
                dfs(dfs, cluster.id);
            }

            if (debug) {
                std::cerr << "rs_bottomup_cleanup_start order_count="
                          << bottomup_order.size() << '\n';
            }
            for (int cluster_id : bottomup_order) {
                solved[cluster_id] = solve_local_route_with_suppressed(
                    t1, t2, plan, cluster_id, true, cluster_id,
                    bottomup_sealed);
                bool sealed =
                    !solved[cluster_id].components.empty() &&
                    solved[cluster_id].components[0].empty();
                bottomup_sealed[cluster_id] = sealed ? 1 : 0;
                if (debug) {
                    int suppressed_direct = 0;
                    for (int child : children[cluster_id]) {
                        if (bottomup_sealed[child]) ++suppressed_direct;
                    }
                    std::cerr << "rs_cluster_solved id=" << cluster_id
                              << " components=" << solved[cluster_id].score
                              << " boundary_atoms="
                              << (solved[cluster_id].components.empty()
                                      ? 0
                                      : solved[cluster_id].components[0].size())
                              << " bottomup=1 sealed=" << (sealed ? 1 : 0)
                              << " suppressed_children="
                              << suppressed_direct << '\n';
                }
            }
            solved[plan.residual_id] = solve_local_route_with_suppressed(
                t1, t2, plan, 0, false, plan.residual_id, bottomup_sealed);
            if (debug) {
                int sealed_top = 0;
                if (std::getenv("RS_DEBUG_BOTTOMUP_MARKERS") != nullptr) {
                    std::cerr << "rs_bottomup_top_sealed ids=";
                }
                for (int cluster_id : top_level) {
                    if (bottomup_sealed[cluster_id]) {
                        ++sealed_top;
                        if (std::getenv("RS_DEBUG_BOTTOMUP_MARKERS") != nullptr) {
                            std::cerr << cluster_id << ',';
                        }
                    }
                }
                if (std::getenv("RS_DEBUG_BOTTOMUP_MARKERS") != nullptr) {
                    std::cerr << '\n';
                }
                std::cerr << "rs_cluster_solved id=" << plan.residual_id
                          << " components=" << solved[plan.residual_id].score
                          << " bottomup=1 sealed_top=" << sealed_top << '\n';
            }
        } else {
            solved[plan.residual_id] =
                solve_local_route(t1, t2, plan, 0, false, plan.residual_id);
            if (debug) {
                std::cerr << "rs_cluster_solved id=" << plan.residual_id
                          << " components=" << solved[plan.residual_id].score
                          << '\n';
            }

            for (int cluster_id : plan.solve_order) {
                solved[cluster_id] =
                    solve_local_route(t1, t2, plan, cluster_id, true, cluster_id);
                if (debug) {
                    std::cerr << "rs_cluster_solved id=" << cluster_id
                              << " components=" << solved[cluster_id].score
                              << " boundary_atoms="
                              << (solved[cluster_id].components.empty()
                                      ? 0
                                      : solved[cluster_id].components[0].size())
                              << '\n';
                }
            }
        }

        auto build_expanded_forest =
            [&](const std::vector<LocalSolution>& local_solved,
                const std::vector<LocalSolution>& plain_solved,
                const std::vector<unsigned char>& nonbridge,
                const std::vector<unsigned char>& suppressed) {
            std::vector<std::vector<int>> global_components;
            std::vector<std::vector<int>> boundary_cache(local_solved.size());
            std::vector<unsigned char> visiting(local_solved.size(), 0);
            int empty_expansions = 0;
            int abstract_count = 0;
            for (size_t index = 0;
                 index < local_solved[plan.residual_id].components.size(); ++index) {
                const std::vector<std::string>& component =
                    local_solved[plan.residual_id].components[index];
                std::vector<int> expanded_component =
                    expand_component_variant(component, local_solved, nonbridge,
                                             visiting, boundary_cache);
                if (expanded_component.empty()) {
                    ++empty_expansions;
                    if (debug) {
                        std::cerr << "rs_cluster_empty_expansion index="
                                  << abstract_count << " atoms=";
                        for (const std::string& atom : component) {
                            std::cerr << atom << ',';
                        }
                        std::cerr << '\n';
                    }
                }
                global_components.push_back(std::move(expanded_component));
                ++abstract_count;
            }
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.id < static_cast<int>(suppressed.size()) &&
                    suppressed[cluster.id]) {
                    continue;
                }
                const bool use_plain =
                    cluster.id < static_cast<int>(nonbridge.size()) &&
                    nonbridge[cluster.id];
                const LocalSolution& solution =
                    use_plain ? plain_solved[cluster.id] : local_solved[cluster.id];
                const size_t begin = use_plain ? 0 : 1;
                for (size_t i = begin; i < solution.components.size(); ++i) {
                    const std::vector<std::string>& component =
                        solution.components[i];
                    std::vector<int> expanded_component =
                        expand_component_variant(component, local_solved, nonbridge,
                                                 visiting, boundary_cache);
                    if (expanded_component.empty()) {
                        ++empty_expansions;
                        if (debug) {
                            std::cerr << "rs_cluster_empty_expansion index="
                                      << abstract_count << " atoms=";
                            for (const std::string& atom : component) {
                                std::cerr << atom << ',';
                            }
                            std::cerr << '\n';
                        }
                    }
                    global_components.push_back(std::move(expanded_component));
                    ++abstract_count;
                }
            }
            if (debug) {
                std::cerr << "rs_cluster_empty_expansions count="
                          << empty_expansions << '\n';
                std::cerr << "rs_cluster_abstract_components count="
                          << abstract_count << '\n';
            }
            return forest_from_components(global_components, t1.n, t1, t2);
        };

        std::vector<LocalSolution> plain_solved(plan.residual_id + 1);
        std::vector<unsigned char> marker_nonbridge(solved.size(), 0);
        std::vector<unsigned char> suppressed_clusters(solved.size(), 0);
        Forest expanded = build_expanded_forest(
            solved, plain_solved, marker_nonbridge, suppressed_clusters);
        if (debug) {
            std::cerr << "rs_cluster_expanded components="
                      << expanded.active_count << '\n';
        }
        publish_valid(expanded, "expanded");

        const bool run_marker_rule =
            !bottomup_cleanup &&
            env_int("RS_CLUSTER_MARKER_RULE", 1) != 0 &&
            expanded.active_count <=
                env_int("RS_CLUSTER_MARKER_RULE_MAX_COMPONENTS", 220);
        if (run_marker_rule) {
            for (int cluster_id : plan.solve_order) {
                plain_solved[cluster_id] =
                    solve_local_route(t1, t2, plan, cluster_id, false, cluster_id);
                if (debug) {
                    std::cerr << "rs_cluster_plain_solved id=" << cluster_id
                              << " components=" << plain_solved[cluster_id].score
                              << '\n';
                }
            }

            auto consider_marker_flags =
                [&](const std::vector<unsigned char>& flags,
                    const char* label) {
                    try {
                        Forest candidate = build_expanded_forest(
                            solved, plain_solved, flags, suppressed_clusters);
                        candidate =
                            exact_output_repair(std::move(candidate), t1, t2);
                        ValidationResult check =
                            validate_forest(candidate, t1, t2);
                        if (debug) {
                            std::cerr << "rs_cluster_marker_variant label="
                                      << label
                                      << " components=" << candidate.active_count
                                      << " valid=" << (check.ok ? 1 : 0) << '\n';
                        }
                        if (check.ok &&
                            candidate.active_count < expanded.active_count) {
                            expanded = std::move(candidate);
                            marker_nonbridge = flags;
                            publish_valid(expanded, "marker");
                            if (debug) {
                                std::cerr
                                    << "rs_cluster_marker_improve components="
                                    << expanded.active_count << '\n';
                            }
                            return true;
                        }
                    } catch (const std::exception& e) {
                        if (debug) {
                            std::cerr << "rs_cluster_marker_variant_failed label="
                                      << label << " reason=" << e.what() << '\n';
                        }
                    }
                    return false;
                };

            if (run_marker_rule) {
                std::vector<unsigned char> all_nonbridge = marker_nonbridge;
                for (const PlannedCluster& cluster : plan.clusters) {
                    all_nonbridge[cluster.id] = 1;
                }
                consider_marker_flags(all_nonbridge, "all_nonbridge");

                bool changed = true;
                while (changed) {
                    changed = false;
                    for (const PlannedCluster& cluster : plan.clusters) {
                        std::vector<unsigned char> trial = marker_nonbridge;
                        trial[cluster.id] = trial[cluster.id] ? 0 : 1;
                        if (consider_marker_flags(trial, "flip")) {
                            changed = true;
                        }
                    }
                }
            }
        }

        if (!bottomup_cleanup &&
            env_int("RS_CLUSTER_OPEN_CHILD_VARIANTS", 1) != 0) {
            std::vector<std::vector<int>> children(plan.residual_id + 1);
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.parent > 0 &&
                    cluster.parent < static_cast<int>(children.size())) {
                    children[cluster.parent].push_back(cluster.id);
                }
            }

            std::vector<int> parent_order;
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.id < static_cast<int>(children.size()) &&
                    !children[cluster.id].empty()) {
                    parent_order.push_back(cluster.id);
                }
            }
            std::sort(parent_order.begin(), parent_order.end(), [&](int a, int b) {
                int sa = plan.clusters[a - 1].original_leaves;
                int sb = plan.clusters[b - 1].original_leaves;
                if (sa != sb) return sa > sb;
                return a < b;
            });

            const int max_variants =
                env_int("RS_CLUSTER_OPEN_CHILD_VARIANT_LIMIT", 12);
            const int max_atoms =
                env_int("RS_CLUSTER_OPEN_CHILD_MAX_ATOMS", 120);
            int tried = 0;
            for (int parent_id : parent_order) {
                if (tried >= max_variants) break;
                const std::vector<int>& open_children = children[parent_id];
                if (open_children.empty()) continue;

                LocalProblem open_problem =
                    build_local_problem_with_open_children(
                        t1, t2, plan, parent_id, true, open_children);
                if (open_problem.solve_n > max_atoms) {
                    if (debug) {
                        std::cerr << "rs_cluster_open_skip parent=" << parent_id
                                  << " atoms=" << open_problem.solve_n << '\n';
                    }
                    continue;
                }
                ++tried;
                try {
                    LocalSolution open_solution = solve_prebuilt_local_route(
                        open_problem, parent_id, true, false, parent_id + 7000);
                    std::vector<LocalSolution> trial_solved = solved;
                    trial_solved[parent_id] = std::move(open_solution);
                    std::vector<unsigned char> trial_suppressed =
                        suppressed_clusters;
                    for (int child_id : open_children) {
                        if (child_id >= 0 &&
                            child_id < static_cast<int>(trial_suppressed.size())) {
                            trial_suppressed[child_id] = 1;
                        }
                    }

                    Forest candidate = build_expanded_forest(
                        trial_solved, plain_solved, marker_nonbridge,
                        trial_suppressed);
                    candidate =
                        exact_output_repair(std::move(candidate), t1, t2);
                    ValidationResult check = validate_forest(candidate, t1, t2);
                    if (debug) {
                        std::cerr << "rs_cluster_open_variant parent="
                                  << parent_id
                                  << " children=" << open_children.size()
                                  << " atoms=" << open_problem.solve_n
                                  << " components=" << candidate.active_count
                                  << " valid=" << (check.ok ? 1 : 0)
                                  << " reason=" << check.message << '\n';
                    }
                    if (check.ok && candidate.active_count < expanded.active_count) {
                        expanded = std::move(candidate);
                        suppressed_clusters = std::move(trial_suppressed);
                        solved = std::move(trial_solved);
                        publish_valid(expanded, "open_children");
                        if (debug) {
                            std::cerr
                                << "rs_cluster_open_improve parent="
                                << parent_id
                                << " components=" << expanded.active_count
                                << '\n';
                        }
                    }
                } catch (const std::exception& e) {
                    if (debug) {
                        std::cerr << "rs_cluster_open_failed parent="
                                  << parent_id << " reason=" << e.what()
                                  << '\n';
                    }
                }
            }
        }

        if (!bottomup_cleanup &&
            env_int("RS_CLUSTER_OPEN_RESIDUAL_VARIANTS", 0) != 0) {
            std::vector<std::vector<int>> residual_children(plan.residual_id + 1);
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.parent >= 0 &&
                    cluster.parent < static_cast<int>(residual_children.size())) {
                    residual_children[cluster.parent].push_back(cluster.id);
                }
            }
            std::vector<int> top_level;
            for (const PlannedCluster& cluster : plan.clusters) {
                if (cluster.parent == 0) top_level.push_back(cluster.id);
            }
            std::sort(top_level.begin(), top_level.end(), [&](int a, int b) {
                int sa = plan.clusters[a - 1].original_leaves;
                int sb = plan.clusters[b - 1].original_leaves;
                if (sa != sb) return sa > sb;
                return a < b;
            });
            const int top_limit =
                env_int("RS_CLUSTER_OPEN_RESIDUAL_TOP_LIMIT", 8);
            if (static_cast<int>(top_level.size()) > top_limit) {
                top_level.resize(top_limit);
            }

            const int max_variants =
                env_int("RS_CLUSTER_OPEN_RESIDUAL_VARIANT_LIMIT", 24);
            const int max_atoms =
                env_int("RS_CLUSTER_OPEN_RESIDUAL_MAX_ATOMS", 220);
            const int max_open =
                env_int("RS_CLUSTER_OPEN_RESIDUAL_MAX_OPEN", 3);
            int tried = 0;

            auto try_residual_open = [&](const std::vector<int>& opened) {
                if (opened.empty() || tried >= max_variants) return;
                LocalProblem open_problem =
                    build_local_problem_with_open_children(
                        t1, t2, plan, 0, false, opened);
                if (open_problem.solve_n > max_atoms) {
                    if (debug) {
                        std::cerr << "rs_cluster_residual_open_skip opened=";
                        for (int id : opened) std::cerr << id << ',';
                        std::cerr << " atoms=" << open_problem.solve_n << '\n';
                    }
                    return;
                }
                ++tried;
                try {
                    LocalSolution open_solution = solve_prebuilt_local_route(
                        open_problem, plan.residual_id, false, true,
                        plan.residual_id + 9000 + tried);
                    std::vector<LocalSolution> trial_solved = solved;
                    trial_solved[plan.residual_id] = std::move(open_solution);
                    std::vector<unsigned char> trial_suppressed =
                        suppressed_clusters;
                    for (int cluster_id : opened) {
                        bool marker_still_present = false;
                        const std::string marker =
                            "C" + std::to_string(cluster_id);
                        for (const std::string& atom : open_problem.atom_names) {
                            if (atom == marker) {
                                marker_still_present = true;
                                break;
                            }
                        }
                        if (marker_still_present) continue;
                        if (cluster_id >= 0 &&
                            cluster_id < static_cast<int>(trial_suppressed.size())) {
                            trial_suppressed[cluster_id] = 1;
                        }
                    }
                    Forest candidate = build_expanded_forest(
                        trial_solved, plain_solved, marker_nonbridge,
                        trial_suppressed);
                    candidate =
                        exact_output_repair(std::move(candidate), t1, t2);
                    ValidationResult check = validate_forest(candidate, t1, t2);
                    if (!check.ok && check.message == "missing leaf") {
                        candidate = add_missing_singletons(
                            std::move(candidate), t1.n, t1, t2);
                        candidate = exact_output_repair(
                            std::move(candidate), t1, t2);
                        check = validate_forest(candidate, t1, t2);
                    }
                    if (debug) {
                        std::cerr << "rs_cluster_residual_open_variant opened=";
                        for (int id : opened) std::cerr << id << ',';
                        std::cerr << " atoms=" << open_problem.solve_n
                                  << " components=" << candidate.active_count
                                  << " valid=" << (check.ok ? 1 : 0)
                                  << " reason=" << check.message << '\n';
                    }
                    if (check.ok && candidate.active_count < expanded.active_count) {
                        expanded = std::move(candidate);
                        solved = std::move(trial_solved);
                        suppressed_clusters = std::move(trial_suppressed);
                        publish_valid(expanded, "residual_open");
                        if (debug) {
                            std::cerr
                                << "rs_cluster_residual_open_improve components="
                                << expanded.active_count << '\n';
                        }
                    }
                } catch (const std::exception& e) {
                    if (debug) {
                        std::cerr << "rs_cluster_residual_open_failed opened=";
                        for (int id : opened) std::cerr << id << ',';
                        std::cerr << " reason=" << e.what() << '\n';
                    }
                }
            };

            auto sorted_unique = [](std::vector<int> ids) {
                std::sort(ids.begin(), ids.end());
                ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                return ids;
            };

            auto open_with_family = [&](const std::vector<int>& roots) {
                std::vector<int> opened = roots;
                for (int root : roots) {
                    if (root <= 0 ||
                        root >= static_cast<int>(residual_children.size())) {
                        continue;
                    }
                    int largest_child = 0;
                    int largest_size = -1;
                    for (int child : residual_children[root]) {
                        opened.push_back(child);
                        int size = plan.clusters[child - 1].original_leaves;
                        if (size > largest_size) {
                            largest_size = size;
                            largest_child = child;
                        }
                    }
                    int cursor = largest_child;
                    while (cursor > 0 &&
                           cursor < static_cast<int>(residual_children.size()) &&
                           !residual_children[cursor].empty()) {
                        int next = 0;
                        int next_size = -1;
                        for (int child : residual_children[cursor]) {
                            int size = plan.clusters[child - 1].original_leaves;
                            if (size > next_size) {
                                next_size = size;
                                next = child;
                            }
                        }
                        if (next <= 0) break;
                        opened.push_back(next);
                        cursor = next;
                    }
                }
                return sorted_unique(std::move(opened));
            };

            for (size_t i = 0; i < top_level.size() && tried < max_variants; ++i) {
                try_residual_open({top_level[i]});
            }
            if (max_open >= 2) {
                for (size_t i = 0; i < top_level.size() && tried < max_variants; ++i) {
                    for (size_t j = i + 1;
                         j < top_level.size() && tried < max_variants; ++j) {
                        try_residual_open({top_level[i], top_level[j]});
                    }
                }
            }
            if (max_open >= 3) {
                for (size_t i = 0; i < top_level.size() && tried < max_variants; ++i) {
                    for (size_t j = i + 1;
                         j < top_level.size() && tried < max_variants; ++j) {
                        for (size_t k = j + 1;
                             k < top_level.size() && tried < max_variants; ++k) {
                            try_residual_open(
                                {top_level[i], top_level[j], top_level[k]});
                        }
                    }
                }
            }
            if (env_int("RS_CLUSTER_OPEN_RESIDUAL_FAMILY", 1) != 0) {
                for (size_t i = 0; i < top_level.size() && tried < max_variants; ++i) {
                    try_residual_open(open_with_family({top_level[i]}));
                }
                if (max_open >= 2) {
                    for (size_t i = 0; i < top_level.size() &&
                                       tried < max_variants; ++i) {
                        for (size_t j = i + 1;
                             j < top_level.size() && tried < max_variants; ++j) {
                            try_residual_open(
                                open_with_family({top_level[i], top_level[j]}));
                        }
                    }
                }
                if (max_open >= 3) {
                    for (size_t i = 0; i < top_level.size() &&
                                       tried < max_variants; ++i) {
                        for (size_t j = i + 1;
                             j < top_level.size() && tried < max_variants; ++j) {
                            for (size_t k = j + 1;
                                 k < top_level.size() && tried < max_variants; ++k) {
                                try_residual_open(open_with_family(
                                    {top_level[i], top_level[j], top_level[k]}));
                            }
                        }
                    }
                }
            }
        }

        Forest final_candidate =
            run_h45_greedy_merge(expanded, t1, t2, 14, false, true);
        final_candidate = exact_output_repair(std::move(final_candidate), t1, t2);
        if (validate_forest(final_candidate, t1, t2).ok &&
            final_candidate.active_count < expanded.active_count) {
            expanded = std::move(final_candidate);
            publish_valid(expanded, "final_greedy");
            if (debug) {
                std::cerr << "rs_cluster_final_greedy_improve components="
                          << expanded.active_count << '\n';
            }
        }

        bool whole_probe_selected = false;
        Forest whole_probe;
        if (run_standalone_whole_probe(t1, t2, expanded, whole_probe, debug)) {
            expanded = std::move(whole_probe);
            whole_probe_selected = true;
            publish_valid(expanded, "whole_probe");
        }

        int final_restart_budget = env_int("RS_FINAL_RESTART_BUDGET", 200000);
        if (final_restart_budget > 0 && !whole_probe_selected) {
            const int final_score_mode = env_int("RS_FINAL_SCORE_MODE", 14);
            const int final_seed_offset = env_int("RS_FINAL_SEED_OFFSET", 6700);
            const int final_extra_splits = env_int("RS_FINAL_EXTRA_SPLITS", 2);
            const int final_seconds = env_int("RS_FINAL_SECONDS", 60);
            const int final_target =
                env_int("RS_FINAL_TARGET", std::max(1, expanded.active_count - 1));
            Forest restart_candidate = run_h45_split_restarts(
                expanded, t1, t2, final_score_mode, final_restart_budget,
                final_seed_offset, final_extra_splits, true, final_target,
                static_cast<double>(final_seconds));
            restart_candidate = exact_output_repair(std::move(restart_candidate), t1, t2);
            if (validate_forest(restart_candidate, t1, t2).ok &&
                restart_candidate.active_count < expanded.active_count) {
                expanded = std::move(restart_candidate);
                publish_valid(expanded, "final_restart");
                if (debug) {
                    std::cerr << "rs_cluster_final_restart_improve components="
                              << expanded.active_count << '\n';
                }
            }
        }
        if (env_int("RS_CLUSTER_BLOCKER_REPAIR", 0) != 0 &&
            expanded.active_count <= env_int("RS_CLUSTER_BLOCKER_MAX_COMPONENTS", 190)) {
            Forest blocker = run_h45_blocker_repair(
                expanded, t1, t2,
                env_int("RS_CLUSTER_BLOCKER_MAX_COMPONENTS", 190),
                env_int("RS_CLUSTER_BLOCKER_MAX_BLOCKERS", 3),
                env_int("RS_CLUSTER_BLOCKER_CANDIDATES", 12),
                env_int("RS_CLUSTER_BLOCKER_RESTART_BUDGET", 256),
                env_int("RS_CLUSTER_BLOCKER_TARGET",
                        std::max(1, expanded.active_count - 1)),
                static_cast<double>(env_int("RS_CLUSTER_BLOCKER_SECONDS", 8)));
            if (validate_forest(blocker, t1, t2).ok &&
                blocker.active_count < expanded.active_count) {
                expanded = std::move(blocker);
                publish_valid(expanded, "blocker");
                if (debug) {
                    std::cerr << "rs_cluster_blocker_improve components="
                              << expanded.active_count << '\n';
                }
            }
        }
        if (env_int("RS_CLUSTER_WINDOW_REPAIR", 0) != 0 &&
            expanded.active_count <= env_int("RS_CLUSTER_WINDOW_MAX_COMPONENTS", 220)) {
            Forest window = run_h45_window_repair(
                expanded, t1, t2,
                env_int("RS_CLUSTER_WINDOW_MAX_COMPONENTS", 220),
                env_int("RS_CLUSTER_WINDOW_CANDIDATES", 10),
                env_int("RS_CLUSTER_WINDOW_SIZE", 18),
                env_int("RS_CLUSTER_WINDOW_RESTART_BUDGET", 256),
                env_int("RS_CLUSTER_WINDOW_TARGET",
                        std::max(1, expanded.active_count - 1)),
                static_cast<double>(env_int("RS_CLUSTER_WINDOW_SECONDS", 10)));
            if (validate_forest(window, t1, t2).ok &&
                window.active_count < expanded.active_count) {
                expanded = std::move(window);
                publish_valid(expanded, "window");
                if (debug) {
                    std::cerr << "rs_cluster_window_improve components="
                              << expanded.active_count << '\n';
                }
            }
        }
        if (env_int("RS_CROSS_RECOMBINE", 0) != 0 &&
            expanded.active_count <= env_int("RS_CROSS_RECOMBINE_MAX_COMPONENTS", 260)) {
            Forest recombined = run_h45_cross_recombine_repair(
                expanded, t1, t2,
                env_int("RS_CROSS_RECOMBINE_MAX_COMPONENTS", 260),
                env_int("RS_CROSS_RECOMBINE_CANDIDATES", 12),
                env_int("RS_CROSS_RECOMBINE_WINDOW", 48),
                env_int("RS_CROSS_RECOMBINE_RESTART_BUDGET", 512),
                env_int("RS_CROSS_RECOMBINE_TARGET",
                        std::max(1, expanded.active_count - 1)),
                static_cast<double>(env_int("RS_CROSS_RECOMBINE_SECONDS", 20)));
            if (validate_forest(recombined, t1, t2).ok &&
                recombined.active_count < expanded.active_count) {
                expanded = std::move(recombined);
                publish_valid(expanded, "cross_recombine");
                if (debug) {
                    std::cerr << "rs_cluster_cross_recombine_improve components="
                              << expanded.active_count << '\n';
                }
            }
        }
        if (env_int("RS_CLUSTER_SPLIT_REPACK", 0) != 0) {
            Forest repacked = run_h45_targeted_split_repack(
                expanded, t1, t2,
                env_int("RS_CLUSTER_SPLIT_REPACK_MAX_COMPONENTS", 800),
                env_int("RS_CLUSTER_SPLIT_REPACK_CANDIDATES", 24),
                env_int("RS_CLUSTER_SPLIT_REPACK_GROUP", 2),
                env_int("RS_CLUSTER_SPLIT_REPACK_RESTART_BUDGET", 128),
                env_int("RS_CLUSTER_SPLIT_REPACK_PASSES", 2),
                env_int("RS_CLUSTER_SPLIT_REPACK_TARGET",
                        std::max(1, expanded.active_count - 1)),
                static_cast<double>(env_int("RS_CLUSTER_SPLIT_REPACK_SECONDS", 20)));
            if (validate_forest(repacked, t1, t2).ok &&
                repacked.active_count < expanded.active_count) {
                expanded = std::move(repacked);
                publish_valid(expanded, "split_repack");
                if (debug) {
                    std::cerr << "rs_cluster_split_repack_improve components="
                              << expanded.active_count << '\n';
                }
            }
        }
        ValidationResult check = validate_forest(expanded, t1, t2);
        if (!check.ok) {
            if (debug) {
                std::cerr << "rs_cluster_invalid reason=" << check.message << '\n';
            }
            return false;
        }
        if (debug) {
            std::cerr << "rs_cluster_final components=" << expanded.active_count << '\n';
        }
        out = std::move(expanded);
        return true;
    } catch (const std::exception& e) {
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
            std::cerr << "rs_cluster_exception " << e.what() << '\n';
        }
        return false;
    }
}
