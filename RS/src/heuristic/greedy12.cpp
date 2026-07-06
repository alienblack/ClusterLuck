#include "heuristic/greedy12.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <random>
#include <unordered_set>
#include <vector>

#include "tree/induced.h"
#include "util/signal_handler.h"

namespace {

constexpr int kLocalWindow = 64;
constexpr int kNeighborWalkLimit = 1000000;
constexpr int kRandomCandidateMultiplier = 2;
constexpr int kMinMergeChecks = 10000;
constexpr int kMergeCheckMultiplier = 80;

struct Greedy12Config {
    int local_window = kLocalWindow;
    int node_pair_cap = 0;
    int neighbor_pair_cap = 0;
    int random_multiplier = kRandomCandidateMultiplier;
    int merge_check_multiplier = kMergeCheckMultiplier;
    int seed = 1;
    int score_mode = 0;
    int score_noise = 0;
};

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

struct Component12 {
    bool active = true;
    std::vector<Leaf> leaves;
    Leaf repr = -1;
    NodeId root1 = -1;
    NodeId root2 = -1;
};

struct Forest12 {
    std::vector<Component12> comps;
    std::vector<CompId> leaf_comp;
    std::vector<CompId> owner1;
    std::vector<CompId> owner2;
    int active_count = 0;
};

struct Candidate12 {
    int score = 0;
    Leaf a = -1;
    Leaf b = -1;
};

struct CandidateGreater12 {
    bool operator()(const Candidate12& x, const Candidate12& y) const {
        return x.score > y.score;
    }
};

using Queue12 = std::priority_queue<Candidate12, std::vector<Candidate12>,
                                    CandidateGreater12>;

uint64_t pair_key(Leaf a, Leaf b) {
    Leaf lo = std::min(a, b);
    Leaf hi = std::max(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(lo)) << 32) |
           static_cast<uint32_t>(hi);
}

std::vector<int> congestion_prefix(const Tree& tree) {
    std::vector<int> subtree(tree.node_count(), 0);
    std::vector<NodeId> order(tree.node_count());
    for (NodeId v = 0; v < tree.node_count(); ++v) order[v] = v;
    std::sort(order.begin(), order.end(), [&](NodeId a, NodeId b) {
        return tree.depth[a] > tree.depth[b];
    });
    for (NodeId v : order) {
        if (tree.is_leaf(v)) {
            subtree[v] = 1;
        } else {
            int count = 0;
            if (tree.left[v] != -1) count += subtree[tree.left[v]];
            if (tree.right[v] != -1) count += subtree[tree.right[v]];
            subtree[v] = count;
        }
    }

    std::vector<int> prefix(tree.node_count(), 0);
    std::vector<NodeId> top_down(tree.node_count());
    for (NodeId v = 0; v < tree.node_count(); ++v) top_down[v] = v;
    std::sort(top_down.begin(), top_down.end(), [&](NodeId a, NodeId b) {
        return tree.depth[a] < tree.depth[b];
    });
    for (NodeId v : top_down) {
        if (v == tree.root) continue;
        long long s = subtree[v];
        long long demand = std::max(1LL, s * (tree.n - s));
        int edge_price = 1;
        while ((1LL << std::min(edge_price, 60)) <= demand &&
               edge_price < 60) {
            ++edge_price;
        }
        prefix[v] = prefix[tree.parent[v]] + edge_price;
    }
    return prefix;
}

int congestion_dist(const Tree& tree, const std::vector<int>& prefix,
                    NodeId a, NodeId b) {
    if (a < 0 || b < 0) return INT_MAX / 4;
    NodeId w = tree.lca(a, b);
    if (w < 0) return INT_MAX / 4;
    return prefix[a] + prefix[b] - 2 * prefix[w];
}

bool paths_allowed(const Tree& tree, const std::vector<CompId>& owner,
                   const std::vector<Leaf>& leaves, NodeId root,
                   CompId allowed_a, CompId allowed_b) {
    for (Leaf leaf : leaves) {
        NodeId v = tree.leaf_node[leaf];
        while (v >= 0 && v != root) {
            CompId current = owner[v];
            if (current != -1 && current != allowed_a && current != allowed_b) {
                return false;
            }
            v = tree.parent[v];
        }
    }
    return true;
}

void mark_paths(const Tree& tree, std::vector<CompId>& owner,
                const std::vector<Leaf>& leaves, NodeId root, CompId id) {
    for (Leaf leaf : leaves) {
        NodeId v = tree.leaf_node[leaf];
        while (v >= 0 && v != root) {
            owner[v] = id;
            v = tree.parent[v];
        }
    }
}

bool rebuild_owners(Forest12& forest, const Tree& t1, const Tree& t2) {
    forest.owner1.assign(t1.node_count(), -1);
    forest.owner2.assign(t2.node_count(), -1);
    int active = 0;
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        Component12& component = forest.comps[id];
        if (!component.active) continue;
        ++active;
        component.root1 = t1.lca_leaves(component.leaves);
        component.root2 = t2.lca_leaves(component.leaves);
        if (component.root1 < 0 || component.root2 < 0) return false;
        if (component.leaves.size() <= 1) continue;
        if (!paths_allowed(t1, forest.owner1, component.leaves,
                           component.root1, id, id)) {
            return false;
        }
        if (!paths_allowed(t2, forest.owner2, component.leaves,
                           component.root2, id, id)) {
            return false;
        }
        mark_paths(t1, forest.owner1, component.leaves, component.root1, id);
        mark_paths(t2, forest.owner2, component.leaves, component.root2, id);
    }
    return active == forest.active_count;
}

Forest12 singleton_forest12(const Tree& t1, const Tree& t2) {
    Forest12 forest;
    forest.comps.reserve(t1.n);
    forest.leaf_comp.assign(t1.n + 1, -1);
    for (Leaf leaf = 1; leaf <= t1.n; ++leaf) {
        Component12 component;
        component.leaves.push_back(leaf);
        component.repr = leaf;
        component.root1 = t1.leaf_node[leaf];
        component.root2 = t2.leaf_node[leaf];
        forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
        forest.comps.push_back(std::move(component));
    }
    forest.active_count = t1.n;
    rebuild_owners(forest, t1, t2);
    return forest;
}

bool same_topology(const Tree& t1, const Tree& t2,
                   const std::vector<Leaf>& leaves) {
    InducedInfo i1 = compute_induced(t1, leaves);
    InducedInfo i2 = compute_induced(t2, leaves);
    return i1.signature == i2.signature;
}

int score_pair(const Forest12& forest, CompId a, CompId b,
               const Tree& t1, const Tree& t2,
               const std::vector<int>& prefix1,
               const std::vector<int>& prefix2,
               const Greedy12Config& config) {
    const Component12& ca = forest.comps[a];
    const Component12& cb = forest.comps[b];
    int d1 = congestion_dist(t1, prefix1, ca.root1, cb.root1);
    int d2 = congestion_dist(t2, prefix2, ca.root2, cb.root2);
    int size = static_cast<int>(ca.leaves.size() + cb.leaves.size());
    int depth_bonus = 0;
    if (ca.root1 >= 0) depth_bonus += t1.depth[ca.root1];
    if (cb.root1 >= 0) depth_bonus += t1.depth[cb.root1];
    if (ca.root2 >= 0) depth_bonus += t2.depth[ca.root2];
    if (cb.root2 >= 0) depth_bonus += t2.depth[cb.root2];
    int score = d1 + d2;
    switch (config.score_mode) {
        case 1:
            score = 2 * d1 + d2;
            break;
        case 2:
            score = d1 + 2 * d2;
            break;
        case 3:
            score = d1 + d2 - 2 * size;
            break;
        case 4:
            score = d1 + d2 - depth_bonus / 4;
            break;
        case 5:
            score = std::max(d1, d2) * 2 + std::min(d1, d2);
            break;
        default:
            break;
    }
    if (config.score_noise > 0) {
        uint64_t h = pair_key(ca.repr, cb.repr);
        h ^= static_cast<uint64_t>(config.seed) * 0x9e3779b97f4a7c15ULL;
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        score += static_cast<int>(h % static_cast<uint64_t>(config.score_noise));
    }
    return score;
}

void push_candidate(Queue12& queue, std::unordered_set<uint64_t>& seen,
                    const Forest12& forest, const Tree& t1, const Tree& t2,
                    const std::vector<int>& prefix1,
                    const std::vector<int>& prefix2,
                    const Greedy12Config& config, CompId a, CompId b) {
    if (a < 0 || b < 0 || a == b) return;
    if (!forest.comps[a].active || !forest.comps[b].active) return;
    Leaf ra = forest.comps[a].repr;
    Leaf rb = forest.comps[b].repr;
    uint64_t key = pair_key(ra, rb);
    if (!seen.insert(key).second) return;
    queue.push({score_pair(forest, a, b, t1, t2, prefix1, prefix2, config), ra, rb});
}

void push_capped_node_pairs(Queue12& queue, std::unordered_set<uint64_t>& seen,
                            const Forest12& forest, const Tree& t1, const Tree& t2,
                            const std::vector<int>& prefix1,
                            const std::vector<int>& prefix2,
                            const Greedy12Config& config,
                            const std::vector<CompId>& ids) {
    if (config.node_pair_cap <= 0 ||
        static_cast<int>(ids.size()) * (static_cast<int>(ids.size()) - 1) / 2 <=
            config.node_pair_cap) {
        for (size_t i = 0; i < ids.size(); ++i) {
            for (size_t j = i + 1; j < ids.size(); ++j) {
                push_candidate(queue, seen, forest, t1, t2, prefix1, prefix2,
                               config, ids[i], ids[j]);
            }
        }
        return;
    }

    std::vector<Candidate12> local;
    local.reserve(static_cast<size_t>(config.node_pair_cap) * 4);
    for (size_t i = 0; i < ids.size(); ++i) {
        for (size_t j = i + 1; j < ids.size(); ++j) {
            CompId a = ids[i];
            CompId b = ids[j];
            if (!forest.comps[a].active || !forest.comps[b].active) continue;
            local.push_back({
                score_pair(forest, a, b, t1, t2, prefix1, prefix2, config),
                forest.comps[a].repr,
                forest.comps[b].repr,
            });
        }
    }
    if (static_cast<int>(local.size()) > config.node_pair_cap) {
        std::nth_element(local.begin(), local.begin() + config.node_pair_cap,
                         local.end(), [](const Candidate12& a,
                                         const Candidate12& b) {
                             return a.score < b.score;
                         });
        local.resize(config.node_pair_cap);
    }
    std::sort(local.begin(), local.end(), [](const Candidate12& a,
                                             const Candidate12& b) {
        if (a.score != b.score) return a.score < b.score;
        if (a.a != b.a) return a.a < b.a;
        return a.b < b.b;
    });
    for (const Candidate12& candidate : local) {
        push_candidate(queue, seen, forest, t1, t2, prefix1, prefix2, config,
                       forest.leaf_comp[candidate.a],
                       forest.leaf_comp[candidate.b]);
    }
}

std::vector<CompId> local_components(const Tree& tree, const Forest12& forest,
                                     NodeId root, int limit) {
    std::vector<CompId> result;
    if (root < 0 || root >= tree.node_count() || limit <= 0) return result;
    std::unordered_set<CompId> seen;
    seen.reserve(static_cast<size_t>(limit) * 2 + 1);
    std::vector<NodeId> stack{root};
    while (!stack.empty() && static_cast<int>(result.size()) < limit) {
        NodeId v = stack.back();
        stack.pop_back();
        if (tree.is_leaf(v)) {
            CompId id = forest.leaf_comp[tree.label[v]];
            if (id >= 0 && forest.comps[id].active && seen.insert(id).second) {
                result.push_back(id);
            }
            continue;
        }
        if (tree.right[v] != -1) stack.push_back(tree.right[v]);
        if (tree.left[v] != -1) stack.push_back(tree.left[v]);
    }
    return result;
}

void seed_candidates(Queue12& queue, std::unordered_set<uint64_t>& seen,
                     const Forest12& forest, const Tree& t1, const Tree& t2,
                     const std::vector<int>& prefix1,
                     const std::vector<int>& prefix2,
                     const Greedy12Config& config,
                     std::mt19937& rng) {
    for (const Tree* tree : {&t1, &t2}) {
        for (NodeId v = 0; v < tree->node_count(); ++v) {
            if (tree->is_leaf(v)) continue;
            std::vector<CompId> ids = local_components(*tree, forest, v,
                                                       config.local_window);
            push_capped_node_pairs(queue, seen, forest, t1, t2,
                                   prefix1, prefix2, config, ids);
        }
    }

    if (t1.n < 2) return;
    unsigned range = static_cast<unsigned>(t1.n);
    unsigned mask = 1;
    while (mask < range - 1) mask = (mask << 1) | 1u;
    auto random_leaf = [&]() {
        unsigned value = 0;
        do {
            value = rng() & mask;
        } while (value >= range);
        return static_cast<Leaf>(value + 1);
    };
    int trials = t1.n * config.random_multiplier;
    for (int t = 0; t < trials; ++t) {
        Leaf la = random_leaf();
        Leaf lb = random_leaf();
        if (la == lb) continue;
        push_candidate(queue, seen, forest, t1, t2, prefix1, prefix2,
                       config, forest.leaf_comp[la], forest.leaf_comp[lb]);
    }
}

bool active_candidate(const Candidate12& candidate, const Forest12& forest,
                      CompId& a, CompId& b) {
    if (candidate.a <= 0 || candidate.b <= 0 ||
        candidate.a >= static_cast<int>(forest.leaf_comp.size()) ||
        candidate.b >= static_cast<int>(forest.leaf_comp.size())) {
        return false;
    }
    a = forest.leaf_comp[candidate.a];
    b = forest.leaf_comp[candidate.b];
    return a >= 0 && b >= 0 && a != b &&
           forest.comps[a].active && forest.comps[b].active;
}

bool pop_active(Queue12& queue, const Forest12& forest, Candidate12& out) {
    while (!queue.empty()) {
        out = queue.top();
        queue.pop();
        CompId a = -1;
        CompId b = -1;
        if (active_candidate(out, forest, a, b)) return true;
    }
    return false;
}

CompId try_merge(Forest12& forest, const Tree& t1, const Tree& t2,
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
    NodeId root1 = t1.lca_leaves(merged);
    NodeId root2 = t2.lca_leaves(merged);
    if (root1 < 0 || root2 < 0) return -1;
    if (!same_topology(t1, t2, merged)) return -1;
    if (!paths_allowed(t1, forest.owner1, merged, root1, a, b)) return -1;
    if (!paths_allowed(t2, forest.owner2, merged, root2, a, b)) return -1;

    Component12 component;
    component.leaves = std::move(merged);
    component.repr = forest.comps[a].repr;
    component.root1 = root1;
    component.root2 = root2;

    forest.comps[a].active = false;
    forest.comps[b].active = false;
    forest.comps[a].leaves.clear();
    forest.comps[b].leaves.clear();

    CompId id = static_cast<CompId>(forest.comps.size());
    for (Leaf leaf : component.leaves) forest.leaf_comp[leaf] = id;
    forest.comps.push_back(std::move(component));
    --forest.active_count;
    mark_paths(t1, forest.owner1, forest.comps[id].leaves, root1, id);
    mark_paths(t2, forest.owner2, forest.comps[id].leaves, root2, id);
    return id;
}

void push_neighbors(Queue12& queue, std::unordered_set<uint64_t>& seen,
                    std::unordered_set<uint64_t>& failed,
                    const Forest12& forest, const Tree& t1, const Tree& t2,
                    const std::vector<int>& prefix1,
                    const std::vector<int>& prefix2,
                    const Greedy12Config& config,
                    CompId focus) {
    if (focus < 0 || !forest.comps[focus].active) return;
    std::unordered_set<CompId> refreshed;
    refreshed.reserve(static_cast<size_t>(config.local_window) * 4 + 1);
    for (const Tree* tree : {&t1, &t2}) {
        std::unordered_set<NodeId> touched;
        for (Leaf leaf : forest.comps[focus].leaves) {
            NodeId v = tree->leaf_node[leaf];
            for (int walk = 0; v >= 0 && walk < kNeighborWalkLimit; ++walk) {
                if (touched.insert(v).second) {
                    std::vector<CompId> ids = local_components(*tree, forest, v,
                                                               config.local_window);
                    std::vector<Candidate12> local;
                    for (CompId other : ids) {
                        if (other == focus || !forest.comps[other].active) continue;
                        if (refreshed.count(other) != 0) continue;
                        local.push_back({
                            score_pair(forest, focus, other, t1, t2,
                                       prefix1, prefix2, config),
                            forest.comps[focus].repr,
                            forest.comps[other].repr,
                        });
                    }
                    int cap = config.neighbor_pair_cap;
                    if (cap > 0 && static_cast<int>(local.size()) > cap) {
                        std::nth_element(local.begin(), local.begin() + cap,
                                         local.end(), [](const Candidate12& a,
                                                         const Candidate12& b) {
                                             return a.score < b.score;
                                         });
                        local.resize(cap);
                    }
                    std::sort(local.begin(), local.end(), [](const Candidate12& a,
                                                             const Candidate12& b) {
                        if (a.score != b.score) return a.score < b.score;
                        return a.b < b.b;
                    });
                    for (const Candidate12& candidate : local) {
                        CompId other = forest.leaf_comp[candidate.b];
                        if (other < 0 || other == focus ||
                            !forest.comps[other].active) {
                            continue;
                        }
                        if (!refreshed.insert(other).second) continue;
                        uint64_t key = pair_key(forest.comps[focus].repr,
                                                forest.comps[other].repr);
                        failed.erase(key);
                        seen.erase(key);
                        push_candidate(queue, seen, forest, t1, t2, prefix1, prefix2,
                                       config, focus, other);
                    }
                }
                v = tree->parent[v];
            }
        }
    }
}

void greedy_merge12(Forest12& forest, const Tree& t1, const Tree& t2,
                    const Greedy12Config& config) {
    std::vector<int> prefix1 = congestion_prefix(t1);
    std::vector<int> prefix2 = congestion_prefix(t2);
    std::mt19937 rng(config.seed);

    Queue12 queue;
    std::unordered_set<uint64_t> seen;
    std::unordered_set<uint64_t> failed;
    seen.reserve(static_cast<size_t>(t1.n) * 8 + 1);
    failed.reserve(static_cast<size_t>(t1.n) * 2 + 1);
    seed_candidates(queue, seen, forest, t1, t2, prefix1, prefix2, config, rng);

    int max_checks = std::max(kMinMergeChecks, t1.n * config.merge_check_multiplier);
    int checks = 0;
    while (!queue.empty() && checks < max_checks) {
        Candidate12 candidate;
        if (!pop_active(queue, forest, candidate)) break;
        CompId a = -1;
        CompId b = -1;
        if (!active_candidate(candidate, forest, a, b)) continue;
        ++checks;
        CompId merged = try_merge(forest, t1, t2, a, b);
        if (merged >= 0) {
            failed.erase(pair_key(candidate.a, candidate.b));
            push_neighbors(queue, seen, failed, forest, t1, t2,
                           prefix1, prefix2, config, merged);
        } else {
            failed.insert(pair_key(candidate.a, candidate.b));
        }
    }
}

Forest to_rs_forest(const Forest12& source) {
    Forest out;
    out.leaf_comp.assign(source.leaf_comp.size(), -1);
    for (const Component12& component : source.comps) {
        if (!component.active) continue;
        Component dst;
        dst.leaves = component.leaves;
        dst.repr = component.repr;
        CompId id = static_cast<CompId>(out.comps.size());
        for (Leaf leaf : dst.leaves) out.leaf_comp[leaf] = id;
        out.comps.push_back(std::move(dst));
        ++out.active_count;
    }
    return out;
}

Forest12 from_rs_forest(const Forest& source, int leaf_count,
                        const Tree& t1, const Tree& t2) {
    Forest12 out;
    out.leaf_comp.assign(leaf_count + 1, -1);
    for (const Component& component : source.comps) {
        if (!component.active) continue;
        Component12 dst;
        dst.leaves = component.leaves;
        dst.repr = component.leaves.empty() ? -1 : component.leaves.front();
        CompId id = static_cast<CompId>(out.comps.size());
        for (Leaf leaf : dst.leaves) out.leaf_comp[leaf] = id;
        out.comps.push_back(std::move(dst));
        ++out.active_count;
    }
    rebuild_owners(out, t1, t2);
    return out;
}

Forest12 split_largest_components(const Forest12& source, int max_components,
                                  int max_split_leaves,
                                  const Tree& t1, const Tree& t2) {
    std::vector<CompId> ids;
    for (CompId id = 0; id < static_cast<CompId>(source.comps.size()); ++id) {
        if (source.comps[id].active && source.comps[id].leaves.size() > 1) {
            ids.push_back(id);
        }
    }
    std::sort(ids.begin(), ids.end(), [&](CompId a, CompId b) {
        if (source.comps[a].leaves.size() != source.comps[b].leaves.size()) {
            return source.comps[a].leaves.size() > source.comps[b].leaves.size();
        }
        return source.comps[a].repr < source.comps[b].repr;
    });
    std::unordered_set<CompId> split;
    int split_leaves = 0;
    for (CompId id : ids) {
        int size = static_cast<int>(source.comps[id].leaves.size());
        if (static_cast<int>(split.size()) >= max_components) break;
        if (max_split_leaves > 0 && split_leaves + size > max_split_leaves) break;
        split.insert(id);
        split_leaves += size;
    }

    Forest12 out;
    out.leaf_comp.assign(t1.n + 1, -1);
    for (CompId id = 0; id < static_cast<CompId>(source.comps.size()); ++id) {
        const Component12& component = source.comps[id];
        if (!component.active) continue;
        if (split.count(id) == 0) {
            Component12 dst;
            dst.leaves = component.leaves;
            dst.repr = component.repr;
            CompId next = static_cast<CompId>(out.comps.size());
            for (Leaf leaf : dst.leaves) out.leaf_comp[leaf] = next;
            out.comps.push_back(std::move(dst));
            ++out.active_count;
            continue;
        }
        for (Leaf leaf : component.leaves) {
            Component12 singleton;
            singleton.leaves.push_back(leaf);
            singleton.repr = leaf;
            singleton.root1 = t1.leaf_node[leaf];
            singleton.root2 = t2.leaf_node[leaf];
            out.leaf_comp[leaf] = static_cast<CompId>(out.comps.size());
            out.comps.push_back(std::move(singleton));
            ++out.active_count;
        }
    }
    rebuild_owners(out, t1, t2);
    return out;
}

Forest run_one_greedy12_variant(const Tree& t1, const Tree& t2,
                                const Greedy12Config& config) {
    Forest12 forest = singleton_forest12(t1, t2);
    if (!rebuild_owners(forest, t1, t2)) return Forest();
    greedy_merge12(forest, t1, t2, config);
    Forest out = to_rs_forest(forest);
    rebuild_forest_data(out, t1, t2);
    return out;
}

double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
}

}  // namespace

bool should_run_greedy12_tail(int leaves, int residual_leaves) {
    if (env_int("RS_GREEDY12_TAIL", 1) == 0) return false;
    int min_leaves = env_int("RS_GREEDY12_MIN_LEAVES", 9000);
    int min_residual_percent =
        env_int("RS_GREEDY12_MIN_RESIDUAL_PERCENT", 0);
    return leaves >= min_leaves &&
           residual_leaves * 100LL >=
               static_cast<long long>(min_residual_percent) * leaves;
}

bool run_greedy12_tail(const Tree& t1, const Tree& t2, Forest& out) {
    const bool debug = std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr;
    const int portfolio = env_int("RS_GREEDY12_PORTFOLIO", 0);
    if (portfolio == 0) {
        Greedy12Config config;
        config.seed = env_int("RS_GREEDY12_SEED", 1);
        config.local_window = env_int("RS_GREEDY12_LOCAL_WINDOW", kLocalWindow);
        config.local_window =
            env_int("RS_GREEDY12_CAPPED_WINDOW", config.local_window);
        config.node_pair_cap = env_int("RS_GREEDY12_NODE_PAIR_CAP", 0);
        config.neighbor_pair_cap = env_int("RS_GREEDY12_NEIGHBOR_PAIR_CAP", 0);
        config.random_multiplier =
            env_int("RS_GREEDY12_RANDOM_MULTIPLIER", kRandomCandidateMultiplier);
        config.merge_check_multiplier =
            env_int("RS_GREEDY12_CHECK_MULTIPLIER", kMergeCheckMultiplier);
        config.score_mode = env_int("RS_GREEDY12_SCORE_MODE", 0);
        config.score_noise = env_int("RS_GREEDY12_SCORE_NOISE", 0);
        out = run_one_greedy12_variant(t1, t2, config);
        return out.active_count > 0 && rebuild_forest_data(out, t1, t2);
    }

    const auto start = std::chrono::steady_clock::now();
    const int max_variants = env_int("RS_GREEDY12_VARIANTS", 4);
    const double max_seconds =
        static_cast<double>(env_int("RS_GREEDY12_SECONDS", 275));
    const int base_seed = env_int("RS_GREEDY12_SEED", 1);
    const int base_window = env_int("RS_GREEDY12_LOCAL_WINDOW", kLocalWindow);
    const int base_checks =
        env_int("RS_GREEDY12_CHECK_MULTIPLIER", kMergeCheckMultiplier);
    const int base_random =
        env_int("RS_GREEDY12_RANDOM_MULTIPLIER", kRandomCandidateMultiplier);

    bool have_best = false;
    Forest best;
    for (int variant = 0; variant < max_variants; ++variant) {
        if (variant > 0 && seconds_since(start) >= max_seconds) break;
        Greedy12Config config;
        config.seed = (variant == 1) ? base_seed : base_seed + variant * 1009;
        config.local_window = base_window;
        config.local_window =
            env_int("RS_GREEDY12_CAPPED_WINDOW", config.local_window);
        config.node_pair_cap = env_int("RS_GREEDY12_NODE_PAIR_CAP", 0);
        config.neighbor_pair_cap = env_int("RS_GREEDY12_NEIGHBOR_PAIR_CAP", 0);
        config.random_multiplier = base_random;
        config.merge_check_multiplier = base_checks;
        config.score_mode = (variant == 0)
            ? env_int("RS_GREEDY12_FIRST_MODE", 0)
            : (variant % 6);
        config.score_noise =
            (variant <= 1) ? 0 : env_int("RS_GREEDY12_SCORE_NOISE", 17);
        if (variant == 1) config.score_mode = 5;
        if (variant == 2) config.score_mode = 4;
        if (variant == 3) {
            config.score_mode = 3;
            config.random_multiplier = std::max(base_random, 3);
            config.merge_check_multiplier = std::max(base_checks, 95);
        }
        if (variant >= 4) {
            config.score_mode = variant % 6;
            config.random_multiplier = std::max(base_random, 3 + (variant % 3));
        }

        Forest candidate = run_one_greedy12_variant(t1, t2, config);
        const bool valid =
            candidate.active_count > 0 && rebuild_forest_data(candidate, t1, t2);
        if (debug) {
            std::cerr << "rs_greedy12_variant variant=" << variant
                      << " valid=" << (valid ? 1 : 0)
                      << " components=" << (valid ? candidate.active_count : -1)
                      << " seed=" << config.seed
                      << " window=" << config.local_window
                      << " mode=" << config.score_mode
                      << " random=" << config.random_multiplier
                      << " checks=" << config.merge_check_multiplier
                      << " seconds=" << seconds_since(start) << '\n';
        }
        if (valid && (!have_best || candidate.active_count < best.active_count)) {
            best = std::move(candidate);
            have_best = true;
            publish_best_output(forest_to_output(best, t1));
        }
    }
    const int destroy_passes = env_int("RS_GREEDY12_DESTROY_PASSES", 0);
    for (int pass = 0; have_best && pass < destroy_passes; ++pass) {
        if (seconds_since(start) >= max_seconds) break;
        const int split_components =
            env_int("RS_GREEDY12_DESTROY_COMPONENTS", 16) * (pass + 1);
        const int split_leaves =
            env_int("RS_GREEDY12_DESTROY_LEAVES", 900) * (pass + 1);
        Forest12 seed = from_rs_forest(best, t1.n, t1, t2);
        seed = split_largest_components(seed, split_components, split_leaves,
                                        t1, t2);
        Greedy12Config config;
        config.seed = base_seed + 7001 + pass * 977;
        config.local_window = base_window;
        config.local_window =
            env_int("RS_GREEDY12_CAPPED_WINDOW", config.local_window);
        config.node_pair_cap = env_int("RS_GREEDY12_NODE_PAIR_CAP", 0);
        config.neighbor_pair_cap = env_int("RS_GREEDY12_NEIGHBOR_PAIR_CAP", 0);
        config.random_multiplier =
            env_int("RS_GREEDY12_DESTROY_RANDOM_MULTIPLIER",
                    std::max(base_random, 4));
        config.merge_check_multiplier =
            env_int("RS_GREEDY12_DESTROY_CHECK_MULTIPLIER",
                    std::max(base_checks, 100));
        config.score_mode = env_int("RS_GREEDY12_DESTROY_SCORE_MODE", 5);
        config.score_noise = env_int("RS_GREEDY12_DESTROY_SCORE_NOISE", 0);
        greedy_merge12(seed, t1, t2, config);
        Forest candidate = to_rs_forest(seed);
        const bool valid =
            candidate.active_count > 0 && rebuild_forest_data(candidate, t1, t2);
        if (debug) {
            std::cerr << "rs_greedy12_destroy pass=" << pass
                      << " valid=" << (valid ? 1 : 0)
                      << " components=" << (valid ? candidate.active_count : -1)
                      << " split_components=" << split_components
                      << " split_leaves=" << split_leaves
                      << " mode=" << config.score_mode
                      << " seconds=" << seconds_since(start) << '\n';
        }
        if (valid && candidate.active_count < best.active_count) {
            best = std::move(candidate);
            publish_best_output(forest_to_output(best, t1));
        }
    }
    if (!have_best) return false;
    out = std::move(best);
    return rebuild_forest_data(out, t1, t2);
}
