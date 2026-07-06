#include "heuristic/highs_pool_repair.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <cstdio>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <unistd.h>

#include "Highs.h"

#include "forest/validator.h"
#include "tree/induced.h"
#include "util/signal_handler.h"

namespace {

int env_int_hp(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

double env_double_hp(const char* name, double fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    double value = std::strtod(raw, &end);
    if (end == raw || *end != '\0' || value <= 0.0) return fallback;
    return value;
}

double elapsed_seconds(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
        .count();
}

std::string key_of(const std::vector<Leaf>& leaves) {
    std::string key;
    for (size_t i = 0; i < leaves.size(); ++i) {
        if (i) key += ',';
        key += std::to_string(leaves[i]);
    }
    return key;
}

int leaf_distance(const Tree& tree, Leaf a, Leaf b) {
    NodeId na = tree.leaf_node[a];
    NodeId nb = tree.leaf_node[b];
    NodeId root = tree.lca(na, nb);
    return tree.depth[na] + tree.depth[nb] - 2 * tree.depth[root];
}

int node_distance(const Tree& tree, NodeId a, NodeId b) {
    NodeId root = tree.lca(a, b);
    return tree.depth[a] + tree.depth[b] - 2 * tree.depth[root];
}

std::vector<Leaf> normalized(std::vector<Leaf> leaves) {
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    return leaves;
}

std::vector<Leaf> set_intersection_vec(const std::vector<Leaf>& a,
                                       const std::vector<Leaf>& b) {
    std::vector<Leaf> out;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(out));
    return out;
}

std::vector<Leaf> set_difference_vec(const std::vector<Leaf>& a,
                                     const std::vector<Leaf>& b) {
    std::vector<Leaf> out;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(out));
    return out;
}

struct Column {
    std::vector<Leaf> leaves;
    std::vector<NodeId> edges1;
    std::vector<NodeId> edges2;
};

struct Region {
    std::vector<Leaf> leaves;
    std::vector<CompId> comp_ids;
    int score = 0;
};

struct Pool {
    const Tree& t1;
    const Tree& t2;
    int max_pool;
    int max_component_leaves;
    std::unordered_map<std::string, int> index;
    std::vector<Column> columns;

    Pool(const Tree& first, const Tree& second, int pool_limit, int component_limit)
        : t1(first),
          t2(second),
          max_pool(pool_limit),
          max_component_leaves(component_limit) {}

    bool add(std::vector<Leaf> leaves) {
        if (static_cast<int>(columns.size()) >= max_pool) return false;
        leaves = normalized(std::move(leaves));
        if (leaves.empty()) return false;
        if (static_cast<int>(leaves.size()) > max_component_leaves) return false;
        std::string key = key_of(leaves);
        if (index.find(key) != index.end()) return false;
        InducedInfo i1;
        InducedInfo i2;
        try {
            i1 = compute_induced(t1, leaves);
            i2 = compute_induced(t2, leaves);
        } catch (...) {
            return false;
        }
        if (i1.signature != i2.signature) return false;
        Column col;
        col.leaves = std::move(leaves);
        col.edges1 = std::move(i1.used_edges);
        col.edges2 = std::move(i2.used_edges);
        index.emplace(key, static_cast<int>(columns.size()));
        columns.push_back(std::move(col));
        return true;
    }
};

int component_score(const std::vector<Leaf>& leaves, const Tree& t1,
                    const Tree& t2) {
    if (leaves.empty()) return 0;
    NodeId r1 = t1.lca_leaves(leaves);
    NodeId r2 = t2.lca_leaves(leaves);
    int depth = (r1 >= 0 ? t1.depth[r1] : 0) + (r2 >= 0 ? t2.depth[r2] : 0);
    return 100000 * (static_cast<int>(leaves.size()) - 1) + 1200 * depth;
}

Forest forest_from_columns(const std::vector<Column>& columns,
                           const std::vector<int>& selected,
                           int leaf_count,
                           const Tree& t1,
                           const Tree& t2) {
    Forest forest;
    forest.leaf_comp.assign(leaf_count + 1, -1);
    forest.comps.reserve(selected.size());
    for (int idx : selected) {
        Component c;
        c.active = true;
        c.leaves = columns[idx].leaves;
        c.repr = c.leaves.empty() ? -1 : c.leaves.front();
        CompId cid = static_cast<CompId>(forest.comps.size());
        for (Leaf x : c.leaves) forest.leaf_comp[x] = cid;
        forest.comps.push_back(std::move(c));
    }
    forest.active_count = static_cast<int>(forest.comps.size());
    if (!rebuild_forest_data(forest, t1, t2)) return Forest{};
    return forest;
}

std::vector<Leaf> collect_limited_leaves(const Tree& tree, NodeId root,
                                         int limit) {
    std::vector<Leaf> leaves;
    std::vector<NodeId> stack{root};
    while (!stack.empty()) {
        NodeId v = stack.back();
        stack.pop_back();
        if (tree.is_leaf(v)) {
            leaves.push_back(tree.label[v]);
            if (static_cast<int>(leaves.size()) > limit) return {};
        } else {
            if (tree.left[v] != -1) stack.push_back(tree.left[v]);
            if (tree.right[v] != -1) stack.push_back(tree.right[v]);
        }
    }
    return normalized(std::move(leaves));
}

std::vector<std::vector<Leaf>> restricted_node_sets(
    const Tree& tree, const std::vector<Leaf>& region_leaves, int max_set,
    int keep_limit) {
    std::vector<char> mark(tree.n + 1, 0);
    for (Leaf x : region_leaves) mark[x] = 1;

    std::vector<NodeId> order;
    std::vector<NodeId> stack{tree.root};
    while (!stack.empty()) {
        NodeId v = stack.back();
        stack.pop_back();
        order.push_back(v);
        if (!tree.is_leaf(v)) {
            stack.push_back(tree.left[v]);
            stack.push_back(tree.right[v]);
        }
    }

    std::vector<std::vector<Leaf>> below(tree.node_count());
    std::vector<std::vector<Leaf>> sets;
    for (int oi = static_cast<int>(order.size()) - 1; oi >= 0; --oi) {
        NodeId v = order[oi];
        std::vector<Leaf> leaves;
        if (tree.is_leaf(v)) {
            Leaf x = tree.label[v];
            if (x >= 1 && x <= tree.n && mark[x]) leaves.push_back(x);
        } else {
            leaves.insert(leaves.end(), below[tree.left[v]].begin(),
                          below[tree.left[v]].end());
            leaves.insert(leaves.end(), below[tree.right[v]].begin(),
                          below[tree.right[v]].end());
            leaves = normalized(std::move(leaves));
        }
        if (static_cast<int>(leaves.size()) <= max_set) below[v] = leaves;
        if (leaves.size() >= 2 && static_cast<int>(leaves.size()) <= max_set) {
            sets.push_back(std::move(leaves));
        }
    }
    std::sort(sets.begin(), sets.end(),
              [&](const auto& a, const auto& b) {
                  if (a.size() != b.size()) return a.size() > b.size();
                  return component_score(a, tree, tree) >
                         component_score(b, tree, tree);
              });
    sets.erase(std::unique(sets.begin(), sets.end()), sets.end());
    if (static_cast<int>(sets.size()) > keep_limit) sets.resize(keep_limit);
    return sets;
}

std::vector<Region> collect_regions(const Forest& best, const Tree& t1,
                                    const Tree& t2) {
    const int region_limit = env_int_hp("RS_HIGHS_REGION_LIMIT", 90);
    const int component_limit =
        env_int_hp("RS_HIGHS_REGION_COMPONENT_LIMIT", 18);
    const int min_leaves = env_int_hp("RS_HIGHS_REGION_MIN_LEAVES", 6);
    const int max_leaves = env_int_hp("RS_HIGHS_REGION_MAX_LEAVES", 42);

    std::vector<CompId> active = active_components(best);
    std::vector<int> leaf_to_pos(t1.n + 1, -1);
    for (int pos = 0; pos < static_cast<int>(active.size()); ++pos) {
        const Component& c = best.comps[active[pos]];
        for (Leaf x : c.leaves) leaf_to_pos[x] = pos;
    }

    std::unordered_map<std::string, Region> by_key;
    auto add_region = [&](std::vector<int> ids, const Tree& guide, NodeId node) {
        ids = normalized(std::move(ids));
        if (ids.size() < 2 || static_cast<int>(ids.size()) > component_limit) {
            return;
        }
        std::vector<Leaf> leaves;
        for (int pos : ids) {
            const auto& comp = best.comps[active[pos]].leaves;
            leaves.insert(leaves.end(), comp.begin(), comp.end());
        }
        leaves = normalized(std::move(leaves));
        if (static_cast<int>(leaves.size()) < min_leaves ||
            static_cast<int>(leaves.size()) > max_leaves) {
            return;
        }
        std::vector<CompId> comp_ids;
        comp_ids.reserve(ids.size());
        for (int pos : ids) comp_ids.push_back(active[pos]);
        std::sort(comp_ids.begin(), comp_ids.end());
        comp_ids.erase(std::unique(comp_ids.begin(), comp_ids.end()),
                       comp_ids.end());

        std::string key = key_of(leaves);
        int score = 4096 * static_cast<int>(ids.size()) -
                    17 * static_cast<int>(leaves.size()) + guide.depth[node];
        auto it = by_key.find(key);
        if (it == by_key.end() || score > it->second.score) {
            by_key[key] = Region{std::move(leaves), std::move(comp_ids), score};
        }
    };

    auto scan_tree = [&](const Tree& tree) {
        std::function<std::vector<int>(NodeId)> dfs = [&](NodeId v) {
            if (tree.is_leaf(v)) {
                int pos = leaf_to_pos[tree.label[v]];
                return pos < 0 ? std::vector<int>{} : std::vector<int>{pos};
            }
            std::vector<int> left = dfs(tree.left[v]);
            std::vector<int> right = dfs(tree.right[v]);
            left.insert(left.end(), right.begin(), right.end());
            left = normalized(std::move(left));
            if (left.size() > static_cast<size_t>(component_limit)) {
                left.resize(component_limit + 1);
                return left;
            }
            add_region(left, tree, v);
            return left;
        };
        dfs(tree.root);
    };
    scan_tree(t1);
    scan_tree(t2);

    for (CompId cid : active) {
        for (const Tree* tree : {&t1, &t2}) {
            NodeId node = tree->lca_leaves(best.comps[cid].leaves);
            for (int step = 0; node >= 0 && step <= 4;
                 ++step, node = tree->parent[node]) {
                std::vector<Leaf> leaves =
                    collect_limited_leaves(*tree, node, max_leaves);
                if (leaves.empty()) continue;
                std::vector<int> ids;
                for (Leaf x : leaves) {
                    int pos = leaf_to_pos[x];
                    if (pos >= 0) ids.push_back(pos);
                }
                add_region(std::move(ids), *tree, node);
            }
        }
    }

    std::vector<Region> regions;
    regions.reserve(by_key.size());
    for (auto& entry : by_key) regions.push_back(std::move(entry.second));
    std::sort(regions.begin(), regions.end(), [](const Region& a, const Region& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
        return a.leaves < b.leaves;
    });
    if (static_cast<int>(regions.size()) > region_limit) regions.resize(region_limit);
    return regions;
}

std::vector<Leaf> union_component_leaves(const Forest& forest,
                                         const std::vector<CompId>& ids) {
    std::vector<Leaf> leaves;
    for (CompId cid : ids) {
        if (cid < 0 || cid >= static_cast<CompId>(forest.comps.size())) continue;
        const Component& comp = forest.comps[cid];
        if (!comp.active) continue;
        leaves.insert(leaves.end(), comp.leaves.begin(), comp.leaves.end());
    }
    return normalized(std::move(leaves));
}

void add_region_union_pool(Pool& pool, const Forest& forest,
                           const Region& region, int& added) {
    const int union_leaf_limit =
        env_int_hp("RS_HIGHS_COMPONENT_UNION_LEAVES", 72);
    const int enum_limit =
        env_int_hp("RS_HIGHS_COMPONENT_UNION_REGION_LIMIT", 18);
    std::vector<CompId> ids = region.comp_ids;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (ids.size() < 2 || static_cast<int>(ids.size()) > enum_limit) return;

    auto add_ids = [&](std::vector<CompId> chosen) {
        std::vector<Leaf> leaves = union_component_leaves(forest, chosen);
        if (leaves.empty() ||
            static_cast<int>(leaves.size()) > union_leaf_limit) {
            return;
        }
        if (pool.add(std::move(leaves))) ++added;
    };

    const int m = static_cast<int>(ids.size());
    for (int a = 0; a < m; ++a) {
        for (int b = a + 1; b < m; ++b) {
            add_ids({ids[a], ids[b]});
        }
    }
    for (int a = 0; a < m; ++a) {
        for (int b = a + 1; b < m; ++b) {
            for (int c = b + 1; c < m; ++c) {
                add_ids({ids[a], ids[b], ids[c]});
            }
        }
    }
    if (m <= enum_limit) {
        for (int a = 0; a < m; ++a) {
            for (int b = a + 1; b < m; ++b) {
                for (int c = b + 1; c < m; ++c) {
                    for (int d = c + 1; d < m; ++d) {
                        add_ids({ids[a], ids[b], ids[c], ids[d]});
                    }
                }
            }
        }
    }
}

void add_region_close_leaf_pool(Pool& pool, const Region& region, int& added) {
    const int distance_limit =
        env_int_hp("RS_HIGHS_REGION_LEAF_DISTANCE", 10);
    const int pair_limit =
        env_int_hp("RS_HIGHS_REGION_CLOSE_PAIR_LIMIT", 1200);
    std::vector<Leaf> leaves = region.leaves;
    if (leaves.size() < 2) return;
    std::vector<std::vector<int>> nearby(leaves.size());
    int pairs_seen = 0;
    for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(leaves.size()); ++j) {
            Leaf a = leaves[i];
            Leaf b = leaves[j];
            int distance =
                std::max(leaf_distance(pool.t1, a, b),
                         leaf_distance(pool.t2, a, b));
            if (distance > distance_limit) continue;
            nearby[i].push_back(j);
            if (pool.add({a, b})) ++added;
            if (++pairs_seen >= pair_limit) return;
        }
    }

    for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
        for (int left = 0; left < static_cast<int>(nearby[i].size()); ++left) {
            for (int right = left + 1;
                 right < static_cast<int>(nearby[i].size()); ++right) {
                int j = nearby[i][left];
                int k = nearby[i][right];
                Leaf a = leaves[i];
                Leaf b = leaves[j];
                Leaf c = leaves[k];
                int distance =
                    std::max(leaf_distance(pool.t1, b, c),
                             leaf_distance(pool.t2, b, c));
                if (distance <= distance_limit && pool.add({a, b, c})) {
                    ++added;
                }
            }
        }
    }
}

void add_region_pool(Pool& pool, const Forest& forest,
                     const std::vector<Region>& regions,
                     std::chrono::steady_clock::time_point start,
                     double deadline) {
    const int local_pool = env_int_hp("RS_HIGHS_LOCAL_POOL_LIMIT", 5000);
    const int max_set = env_int_hp("RS_HIGHS_NODE_SET_MAX_LEAVES", 42);
    const int keep = env_int_hp("RS_HIGHS_NODE_SET_KEEP", 120);
    const int clade_per_region =
        env_int_hp("RS_HIGHS_CLADE_PER_REGION_LIMIT", 80);

    for (const Region& region : regions) {
        if (static_cast<int>(pool.columns.size()) >= local_pool ||
            elapsed_seconds(start) >= deadline) {
            break;
        }
        int added = 0;
        add_region_union_pool(pool, forest, region, added);
        add_region_close_leaf_pool(pool, region, added);
        auto first_sets = restricted_node_sets(pool.t1, region.leaves, max_set, keep);
        auto second_sets = restricted_node_sets(pool.t2, region.leaves, max_set, keep);
        for (const auto& leaves : first_sets) {
            if (pool.add(leaves)) ++added;
            if (added >= clade_per_region) break;
        }
        for (const auto& leaves : second_sets) {
            if (pool.add(leaves)) ++added;
            if (added >= clade_per_region) break;
        }
        for (const auto& a : first_sets) {
            for (const auto& b : second_sets) {
                pool.add(set_intersection_vec(a, b));
                pool.add(set_difference_vec(a, b));
                pool.add(set_difference_vec(b, a));
                if (++added >= clade_per_region ||
                    static_cast<int>(pool.columns.size()) >= local_pool ||
                    elapsed_seconds(start) >= deadline) {
                    break;
                }
            }
            if (added >= clade_per_region ||
                static_cast<int>(pool.columns.size()) >= local_pool ||
                elapsed_seconds(start) >= deadline) {
                break;
            }
        }
    }
}

std::vector<std::vector<int>> build_adjacency(const Tree& t1, const Tree& t2,
                                              int distance_limit) {
    std::vector<std::vector<int>> adjacency(t1.n + 1);
    for (Leaf a = 1; a <= t1.n; ++a) {
        for (Leaf b = a + 1; b <= t1.n; ++b) {
            int d = std::max(leaf_distance(t1, a, b), leaf_distance(t2, a, b));
            if (d <= distance_limit) {
                adjacency[a].push_back(b);
                adjacency[b].push_back(a);
            }
        }
    }
    for (auto& row : adjacency) std::sort(row.begin(), row.end());
    return adjacency;
}

bool adjacent_to_all(const std::vector<std::vector<int>>& adjacency, Leaf x,
                     const std::vector<Leaf>& state) {
    const auto& row = adjacency[x];
    for (Leaf y : state) {
        if (!std::binary_search(row.begin(), row.end(), y)) return false;
    }
    return true;
}

void add_bucket_clique_pool(Pool& pool,
                            std::chrono::steady_clock::time_point start,
                            double deadline) {
    const int n = pool.t1.n;
    const int distance = env_int_hp("RS_HIGHS_BUCKET_DISTANCE", 22);
    const int max_component_size =
        env_int_hp("RS_HIGHS_BUCKET_MAX_COMPONENT_SIZE", 16);
    const int selected_limit = env_int_hp("RS_HIGHS_BUCKET_SEED_LIMIT", 6500);
    const int add_limit = env_int_hp("RS_HIGHS_BUCKET_ADD_LIMIT", 2800);
    const int beam_width = env_int_hp("RS_HIGHS_BUCKET_BEAM_WIDTH", 20);
    const int option_limit = env_int_hp("RS_HIGHS_BUCKET_OPTION_LIMIT", 40);
    const int bucket_step = std::max(1, env_int_hp("RS_HIGHS_BUCKET_STEP", 250));
    const int bucket_half = std::max(0, env_int_hp("RS_HIGHS_BUCKET_HALF_WIDTH", 20));

    auto adjacency = build_adjacency(pool.t1, pool.t2, distance);
    std::vector<std::vector<Leaf>> seeds;
    for (Leaf a = 1; a <= n; ++a) {
        for (Leaf b : adjacency[a]) {
            if (a < b) seeds.push_back({a, b});
        }
    }
    std::sort(seeds.begin(), seeds.end(), [&](const auto& a, const auto& b) {
        int sa = component_score(a, pool.t1, pool.t2);
        int sb = component_score(b, pool.t1, pool.t2);
        if (sa != sb) return sa > sb;
        return a < b;
    });

    std::vector<std::vector<Leaf>> selected;
    std::vector<char> used(seeds.size(), 0);
    auto add_seed_index = [&](int idx) {
        if (idx < 0 || idx >= static_cast<int>(seeds.size()) || used[idx]) return;
        used[idx] = 1;
        selected.push_back(seeds[idx]);
    };
    int bucket_count = (static_cast<int>(seeds.size()) + bucket_step - 1) / bucket_step;
    for (int delta = 0;
         delta <= bucket_half && static_cast<int>(selected.size()) < selected_limit;
         ++delta) {
        for (int bucket = 0;
             bucket < bucket_count && static_cast<int>(selected.size()) < selected_limit;
             ++bucket) {
            int center = bucket * bucket_step;
            add_seed_index(center + delta);
            if (delta) add_seed_index(center - delta);
        }
    }

    int accepted = 0;
    std::unordered_set<std::string> global_seen;
    for (const auto& seed : selected) {
        if (accepted >= add_limit || static_cast<int>(pool.columns.size()) >= pool.max_pool ||
            elapsed_seconds(start) >= deadline) {
            break;
        }
        std::vector<std::vector<Leaf>> beam{seed};
        std::unordered_set<std::string> local_seen;
        local_seen.insert(key_of(seed));
        pool.add(seed);

        while (!beam.empty() && accepted < add_limit &&
               static_cast<int>(pool.columns.size()) < pool.max_pool &&
               elapsed_seconds(start) < deadline) {
            std::vector<std::vector<Leaf>> next_states;
            for (const auto& state : beam) {
                if (static_cast<int>(state.size()) >= max_component_size) continue;
                std::vector<std::pair<int, Leaf>> options;
                for (Leaf x = 1; x <= n; ++x) {
                    if (std::binary_search(state.begin(), state.end(), x)) continue;
                    if (!adjacent_to_all(adjacency, x, state)) continue;
                    std::vector<Leaf> grown = state;
                    grown.push_back(x);
                    grown = normalized(std::move(grown));
                    options.push_back({-component_score(grown, pool.t1, pool.t2), x});
                }
                std::sort(options.begin(), options.end());
                if (static_cast<int>(options.size()) > option_limit) {
                    options.resize(option_limit);
                }
                for (auto [score, x] : options) {
                    (void)score;
                    std::vector<Leaf> grown = state;
                    grown.push_back(x);
                    grown = normalized(std::move(grown));
                    std::string key = key_of(grown);
                    if (!local_seen.insert(key).second) continue;
                    if (global_seen.insert(key).second && pool.add(grown)) {
                        ++accepted;
                        next_states.push_back(std::move(grown));
                    }
                    if (accepted >= add_limit ||
                        static_cast<int>(pool.columns.size()) >= pool.max_pool) {
                        break;
                    }
                }
            }
            std::sort(next_states.begin(), next_states.end(),
                      [&](const auto& a, const auto& b) {
                          int sa = component_score(a, pool.t1, pool.t2);
                          int sb = component_score(b, pool.t1, pool.t2);
                          if (sa != sb) return sa > sb;
                          if (a.size() != b.size()) return a.size() > b.size();
                          return a < b;
                      });
            next_states.erase(std::unique(next_states.begin(), next_states.end()),
                              next_states.end());
            if (static_cast<int>(next_states.size()) > beam_width) {
                next_states.resize(beam_width);
            }
            beam = std::move(next_states);
        }
    }
}

void add_edit_pool(Pool& pool, std::chrono::steady_clock::time_point start,
                   double deadline) {
    const int n = pool.t1.n;
    const int distance = env_int_hp("RS_HIGHS_EDIT_DISTANCE", 22);
    const int min_base = env_int_hp("RS_HIGHS_EDIT_MIN_BASE_SIZE", 4);
    const int max_base = env_int_hp("RS_HIGHS_EDIT_MAX_BASE_SIZE", 13);
    const int max_remove = env_int_hp("RS_HIGHS_EDIT_MAX_REMOVE", 3);
    const int option_limit = env_int_hp("RS_HIGHS_EDIT_OPTION_LIMIT", 80);
    const int base_limit = env_int_hp("RS_HIGHS_EDIT_BASE_LIMIT", 4400);
    const int add_limit = env_int_hp("RS_HIGHS_EDIT_ADD_LIMIT", 2200);

    auto adjacency = build_adjacency(pool.t1, pool.t2, distance);
    std::vector<std::vector<Leaf>> bases;
    for (const Column& col : pool.columns) {
        if (static_cast<int>(col.leaves.size()) >= min_base &&
            static_cast<int>(col.leaves.size()) <= max_base) {
            bases.push_back(col.leaves);
        }
    }
    std::sort(bases.begin(), bases.end(), [&](const auto& a, const auto& b) {
        int da = std::abs(static_cast<int>(a.size()) - 10);
        int db = std::abs(static_cast<int>(b.size()) - 10);
        if (da != db) return da < db;
        return component_score(a, pool.t1, pool.t2) >
               component_score(b, pool.t1, pool.t2);
    });
    bases.erase(std::unique(bases.begin(), bases.end()), bases.end());
    if (static_cast<int>(bases.size()) > base_limit) bases.resize(base_limit);

    int accepted = 0;
    std::unordered_set<std::string> generated;
    for (const auto& base : bases) {
        if (accepted >= add_limit || static_cast<int>(pool.columns.size()) >= pool.max_pool ||
            elapsed_seconds(start) >= deadline) {
            break;
        }
        std::vector<char> in_base(n + 1, 0);
        for (Leaf x : base) in_base[x] = 1;
        std::vector<std::pair<int, Leaf>> options;
        for (Leaf x = 1; x <= n; ++x) {
            if (in_base[x]) continue;
            int conflicts = 0;
            int total_distance = 0;
            for (Leaf y : base) {
                int d = std::max(leaf_distance(pool.t1, x, y),
                                 leaf_distance(pool.t2, x, y));
                total_distance += d;
                if (d > distance) ++conflicts;
            }
            if (conflicts <= max_remove) options.push_back({total_distance, x});
        }
        std::sort(options.begin(), options.end());
        if (static_cast<int>(options.size()) > option_limit) {
            options.resize(option_limit);
        }
        for (int remove_count = 1; remove_count <= max_remove; ++remove_count) {
            for (int ri = 0; ri < static_cast<int>(base.size()); ++ri) {
                std::vector<Leaf> kept;
                for (int bi = 0; bi < static_cast<int>(base.size()); ++bi) {
                    if (bi != ri) kept.push_back(base[bi]);
                }
                for (auto [unused, x] : options) {
                    (void)unused;
                    if (!adjacent_to_all(adjacency, x, kept)) continue;
                    std::vector<Leaf> candidate = kept;
                    candidate.push_back(x);
                    candidate = normalized(std::move(candidate));
                    std::string key = key_of(candidate);
                    if (!generated.insert(key).second) continue;
                    if (pool.add(candidate)) ++accepted;
                    if (accepted >= add_limit ||
                        static_cast<int>(pool.columns.size()) >= pool.max_pool ||
                        elapsed_seconds(start) >= deadline) {
                        return;
                    }
                }
                if (remove_count == 1) break;
            }
        }
    }
}

bool solve_pool(Pool& pool, int incumbent, double seconds,
                Forest& improved) {
    if (pool.columns.empty() || incumbent <= 1) return false;
    const int leaf_count = pool.t1.n;
    std::vector<std::vector<int>> by_leaf(leaf_count + 1);
    std::vector<std::vector<int>> by_edge1(pool.t1.node_count());
    std::vector<std::vector<int>> by_edge2(pool.t2.node_count());
    for (int i = 0; i < static_cast<int>(pool.columns.size()); ++i) {
        for (Leaf x : pool.columns[i].leaves) by_leaf[x].push_back(i);
        for (NodeId e : pool.columns[i].edges1) by_edge1[e].push_back(i);
        for (NodeId e : pool.columns[i].edges2) by_edge2[e].push_back(i);
    }

    std::string model_path = "/tmp/rs_highs_model_" +
                             std::to_string(static_cast<long long>(getpid())) +
                             ".lp";
    std::ofstream model(model_path);
    if (!model) return false;
    auto write_terms = [&](const std::vector<int>& cols) {
        if (cols.empty()) {
            model << " 0\n";
            return;
        }
        int width = 1;
        model << ' ';
        for (int col : cols) {
            std::string term = "+ x" + std::to_string(col) + " ";
            if (width + static_cast<int>(term.size()) > 220) {
                model << "\n  ";
                width = 2;
            }
            model << term;
            width += static_cast<int>(term.size());
        }
        model << '\n';
    };

    std::vector<int> all(pool.columns.size());
    std::iota(all.begin(), all.end(), 0);
    model << "Minimize\n obj:\n";
    write_terms(all);
    model << "Subject To\n";
    for (Leaf x = 1; x <= leaf_count; ++x) {
        model << " leaf_" << x << ":\n";
        write_terms(by_leaf[x]);
        model << " = 1\n";
    }
    for (int e = 0; e < static_cast<int>(by_edge1.size()); ++e) {
        if (by_edge1[e].size() <= 1) continue;
        model << " tree1_edge_" << e << ":\n";
        write_terms(by_edge1[e]);
        model << " <= 1\n";
    }
    for (int e = 0; e < static_cast<int>(by_edge2.size()); ++e) {
        if (by_edge2[e].size() <= 1) continue;
        model << " tree2_edge_" << e << ":\n";
        write_terms(by_edge2[e]);
        model << " <= 1\n";
    }
    model << " better_than_incumbent:\n";
    write_terms(all);
    model << " <= " << (incumbent - 1) << "\n";
    model << "Binary\n";
    for (int i = 0; i < static_cast<int>(pool.columns.size()); ++i) {
        model << " x" << i << '\n';
    }
    model << "End\n";
    model.close();

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("time_limit", std::max(1.0, seconds));
    highs.setOptionValue("mip_rel_gap", 0.0);
    highs.setOptionValue("random_seed", 1);
    HighsStatus read_status = highs.readModel(model_path);
    std::remove(model_path.c_str());
    if (read_status == HighsStatus::kError) return false;
    HighsStatus run_status = highs.run();
    if (run_status == HighsStatus::kError) return false;

    std::cerr << "rs_highs_status="
              << highs.modelStatusToString(highs.getModelStatus())
              << " objective=" << highs.getInfo().objective_function_value
              << '\n';
    const HighsSolution& solution = highs.getSolution();
    if (!solution.value_valid ||
        solution.col_value.size() != pool.columns.size()) {
        highs.resetGlobalScheduler(true);
        return false;
    }
    std::vector<int> selected;
    for (size_t i = 0; i < solution.col_value.size(); ++i) {
        if (solution.col_value[i] > 0.5) selected.push_back(static_cast<int>(i));
    }
    if (selected.empty() || static_cast<int>(selected.size()) >= incumbent) {
        highs.resetGlobalScheduler(true);
        return false;
    }
    Forest candidate =
        forest_from_columns(pool.columns, selected, leaf_count, pool.t1, pool.t2);
    ValidationResult check = validate_forest(candidate, pool.t1, pool.t2);
    if (!check.ok || candidate.active_count != static_cast<int>(selected.size()) ||
        candidate.active_count >= incumbent) {
        highs.resetGlobalScheduler(true);
        return false;
    }
    improved = std::move(candidate);
    highs.resetGlobalScheduler(true);
    return true;
}

void add_forest_to_pool(Pool& pool, const Forest& forest) {
    for (const Component& c : forest.comps) {
        if (c.active) pool.add(c.leaves);
    }
}

}  // namespace

Forest run_highs_pool_repair(
    Forest best, const Instance& instance,
    std::chrono::steady_clock::time_point program_start) {
    if (env_int_hp("RS_HIGHS_POOL_REPAIR", 1) == 0) return best;
    if (instance.n > env_int_hp("RS_HIGHS_MAX_LEAVES", 500)) return best;

    const double total_seconds = env_double_hp("RS_HIGHS_TOTAL_SECONDS", 1700.0);
    const double guard_seconds = env_double_hp("RS_HIGHS_GUARD_SECONDS", 90.0);
    const int max_passes = env_int_hp("RS_HIGHS_ILP_PASSES", 2);
    const int max_pool = env_int_hp("RS_HIGHS_MAX_POOL", 12000);
    const int max_component_leaves =
        env_int_hp("RS_HIGHS_MAX_COMPONENT_LEAVES", instance.n);
    auto remaining = [&]() {
        return total_seconds - guard_seconds - elapsed_seconds(program_start);
    };
    if (remaining() < env_double_hp("RS_HIGHS_MIN_SECONDS", 30.0)) return best;

    Pool pool(instance.t1, instance.t2, max_pool, max_component_leaves);
    for (Leaf x = 1; x <= instance.n; ++x) pool.add({x});
    add_forest_to_pool(pool, best);

    for (int pass = 0; pass < max_passes && remaining() > 20.0; ++pass) {
        int before_pool = static_cast<int>(pool.columns.size());
        double pool_budget = std::min(
            remaining() * 0.35,
            env_double_hp(pass == 0 ? "RS_HIGHS_POOL_SECONDS"
                                    : "RS_HIGHS_REENRICH_SECONDS",
                          pass == 0 ? 300.0 : 120.0));
        double pool_deadline = elapsed_seconds(program_start) + pool_budget;
        std::vector<Region> regions =
            collect_regions(best, instance.t1, instance.t2);
        add_region_pool(pool, best, regions, program_start, pool_deadline);
        if (elapsed_seconds(program_start) < pool_deadline &&
            static_cast<int>(pool.columns.size()) < max_pool) {
            add_bucket_clique_pool(pool, program_start, pool_deadline);
        }
        if (elapsed_seconds(program_start) < pool_deadline &&
            static_cast<int>(pool.columns.size()) < max_pool) {
            add_edit_pool(pool, program_start, pool_deadline);
        }
        std::cerr << "rs_highs_pool pass=" << pass
                  << " before=" << before_pool
                  << " after=" << pool.columns.size()
                  << " regions=" << regions.size()
                  << " best=" << best.active_count
                  << " remaining=" << remaining() << '\n';

        double highs_seconds = std::min(
            remaining(),
            env_double_hp("RS_HIGHS_SECONDS", pass == 0 ? 900.0 : 360.0));
        if (highs_seconds < env_double_hp("RS_HIGHS_MIN_SOLVE_SECONDS", 20.0)) {
            break;
        }
        Forest improved;
        if (!solve_pool(pool, best.active_count, highs_seconds, improved)) {
            std::cerr << "rs_highs_no_improvement pass=" << pass
                      << " best=" << best.active_count << '\n';
            break;
        }
        std::cerr << "rs_highs_improvement pass=" << pass
                  << " " << best.active_count << " -> "
                  << improved.active_count << '\n';
        best = std::move(improved);
        publish_best_output(forest_to_output(best, instance.t1));
        add_forest_to_pool(pool, best);
    }
    return best;
}
