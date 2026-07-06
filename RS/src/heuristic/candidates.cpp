#include "heuristic/candidates.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace {

using Pair = std::pair<CompId, CompId>;

Pair ordered_pair(CompId a, CompId b) {
    if (b < a) std::swap(a, b);
    return {a, b};
}

void add_score(std::map<Pair, int>& scores, CompId a, CompId b, int score) {
    if (a == b || a < 0 || b < 0) return;
    scores[ordered_pair(a, b)] += score;
}

CompId subtree_component(const Tree& tree, const Forest& forest, NodeId v,
                         std::map<Pair, int>& scores) {
    if (tree.is_leaf(v)) return forest.leaf_comp[tree.label[v]];

    CompId a = subtree_component(tree, forest, tree.left[v], scores);
    CompId b = subtree_component(tree, forest, tree.right[v], scores);
    if (a == b) return a;

    if (a >= 0 && b >= 0 && forest.comps[a].active && forest.comps[b].active) {
        add_score(scores, a, b, 5000);
    }
    return -2;
}

NodeId comp_root(const Tree& tree, const Component& component) {
    return tree.lca_leaves(component.leaves);
}

void add_window_scores(const Tree& tree, const Forest& forest, int base_score,
                       int levels, int max_bucket_size,
                       std::map<Pair, int>& scores) {
    std::map<NodeId, std::vector<CompId>> buckets;
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        if (!forest.comps[id].active) continue;
        NodeId v = comp_root(tree, forest.comps[id]);
        for (int level = 0; level < levels && v != -1; ++level) {
            buckets[v].push_back(id);
            v = tree.parent[v];
        }
    }

    for (const auto& item : buckets) {
        const std::vector<CompId>& ids = item.second;
        if (ids.size() < 2 || ids.size() > static_cast<size_t>(max_bucket_size)) continue;
        int score = base_score - static_cast<int>(ids.size());
        for (size_t i = 0; i < ids.size(); ++i) {
            for (size_t j = i + 1; j < ids.size(); ++j) {
                add_score(scores, ids[i], ids[j], score);
            }
        }
    }
}

void collect_local_components(const Tree& tree, const Forest& forest, NodeId root,
                              int limit, std::vector<CompId>& out) {
    out.clear();
    std::vector<NodeId> stack{root};
    std::vector<char> seen(forest.comps.size(), 0);
    while (!stack.empty() && static_cast<int>(out.size()) < limit) {
        NodeId v = stack.back();
        stack.pop_back();
        if (tree.is_leaf(v)) {
            CompId id = forest.leaf_comp[tree.label[v]];
            if (id >= 0 && forest.comps[id].active && !seen[id]) {
                seen[id] = 1;
                out.push_back(id);
            }
        } else {
            stack.push_back(tree.right[v]);
            stack.push_back(tree.left[v]);
        }
    }
}

int merged_size(const Forest& forest, const Pair& p) {
    return static_cast<int>(forest.comps[p.first].leaves.size() +
                            forest.comps[p.second].leaves.size());
}

uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

int pair_jitter(CompId a, CompId b, int seed, int span) {
    if (span <= 0) return 0;
    uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) ^
                   static_cast<uint32_t>(b) ^
                   (static_cast<uint64_t>(static_cast<uint32_t>(seed)) << 1);
    return static_cast<int>(splitmix64(key) % static_cast<uint64_t>(span));
}

int component_distance_score(const Forest& forest, const Tree& t1, const Tree& t2,
                             const Pair& p) {
    const Component& a = forest.comps[p.first];
    const Component& b = forest.comps[p.second];
    if (a.root1 < 0 || a.root2 < 0 || b.root1 < 0 || b.root2 < 0) return 0;
    NodeId l1 = t1.lca(a.root1, b.root1);
    NodeId l2 = t2.lca(a.root2, b.root2);
    int d1 = t1.depth[a.root1] + t1.depth[b.root1] - 2 * t1.depth[l1];
    int d2 = t2.depth[a.root2] + t2.depth[b.root2] - 2 * t2.depth[l2];
    int max_distance = std::max(d1, d2);
    int total_distance = d1 + d2;
    int combined_size = std::min(255, merged_size(forest, p));
    return 1000000 - 100 * (32 * max_distance + total_distance) + combined_size;
}

}  // namespace

const char* candidate_mode_name(CandidateMode mode) {
    switch (mode) {
        case CandidateMode::SiblingAllPairs:
            return "sibling_all_pairs";
        case CandidateMode::WindowScored:
            return "window_scored";
        case CandidateMode::WindowTight:
            return "window_tight";
        case CandidateMode::WindowWide:
            return "window_wide";
        case CandidateMode::DistanceScored:
            return "distance_scored";
        case CandidateMode::LocalDistance:
            return "local_distance";
    }
    return "unknown";
}

std::vector<MergeCandidate> generate_candidates(const Forest& forest, const Tree& t1,
                                                const Tree& t2, int all_pair_cap,
                                                CandidateMode mode, int jitter_seed,
                                                int jitter_span) {
    std::map<Pair, int> scores;

    if (mode == CandidateMode::DistanceScored || mode == CandidateMode::LocalDistance) {
        auto ids = active_components(forest);
        if (mode == CandidateMode::DistanceScored) {
            for (size_t i = 0; i < ids.size(); ++i) {
                for (size_t j = i + 1; j < ids.size(); ++j) {
                    Pair p = ordered_pair(ids[i], ids[j]);
                    scores[p] = component_distance_score(forest, t1, t2, p);
                }
            }
        } else {
            std::vector<CompId> local;
            for (const Tree* tree : {&t1, &t2}) {
                for (NodeId v = 0; v < tree->node_count(); ++v) {
                    if (tree->is_leaf(v)) continue;
                    collect_local_components(*tree, forest, v, 24, local);
                    for (size_t i = 0; i < local.size(); ++i) {
                        for (size_t j = i + 1; j < local.size(); ++j) {
                            Pair p = ordered_pair(local[i], local[j]);
                            scores[p] = std::max(scores[p],
                                                 component_distance_score(forest, t1, t2, p));
                        }
                    }
                }
            }
        }
        std::vector<MergeCandidate> candidates;
        candidates.reserve(scores.size());
        for (const auto& item : scores) {
            int jitter = pair_jitter(item.first.first, item.first.second,
                                     jitter_seed, jitter_span);
            candidates.push_back({item.first.first, item.first.second,
                                  item.second + jitter});
        }
        std::sort(candidates.begin(), candidates.end(), [&](const auto& x, const auto& y) {
            if (x.score != y.score) return x.score > y.score;
            Pair px{x.a, x.b};
            Pair py{y.a, y.b};
            int sx = merged_size(forest, px);
            int sy = merged_size(forest, py);
            if (sx != sy) return sx > sy;
            return px < py;
        });
        return candidates;
    }

    subtree_component(t1, forest, t1.root, scores);
    subtree_component(t2, forest, t2.root, scores);
    if (mode == CandidateMode::WindowScored) {
        add_window_scores(t1, forest, 900, 5, 80, scores);
        add_window_scores(t2, forest, 900, 5, 80, scores);
    } else if (mode == CandidateMode::WindowTight) {
        add_window_scores(t1, forest, 1300, 3, 40, scores);
        add_window_scores(t2, forest, 1300, 3, 40, scores);
    } else if (mode == CandidateMode::WindowWide) {
        add_window_scores(t1, forest, 650, 8, 140, scores);
        add_window_scores(t2, forest, 650, 8, 140, scores);
    }

    if (forest.active_count <= all_pair_cap) {
        auto ids = active_components(forest);
        for (size_t i = 0; i < ids.size(); ++i) {
            for (size_t j = i + 1; j < ids.size(); ++j) {
                add_score(scores, ids[i], ids[j], 1);
            }
        }
    }

    std::vector<MergeCandidate> candidates;
    candidates.reserve(scores.size());
    for (const auto& item : scores) {
        int jitter = pair_jitter(item.first.first, item.first.second, jitter_seed, jitter_span);
        candidates.push_back({item.first.first, item.first.second, item.second + jitter});
    }
    std::sort(candidates.begin(), candidates.end(), [&](const auto& x, const auto& y) {
        if (x.score != y.score) return x.score > y.score;
        Pair px{x.a, x.b};
        Pair py{y.a, y.b};
        if (mode == CandidateMode::SiblingAllPairs) return px < py;
        int sx = merged_size(forest, px);
        int sy = merged_size(forest, py);
        if (sx != sy) return sx > sy;
        return px < py;
    });
    return candidates;
}
