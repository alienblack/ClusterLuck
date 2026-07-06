#include "heuristic/structured_cut_repair.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_set>
#include <vector>

#include "forest/validator.h"
#include "heuristic/greedy_merge.h"
#include "tree/induced.h"
#include "util/signal_handler.h"

namespace {

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

struct ClockLimit {
    std::chrono::steady_clock::time_point start;
    double max_seconds = 0.0;

    bool ok() const {
        if (max_seconds <= 0.0) return true;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        return elapsed < max_seconds;
    }

    double remaining() const {
        if (max_seconds <= 0.0) return 0.0;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        return std::max(0.1, max_seconds - elapsed);
    }
};

struct CutCandidate {
    CompId component = -1;
    int guide = 0;
    NodeId node = -1;
    int comp_size = 0;
    int block_size = 0;
    int remainder_size = 0;
    int depth = 0;
    int score = 0;
    std::vector<Leaf> block;
    std::vector<Leaf> remainder;
    bool remainder_exact = false;
};

std::vector<std::vector<Leaf>> leaves_below(const Tree& tree) {
    std::vector<std::vector<Leaf>> below(tree.node_count());
    std::vector<NodeId> order;
    order.reserve(tree.node_count());
    if (tree.root >= 0) order.push_back(tree.root);
    for (size_t i = 0; i < order.size(); ++i) {
        NodeId v = order[i];
        if (tree.left[v] >= 0) order.push_back(tree.left[v]);
        if (tree.right[v] >= 0) order.push_back(tree.right[v]);
    }
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        NodeId v = *it;
        if (tree.is_leaf(v)) {
            if (tree.label[v] > 0) below[v].push_back(tree.label[v]);
        } else {
            if (tree.left[v] >= 0) {
                below[v].insert(below[v].end(), below[tree.left[v]].begin(),
                                below[tree.left[v]].end());
            }
            if (tree.right[v] >= 0) {
                below[v].insert(below[v].end(), below[tree.right[v]].begin(),
                                below[tree.right[v]].end());
            }
            std::sort(below[v].begin(), below[v].end());
        }
    }
    return below;
}

std::vector<Leaf> intersect_component(const std::vector<Leaf>& subtree,
                                      const std::vector<unsigned char>& in_comp) {
    std::vector<Leaf> out;
    for (Leaf x : subtree) {
        if (x > 0 && x < static_cast<int>(in_comp.size()) && in_comp[x]) {
            out.push_back(x);
        }
    }
    return out;
}

std::vector<Leaf> subtract_sorted(const std::vector<Leaf>& all,
                                  const std::vector<Leaf>& remove) {
    std::vector<Leaf> out;
    std::set_difference(all.begin(), all.end(), remove.begin(), remove.end(),
                        std::back_inserter(out));
    return out;
}

bool exact_component_valid(const std::vector<Leaf>& leaves, const Tree& t1,
                           const Tree& t2) {
    if (leaves.empty()) return false;
    if (leaves.size() == 1) return true;
    InducedInfo i1 = compute_induced(t1, leaves);
    InducedInfo i2 = compute_induced(t2, leaves);
    return i1.signature == i2.signature;
}

std::string candidate_key(CompId component, const std::vector<Leaf>& block) {
    std::string key = std::to_string(component);
    key.push_back(':');
    for (Leaf x : block) {
        key += std::to_string(x);
        key.push_back(',');
    }
    return key;
}

void add_candidates_from_guide(
    const Forest& forest, const Tree& guide, const Tree& t1, const Tree& t2,
    const std::vector<std::vector<Leaf>>& guide_leaves, CompId cid, int guide_id,
    int min_block_leaves, int max_block_leaves, int per_component_limit,
    std::unordered_set<std::string>& seen, std::vector<CutCandidate>& out) {
    const Component& comp = forest.comps[cid];
    const int comp_size = static_cast<int>(comp.leaves.size());
    if (comp_size <= min_block_leaves) return;

    std::vector<unsigned char> in_comp(guide.n + 1, 0);
    for (Leaf x : comp.leaves) {
        if (x > 0 && x <= guide.n) in_comp[x] = 1;
    }

    std::vector<CutCandidate> local;
    local.reserve(64);
    const int max_fraction_block =
        std::max(min_block_leaves, (comp_size * 4) / 5);
    const int effective_max_block =
        max_block_leaves > 0 ? std::min(max_block_leaves, max_fraction_block)
                             : max_fraction_block;

    for (NodeId node = 0; node < guide.node_count(); ++node) {
        if (guide.is_leaf(node) || node == guide.root) continue;
        std::vector<Leaf> block = intersect_component(guide_leaves[node], in_comp);
        int block_size = static_cast<int>(block.size());
        if (block_size < min_block_leaves || block_size >= comp_size ||
            block_size > effective_max_block) {
            continue;
        }
        std::sort(block.begin(), block.end());
        block.erase(std::unique(block.begin(), block.end()), block.end());
        std::string key = candidate_key(cid, block);
        if (!seen.insert(key).second) continue;
        if (!exact_component_valid(block, t1, t2)) continue;

        std::vector<Leaf> remainder = subtract_sorted(comp.leaves, block);
        if (remainder.empty()) continue;
        bool remainder_exact = exact_component_valid(remainder, t1, t2);
        if (!remainder_exact &&
            env_int("RS_STRUCTURED_CUT_LEFTOVER_SINGLETONS", 1) == 0) {
            continue;
        }

        int smaller = std::min(block_size, static_cast<int>(remainder.size()));
        int larger = std::max(block_size, static_cast<int>(remainder.size()));
        CutCandidate cand;
        cand.component = cid;
        cand.guide = guide_id;
        cand.node = node;
        cand.comp_size = comp_size;
        cand.block_size = block_size;
        cand.remainder_size = static_cast<int>(remainder.size());
        cand.depth = guide.depth[node];
        cand.remainder_exact = remainder_exact;
        cand.score = 100000 * smaller / std::max(1, larger) +
                     100 * block_size + cand.depth +
                     (remainder_exact ? 25000 : 0);
        cand.block = std::move(block);
        cand.remainder = std::move(remainder);
        local.push_back(std::move(cand));
    }

    std::sort(local.begin(), local.end(), [](const CutCandidate& a,
                                             const CutCandidate& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.block_size != b.block_size) return a.block_size > b.block_size;
        return a.depth > b.depth;
    });
    if (per_component_limit > 0 &&
        static_cast<int>(local.size()) > per_component_limit) {
        local.resize(per_component_limit);
    }
    out.insert(out.end(), std::make_move_iterator(local.begin()),
               std::make_move_iterator(local.end()));
}

std::vector<CutCandidate> generate_cut_candidates(
    const Forest& forest, const Tree& t1, const Tree& t2, int max_components,
    int max_candidates, int min_component_leaves, int min_block_leaves,
    int max_block_leaves) {
    std::vector<CompId> ids = active_components(forest);
    std::sort(ids.begin(), ids.end(), [&](CompId a, CompId b) {
        return forest.comps[a].leaves.size() > forest.comps[b].leaves.size();
    });
    if (max_components > 0 && static_cast<int>(ids.size()) > max_components) {
        ids.resize(max_components);
    }

    std::vector<std::vector<Leaf>> below1 = leaves_below(t1);
    std::vector<std::vector<Leaf>> below2 = leaves_below(t2);
    std::unordered_set<std::string> seen;
    std::vector<CutCandidate> candidates;
    int per_component_limit = env_int("RS_STRUCTURED_CUT_PER_COMPONENT", 12);

    for (CompId cid : ids) {
        const Component& comp = forest.comps[cid];
        if (!comp.active ||
            static_cast<int>(comp.leaves.size()) < min_component_leaves) {
            continue;
        }
        add_candidates_from_guide(forest, t1, t1, t2, below1, cid, 1,
                                  min_block_leaves, max_block_leaves,
                                  per_component_limit, seen, candidates);
        add_candidates_from_guide(forest, t2, t1, t2, below2, cid, 2,
                                  min_block_leaves, max_block_leaves,
                                  per_component_limit, seen, candidates);
    }

    std::sort(candidates.begin(), candidates.end(), [](const CutCandidate& a,
                                                       const CutCandidate& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.comp_size != b.comp_size) return a.comp_size > b.comp_size;
        if (a.block_size != b.block_size) return a.block_size > b.block_size;
        return a.node < b.node;
    });
    if (max_candidates > 0 &&
        static_cast<int>(candidates.size()) > max_candidates) {
        candidates.resize(max_candidates);
    }
    return candidates;
}

void append_component(Forest& forest, const std::vector<Leaf>& leaves) {
    Component comp;
    comp.active = true;
    comp.leaves = leaves;
    comp.repr = leaves.empty() ? -1 : leaves.front();
    forest.comps.push_back(std::move(comp));
    ++forest.active_count;
}

bool build_cut_trial(const Forest& base, const CutCandidate& cand,
                     const Tree& t1, const Tree& t2, Forest& trial) {
    trial = Forest{};
    trial.leaf_comp.assign(t1.n + 1, -1);
    trial.comps.reserve(base.active_count + cand.remainder_size + 2);

    for (CompId id = 0; id < static_cast<CompId>(base.comps.size()); ++id) {
        const Component& comp = base.comps[id];
        if (!comp.active || id == cand.component) continue;
        append_component(trial, comp.leaves);
    }
    append_component(trial, cand.block);
    if (cand.remainder_exact) {
        append_component(trial, cand.remainder);
    } else {
        for (Leaf x : cand.remainder) {
            append_component(trial, std::vector<Leaf>{x});
        }
    }
    return rebuild_forest_data(trial, t1, t2);
}

bool build_cut_trial_pair(const Forest& base, const CutCandidate& a,
                          const CutCandidate& b, const Tree& t1,
                          const Tree& t2, Forest& trial) {
    if (a.component == b.component) return false;

    trial = Forest{};
    trial.leaf_comp.assign(t1.n + 1, -1);
    trial.comps.reserve(base.active_count + a.remainder_size + b.remainder_size + 4);

    for (CompId id = 0; id < static_cast<CompId>(base.comps.size()); ++id) {
        const Component& comp = base.comps[id];
        if (!comp.active || id == a.component || id == b.component) continue;
        append_component(trial, comp.leaves);
    }

    auto append_split = [&](const CutCandidate& cand) {
        append_component(trial, cand.block);
        if (cand.remainder_exact) {
            append_component(trial, cand.remainder);
        } else {
            for (Leaf x : cand.remainder) {
                append_component(trial, std::vector<Leaf>{x});
            }
        }
    };
    append_split(a);
    append_split(b);
    return rebuild_forest_data(trial, t1, t2);
}

Forest repair_trial(Forest trial, const Tree& t1, const Tree& t2, int mode,
                    int restart_budget, int seed, int target_count,
                    const ClockLimit& limit) {
    trial = run_h45_greedy_merge(std::move(trial), t1, t2, mode, true, true);
    trial = exact_output_repair(std::move(trial), t1, t2);
    if (!validate_forest(trial, t1, t2).ok) return Forest{};
    if (restart_budget > 0 && limit.ok()) {
        Forest restarted = run_h45_split_restarts(
            trial, t1, t2, mode, restart_budget, seed, 1, true, target_count,
            limit.remaining());
        restarted = exact_output_repair(std::move(restarted), t1, t2);
        if (validate_forest(restarted, t1, t2).ok &&
            restarted.active_count < trial.active_count) {
            trial = std::move(restarted);
        }
    }
    return trial;
}

}  // namespace

Forest run_structured_cut_repair(Forest initial, const Tree& t1,
                                 const Tree& t2, int max_leaves,
                                 int max_components, int max_candidates,
                                 int min_component_leaves,
                                 int min_block_leaves,
                                 int max_block_leaves,
                                 int repair_restart_budget,
                                 int target_count,
                                 double max_seconds) {
    if (t1.n > max_leaves || max_candidates == 0 ||
        min_component_leaves < 2 || min_block_leaves < 2) {
        return initial;
    }

    ClockLimit limit{std::chrono::steady_clock::now(), max_seconds};
    const bool debug = std::getenv("RS_STRUCTURED_CUT_DEBUG") != nullptr;

    if (!rebuild_forest_data(initial, t1, t2)) return initial;
    Forest best = exact_output_repair(std::move(initial), t1, t2);
    if (!validate_forest(best, t1, t2).ok) return best;

    int passes = env_int("RS_STRUCTURED_CUT_PASSES", 2);
    int considered = 0;
    int valid_trials = 0;
    int improved = 0;
    int start_count = best.active_count;

    for (int pass = 0; pass < passes && limit.ok(); ++pass) {
        std::vector<CutCandidate> candidates = generate_cut_candidates(
            best, t1, t2, max_components, max_candidates, min_component_leaves,
            min_block_leaves, max_block_leaves);
        if (debug) {
            std::cerr << "rs_structured_cut_pass pass=" << pass
                      << " components=" << best.active_count
                      << " candidates=" << candidates.size() << '\n';
        }
        bool pass_improved = false;
        if (env_int("RS_STRUCTURED_CUT_SINGLE_ROUTE", 1) != 0) {
            for (const CutCandidate& cand : candidates) {
                if (!limit.ok()) break;
                ++considered;
                Forest trial;
                if (!build_cut_trial(best, cand, t1, t2, trial)) continue;
                ++valid_trials;

                int modes[] = {
                    env_int("RS_STRUCTURED_CUT_MODE1", 14),
                    env_int("RS_STRUCTURED_CUT_MODE2", 2),
                    env_int("RS_STRUCTURED_CUT_MODE3", 20),
                    env_int("RS_STRUCTURED_CUT_MODE4", 1),
                };
                for (int mode : modes) {
                    if (!limit.ok()) break;
                    Forest repaired = repair_trial(
                        trial, t1, t2, mode, repair_restart_budget,
                        91000 + pass * 10007 + considered * 257 + mode * 13,
                        target_count, limit);
                    if (repaired.active_count == 0) continue;
                    if (repaired.active_count < best.active_count) {
                        best = std::move(repaired);
                        publish_best_output(forest_to_output(best, t1));
                        ++improved;
                        pass_improved = true;
                        if (debug) {
                            std::cerr << "rs_structured_cut_improve components="
                                      << best.active_count
                                      << " pass=" << pass
                                      << " guide=" << cand.guide
                                      << " comp_size=" << cand.comp_size
                                      << " block=" << cand.block_size
                                      << " rem=" << cand.remainder_size
                                      << " rem_exact="
                                      << (cand.remainder_exact ? 1 : 0)
                                      << " mode=" << mode << '\n';
                        }
                        if (target_count > 0 && best.active_count <= target_count) {
                            return best;
                        }
                        break;
                    }
                }
                if (pass_improved &&
                    env_int("RS_STRUCTURED_CUT_REGEN_AFTER_IMPROVE", 1) != 0) {
                    break;
                }
            }
        }
        if (!pass_improved &&
            env_int("RS_STRUCTURED_CUT_PAIR_ROUTE", 1) != 0 && limit.ok()) {
            const int pair_limit =
                env_int("RS_STRUCTURED_CUT_PAIR_COMBO_LIMIT", 1200);
            int pair_considered = 0;
            for (int i = 0; i < static_cast<int>(candidates.size()) &&
                            pair_considered < pair_limit && limit.ok(); ++i) {
                for (int j = i + 1; j < static_cast<int>(candidates.size()) &&
                                pair_considered < pair_limit && limit.ok(); ++j) {
                    if (candidates[i].component == candidates[j].component) continue;
                    ++pair_considered;
                    Forest trial;
                    if (!build_cut_trial_pair(best, candidates[i], candidates[j],
                                              t1, t2, trial)) {
                        continue;
                    }
                    ++valid_trials;
                    int modes[] = {
                        env_int("RS_STRUCTURED_CUT_PAIR_MODE1", 2),
                        env_int("RS_STRUCTURED_CUT_PAIR_MODE2", 14),
                        env_int("RS_STRUCTURED_CUT_PAIR_MODE3", 20),
                    };
                    for (int mode : modes) {
                        if (!limit.ok()) break;
                        Forest repaired = repair_trial(
                            trial, t1, t2, mode, repair_restart_budget,
                            191000 + pass * 10007 + pair_considered * 257 +
                                mode * 13,
                            target_count, limit);
                        if (repaired.active_count == 0) continue;
                        if (repaired.active_count < best.active_count) {
                            best = std::move(repaired);
                            publish_best_output(forest_to_output(best, t1));
                            ++improved;
                            pass_improved = true;
                            if (debug) {
                                std::cerr
                                    << "rs_structured_cut_pair_improve components="
                                    << best.active_count
                                    << " pass=" << pass
                                    << " pair_considered=" << pair_considered
                                    << " a_comp_size=" << candidates[i].comp_size
                                    << " a_block=" << candidates[i].block_size
                                    << " b_comp_size=" << candidates[j].comp_size
                                    << " b_block=" << candidates[j].block_size
                                    << " mode=" << mode << '\n';
                            }
                            if (target_count > 0 &&
                                best.active_count <= target_count) {
                                return best;
                            }
                            break;
                        }
                    }
                    if (pass_improved &&
                        env_int("RS_STRUCTURED_CUT_REGEN_AFTER_IMPROVE", 1) != 0) {
                        break;
                    }
                }
            }
        }
        if (!pass_improved &&
            env_int("RS_STRUCTURED_CUT_STOP_ON_STALE_PASS", 1) != 0) {
            break;
        }
    }

    if (debug) {
        std::cerr << "rs_structured_cut_done start=" << start_count
                  << " components=" << best.active_count
                  << " considered=" << considered
                  << " valid_trials=" << valid_trials
                  << " improved=" << improved << '\n';
    }
    return best;
}
