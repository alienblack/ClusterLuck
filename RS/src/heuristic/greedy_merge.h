#pragma once

#include "forest/forest.h"

struct GreedyConfig {
    int all_pair_cap = 1800;
    int repair_top_components = 18;
    int repair_pair_top_components = 10;
    int repair_pair_limit = 24;
    int conflict_target_limit = 36;
    int conflict_max_split_components = 4;
    int repair_passes = 2;
    bool publish_improvements = true;
    bool enable_distance_routes = false;
    bool enable_common_cherry_contraction = false;
};

Forest run_greedy_merge(Forest forest, const Tree& t1, const Tree& t2,
                        const GreedyConfig& config);

Forest run_h45_greedy_merge(Forest forest, const Tree& t1, const Tree& t2,
                            int score_mode, bool use_lookahead,
                            bool reconsider_failed);

Forest exact_output_repair(Forest forest, const Tree& t1, const Tree& t2);

Forest run_h45_split_restarts(Forest initial, const Tree& t1, const Tree& t2,
                              int score_mode, int attempts,
                              int seed_offset, int extra_splits,
                              bool use_fpt_conflict = false,
                              int target_count = 0,
                              double max_seconds = 0.0);

Forest run_h45_blocker_repair(Forest initial, const Tree& t1, const Tree& t2,
                              int max_components, int max_blockers,
                              int candidate_limit, int restart_budget,
                              int target_count = 0,
                              double max_seconds = 0.0);

Forest run_h45_window_repair(Forest initial, const Tree& t1, const Tree& t2,
                             int max_components, int candidate_limit,
                             int window_limit, int restart_budget,
                             int target_count = 0,
                             double max_seconds = 0.0);

Forest run_h45_local_exchange_repair(Forest initial, const Tree& t1,
                                     const Tree& t2, int max_components,
                                     int candidate_limit, int window_limit,
                                     int max_blockers, int restart_budget,
                                     int target_count = 0,
                                     double max_seconds = 0.0);

Forest run_h45_cross_recombine_repair(Forest initial, const Tree& t1,
                                      const Tree& t2, int max_components,
                                      int candidate_limit, int window_limit,
                                      int restart_budget,
                                      int target_count = 0,
                                      double max_seconds = 0.0);

Forest run_h45_targeted_split_repack(Forest initial, const Tree& t1,
                                     const Tree& t2, int max_components,
                                     int candidate_limit, int group_width,
                                     int restart_budget, int passes = 1,
                                     int target_count = 0,
                                     double max_seconds = 0.0);

Forest run_split_restarts(Forest initial, const Tree& t1, const Tree& t2,
                          const GreedyConfig& config, int attempts,
                          int seed_offset, int extra_splits);
