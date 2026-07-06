#pragma once

#include "forest/forest.h"

Forest run_near_cluster_route(Forest initial, const Tree& t1, const Tree& t2,
                              int max_leaves, int min_island_size,
                              int max_spill, int max_trim,
                              int max_candidates, int max_islands,
                              int restart_budget, int target_count = 0);

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
                                       int target_count = 0,
                                       double max_seconds = 0.0);
