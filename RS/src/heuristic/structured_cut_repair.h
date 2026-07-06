#pragma once

#include "forest/forest.h"

Forest run_structured_cut_repair(Forest initial, const Tree& t1,
                                 const Tree& t2, int max_leaves,
                                 int max_components, int max_candidates,
                                 int min_component_leaves,
                                 int min_block_leaves,
                                 int max_block_leaves,
                                 int repair_restart_budget,
                                 int target_count,
                                 double max_seconds);
