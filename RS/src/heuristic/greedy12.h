#pragma once

#include "cluster/cluster_points.h"
#include "forest/forest.h"
#include "tree/tree.h"

bool should_run_greedy12_tail(int leaves, int residual_leaves);
bool run_greedy12_tail(const Tree& t1, const Tree& t2, Forest& out);

