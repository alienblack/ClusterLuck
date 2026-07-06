#pragma once

#include <chrono>

#include "forest/forest.h"
#include "io/pace_io.h"

Forest run_highs_pool_repair(
    Forest best, const Instance& instance,
    std::chrono::steady_clock::time_point program_start);
