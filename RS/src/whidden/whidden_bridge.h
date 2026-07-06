#pragma once

#include <string>
#include <vector>

bool whidden_seed_components(const std::string& first_newick,
                             const std::string& second_newick,
                             int leaf_count,
                             int greedy_upper_bound,
                             int cluster_tune,
                             int split_threshold,
                             int random_seed,
                             std::vector<std::vector<int>>& components);
