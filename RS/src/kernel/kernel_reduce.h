#pragma once

#include <optional>
#include <unordered_map>
#include <vector>

#include "forest/forest.h"
#include "io/pace_io.h"

struct KernelMap {
    struct ChainRecord {
        std::vector<int> anchor1;
        std::vector<int> anchor2;
        std::vector<int> anchor3;
        std::vector<int> deleted;
    };

    std::unordered_map<int, std::vector<int>> active_payload;
    std::vector<std::vector<int>> forced;
    std::vector<ChainRecord> chains;
};

struct KernelReducedInstance {
    Instance reduced;
    KernelMap map;
    int original_n = 0;
    int subtree_deleted = 0;
    int chain_deleted = 0;
    int param_reduction = 0;
};

std::optional<KernelReducedInstance> kernel_reduce_instance(
    const Instance& instance, bool verbose = false, int protected_taxon = -1);

Forest expand_kernel_reduced_forest(const Forest& reduced_forest,
                                    const KernelMap& map,
                                    const Instance& original);
