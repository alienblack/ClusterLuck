#pragma once

#include <string>
#include <vector>

#include "tree/tree.h"

struct InducedInfo {
    NodeId root = -1;
    std::string signature;
    std::string newick;
    std::vector<NodeId> used_edges;
};

InducedInfo compute_induced(const Tree& tree, const std::vector<Leaf>& leaves);

