#pragma once

#include <string>
#include <vector>

#include "tree/tree.h"

struct Component {
    bool active = true;
    std::vector<Leaf> leaves;
    Leaf repr = -1;
    NodeId root1 = -1;
    NodeId root2 = -1;
    std::string signature1;
    std::string signature2;
    std::vector<NodeId> edges1;
    std::vector<NodeId> edges2;
};

struct Forest {
    std::vector<Component> comps;
    std::vector<CompId> leaf_comp;
    std::vector<CompId> owner1;
    std::vector<CompId> owner2;
    int active_count = 0;
};

Forest make_singleton_forest(int n);
bool rebuild_forest_data(Forest& forest, const Tree& t1, const Tree& t2);
std::string forest_to_output(const Forest& forest, const Tree& t1);
std::vector<CompId> active_components(const Forest& forest);
void merge_components(Forest& forest, CompId a, CompId b, std::vector<Leaf> merged_leaves);
