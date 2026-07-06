#pragma once

#include <string>
#include <vector>

#include "types.h"

struct Tree {
    int n = 0;
    NodeId root = -1;

    std::vector<NodeId> parent;
    std::vector<NodeId> left;
    std::vector<NodeId> right;
    std::vector<Leaf> label;
    std::vector<NodeId> leaf_node;

    std::vector<int> depth;
    std::vector<int> tin;
    std::vector<int> tout;
    std::vector<std::vector<NodeId>> up;

    NodeId add_node(Leaf leaf_label);
    bool is_leaf(NodeId v) const;
    int node_count() const;
    void prepare(int leaf_count);
    bool is_ancestor(NodeId a, NodeId b) const;
    NodeId lca(NodeId a, NodeId b) const;
    NodeId lca_leaves(const std::vector<Leaf>& leaves) const;

private:
    int timer_ = 0;
    void dfs_prepare(NodeId v, NodeId p);
};
