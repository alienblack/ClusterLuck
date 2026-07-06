#include "tree/tree.h"

#include <algorithm>
#include <stdexcept>

NodeId Tree::add_node(Leaf leaf_label) {
    NodeId id = static_cast<NodeId>(parent.size());
    parent.push_back(-1);
    left.push_back(-1);
    right.push_back(-1);
    label.push_back(leaf_label);
    return id;
}

bool Tree::is_leaf(NodeId v) const {
    return left[v] == -1 && right[v] == -1;
}

int Tree::node_count() const {
    return static_cast<int>(parent.size());
}

void Tree::prepare(int leaf_count) {
    n = leaf_count;
    leaf_node.assign(n + 1, -1);
    for (NodeId v = 0; v < node_count(); ++v) {
        if (label[v] > 0) {
            if (label[v] > n) {
                throw std::runtime_error("leaf label exceeds declared n");
            }
            if (leaf_node[label[v]] != -1) {
                throw std::runtime_error("duplicate leaf label");
            }
            leaf_node[label[v]] = v;
        }
    }
    for (Leaf x = 1; x <= n; ++x) {
        if (leaf_node[x] == -1) {
            throw std::runtime_error("missing leaf label");
        }
    }

    int log = 1;
    while ((1 << log) <= std::max(1, node_count())) ++log;
    depth.assign(node_count(), 0);
    tin.assign(node_count(), 0);
    tout.assign(node_count(), 0);
    up.assign(log, std::vector<NodeId>(node_count(), -1));
    timer_ = 0;
    dfs_prepare(root, root);
}

void Tree::dfs_prepare(NodeId v, NodeId p) {
    parent[v] = (v == root ? -1 : p);
    tin[v] = timer_++;
    up[0][v] = (v == root ? v : p);
    for (int j = 1; j < static_cast<int>(up.size()); ++j) {
        up[j][v] = up[j - 1][up[j - 1][v]];
    }
    if (left[v] != -1) {
        depth[left[v]] = depth[v] + 1;
        dfs_prepare(left[v], v);
    }
    if (right[v] != -1) {
        depth[right[v]] = depth[v] + 1;
        dfs_prepare(right[v], v);
    }
    tout[v] = timer_;
}

bool Tree::is_ancestor(NodeId a, NodeId b) const {
    return tin[a] <= tin[b] && tout[b] <= tout[a];
}

NodeId Tree::lca(NodeId a, NodeId b) const {
    if (is_ancestor(a, b)) return a;
    if (is_ancestor(b, a)) return b;
    for (int j = static_cast<int>(up.size()) - 1; j >= 0; --j) {
        NodeId next = up[j][a];
        if (!is_ancestor(next, b)) a = next;
    }
    return up[0][a];
}

NodeId Tree::lca_leaves(const std::vector<Leaf>& leaves) const {
    if (leaves.empty()) return -1;
    NodeId cur = leaf_node[leaves[0]];
    for (size_t i = 1; i < leaves.size(); ++i) {
        cur = lca(cur, leaf_node[leaves[i]]);
    }
    return cur;
}
