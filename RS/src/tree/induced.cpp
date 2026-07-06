#include "tree/induced.h"

#include <algorithm>
#include <stdexcept>

namespace {

struct Builder {
    const Tree& t;
    std::vector<int> selected_count;

    explicit Builder(const Tree& tree) : t(tree), selected_count(tree.node_count(), 0) {}

    void mark_path(NodeId leaf, NodeId root) {
        NodeId v = leaf;
        while (true) {
            ++selected_count[v];
            if (v == root) break;
            v = t.parent[v];
            if (v == -1) throw std::runtime_error("bad induced root");
        }
    }

    std::pair<std::string, std::string> build(NodeId v) {
        if (selected_count[v] == 0) return {"", ""};
        if (t.is_leaf(v)) {
            std::string s = std::to_string(t.label[v]);
            return {s, s};
        }

        auto left_part = build(t.left[v]);
        auto right_part = build(t.right[v]);
        bool has_left = !left_part.first.empty();
        bool has_right = !right_part.first.empty();

        if (has_left && has_right) {
            std::string a = left_part.first;
            std::string b = right_part.first;
            if (b < a) std::swap(a, b);
            std::string sig = "(" + a + "," + b + ")";
            std::string nw = "(" + left_part.second + "," + right_part.second + ")";
            return {sig, nw};
        }
        if (has_left) return left_part;
        if (has_right) return right_part;
        return {"", ""};
    }
};

}  // namespace

InducedInfo compute_induced(const Tree& tree, const std::vector<Leaf>& leaves) {
    if (leaves.empty()) throw std::runtime_error("empty component");

    InducedInfo info;
    info.root = tree.lca_leaves(leaves);

    Builder builder(tree);
    for (Leaf x : leaves) {
        NodeId leaf = tree.leaf_node[x];
        builder.mark_path(leaf, info.root);
        NodeId v = leaf;
        while (v != info.root) {
            info.used_edges.push_back(v);
            v = tree.parent[v];
        }
    }

    std::sort(info.used_edges.begin(), info.used_edges.end());
    info.used_edges.erase(std::unique(info.used_edges.begin(), info.used_edges.end()),
                          info.used_edges.end());

    auto built = builder.build(info.root);
    info.signature = built.first;
    info.newick = built.second;
    if (info.signature.empty()) throw std::runtime_error("failed induced tree");
    return info;
}

