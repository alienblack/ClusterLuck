#include "lowerbound/lower_bound.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

namespace {

std::vector<std::pair<Leaf, Leaf>> enumerate_leaf_cherries(const Tree& tree) {
    std::vector<std::pair<Leaf, Leaf>> cherries;
    for (NodeId v = 0; v < tree.node_count(); ++v) {
        NodeId a = tree.left[v];
        NodeId b = tree.right[v];
        if (a == -1 || b == -1) continue;
        Leaf la = tree.label[a];
        Leaf lb = tree.label[b];
        if (la < 1 || lb < 1) continue;
        if (la > lb) std::swap(la, lb);
        cherries.push_back({la, lb});
    }
    return cherries;
}

uint64_t cherry_key(Leaf a, Leaf b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<unsigned>(a)) << 32) |
           static_cast<unsigned>(b);
}

int subtree_node_count(const Tree& tree, NodeId v) {
    return tree.tout[v] - tree.tin[v];
}

int pack_conflicting_cherry_subtrees(const Tree& source, const Tree& other) {
    std::vector<std::pair<Leaf, Leaf>> source_cherries =
        enumerate_leaf_cherries(source);
    std::vector<std::pair<Leaf, Leaf>> other_cherries =
        enumerate_leaf_cherries(other);

    std::unordered_set<uint64_t> other_set;
    other_set.reserve(other_cherries.size() * 2 + 1);
    for (const auto& p : other_cherries) {
        other_set.insert(cherry_key(p.first, p.second));
    }

    std::vector<NodeId> roots;
    std::unordered_set<NodeId> seen_roots;
    seen_roots.reserve(source_cherries.size() * 2 + 1);
    for (const auto& p : source_cherries) {
        if (other_set.count(cherry_key(p.first, p.second)) != 0) continue;
        NodeId a = other.leaf_node[p.first];
        NodeId b = other.leaf_node[p.second];
        NodeId root = other.lca(a, b);
        if (seen_roots.insert(root).second) roots.push_back(root);
    }

    std::sort(roots.begin(), roots.end(), [&](NodeId a, NodeId b) {
        int size_a = subtree_node_count(other, a);
        int size_b = subtree_node_count(other, b);
        if (size_a != size_b) return size_a < size_b;
        if (other.depth[a] != other.depth[b]) {
            return other.depth[a] > other.depth[b];
        }
        return a < b;
    });

    std::vector<NodeId> chosen;
    int packed = 0;
    for (NodeId root : roots) {
        bool overlaps = false;
        for (NodeId prev : chosen) {
            if (other.is_ancestor(root, prev) || other.is_ancestor(prev, root)) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;
        chosen.push_back(root);
        ++packed;
    }
    return packed;
}

int compute_cherry_packing_lower_bound(const Tree& t1, const Tree& t2) {
    return 1 + std::max(pack_conflicting_cherry_subtrees(t1, t2),
                        pack_conflicting_cherry_subtrees(t2, t1));
}

}  // namespace

LowerBoundCert compute_lower_bound_certificate(const Instance& instance) {
    LowerBoundCert cert;
    if (!instance.has_lower_bound_rule) return cert;
    cert.enabled = true;
    cert.lower_bound =
        std::max(1, compute_cherry_packing_lower_bound(instance.t1, instance.t2));
    cert.target_components =
        static_cast<int>(std::floor(instance.lower_bound_factor *
                                    static_cast<double>(cert.lower_bound) +
                                    1e-12)) +
        instance.lower_bound_slack;
    cert.target_components = std::max(1, cert.target_components);
    return cert;
}

bool lower_bound_target_met(int components, const LowerBoundCert& cert) {
    return cert.enabled && components <= cert.target_components;
}
