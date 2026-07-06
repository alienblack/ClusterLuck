#include "forest/forest.h"

#include <algorithm>
#include <stdexcept>

#include "tree/induced.h"

Forest make_singleton_forest(int n) {
    Forest forest;
    forest.comps.reserve(n);
    forest.leaf_comp.assign(n + 1, -1);
    for (Leaf x = 1; x <= n; ++x) {
        Component c;
        c.leaves.push_back(x);
        c.repr = x;
        forest.leaf_comp[x] = static_cast<CompId>(forest.comps.size());
        forest.comps.push_back(std::move(c));
    }
    forest.active_count = n;
    return forest;
}

bool rebuild_forest_data(Forest& forest, const Tree& t1, const Tree& t2) {
    forest.owner1.assign(t1.node_count(), -1);
    forest.owner2.assign(t2.node_count(), -1);

    int active = 0;
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        Component& c = forest.comps[id];
        if (!c.active) continue;
        if (c.repr < 1 && !c.leaves.empty()) c.repr = c.leaves.front();
        ++active;

        InducedInfo i1 = compute_induced(t1, c.leaves);
        InducedInfo i2 = compute_induced(t2, c.leaves);
        c.root1 = i1.root;
        c.root2 = i2.root;
        c.signature1 = std::move(i1.signature);
        c.signature2 = std::move(i2.signature);
        c.edges1 = std::move(i1.used_edges);
        c.edges2 = std::move(i2.used_edges);

        if (c.signature1 != c.signature2) return false;
        for (NodeId e : c.edges1) {
            if (forest.owner1[e] != -1) return false;
            forest.owner1[e] = id;
        }
        for (NodeId e : c.edges2) {
            if (forest.owner2[e] != -1) return false;
            forest.owner2[e] = id;
        }
    }
    return active == forest.active_count;
}

std::vector<CompId> active_components(const Forest& forest) {
    std::vector<CompId> ids;
    ids.reserve(forest.active_count);
    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        if (forest.comps[id].active) ids.push_back(id);
    }
    return ids;
}

std::string forest_to_output(const Forest& forest, const Tree& t1) {
    std::string out;
    for (const Component& c : forest.comps) {
        if (!c.active) continue;
        InducedInfo induced = compute_induced(t1, c.leaves);
        out += induced.newick;
        out += ";\n";
    }
    return out;
}

void merge_components(Forest& forest, CompId a, CompId b, std::vector<Leaf> merged_leaves) {
    if (a == b || !forest.comps[a].active || !forest.comps[b].active) {
        throw std::runtime_error("invalid merge ids");
    }
    std::sort(merged_leaves.begin(), merged_leaves.end());
    merged_leaves.erase(std::unique(merged_leaves.begin(), merged_leaves.end()),
                        merged_leaves.end());

    forest.comps[a].leaves = std::move(merged_leaves);
    forest.comps[b].active = false;
    for (Leaf x : forest.comps[a].leaves) forest.leaf_comp[x] = a;
    --forest.active_count;
}
