#include "forest/validator.h"

#include <algorithm>
#include <vector>

#include "tree/induced.h"

namespace {

bool claim_edges(const std::vector<NodeId>& edges, CompId id, std::vector<CompId>& owner) {
    for (NodeId edge : edges) {
        if (owner[edge] != -1 && owner[edge] != id) return false;
        owner[edge] = id;
    }
    return true;
}

}  // namespace

ValidationResult validate_forest(const Forest& forest, const Tree& t1, const Tree& t2) {
    ValidationResult result;
    int n = t1.n;
    if (t2.n != n) return {false, "tree leaf counts differ"};

    std::vector<int> seen(n + 1, 0);
    int active = 0;
    std::vector<CompId> owner1(t1.node_count(), -1);
    std::vector<CompId> owner2(t2.node_count(), -1);

    for (CompId id = 0; id < static_cast<CompId>(forest.comps.size()); ++id) {
        const Component& c = forest.comps[id];
        if (!c.active) continue;
        ++active;
        if (c.leaves.empty()) return {false, "empty active component"};
        for (Leaf x : c.leaves) {
            if (x < 1 || x > n) return {false, "leaf out of range"};
            if (++seen[x] != 1) return {false, "duplicate leaf"};
        }

        InducedInfo i1 = compute_induced(t1, c.leaves);
        InducedInfo i2 = compute_induced(t2, c.leaves);
        if (i1.signature != i2.signature) return {false, "component topology mismatch"};
        if (!claim_edges(i1.used_edges, id, owner1)) return {false, "T1 edge overlap"};
        if (!claim_edges(i2.used_edges, id, owner2)) return {false, "T2 edge overlap"};
    }

    if (active != forest.active_count) return {false, "active count mismatch"};
    for (Leaf x = 1; x <= n; ++x) {
        if (seen[x] != 1) return {false, "missing leaf"};
    }
    result.ok = true;
    result.message = "ok";
    return result;
}

