#include "kernel/kernel_reduce.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "forest/validator.h"
#include "io/newick.h"

namespace {

struct MutableTree {
    struct Node {
        int parent = -1;
        int left = -1;
        int right = -1;
        int label = 0;
        bool alive = true;
    };

    int n = 0;
    int root = -1;
    std::vector<Node> nodes;
    std::vector<int> leaf_node;

    static MutableTree from_tree(const Tree& t) {
        MutableTree mt;
        mt.n = t.n;
        mt.root = t.root;
        mt.nodes.resize(t.node_count());
        mt.leaf_node.assign(t.n + 1, -1);
        for (int v = 0; v < t.node_count(); ++v) {
            mt.nodes[v].parent = t.parent[v];
            mt.nodes[v].left = t.left[v];
            mt.nodes[v].right = t.right[v];
            mt.nodes[v].label = t.label[v];
            if (t.label[v] > 0) mt.leaf_node[t.label[v]] = v;
        }
        return mt;
    }

    bool is_leaf(int v) const {
        return v >= 0 && nodes[v].alive && nodes[v].left == -1 &&
               nodes[v].right == -1;
    }

    int other_child(int p, int child) const {
        if (p < 0) return -1;
        if (nodes[p].left == child) return nodes[p].right;
        if (nodes[p].right == child) return nodes[p].left;
        return -1;
    }

    int leaf_node_of(int taxon) const {
        if (taxon <= 0 || taxon >= static_cast<int>(leaf_node.size())) return -1;
        int v = leaf_node[taxon];
        if (v < 0 || !nodes[v].alive) return -1;
        return v;
    }

    int part_of_cherry(int taxon) const {
        int v = leaf_node_of(taxon);
        if (v < 0) return -1;
        int p = nodes[v].parent;
        if (p < 0 || !nodes[p].alive) return -1;
        int sibling = other_child(p, v);
        if (!is_leaf(sibling)) return -1;
        return nodes[sibling].label;
    }

    int number_leaf_children(int v) const {
        if (v < 0 || !nodes[v].alive) return 0;
        return (is_leaf(nodes[v].left) ? 1 : 0) +
               (is_leaf(nodes[v].right) ? 1 : 0);
    }

    int get_non_leaf_child(int v) const {
        if (v < 0 || !nodes[v].alive) return -1;
        int l = nodes[v].left;
        if (l >= 0 && nodes[l].alive && !is_leaf(l)) return l;
        int r = nodes[v].right;
        if (r >= 0 && nodes[r].alive && !is_leaf(r)) return r;
        return -1;
    }

    int get_leaf_child_label(int v) const {
        if (v < 0 || !nodes[v].alive) return -1;
        if (is_leaf(nodes[v].left)) return nodes[nodes[v].left].label;
        if (is_leaf(nodes[v].right)) return nodes[nodes[v].right].label;
        return -1;
    }

    bool has_leaf_child(int v, int taxon) const {
        if (v < 0 || !nodes[v].alive) return false;
        int l = nodes[v].left;
        int r = nodes[v].right;
        return (is_leaf(l) && nodes[l].label == taxon) ||
               (is_leaf(r) && nodes[r].label == taxon);
    }

    void delete_taxon(int taxon) {
        int v = leaf_node_of(taxon);
        if (v < 0) return;
        int p = nodes[v].parent;
        if (p < 0) {
            nodes[v].alive = false;
            leaf_node[taxon] = -1;
            root = -1;
            return;
        }
        int sibling = other_child(p, v);
        int gp = nodes[p].parent;
        if (gp < 0) {
            root = sibling;
            nodes[sibling].parent = -1;
        } else {
            if (nodes[gp].left == p) {
                nodes[gp].left = sibling;
            } else if (nodes[gp].right == p) {
                nodes[gp].right = sibling;
            } else {
                throw std::runtime_error("broken parent pointer in kernel delete");
            }
            nodes[sibling].parent = gp;
        }
        nodes[v].alive = false;
        nodes[p].alive = false;
        leaf_node[taxon] = -1;
    }

    std::string newick_rec(int v) const {
        if (v < 0 || !nodes[v].alive) return "";
        if (is_leaf(v)) return std::to_string(nodes[v].label);
        return "(" + newick_rec(nodes[v].left) + "," +
               newick_rec(nodes[v].right) + ")";
    }

    std::string newick() const {
        return newick_rec(root) + ";";
    }
};

struct KernelStats {
    int leaves_left = 0;
    int subtree_deleted = 0;
    int chain_deleted = 0;
    int param_reduction = 0;
    std::vector<char> active;
    std::vector<std::vector<int>> forced_components;
    std::vector<KernelMap::ChainRecord> chain_records;
};

std::string join_payload(const std::vector<int>& payload) {
    std::ostringstream out;
    for (size_t i = 0; i < payload.size(); ++i) {
        if (i) out << ",";
        out << payload[i];
    }
    return out.str();
}

KernelStats reduce_kernel(MutableTree& t1, MutableTree& t2,
                          std::vector<std::vector<int>>& payloads,
                          bool verbose, int protected_taxon) {
    const int n = static_cast<int>(payloads.size()) - 1;
    std::vector<char> active(n + 1, true);
    int leaves_left = n;
    KernelStats stats;

    bool finished = false;
    while (!finished) {
        finished = true;
        if (leaves_left <= 2) break;

        for (int x = 1; x <= n; ++x) {
            if (!active[x]) continue;
            if (x == protected_taxon) continue;
            int sibling_t1 = t1.part_of_cherry(x);
            int sibling_t2 = t2.part_of_cherry(x);
            if (sibling_t1 < 0 || sibling_t2 < 0) continue;
            if (sibling_t1 != sibling_t2) continue;

            int y = sibling_t1;
            if (y == protected_taxon) continue;
            active[y] = false;
            --leaves_left;
            ++stats.subtree_deleted;
            payloads[x].insert(payloads[x].end(), payloads[y].begin(),
                               payloads[y].end());
            std::sort(payloads[x].begin(), payloads[x].end());
            t1.delete_taxon(y);
            t2.delete_taxon(y);
            if (verbose) {
                std::cerr << "// Reduced common cherry (" << x << "," << y
                          << ") into single taxon " << x << "\n";
            }
            finished = false;
            goto next_reduce_pass;
        }

        for (int x = 1; x <= n; ++x) {
            if (!active[x]) continue;
            if (x == protected_taxon) continue;
            int chain_base[4] = {0, x, 0, 0};
            int at_t1 = t1.nodes[t1.leaf_node_of(x)].parent;
            int at_t2 = t2.nodes[t2.leaf_node_of(x)].parent;

            if (t1.number_leaf_children(at_t1) != 1) continue;
            if (t2.number_leaf_children(at_t2) != 1) continue;
            at_t1 = t1.get_non_leaf_child(at_t1);
            at_t2 = t2.get_non_leaf_child(at_t2);
            if (t1.number_leaf_children(at_t1) != 1) continue;
            if (t2.number_leaf_children(at_t2) != 1) continue;
            int c1 = t1.get_leaf_child_label(at_t1);
            int c2 = t2.get_leaf_child_label(at_t2);
            if (c1 != c2) continue;
            if (c1 == protected_taxon) continue;
            chain_base[2] = c1;
            at_t1 = t1.get_non_leaf_child(at_t1);
            at_t2 = t2.get_non_leaf_child(at_t2);

            int leaf_count1 = t1.number_leaf_children(at_t1);
            int leaf_count2 = t2.number_leaf_children(at_t2);
            if (leaf_count1 != 1 && leaf_count2 != 1) continue;

            if (leaf_count1 == 1 && leaf_count2 == 1) {
                c1 = t1.get_leaf_child_label(at_t1);
                c2 = t2.get_leaf_child_label(at_t2);
                if (c1 != c2) continue;
                if (c1 == protected_taxon) continue;
                chain_base[3] = c1;
                at_t1 = t1.get_non_leaf_child(at_t1);
                at_t2 = t2.get_non_leaf_child(at_t2);
            } else if (leaf_count1 == 1 && leaf_count2 == 2) {
                c1 = t1.get_leaf_child_label(at_t1);
                if (c1 == protected_taxon) continue;
                if (!t2.has_leaf_child(at_t2, c1)) continue;
                chain_base[3] = c1;
                at_t1 = t1.get_non_leaf_child(at_t1);
            } else if (leaf_count1 == 2 && leaf_count2 == 1) {
                c2 = t2.get_leaf_child_label(at_t2);
                if (c2 == protected_taxon) continue;
                if (!t1.has_leaf_child(at_t1, c2)) continue;
                chain_base[3] = c2;
                at_t2 = t2.get_non_leaf_child(at_t2);
            } else {
                continue;
            }

            if (at_t1 < 0 || at_t2 < 0) continue;
            for (int y = 1; y <= n; ++y) {
                if (!active[y]) continue;
                if (y == protected_taxon) continue;
                if (!t1.has_leaf_child(at_t1, y) ||
                    !t2.has_leaf_child(at_t2, y)) {
                    continue;
                }
                active[y] = false;
                --leaves_left;
                ++stats.chain_deleted;
                stats.chain_records.push_back({payloads[chain_base[1]],
                                               payloads[chain_base[2]],
                                               payloads[chain_base[3]],
                                               payloads[y]});
                if (verbose) {
                    std::cerr << "// Chain reduction deleted taxon " << y
                              << " representing: " << join_payload(payloads[y])
                              << "\n";
                }
                t1.delete_taxon(y);
                t2.delete_taxon(y);
                finished = false;
                goto next_reduce_pass;
            }
        }

        for (int p = 1; p <= n; ++p) {
            if (!active[p]) continue;
            if (p == protected_taxon) continue;
            int sibling_t1 = t1.part_of_cherry(p);
            int sibling_t2 = t2.part_of_cherry(p);
            if (sibling_t1 < 0 || sibling_t2 < 0) continue;
            if (sibling_t1 == protected_taxon ||
                sibling_t2 == protected_taxon) {
                continue;
            }

            int q_in_t1 = sibling_t1;
            int q_node_t2 = t2.leaf_node_of(q_in_t1);
            int r_in_t2 = sibling_t2;
            int r_node_t2 = t2.leaf_node_of(r_in_t2);
            if (q_node_t2 >= 0 && r_node_t2 >= 0) {
                int r_parent = t2.nodes[r_node_t2].parent;
                int r_grandparent =
                    r_parent >= 0 ? t2.nodes[r_parent].parent : -1;
                int q_parent = t2.nodes[q_node_t2].parent;
                if (r_grandparent == q_parent) {
                    int taxon = r_in_t2;
                    if (taxon == protected_taxon) continue;
                    active[taxon] = false;
                    --leaves_left;
                    ++stats.param_reduction;
                    stats.forced_components.push_back(payloads[taxon]);
                    if (verbose) {
                        std::cerr << "// 3-2 chain reduction deleted taxon "
                                  << taxon << " representing: "
                                  << join_payload(payloads[taxon]) << "\n";
                    }
                    t1.delete_taxon(taxon);
                    t2.delete_taxon(taxon);
                    finished = false;
                    goto next_reduce_pass;
                }
            }

            int q_in_t2 = sibling_t2;
            int q_node_t1 = t1.leaf_node_of(q_in_t2);
            int r_in_t1 = sibling_t1;
            int r_node_t1 = t1.leaf_node_of(r_in_t1);
            if (q_node_t1 >= 0 && r_node_t1 >= 0) {
                int r_parent = t1.nodes[r_node_t1].parent;
                int r_grandparent =
                    r_parent >= 0 ? t1.nodes[r_parent].parent : -1;
                int q_parent = t1.nodes[q_node_t1].parent;
                if (r_grandparent == q_parent) {
                    int taxon = r_in_t1;
                    if (taxon == protected_taxon) continue;
                    active[taxon] = false;
                    --leaves_left;
                    ++stats.param_reduction;
                    stats.forced_components.push_back(payloads[taxon]);
                    if (verbose) {
                        std::cerr << "// 3-2 chain reduction deleted taxon "
                                  << taxon << " representing: "
                                  << join_payload(payloads[taxon]) << "\n";
                    }
                    t1.delete_taxon(taxon);
                    t2.delete_taxon(taxon);
                    finished = false;
                    goto next_reduce_pass;
                }
            }
        }

    next_reduce_pass:
        continue;
    }

    stats.leaves_left = leaves_left;
    stats.active = std::move(active);
    return stats;
}

Forest forest_from_components(std::vector<std::vector<int>> components,
                              const Instance& instance) {
    Forest forest;
    forest.leaf_comp.assign(instance.n + 1, -1);
    std::vector<unsigned char> seen(instance.n + 1, 0);
    for (auto& leaves : components) {
        std::sort(leaves.begin(), leaves.end());
        leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
        if (leaves.empty()) continue;
        Component component;
        component.leaves = std::move(leaves);
        component.repr = component.leaves.front();
        for (int leaf : component.leaves) {
            if (leaf < 1 || leaf > instance.n || seen[leaf]) {
                throw std::runtime_error("expanded forest does not cover once");
            }
            seen[leaf] = 1;
        }
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    for (int leaf = 1; leaf <= instance.n; ++leaf) {
        if (seen[leaf]) continue;
        Component component;
        component.leaves = {leaf};
        component.repr = leaf;
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    if (!rebuild_forest_data(forest, instance.t1, instance.t2)) {
        throw std::runtime_error("expanded forest rebuild failed");
    }
    return forest;
}

int component_containing_any(const std::vector<std::vector<int>>& components,
                             const std::vector<int>& payload) {
    for (size_t id = 0; id < components.size(); ++id) {
        for (int leaf : payload) {
            if (std::find(components[id].begin(), components[id].end(), leaf) !=
                components[id].end()) {
                return static_cast<int>(id);
            }
        }
    }
    return -1;
}

}  // namespace

std::optional<KernelReducedInstance> kernel_reduce_instance(
    const Instance& instance, bool verbose, int protected_taxon) {
    MutableTree t1 = MutableTree::from_tree(instance.t1);
    MutableTree t2 = MutableTree::from_tree(instance.t2);
    std::vector<std::vector<int>> payloads(instance.n + 1);
    for (int x = 1; x <= instance.n; ++x) payloads[x] = {x};

    KernelStats stats = reduce_kernel(t1, t2, payloads, verbose,
                                      protected_taxon);
    if (stats.leaves_left <= 0 || stats.leaves_left >= instance.n) {
        return std::nullopt;
    }

    std::vector<int> active;
    active.reserve(stats.leaves_left);
    for (int x = 1; x <= instance.n; ++x) {
        if (x < static_cast<int>(stats.active.size()) && stats.active[x]) {
            active.push_back(x);
        }
    }
    std::unordered_map<int, int> compact;
    for (int old : active) {
        compact[old] = static_cast<int>(compact.size()) + 1;
    }

    auto compact_newick = [&compact](const std::string& text) {
        std::string out;
        for (size_t i = 0; i < text.size();) {
            if (!std::isdigit(static_cast<unsigned char>(text[i]))) {
                out.push_back(text[i++]);
                continue;
            }
            int value = 0;
            while (i < text.size() &&
                   std::isdigit(static_cast<unsigned char>(text[i]))) {
                value = value * 10 + (text[i] - '0');
                ++i;
            }
            auto found = compact.find(value);
            if (found == compact.end()) {
                throw std::runtime_error("active taxon missing from compact map");
            }
            out += std::to_string(found->second);
        }
        return out;
    };

    KernelReducedInstance result;
    result.original_n = instance.n;
    result.subtree_deleted = stats.subtree_deleted;
    result.chain_deleted = stats.chain_deleted;
    result.param_reduction = stats.param_reduction;
    result.reduced.n = static_cast<int>(active.size());
    result.reduced.t1 =
        parse_newick_tree_h45_order(compact_newick(t1.newick()), result.reduced.n);
    result.reduced.t2 =
        parse_newick_tree_h45_order(compact_newick(t2.newick()), result.reduced.n);

    for (int old : active) {
        int reduced_label = compact[old];
        result.map.active_payload[reduced_label] = payloads[old];
    }
    result.map.forced = std::move(stats.forced_components);
    result.map.chains = std::move(stats.chain_records);

    return result;
}

Forest expand_kernel_reduced_forest(const Forest& reduced_forest,
                                    const KernelMap& map,
                                    const Instance& original) {
    std::vector<std::vector<int>> components = map.forced;
    for (const Component& component : reduced_forest.comps) {
        if (!component.active) continue;
        std::vector<int> expanded;
        for (int label : component.leaves) {
            auto found = map.active_payload.find(label);
            if (found == map.active_payload.end()) {
                throw std::runtime_error("missing reduced label in kernel map");
            }
            expanded.insert(expanded.end(), found->second.begin(),
                            found->second.end());
        }
        if (!expanded.empty()) components.push_back(std::move(expanded));
    }
    for (const auto& chain : map.chains) {
        int c1 = component_containing_any(components, chain.anchor1);
        int c2 = component_containing_any(components, chain.anchor2);
        int c3 = component_containing_any(components, chain.anchor3);
        if (c1 >= 0 && c1 == c2 && c1 == c3) {
            components[c1].insert(components[c1].end(), chain.deleted.begin(),
                                  chain.deleted.end());
        } else {
            components.push_back(chain.deleted);
        }
    }

    Forest forest = forest_from_components(std::move(components), original);
    ValidationResult check = validate_forest(forest, original.t1, original.t2);
    if (!check.ok) {
        throw std::runtime_error("expanded kernel forest invalid: " + check.message);
    }
    return forest;
}
