#include "cluster/cluster_points.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct ClusterContext {
    ClusterContext(const Tree& first, const Tree& second) : t1(first), t2(second) {}

    const Tree& t1;
    const Tree& t2;

    std::vector<int> count1;
    std::vector<int> count2;
    std::vector<NodeId> twin1;
    std::vector<NodeId> twin2;
};

struct ClusterToken {
    NodeId root = -1;
    int atomic_count = 1;
};

NodeId sibling_of(const Tree& tree, NodeId v) {
    NodeId p = tree.parent[v];
    if (p < 0) return -1;
    return tree.left[p] == v ? tree.right[p] : tree.left[p];
}

int subtree_count(const std::vector<int>& counts, NodeId v) {
    return v < 0 ? 0 : counts[v];
}

int atomic_count(const Tree& tree, const std::vector<int>& counts, NodeId v,
                 const std::vector<ClusterToken>& tokens) {
    if (v < 0) return 0;
    int total = subtree_count(counts, v);
    for (size_t i = 0; i < tokens.size(); ++i) {
        const ClusterToken& token = tokens[i];
        if (token.root < 0 || !tree.is_ancestor(v, token.root)) continue;

        bool covered_by_larger_token = false;
        for (size_t j = 0; j < tokens.size(); ++j) {
            if (i == j) continue;
            const ClusterToken& other = tokens[j];
            if (other.root >= 0 &&
                other.root != token.root &&
                tree.is_ancestor(v, other.root) &&
                tree.is_ancestor(other.root, token.root)) {
                covered_by_larger_token = true;
                break;
            }
        }
        if (!covered_by_larger_token) {
            total -= subtree_count(counts, token.root) - token.atomic_count;
        }
    }
    return total;
}

int fill_counts(const Tree& tree, NodeId v, std::vector<int>& counts) {
    if (tree.is_leaf(v)) {
        counts[v] = 1;
    } else {
        counts[v] = fill_counts(tree, tree.left[v], counts) +
                    fill_counts(tree, tree.right[v], counts);
    }
    return counts[v];
}

NodeId fill_twin_lcas(const Tree& source, const Tree& other, NodeId v,
                      std::vector<NodeId>& twins) {
    if (source.is_leaf(v)) {
        twins[v] = other.leaf_node[source.label[v]];
    } else {
        NodeId a = fill_twin_lcas(source, other, source.left[v], twins);
        NodeId b = fill_twin_lcas(source, other, source.right[v], twins);
        twins[v] = other.lca(a, b);
    }
    return twins[v];
}

bool is_cluster_point_raw(const ClusterContext& ctx, NodeId v) {
    if (v < 0 || ctx.t1.is_leaf(v) || ctx.t1.parent[v] < 0) return false;
    NodeId twin = ctx.twin1[v];
    if (twin < 0) return false;
    NodeId reverse = ctx.twin2[twin];
    if (reverse < 0) return false;
    if (ctx.t1.depth[v] > ctx.t1.depth[reverse]) return false;

    int clustered_children = 0;
    for (NodeId child : {ctx.t1.left[v], ctx.t1.right[v]}) {
        NodeId child_twin = ctx.twin1[child];
        if (child_twin < 0) continue;
        NodeId child_reverse = ctx.twin2[child_twin];
        if (child_reverse >= 0 &&
            ctx.t1.depth[child] <= ctx.t1.depth[child_reverse]) {
            ++clustered_children;
        }
    }
    return clustered_children < 2;
}

bool auto_skip_oversplit(ClusterCandidate& candidate) {
    const char* disable = std::getenv("RS_DISABLE_CLUSTER_AUTOSKIP");
    if (disable != nullptr && *disable != '\0' && std::string(disable) != "0") {
        return false;
    }
    if (candidate.child_first <= 0 || candidate.child_second <= 0 ||
        candidate.parent_first <= 0 || candidate.parent_second <= 0) {
        return false;
    }

    int child_min = std::min(candidate.child_first, candidate.child_second);
    int child_max = std::max(candidate.child_first, candidate.child_second);
    int small_parent = std::min(candidate.parent_first, candidate.parent_second);
    int large_parent = std::max(candidate.parent_first, candidate.parent_second);

    if (candidate.child_first != candidate.child_second ||
        candidate.parent_first != candidate.parent_second) {
        int pendant_spill = small_parent - child_min;
        if (child_min >= 10 && child_max <= 40 &&
            child_max - child_min <= 2 &&
            pendant_spill >= 1 && pendant_spill <= 2 &&
            large_parent >= std::max(200, 8 * child_max) &&
            large_parent >= 4 * small_parent) {
            candidate.skipped = true;
            candidate.skip_rule = "asymmetric_pendant";
            return true;
        }
        return false;
    }

    if (candidate.parent_first <= candidate.child_first) return false;
    int spill = candidate.parent_first - candidate.child_first;
    if (candidate.child_first >= 60 &&
        candidate.child_first * 100 >= 65 * candidate.parent_first &&
        spill >= 10 &&
        spill * 4 <= candidate.parent_first &&
        spill <= 25) {
        candidate.skipped = true;
        candidate.skip_rule = "over_split_parent";
        return true;
    }
    return false;
}

void collect_candidates_postorder(const ClusterContext& ctx, NodeId v,
                                  std::vector<ClusterCandidate>& out,
                                  int& original_id,
                                  std::vector<ClusterToken>& tokens1,
                                  std::vector<ClusterToken>& tokens2) {
    if (ctx.t1.is_leaf(v)) return;
    collect_candidates_postorder(ctx, ctx.t1.left[v], out, original_id,
                                 tokens1, tokens2);
    collect_candidates_postorder(ctx, ctx.t1.right[v], out, original_id,
                                 tokens1, tokens2);

    if (!is_cluster_point_raw(ctx, v)) return;
    NodeId twin = ctx.twin1[v];
    if (twin < 0 || ctx.t2.parent[twin] < 0) return;

    ClusterCandidate candidate;
    candidate.original_id = ++original_id;
    candidate.node1 = v;
    candidate.twin2 = twin;
    candidate.child_first = atomic_count(ctx.t1, ctx.count1, v, tokens1);
    candidate.original_child_first = subtree_count(ctx.count1, v);
    candidate.parent_first =
        atomic_count(ctx.t1, ctx.count1, ctx.t1.parent[v], tokens1);
    candidate.sibling_first =
        atomic_count(ctx.t1, ctx.count1, sibling_of(ctx.t1, v), tokens1);
    candidate.child_second = atomic_count(ctx.t2, ctx.count2, twin, tokens2);
    candidate.parent_second =
        atomic_count(ctx.t2, ctx.count2, ctx.t2.parent[twin], tokens2);
    candidate.sibling_second =
        atomic_count(ctx.t2, ctx.count2, sibling_of(ctx.t2, twin), tokens2);
    auto_skip_oversplit(candidate);
    if (!candidate.skipped) {
        tokens1.push_back({v, 1});
        tokens2.push_back({twin, 1});
    }
    out.push_back(candidate);
}

}  // namespace

std::vector<ClusterCandidate> find_cluster_candidates(const Tree& t1,
                                                      const Tree& t2) {
    ClusterContext ctx(t1, t2);
    ctx.count1.assign(t1.node_count(), 0);
    ctx.count2.assign(t2.node_count(), 0);
    ctx.twin1.assign(t1.node_count(), -1);
    ctx.twin2.assign(t2.node_count(), -1);
    fill_counts(t1, t1.root, ctx.count1);
    fill_counts(t2, t2.root, ctx.count2);
    fill_twin_lcas(t1, t2, t1.root, ctx.twin1);
    fill_twin_lcas(t2, t1, t2.root, ctx.twin2);

    std::vector<ClusterCandidate> candidates;
    int original_id = 0;
    std::vector<ClusterToken> tokens1;
    std::vector<ClusterToken> tokens2;
    collect_candidates_postorder(ctx, t1.root, candidates, original_id,
                                 tokens1, tokens2);
    return candidates;
}

void dump_cluster_candidates(const std::vector<ClusterCandidate>& candidates,
                             std::ostream& out) {
    const bool debug_nodes = std::getenv("RS_DUMP_CLUSTER_NODE_IDS") != nullptr;
    int retained = 0;
    std::vector<int> retained_sizes;
    for (const ClusterCandidate& candidate : candidates) {
        out << "cluster_candidate original_id=" << candidate.original_id
            << " child_first=" << candidate.child_first
            << " parent_first=" << candidate.parent_first
            << " sibling_first=" << candidate.sibling_first
            << " child_second=" << candidate.child_second
            << " parent_second=" << candidate.parent_second
            << " sibling_second=" << candidate.sibling_second << '\n';
        if (debug_nodes) {
            out << "rs_cluster_node original_id=" << candidate.original_id
                << " node1=" << candidate.node1
                << " twin2=" << candidate.twin2 << '\n';
        }
        if (candidate.skipped) {
            out << "cluster_candidate_auto_skip original_id="
                << candidate.original_id
                << " child_leaves=" << candidate.child_first;
            if (std::string(candidate.skip_rule) == "asymmetric_pendant") {
                out << " parent_first=" << candidate.parent_first
                    << " parent_second=" << candidate.parent_second;
            } else {
                out << " parent_leaves=" << candidate.parent_first
                    << " spill=" << candidate.parent_first - candidate.child_first;
            }
            out << " rule=" << candidate.skip_rule << '\n';
        } else {
            ++retained;
            retained_sizes.push_back(candidate.child_first);
        }
    }
    std::sort(retained_sizes.begin(), retained_sizes.end(), std::greater<int>());
    out << "rs_cluster_summary candidates=" << candidates.size()
        << " retained=" << retained
        << " retained_sizes=";
    for (size_t i = 0; i < retained_sizes.size(); ++i) {
        if (i) out << ',';
        out << retained_sizes[i];
    }
    out << '\n';
}

ClusterPlan build_cluster_plan(const Tree& t1,
                               const std::vector<ClusterCandidate>& candidates) {
    ClusterPlan plan;
    for (const ClusterCandidate& candidate : candidates) {
        if (candidate.skipped) continue;
        PlannedCluster item;
        item.id = static_cast<int>(plan.clusters.size()) + 1;
        item.original_id = candidate.original_id;
        item.node1 = candidate.node1;
        item.node2 = candidate.twin2;
        item.local_leaves = candidate.child_first;
        item.original_leaves = candidate.original_child_first;
        plan.clusters.push_back(item);
    }

    for (PlannedCluster& child : plan.clusters) {
        int best_parent = 0;
        int best_original_leaves = t1.n + 1;
        for (const PlannedCluster& parent : plan.clusters) {
            if (parent.id == child.id) continue;
            if (parent.original_leaves <= child.original_leaves) continue;
            if (!t1.is_ancestor(parent.node1, child.node1)) continue;
            if (parent.original_leaves < best_original_leaves) {
                best_parent = parent.id;
                best_original_leaves = parent.original_leaves;
            }
        }
        child.parent = best_parent;
    }

    plan.residual_id = static_cast<int>(plan.clusters.size()) + 1;
    plan.residual_leaves = t1.n;
    for (const PlannedCluster& cluster : plan.clusters) {
        if (cluster.parent != 0) continue;
        plan.residual_leaves -= cluster.original_leaves - 1;
    }

    for (const PlannedCluster& cluster : plan.clusters) {
        plan.solve_order.push_back(cluster.id);
    }
    std::sort(plan.solve_order.begin(), plan.solve_order.end(), [&](int a, int b) {
        const PlannedCluster& ca = plan.clusters[a - 1];
        const PlannedCluster& cb = plan.clusters[b - 1];
        if (ca.local_leaves != cb.local_leaves) return ca.local_leaves > cb.local_leaves;
        return ca.id < cb.id;
    });
    return plan;
}

void dump_cluster_plan(const ClusterPlan& plan, std::ostream& out) {
    out << "rs_cluster_plan proper_clusters=" << plan.clusters.size()
        << " residual_id=" << plan.residual_id
        << " residual_leaves=" << plan.residual_leaves
        << " order_count=" << plan.solve_order.size() << '\n';
    for (int id : plan.solve_order) {
        const PlannedCluster& cluster = plan.clusters[id - 1];
        out << "rs_cluster_plan_item id=" << cluster.id
            << " original_id=" << cluster.original_id
            << " local_leaves=" << cluster.local_leaves
            << " original_leaves=" << cluster.original_leaves
            << " parent=" << cluster.parent << '\n';
    }
}

namespace {

NodeId root_for_cluster(const ClusterPlan& plan, int cluster_id, bool first_tree) {
    if (cluster_id <= 0 || cluster_id > static_cast<int>(plan.clusters.size())) {
        return -1;
    }
    const PlannedCluster& cluster = plan.clusters[cluster_id - 1];
    return first_tree ? cluster.node1 : cluster.node2;
}

int replacement_cluster_at(const ClusterPlan& plan, NodeId node, int parent_cluster,
                           bool first_tree) {
    for (const PlannedCluster& cluster : plan.clusters) {
        if (cluster.parent != parent_cluster) continue;
        NodeId root = first_tree ? cluster.node1 : cluster.node2;
        if (root == node) return cluster.id;
    }
    return 0;
}

std::string leaf_atom_name(Leaf leaf, const std::vector<int>& original_to_internal) {
    if (leaf > 0 && leaf < static_cast<int>(original_to_internal.size()) &&
        original_to_internal[leaf] >= 0) {
        return std::to_string(original_to_internal[leaf]);
    }
    return std::to_string(leaf);
}

std::string abstract_subtree(const Tree& tree, const ClusterPlan& plan,
                             NodeId node, int parent_cluster, bool first_tree,
                             const std::vector<int>* original_to_internal = nullptr) {
    int replacement = replacement_cluster_at(plan, node, parent_cluster, first_tree);
    if (replacement != 0) return "C" + std::to_string(replacement);
    if (tree.is_leaf(node)) {
        return original_to_internal == nullptr
                   ? std::to_string(tree.label[node])
                   : leaf_atom_name(tree.label[node], *original_to_internal);
    }

    std::vector<std::string> children;
    std::vector<std::string> appended_clusters;
    for (NodeId child : {tree.left[node], tree.right[node]}) {
        int child_replacement =
            replacement_cluster_at(plan, child, parent_cluster, first_tree);
        if (child_replacement != 0) {
            appended_clusters.push_back("C" + std::to_string(child_replacement));
            continue;
        }
        children.push_back(abstract_subtree(tree, plan, child, parent_cluster,
                                            first_tree, original_to_internal));
    }
    children.insert(children.end(), appended_clusters.begin(),
                    appended_clusters.end());
    std::string out = "(";
    for (size_t i = 0; i < children.size(); ++i) {
        if (i) out.push_back(',');
        out += children[i];
    }
    out.push_back(')');
    return out;
}

bool is_open_cluster(int cluster_id, const std::vector<int>& open_children) {
    return std::find(open_children.begin(), open_children.end(), cluster_id) !=
           open_children.end();
}

bool is_suppressed_cluster(int cluster_id,
                           const std::vector<unsigned char>& suppressed_clusters) {
    return cluster_id > 0 &&
           cluster_id < static_cast<int>(suppressed_clusters.size()) &&
           suppressed_clusters[cluster_id] != 0;
}

std::string abstract_subtree_open(const Tree& tree, const ClusterPlan& plan,
                                  NodeId node, int parent_cluster,
                                  bool first_tree,
                                  const std::vector<int>& open_children,
                                  const std::vector<unsigned char>& suppressed_clusters,
                                  const std::vector<int>* original_to_internal) {
    int effective_parent = parent_cluster;
    int replacement = replacement_cluster_at(plan, node, parent_cluster, first_tree);
    if (replacement != 0) {
        if (is_suppressed_cluster(replacement, suppressed_clusters)) return "";
        if (!is_open_cluster(replacement, open_children)) {
            return "C" + std::to_string(replacement);
        }
        effective_parent = replacement;
    }

    if (tree.is_leaf(node)) {
        return original_to_internal == nullptr
                   ? std::to_string(tree.label[node])
                   : leaf_atom_name(tree.label[node], *original_to_internal);
    }

    std::vector<std::string> children;
    std::vector<std::string> appended_clusters;
    for (NodeId child : {tree.left[node], tree.right[node]}) {
        int child_replacement =
            replacement_cluster_at(plan, child, effective_parent, first_tree);
        if (is_suppressed_cluster(child_replacement, suppressed_clusters)) {
            continue;
        }
        if (child_replacement != 0 &&
            !is_open_cluster(child_replacement, open_children)) {
            appended_clusters.push_back("C" + std::to_string(child_replacement));
            continue;
        }
        children.push_back(abstract_subtree_open(
            tree, plan, child, effective_parent, first_tree, open_children,
            suppressed_clusters, original_to_internal));
    }
    children.insert(children.end(), appended_clusters.begin(),
                    appended_clusters.end());
    children.erase(std::remove(children.begin(), children.end(), ""),
                   children.end());
    if (children.empty()) return "";
    if (children.size() == 1) return children.front();
    std::string out = "(";
    for (size_t i = 0; i < children.size(); ++i) {
        if (i) out.push_back(',');
        out += children[i];
    }
    out.push_back(')');
    return out;
}

void assign_internal_labels_preorder(const Tree& tree, NodeId node,
                                     std::vector<int>& original_to_internal,
                                     std::vector<int>& internal_to_original) {
    if (node < 0) return;
    if (tree.is_leaf(node)) {
        Leaf leaf = tree.label[node];
        if (leaf > 0 && leaf < static_cast<int>(original_to_internal.size()) &&
            original_to_internal[leaf] < 0) {
            int internal = static_cast<int>(internal_to_original.size());
            original_to_internal[leaf] = internal;
            internal_to_original.push_back(leaf);
        }
        return;
    }
    assign_internal_labels_preorder(tree, tree.left[node],
                                    original_to_internal, internal_to_original);
    assign_internal_labels_preorder(tree, tree.right[node],
                                    original_to_internal, internal_to_original);
}

void h67_label_maps(const Tree& t1, std::vector<int>& original_to_internal,
                    std::vector<int>& internal_to_original) {
    original_to_internal.assign(t1.n + 1, -1);
    internal_to_original.clear();
    internal_to_original.reserve(t1.n);
    assign_internal_labels_preorder(t1, t1.root, original_to_internal,
                                    internal_to_original);
}

}  // namespace

std::string abstract_newick(const Tree& tree, const ClusterPlan& plan,
                            int cluster_id, bool first_tree) {
    NodeId root = cluster_id == 0 ? tree.root
                                  : root_for_cluster(plan, cluster_id, first_tree);
    if (root < 0) return ";";
    return abstract_subtree(tree, plan, root, cluster_id, first_tree) + ";";
}

void dump_cluster_newicks(const Tree& t1, const Tree& t2,
                          const ClusterPlan& plan, std::ostream& out) {
    out << "rs_cluster_residual_first=" << abstract_newick(t1, plan, 0, true)
        << '\n';
    out << "rs_cluster_residual_second=" << abstract_newick(t2, plan, 0, false)
        << '\n';
    for (const PlannedCluster& cluster : plan.clusters) {
        out << "rs_cluster_local id=" << cluster.id
            << " first=" << abstract_newick(t1, plan, cluster.id, true)
            << " second=" << abstract_newick(t2, plan, cluster.id, false)
            << '\n';
    }
}

namespace {

std::string strip_semicolon(std::string text) {
    if (!text.empty() && text.back() == ';') text.pop_back();
    return text;
}

std::string numeric_atom_newick(const std::string& abstract_newick,
                                std::map<std::string, int>& atom_id,
                                std::vector<std::string>& atom_names) {
    std::string out;
    for (size_t i = 0; i < abstract_newick.size();) {
        char ch = abstract_newick[i];
        if (ch == 'C') {
            size_t j = i + 1;
            while (j < abstract_newick.size() &&
                   std::isdigit(static_cast<unsigned char>(abstract_newick[j]))) {
                ++j;
            }
            std::string atom = abstract_newick.substr(i, j - i);
            auto [it, inserted] =
                atom_id.insert({atom, static_cast<int>(atom_id.size()) + 1});
            if (inserted) atom_names.push_back(atom);
            out += std::to_string(it->second);
            i = j;
        } else if (std::isdigit(static_cast<unsigned char>(ch))) {
            size_t j = i + 1;
            while (j < abstract_newick.size() &&
                   std::isdigit(static_cast<unsigned char>(abstract_newick[j]))) {
                ++j;
            }
            std::string atom = abstract_newick.substr(i, j - i);
            auto [it, inserted] =
                atom_id.insert({atom, static_cast<int>(atom_id.size()) + 1});
            if (inserted) atom_names.push_back(atom);
            out += std::to_string(it->second);
            i = j;
        } else {
            out.push_back(ch);
            ++i;
        }
    }
    return out;
}

std::vector<std::string> h67_atom_order_open(
    const Tree& tree, const ClusterPlan& plan, int cluster_id, bool first_tree,
    const std::vector<int>& original_to_internal,
    const std::vector<int>& open_children,
    const std::vector<unsigned char>& suppressed_clusters) {
    std::vector<std::string> atoms;
    NodeId root = cluster_id == 0 ? tree.root
                                  : root_for_cluster(plan, cluster_id, first_tree);
    if (root < 0) return atoms;

    std::vector<std::pair<NodeId, int>> stack{{root, cluster_id}};
    while (!stack.empty()) {
        auto [node, parent_cluster] = stack.back();
        stack.pop_back();

        int effective_parent = parent_cluster;
        int replacement =
            replacement_cluster_at(plan, node, parent_cluster, first_tree);
        if (replacement != 0) {
            if (is_suppressed_cluster(replacement, suppressed_clusters)) {
                continue;
            }
            if (!is_open_cluster(replacement, open_children)) {
                atoms.push_back("C" + std::to_string(replacement));
                continue;
            }
            effective_parent = replacement;
        }

        if (tree.is_leaf(node)) {
            atoms.push_back(leaf_atom_name(tree.label[node], original_to_internal));
            continue;
        }
        std::vector<NodeId> ordinary;
        std::vector<std::string> appended_clusters;
        for (NodeId child : {tree.left[node], tree.right[node]}) {
            int child_replacement =
                replacement_cluster_at(plan, child, effective_parent, first_tree);
            if (is_suppressed_cluster(child_replacement, suppressed_clusters)) {
                continue;
            }
            if (child_replacement != 0 &&
                !is_open_cluster(child_replacement, open_children)) {
                appended_clusters.push_back("C" + std::to_string(child_replacement));
            } else {
                ordinary.push_back(child);
            }
        }
        for (NodeId child : ordinary) stack.push_back({child, effective_parent});
        for (auto it = appended_clusters.rbegin(); it != appended_clusters.rend();
             ++it) {
            atoms.push_back(*it);
        }
    }
    return atoms;
}

}  // namespace

LocalProblem build_local_problem(const Tree& t1, const Tree& t2,
                                 const ClusterPlan& plan, int cluster_id,
                                 bool attach_boundary_marker) {
    return build_local_problem_with_open_children(
        t1, t2, plan, cluster_id, attach_boundary_marker, {});
}

LocalProblem build_local_problem_with_open_children(
    const Tree& t1, const Tree& t2, const ClusterPlan& plan, int cluster_id,
    bool attach_boundary_marker, const std::vector<int>& open_children) {
    return build_local_problem_with_open_children_and_suppressed(
        t1, t2, plan, cluster_id, attach_boundary_marker, open_children, {});
}

LocalProblem build_local_problem_with_open_children_and_suppressed(
    const Tree& t1, const Tree& t2, const ClusterPlan& plan, int cluster_id,
    bool attach_boundary_marker, const std::vector<int>& open_children,
    const std::vector<unsigned char>& suppressed_clusters) {
    std::map<std::string, int> atom_id;
    LocalProblem problem;
    problem.atom_names.push_back("");
    std::vector<int> original_to_internal;
    h67_label_maps(t1, original_to_internal, problem.internal_to_original);

    for (const std::string& atom :
         h67_atom_order_open(t1, plan, cluster_id, true, original_to_internal,
                             open_children, suppressed_clusters)) {
        if (atom_id.find(atom) != atom_id.end()) continue;
        atom_id.insert({atom, static_cast<int>(atom_id.size()) + 1});
        problem.atom_names.push_back(atom);
    }

    NodeId first_root = cluster_id == 0 ? t1.root
                                        : root_for_cluster(plan, cluster_id, true);
    NodeId second_root = cluster_id == 0 ? t2.root
                                         : root_for_cluster(plan, cluster_id, false);
    std::string first =
        first_root < 0
            ? ";"
            : abstract_subtree_open(t1, plan, first_root, cluster_id, true,
                                    open_children, suppressed_clusters,
                                    &original_to_internal) +
                  ";";
    std::string second =
        second_root < 0
            ? ";"
            : abstract_subtree_open(t2, plan, second_root, cluster_id, false,
                                    open_children, suppressed_clusters,
                                    &original_to_internal) +
                  ";";
    problem.first_newick = numeric_atom_newick(first, atom_id, problem.atom_names);
    problem.second_newick = numeric_atom_newick(second, atom_id, problem.atom_names);
    problem.solve_n = static_cast<int>(problem.atom_names.size()) - 1;

    if (attach_boundary_marker) {
        int boundary = problem.solve_n + 1;
        problem.first_newick = "(" + strip_semicolon(problem.first_newick) +
                               "," + std::to_string(boundary) + ");";
        problem.second_newick = "(" + strip_semicolon(problem.second_newick) +
                                "," + std::to_string(boundary) + ");";
        problem.atom_names.push_back("__BOUNDARY__");
        problem.solve_n = boundary;
    }
    return problem;
}

void dump_local_problem(const LocalProblem& problem, int cluster_id,
                        std::ostream& out) {
    out << "rs_local_problem id=" << cluster_id
        << " solve_n=" << problem.solve_n
        << " atoms=" << problem.atom_names.size() - 1 << '\n';
    out << "rs_local_first id=" << cluster_id << " newick="
        << problem.first_newick << '\n';
    out << "rs_local_second id=" << cluster_id << " newick="
        << problem.second_newick << '\n';
    out << "rs_local_atom_map id=" << cluster_id << " map=";
    for (size_t i = 1; i < problem.atom_names.size(); ++i) {
        if (i > 1) out << ' ';
        out << i << ':' << problem.atom_names[i];
    }
    out << '\n';
}
