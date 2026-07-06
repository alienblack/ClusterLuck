#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "tree/tree.h"

struct ClusterCandidate {
    int original_id = 0;
    NodeId node1 = -1;
    NodeId twin2 = -1;
    bool skipped = false;
    const char* skip_rule = nullptr;

    int child_first = 0;
    int parent_first = 0;
    int sibling_first = 0;
    int child_second = 0;
    int parent_second = 0;
    int sibling_second = 0;

    int original_child_first = 0;
};

struct PlannedCluster {
    int id = 0;
    int original_id = 0;
    NodeId node1 = -1;
    NodeId node2 = -1;
    int local_leaves = 0;
    int original_leaves = 0;
    int parent = 0;
};

struct ClusterPlan {
    std::vector<PlannedCluster> clusters;
    int residual_id = 0;
    int residual_leaves = 0;
    std::vector<int> solve_order;
};

struct LocalProblem {
    std::string first_newick;
    std::string second_newick;
    std::vector<std::string> atom_names;
    std::vector<int> internal_to_original;
    int solve_n = 0;
};

std::vector<ClusterCandidate> find_cluster_candidates(const Tree& t1,
                                                      const Tree& t2);

ClusterPlan build_cluster_plan(const Tree& t1,
                               const std::vector<ClusterCandidate>& candidates);

void dump_cluster_candidates(const std::vector<ClusterCandidate>& candidates,
                             std::ostream& out);

void dump_cluster_plan(const ClusterPlan& plan, std::ostream& out);

std::string abstract_newick(const Tree& tree, const ClusterPlan& plan,
                            int cluster_id, bool first_tree);

void dump_cluster_newicks(const Tree& t1, const Tree& t2,
                          const ClusterPlan& plan, std::ostream& out);

LocalProblem build_local_problem(const Tree& t1, const Tree& t2,
                                 const ClusterPlan& plan, int cluster_id,
                                 bool attach_boundary_marker);

LocalProblem build_local_problem_with_open_children(
    const Tree& t1, const Tree& t2, const ClusterPlan& plan, int cluster_id,
    bool attach_boundary_marker, const std::vector<int>& open_children);

LocalProblem build_local_problem_with_open_children_and_suppressed(
    const Tree& t1, const Tree& t2, const ClusterPlan& plan, int cluster_id,
    bool attach_boundary_marker, const std::vector<int>& open_children,
    const std::vector<unsigned char>& suppressed_clusters);

void dump_local_problem(const LocalProblem& problem, int cluster_id,
                        std::ostream& out);
