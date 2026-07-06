#include <exception>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cluster/cluster_points.h"
#include "cluster/cluster_solver.h"
#include "forest/forest.h"
#include "forest/validator.h"
#include "heuristic/greedy12.h"
#include "heuristic/greedy_merge.h"
#include "heuristic/highs_pool_repair.h"
#include "heuristic/near_cluster_route.h"
#include "heuristic/structured_cut_repair.h"
#include "io/pace_io.h"
#include "kernel/kernel_reduce.h"
#include "util/signal_handler.h"

namespace {

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    char* end = nullptr;
    long value = std::strtol(raw, &end, 10);
    if (end == raw || *end != '\0') return fallback;
    return static_cast<int>(value);
}

void set_env_default(const char* name, const char* value) {
    const char* raw = std::getenv(name);
    if (raw != nullptr && *raw != '\0') return;
    setenv(name, value, 0);
}

void set_env_default_int(const char* name, int value) {
    std::ostringstream out;
    out << value;
    set_env_default(name, out.str().c_str());
}

void set_env_force(const char* name, const char* value) {
    setenv(name, value, 1);
}

double average_leaf_depth(const Tree& tree) {
    if (tree.n <= 0) return 0.0;
    long long total = 0;
    for (int leaf = 1; leaf <= tree.n; ++leaf) {
        NodeId node = tree.leaf_node[leaf];
        if (node >= 0 && node < static_cast<NodeId>(tree.depth.size())) {
            total += tree.depth[node];
        }
    }
    return static_cast<double>(total) / static_cast<double>(tree.n);
}

int max_cluster_leaves(const ClusterPlan& plan) {
    int largest = 0;
    for (const PlannedCluster& cluster : plan.clusters) {
        largest = std::max(largest, cluster.original_leaves);
    }
    return largest;
}

void keep_and_publish_if_better(Forest candidate, const Instance& instance,
                                Forest& best);
Forest run_fast_restart_chain(Forest best, const Instance& instance);
Forest run_direct_adaptive_chain(Forest best, const Instance& instance);
Forest run_no_cluster_fallback(const Instance& instance);
Forest run_final_anytime_tail(
    Forest best, const Instance& instance,
    std::chrono::steady_clock::time_point program_start);

struct PublishGuard {
    explicit PublishGuard(bool enabled) {
        set_publish_enabled(enabled);
    }
    ~PublishGuard() {
        set_publish_enabled(true);
    }
};

bool should_run_big_101_route(const Instance& instance,
                              const ClusterPlan& plan) {
    if (env_int("RS_AUTO_101_ROUTE", 1) == 0) return false;
    if (instance.n <= env_int("RS_101_MIN_LEAVES", 500)) return false;
    return plan.clusters.empty() && plan.residual_leaves == instance.n;
}

bool should_run_102_route(const Instance& instance, const ClusterPlan& plan,
                          int& largest_cluster_out) {
    if (env_int("RS_AUTO_102_ROUTE", 1) == 0) return false;
    if (plan.clusters.empty()) return false;
    largest_cluster_out = max_cluster_leaves(plan);
    if (largest_cluster_out > env_int("RS_102_MAX_CLUSTER_LEAVES", 25)) {
        return false;
    }
    const int min_residue_percent = env_int("RS_102_MIN_RESIDUE_PERCENT", 75);
    return plan.residual_leaves * 100LL >=
           static_cast<long long>(instance.n) * min_residue_percent;
}

void configure_big_101_defaults() {
    set_env_default("RS_EARLY_GREEDY_BASELINE", "0");
    set_env_default("RS_DIRECT_ADAPTIVE_CHAIN", "1");
    set_env_default("RS_DIRECT_ADAPTIVE_MIN_LEAVES", "800");
    set_env_default("RS_DIRECT_ADAPTIVE_MIN_COMPONENTS", "600");
    set_env_default("RS_DIRECT_ADAPTIVE_INITIAL_BUDGET", "17");
    set_env_default("RS_DIRECT_ADAPTIVE_SECONDS", "220");
    set_env_default("RS_DIRECT_ADAPTIVE_STEPS", "16");
    set_env_default("RS_DIRECT_ADAPTIVE_STALE_LIMIT", "3");
    set_env_default("RS_DIRECT_PORTFOLIO", "1");
    set_env_default("RS_DIRECT_SPLIT_REPACK", "0");
    set_env_default("RS_FAST_RESTART_CHAIN", "0");
}

void configure_101_mode2_defaults() {
    configure_big_101_defaults();
    set_env_default("RS_DIRECT_ADAPTIVE_SCORE_MODE", "2");
}

void configure_102_defaults(const Instance& instance) {
    set_env_default("RS_DIRECT_ADAPTIVE_CHAIN", "0");
    set_env_default("RS_DIRECT_RESTART_BUDGET", "128");
    set_env_default("RS_DIRECT_RESTART_SCORE_MODE", "14");
    set_env_default("RS_DIRECT_EXTRA_SPLITS", "1");
    set_env_default("RS_DIRECT_PORTFOLIO", "0");
    set_env_default("RS_DIRECT_SPLIT_REPACK", "0");

    set_env_default("RS_NEAR_CLUSTER_ROUTE", "1");
    set_env_default_int("RS_NEAR_CLUSTER_MAX_LEAVES", instance.n);
    set_env_default_int("RS_CONFLICT_POLLUTER_MAX_LEAVES", instance.n);
    set_env_default_int("RS_NEAR_CLADE_POLLUTER_MAX_LEAVES", instance.n);

    set_env_default("RS_CONFLICT_POLLUTER_ROUTE", "1");
    set_env_default("RS_CONFLICT_POLLUTER_PAIR_LIMIT", "400");
    set_env_default("RS_CONFLICT_POLLUTER_MAX_REGION_COMPONENTS", "150");
    set_env_default("RS_CONFLICT_POLLUTER_MAX_REGION_LEAVES", "600");
    set_env_default("RS_CONFLICT_POLLUTER_LOCAL_RESTART_BUDGET", "48");
    set_env_default("RS_CONFLICT_POLLUTER_LOCAL_SLACK", "3");
    set_env_default("RS_CONFLICT_POLLUTER_FINAL_RESTART_BUDGET", "96");

    if (std::getenv("RS_NEAR_CLADE_POLLUTER_ROUTE") == nullptr) {
        const int max_soft_window_leaves =
            env_int("RS_102_SOFT_WINDOW_MAX_LEAVES", 2200);
        setenv("RS_NEAR_CLADE_POLLUTER_ROUTE",
               instance.n <= max_soft_window_leaves ? "1" : "0", 0);
    }
    set_env_default("RS_NEAR_CLADE_POLLUTER_MIN_COMMON", "5");
    set_env_default("RS_NEAR_CLADE_POLLUTER_MAX_SPILL", "14");
    set_env_default("RS_NEAR_CLADE_POLLUTER_CANDIDATES", "800");
    set_env_default("RS_NEAR_CLADE_POLLUTER_MAX_REGION_COMPONENTS", "22");
    set_env_default("RS_NEAR_CLADE_POLLUTER_MAX_REGION_LEAVES", "100");
    set_env_default("RS_NEAR_CLADE_POLLUTER_LOCAL_RESTART_BUDGET", "48");
    set_env_default("RS_NEAR_CLADE_POLLUTER_LOCAL_SLACK", "4");
    set_env_default("RS_NEAR_CLADE_POLLUTER_FINAL_RESTART_BUDGET", "96");

    set_env_default("RS_NEAR_CLUSTER_SEED_ROUTE", "0");
    set_env_default("RS_DIRECT_WINDOW_REPAIR", "1");
    set_env_default("RS_DIRECT_WINDOW_SIZE", "24");
    set_env_default("RS_DIRECT_WINDOW_CANDIDATES", "14");
    set_env_default("RS_DIRECT_WINDOW_RESTART_BUDGET", "256");
    set_env_default("RS_DIRECT_WINDOW_SECONDS", "8");
    set_env_default("RS_LOCAL_EXCHANGE_REPAIR", "1");
    set_env_default("RS_LOCAL_EXCHANGE_SECONDS", "8");

    set_env_default("RS_FAST_RESTART_CHAIN", "1");
    set_env_default("RS_FAST_RESTART_CHAIN_SECONDS", "165");
}

void configure_giant_fast_mode(int leaves) {
    if (env_int("RS_GIANT_FAST_MODE", 1) == 0) return;
    const int min_leaves = env_int("RS_GIANT_FAST_MIN_LEAVES", 9000);
    if (leaves < min_leaves) return;

    set_env_default("RS_EARLY_GREEDY_BASELINE", "0");
    set_env_default("RS_LOCAL_CLUSTER_RESTART1_BUDGET", "0");
    set_env_default("RS_LOCAL_CLUSTER_RESTART2_BUDGET", "0");
    set_env_default("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "0");
    set_env_default("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", "0");
    set_env_default("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", "0");
    set_env_default("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "0");
    set_env_default("RS_CLUSTER_OPEN_CHILD_VARIANTS", "0");
    set_env_default("RS_WHOLE_PROBE", "0");
    set_env_default("RS_FINAL_RESTART_BUDGET", "0");
}

enum class KernelClusterRoute {
    None,
    HugeClustered,
    ClusterRichNormal,
    Medium087,
    Medium068,
};

const char* kernel_cluster_route_name(KernelClusterRoute route) {
    switch (route) {
        case KernelClusterRoute::HugeClustered:
            return "huge_clustered";
        case KernelClusterRoute::ClusterRichNormal:
            return "cluster_rich_normal";
        case KernelClusterRoute::Medium087:
            return "medium_087";
        case KernelClusterRoute::Medium068:
            return "medium_068";
        case KernelClusterRoute::None:
            return "none";
    }
    return "none";
}

KernelClusterRoute classify_kernel_cluster_route(
    const Instance& reduced, const ClusterPlan& plan) {
    if (env_int("RS_KERNEL_CLUSTER_ROUTES", 1) == 0) {
        return KernelClusterRoute::None;
    }

    const int compact_n = reduced.n;
    const int proper_clusters = static_cast<int>(plan.clusters.size());
    const int residual_leaves = plan.residual_leaves;

    if (compact_n >= env_int("RS_KERNEL_HUGE_MIN_COMPACT", 6000) &&
        proper_clusters >= env_int("RS_KERNEL_HUGE_MIN_CLUSTERS", 500)) {
        return KernelClusterRoute::HugeClustered;
    }
    if (compact_n >= env_int("RS_KERNEL_087_MIN_COMPACT", 3500) &&
        proper_clusters >= env_int("RS_KERNEL_087_MIN_CLUSTERS", 150) &&
        residual_leaves <= env_int("RS_KERNEL_087_MAX_RESIDUAL", 500)) {
        return KernelClusterRoute::Medium087;
    }
    if (compact_n >= env_int("RS_KERNEL_RICH_MIN_COMPACT", 1500) &&
        proper_clusters >= env_int("RS_KERNEL_RICH_MIN_CLUSTERS", 50) &&
        residual_leaves * 100 <=
            compact_n * env_int("RS_KERNEL_RICH_MAX_RESIDUAL_PERCENT", 55)) {
        return KernelClusterRoute::ClusterRichNormal;
    }
    if (compact_n >= env_int("RS_KERNEL_068_MIN_COMPACT", 1800) &&
        compact_n < env_int("RS_KERNEL_068_MAX_COMPACT", 3500) &&
        proper_clusters >= env_int("RS_KERNEL_068_MIN_CLUSTERS", 30)) {
        return KernelClusterRoute::Medium068;
    }
    return KernelClusterRoute::None;
}

void configure_kernel_cluster_common() {
    set_env_force("RS_EARLY_GREEDY_BASELINE", "0");
    set_env_force("RS_WHOLE_PROBE", "0");
    set_env_force("RS_FINAL_RESTART_BUDGET", "0");
    set_env_force("RS_GREEDY12_RETURN", "0");
    set_env_force("RS_DIRECT_ADAPTIVE_CHAIN", "0");
    set_env_force("RS_FAST_RESTART_CHAIN", "0");
}

void configure_kernel_huge_clustered() {
    configure_kernel_cluster_common();
    set_env_force("RS_BOTTOMUP_CLEANUP", "1");
    set_env_force("RS_LOCAL_CHERRY_KERNEL", "1");
    set_env_force("RS_LOCAL_RESTART1_MODE", "2");
    set_env_force("RS_LOCAL_RESTART2_MODE", "2");
    set_env_force("RS_LOCAL_CLUSTER_RESTART1_BUDGET", "256");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", "256");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_BUDGET", "0");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", "0");
    set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANTS", "0");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART", "1");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_MIN_SOLVE_N", "1000");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_THRESHOLD", "25");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_BUDGET", "128");
}

void configure_kernel_cluster_rich_normal() {
    configure_kernel_cluster_common();
    set_env_force("RS_BOTTOMUP_CLEANUP", "1");
    set_env_force("RS_LOCAL_FULL_KERNEL", "1");
    set_env_force("RS_LOCAL_FULL_KERNEL_MIN_N", "2");
    set_env_force("RS_LOCAL_CHERRY_KERNEL", "0");
    set_env_force("RS_LOCAL_COMPACT_BUDGETS", "0");
    set_env_force("RS_LOCAL_RESTART1_MODE", "2");
    set_env_force("RS_LOCAL_RESTART2_MODE", "2");
    set_env_force("RS_LOCAL_CLUSTER_RESTART1_BUDGET", "256");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", "256");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_BUDGET", "0");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", "0");
    set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANTS", "0");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART", "0");
}

void configure_kernel_087_route() {
    configure_kernel_cluster_common();
    set_env_force("RS_BOTTOMUP_CLEANUP", "1");
    set_env_force("RS_LOCAL_CHERRY_KERNEL", "1");
    set_env_force("RS_LOCAL_RESTART1_MODE", "2");
    set_env_force("RS_LOCAL_RESTART2_MODE", "0");
    set_env_force("RS_LOCAL_CLUSTER_RESTART1_BUDGET", "750");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", "750");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_BUDGET", "128");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", "128");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_MIN_LEAVES", "80");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_MIN_LEAVES", "80");
    set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANTS", "0");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART", "1");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_MIN_SOLVE_N", "1000");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_THRESHOLD", "25");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_BUDGET", "128");
}

void configure_kernel_068_route() {
    configure_kernel_cluster_common();
    set_env_force("RS_BOTTOMUP_CLEANUP", "1");
    set_env_force("RS_LOCAL_CHERRY_KERNEL", "1");
    set_env_force("RS_LOCAL_RESTART1_MODE", "2");
    set_env_force("RS_LOCAL_RESTART2_MODE", "0");
    set_env_force("RS_LOCAL_CLUSTER_RESTART1_BUDGET", "750");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART1_BUDGET", "750");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_BUDGET", "128");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_BUDGET", "128");
    set_env_force("RS_LOCAL_CLUSTER_RESTART2_MIN_LEAVES", "80");
    set_env_force("RS_LOCAL_RESIDUAL_RESTART2_MIN_LEAVES", "80");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANTS", "1");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANT_LIMIT", "12");
    set_env_force("RS_CLUSTER_OPEN_CHILD_MAX_ATOMS", "120");
    if (env_int("RS_KERNEL_068_PORTFOLIO", 0) != 0) {
        set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "64");
        set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "128");
        set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_MIN_LEAVES", "80");
        set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_MIN_LEAVES", "80");
    } else {
        set_env_force("RS_LOCAL_CLUSTER_PORTFOLIO_BUDGET", "0");
        set_env_force("RS_LOCAL_RESIDUAL_PORTFOLIO_BUDGET", "0");
    }
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART", "1");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_MIN_SOLVE_N", "1000");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_THRESHOLD", "25");
    set_env_force("RS_CHILD_MARKER_LIGHT_RESTART_BUDGET", "128");
}

void configure_kernel_first_default_route() {
    configure_kernel_cluster_common();
    set_env_force("RS_BOTTOMUP_CLEANUP", "1");
    set_env_force("RS_LOCAL_FULL_KERNEL", "1");
    set_env_force("RS_LOCAL_FULL_KERNEL_MIN_N", "2");
    set_env_force("RS_LOCAL_CHERRY_KERNEL", "0");
    set_env_force("RS_LOCAL_COMPACT_BUDGETS", "1");
    set_env_force("RS_LOCAL_COMPACT_BUDGET_THRESHOLD", "150");
    set_env_force(
        "RS_LOCAL_COMPACT_BIG_RESTART1",
        std::to_string(env_int("RS_KERNEL_FIRST_BIG_RESTART1", 256)).c_str());
    set_env_force("RS_LOCAL_COMPACT_BIG_RESTART2", "0");
    set_env_force("RS_LOCAL_COMPACT_SMALL_RESTART1", "256");
    set_env_force("RS_LOCAL_COMPACT_SMALL_RESTART2", "256");
    set_env_force("RS_LOCAL_COMPACT_PORTFOLIO_BUDGET", "0");
    set_env_force("RS_LOCAL_RESTART1_MODE", "2");
    set_env_force("RS_LOCAL_RESTART2_MODE", "2");
    set_env_force("RS_CLUSTER_OPEN_CHILD_VARIANTS", "0");
    set_env_force("RS_WHOLE_PROBE", "0");
    set_env_force("RS_FINAL_RESTART_BUDGET", "0");
}

void configure_kernel_cluster_route(KernelClusterRoute route) {
    switch (route) {
        case KernelClusterRoute::HugeClustered:
            configure_kernel_huge_clustered();
            break;
        case KernelClusterRoute::ClusterRichNormal:
            configure_kernel_cluster_rich_normal();
            break;
        case KernelClusterRoute::Medium087:
            configure_kernel_087_route();
            break;
        case KernelClusterRoute::Medium068:
            configure_kernel_068_route();
            break;
        case KernelClusterRoute::None:
            break;
    }
}

Forest run_102_route(const Instance& instance, const ClusterPlan& plan,
                     Forest best) {
    configure_102_defaults(instance);
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_NEAR_CLUSTER_DEBUG") != nullptr) {
        std::cerr << "rs_102_route_start leaves=" << instance.n
                  << " avg_depth="
                  << 0.5 * (average_leaf_depth(instance.t1) +
                            average_leaf_depth(instance.t2))
                  << " proper_clusters=" << plan.clusters.size()
                  << " largest_cluster=" << max_cluster_leaves(plan)
                  << " residual_leaves=" << plan.residual_leaves << '\n';
    }

    if (should_run_greedy12_tail(instance.n, plan.residual_leaves)) {
        Forest greedy12;
        if (run_greedy12_tail(instance.t1, instance.t2, greedy12)) {
            keep_and_publish_if_better(std::move(greedy12), instance, best);
        }
    }

    Forest shallow = best;
    if (instance.n <= env_int("RS_102_FULL_FALLBACK_MAX_LEAVES", 2500)) {
        Forest fallback = run_no_cluster_fallback(instance);
        keep_and_publish_if_better(fallback, instance, best);
        shallow = best;
    }
    if (shallow.active_count == instance.n &&
        instance.n <= env_int("RS_102_DIRECT_GREEDY_MAX_LEAVES", 5000)) {
        shallow = run_h45_greedy_merge(
            std::move(shallow), instance.t1, instance.t2,
            env_int("RS_DIRECT_SCORE_MODE", 2), false, false);
        keep_and_publish_if_better(shallow, instance, best);
    }

    if (shallow.active_count < instance.n) {
        shallow = run_h45_split_restarts(
            shallow, instance.t1, instance.t2,
            env_int("RS_DIRECT_RESTART_SCORE_MODE", 14),
            env_int("RS_DIRECT_RESTART_BUDGET", 128),
            env_int("RS_DIRECT_RESTART_SEED_OFFSET", 0),
            env_int("RS_DIRECT_EXTRA_SPLITS", 1),
            true,
            env_int("RS_102_TARGET", std::max(1, shallow.active_count - 1)),
            static_cast<double>(env_int("RS_102_RESTART_SECONDS", 30)));
        keep_and_publish_if_better(shallow, instance, best);
    }

    shallow = run_near_cluster_route(
        shallow, instance.t1, instance.t2,
        env_int("RS_NEAR_CLUSTER_MAX_LEAVES", instance.n),
        env_int("RS_NEAR_CLUSTER_MIN_SIZE", 6),
        env_int("RS_NEAR_CLUSTER_MAX_SPILL", 8),
        env_int("RS_NEAR_CLUSTER_MAX_TRIM", 2),
        env_int("RS_NEAR_CLUSTER_MAX_CANDIDATES", 2000),
        env_int("RS_NEAR_CLUSTER_MAX_ISLANDS", 32),
        env_int("RS_NEAR_CLUSTER_RESTART_BUDGET", 256),
        env_int("RS_NEAR_CLUSTER_TARGET", std::max(1, best.active_count - 1)));
    keep_and_publish_if_better(std::move(shallow), instance, best);

    if (env_int("RS_FAST_RESTART_CHAIN", 1) != 0) {
        best = run_fast_restart_chain(std::move(best), instance);
    }

    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_NEAR_CLUSTER_DEBUG") != nullptr) {
        std::cerr << "rs_102_route_done components=" << best.active_count << '\n';
    }
    return best;
}

bool run_kernel_101_route(const Instance& original, Forest& best) {
    auto kernel = kernel_reduce_instance(
        original, std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr);
    if (!kernel.has_value()) return false;

    const int mode2_max_leaves = env_int("RS_101_MODE2_MAX_LEAVES", 1500);
    if (kernel->reduced.n <= mode2_max_leaves) {
        configure_101_mode2_defaults();
    } else {
        configure_big_101_defaults();
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_101_start original_leaves=" << original.n
                  << " reduced_leaves=" << kernel->reduced.n
                  << " mode="
                  << (kernel->reduced.n <= mode2_max_leaves ? "mode2" : "mixed")
                  << " subtree=" << kernel->subtree_deleted
                  << " chain=" << kernel->chain_deleted
                  << " three_two=" << kernel->param_reduction << '\n';
    }

    Forest reduced_best;
    {
        PublishGuard guard(false);
        reduced_best = run_no_cluster_fallback(kernel->reduced);
    }
    Forest expanded =
        expand_kernel_reduced_forest(reduced_best, kernel->map, original);
    if (expanded.active_count < best.active_count) {
        best = std::move(expanded);
        publish_best_output(forest_to_output(best, original.t1));
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_101_done reduced_components="
                  << reduced_best.active_count
                  << " expanded_components=" << best.active_count << '\n';
    }
    return true;
}

bool run_kernel_102_route(const Instance& original, Forest& best) {
    auto kernel = kernel_reduce_instance(
        original, std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr);
    if (!kernel.has_value()) return false;
    configure_kernel_first_default_route();

    std::vector<ClusterCandidate> reduced_candidates =
        find_cluster_candidates(kernel->reduced.t1, kernel->reduced.t2);
    ClusterPlan reduced_plan =
        build_cluster_plan(kernel->reduced.t1, reduced_candidates);

    Forest reduced_best = make_singleton_forest(kernel->reduced.n);
    if (!rebuild_forest_data(reduced_best, kernel->reduced.t1,
                             kernel->reduced.t2)) {
        throw std::runtime_error("kernel reduced singleton invalid");
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_102_start original_leaves=" << original.n
                  << " reduced_leaves=" << kernel->reduced.n
                  << " proper_clusters=" << reduced_plan.clusters.size()
                  << " residual_leaves=" << reduced_plan.residual_leaves
                  << " subtree=" << kernel->subtree_deleted
                  << " chain=" << kernel->chain_deleted
                  << " three_two=" << kernel->param_reduction << '\n';
    }

    {
        PublishGuard guard(false);
        reduced_best =
            run_102_route(kernel->reduced, reduced_plan, std::move(reduced_best));
    }
    Forest expanded =
        expand_kernel_reduced_forest(reduced_best, kernel->map, original);
    if (expanded.active_count < best.active_count) {
        best = std::move(expanded);
        publish_best_output(forest_to_output(best, original.t1));
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_102_done reduced_components="
                  << reduced_best.active_count
                  << " expanded_components=" << best.active_count << '\n';
    }
    return true;
}

bool run_kernel_cluster_route(const Instance& original, Forest& best) {
    if (env_int("RS_KERNEL_CLUSTER_ROUTES", 1) == 0) return false;
    if (original.n < env_int("RS_KERNEL_068_MIN_COMPACT", 1800)) return false;

    auto kernel = kernel_reduce_instance(
        original, std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr);
    if (!kernel.has_value()) return false;

    std::vector<ClusterCandidate> reduced_candidates =
        find_cluster_candidates(kernel->reduced.t1, kernel->reduced.t2);
    ClusterPlan reduced_plan =
        build_cluster_plan(kernel->reduced.t1, reduced_candidates);
    KernelClusterRoute route =
        classify_kernel_cluster_route(kernel->reduced, reduced_plan);
    if (route == KernelClusterRoute::None) return false;

    configure_kernel_cluster_route(route);
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_cluster_start kind="
                  << kernel_cluster_route_name(route)
                  << " original_leaves=" << original.n
                  << " reduced_leaves=" << kernel->reduced.n
                  << " proper_clusters=" << reduced_plan.clusters.size()
                  << " residual_leaves=" << reduced_plan.residual_leaves
                  << " largest_cluster=" << max_cluster_leaves(reduced_plan)
                  << " subtree=" << kernel->subtree_deleted
                  << " chain=" << kernel->chain_deleted
                  << " three_two=" << kernel->param_reduction << '\n';
    }

    Forest reduced_best;
    {
        PublishGuard guard(false);
        if (!run_cluster_solver(kernel->reduced.t1, kernel->reduced.t2,
                                reduced_best)) {
            return false;
        }
    }
    reduced_best = exact_output_repair(std::move(reduced_best),
                                       kernel->reduced.t1,
                                       kernel->reduced.t2);
    ValidationResult reduced_check =
        validate_forest(reduced_best, kernel->reduced.t1, kernel->reduced.t2);
    if (!reduced_check.ok) {
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
            std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
            std::cerr << "rs_kernel_cluster_invalid kind="
                      << kernel_cluster_route_name(route)
                      << " message=" << reduced_check.message << '\n';
        }
        return false;
    }

    Forest expanded =
        expand_kernel_reduced_forest(reduced_best, kernel->map, original);
    expanded = exact_output_repair(std::move(expanded), original.t1, original.t2);
    ValidationResult expanded_check =
        validate_forest(expanded, original.t1, original.t2);
    if (!expanded_check.ok) {
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
            std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
            std::cerr << "rs_kernel_cluster_expand_invalid kind="
                      << kernel_cluster_route_name(route)
                      << " message=" << expanded_check.message << '\n';
        }
        return false;
    }
    if (expanded.active_count < best.active_count) {
        best = std::move(expanded);
        publish_best_output(forest_to_output(best, original.t1));
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_cluster_done kind="
                  << kernel_cluster_route_name(route)
                  << " reduced_components=" << reduced_best.active_count
                  << " expanded_components=" << best.active_count << '\n';
    }
    return true;
}

bool run_kernel_first_default_route(const Instance& original, Forest& best) {
    if (env_int("RS_KERNEL_FIRST_ROUTE", 1) == 0) return false;
    if (original.n <= env_int("RS_KERNEL_FIRST_MIN_LEAVES", 2000)) {
        return false;
    }

    auto kernel = kernel_reduce_instance(
        original, std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr);
    if (!kernel.has_value()) return false;

    std::vector<ClusterCandidate> reduced_candidates =
        find_cluster_candidates(kernel->reduced.t1, kernel->reduced.t2);
    ClusterPlan reduced_plan =
        build_cluster_plan(kernel->reduced.t1, reduced_candidates);

    configure_kernel_first_default_route();

    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_first_start original_leaves=" << original.n
                  << " reduced_leaves=" << kernel->reduced.n
                  << " proper_clusters=" << reduced_plan.clusters.size()
                  << " residual_leaves=" << reduced_plan.residual_leaves
                  << " largest_cluster=" << max_cluster_leaves(reduced_plan)
                  << " subtree=" << kernel->subtree_deleted
                  << " chain=" << kernel->chain_deleted
                  << " three_two=" << kernel->param_reduction << '\n';
    }

    Forest reduced_best;
    {
        PublishGuard guard(false);
        if (reduced_plan.clusters.empty()) {
            reduced_best = run_no_cluster_fallback(kernel->reduced);
        } else if (!run_cluster_solver(kernel->reduced.t1, kernel->reduced.t2,
                                       reduced_best)) {
            reduced_best = run_no_cluster_fallback(kernel->reduced);
        }
    }
    reduced_best = exact_output_repair(std::move(reduced_best),
                                       kernel->reduced.t1,
                                       kernel->reduced.t2);
    ValidationResult reduced_check =
        validate_forest(reduced_best, kernel->reduced.t1, kernel->reduced.t2);
    if (!reduced_check.ok) {
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
            std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
            std::cerr << "rs_kernel_first_invalid message="
                      << reduced_check.message << '\n';
        }
        return false;
    }

    Forest expanded =
        expand_kernel_reduced_forest(reduced_best, kernel->map, original);
    expanded = exact_output_repair(std::move(expanded), original.t1, original.t2);
    ValidationResult expanded_check =
        validate_forest(expanded, original.t1, original.t2);
    if (!expanded_check.ok) {
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
            std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
            std::cerr << "rs_kernel_first_expand_invalid message="
                      << expanded_check.message << '\n';
        }
        return false;
    }
    if (expanded.active_count < best.active_count) {
        best = std::move(expanded);
        publish_best_output(forest_to_output(best, original.t1));
    }
    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr ||
        std::getenv("RS_DEBUG_KERNEL_ROUTE") != nullptr) {
        std::cerr << "rs_kernel_first_done reduced_components="
                  << reduced_best.active_count
                  << " expanded_components=" << best.active_count << '\n';
    }
    return true;
}

std::vector<int> parse_leaves_from_line(const std::string& line) {
    std::vector<int> leaves;
    int value = 0;
    bool in_number = false;
    for (char ch : line) {
        if (ch >= '0' && ch <= '9') {
            value = value * 10 + (ch - '0');
            in_number = true;
        } else if (in_number) {
            leaves.push_back(value);
            value = 0;
            in_number = false;
        }
    }
    if (in_number) leaves.push_back(value);
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    return leaves;
}

bool load_seed_forest(const char* path, int n, Forest& forest) {
    std::ifstream input(path);
    if (!input) return false;
    forest = Forest{};
    forest.leaf_comp.assign(n + 1, -1);
    std::string line;
    while (std::getline(input, line)) {
        std::vector<int> leaves = parse_leaves_from_line(line);
        if (leaves.empty()) continue;
        Component component;
        component.leaves.reserve(leaves.size());
        for (int leaf : leaves) {
            if (leaf < 1 || leaf > n || forest.leaf_comp[leaf] != -1) {
                return false;
            }
            forest.leaf_comp[leaf] = static_cast<CompId>(forest.comps.size());
            component.leaves.push_back(leaf);
        }
        component.repr = component.leaves.front();
        forest.comps.push_back(std::move(component));
        ++forest.active_count;
    }
    for (int leaf = 1; leaf <= n; ++leaf) {
        if (forest.leaf_comp[leaf] == -1) return false;
    }
    return true;
}

void keep_if_better(Forest candidate, const Instance& instance, Forest& best) {
    candidate = exact_output_repair(std::move(candidate), instance.t1, instance.t2);
    ValidationResult check = validate_forest(candidate, instance.t1, instance.t2);
    if (check.ok && candidate.active_count < best.active_count) {
        best = std::move(candidate);
    }
}

void keep_and_publish_if_better(Forest candidate, const Instance& instance,
                                Forest& best) {
    candidate = exact_output_repair(std::move(candidate), instance.t1, instance.t2);
    ValidationResult check = validate_forest(candidate, instance.t1, instance.t2);
    if (check.ok && candidate.active_count < best.active_count) {
        best = std::move(candidate);
        publish_best_output(forest_to_output(best, instance.t1));
    }
}

Forest run_fast_restart_chain(Forest best, const Instance& instance) {
    if (env_int("RS_FAST_RESTART_CHAIN", 0) == 0) return best;

    set_env_default(
        "RS_H45_MAX_MERGE_CHECKS",
        std::to_string(env_int("RS_FAST_RESTART_MAX_CHECKS", 8000)).c_str());

    const int total_seconds = env_int("RS_FAST_RESTART_CHAIN_SECONDS", 165);
    const auto start_time = std::chrono::steady_clock::now();
    auto remaining_seconds = [&]() {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        return std::max(0.0, static_cast<double>(total_seconds) - elapsed);
    };

    auto run_stage = [&](const char* label, int score_mode, int seconds,
                         int attempts, int seed_offset, int extra_splits,
                         int target_drop, bool use_fpt_conflict) {
        double remaining = remaining_seconds();
        double budget_seconds = std::min<double>(seconds, remaining);
        if (budget_seconds < 1.0) return;
        int target = std::max(1, best.active_count - target_drop);
        Forest candidate = run_h45_split_restarts(
            best, instance.t1, instance.t2, score_mode, attempts, seed_offset,
            extra_splits, use_fpt_conflict, target, budget_seconds);
        candidate =
            exact_output_repair(std::move(candidate), instance.t1, instance.t2);
        ValidationResult check =
            validate_forest(candidate, instance.t1, instance.t2);
        if (check.ok && candidate.active_count < best.active_count) {
            best = std::move(candidate);
            publish_best_output(forest_to_output(best, instance.t1));
        }
        if (std::getenv("RS_FAST_RESTART_DEBUG") != nullptr) {
            std::cerr << "fast_restart_stage label=" << label
                      << " best=" << best.active_count
                      << " remaining=" << remaining_seconds() << '\n';
        }
    };

    run_stage("coarse_mode20", env_int("RS_FAST_STAGE1_MODE", 20),
              env_int("RS_FAST_STAGE1_SECONDS", 18),
              env_int("RS_FAST_STAGE1_ATTEMPTS", 3000),
              env_int("RS_FAST_STAGE1_SEED_OFFSET", 6700),
              env_int("RS_FAST_STAGE1_EXTRA_SPLITS", 4),
              env_int("RS_FAST_STAGE1_TARGET_DROP", 5),
              env_int("RS_FAST_STAGE1_FPT", 1) != 0);
    run_stage("bulk_mode2", env_int("RS_FAST_STAGE2_MODE", 2),
              env_int("RS_FAST_STAGE2_SECONDS", 75),
              env_int("RS_FAST_STAGE2_ATTEMPTS", 6000),
              env_int("RS_FAST_STAGE2_SEED_OFFSET", 13000),
              env_int("RS_FAST_STAGE2_EXTRA_SPLITS", 4),
              env_int("RS_FAST_STAGE2_TARGET_DROP", 15),
              env_int("RS_FAST_STAGE2_FPT", 0) != 0);
    run_stage("late_mode1", env_int("RS_FAST_STAGE3_MODE", 1),
              env_int("RS_FAST_STAGE3_SECONDS", 65),
              env_int("RS_FAST_STAGE3_ATTEMPTS", 5000),
              env_int("RS_FAST_STAGE3_SEED_OFFSET", 6700),
              env_int("RS_FAST_STAGE3_EXTRA_SPLITS", 4),
              env_int("RS_FAST_STAGE3_TARGET_DROP", 3),
              env_int("RS_FAST_STAGE3_FPT", 1) != 0);
    run_stage("tail_mode1", env_int("RS_FAST_STAGE4_MODE", 1),
              env_int("RS_FAST_STAGE4_SECONDS", 0),
              env_int("RS_FAST_STAGE4_ATTEMPTS", 5000),
              env_int("RS_FAST_STAGE4_SEED_OFFSET", 6700),
              env_int("RS_FAST_STAGE4_EXTRA_SPLITS", 4),
              env_int("RS_FAST_STAGE4_TARGET_DROP", 2),
              env_int("RS_FAST_STAGE4_FPT", 0) != 0);

    return best;
}

Forest run_direct_adaptive_chain(Forest best, const Instance& instance) {
    if (env_int("RS_DIRECT_ADAPTIVE_CHAIN", 1) == 0) return best;
    if (best.active_count <= env_int("RS_DIRECT_ADAPTIVE_MIN_COMPONENTS", 600)) {
        return best;
    }

    const int total_seconds = env_int("RS_DIRECT_ADAPTIVE_SECONDS", 220);
    const int max_steps = env_int("RS_DIRECT_ADAPTIVE_STEPS", 16);
    const auto start_time = std::chrono::steady_clock::now();
    auto remaining_seconds = [&]() {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        return std::max(0.0, static_cast<double>(total_seconds) - elapsed);
    };

    int stale = 0;
    for (int step = 1; step <= max_steps; ++step) {
        double remaining = remaining_seconds();
        if (remaining < 8.0) break;

        int extra_splits = 2;
        int seconds = 35;
        int attempts = 96;
        int score_mode = 0;
        if (step == 3) {
            extra_splits = 4;
            seconds = 30;
            attempts = 160;
        } else if (step <= 5 && step > 3) {
            extra_splits = 5;
            seconds = 30;
            attempts = 192;
            score_mode = 20;
        } else if (step <= 7 && step > 5) {
            extra_splits = 6;
            seconds = 30;
            attempts = 224;
            score_mode = 20;
        } else if (step <= 9 && step > 7) {
            extra_splits = 7;
            attempts = 256;
            score_mode = 20;
        } else if (step > 9) {
            extra_splits = 8 + std::min(3, step - 10);
            attempts = 288;
            score_mode = 20;
        }

        extra_splits = env_int("RS_DIRECT_ADAPTIVE_EXTRA_SPLITS", extra_splits);
        seconds = env_int("RS_DIRECT_ADAPTIVE_STAGE_SECONDS", seconds);
        attempts = env_int("RS_DIRECT_ADAPTIVE_ATTEMPTS", attempts);
        score_mode = env_int("RS_DIRECT_ADAPTIVE_SCORE_MODE", score_mode);

        double budget_seconds = std::min<double>(seconds, remaining);
        int target = std::max(1, best.active_count - 1);
        Forest candidate = run_h45_split_restarts(
            best, instance.t1, instance.t2, score_mode, attempts,
            env_int("RS_DIRECT_ADAPTIVE_SEED_OFFSET", 0), extra_splits,
            env_int("RS_DIRECT_ADAPTIVE_FPT", 0) != 0, target, budget_seconds);
        candidate =
            exact_output_repair(std::move(candidate), instance.t1, instance.t2);
        ValidationResult check =
            validate_forest(candidate, instance.t1, instance.t2);
        bool improved = check.ok && candidate.active_count < best.active_count;
        if (improved) {
            best = std::move(candidate);
            publish_best_output(forest_to_output(best, instance.t1));
            stale = 0;
        } else {
            ++stale;
        }
        if (std::getenv("RS_DIRECT_ADAPTIVE_DEBUG") != nullptr) {
            std::cerr << "rs_direct_adaptive_stage step=" << step
                      << " best=" << best.active_count
                      << " improved=" << (improved ? 1 : 0)
                      << " mode=" << score_mode
                      << " extra=" << extra_splits
                      << " remaining=" << remaining_seconds() << '\n';
        }
        if (stale >= env_int("RS_DIRECT_ADAPTIVE_STALE_LIMIT", 3)) break;
    }

    return best;
}

Forest run_no_cluster_fallback(const Instance& instance) {
    const int score_mode = env_int("RS_DIRECT_SCORE_MODE", 2);
    Forest singleton = make_singleton_forest(instance.n);
    if (!rebuild_forest_data(singleton, instance.t1, instance.t2)) {
        throw std::runtime_error("direct fallback singleton invalid");
    }

    Forest raw = run_h45_greedy_merge(std::move(singleton), instance.t1,
                                      instance.t2, score_mode, false, false);
    Forest best = exact_output_repair(raw, instance.t1, instance.t2);
    if (!validate_forest(best, instance.t1, instance.t2).ok) {
        throw std::runtime_error("direct fallback greedy invalid");
    }

    int direct_restart_budget = env_int("RS_DIRECT_RESTART_BUDGET", 2048);
    if (env_int("RS_DIRECT_ADAPTIVE_CHAIN", 1) != 0 &&
        instance.n >= env_int("RS_DIRECT_ADAPTIVE_MIN_LEAVES", 800)) {
        direct_restart_budget =
            env_int("RS_DIRECT_ADAPTIVE_INITIAL_BUDGET", 17);
    }
    Forest restarted = run_h45_split_restarts(
        best, instance.t1, instance.t2,
        env_int("RS_DIRECT_RESTART_SCORE_MODE", 0),
        direct_restart_budget,
        env_int("RS_DIRECT_RESTART_SEED_OFFSET", 0),
        env_int("RS_DIRECT_EXTRA_SPLITS", 1));
    restarted = exact_output_repair(std::move(restarted),
                                    instance.t1, instance.t2);
    ValidationResult restart_check =
        validate_forest(restarted, instance.t1, instance.t2);
    if (restart_check.ok && restarted.active_count < best.active_count) {
        best = std::move(restarted);
    }
    if (env_int("RS_DIRECT_ADAPTIVE_CHAIN", 1) != 0 &&
        instance.n >= env_int("RS_DIRECT_ADAPTIVE_MIN_LEAVES", 800)) {
        best = run_direct_adaptive_chain(std::move(best), instance);
    }

    if (env_int("RS_DIRECT_PORTFOLIO", 1) != 0 &&
        best.active_count <= env_int("RS_DIRECT_PORTFOLIO_MAX_COMPONENTS", 260)) {
        if (best.active_count > env_int("RS_DIRECT_SMALL_PORTFOLIO_MAX", 180)) {
            Forest mode1_seed = make_singleton_forest(instance.n);
            if (!rebuild_forest_data(mode1_seed, instance.t1, instance.t2)) {
                throw std::runtime_error("direct portfolio singleton invalid");
            }
            Forest mode1_raw =
                run_h45_greedy_merge(std::move(mode1_seed), instance.t1,
                                     instance.t2, 1, false, false);
            Forest mode1 = exact_output_repair(mode1_raw, instance.t1, instance.t2);
            mode1 = run_h45_split_restarts(
                mode1, instance.t1, instance.t2,
                env_int("RS_DIRECT_PORTFOLIO_MODE1_RESTART_MODE", 1),
                env_int("RS_DIRECT_PORTFOLIO_MODE1_BUDGET", 1024),
                env_int("RS_DIRECT_PORTFOLIO_MODE1_SEED_OFFSET", 0),
                env_int("RS_DIRECT_PORTFOLIO_MODE1_EXTRA_SPLITS", 1));
            keep_if_better(std::move(mode1), instance, best);
        } else {
            Forest extra2 = run_h45_split_restarts(
                best, instance.t1, instance.t2,
                env_int("RS_DIRECT_RESTART_SCORE_MODE", 0),
                env_int("RS_DIRECT_PORTFOLIO_EXTRA2_BUDGET", 2048),
                env_int("RS_DIRECT_PORTFOLIO_EXTRA2_SEED_OFFSET", 0),
                env_int("RS_DIRECT_PORTFOLIO_EXTRA2_SPLITS", 2));
            keep_if_better(std::move(extra2), instance, best);

            Forest blocker = run_h45_blocker_repair(
                best, instance.t1, instance.t2,
                env_int("RS_DIRECT_BLOCKER_MAX_COMPONENTS", 180),
                env_int("RS_DIRECT_BLOCKER_MAX_BLOCKERS", 3),
                env_int("RS_DIRECT_BLOCKER_CANDIDATES", 12),
                env_int("RS_DIRECT_BLOCKER_RESTART_BUDGET", 256),
                env_int("RS_DIRECT_BLOCKER_TARGET", std::max(1, best.active_count - 1)),
                static_cast<double>(env_int("RS_DIRECT_BLOCKER_SECONDS", 6)));
            keep_if_better(std::move(blocker), instance, best);
        }
        if (env_int("RS_DIRECT_WINDOW_REPAIR", 0) != 0) {
            Forest window = run_h45_window_repair(
                best, instance.t1, instance.t2,
                env_int("RS_DIRECT_WINDOW_MAX_COMPONENTS", 220),
                env_int("RS_DIRECT_WINDOW_CANDIDATES", 10),
                env_int("RS_DIRECT_WINDOW_SIZE", 18),
                env_int("RS_DIRECT_WINDOW_RESTART_BUDGET", 256),
                env_int("RS_DIRECT_WINDOW_TARGET",
                        std::max(1, best.active_count - 1)),
                static_cast<double>(env_int("RS_DIRECT_WINDOW_SECONDS", 10)));
            keep_if_better(std::move(window), instance, best);
        }
        if (env_int("RS_LOCAL_EXCHANGE_REPAIR", 0) != 0) {
            Forest exchange = run_h45_local_exchange_repair(
                best, instance.t1, instance.t2,
                env_int("RS_LOCAL_EXCHANGE_MAX_COMPONENTS", 260),
                env_int("RS_LOCAL_EXCHANGE_CANDIDATES", 16),
                env_int("RS_LOCAL_EXCHANGE_WINDOW", 24),
                env_int("RS_LOCAL_EXCHANGE_MAX_BLOCKERS", 8),
                env_int("RS_LOCAL_EXCHANGE_RESTART_BUDGET", 256),
                env_int("RS_LOCAL_EXCHANGE_TARGET",
                        std::max(1, best.active_count - 1)),
                static_cast<double>(env_int("RS_LOCAL_EXCHANGE_SECONDS", 12)));
            keep_if_better(std::move(exchange), instance, best);
        }
    }
    if (env_int("RS_DIRECT_SPLIT_REPACK", 0) != 0) {
        Forest repacked = run_h45_targeted_split_repack(
            best, instance.t1, instance.t2,
            env_int("RS_DIRECT_SPLIT_REPACK_MAX_COMPONENTS", 400),
            env_int("RS_DIRECT_SPLIT_REPACK_CANDIDATES", 16),
            env_int("RS_DIRECT_SPLIT_REPACK_GROUP", 2),
            env_int("RS_DIRECT_SPLIT_REPACK_RESTART_BUDGET", 128),
            env_int("RS_DIRECT_SPLIT_REPACK_PASSES", 2),
            env_int("RS_DIRECT_SPLIT_REPACK_TARGET",
                    std::max(1, best.active_count - 1)),
            static_cast<double>(env_int("RS_DIRECT_SPLIT_REPACK_SECONDS", 20)));
        keep_if_better(std::move(repacked), instance, best);
    }
    if (env_int("RS_NEAR_CLUSTER_ROUTE", 0) != 0) {
        Forest near_cluster = run_near_cluster_route(
            best, instance.t1, instance.t2,
            env_int("RS_NEAR_CLUSTER_MAX_LEAVES", 700),
            env_int("RS_NEAR_CLUSTER_MIN_SIZE", 6),
            env_int("RS_NEAR_CLUSTER_MAX_SPILL", 8),
            env_int("RS_NEAR_CLUSTER_MAX_TRIM", 2),
            env_int("RS_NEAR_CLUSTER_MAX_CANDIDATES", 2000),
            env_int("RS_NEAR_CLUSTER_MAX_ISLANDS", 32),
            env_int("RS_NEAR_CLUSTER_RESTART_BUDGET", 256),
            env_int("RS_NEAR_CLUSTER_TARGET",
                    std::max(1, best.active_count - 1)));
        keep_if_better(std::move(near_cluster), instance, best);
    }
    return best;
}

Forest run_final_anytime_tail(
    Forest best, const Instance& instance,
    std::chrono::steady_clock::time_point program_start) {
    if (env_int("RS_FINAL_ANYTIME_TAIL", 1) == 0) return best;

    publish_best_output(forest_to_output(best, instance.t1));

    const int total_seconds = env_int("RS_TOTAL_SECONDS", 296);
    const int guard_seconds = env_int("RS_FINAL_TAIL_GUARD_SECONDS", 3);
    const double stop_at =
        static_cast<double>(std::max(1, total_seconds - guard_seconds));
    auto elapsed_seconds = [&]() {
        return std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - program_start)
            .count();
    };
    auto remaining_seconds = [&]() {
        return std::max(0.0, stop_at - elapsed_seconds());
    };

    const int min_seconds = env_int("RS_FINAL_TAIL_MIN_SECONDS", 4);
    int round = 0;
    while (remaining_seconds() >= static_cast<double>(min_seconds)) {
        bool improved_this_round = false;

        if (env_int("RS_FINAL_TAIL_SPLIT_REPACK", 1) != 0 &&
            remaining_seconds() >= static_cast<double>(min_seconds) &&
            best.active_count <=
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_MAX_COMPONENTS", 12000)) {
            const double budget = std::min<double>(
                remaining_seconds(),
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_SECONDS", 18));
            Forest repacked = run_h45_targeted_split_repack(
                best, instance.t1, instance.t2,
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_MAX_COMPONENTS", 12000),
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_CANDIDATES", 24),
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_GROUP", 2),
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_RESTART_BUDGET", 192),
                env_int("RS_FINAL_TAIL_SPLIT_REPACK_PASSES", 1),
                env_int("RS_FINAL_TAIL_TARGET",
                        std::max(1, best.active_count - 1)),
                budget);
            repacked =
                exact_output_repair(std::move(repacked), instance.t1, instance.t2);
            ValidationResult check =
                validate_forest(repacked, instance.t1, instance.t2);
            if (check.ok && repacked.active_count < best.active_count) {
                best = std::move(repacked);
                improved_this_round = true;
                publish_best_output(forest_to_output(best, instance.t1));
            }
        }

        if (remaining_seconds() < static_cast<double>(min_seconds)) break;

        const int modes[] = {14, 2, 0, 1, 20};
        const int mode = modes[round % 5];
        const double budget = std::min<double>(
            remaining_seconds(), env_int("RS_FINAL_TAIL_RESTART_SECONDS", 22));
        Forest restarted = run_h45_split_restarts(
            best, instance.t1, instance.t2,
            env_int("RS_FINAL_TAIL_RESTART_MODE", mode),
            env_int("RS_FINAL_TAIL_RESTART_ATTEMPTS", 200000),
            env_int("RS_FINAL_TAIL_SEED_OFFSET", 90000) + round * 4099,
            env_int("RS_FINAL_TAIL_EXTRA_SPLITS", 2 + (round % 3)),
            env_int("RS_FINAL_TAIL_FPT", 1) != 0,
            env_int("RS_FINAL_TAIL_TARGET", std::max(1, best.active_count - 1)),
            budget);
        restarted =
            exact_output_repair(std::move(restarted), instance.t1, instance.t2);
        ValidationResult restart_check =
            validate_forest(restarted, instance.t1, instance.t2);
        if (restart_check.ok && restarted.active_count < best.active_count) {
            best = std::move(restarted);
            improved_this_round = true;
            publish_best_output(forest_to_output(best, instance.t1));
        }

        if (std::getenv("RS_FINAL_TAIL_DEBUG") != nullptr) {
            std::cerr << "rs_final_tail_round round=" << round
                      << " best=" << best.active_count
                      << " improved=" << (improved_this_round ? 1 : 0)
                      << " remaining=" << remaining_seconds() << '\n';
        }

        ++round;
        if (round >= env_int("RS_FINAL_TAIL_MAX_ROUNDS", 1000)) break;
    }

    return best;
}

}  // namespace

int main() {
    try {
        const auto program_start = std::chrono::steady_clock::now();
        install_signal_handler();

        Instance instance = read_instance(std::cin);
        configure_giant_fast_mode(instance.n);
        if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr &&
            env_int("RS_GIANT_FAST_MODE", 1) != 0 &&
            instance.n >= env_int("RS_GIANT_FAST_MIN_LEAVES", 10000)) {
            std::cerr << "rs_giant_fast_mode leaves=" << instance.n << '\n';
        }
        if (std::getenv("RS_DUMP_CLUSTER_CANDIDATES") != nullptr) {
            std::vector<ClusterCandidate> candidates =
                find_cluster_candidates(instance.t1, instance.t2);
            dump_cluster_candidates(candidates, std::cerr);
            if (std::getenv("RS_DUMP_CLUSTER_PLAN") != nullptr) {
                ClusterPlan plan = build_cluster_plan(instance.t1, candidates);
                dump_cluster_plan(plan, std::cerr);
                if (std::getenv("RS_DUMP_CLUSTER_NEWICKS") != nullptr) {
                    dump_cluster_newicks(instance.t1, instance.t2, plan, std::cerr);
                }
                if (std::getenv("RS_DUMP_LOCAL_PROBLEMS") != nullptr) {
                    dump_local_problem(
                        build_local_problem(instance.t1, instance.t2, plan, 0, false),
                        plan.residual_id, std::cerr);
                    for (int id : plan.solve_order) {
                        dump_local_problem(
                            build_local_problem(instance.t1, instance.t2, plan, id, true),
                            id, std::cerr);
                    }
                }
            }
            if (std::getenv("RS_DUMP_CLUSTER_CANDIDATES_EXIT") != nullptr) {
                return 0;
            }
        }
        Forest best = make_singleton_forest(instance.n);
        if (!rebuild_forest_data(best, instance.t1, instance.t2)) {
            std::cerr << "failed to initialize singleton forest\n";
            return 1;
        }
        publish_best_output(forest_to_output(best, instance.t1));

        auto finish = [&]() -> int {
            best = run_final_anytime_tail(std::move(best), instance, program_start);
            best = run_highs_pool_repair(std::move(best), instance, program_start);
            std::cout << forest_to_output(best, instance.t1);
            return 0;
        };

        ValidationResult singleton_check = validate_forest(best, instance.t1, instance.t2);
        if (!singleton_check.ok) {
            std::cerr << "singleton invalid: " << singleton_check.message << "\n";
            std::cout << forest_to_output(best, instance.t1);
            return 1;
        }

        std::vector<ClusterCandidate> route_candidates =
            find_cluster_candidates(instance.t1, instance.t2);
        ClusterPlan route_plan = build_cluster_plan(instance.t1, route_candidates);

        if (should_run_big_101_route(instance, route_plan)) {
            if (env_int("RS_KERNEL_101_ROUTE", 1) != 0 &&
                run_kernel_101_route(instance, best)) {
                return finish();
            }
            const int mode2_max_leaves = env_int("RS_101_MODE2_MAX_LEAVES", 1500);
            if (instance.n <= mode2_max_leaves) {
                configure_101_mode2_defaults();
            } else {
                configure_big_101_defaults();
            }
            if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                std::cerr << "rs_101_route_start leaves=" << instance.n
                          << " mode="
                          << (instance.n <= mode2_max_leaves ? "mode2" : "mixed")
                          << " proper_clusters=" << route_plan.clusters.size()
                          << " residual_leaves=" << route_plan.residual_leaves
                          << '\n';
            }
            Forest fallback = run_no_cluster_fallback(instance);
            ValidationResult fallback_check =
                validate_forest(fallback, instance.t1, instance.t2);
            if (fallback_check.ok && fallback.active_count < best.active_count) {
                best = std::move(fallback);
                publish_best_output(forest_to_output(best, instance.t1));
            }
            if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                std::cerr << "rs_101_route_done components="
                          << best.active_count << '\n';
            }
            return finish();
        }

        bool skip_early_for_direct_adaptive = false;
        if (env_int("RS_DIRECT_ADAPTIVE_SKIP_EARLY", 1) != 0 &&
            env_int("RS_DIRECT_ADAPTIVE_CHAIN", 1) != 0 &&
            instance.n >= env_int("RS_DIRECT_ADAPTIVE_MIN_LEAVES", 800)) {
            skip_early_for_direct_adaptive = route_plan.clusters.empty();
        }

        if (env_int("RS_EARLY_GREEDY_BASELINE", 1) != 0 &&
            !skip_early_for_direct_adaptive) {
            Forest early = make_singleton_forest(instance.n);
            if (!rebuild_forest_data(early, instance.t1, instance.t2)) {
                std::cerr << "early greedy singleton invalid\n";
                return 1;
            }
            early = run_h45_greedy_merge(
                std::move(early), instance.t1, instance.t2,
                env_int("RS_EARLY_GREEDY_MODE", 2), false, false);
            if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                std::cerr << "rs_early_greedy_baseline components="
                          << early.active_count << '\n';
            }
            keep_and_publish_if_better(std::move(early), instance, best);
        }

        if (const char* seed_path = std::getenv("RS_SEED_FOREST_PATH")) {
            Forest seed;
            if (!load_seed_forest(seed_path, instance.n, seed) ||
                !rebuild_forest_data(seed, instance.t1, instance.t2)) {
                std::cerr << "seed forest load failed\n";
                return 1;
            }
            ValidationResult seed_check =
                validate_forest(seed, instance.t1, instance.t2);
            if (!seed_check.ok) {
                std::cerr << "seed forest invalid: " << seed_check.message << "\n";
                return 1;
            }
            publish_best_output(forest_to_output(seed, instance.t1));
            Forest candidate = run_h45_split_restarts(
                seed, instance.t1, instance.t2,
                env_int("RS_FINAL_SCORE_MODE", 14),
                env_int("RS_FINAL_RESTART_BUDGET", 200000),
                env_int("RS_FINAL_SEED_OFFSET", 6700),
                env_int("RS_FINAL_EXTRA_SPLITS", 2),
                true,
                env_int("RS_FINAL_TARGET", std::max(1, seed.active_count - 1)),
                static_cast<double>(env_int("RS_FINAL_SECONDS", 60)));
            candidate = exact_output_repair(std::move(candidate),
                                            instance.t1, instance.t2);
            ValidationResult candidate_check =
                validate_forest(candidate, instance.t1, instance.t2);
            if (candidate_check.ok && candidate.active_count < seed.active_count) {
                seed = std::move(candidate);
            }
            if (env_int("RS_LOCAL_EXCHANGE_REPAIR", 0) != 0) {
                Forest exchange = run_h45_local_exchange_repair(
                    seed, instance.t1, instance.t2,
                    env_int("RS_LOCAL_EXCHANGE_MAX_COMPONENTS", 260),
                    env_int("RS_LOCAL_EXCHANGE_CANDIDATES", 16),
                    env_int("RS_LOCAL_EXCHANGE_WINDOW", 24),
                    env_int("RS_LOCAL_EXCHANGE_MAX_BLOCKERS", 8),
                    env_int("RS_LOCAL_EXCHANGE_RESTART_BUDGET", 256),
                    env_int("RS_LOCAL_EXCHANGE_TARGET",
                            std::max(1, seed.active_count - 1)),
                    static_cast<double>(env_int("RS_LOCAL_EXCHANGE_SECONDS", 12)));
                exchange = exact_output_repair(std::move(exchange),
                                               instance.t1, instance.t2);
                ValidationResult exchange_check =
                    validate_forest(exchange, instance.t1, instance.t2);
                if (exchange_check.ok && exchange.active_count < seed.active_count) {
                    seed = std::move(exchange);
                }
            }
            std::cerr << "seed_restart_final components="
                      << seed.active_count << '\n';
            best = std::move(seed);
            publish_best_output(forest_to_output(best, instance.t1));
            return finish();
        }

        if (std::getenv("RS_DISABLE_CLUSTER_SOLVER") == nullptr) {
            const ClusterPlan& plan = route_plan;
            int route_102_largest_cluster = 0;
            if (should_run_102_route(instance, plan, route_102_largest_cluster)) {
                if (env_int("RS_KERNEL_102_ROUTE", 1) != 0 &&
                    run_kernel_102_route(instance, best)) {
                    return finish();
                }
                best = run_102_route(instance, plan, std::move(best));
                return finish();
            }
            if (run_kernel_cluster_route(instance, best)) {
                return finish();
            }
            if (run_kernel_first_default_route(instance, best)) {
                return finish();
            }
            bool greedy12_ran = false;
            bool greedy12_improved = false;
            if (should_run_greedy12_tail(instance.n, plan.residual_leaves)) {
                greedy12_ran = true;
                if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                    std::cerr << "rs_greedy12_start leaves=" << instance.n
                              << " residual_leaves=" << plan.residual_leaves
                              << '\n';
                }
                Forest greedy12;
                if (run_greedy12_tail(instance.t1, instance.t2, greedy12)) {
                    greedy12 = exact_output_repair(std::move(greedy12),
                                                   instance.t1, instance.t2);
                    ValidationResult greedy12_check =
                        validate_forest(greedy12, instance.t1, instance.t2);
                    if (std::getenv("RS_DEBUG_CLUSTER_SOLVER") != nullptr) {
                        std::cerr << "rs_greedy12_done valid="
                                  << (greedy12_check.ok ? 1 : 0)
                                  << " components=" << greedy12.active_count
                                  << '\n';
                    }
                    if (greedy12_check.ok &&
                        greedy12.active_count < best.active_count) {
                        best = std::move(greedy12);
                        greedy12_improved = true;
                        publish_best_output(forest_to_output(best, instance.t1));
                        if (env_int("RS_GREEDY12_RETURN", 1) != 0) {
                            return finish();
                        }
                    }
                }
            }
            if (greedy12_ran && greedy12_improved &&
                env_int("RS_GREEDY12_RETURN", 1) != 0) {
                return finish();
            }
            if (plan.clusters.empty()) {
                Forest fallback = run_no_cluster_fallback(instance);
                ValidationResult fallback_check =
                    validate_forest(fallback, instance.t1, instance.t2);
                if (fallback_check.ok && fallback.active_count < best.active_count) {
                    best = std::move(fallback);
                    publish_best_output(forest_to_output(best, instance.t1));
                }
                return finish();
            }

            Forest clustered;
            if (run_cluster_solver(instance.t1, instance.t2, clustered)) {
                ValidationResult cluster_check =
                    validate_forest(clustered, instance.t1, instance.t2);
                if (cluster_check.ok && clustered.active_count < best.active_count) {
                    best = std::move(clustered);
                    publish_best_output(forest_to_output(best, instance.t1));
                    return finish();
                }
            }
        }

        GreedyConfig config;
        Forest candidate = run_greedy_merge(best, instance.t1, instance.t2, config);
        ValidationResult check = validate_forest(candidate, instance.t1, instance.t2);
        if (check.ok && candidate.active_count < best.active_count) {
            best = std::move(candidate);
            publish_best_output(forest_to_output(best, instance.t1));
        }

        return finish();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
