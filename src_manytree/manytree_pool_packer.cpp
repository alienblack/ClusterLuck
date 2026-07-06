#include <algorithm>
#include <array>
#include <chrono>
#include <csignal>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>

#ifdef USE_HIGHS
#include "Highs.h"
#endif
#ifdef USE_CADICAL
#include "cadical.hpp"
#endif

using namespace std;

static string g_signal_output;

static void sigterm_emit_best(int) {
    _Exit(124);
}

#ifndef TLE_TWO_TREE_DEFAULT
#define TLE_TWO_TREE_DEFAULT 0
#endif
#ifndef AUTO_TARGET_LOOP_DEFAULT
#define AUTO_TARGET_LOOP_DEFAULT 0
#endif
#ifndef STRICT_EXACT_DEFAULT
#define STRICT_EXACT_DEFAULT 0
#endif
#ifndef ALLOW_STABLE_EXCHANGE_DEFAULT
#define ALLOW_STABLE_EXCHANGE_DEFAULT 0
#endif
#ifndef ALLOW_EMPIRICAL_EXACT_DEFAULT
#define ALLOW_EMPIRICAL_EXACT_DEFAULT 0
#endif
#ifndef DEFAULT_TOTAL_SECONDS
#define DEFAULT_TOTAL_SECONDS 0.0
#endif
#ifndef DEFAULT_OUTPUT_SAFETY_SECONDS
#define DEFAULT_OUTPUT_SAFETY_SECONDS 8.0
#endif
#ifndef DEFAULT_GLOBAL_SEARCH_SECONDS
#define DEFAULT_GLOBAL_SEARCH_SECONDS 0.0
#endif
#ifndef DEFAULT_GLOBAL_BRANCH_SECONDS
#define DEFAULT_GLOBAL_BRANCH_SECONDS 0.0
#endif
#ifndef DEFAULT_HIGH_DIVERGENCE_PACK_SECONDS
#define DEFAULT_HIGH_DIVERGENCE_PACK_SECONDS 0.0
#endif
#ifndef DEFAULT_EXACT_EXCHANGE_SECONDS
#define DEFAULT_EXACT_EXCHANGE_SECONDS 0.0
#endif
#ifndef DEFAULT_PAIR_AUGMENT_SECONDS
#define DEFAULT_PAIR_AUGMENT_SECONDS 0.0
#endif
#ifndef DEFAULT_PAIR_AUGMENT_ATTEMPT_SECONDS
#define DEFAULT_PAIR_AUGMENT_ATTEMPT_SECONDS 1.5
#endif
#ifndef DEFAULT_INCUMBENT_REPACK_SECONDS
#define DEFAULT_INCUMBENT_REPACK_SECONDS 0.0
#endif
#ifndef DEFAULT_INCUMBENT_REPACK_ATTEMPT_SECONDS
#define DEFAULT_INCUMBENT_REPACK_ATTEMPT_SECONDS 2.0
#endif
#ifndef DEFAULT_SMALL_SAVINGS_EXACT_SECONDS
#define DEFAULT_SMALL_SAVINGS_EXACT_SECONDS 180.0
#endif
#ifndef DEFAULT_SMALL_SAVINGS_PAIR_NODE_LIMIT
#define DEFAULT_SMALL_SAVINGS_PAIR_NODE_LIMIT 80000000LL
#endif
#ifndef DEFAULT_SMALL_SAVINGS_SPECIAL_NODE_LIMIT
#define DEFAULT_SMALL_SAVINGS_SPECIAL_NODE_LIMIT 5000000LL
#endif
#ifndef DEFAULT_SMALL_SAVINGS_MAX_TARGET
#define DEFAULT_SMALL_SAVINGS_MAX_TARGET 13
#endif
#ifndef DEFAULT_SMALL_SAVINGS_COMPONENT_LIMIT
#define DEFAULT_SMALL_SAVINGS_COMPONENT_LIMIT 100000
#endif
#ifndef DEFAULT_REFERENCE_CUT_SECONDS
#define DEFAULT_REFERENCE_CUT_SECONDS 0.0
#endif
#ifndef DEFAULT_EXCHANGE_SECONDS
#define DEFAULT_EXCHANGE_SECONDS 60.0
#endif
#ifndef DEFAULT_EXCHANGE_MAX_COMPONENT_SIZE
#define DEFAULT_EXCHANGE_MAX_COMPONENT_SIZE 8
#endif
#ifndef DEFAULT_EXCHANGE_MAX_REGION_COMPONENTS
#define DEFAULT_EXCHANGE_MAX_REGION_COMPONENTS 10
#endif
#ifndef DEFAULT_EXCHANGE_MAX_CANDIDATES
#define DEFAULT_EXCHANGE_MAX_CANDIDATES 500000
#endif
#ifndef DEFAULT_EXCHANGE_MAX_REGIONS
#define DEFAULT_EXCHANGE_MAX_REGIONS 200000
#endif
#ifndef DEFAULT_EXCHANGE_NODE_LIMIT
#define DEFAULT_EXCHANGE_NODE_LIMIT 2000000
#endif
#ifndef DEFAULT_POOL_SECONDS
#define DEFAULT_POOL_SECONDS 260.0
#endif
#ifndef DEFAULT_SWAP_SECONDS
#define DEFAULT_SWAP_SECONDS 240.0
#endif
#ifndef DEFAULT_MEDIUM_SWAP_SECONDS
#define DEFAULT_MEDIUM_SWAP_SECONDS 180.0
#endif

struct TreeNode { int label = -1, parent = -1; vector<int> child; };
struct Tree {
    vector<TreeNode> nodes;
    int root = -1, leaf_count = 0, LOG = 1;
    vector<int> leaf_node, depth;
    vector<vector<int>> up;
    int add_node(int label) { nodes.push_back(TreeNode{label, -1, {}}); return (int)nodes.size() - 1; }
};
struct Instance { int tree_count = 0, leaf_count = 0; vector<Tree> trees; };

struct Options {
    string instance_path, output_components, output_newick;
    string check_components;
    string dump_packing_lp, dump_packing_blocks;
    string dump_high_divergence_lp, dump_high_divergence_blocks;
    int target_components = -1;
    bool has_target = false;
    bool tle_two_trees = TLE_TWO_TREE_DEFAULT;
    bool auto_target_loop = AUTO_TARGET_LOOP_DEFAULT;
    bool strict_exact = STRICT_EXACT_DEFAULT;
    bool exact_branch_decision = false;
    bool exact_pair_lb_decision = false;
    bool direct_lazy_decision = false;
    bool reference_cut_decision = false;
    bool savings_lb_decision = false;
    bool edge_owner_decision = false;
    bool high_divergence_pool_decision = false;
    bool seeded_pair_decomposition_decision = false;
    bool exchange_certify = false;
    bool allow_stable_exchange = ALLOW_STABLE_EXCHANGE_DEFAULT;
    bool allow_empirical_exact = ALLOW_EMPIRICAL_EXACT_DEFAULT;
    bool full_sat_before_lazy = false;
    string exchange_incumbent;
    int seed = 3403;
    int max_pool = 900000;
    int subtree_samples = 500;
    int beam_width = 45000;
    int pack_pool_limit = 900000;
    int always_keep_size = 5;
    int anchor_limit = 3000;
    int swap_anchor_limit = 2;
    int swap_candidate_limit = 50000;
    int swap_max_old_components = 7;
    int swap_branch_width = 500;
    int swap_small_keep_size = 3;
    int swap_anchor_slack = 10;
    int swap_extra_leaf_limit = 14;
    int medium_swap_candidate_limit = 120000;
    int medium_swap_max_old_components = 8;
    int medium_swap_max_replacement_size = 14;
    int local_guided_components = 62;
    int local_seed_triples_limit = 1000;
    int local_variant_seed_limit = 12000;
    int local_trials = 7000;
    int local_random_extra_components = 8;
    int local_random_anchor_components = 1;
    int local_guided_growth_trials = 0;
    int local_guided_choice_window = 6;
    int local_guided_future_probe = 48;
    int local_max_pool = 45000;
    int local_max_component_size = 40;
    int local_top_per_leaf = 24;
    int local_column_max_size = 8;
    int local_column_max_edges = 80;
    int global_greedy_starts = 100;
    int global_remove_count = 6;
    int global_remove_jitter = 8;
    int global_candidate_limit = 250000;
    int global_highs_top_per_leaf = 160;
    int global_highs_column_max_size = 8;
    int global_highs_column_max_edges = 120;
    int global_branch_top_per_leaf = 80;
    int global_branch_column_max_size = 8;
    int global_branch_column_max_edges = 120;
    int global_branch_width = 160;
    int high_divergence_max_component_size = 8;
    int high_divergence_column_max_edges = 180;
    int high_divergence_max_components = 200000;
    int high_divergence_sat_candidate_limit = 12000;
    int high_divergence_top_per_leaf = 220;
    int seeded_pair_decision_seed_limit = 0;
    int exact_exchange_max_region_components = 14;
    int exact_exchange_max_region_leaves = 42;
    int exact_exchange_max_candidates = 4500;
    int exact_exchange_top_per_leaf = 160;
    int exact_exchange_column_max_size = 8;
    int exact_exchange_column_max_edges = 180;
    int exact_exchange_seed_limit = 12000;
    int pair_augment_candidate_limit = 6000;
    int pair_augment_bundle_size = 5;
    int pair_augment_max_region_components = 20;
    int pair_augment_column_max_size = 8;
    int pair_augment_column_max_edges = 120;
    int pair_augment_top_per_leaf = 80;
    int incumbent_repack_candidate_limit = 8000;
    int incumbent_repack_max_region_components = 28;
    int incumbent_repack_column_max_size = 8;
    int incumbent_repack_column_max_edges = 140;
    int incumbent_repack_top_per_leaf = 120;
    int incumbent_repack_branch_width = 220;
    int direct_pair_leaf_limit = 64;
    int direct_lazy_preload_side = 0;
    double pool_seconds = DEFAULT_POOL_SECONDS;
    double swap_seconds = DEFAULT_SWAP_SECONDS;
    double medium_swap_seconds = DEFAULT_MEDIUM_SWAP_SECONDS;
    double local_seconds = 0.0;
    double local_round_seconds = 8.0;
    double global_search_seconds = DEFAULT_GLOBAL_SEARCH_SECONDS;
    double global_highs_seconds = 0.0;
    double global_branch_seconds = DEFAULT_GLOBAL_BRANCH_SECONDS;
    double high_divergence_pack_seconds = DEFAULT_HIGH_DIVERGENCE_PACK_SECONDS;
    double high_divergence_decision_seconds = 0.0;
    double seeded_pair_decision_seconds = 0.0;
    double exact_exchange_seconds = DEFAULT_EXACT_EXCHANGE_SECONDS;
    double exact_exchange_attempt_seconds = 2.5;
    double pair_augment_seconds = DEFAULT_PAIR_AUGMENT_SECONDS;
    double pair_augment_attempt_seconds = DEFAULT_PAIR_AUGMENT_ATTEMPT_SECONDS;
    double incumbent_repack_seconds = DEFAULT_INCUMBENT_REPACK_SECONDS;
    double incumbent_repack_attempt_seconds = DEFAULT_INCUMBENT_REPACK_ATTEMPT_SECONDS;
    double small_savings_exact_seconds = DEFAULT_SMALL_SAVINGS_EXACT_SECONDS;
    double direct_pair_seconds = 120.0;
    double exact_branch_seconds = 0.0;
    double exact_pair_lb_seconds = 0.0;
    double reference_cut_seconds = DEFAULT_REFERENCE_CUT_SECONDS;
    double savings_lb_seconds = 0.0;
    long long exact_branch_node_limit = 2000000;
    int exact_pair_lb_max_pairs = 0;
    int savings_lb_max_depth = 64;
    int savings_lb_max_components_per_triple = 250000;
    int exchange_max_component_size = DEFAULT_EXCHANGE_MAX_COMPONENT_SIZE;
    int exchange_max_component_edges = 160;
    int exchange_max_region_components = DEFAULT_EXCHANGE_MAX_REGION_COMPONENTS;
    int exchange_max_candidates = DEFAULT_EXCHANGE_MAX_CANDIDATES;
    int exchange_max_regions = DEFAULT_EXCHANGE_MAX_REGIONS;
    long long savings_lb_pair_node_limit = 200000;
    long long small_savings_pair_node_limit = DEFAULT_SMALL_SAVINGS_PAIR_NODE_LIMIT;
    long long small_savings_special_node_limit = DEFAULT_SMALL_SAVINGS_SPECIAL_NODE_LIMIT;
    int small_savings_component_limit = DEFAULT_SMALL_SAVINGS_COMPONENT_LIMIT;
    int small_savings_max_target = DEFAULT_SMALL_SAVINGS_MAX_TARGET;
    long long exchange_node_limit = DEFAULT_EXCHANGE_NODE_LIMIT;
    double exchange_seconds = DEFAULT_EXCHANGE_SECONDS;
    double total_seconds = DEFAULT_TOTAL_SECONDS;
    double output_safety_seconds = DEFAULT_OUTPUT_SAFETY_SECONDS;
};

static void skip_ws(const string &s, int &p) { while (p < (int)s.size() && isspace((unsigned char)s[p])) p++; }
static void skip_label_or_length(const string &s, int &p) {
    skip_ws(s, p);
    while (p < (int)s.size() && s[p] != ',' && s[p] != ')' && s[p] != ';') p++;
}
static int parse_subtree(Tree &tree, const string &s, int &p) {
    skip_ws(s, p);
    if (p >= (int)s.size()) throw runtime_error("unexpected end of Newick");
    if (s[p] == '(') {
        p++;
        int node = tree.add_node(-1);
        while (true) {
            int child = parse_subtree(tree, s, p);
            tree.nodes[child].parent = node;
            tree.nodes[node].child.push_back(child);
            skip_ws(s, p);
            if (p >= (int)s.size()) throw runtime_error("unterminated internal node");
            if (s[p] == ',') { p++; continue; }
            if (s[p] == ')') { p++; skip_label_or_length(s, p); return node; }
            throw runtime_error("unexpected character inside internal node");
        }
    }
    if (!isdigit((unsigned char)s[p])) throw runtime_error("expected numeric leaf");
    int label = 0;
    while (p < (int)s.size() && isdigit((unsigned char)s[p])) label = label * 10 + (s[p++] - '0');
    int node = tree.add_node(label);
    skip_label_or_length(s, p);
    return node;
}
static void index_tree(Tree &tree, int node, int parent, int d) {
    tree.nodes[node].parent = parent;
    tree.depth[node] = d;
    if (tree.nodes[node].label >= 1) tree.leaf_node[tree.nodes[node].label] = node;
    for (int child : tree.nodes[node].child) index_tree(tree, child, node, d + 1);
}
static void build_lca(Tree &tree) {
    tree.LOG = 1;
    while ((1 << tree.LOG) <= max(1, (int)tree.nodes.size())) tree.LOG++;
    tree.up.assign(tree.LOG, vector<int>(tree.nodes.size(), -1));
    for (int i = 0; i < (int)tree.nodes.size(); i++) tree.up[0][i] = tree.nodes[i].parent;
    for (int k = 1; k < tree.LOG; k++) {
        for (int i = 0; i < (int)tree.nodes.size(); i++) {
            int p = tree.up[k - 1][i];
            tree.up[k][i] = p < 0 ? -1 : tree.up[k - 1][p];
        }
    }
}
static int lca_node(const Tree &tree, int a, int b) {
    if (tree.depth[a] < tree.depth[b]) swap(a, b);
    int diff = tree.depth[a] - tree.depth[b];
    for (int k = 0; k < tree.LOG; k++) if (diff & (1 << k)) a = tree.up[k][a];
    if (a == b) return a;
    for (int k = tree.LOG - 1; k >= 0; k--) if (tree.up[k][a] != tree.up[k][b]) {
        a = tree.up[k][a]; b = tree.up[k][b];
    }
    return tree.nodes[a].parent;
}
static int lca_leaf_set(const Tree &tree, const vector<int> &leaves) {
    int root = tree.leaf_node[leaves[0]];
    for (int i = 1; i < (int)leaves.size(); i++) root = lca_node(tree, root, tree.leaf_node[leaves[i]]);
    return root;
}
static Tree parse_tree(const string &text, int n) {
    Tree tree; tree.leaf_count = n; tree.leaf_node.assign(n + 1, -1);
    int p = 0; tree.root = parse_subtree(tree, text, p);
    tree.depth.assign(tree.nodes.size(), 0);
    index_tree(tree, tree.root, -1, 0);
    build_lca(tree);
    return tree;
}
static Instance read_instance_stream(istream &in) {
    Instance inst; vector<string> texts; string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.size() >= 2 && line[1] == 'p') {
                string tag; stringstream ss(line); ss >> tag >> inst.tree_count >> inst.leaf_count;
            }
            continue;
        }
        texts.push_back(line);
    }
    if ((int)texts.size() != inst.tree_count) throw runtime_error("tree count mismatch");
    for (const string &text : texts) inst.trees.push_back(parse_tree(text, inst.leaf_count));
    return inst;
}
static Instance read_instance_path(const string &path) {
    if (path.empty() || path == "-") return read_instance_stream(cin);
    ifstream in(path);
    if (!in) throw runtime_error("could not open instance");
    return read_instance_stream(in);
}

static string combine_sig(string a, string b) { if (b < a) swap(a, b); return "(" + a + "," + b + ")"; }
static string restriction_sig(const Tree &tree, int node, const vector<unsigned char> &mark) {
    const TreeNode &tn = tree.nodes[node];
    if (tn.label >= 1) return mark[tn.label] ? to_string(tn.label) : "";
    vector<string> kids;
    for (int child : tn.child) {
        string s = restriction_sig(tree, child, mark);
        if (!s.empty()) kids.push_back(std::move(s));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return combine_sig(kids[0], kids[1]);
}
static bool valid_topology(const Instance &inst, const vector<int> &leaves) {
    if (leaves.empty()) return false;
    vector<unsigned char> mark(inst.leaf_count + 1, 0);
    for (int x : leaves) mark[x] = 1;
    string ref;
    for (int t = 0; t < inst.tree_count; t++) {
        int root = lca_leaf_set(inst.trees[t], leaves);
        string sig = restriction_sig(inst.trees[t], root, mark);
        if (sig.empty()) return false;
        if (t == 0) ref = std::move(sig);
        else if (sig != ref) return false;
    }
    return true;
}
static vector<int> tree_component_edges(const Tree &tree, const vector<int> &leaves) {
    vector<int> edges;
    if (leaves.size() <= 1) return edges;
    int root = lca_leaf_set(tree, leaves);
    for (int label : leaves) {
        int node = tree.leaf_node[label];
        while (node != root) { edges.push_back(node); node = tree.nodes[node].parent; }
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    return edges;
}
static vector<int> component_edges_all(const Instance &inst, const vector<int> &leaves, const vector<int> &edge_offset) {
    vector<int> out;
    for (int t = 0; t < inst.tree_count; t++) {
        vector<int> local = tree_component_edges(inst.trees[t], leaves);
        for (int e : local) out.push_back(edge_offset[t] + e);
    }
    sort(out.begin(), out.end());
    return out;
}
static int triple_sig(const Tree &tree, int a, int b, int c) {
    int ab = lca_node(tree, tree.leaf_node[a], tree.leaf_node[b]);
    int ac = lca_node(tree, tree.leaf_node[a], tree.leaf_node[c]);
    int bc = lca_node(tree, tree.leaf_node[b], tree.leaf_node[c]);
    if (tree.depth[ab] >= tree.depth[ac] && tree.depth[ab] >= tree.depth[bc]) return 0;
    if (tree.depth[ac] >= tree.depth[ab] && tree.depth[ac] >= tree.depth[bc]) return 1;
    return 2;
}
static uint32_t triple_key(int a, int b, int c, int n) {
    return (uint32_t)((a * (n + 1) + b) * (n + 1) + c);
}
static unordered_set<uint32_t> collect_valid_triples(const Instance &inst) {
    unordered_set<uint32_t> valid;
    int n = inst.leaf_count;
    for (int a = 1; a <= n; a++) for (int b = a + 1; b <= n; b++) for (int c = b + 1; c <= n; c++) {
        int sig = triple_sig(inst.trees[0], a, b, c);
        bool ok = true;
        for (int t = 1; t < inst.tree_count; t++) if (triple_sig(inst.trees[t], a, b, c) != sig) { ok = false; break; }
        if (ok) valid.insert(triple_key(a, b, c, n));
    }
    return valid;
}
static bool triple_compatible(const vector<int> &block, const unordered_set<uint32_t> &valid, int n) {
    if (block.size() <= 2) return true;
    for (int i = 0; i < (int)block.size(); i++) for (int j = i + 1; j < (int)block.size(); j++) for (int k = j + 1; k < (int)block.size(); k++) {
        if (!valid.count(triple_key(block[i], block[j], block[k], n))) return false;
    }
    return true;
}
static string key_of(const vector<int> &v) {
    string s;
    for (int x : v) { s += to_string(x); s += ','; }
    return s;
}

struct Component { vector<int> leaves, edges; double score = 0.0; };
static vector<vector<int>> local_branch_cover(const vector<Component> &columns,
        const vector<int> &selected_leaves, int max_components, double seconds,
        int branch_width);
static vector<Component> reduce_repair_columns_preserve_pairs(const vector<Component> &columns,
        const unordered_map<int, int> &leaf_to_local, int local_leaf_count, int per_leaf);
struct Pool {
    vector<Component> comps;
    unordered_map<string, int> id;
};
static bool add_component(Pool &pool, const Instance &inst, const unordered_set<uint32_t> &valid,
        const vector<int> &edge_offset, vector<int> block) {
    sort(block.begin(), block.end());
    block.erase(unique(block.begin(), block.end()), block.end());
    if (block.empty()) return false;
    string key = key_of(block);
    if (pool.id.count(key)) return false;
    if (!triple_compatible(block, valid, inst.leaf_count)) return false;
    if (!valid_topology(inst, block)) return false;
    vector<int> edges = component_edges_all(inst, block, edge_offset);
    double score = (double)(block.size() - 1) / pow(max(1, (int)edges.size()), 0.25);
    int idx = (int)pool.comps.size();
    pool.id[key] = idx;
    pool.comps.push_back(Component{std::move(block), std::move(edges), score});
    return true;
}
static void collect_leaf_sets(const Tree &tree, int node, vector<vector<int>> &sets) {
    if (tree.nodes[node].label >= 1) { sets.push_back({tree.nodes[node].label}); return; }
    vector<int> leaves;
    for (int child : tree.nodes[node].child) {
        collect_leaf_sets(tree, child, sets);
        vector<int> child_leaves = sets.back();
        leaves.insert(leaves.end(), child_leaves.begin(), child_leaves.end());
    }
    sort(leaves.begin(), leaves.end());
    if (leaves.size() >= 2) sets.push_back(leaves);
}
static vector<int> random_sample(const vector<int> &labels, int k, mt19937_64 &rng) {
    vector<int> tmp = labels;
    shuffle(tmp.begin(), tmp.end(), rng);
    tmp.resize(min(k, (int)tmp.size()));
    sort(tmp.begin(), tmp.end());
    return tmp;
}

using Clock = chrono::steady_clock;
#ifdef USE_CADICAL
struct DeadlineTerminator : public CaDiCaL::Terminator {
    Clock::time_point deadline;
    explicit DeadlineTerminator(Clock::time_point deadline_) : deadline(deadline_) {}
    bool terminate() override { return Clock::now() >= deadline; }
};
#endif
static double elapsed(Clock::time_point start) { return chrono::duration<double>(Clock::now() - start).count(); }

[[noreturn]] static void wait_for_timeout(const string &reason) {
    cerr << "STRICT_TLE reason=" << reason << "\n";
    while (true) this_thread::sleep_for(chrono::hours(24));
}

static Pool build_pool(const Instance &inst, const Options &opt, const vector<int> &edge_offset) {
    auto start = Clock::now();
    mt19937_64 rng(opt.seed);
    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    Pool pool;
    int n = inst.leaf_count;
    int pool_limit = opt.max_pool;
    double pool_seconds = opt.pool_seconds;
    if (opt.strict_exact && !opt.has_target &&
            inst.tree_count >= 20 && n >= 180 &&
            (int)valid.size() <= max(32, n / 2)) {
        pool_limit = min(pool_limit, 60000);
        pool_seconds = min(pool_seconds, 4.0);
        cerr << "LOW_TRIPLE_POOL_PROFILE valid_triples=" << valid.size()
             << " max_pool=" << pool_limit
             << " pool_seconds=" << pool_seconds << "\n";
    }
    vector<int> labels(n);
    for (int i = 0; i < n; i++) labels[i] = i + 1;
    for (int a = 1; a <= n; a++) for (int b = a + 1; b <= n; b++) add_component(pool, inst, valid, edge_offset, vector<int>{a, b});
    for (int a = 1; a <= n; a++) for (int b = a + 1; b <= n; b++) for (int c = b + 1; c <= n; c++) {
        if (valid.count(triple_key(a, b, c, n))) add_component(pool, inst, valid, edge_offset, vector<int>{a, b, c});
    }
    for (const Tree &tree : inst.trees) {
        vector<vector<int>> sets;
        collect_leaf_sets(tree, tree.root, sets);
        for (const vector<int> &s : sets) {
            if ((int)s.size() <= n) add_component(pool, inst, valid, edge_offset, s);
            if (s.size() > 3) {
                for (int r = 0; r < opt.subtree_samples && (int)pool.comps.size() < pool_limit; r++) {
                    uniform_int_distribution<int> dist(3, (int)s.size());
                    add_component(pool, inst, valid, edge_offset, random_sample(s, dist(rng), rng));
                }
            }
        }
    }
    uniform_real_distribution<double> prob(0.0, 1.0);
    while ((int)pool.comps.size() < pool_limit && elapsed(start) < pool_seconds) {
        int start_size = prob(rng) < 0.75 ? 2 : 3;
        vector<int> current = random_sample(labels, start_size, rng);
        if (!triple_compatible(current, valid, n)) continue;
        vector<int> order = labels;
        shuffle(order.begin(), order.end(), rng);
        for (int leaf : order) {
            if ((int)current.size() >= n) break;
            if (binary_search(current.begin(), current.end(), leaf)) continue;
            vector<int> cand = current;
            cand.push_back(leaf);
            sort(cand.begin(), cand.end());
            if (triple_compatible(cand, valid, n)) {
                current = cand;
                if (prob(rng) < 0.92) add_component(pool, inst, valid, edge_offset, current);
            }
        }
        add_component(pool, inst, valid, edge_offset, current);
    }
    cerr << "POOL components=" << pool.comps.size() << " seconds=" << elapsed(start) << "\n";
    if ((int)pool.comps.size() > opt.pack_pool_limit) {
        vector<int> order(pool.comps.size());
        for (int i = 0; i < (int)order.size(); i++) order[i] = i;
        sort(order.begin(), order.end(), [&](int a, int b) {
            bool keepa = (int)pool.comps[a].leaves.size() <= opt.always_keep_size;
            bool keepb = (int)pool.comps[b].leaves.size() <= opt.always_keep_size;
            if (keepa != keepb) return keepa > keepb;
            if (pool.comps[a].score != pool.comps[b].score) return pool.comps[a].score > pool.comps[b].score;
            return pool.comps[a].leaves.size() > pool.comps[b].leaves.size();
        });
        Pool pruned;
        for (int i = 0; i < opt.pack_pool_limit; i++) {
            Component c = pool.comps[order[i]];
            pruned.id[key_of(c.leaves)] = (int)pruned.comps.size();
            pruned.comps.push_back(std::move(c));
        }
        pool = std::move(pruned);
        cerr << "PRUNE after=" << pool.comps.size() << "\n";
    }
    return pool;
}

struct Column { vector<int> leaves, edges; vector<unsigned long long> mask, leaf_mask; double score; };
enum class ExactDecision { No, Yes, Unknown };
static void set_bit(vector<unsigned long long> &bits, int pos) { bits[pos >> 6] |= 1ULL << (pos & 63); }
static bool intersects(const vector<unsigned long long> &a, const vector<unsigned long long> &b) {
    for (int i = 0; i < (int)a.size(); i++) if (a[i] & b[i]) return true;
    return false;
}
static void bit_or_into(vector<unsigned long long> &a, const vector<unsigned long long> &b) {
    for (int i = 0; i < (int)a.size(); i++) a[i] |= b[i];
}
static int bit_count(const vector<unsigned long long> &a) {
    int total = 0; for (auto x : a) total += __builtin_popcountll(x); return total;
}
static vector<Column> encode_columns(const Instance &inst, const Pool &pool, int total_bits) {
    int words = (total_bits + 63) >> 6;
    vector<Column> cols;
    for (const Component &c : pool.comps) {
        Column col; col.leaves = c.leaves; col.edges = c.edges; col.score = c.score;
        col.mask.assign(words, 0); col.leaf_mask.assign(words, 0);
        for (int leaf : c.leaves) { set_bit(col.mask, leaf - 1); set_bit(col.leaf_mask, leaf - 1); }
        for (int edge : c.edges) set_bit(col.mask, inst.leaf_count + edge);
        cols.push_back(std::move(col));
    }
    sort(cols.begin(), cols.end(), [](const Column &a, const Column &b) {
        if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
        return a.score > b.score;
    });
    return cols;
}
static vector<int> greedy_with_anchor(const vector<Column> &cols, int anchor, int n) {
    (void)n;
    vector<int> selected{anchor};
    vector<unsigned long long> used = cols[anchor].mask;
    vector<int> order(cols.size());
    for (int i = 0; i < (int)order.size(); i++) order[i] = i;
    sort(order.begin(), order.end(), [&](int a, int b) {
        if (cols[a].score != cols[b].score) return cols[a].score > cols[b].score;
        return cols[a].leaves.size() > cols[b].leaves.size();
    });
    for (int idx : order) {
        if (idx == anchor) continue;
        if (!intersects(used, cols[idx].mask)) {
            selected.push_back(idx);
            bit_or_into(used, cols[idx].mask);
        }
    }
    return selected;
}
static int component_count(const vector<Column> &cols, const vector<int> &selected, int n) {
    vector<unsigned long long> leaf(((n + 63) >> 6), 0);
    for (int idx : selected) {
        for (int x : cols[idx].leaves) leaf[(x - 1) >> 6] |= 1ULL << ((x - 1) & 63);
    }
    return (int)selected.size() + n - bit_count(leaf);
}

static int column_weight(const Column &col) {
    return max(0, (int)col.leaves.size() - 1);
}

static bool mask_has_bit(const vector<unsigned long long> &bits, int pos) {
    return (bits[pos >> 6] >> (pos & 63)) & 1ULL;
}

static void dump_packing_lp_files(const vector<Column> &cols, const Options &opt) {
    if (opt.dump_packing_lp.empty()) return;
    ofstream lp(opt.dump_packing_lp);
    if (!lp) throw runtime_error("could not open --dump-packing-lp");
    lp << "Maximize\n obj:";
    for (int i = 0; i < (int)cols.size(); i++) {
        int weight = column_weight(cols[i]);
        if (weight > 0) lp << " + " << weight << " x" << i;
        if (i % 1000 == 999) lp << "\n";
    }
    lp << "\nSubject To\n";
    int row = 0;
    if (!cols.empty()) {
        vector<vector<int>> by_bit(cols[0].mask.size() * 64);
        for (int i = 0; i < (int)cols.size(); i++) {
            for (int w = 0; w < (int)cols[i].mask.size(); w++) {
                unsigned long long bits = cols[i].mask[w];
                while (bits) {
                    int bit = __builtin_ctzll(bits);
                    int global_bit = (w << 6) + bit;
                    if (global_bit < (int)by_bit.size()) by_bit[global_bit].push_back(i);
                    bits &= bits - 1;
                }
            }
        }
        for (const vector<int> &bucket : by_bit) {
            if (bucket.empty()) continue;
            lp << " c" << row++ << ":";
            for (int j = 0; j < (int)bucket.size(); j++) {
                lp << " + x" << bucket[j];
                if (j % 1000 == 999) lp << "\n";
            }
            lp << " <= 1\n";
        }
    }
    lp << "Binary\n";
    for (int i = 0; i < (int)cols.size(); i++) lp << " x" << i << "\n";
    lp << "End\n";

    if (!opt.dump_packing_blocks.empty()) {
        ofstream blocks(opt.dump_packing_blocks);
        if (!blocks) throw runtime_error("could not open --dump-packing-blocks");
        for (int i = 0; i < (int)cols.size(); i++) {
            blocks << i;
            for (int leaf : cols[i].leaves) blocks << "," << leaf;
            blocks << "\n";
        }
    }
    cerr << "DUMP_PACKING_LP vars=" << cols.size()
         << " rows=" << row
         << " lp=" << opt.dump_packing_lp
         << " blocks=" << opt.dump_packing_blocks << "\n";
}

static int mask_bit_count(const vector<unsigned long long> &mask) {
    int total = 0;
    for (unsigned long long word : mask) total += __builtin_popcountll(word);
    return total;
}

static void rebuild_state(const vector<Column> &cols, const vector<int> &input,
        vector<int> &selected, vector<unsigned long long> &used, int &total) {
    selected.clear();
    used.assign(cols.empty() ? 0 : cols[0].mask.size(), 0);
    total = 0;
    for (int idx : input) {
        if (idx < 0 || idx >= (int)cols.size()) continue;
        if (intersects(used, cols[idx].mask)) continue;
        selected.push_back(idx);
        bit_or_into(used, cols[idx].mask);
        total += column_weight(cols[idx]);
    }
}

static vector<int> greedy_lns_solution(const vector<Column> &cols, const vector<int> &candidates,
        int seed, vector<unsigned long long> &used, int &total) {
    mt19937_64 rng(seed);
    uniform_real_distribution<double> noise(0.0, 1.0);
    vector<pair<double, int>> order;
    order.reserve(candidates.size());
    for (int idx : candidates) {
        int footprint = max(1, mask_bit_count(cols[idx].mask));
        double score = column_weight(cols[idx]) / pow((double)footprint, 0.25) + noise(rng);
        order.push_back({score, idx});
    }
    sort(order.begin(), order.end(), greater<pair<double, int>>());
    vector<int> selected;
    used.assign(cols.empty() ? 0 : cols[0].mask.size(), 0);
    total = 0;
    for (auto &item : order) {
        int idx = item.second;
        if (intersects(used, cols[idx].mask)) continue;
        selected.push_back(idx);
        bit_or_into(used, cols[idx].mask);
        total += column_weight(cols[idx]);
    }
    return selected;
}

static bool improve_one_exchange_lns(const vector<Column> &cols, const vector<int> &candidates,
        vector<int> &selected, vector<unsigned long long> &used, int &total) {
    if (cols.empty()) return false;
    vector<int> owner(cols[0].mask.size() * 64, -1);
    vector<unsigned char> selected_flag(cols.size(), 0);
    for (int pos = 0; pos < (int)selected.size(); pos++) {
        int idx = selected[pos];
        selected_flag[idx] = 1;
        for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
            unsigned long long bits = cols[idx].mask[w];
            while (bits) {
                unsigned long long low = bits & -bits;
                owner[(w << 6) + __builtin_ctzll(bits)] = pos;
                bits ^= low;
            }
        }
    }
    vector<int> touched;
    vector<int> mark(selected.size(), 0);
    int stamp = 0;
    for (int cand : candidates) {
        if (selected_flag[cand]) continue;
        stamp++;
        touched.clear();
        bool any = false;
        for (int w = 0; w < (int)cols[cand].mask.size(); w++) {
            unsigned long long bits = cols[cand].mask[w] & used[w];
            if (bits) any = true;
            while (bits) {
                unsigned long long low = bits & -bits;
                int bit = (w << 6) + __builtin_ctzll(bits);
                int pos = owner[bit];
                if (pos >= 0 && mark[pos] != stamp) {
                    mark[pos] = stamp;
                    touched.push_back(pos);
                }
                bits ^= low;
            }
        }
        if (!any) continue;
        int lost = 0;
        for (int pos : touched) lost += column_weight(cols[selected[pos]]);
        if (column_weight(cols[cand]) <= lost) continue;
        vector<int> next;
        next.reserve(selected.size() - touched.size() + 1);
        vector<unsigned char> remove(selected.size(), 0);
        for (int pos : touched) remove[pos] = 1;
        for (int pos = 0; pos < (int)selected.size(); pos++) if (!remove[pos]) next.push_back(selected[pos]);
        next.push_back(cand);
        rebuild_state(cols, next, selected, used, total);
        return true;
    }
    return false;
}

static vector<int> lns_refill_solution(const vector<Column> &cols, const vector<int> &candidates,
        const vector<int> &base, int seed, int remove_count, vector<unsigned long long> &used, int &total) {
    mt19937_64 rng(seed);
    vector<int> selected = base;
    shuffle(selected.begin(), selected.end(), rng);
    if (remove_count < (int)selected.size()) selected.resize((int)selected.size() - remove_count);
    else selected.clear();
    vector<int> kept = selected;
    rebuild_state(cols, kept, selected, used, total);
    vector<unsigned char> selected_flag(cols.size(), 0);
    for (int idx : selected) selected_flag[idx] = 1;
    uniform_real_distribution<double> noise(0.0, 1.0);
    double alpha_values[] = {0.12, 0.18, 0.25, 0.35, 0.50};
    double noise_values[] = {0.05, 0.2, 0.6, 1.2};
    double alpha = alpha_values[seed % 5];
    double noise_scale = noise_values[seed % 4];
    vector<pair<double, int>> order;
    for (int idx : candidates) {
        if (selected_flag[idx]) continue;
        int footprint = max(1, mask_bit_count(cols[idx].mask));
        double score = column_weight(cols[idx]) / pow((double)footprint, alpha) + noise(rng) * noise_scale;
        order.push_back({score, idx});
    }
    sort(order.begin(), order.end(), greater<pair<double, int>>());
    for (auto &item : order) {
        int idx = item.second;
        if (intersects(used, cols[idx].mask)) continue;
        selected.push_back(idx);
        bit_or_into(used, cols[idx].mask);
        total += column_weight(cols[idx]);
    }
    return selected;
}

static vector<int> global_lns_pack(const vector<Column> &cols, int n, const Options &opt) {
    if (opt.global_search_seconds <= 0.0 || cols.empty()) return {};
    vector<int> candidates(cols.size());
    for (int i = 0; i < (int)cols.size(); i++) candidates[i] = i;
    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        double sa = wa / pow((double)max(1, mask_bit_count(cols[a].mask)), 0.25);
        double sb = wb / pow((double)max(1, mask_bit_count(cols[b].mask)), 0.25);
        if (sa != sb) return sa > sb;
        if (cols[a].leaves.size() != cols[b].leaves.size()) return cols[a].leaves.size() > cols[b].leaves.size();
        return cols[a].leaves < cols[b].leaves;
    });
    if (opt.global_candidate_limit > 0 && (int)candidates.size() > opt.global_candidate_limit) {
        vector<int> keep;
        vector<unsigned char> used(cols.size(), 0);
        for (int idx : candidates) {
            if ((int)cols[idx].leaves.size() <= opt.always_keep_size) {
                keep.push_back(idx);
                used[idx] = 1;
            }
        }
        for (int idx : candidates) {
            if ((int)keep.size() >= opt.global_candidate_limit) break;
            if (used[idx]) continue;
            keep.push_back(idx);
        }
        candidates.swap(keep);
    }
    auto deadline = Clock::now() + chrono::milliseconds((int)(opt.global_search_seconds * 1000));
    vector<int> best;
    int best_total = -1;
    int best_count = n;
    int seed = 1;
    while (Clock::now() < deadline) {
        vector<unsigned long long> used;
        int total = 0;
        vector<int> selected;
        if (seed <= opt.global_greedy_starts || best.empty()) {
            selected = greedy_lns_solution(cols, candidates, opt.seed * 100000 + seed, used, total);
        } else {
            int remove_count = opt.global_remove_count + (seed % max(1, opt.global_remove_jitter));
            selected = lns_refill_solution(cols, candidates, best, opt.seed * 100000 + seed, remove_count, used, total);
        }
        bool changed = true;
        while (changed && Clock::now() < deadline) {
            changed = improve_one_exchange_lns(cols, candidates, selected, used, total);
        }
        int count = component_count(cols, selected, n);
        if (count < best_count || total > best_total) {
            best = selected;
            best_total = total;
            best_count = count;
            cerr << "GLOBAL_LNS_BEST components=" << best_count
                 << " savings=" << best_total
                 << " selected=" << best.size()
                 << " seed=" << seed << "\n";
            if (best_count <= opt.target_components) break;
        }
        seed++;
    }
    return best;
}

static vector<int> anchored_pack(const vector<Column> &cols, int n, const Options &opt) {
    vector<int> anchors;
    for (int i = 0; i < (int)cols.size(); i++) if ((int)cols[i].leaves.size() >= 30) anchors.push_back(i);
    if ((int)anchors.size() > opt.anchor_limit) anchors.resize(opt.anchor_limit);
    vector<int> best;
    int best_count = n;
    for (int pos = 0; pos < (int)anchors.size(); pos++) {
        vector<int> candidate = greedy_with_anchor(cols, anchors[pos], n);
        int count = component_count(cols, candidate, n);
        if (count < best_count) {
            best_count = count; best = candidate;
            cerr << "ANCHOR_BEST components=" << best_count << " anchor=" << (pos + 1)
                 << " size=" << cols[anchors[pos]].leaves.size() << "\n";
            if (best_count <= opt.target_components) break;
        }
    }
    return best;
}

static vector<vector<int>> selected_partition(const vector<Column> &cols, const vector<int> &selected, int n) {
    vector<unsigned char> used(n + 1, 0);
    vector<vector<int>> out;
    for (int idx : selected) { out.push_back(cols[idx].leaves); for (int x : cols[idx].leaves) used[x] = 1; }
    for (int x = 1; x <= n; x++) if (!used[x]) out.push_back(vector<int>{x});
    sort(out.begin(), out.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return out;
}

static bool validate_partition(const Instance &inst, const vector<vector<int>> &partition,
        const vector<int> &edge_offset) {
    vector<unsigned char> seen(inst.leaf_count + 1, 0);
    unordered_set<int> used_edges;
    for (const vector<int> &block : partition) {
        if (!valid_topology(inst, block)) return false;
        for (int x : block) {
            if (seen[x]) return false;
            seen[x] = 1;
        }
        vector<int> edges = component_edges_all(inst, block, edge_offset);
        for (int e : edges) {
            if (!used_edges.insert(e).second) return false;
        }
    }
    for (int x = 1; x <= inst.leaf_count; x++) if (!seen[x]) return false;
    return true;
}

static vector<vector<int>> singleton_merge_repair(const Instance &inst, vector<vector<int>> partition,
        int target, const vector<int> &edge_offset) {
    bool changed = true;
    while ((int)partition.size() > target && changed) {
        changed = false;
        for (int si = 0; si < (int)partition.size() && !changed; si++) {
            if (partition[si].size() != 1) continue;
            for (int bi = 0; bi < (int)partition.size() && !changed; bi++) {
                if (bi == si || partition[bi].empty()) continue;
                vector<vector<int>> candidate;
                candidate.reserve(partition.size() - 1);
                vector<int> merged = partition[bi];
                merged.push_back(partition[si][0]);
                sort(merged.begin(), merged.end());
                for (int i = 0; i < (int)partition.size(); i++) {
                    if (i == si) continue;
                    if (i == bi) candidate.push_back(merged);
                    else candidate.push_back(partition[i]);
                }
                sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
                    if (a[0] != b[0]) return a[0] < b[0];
                    if (a.size() != b.size()) return a.size() < b.size();
                    return a < b;
                });
                if (validate_partition(inst, candidate, edge_offset)) {
                    cerr << "MERGE_SINGLETON leaf=" << partition[si][0]
                         << " into_size=" << partition[bi].size()
                         << " components=" << candidate.size() << "\n";
                    partition = std::move(candidate);
                    changed = true;
                }
            }
        }
    }
    return partition;
}

static vector<vector<int>> pair_packing_repair(const Instance &inst,
        const vector<vector<int>> &partition, int target_components,
        const vector<int> &edge_offset, double seconds = 8.0) {
    int n = inst.leaf_count;
    int target_savings = n - target_components;
    int current_savings = n - (int)partition.size();
    if (target_savings <= current_savings || target_savings <= 0) return {};
    if (target_savings > 12 || n > 240 || seconds <= 0.0) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    int edge_words = (edge_offset.back() + 63) >> 6;
    struct PairChoice {
        int a = 0, b = 0;
        vector<unsigned long long> edges;
        int edge_count = 0;
    };
    vector<PairChoice> pairs;
    vector<vector<int>> by_leaf(n + 1);
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            vector<int> leaves{a, b};
            if (!valid_topology(inst, leaves)) continue;
            vector<int> edges = component_edges_all(inst, leaves, edge_offset);
            PairChoice pc;
            pc.a = a;
            pc.b = b;
            pc.edge_count = (int)edges.size();
            pc.edges.assign(edge_words, 0);
            for (int edge : edges) set_bit(pc.edges, edge);
            int idx = (int)pairs.size();
            pairs.push_back(std::move(pc));
            by_leaf[a].push_back(idx);
            by_leaf[b].push_back(idx);
        }
    }
    for (int leaf = 1; leaf <= n; leaf++) {
        sort(by_leaf[leaf].begin(), by_leaf[leaf].end(), [&](int x, int y) {
            if (pairs[x].edge_count != pairs[y].edge_count) {
                return pairs[x].edge_count < pairs[y].edge_count;
            }
            if (pairs[x].a != pairs[y].a) return pairs[x].a < pairs[y].a;
            return pairs[x].b < pairs[y].b;
        });
    }
    vector<unsigned char> used_leaf(n + 1, 0);
    vector<unsigned long long> used_edges(edge_words, 0);
    vector<int> chosen, answer;
    long long nodes = 0;
    const long long node_limit = 2000000;
    function<bool()> dfs = [&]() {
        if (Clock::now() >= deadline || nodes++ >= node_limit) return false;
        if ((int)chosen.size() >= target_savings) {
            answer = chosen;
            return true;
        }
        int remaining = 0;
        for (int leaf = 1; leaf <= n; leaf++) if (!used_leaf[leaf]) remaining++;
        if ((int)chosen.size() + remaining / 2 < target_savings) return false;
        int pivot = -1;
        int best_options = INT_MAX;
        for (int leaf = 1; leaf <= n; leaf++) {
            if (used_leaf[leaf]) continue;
            int feasible = 0;
            for (int idx : by_leaf[leaf]) {
                int other = pairs[idx].a == leaf ? pairs[idx].b : pairs[idx].a;
                if (used_leaf[other]) continue;
                if (intersects(used_edges, pairs[idx].edges)) continue;
                feasible++;
            }
            if (feasible < best_options) {
                best_options = feasible;
                pivot = leaf;
                if (feasible == 0) break;
            }
        }
        if (pivot < 0) return false;
        used_leaf[pivot] = 1;
        if (dfs()) return true;
        used_leaf[pivot] = 0;
        for (int idx : by_leaf[pivot]) {
            int other = pairs[idx].a == pivot ? pairs[idx].b : pairs[idx].a;
            if (used_leaf[pivot] || used_leaf[other]) continue;
            if (intersects(used_edges, pairs[idx].edges)) continue;
            vector<unsigned long long> old_edges = used_edges;
            used_leaf[pivot] = used_leaf[other] = 1;
            bit_or_into(used_edges, pairs[idx].edges);
            chosen.push_back(idx);
            if (dfs()) return true;
            chosen.pop_back();
            used_edges = std::move(old_edges);
            used_leaf[pivot] = used_leaf[other] = 0;
        }
        return false;
    };
    bool hit = dfs();
    cerr << "PAIR_PACK_REPAIR target_savings=" << target_savings
         << " current_savings=" << current_savings
         << " hit=" << (hit ? 1 : 0)
         << " pairs=" << pairs.size()
         << " nodes=" << nodes
         << " seconds=" << elapsed(start) << "\n";
    if (!hit) return {};
    vector<unsigned char> covered(n + 1, 0);
    vector<vector<int>> out;
    for (int idx : answer) {
        out.push_back(vector<int>{pairs[idx].a, pairs[idx].b});
        covered[pairs[idx].a] = covered[pairs[idx].b] = 1;
    }
    for (int leaf = 1; leaf <= n; leaf++) {
        if (!covered[leaf]) out.push_back(vector<int>{leaf});
    }
    sort(out.begin(), out.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    if (!validate_partition(inst, out, edge_offset)) return {};
    if ((int)out.size() >= (int)partition.size()) return {};
    return out;
}

static vector<vector<int>> small_savings_pack_repair(const Instance &inst,
        const vector<Column> &cols, int target_components,
        const vector<int> &edge_offset) {
    int target_savings = inst.leaf_count - target_components;
    if (target_savings <= 0 || target_savings > 12 || cols.empty()) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds(30000);
    int words = (int)cols[0].mask.size();
    vector<int> pairs, special;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 3 && sz <= min(8, target_savings + 1)) special.push_back(i);
    }
    auto better = [&](int a, int b) {
        int sa = (int)cols[a].leaves.size() - 1;
        int sb = (int)cols[b].leaves.size() - 1;
        if (sa != sb) return sa > sb;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    sort(pairs.begin(), pairs.end(), better);
    sort(special.begin(), special.end(), better);
    if ((int)special.size() > 5000) special.resize(5000);

    vector<int> chosen_special, chosen_pairs, answer_special, answer_pairs;
    vector<unsigned long long> used(words, 0);
    long long special_nodes = 0, pair_nodes = 0;

    auto rebuild_used = [&]() {
        vector<unsigned long long> rebuilt(words, 0);
        for (int idx : chosen_special) bit_or_into(rebuilt, cols[idx].mask);
        for (int idx : chosen_pairs) bit_or_into(rebuilt, cols[idx].mask);
        return rebuilt;
    };

    function<bool(int, int, vector<unsigned long long>&)> fill_pairs =
        [&](int start_idx, int need, vector<unsigned long long> &pair_used) {
            if (need <= 0) return true;
            if (Clock::now() >= deadline || pair_nodes++ > 3000000) return false;
            if ((int)pairs.size() - start_idx < need) return false;
            for (int pos = start_idx; pos < (int)pairs.size(); pos++) {
                int idx = pairs[pos];
                if (intersects(pair_used, cols[idx].mask)) continue;
                bit_or_into(pair_used, cols[idx].mask);
                chosen_pairs.push_back(idx);
                if (fill_pairs(pos + 1, need - 1, pair_used)) return true;
                chosen_pairs.pop_back();
                pair_used = rebuild_used();
            }
            return false;
        };

    function<bool(int, int, vector<unsigned long long>&)> dfs_special =
        [&](int start_idx, int savings, vector<unsigned long long> &cur_used) {
            if (Clock::now() >= deadline || special_nodes++ > 300000) return false;
            int remaining = target_savings - savings;
            if (remaining <= 0) {
                answer_special = chosen_special;
                answer_pairs = chosen_pairs;
                return true;
            }
            chosen_pairs.clear();
            vector<unsigned long long> pair_used = cur_used;
            if (fill_pairs(0, remaining, pair_used)) {
                answer_special = chosen_special;
                answer_pairs = chosen_pairs;
                return true;
            }
            chosen_pairs.clear();
            for (int pos = start_idx; pos < (int)special.size(); pos++) {
                int idx = special[pos];
                int add = (int)cols[idx].leaves.size() - 1;
                if (add <= 0 || add > remaining) continue;
                if (intersects(cur_used, cols[idx].mask)) continue;
                vector<unsigned long long> next_used = cur_used;
                bit_or_into(next_used, cols[idx].mask);
                chosen_special.push_back(idx);
                if (dfs_special(pos + 1, savings + add, next_used)) return true;
                chosen_special.pop_back();
            }
            return false;
        };

    bool hit = dfs_special(0, 0, used);
    cerr << "SMALL_SAVINGS_PACK target_savings=" << target_savings
         << " hit=" << (hit ? 1 : 0)
         << " special=" << special.size()
         << " pairs=" << pairs.size()
         << " special_nodes=" << special_nodes
         << " pair_nodes=" << pair_nodes
         << " seconds=" << elapsed(start) << "\n";
    if (!hit) return {};

    vector<unsigned char> covered(inst.leaf_count + 1, 0);
    vector<vector<int>> out;
    for (int idx : answer_special) {
        out.push_back(cols[idx].leaves);
        for (int leaf : cols[idx].leaves) covered[leaf] = 1;
    }
    for (int idx : answer_pairs) {
        out.push_back(cols[idx].leaves);
        for (int leaf : cols[idx].leaves) covered[leaf] = 1;
    }
    for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
        if (!covered[leaf]) out.push_back(vector<int>{leaf});
    }
    sort(out.begin(), out.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    if ((int)out.size() <= target_components && validate_partition(inst, out, edge_offset)) return out;
    return {};
}

struct SmallSavingsExactResult {
    vector<vector<int>> partition;
    bool proved_infeasible = false;
};

static SmallSavingsExactResult small_savings_exhaustive_pack(
        const Instance &inst, int target_components,
        const vector<int> &edge_offset, double seconds = 45.0,
        long long pair_node_limit = DEFAULT_SMALL_SAVINGS_PAIR_NODE_LIMIT,
        long long special_node_limit = DEFAULT_SMALL_SAVINGS_SPECIAL_NODE_LIMIT,
        int max_target_savings = DEFAULT_SMALL_SAVINGS_MAX_TARGET,
        int component_limit = DEFAULT_SMALL_SAVINGS_COMPONENT_LIMIT) {
    SmallSavingsExactResult result;
    int target_savings = inst.leaf_count - target_components;
    if (target_savings <= 0 || target_savings > max_target_savings) return result;
    auto start = Clock::now();
    cerr << "SMALL_SAVINGS_EXACT_START target_savings=" << target_savings
         << " seconds=" << seconds
         << " pair_node_limit=" << pair_node_limit
         << " special_node_limit=" << special_node_limit
         << " component_limit=" << component_limit << "\n";
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    int total_bits = inst.leaf_count + edge_offset.back();
    int words = (total_bits + 63) >> 6;
    struct MiniComponent {
        vector<int> leaves;
        vector<unsigned long long> mask;
    };
    vector<MiniComponent> pairs, special;
    unordered_set<string> seen;
    bool exhausted = true;

    auto add_component_exact = [&](const vector<int> &leaves, vector<MiniComponent> &where) {
        string key = key_of(leaves);
        if (!seen.insert(key).second) return;
        if (!valid_topology(inst, leaves)) return;
        MiniComponent c;
        c.leaves = leaves;
        c.mask.assign(words, 0);
        for (int leaf : leaves) set_bit(c.mask, leaf - 1);
        vector<int> edges = component_edges_all(inst, leaves, edge_offset);
        for (int edge : edges) set_bit(c.mask, inst.leaf_count + edge);
        where.push_back(std::move(c));
    };

    for (int a = 1; a <= inst.leaf_count; a++) {
        for (int b = a + 1; b <= inst.leaf_count; b++) {
            vector<int> leaves{a, b};
            add_component_exact(leaves, pairs);
        }
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    cerr << "SMALL_SAVINGS_EXACT_TRIPLES target_savings=" << target_savings
         << " valid=" << valid.size()
         << " seconds=" << elapsed(start) << "\n";
    vector<int> current;
    function<void(int)> grow = [&](int next_leaf) {
        if (!exhausted) return;
        if (Clock::now() >= deadline) {
            exhausted = false;
            return;
        }
        if ((int)current.size() >= 3) {
            add_component_exact(current, special);
            if (target_savings >= 10 && component_limit > 0 &&
                    (int)special.size() > component_limit) {
                exhausted = false;
                return;
            }
        }
        if ((int)current.size() >= target_savings + 1) return;
        for (int leaf = next_leaf; leaf <= inst.leaf_count; leaf++) {
            current.push_back(leaf);
            if (triple_compatible(current, valid, inst.leaf_count)) grow(leaf + 1);
            current.pop_back();
            if (!exhausted) return;
        }
    };
    for (int a = 1; a <= inst.leaf_count && exhausted; a++) {
        for (int b = a + 1; b <= inst.leaf_count && exhausted; b++) {
            for (int c = b + 1; c <= inst.leaf_count; c++) {
                if (!valid.count(triple_key(a, b, c, inst.leaf_count))) continue;
                current = vector<int>{a, b, c};
                grow(c + 1);
                if (!exhausted) break;
            }
        }
    }

    auto better_pair = [](const MiniComponent &a, const MiniComponent &b) {
        if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
        return a.leaves < b.leaves;
    };
    sort(pairs.begin(), pairs.end(), better_pair);
    sort(special.begin(), special.end(), better_pair);

    if (special.empty()) {
        struct PairOnlyResult {
            bool hit = false;
            bool proved_infeasible = false;
            vector<int> chosen;
            long long nodes = 0;
        };
        PairOnlyResult pair_only;
        vector<vector<int>> by_leaf(inst.leaf_count + 1);
        for (int i = 0; i < (int)pairs.size(); i++) {
            if (pairs[i].leaves.size() != 2) continue;
            by_leaf[pairs[i].leaves[0]].push_back(i);
            by_leaf[pairs[i].leaves[1]].push_back(i);
        }
        for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
            sort(by_leaf[leaf].begin(), by_leaf[leaf].end(), [&](int a, int b) {
                if (pairs[a].leaves != pairs[b].leaves) return pairs[a].leaves < pairs[b].leaves;
                return a < b;
            });
        }
        vector<unsigned char> used_leaf(inst.leaf_count + 1, 0);
        vector<unsigned long long> used_mask(words, 0);
        vector<int> chosen;
        bool pair_search_exhausted = true;
        function<bool()> dfs_pairs_exact = [&]() {
            if (Clock::now() >= deadline || pair_only.nodes++ >= pair_node_limit) {
                pair_search_exhausted = false;
                return false;
            }
            if ((int)chosen.size() >= target_savings) {
                pair_only.hit = true;
                pair_only.chosen = chosen;
                return true;
            }
            int remaining_leaves = 0;
            for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
                if (!used_leaf[leaf]) remaining_leaves++;
            }
            if ((int)chosen.size() + remaining_leaves / 2 < target_savings) return false;

            int pivot = -1;
            int best_options = INT_MAX;
            for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
                if (used_leaf[leaf]) continue;
                int feasible = 0;
                for (int idx : by_leaf[leaf]) {
                    int a = pairs[idx].leaves[0], b = pairs[idx].leaves[1];
                    int other = (a == leaf) ? b : a;
                    if (used_leaf[other]) continue;
                    if (intersects(used_mask, pairs[idx].mask)) continue;
                    feasible++;
                }
                if (feasible < best_options) {
                    best_options = feasible;
                    pivot = leaf;
                    if (feasible == 0) break;
                }
            }
            if (pivot < 0) return false;

            used_leaf[pivot] = 1;
            if (dfs_pairs_exact()) return true;
            used_leaf[pivot] = 0;
            if (!pair_search_exhausted) return false;

            for (int idx : by_leaf[pivot]) {
                int a = pairs[idx].leaves[0], b = pairs[idx].leaves[1];
                int other = (a == pivot) ? b : a;
                if (used_leaf[pivot] || used_leaf[other]) continue;
                if (intersects(used_mask, pairs[idx].mask)) continue;
                vector<unsigned long long> old_mask = used_mask;
                used_leaf[pivot] = used_leaf[other] = 1;
                bit_or_into(used_mask, pairs[idx].mask);
                chosen.push_back(idx);
                if (dfs_pairs_exact()) return true;
                chosen.pop_back();
                used_mask = std::move(old_mask);
                used_leaf[pivot] = used_leaf[other] = 0;
                if (!pair_search_exhausted) return false;
            }
            return false;
        };
        dfs_pairs_exact();
        pair_only.proved_infeasible = !pair_only.hit && pair_search_exhausted;
        cerr << "SMALL_SAVINGS_PAIR_ONLY_EXACT target_savings=" << target_savings
             << " hit=" << (pair_only.hit ? 1 : 0)
             << " proved_infeasible=" << (pair_only.proved_infeasible ? 1 : 0)
             << " pairs=" << pairs.size()
             << " nodes=" << pair_only.nodes
             << " exhausted=" << (pair_search_exhausted ? 1 : 0)
             << " seconds=" << elapsed(start) << "\n";
        if (pair_only.hit) {
            vector<unsigned char> covered(inst.leaf_count + 1, 0);
            for (int idx : pair_only.chosen) {
                result.partition.push_back(pairs[idx].leaves);
                for (int leaf : pairs[idx].leaves) covered[leaf] = 1;
            }
            for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
                if (!covered[leaf]) result.partition.push_back(vector<int>{leaf});
            }
            sort(result.partition.begin(), result.partition.end(), [](const vector<int> &a, const vector<int> &b) {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            if ((int)result.partition.size() > target_components ||
                    !validate_partition(inst, result.partition, edge_offset)) {
                result.partition.clear();
            }
            return result;
        }
        if (pair_only.proved_infeasible) {
            result.proved_infeasible = true;
            return result;
        }
    }

#ifdef USE_CADICAL
    int component_sat_limit = target_savings >= 10 ? 25000 : 10000;
    int component_sat_ms = target_savings >= 10 ? 30000 : 12000;
    int component_sat_conflicts = target_savings >= 10 ? 2500000 : 900000;
    int component_sat_decisions = target_savings >= 10 ? 12000000 : 5000000;
    if ((int)pairs.size() + (int)special.size() <= component_sat_limit) {
        vector<const MiniComponent*> all_components;
        for (const MiniComponent &c : special) all_components.push_back(&c);
        for (const MiniComponent &c : pairs) all_components.push_back(&c);
        CaDiCaL::Solver solver;
        solver.set("quiet", 1);
        solver.set("factor", 0);
        solver.limit("conflicts", component_sat_conflicts);
        solver.limit("decisions", component_sat_decisions);
        Clock::time_point sat_deadline = min(deadline, Clock::now() + chrono::milliseconds(component_sat_ms));
        DeadlineTerminator terminator(sat_deadline);
        solver.connect_terminator(&terminator);
        int var_count = (int)all_components.size();
        long long sat_clauses = 0;
        auto sat_add_clause = [&](const vector<int> &lits) {
            for (int lit : lits) solver.add(lit);
            solver.add(0);
            sat_clauses++;
        };
        auto sat_add_at_most_one = [&](const vector<int> &vars) {
            int m = (int)vars.size();
            if (m <= 1) return;
            vector<int> seq(m, 0);
            for (int i = 1; i < m; i++) seq[i] = ++var_count;
            sat_add_clause({-vars[0], seq[1]});
            for (int i = 2; i < m; i++) {
                sat_add_clause({-seq[i - 1], seq[i]});
                sat_add_clause({-vars[i - 1], -seq[i - 1]});
                sat_add_clause({-vars[i - 1], seq[i]});
            }
            sat_add_clause({-vars[m - 1], -seq[m - 1]});
        };
        auto sat_add_weighted_at_least_k = [&](const vector<int> &vars, const vector<int> &weights, int need) {
            int m = (int)vars.size();
            if (need <= 0) return;
            int total_weight = 0;
            for (int w : weights) total_weight += w;
            if (total_weight < need) {
                sat_add_clause({});
                return;
            }
            vector<vector<int>> dp(m, vector<int>(need + 1, 0));
            for (int i = 0; i < m; i++) {
                int prefix = 0;
                for (int p = 0; p <= i; p++) prefix += weights[p];
                for (int s = 1; s <= min(need, prefix); s++) dp[i][s] = ++var_count;
            }
            for (int i = 0; i < m; i++) {
                int w = min(weights[i], need);
                for (int s = 1; s <= need; s++) {
                    if (!dp[i][s]) continue;
                    int same = (i > 0 ? dp[i - 1][s] : 0);
                    int prev_need = max(0, s - w);
                    int take_prev = (i > 0 && prev_need > 0 ? dp[i - 1][prev_need] : 0);
                    if (same) sat_add_clause({-same, dp[i][s]});
                    if (prev_need == 0) {
                        sat_add_clause({-vars[i], dp[i][s]});
                    } else if (take_prev) {
                        sat_add_clause({-vars[i], -take_prev, dp[i][s]});
                    }
                    vector<int> rev_choice{-dp[i][s]};
                    if (same) rev_choice.push_back(same);
                    rev_choice.push_back(vars[i]);
                    sat_add_clause(rev_choice);
                    if (prev_need > 0) {
                        vector<int> rev_prev{-dp[i][s]};
                        if (same) rev_prev.push_back(same);
                        if (take_prev) rev_prev.push_back(take_prev);
                        sat_add_clause(rev_prev);
                    }
                }
            }
            sat_add_clause({dp[m - 1][need]});
        };
        int total_bits = inst.leaf_count + edge_offset.back();
        vector<vector<int>> by_bit(total_bits);
        vector<int> vars;
        vector<int> weights;
        for (int i = 0; i < (int)all_components.size(); i++) {
            int var = i + 1;
            const MiniComponent &c = *all_components[i];
            vars.push_back(var);
            weights.push_back((int)c.leaves.size() - 1);
            for (int w = 0; w < (int)c.mask.size(); w++) {
                unsigned long long x = c.mask[w];
                while (x) {
                    int bit = __builtin_ctzll(x);
                    int global_bit = (w << 6) + bit;
                    if (global_bit < total_bits) by_bit[global_bit].push_back(var);
                    x &= x - 1;
                }
            }
        }
        long long estimated_amo_vars = 0;
        long long estimated_amo_clauses = 0;
        for (const vector<int> &bucket : by_bit) {
            int m = (int)bucket.size();
            if (m <= 1) continue;
            estimated_amo_vars += m - 1;
            estimated_amo_clauses += 3LL * m - 4;
        }
        long long estimated_weight_vars = 1LL * (int)all_components.size() * target_savings;
        long long estimated_weight_clauses = 4LL * (int)all_components.size() * target_savings + 1;
        long long estimated_vars = (int)all_components.size() +
            estimated_amo_vars + estimated_weight_vars;
        long long estimated_clauses = estimated_amo_clauses + estimated_weight_clauses;
        if (estimated_vars > 8000000LL || estimated_clauses > 25000000LL) {
            cerr << "SMALL_SAVINGS_COMPONENT_SAT_SKIPPED reason=model_too_large"
                 << " target_savings=" << target_savings
                 << " components=" << all_components.size()
                 << " estimated_vars=" << estimated_vars
                 << " estimated_clauses=" << estimated_clauses << "\n";
        } else {
        for (const vector<int> &bucket : by_bit) sat_add_at_most_one(bucket);
        sat_add_weighted_at_least_k(vars, weights, target_savings);
        cerr << "SMALL_SAVINGS_COMPONENT_SAT_SOLVE target_savings=" << target_savings
             << " components=" << all_components.size()
             << " vars=" << var_count
             << " clauses=" << sat_clauses
             << " build_seconds=" << elapsed(start) << "\n";
        int sat_result = solver.solve();
        cerr << "SMALL_SAVINGS_COMPONENT_SAT target_savings=" << target_savings
             << " components=" << all_components.size()
             << " result=" << sat_result
             << " clauses=" << sat_clauses
             << " seconds=" << elapsed(start) << "\n";
        if (sat_result == 20) {
            result.proved_infeasible = true;
            return result;
        }
        if (sat_result == 10) {
            vector<unsigned char> covered(inst.leaf_count + 1, 0);
            int savings = 0;
            for (int i = 0; i < (int)all_components.size(); i++) {
                if (solver.val(i + 1) <= 0) continue;
                const vector<int> &leaves = all_components[i]->leaves;
                result.partition.push_back(leaves);
                savings += (int)leaves.size() - 1;
                for (int leaf : leaves) covered[leaf] = 1;
            }
            for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
                if (!covered[leaf]) result.partition.push_back(vector<int>{leaf});
            }
            sort(result.partition.begin(), result.partition.end(), [](const vector<int> &a, const vector<int> &b) {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            if (savings >= target_savings &&
                    (int)result.partition.size() <= target_components &&
                    validate_partition(inst, result.partition, edge_offset)) {
                return result;
            }
            result.partition.clear();
        }
        }
    }
#endif

    vector<int> chosen_special, chosen_pairs, answer_special, answer_pairs;
    long long pair_nodes = 0, special_nodes = 0;
    bool search_exhausted = true;
    struct PairFillCacheEntry {
        unsigned char state = 0; // 1 = hit, 2 = proved infeasible
        vector<int> chosen;
    };
    unordered_map<string, PairFillCacheEntry> pair_sat_cache;
    unordered_map<string, PairFillCacheEntry> pair_clique_cache;
    long long pair_sat_cache_hits = 0, pair_sat_cache_stores = 0;
    long long pair_clique_cache_hits = 0, pair_clique_cache_stores = 0;
    int pair_bit_words = ((int)pairs.size() + 63) >> 6;
    vector<unsigned char> special_pair_block_ready(special.size(), 0);
    vector<vector<unsigned long long>> special_pair_block(special.size());

    auto pair_cache_key_bits = [&](int need, const vector<unsigned long long> &feasible_bits) {
        string key;
        key.reserve(sizeof(int) + feasible_bits.size() * sizeof(unsigned long long));
        key.append(reinterpret_cast<const char*>(&need), sizeof(int));
        for (unsigned long long word : feasible_bits) {
            key.append(reinterpret_cast<const char*>(&word), sizeof(unsigned long long));
        }
        return key;
    };

    auto bitset_count = [&](const vector<unsigned long long> &bits) {
        int count = 0;
        for (unsigned long long word : bits) count += __builtin_popcountll(word);
        return count;
    };

    auto bitset_to_indices = [&](const vector<unsigned long long> &bits) {
        vector<int> indices;
        indices.reserve(bitset_count(bits));
        for (int w = 0; w < (int)bits.size(); w++) {
            unsigned long long x = bits[w];
            while (x) {
                int bit = __builtin_ctzll(x);
                int idx = (w << 6) + bit;
                if (idx < (int)pairs.size()) indices.push_back(idx);
                x &= x - 1;
            }
        }
        return indices;
    };

    auto blocked_pairs_for_special = [&](int special_idx) -> const vector<unsigned long long>& {
        if (!special_pair_block_ready[special_idx]) {
            vector<unsigned long long> blocked(pair_bit_words, 0);
            for (int i = 0; i < (int)pairs.size(); i++) {
                if (intersects(special[special_idx].mask, pairs[i].mask)) {
                    blocked[i >> 6] |= 1ULL << (i & 63);
                }
            }
            special_pair_block[special_idx] = std::move(blocked);
            special_pair_block_ready[special_idx] = 1;
        }
        return special_pair_block[special_idx];
    };

    auto rebuild_used = [&]() {
        vector<unsigned long long> rebuilt(words, 0);
        for (int idx : chosen_special) bit_or_into(rebuilt, special[idx].mask);
        for (int idx : chosen_pairs) bit_or_into(rebuilt, pairs[idx].mask);
        return rebuilt;
    };

    auto fill_pairs_sat = [&](int need, const vector<unsigned long long> &feasible_bits) {
        struct SatPackResult {
            bool hit = false;
            bool proved_infeasible = false;
            vector<int> chosen;
        };
        SatPackResult res;
        if (need <= 0) {
            res.hit = true;
            return res;
        }
#ifdef USE_CADICAL
        string cache_key = pair_cache_key_bits(need, feasible_bits);
        auto cache_it = pair_sat_cache.find(cache_key);
        if (cache_it != pair_sat_cache.end()) {
            pair_sat_cache_hits++;
            if (cache_it->second.state == 1) {
                res.hit = true;
                res.chosen = cache_it->second.chosen;
            } else if (cache_it->second.state == 2) {
                res.proved_infeasible = true;
            }
            return res;
        }
        int feasible_count = bitset_count(feasible_bits);
        if (feasible_count < need) {
            res.proved_infeasible = true;
            pair_sat_cache.emplace(std::move(cache_key), PairFillCacheEntry{2, {}});
            pair_sat_cache_stores++;
            return res;
        }
        if (feasible_count > 6000) {
            return res;
        }
        vector<int> feasible = bitset_to_indices(feasible_bits);
        CaDiCaL::Solver solver;
        solver.set("quiet", 1);
        solver.set("factor", 0);
        solver.limit("conflicts", 500000);
        solver.limit("decisions", 3000000);
        Clock::time_point pair_sat_deadline = min(deadline, Clock::now() + chrono::milliseconds(6000));
        DeadlineTerminator pair_sat_terminator(pair_sat_deadline);
        solver.connect_terminator(&pair_sat_terminator);
        int var_count = (int)feasible.size();
        long long clauses = 0;
        auto add_clause = [&](const vector<int> &lits) {
            for (int lit : lits) solver.add(lit);
            solver.add(0);
            clauses++;
        };
        auto add_at_most_k = [&](const vector<int> &vars, int max_true) {
            int m = (int)vars.size();
            if (max_true < 0) {
                add_clause({});
                return;
            }
            if (max_true >= m) return;
            if (max_true == 0) {
                for (int var : vars) add_clause({-var});
                return;
            }
            vector<vector<int>> seq(m, vector<int>(max_true + 1, 0));
            for (int i = 1; i < m; i++) {
                for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
            }
            for (int i = 1; i < m; i++) add_clause({-vars[i - 1], seq[i][1]});
            for (int i = 2; i < m; i++) {
                for (int j = 1; j <= max_true; j++) add_clause({-seq[i - 1][j], seq[i][j]});
            }
            for (int i = 2; i <= m; i++) add_clause({-vars[i - 1], -seq[i - 1][max_true]});
            for (int i = 2; i < m; i++) {
                for (int j = 2; j <= max_true; j++) {
                    add_clause({-vars[i - 1], -seq[i - 1][j - 1], seq[i][j]});
                }
            }
        };
        auto add_at_least_k = [&](const vector<int> &vars, int min_true) {
            int m = (int)vars.size();
            if (min_true <= 0) return;
            if (min_true > m) {
                add_clause({});
                return;
            }
            vector<vector<int>> seq(m, vector<int>(min_true + 1, 0));
            for (int i = 0; i < m; i++) {
                for (int j = 1; j <= min(min_true, i + 1); j++) {
                    seq[i][j] = ++var_count;
                }
            }
            add_clause({-vars[0], seq[0][1]});
            add_clause({-seq[0][1], vars[0]});
            for (int i = 1; i < m; i++) {
                add_clause({-vars[i], seq[i][1]});
                add_clause({-seq[i - 1][1], seq[i][1]});
                add_clause({-seq[i][1], vars[i], seq[i - 1][1]});
                for (int j = 2; j <= min_true; j++) {
                    if (!seq[i][j]) continue;
                    int prev_same = seq[i - 1][j];
                    int prev_less = seq[i - 1][j - 1];
                    if (prev_same) add_clause({-prev_same, seq[i][j]});
                    if (prev_less) {
                        add_clause({-vars[i], -prev_less, seq[i][j]});
                    }
                    vector<int> reverse1{-seq[i][j]};
                    if (prev_same) reverse1.push_back(prev_same);
                    if (prev_less) reverse1.push_back(prev_less);
                    add_clause(reverse1);
                    vector<int> reverse2{-seq[i][j]};
                    if (prev_same) reverse2.push_back(prev_same);
                    reverse2.push_back(vars[i]);
                    add_clause(reverse2);
                }
            }
            add_clause({seq[m - 1][min_true]});
        };
        int total_bits = inst.leaf_count + edge_offset.back();
        vector<vector<int>> by_bit(total_bits);
        vector<int> vars;
        vars.reserve(feasible.size());
        for (int pos = 0; pos < (int)feasible.size(); pos++) {
            int var = pos + 1;
            vars.push_back(var);
            const vector<unsigned long long> &mask = pairs[feasible[pos]].mask;
            for (int w = 0; w < (int)mask.size(); w++) {
                unsigned long long x = mask[w];
                while (x) {
                    int bit = __builtin_ctzll(x);
                    int global_bit = (w << 6) + bit;
                    if (global_bit < total_bits) by_bit[global_bit].push_back(var);
                    x &= x - 1;
                }
            }
        }
        for (const vector<int> &bucket : by_bit) {
            if (bucket.size() <= 1) continue;
            add_at_most_k(bucket, 1);
        }
        add_at_least_k(vars, need);
        cerr << "SMALL_SAVINGS_PAIR_SAT_SOLVE need=" << need
             << " feasible=" << feasible_count
             << " vars=" << var_count
             << " clauses=" << clauses << "\n";
        int sat = solver.solve();
        if (sat == 20) {
            res.proved_infeasible = true;
        } else if (sat == 10) {
            for (int pos = 0; pos < (int)feasible.size(); pos++) {
                if (solver.val(pos + 1) > 0) res.chosen.push_back(feasible[pos]);
            }
            if ((int)res.chosen.size() >= need) {
                vector<int> filtered;
                vector<unsigned long long> check_used(words, 0);
                for (int idx : res.chosen) {
                    if (intersects(check_used, pairs[idx].mask)) continue;
                    bit_or_into(check_used, pairs[idx].mask);
                    filtered.push_back(idx);
                    if ((int)filtered.size() == need) break;
                }
                if ((int)filtered.size() == need) {
                    res.chosen = std::move(filtered);
                    res.hit = true;
                } else {
                    res.chosen.clear();
                }
            } else {
                res.chosen.clear();
            }
        } else {
            res.hit = false;
        }
        cerr << "SMALL_SAVINGS_PAIR_SAT need=" << need
             << " feasible=" << feasible_count
             << " result=" << sat
             << " chosen=" << res.chosen.size()
             << " hit=" << (res.hit ? 1 : 0)
             << " clauses=" << clauses << "\n";
        if (res.hit) {
            pair_sat_cache.emplace(std::move(cache_key), PairFillCacheEntry{1, res.chosen});
            pair_sat_cache_stores++;
        } else if (res.proved_infeasible) {
            pair_sat_cache.emplace(std::move(cache_key), PairFillCacheEntry{2, {}});
            pair_sat_cache_stores++;
        }
#else
        (void)feasible_bits;
#endif
        return res;
    };

    auto fill_pairs_clique = [&](int need, const vector<unsigned long long> &feasible_bits) {
        struct CliquePackResult {
            bool hit = false;
            bool proved_infeasible = false;
            vector<int> chosen;
            long long nodes = 0;
        };
        CliquePackResult res;
        if (need <= 0) {
            res.hit = true;
            return res;
        }
        string cache_key = pair_cache_key_bits(need, feasible_bits);
        auto cache_it = pair_clique_cache.find(cache_key);
        if (cache_it != pair_clique_cache.end()) {
            pair_clique_cache_hits++;
            if (cache_it->second.state == 1) {
                res.hit = true;
                res.chosen = cache_it->second.chosen;
            } else if (cache_it->second.state == 2) {
                res.proved_infeasible = true;
            }
            return res;
        }
        int feasible_count = bitset_count(feasible_bits);
        if (feasible_count < need) {
            res.proved_infeasible = true;
            pair_clique_cache.emplace(std::move(cache_key), PairFillCacheEntry{2, {}});
            pair_clique_cache_stores++;
            return res;
        }
        vector<int> feasible = bitset_to_indices(feasible_bits);
        int m = (int)feasible.size();
        int bit_words = (m + 63) >> 6;
        vector<vector<unsigned long long>> compat(m, vector<unsigned long long>(bit_words, 0));
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                if (intersects(pairs[feasible[i]].mask, pairs[feasible[j]].mask)) continue;
                compat[i][j >> 6] |= 1ULL << (j & 63);
                compat[j][i >> 6] |= 1ULL << (i & 63);
            }
        }

        auto is_compat = [&](int a, int b) {
            return ((compat[a][b >> 6] >> (b & 63)) & 1ULL) != 0;
        };
        auto color_sort = [&](const vector<int> &input, vector<int> &ordered, vector<int> &colors) {
            vector<int> remaining = input;
            ordered.clear();
            colors.clear();
            int color = 0;
            while (!remaining.empty()) {
                color++;
                vector<int> next_remaining;
                vector<int> class_vertices;
                for (int v : remaining) {
                    bool conflicts_class = false;
                    for (int u : class_vertices) {
                        if (is_compat(v, u)) {
                            conflicts_class = true;
                            break;
                        }
                    }
                    if (conflicts_class) {
                        next_remaining.push_back(v);
                    } else {
                        class_vertices.push_back(v);
                    }
                }
                for (int v : class_vertices) {
                    ordered.push_back(v);
                    colors.push_back(color);
                }
                remaining.swap(next_remaining);
            }
        };

        vector<int> chosen_local;
        bool search_exhausted_local = true;
        function<bool(vector<int>, vector<int>, int)> expand =
            [&](vector<int> cand, vector<int> colors, int depth) -> bool {
                if (Clock::now() >= deadline || res.nodes++ >= pair_node_limit) {
                    search_exhausted_local = false;
                    return false;
                }
                while (!cand.empty()) {
                    if (depth + colors.back() < need) return false;
                    int v = cand.back();
                    cand.pop_back();
                    colors.pop_back();
                    chosen_local.push_back(v);
                    if (depth + 1 >= need) {
                        res.hit = true;
                        for (int local_idx : chosen_local) res.chosen.push_back(feasible[local_idx]);
                        return true;
                    }
                    vector<int> next;
                    next.reserve(cand.size());
                    for (int u : cand) {
                        if (is_compat(v, u)) next.push_back(u);
                    }
                    if (!next.empty()) {
                        vector<int> ordered, next_colors;
                        color_sort(next, ordered, next_colors);
                        if (expand(std::move(ordered), std::move(next_colors), depth + 1)) return true;
                        if (!search_exhausted_local) return false;
                    }
                    chosen_local.pop_back();
                }
                return false;
            };

        vector<int> initial(m);
        for (int i = 0; i < m; i++) initial[i] = i;
        sort(initial.begin(), initial.end(), [&](int a, int b) {
            int da = 0, db = 0;
            for (unsigned long long x : compat[a]) da += __builtin_popcountll(x);
            for (unsigned long long x : compat[b]) db += __builtin_popcountll(x);
            if (da != db) return da > db;
            return pairs[feasible[a]].leaves < pairs[feasible[b]].leaves;
        });
        vector<int> ordered, colors;
        color_sort(initial, ordered, colors);
        expand(std::move(ordered), std::move(colors), 0);
        res.proved_infeasible = !res.hit && search_exhausted_local;
        cerr << "SMALL_SAVINGS_PAIR_CLIQUE need=" << need
             << " feasible=" << feasible_count
             << " hit=" << (res.hit ? 1 : 0)
             << " proved_infeasible=" << (res.proved_infeasible ? 1 : 0)
             << " nodes=" << res.nodes
             << " exhausted=" << (search_exhausted_local ? 1 : 0)
             << " seconds=" << elapsed(start) << "\n";
        if (res.hit) {
            pair_clique_cache.emplace(std::move(cache_key), PairFillCacheEntry{1, res.chosen});
            pair_clique_cache_stores++;
        } else if (res.proved_infeasible) {
            pair_clique_cache.emplace(std::move(cache_key), PairFillCacheEntry{2, {}});
            pair_clique_cache_stores++;
        }
        return res;
    };

    function<bool(int, int, vector<unsigned long long>&)> fill_pairs =
        [&](int start_idx, int need, vector<unsigned long long> &used) {
            if (need <= 0) return true;
            if (Clock::now() >= deadline || pair_nodes++ > pair_node_limit) {
                search_exhausted = false;
                return false;
            }
            if ((int)pairs.size() - start_idx < need) return false;
            for (int pos = start_idx; pos < (int)pairs.size(); pos++) {
                if (intersects(used, pairs[pos].mask)) continue;
                bit_or_into(used, pairs[pos].mask);
                chosen_pairs.push_back(pos);
                if (fill_pairs(pos + 1, need - 1, used)) return true;
                chosen_pairs.pop_back();
                used = rebuild_used();
            }
            return false;
        };

    function<bool(int, int, vector<unsigned long long>&, vector<unsigned long long>&)> dfs_special =
        [&](int start_idx, int savings, vector<unsigned long long> &used,
                vector<unsigned long long> &pair_feasible_bits) {
            if (Clock::now() >= deadline || special_nodes++ > special_node_limit) {
                search_exhausted = false;
                return false;
            }
            int remaining = target_savings - savings;
            if (remaining <= 0) {
                answer_special = chosen_special;
                answer_pairs = chosen_pairs;
                return true;
            }
            chosen_pairs.clear();
            vector<unsigned long long> pair_used = used;
            int feasible_pair_count = bitset_count(pair_feasible_bits);
            bool pair_fill_closed = false;
            if (feasible_pair_count <= 1000) {
                auto pair_clique = fill_pairs_clique(remaining, pair_feasible_bits);
                if (pair_clique.hit) {
                    chosen_pairs = pair_clique.chosen;
                    answer_special = chosen_special;
                    answer_pairs = chosen_pairs;
                    return true;
                }
                pair_fill_closed = pair_clique.proved_infeasible;
            }
            if (!pair_fill_closed) {
                auto pair_sat = fill_pairs_sat(remaining, pair_feasible_bits);
                if (pair_sat.hit) {
                    chosen_pairs = pair_sat.chosen;
                    answer_special = chosen_special;
                    answer_pairs = chosen_pairs;
                    return true;
                }
                pair_fill_closed = pair_sat.proved_infeasible;
            }
            if (!pair_fill_closed) {
                auto pair_clique = fill_pairs_clique(remaining, pair_feasible_bits);
                if (pair_clique.hit) {
                    chosen_pairs = pair_clique.chosen;
                    answer_special = chosen_special;
                    answer_pairs = chosen_pairs;
                    return true;
                }
                pair_fill_closed = pair_clique.proved_infeasible;
            }
            if (!pair_fill_closed) {
                if (fill_pairs(0, remaining, pair_used)) {
                    answer_special = chosen_special;
                    answer_pairs = chosen_pairs;
                    return true;
                }
                chosen_pairs.clear();
                if (!search_exhausted) return false;
            }
            chosen_pairs.clear();
            for (int pos = start_idx; pos < (int)special.size(); pos++) {
                int add = (int)special[pos].leaves.size() - 1;
                if (add <= 0 || add > remaining) continue;
                if (intersects(used, special[pos].mask)) continue;
                vector<unsigned long long> next_used = used;
                bit_or_into(next_used, special[pos].mask);
                vector<unsigned long long> next_pair_feasible = pair_feasible_bits;
                const vector<unsigned long long> &blocked = blocked_pairs_for_special(pos);
                for (int w = 0; w < pair_bit_words; w++) {
                    next_pair_feasible[w] &= ~blocked[w];
                }
                chosen_special.push_back(pos);
                if (dfs_special(pos + 1, savings + add, next_used, next_pair_feasible)) return true;
                chosen_special.pop_back();
                if (!search_exhausted) return false;
            }
            return false;
        };

    vector<unsigned long long> used(words, 0);
    vector<unsigned long long> initial_pair_feasible(pair_bit_words, ~0ULL);
    if (!initial_pair_feasible.empty() && ((int)pairs.size() & 63)) {
        initial_pair_feasible.back() =
            (1ULL << ((int)pairs.size() & 63)) - 1ULL;
    }
    bool hit = exhausted && dfs_special(0, 0, used, initial_pair_feasible);
    result.proved_infeasible = exhausted && search_exhausted && !hit;
    cerr << "SMALL_SAVINGS_EXACT target_savings=" << target_savings
         << " hit=" << (hit ? 1 : 0)
         << " proved_infeasible=" << (result.proved_infeasible ? 1 : 0)
         << " special=" << special.size()
         << " pairs=" << pairs.size()
         << " special_nodes=" << special_nodes
         << " pair_nodes=" << pair_nodes
         << " pair_sat_cache_hits=" << pair_sat_cache_hits
         << " pair_sat_cache_stores=" << pair_sat_cache_stores
         << " pair_clique_cache_hits=" << pair_clique_cache_hits
         << " pair_clique_cache_stores=" << pair_clique_cache_stores
         << " exhausted=" << (exhausted ? 1 : 0)
         << " search_exhausted=" << (search_exhausted ? 1 : 0)
         << " seconds=" << elapsed(start) << "\n";
    if (!hit) return result;

    vector<unsigned char> covered(inst.leaf_count + 1, 0);
    for (int idx : answer_special) {
        result.partition.push_back(special[idx].leaves);
        for (int leaf : special[idx].leaves) covered[leaf] = 1;
    }
    for (int idx : answer_pairs) {
        result.partition.push_back(pairs[idx].leaves);
        for (int leaf : pairs[idx].leaves) covered[leaf] = 1;
    }
    for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
        if (!covered[leaf]) result.partition.push_back(vector<int>{leaf});
    }
    sort(result.partition.begin(), result.partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    if ((int)result.partition.size() > target_components ||
            !validate_partition(inst, result.partition, edge_offset)) {
        result.partition.clear();
    }
    return result;
}

struct SavingsLbResult {
    bool proved = false;
    bool unknown = false;
    int component_lb = 1;
    double upper_savings = 0.0;
    long long calls = 0;
    long long triple_scans = 0;
    long long component_enums = 0;
    long long pair_nodes = 0;
    string reason;
};

struct SavingsComponentChoice {
    vector<int> leaves;
    vector<unsigned long long> edges;
    vector<unsigned long long> leaf_mask;
    int saving = 0;
};

struct LbClusterDfsResult {
    vector<int> leaves;
    string signature;
    bool nonempty = false;
};

static LbClusterDfsResult lb_collect_active_clusters_dfs(
        const Tree &tree, int node, const vector<unsigned char> &active,
        map<vector<int>, string> &clusters) {
    int label = tree.nodes[node].label;
    if (label >= 1) {
        if (!active[label]) return {};
        return LbClusterDfsResult{{label}, to_string(label), true};
    }
    vector<LbClusterDfsResult> kids;
    for (int child : tree.nodes[node].child) {
        LbClusterDfsResult got =
            lb_collect_active_clusters_dfs(tree, child, active, clusters);
        if (got.nonempty) kids.push_back(std::move(got));
    }
    if (kids.empty()) return {};
    if (kids.size() == 1) return kids[0];
    LbClusterDfsResult out;
    out.nonempty = true;
    out.leaves.reserve(kids[0].leaves.size() + kids[1].leaves.size());
    merge(kids[0].leaves.begin(), kids[0].leaves.end(),
          kids[1].leaves.begin(), kids[1].leaves.end(),
          back_inserter(out.leaves));
    out.signature = combine_sig(kids[0].signature, kids[1].signature);
    if (out.leaves.size() >= 2) clusters[out.leaves] = out.signature;
    return out;
}

static map<vector<int>, string> lb_collect_active_clusters(
        const Tree &tree, const vector<unsigned char> &active) {
    map<vector<int>, string> clusters;
    lb_collect_active_clusters_dfs(tree, tree.root, active, clusters);
    return clusters;
}

static string lb_restricted_newick_renumber(const Tree &tree, int node,
        const vector<unsigned char> &active, const vector<int> &old_to_new) {
    int label = tree.nodes[node].label;
    if (label >= 1) {
        if (!active[label]) return "";
        return to_string(old_to_new[label]);
    }
    vector<string> kids;
    for (int child : tree.nodes[node].child) {
        string s = lb_restricted_newick_renumber(tree, child, active, old_to_new);
        if (!s.empty()) kids.push_back(std::move(s));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return "(" + kids[0] + "," + kids[1] + ")";
}

static Instance lb_common_pendant_reduced_instance(const Instance &inst, int &removed) {
    removed = 0;
    vector<unsigned char> active(inst.leaf_count + 1, 1);
    active[0] = 0;
    while (true) {
        vector<map<vector<int>, string>> all;
        all.reserve(inst.tree_count);
        for (const Tree &tree : inst.trees) all.push_back(lb_collect_active_clusters(tree, active));
        int active_count = 0;
        for (int label = 1; label <= inst.leaf_count; label++) if (active[label]) active_count++;
        vector<int> best;
        for (const auto &entry : all[0]) {
            const vector<int> &cluster = entry.first;
            if ((int)cluster.size() < 2 || (int)cluster.size() >= active_count) continue;
            const string &sig = entry.second;
            bool ok = true;
            for (int t = 1; t < inst.tree_count; t++) {
                auto found = all[t].find(cluster);
                if (found == all[t].end() || found->second != sig) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;
            if (cluster.size() > best.size() ||
                    (cluster.size() == best.size() && cluster < best)) {
                best = cluster;
            }
        }
        if (best.empty()) break;
        for (int i = 1; i < (int)best.size(); i++) {
            if (active[best[i]]) {
                active[best[i]] = 0;
                removed++;
            }
        }
    }
    if (removed == 0) return inst;
    vector<int> old_to_new(inst.leaf_count + 1, 0);
    int new_label = 0;
    for (int label = 1; label <= inst.leaf_count; label++) {
        if (active[label]) old_to_new[label] = ++new_label;
    }
    vector<string> newicks;
    newicks.reserve(inst.tree_count);
    for (const Tree &tree : inst.trees) {
        string nw = lb_restricted_newick_renumber(tree, tree.root, active, old_to_new);
        if (nw.empty()) nw = "1";
        newicks.push_back(nw);
    }
    Instance reduced;
    reduced.tree_count = inst.tree_count;
    reduced.leaf_count = new_label;
    for (const string &nw : newicks) reduced.trees.push_back(parse_tree(nw, new_label));
    return reduced;
}

struct SavingsLbProver {
    const Instance &inst;
    const vector<int> &edge_offset;
    int n;
    int edge_words;
    int needed;
    int max_depth;
    int max_components_per_triple;
    long long pair_node_limit;
    Clock::time_point deadline;
    unordered_set<uint32_t> valid_triples;
    vector<array<int, 3>> valid_triple_list;
    vector<vector<unsigned long long>> pair_edges;
    vector<int> leaf_component_scarcity;
    int leaf_words;
    vector<vector<unsigned long long>> edge_leaf_masks;
    unordered_map<string, double> memo_leaf_upper;
    unordered_map<int, vector<SavingsComponentChoice>> leaf_component_cache;
    SavingsLbResult result;
    bool pair_pack_aborted = false;

    SavingsLbProver(const Instance &inst_, const vector<int> &edge_offset_,
            int needed_, double seconds, int max_depth_,
            int max_components_per_triple_, long long pair_node_limit_)
        : inst(inst_), edge_offset(edge_offset_), n(inst_.leaf_count),
          edge_words((edge_offset_.back() + 63) >> 6), needed(needed_),
          max_depth(max_depth_),
          max_components_per_triple(max_components_per_triple_),
          pair_node_limit(pair_node_limit_),
          deadline(Clock::now() + chrono::milliseconds((int)(seconds * 1000))) {
        leaf_words = (n + 64) >> 6;
        valid_triples = collect_valid_triples(inst);
        valid_triple_list.reserve(valid_triples.size());
        for (int a = 1; a <= n; a++) {
            for (int b = a + 1; b <= n; b++) {
                for (int c = b + 1; c <= n; c++) {
                    if (valid_triples.count(triple_key(a, b, c, n))) {
                        valid_triple_list.push_back({a, b, c});
                    }
                }
            }
        }
        pair_edges.assign((n + 1) * (n + 1), {});
        for (int a = 1; a <= n; a++) {
            for (int b = a + 1; b <= n; b++) {
                vector<int> leaves{a, b};
                pair_edges[pair_index(a, b)] = edge_mask(leaves);
            }
        }
        edge_leaf_masks.assign(edge_offset.back(), vector<unsigned long long>(leaf_words, 0));
        for (int t = 0; t < inst.tree_count; t++) {
            const Tree &tree = inst.trees[t];
            function<vector<unsigned long long>(int)> dfs_mask = [&](int node) {
                vector<unsigned long long> mask(leaf_words, 0);
                int label = tree.nodes[node].label;
                if (label >= 1) {
                    mask[label >> 6] |= 1ULL << (label & 63);
                } else {
                    for (int child : tree.nodes[node].child) {
                        vector<unsigned long long> child_mask = dfs_mask(child);
                        bit_or_into(mask, child_mask);
                    }
                }
                edge_leaf_masks[edge_offset[t] + node] = mask;
                return mask;
            };
            dfs_mask(tree.root);
        }
        leaf_component_scarcity.assign(n + 1, 0);
        vector<int> current;
        function<void(int)> dfs_count = [&](int start) {
            if ((int)current.size() >= 2) {
                for (int leaf : current) {
                    if (leaf_component_scarcity[leaf] < INT_MAX / 4) {
                        leaf_component_scarcity[leaf]++;
                    }
                }
            }
            if ((int)current.size() >= min(7, needed + 1)) return;
            for (int leaf = start; leaf <= n; leaf++) {
                if (!can_add_leaf(current, leaf, unordered_set<uint32_t>{})) continue;
                current.push_back(leaf);
                dfs_count(leaf + 1);
                current.pop_back();
            }
        };
        dfs_count(1);
    }

    int pair_index(int a, int b) const {
        if (b < a) swap(a, b);
        return a * (n + 1) + b;
    }

    bool expired() const {
        return Clock::now() >= deadline;
    }

    void mark_unknown(const string &why) {
        result.unknown = true;
        if (result.reason.empty()) result.reason = why;
    }

    vector<unsigned long long> edge_mask(const vector<int> &leaves) const {
        vector<unsigned long long> mask(edge_words, 0);
        vector<int> edges = component_edges_all(inst, leaves, edge_offset);
        for (int edge : edges) set_bit(mask, edge);
        return mask;
    }

    bool edge_disjoint(const vector<unsigned long long> &a,
            const vector<unsigned long long> &b) const {
        return !intersects(a, b);
    }

    vector<unsigned long long> edge_union(vector<unsigned long long> a,
            const vector<unsigned long long> &b) const {
        bit_or_into(a, b);
        return a;
    }

    string leaf_branch_key(const vector<int> &active,
            const vector<unsigned long long> &used_edges,
            int need, int depth) const {
        vector<unsigned long long> active_bits(leaf_words, 0);
        for (int leaf : active) active_bits[leaf >> 6] |= 1ULL << (leaf & 63);
        string key;
        key.reserve(sizeof(int) * 2 +
            sizeof(unsigned long long) * (active_bits.size() + used_edges.size()));
        auto append_int = [&](int value) {
            key.append(reinterpret_cast<const char *>(&value), sizeof(value));
        };
        auto append_word = [&](unsigned long long value) {
            key.append(reinterpret_cast<const char *>(&value), sizeof(value));
        };
        append_int(need);
        append_int(max_depth - depth);
        for (unsigned long long word : active_bits) append_word(word);
        for (unsigned long long word : used_edges) append_word(word);
        return key;
    }

    vector<unsigned long long> active_leaf_bits(const vector<int> &active) const {
        vector<unsigned long long> bits(leaf_words, 0);
        for (int leaf : active) bits[leaf >> 6] |= 1ULL << (leaf & 63);
        return bits;
    }

    int active_intersection_count(const vector<unsigned long long> &a,
            const vector<unsigned long long> &b) const {
        int total = 0;
        for (int i = 0; i < (int)a.size(); i++) {
            total += __builtin_popcountll(a[i] & b[i]);
        }
        return total;
    }

    bool forbidden_contains(const unordered_set<uint32_t> &forbidden,
            int a, int b, int c) const {
        if (b < a) swap(a, b);
        if (c < b) swap(b, c);
        if (b < a) swap(a, b);
        return forbidden.count(triple_key(a, b, c, n));
    }

    bool triple_allowed(const unordered_set<uint32_t> &forbidden,
            int a, int b, int c) const {
        if (b < a) swap(a, b);
        if (c < b) swap(b, c);
        if (b < a) swap(a, b);
        uint32_t key = triple_key(a, b, c, n);
        return !forbidden.count(key) && valid_triples.count(key);
    }

    bool can_add_leaf(const vector<int> &current, int leaf,
            const unordered_set<uint32_t> &forbidden) const {
        if (current.size() < 2) return true;
        for (int i = 0; i < (int)current.size(); i++) {
            for (int j = i + 1; j < (int)current.size(); j++) {
                if (!triple_allowed(forbidden, current[i], current[j], leaf)) return false;
            }
        }
        return true;
    }

    array<int, 3> find_valid_triple(const vector<int> &active,
            const unordered_set<uint32_t> &forbidden) {
        result.triple_scans++;
        vector<unsigned char> is_active(n + 1, 0);
        for (int leaf : active) is_active[leaf] = 1;
        for (auto tri : valid_triple_list) {
            if (!is_active[tri[0]] || !is_active[tri[1]] || !is_active[tri[2]]) continue;
            if (forbidden.count(triple_key(tri[0], tri[1], tri[2], n))) continue;
            return tri;
        }
        return {-1, -1, -1};
    }

    int pair_pack_dfs(vector<unsigned char> &active,
            const vector<unsigned long long> &used_edges, int current,
            int &best, int cap) {
        if (expired()) {
            mark_unknown("timeout_pair_packing");
            return best;
        }
        result.pair_nodes++;
        if (result.pair_nodes > pair_node_limit) {
            pair_pack_aborted = true;
            return best;
        }
        int remaining = 0;
        for (int i = 1; i <= n; i++) if (active[i]) remaining++;
        if (current + remaining / 2 <= best) return best;
        if (best >= cap) return best;

        int pivot = -1;
        int best_options = INT_MAX;
        for (int a = 1; a <= n; a++) {
            if (!active[a]) continue;
            int options = 0;
            for (int b = 1; b <= n; b++) {
                if (a == b || !active[b]) continue;
                const vector<unsigned long long> &edges = pair_edges[pair_index(a, b)];
                if (edge_disjoint(used_edges, edges)) options++;
            }
            if (options < best_options) {
                best_options = options;
                pivot = a;
                if (options == 0) break;
            }
        }
        if (pivot < 0) {
            best = max(best, current);
            return best;
        }

        active[pivot] = 0;
        pair_pack_dfs(active, used_edges, current, best, cap);
        if (best >= cap) {
            active[pivot] = 1;
            return best;
        }
        vector<int> partners;
        for (int b = 1; b <= n; b++) {
            if (!active[b]) continue;
            const vector<unsigned long long> &edges = pair_edges[pair_index(pivot, b)];
            if (edge_disjoint(used_edges, edges)) partners.push_back(b);
        }
        sort(partners.begin(), partners.end(), [&](int x, int y) {
            int dx = bit_count(pair_edges[pair_index(pivot, x)]);
            int dy = bit_count(pair_edges[pair_index(pivot, y)]);
            if (dx != dy) return dx < dy;
            return x < y;
        });
        for (int b : partners) {
            if (!active[b]) continue;
            active[b] = 0;
            pair_pack_dfs(active, used_edges, current + 1, best, cap);
            active[b] = 1;
            if (best >= cap || result.unknown || pair_pack_aborted) break;
        }
        active[pivot] = 1;
        return best;
    }

    double pair_packing_upper(const vector<int> &active,
            const vector<unsigned long long> &used_edges, int cap) {
        int trivial = min(cap, (int)active.size() / 2);
        if (trivial <= 1) return trivial;
        pair_pack_aborted = false;
        vector<unsigned char> active_flag(n + 1, 0);
        for (int leaf : active) active_flag[leaf] = 1;
        int best = 0;
        pair_pack_dfs(active_flag, used_edges, 0, best, cap);
        if (pair_pack_aborted || (result.unknown && result.reason == "timeout_pair_packing")) {
            return trivial;
        }
        return min(cap, best);
    }

    vector<SavingsComponentChoice> enumerate_components_containing(
            const vector<int> &active,
            const vector<unsigned long long> &used_edges,
            const unordered_set<uint32_t> &forbidden,
            const array<int, 3> &base,
            int limit_size,
            bool &blocked) {
        result.component_enums++;
        blocked = false;
        vector<unsigned char> in_base(n + 1, 0);
        for (int x : base) in_base[x] = 1;
        vector<int> others;
        for (int leaf : active) if (!in_base[leaf]) others.push_back(leaf);
        vector<SavingsComponentChoice> out;
        vector<int> current{base[0], base[1], base[2]};
        sort(current.begin(), current.end());

        function<void(int)> dfs = [&](int pos) {
            if (blocked || result.unknown) return;
            if (expired()) {
                mark_unknown("timeout_component_enum");
                return;
            }
            if ((int)current.size() >= limit_size) {
                blocked = true;
                return;
            }
            vector<unsigned long long> edges = edge_mask(current);
            if (!edge_disjoint(used_edges, edges)) return;
            vector<unsigned long long> leaf_bits((n + 64) >> 6, 0);
            for (int leaf : current) leaf_bits[leaf >> 6] |= 1ULL << (leaf & 63);
            out.push_back(SavingsComponentChoice{
                current, std::move(edges), std::move(leaf_bits), (int)current.size() - 1});
            if ((int)out.size() > max_components_per_triple) {
                mark_unknown("component_candidate_limit");
                return;
            }
            for (int i = pos; i < (int)others.size(); i++) {
                int leaf = others[i];
                if (!can_add_leaf(current, leaf, forbidden)) continue;
                current.push_back(leaf);
                sort(current.begin(), current.end());
                dfs(i + 1);
                current.erase(find(current.begin(), current.end(), leaf));
                if (blocked || result.unknown) return;
            }
        };
        dfs(0);
        return out;
    }

    double upper_savings(const vector<int> &active,
            const vector<unsigned long long> &used_edges,
            const unordered_set<uint32_t> &forbidden,
            int need, int depth) {
        if (need <= 0 || active.empty()) return 0.0;
        if (expired()) {
            mark_unknown("timeout");
            return need;
        }
        if (depth > max_depth) {
            mark_unknown("max_depth");
            return need;
        }
        result.calls++;
        int trivial_upper = (int)active.size() - 1;
        if (trivial_upper < need) return (double)trivial_upper;

        array<int, 3> tri = find_valid_triple(active, forbidden);
        if (tri[0] < 0) return pair_packing_upper(active, used_edges, need);

        double best = 0.0;
        bool blocked = false;
        vector<SavingsComponentChoice> choices = enumerate_components_containing(
            active, used_edges, forbidden, tri, need + 1, blocked);
        if (blocked || result.unknown) return need;

        vector<unsigned char> active_flag(n + 1, 0);
        for (int leaf : active) active_flag[leaf] = 1;
        for (const SavingsComponentChoice &choice : choices) {
            if (choice.saving >= need) return need;
            vector<int> remaining;
            remaining.reserve(active.size() - choice.leaves.size());
            for (int leaf : active) active_flag[leaf] = 1;
            for (int leaf : choice.leaves) active_flag[leaf] = 0;
            for (int leaf : active) if (active_flag[leaf]) remaining.push_back(leaf);
            vector<unsigned long long> next_edges = edge_union(used_edges, choice.edges);
            double child = upper_savings(
                remaining, next_edges, forbidden, need - choice.saving, depth + 1);
            double candidate = choice.saving + child;
            if (candidate > best) best = candidate;
            if (best >= need) return need;
        }

        unordered_set<uint32_t> next_forbidden = forbidden;
        next_forbidden.insert(triple_key(tri[0], tri[1], tri[2], n));
        double forbid_upper = upper_savings(active, used_edges, next_forbidden, need, depth + 1);
        if (forbid_upper > best) best = forbid_upper;
        return best;
    }

    int choose_pivot_leaf(const vector<int> &active) const {
        vector<unsigned char> is_active(n + 1, 0);
        for (int leaf : active) is_active[leaf] = 1;
        int best_leaf = active.empty() ? -1 : active[0];
        int best_score = INT_MAX;
        int best_scarcity = INT_MAX;
        for (int leaf : active) {
            int triples = 0;
            for (auto tri : valid_triple_list) {
                if (tri[0] != leaf && tri[1] != leaf && tri[2] != leaf) continue;
                if (is_active[tri[0]] && is_active[tri[1]] && is_active[tri[2]]) triples++;
            }
            int scarcity = leaf_component_scarcity.empty() ? triples : leaf_component_scarcity[leaf];
            if (scarcity < best_scarcity ||
                    (scarcity == best_scarcity && triples < best_score) ||
                    (scarcity == best_scarcity && triples == best_score && leaf < best_leaf)) {
                best_scarcity = scarcity;
                best_score = triples;
                best_leaf = leaf;
            }
        }
        return best_leaf;
    }

    vector<SavingsComponentChoice> &all_components_containing_leaf(
            int pivot, int limit_size) {
        auto found = leaf_component_cache.find(pivot);
        if (found != leaf_component_cache.end()) return found->second;
        result.component_enums++;
        vector<SavingsComponentChoice> cached;
        vector<int> others;
        for (int leaf = 1; leaf <= n; leaf++) if (leaf != pivot) others.push_back(leaf);
        vector<int> current{pivot};

        function<void(int)> dfs = [&](int pos) {
            if (result.unknown) return;
            if (expired()) {
                mark_unknown("timeout_leaf_component_enum");
                return;
            }
            if ((int)current.size() > 1) {
                vector<unsigned long long> edges = edge_mask(current);
                if (!valid_topology(inst, current)) return;
                vector<unsigned long long> leaf_bits((n + 64) >> 6, 0);
                for (int leaf : current) leaf_bits[leaf >> 6] |= 1ULL << (leaf & 63);
                cached.push_back(SavingsComponentChoice{
                    current, std::move(edges), std::move(leaf_bits), (int)current.size() - 1});
                if ((int)cached.size() > max_components_per_triple) {
                    mark_unknown("leaf_component_candidate_limit");
                    return;
                }
            }
            if ((int)current.size() >= limit_size) return;
            for (int i = pos; i < (int)others.size(); i++) {
                int leaf = others[i];
                if (!can_add_leaf(current, leaf, unordered_set<uint32_t>{})) continue;
                current.push_back(leaf);
                sort(current.begin(), current.end());
                dfs(i + 1);
                current.erase(find(current.begin(), current.end(), leaf));
                if (result.unknown) return;
            }
        };
        dfs(0);
        sort(cached.begin(), cached.end(), [](const SavingsComponentChoice &a,
                    const SavingsComponentChoice &b) {
            if (a.saving != b.saving) return a.saving > b.saving;
            if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
            return a.leaves < b.leaves;
        });
        cerr << "SAVINGS_LEAF_CACHE leaf=" << pivot
             << " components=" << cached.size() << "\n";
        auto inserted = leaf_component_cache.emplace(pivot, std::move(cached));
        return inserted.first->second;
    }

    vector<SavingsComponentChoice> enumerate_components_containing_leaf(
            const vector<int> &active,
            const vector<unsigned long long> &used_edges,
            int pivot,
            int limit_size) {
        vector<SavingsComponentChoice> &cached =
            all_components_containing_leaf(pivot, needed + 1);
        if (result.unknown) return {};
        vector<unsigned long long> active_bits((n + 64) >> 6, 0);
        for (int leaf : active) active_bits[leaf >> 6] |= 1ULL << (leaf & 63);
        vector<SavingsComponentChoice> out;
        for (const SavingsComponentChoice &choice : cached) {
            if ((int)choice.leaves.size() > limit_size) continue;
            bool subset = true;
            for (int w = 0; w < (int)choice.leaf_mask.size(); w++) {
                if (choice.leaf_mask[w] & ~active_bits[w]) {
                    subset = false;
                    break;
                }
            }
            if (!subset) continue;
            if (!edge_disjoint(used_edges, choice.edges)) continue;
            out.push_back(choice);
        }
        return out;
    }

    int edge_capacity_upper(const vector<unsigned long long> &used_edges,
            const vector<unsigned long long> &active_bits,
            int active_count) const {
        if (active_count <= 1) return 0;
        int upper = max(0, active_count - 1);
        for (int t = 0; t < inst.tree_count; t++) {
            int free_edges = 0;
            for (int edge = edge_offset[t]; edge < edge_offset[t + 1]; edge++) {
                int word = edge >> 6;
                int bit = edge & 63;
                if ((used_edges[word] & (1ULL << bit)) == 0) free_edges++;
                else continue;
                int below = active_intersection_count(edge_leaf_masks[edge], active_bits);
                if (below <= 0 || below >= active_count) {
                    free_edges--;
                }
            }
            upper = min(upper, free_edges / 2);
        }
        return upper;
    }

    int free_pair_connectivity_upper(const vector<int> &active,
            const vector<unsigned long long> &used_edges) const {
        int m = (int)active.size();
        if (m <= 1) return 0;
        vector<int> parent(m);
        iota(parent.begin(), parent.end(), 0);
        function<int(int)> find_root = [&](int x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        auto unite = [&](int a, int b) {
            int ra = find_root(a), rb = find_root(b);
            if (ra != rb) parent[rb] = ra;
        };
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                const vector<unsigned long long> &edges =
                    pair_edges[pair_index(active[i], active[j])];
                if (edge_disjoint(used_edges, edges)) unite(i, j);
            }
        }
        int components = 0;
        for (int i = 0; i < m; i++) {
            if (find_root(i) == i) components++;
        }
        return m - components;
    }

    int combined_packing_upper(const vector<int> &active,
            const vector<unsigned long long> &used_edges) const {
        if (active.size() <= 1) return 0;
        vector<unsigned long long> active_bits_for_bound = active_leaf_bits(active);
        int edge_upper = edge_capacity_upper(
            used_edges, active_bits_for_bound, (int)active.size());
        int connectivity_upper = free_pair_connectivity_upper(active, used_edges);
        return min(edge_upper, connectivity_upper);
    }

    double upper_savings_leaf_branch(const vector<int> &active,
            const vector<unsigned long long> &used_edges,
            int need, int depth) {
        if (need <= 0 || active.empty()) return 0.0;
        if (expired()) {
            mark_unknown("timeout");
            return need;
        }
        if (depth > max_depth) {
            mark_unknown("max_depth");
            return need;
        }
        result.calls++;
        int trivial_upper = (int)active.size() - 1;
        if (trivial_upper < need) return (double)trivial_upper;
        if ((int)active.size() <= 2) return min(need, (int)active.size() - 1);
        int packing_upper = combined_packing_upper(active, used_edges);
        if (packing_upper < need) return (double)packing_upper;
        string memo_key = leaf_branch_key(active, used_edges, need, depth);
        auto memo_found = memo_leaf_upper.find(memo_key);
        if (memo_found != memo_leaf_upper.end()) return memo_found->second;

        int pivot = choose_pivot_leaf(active);
        vector<int> singleton_remaining;
        singleton_remaining.reserve(active.size() - 1);
        for (int leaf : active) if (leaf != pivot) singleton_remaining.push_back(leaf);
        double best = 0.0;
        int singleton_cap = combined_packing_upper(singleton_remaining, used_edges);
        if (singleton_cap < need) {
            best = singleton_cap;
        } else {
            best = upper_savings_leaf_branch(
                singleton_remaining, used_edges, need, depth + 1);
        }
        if (best >= need) return need;

        vector<SavingsComponentChoice> choices = enumerate_components_containing_leaf(
            active, used_edges, pivot, need + 1);
        if (result.unknown) return need;

        vector<unsigned char> active_flag(n + 1, 0);
        for (int leaf : active) active_flag[leaf] = 1;
        for (const SavingsComponentChoice &choice : choices) {
            if (choice.saving >= need) return need;
            for (int leaf : active) active_flag[leaf] = 1;
            for (int leaf : choice.leaves) active_flag[leaf] = 0;
            vector<int> remaining;
            remaining.reserve(active.size() - choice.leaves.size());
            for (int leaf : active) if (active_flag[leaf]) remaining.push_back(leaf);
            vector<unsigned long long> next_edges = edge_union(used_edges, choice.edges);
            int child_cap = combined_packing_upper(remaining, next_edges);
            int branch_cap = choice.saving + child_cap;
            if (branch_cap <= best + 1e-9) continue;
            if (branch_cap < need) {
                best = max(best, (double)branch_cap);
                continue;
            }
            double child = upper_savings_leaf_branch(
                remaining, next_edges, need - choice.saving, depth + 1);
            double candidate = choice.saving + child;
            if (candidate > best) best = candidate;
            if (best >= need) {
                memo_leaf_upper.emplace(std::move(memo_key), (double)need);
                return need;
            }
        }
        memo_leaf_upper.emplace(std::move(memo_key), best);
        return best;
    }

    SavingsLbResult run() {
        vector<int> active(n);
        for (int i = 0; i < n; i++) active[i] = i + 1;
        vector<unsigned long long> used_edges(edge_words, 0);
        unordered_set<uint32_t> forbidden;
        (void)forbidden;
        result.upper_savings = upper_savings_leaf_branch(active, used_edges, needed, 0);
        if (!result.unknown && result.upper_savings < needed - 1e-8) {
            result.proved = true;
            int max_savings = (int)floor(result.upper_savings + 1e-8);
            result.component_lb = n - max_savings;
        }
        return result;
    }
};

static SavingsLbResult savings_lb_rule_out(const Instance &inst,
        int target_components, double seconds, int max_depth,
        int max_components_per_triple, long long pair_node_limit) {
    SavingsLbResult result;
    if (target_components < 1) {
        result.proved = true;
        result.component_lb = 1;
        return result;
    }
    int needed = inst.leaf_count - target_components;
    if (seconds <= 0.0 || needed <= 0) {
        result.unknown = true;
        result.reason = needed <= 0 ? "target_not_below_n" : "no_time";
        return result;
    }
    int reduced_removed = 0;
    Instance reduced = lb_common_pendant_reduced_instance(inst, reduced_removed);
    const Instance *work = &inst;
    if (reduced_removed > 0) {
        work = &reduced;
        needed = work->leaf_count - target_components;
        cerr << "SAVINGS_LB_REDUCED original_n=" << inst.leaf_count
             << " reduced_n=" << work->leaf_count
             << " pendant_removed=" << reduced_removed
             << " target=" << target_components
             << " needed_savings=" << needed << "\n";
        if (needed <= 0) {
            result.unknown = true;
            result.reason = "target_not_below_reduced_n";
            return result;
        }
    }
    vector<int> edge_offset(inst.tree_count + 1, 0);
    for (int t = 0; t < work->tree_count; t++) {
        edge_offset[t + 1] = edge_offset[t] + (int)work->trees[t].nodes.size();
    }
    SavingsLbProver prover(
        *work, edge_offset, needed, seconds, max_depth,
        max_components_per_triple, pair_node_limit);
    result = prover.run();
    cerr << "SAVINGS_LB_STATUS target=" << target_components
         << " needed_savings=" << needed
         << " result=" << (result.proved ? "NO" : result.unknown ? "UNKNOWN" : "WEAK")
         << " component_lb=" << result.component_lb
         << " upper_savings=" << result.upper_savings
         << " calls=" << result.calls
         << " triple_scans=" << result.triple_scans
         << " component_enums=" << result.component_enums
         << " pair_nodes=" << result.pair_nodes
         << " reason=" << (result.reason.empty() ? "none" : result.reason)
         << "\n";
    return result;
}

#if defined(USE_HIGHS) || defined(USE_CADICAL)
struct ClauseKey {
    int a, b, c;
    bool operator==(const ClauseKey &other) const {
        return a == other.a && b == other.b && c == other.c;
    }
};
struct ClauseKeyHash {
    size_t operator()(const ClauseKey &key) const {
        uint64_t x = (uint32_t)key.a;
        x = x * 1000003ULL ^ (uint32_t)key.b;
        x = x * 1000003ULL ^ (uint32_t)key.c;
        return (size_t)x;
    }
};
#endif

#ifdef USE_CADICAL
class TimeTerminator : public CaDiCaL::Terminator {
    Clock::time_point deadline;
public:
    explicit TimeTerminator(Clock::time_point deadline_) : deadline(deadline_) {}
    bool terminate() override { return Clock::now() >= deadline; }
};

static vector<vector<int>> direct_pair_cadical_target(const Instance &inst, int target_components,
        double seconds, bool *proved_unsat = nullptr) {
    if (seconds <= 0.0) return {};
    if (proved_unsat) *proved_unsat = false;
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    int n = inst.leaf_count;
    vector<vector<int>> pair_var(n + 1, vector<int>(n + 1, -1));
    int var_count = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            pair_var[a][b] = pair_var[b][a] = ++var_count;
        }
    }
    vector<int> rep_var(n + 1, -1);
    for (int i = 1; i <= n; i++) rep_var[i] = ++var_count;

    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    TimeTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    unordered_set<ClauseKey, ClauseKeyHash> seen_clauses;
    int clause_count = 0;

    auto add_clause = [&](vector<int> lits) {
        sort(lits.begin(), lits.end(), [](int x, int y) {
            if (abs(x) != abs(y)) return abs(x) < abs(y);
            return x < y;
        });
        lits.erase(unique(lits.begin(), lits.end()), lits.end());
        for (int i = 0; i < (int)lits.size(); i++) {
            for (int j = i + 1; j < (int)lits.size(); j++) {
                if (lits[i] == -lits[j]) return;
            }
        }
        if (lits.size() == 3) {
            ClauseKey key{lits[0], lits[1], lits[2]};
            if (!seen_clauses.insert(key).second) return;
        }
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clause_count++;
    };

    auto add_at_most_k = [&](const vector<int> &vars, int max_true) {
        int m = (int)vars.size();
        if (max_true < 0) {
            add_clause({});
            return;
        }
        if (max_true >= m) return;
        if (max_true == 0) {
            for (int var : vars) add_clause({-var});
            return;
        }
        vector<vector<int>> seq(m, vector<int>(max_true + 1, 0));
        for (int i = 1; i < m; i++) {
            for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
        }
        for (int i = 1; i < m; i++) add_clause({-vars[i - 1], seq[i][1]});
        for (int i = 2; i < m; i++) {
            for (int j = 1; j <= max_true; j++) add_clause({-seq[i - 1][j], seq[i][j]});
        }
        for (int i = 2; i <= m; i++) add_clause({-vars[i - 1], -seq[i - 1][max_true]});
        for (int i = 2; i < m; i++) {
            for (int j = 2; j <= max_true; j++) add_clause({-vars[i - 1], -seq[i - 1][j - 1], seq[i][j]});
        }
    };

    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                int ab = pair_var[a][b], ac = pair_var[a][c], bc = pair_var[b][c];
                add_clause({-ab, -bc, ac});
                add_clause({-ab, -ac, bc});
                add_clause({-ac, -bc, ab});
            }
        }
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    int conflict_triples = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                if (valid.count(triple_key(a, b, c, n))) continue;
                conflict_triples++;
                add_clause({-pair_var[a][b], -pair_var[a][c], -pair_var[b][c]});
            }
        }
    }

    add_clause({rep_var[1]});
    for (int i = 2; i <= n; i++) {
        vector<int> rep_or_previous{rep_var[i]};
        for (int j = 1; j < i; j++) {
            int same = pair_var[j][i];
            rep_or_previous.push_back(same);
            add_clause({-rep_var[i], -same});
        }
        add_clause(rep_or_previous);
    }
    vector<int> reps;
    for (int i = 1; i <= n; i++) reps.push_back(rep_var[i]);
    add_at_most_k(reps, target_components);

    int edge_clauses = 0;
    vector<int> all_leaves(n);
    for (int i = 0; i < n; i++) all_leaves[i] = i + 1;
    for (const Tree &tree : inst.trees) {
        for (int node = 0; node < (int)tree.nodes.size(); node++) {
            if (node == tree.root) continue;
            vector<int> side;
            function<void(int)> collect = [&](int u) {
                int label = tree.nodes[u].label;
                if (label >= 1) {
                    side.push_back(label);
                    return;
                }
                for (int child : tree.nodes[u].child) collect(child);
            };
            collect(node);
            sort(side.begin(), side.end());
            if (side.empty() || (int)side.size() == n) continue;
            vector<unsigned char> in_side(n + 1, 0);
            for (int x : side) in_side[x] = 1;
            vector<int> outside;
            for (int x : all_leaves) if (!in_side[x]) outside.push_back(x);
            if (side.size() > outside.size()) side.swap(outside);
            for (int ai = 0; ai < (int)side.size(); ai++) {
                int a = side[ai];
                for (int ci = ai + 1; ci < (int)side.size(); ci++) {
                    int c = side[ci];
                    int ac = pair_var[a][c];
                    for (int bi = 0; bi < (int)outside.size(); bi++) {
                        int b = outside[bi];
                        int ab = pair_var[a][b];
                        int cb = pair_var[c][b];
                        for (int di = bi + 1; di < (int)outside.size(); di++) {
                            int d = outside[di];
                            int cd = pair_var[c][d];
                            int ad = pair_var[a][d];
                            int before = clause_count;
                            add_clause({-ab, -cd, ac});
                            if (clause_count != before) edge_clauses++;
                            before = clause_count;
                            add_clause({-ad, -cb, ac});
                            if (clause_count != before) edge_clauses++;
                        }
                    }
                }
            }
        }
    }

    cerr << "DIRECT_SAT_MODEL vars=" << var_count
         << " clauses=" << clause_count
         << " conflict_triples=" << conflict_triples
         << " edge_clauses=" << edge_clauses
         << " build_seconds=" << elapsed(start) << "\n";
    int result = solver.solve();
    cerr << "DIRECT_SAT_STATUS result=" << result
         << " seconds=" << elapsed(start) << "\n";
    solver.disconnect_terminator();
    if (proved_unsat && result == 20) *proved_unsat = true;
    if (result != 10) return {};

    vector<int> parent(n + 1);
    for (int i = 1; i <= n; i++) parent[i] = i;
    function<int(int)> find_root = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find_root(a), rb = find_root(b);
        if (ra != rb) parent[rb] = ra;
    };
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            if (solver.val(pair_var[a][b]) > 0) unite(a, b);
        }
    }
    unordered_map<int, vector<int>> groups;
    for (int leaf = 1; leaf <= n; leaf++) groups[find_root(leaf)].push_back(leaf);
    vector<vector<int>> partition;
    for (auto &kv : groups) {
        sort(kv.second.begin(), kv.second.end());
        partition.push_back(kv.second);
    }
    sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return partition;
}

static bool direct_pair_cadical_rule_out(const Instance &inst, int target_components, double seconds) {
    if (target_components < 1) return true;
    bool proved_unsat = false;
    (void)direct_pair_cadical_target(inst, target_components, seconds, &proved_unsat);
    return proved_unsat;
}

static vector<vector<int>> edge_owner_cadical_target(const Instance &inst,
        int target_components, double seconds, bool *proved_unsat = nullptr) {
    if (proved_unsat) *proved_unsat = false;
    if (seconds <= 0.0 || target_components < 1) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    int n = inst.leaf_count;
    if (inst.tree_count >= 20 && n >= 90 && n <= 100 &&
            target_components >= 80 && target_components <= 85) {
        cerr << "DIRECT_LAZY_SAT_SKIPPED reason=memory_risk_band"
             << " t=" << inst.tree_count
             << " n=" << n
             << " target=" << target_components << "\n";
        return {};
    }
    vector<vector<int>> pair_var(n + 1, vector<int>(n + 1, -1));
    int var_count = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            pair_var[a][b] = pair_var[b][a] = ++var_count;
        }
    }
    vector<int> rep_var(n + 1, -1);
    for (int i = 1; i <= n; i++) rep_var[i] = ++var_count;

    struct EdgeSide {
        int tree_index = -1;
        int edge = -1;
        vector<int> side;
        vector<int> outside;
    };
    vector<EdgeSide> edge_sides;
    for (int tree_index = 0; tree_index < inst.tree_count; tree_index++) {
        const Tree &tree = inst.trees[tree_index];
        for (int edge = 0; edge < (int)tree.nodes.size(); edge++) {
            if (edge == tree.root) continue;
            vector<int> side;
            function<void(int)> collect = [&](int u) {
                int label = tree.nodes[u].label;
                if (label >= 1) {
                    side.push_back(label);
                    return;
                }
                for (int child : tree.nodes[u].child) collect(child);
            };
            collect(edge);
            sort(side.begin(), side.end());
            if (side.empty() || (int)side.size() == n) continue;
            vector<unsigned char> in_side(n + 1, 0);
            for (int x : side) in_side[x] = 1;
            vector<int> outside;
            for (int x = 1; x <= n; x++) if (!in_side[x]) outside.push_back(x);
            if (side.size() > outside.size()) side.swap(outside);
            if (side.empty() || outside.empty()) continue;
            edge_sides.push_back(EdgeSide{tree_index, edge, std::move(side), std::move(outside)});
        }
    }
    int edge_count = (int)edge_sides.size();
    vector<vector<int>> cross_var(edge_count, vector<int>(n + 1, -1));
    for (int e = 0; e < edge_count; e++) {
        for (int r = 1; r <= n; r++) cross_var[e][r] = ++var_count;
    }

    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factorcheck", 0);
    TimeTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    unordered_set<ClauseKey, ClauseKeyHash> seen_ternary;
    long long clause_count = 0;

    auto add_clause = [&](vector<int> lits) {
        sort(lits.begin(), lits.end(), [](int x, int y) {
            if (abs(x) != abs(y)) return abs(x) < abs(y);
            return x < y;
        });
        lits.erase(unique(lits.begin(), lits.end()), lits.end());
        for (int i = 0; i < (int)lits.size(); i++) {
            for (int j = i + 1; j < (int)lits.size(); j++) {
                if (lits[i] == -lits[j]) return;
            }
        }
        if (lits.size() == 3) {
            ClauseKey key{lits[0], lits[1], lits[2]};
            if (!seen_ternary.insert(key).second) return;
        }
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clause_count++;
    };

    auto add_at_most_k = [&](const vector<int> &vars, int max_true) {
        int m = (int)vars.size();
        if (max_true < 0) {
            add_clause({});
            return;
        }
        if (max_true >= m) return;
        if (max_true == 0) {
            for (int var : vars) add_clause({-var});
            return;
        }
        vector<vector<int>> seq(m, vector<int>(max_true + 1, 0));
        for (int i = 1; i < m; i++) {
            for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
        }
        for (int i = 1; i < m; i++) add_clause({-vars[i - 1], seq[i][1]});
        for (int i = 2; i < m; i++) {
            for (int j = 1; j <= max_true; j++) add_clause({-seq[i - 1][j], seq[i][j]});
        }
        for (int i = 2; i <= m; i++) add_clause({-vars[i - 1], -seq[i - 1][max_true]});
        for (int i = 2; i < m; i++) {
            for (int j = 2; j <= max_true; j++) {
                add_clause({-vars[i - 1], -seq[i - 1][j - 1], seq[i][j]});
            }
        }
    };

    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                int ab = pair_var[a][b], ac = pair_var[a][c], bc = pair_var[b][c];
                add_clause({-ab, -bc, ac});
                add_clause({-ab, -ac, bc});
                add_clause({-ac, -bc, ab});
            }
        }
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    int conflict_triples = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                if (valid.count(triple_key(a, b, c, n))) continue;
                conflict_triples++;
                add_clause({-pair_var[a][b], -pair_var[a][c], -pair_var[b][c]});
            }
        }
    }

    add_clause({rep_var[1]});
    for (int i = 2; i <= n; i++) {
        vector<int> rep_or_previous{rep_var[i]};
        for (int j = 1; j < i; j++) {
            int same = pair_var[j][i];
            rep_or_previous.push_back(same);
            add_clause({-rep_var[i], -same});
        }
        add_clause(rep_or_previous);
    }
    vector<int> reps;
    for (int i = 1; i <= n; i++) reps.push_back(rep_var[i]);
    add_at_most_k(reps, target_components);

    long long owner_clauses = 0;
    long long at_most_edge_clauses = 0;
    for (int e = 0; e < edge_count; e++) {
        if (Clock::now() >= deadline) break;
        const EdgeSide &es = edge_sides[e];
        for (int r = 1; r <= n; r++) {
            add_clause({-cross_var[e][r], rep_var[r]});
            for (int a : es.side) {
                int ra = (r == a) ? 0 : pair_var[r][a];
                for (int b : es.outside) {
                    vector<int> clause{-rep_var[r], -pair_var[a][b], cross_var[e][r]};
                    if (ra) clause.push_back(-ra);
                    add_clause(std::move(clause));
                    owner_clauses++;
                }
            }
        }
        for (int r1 = 1; r1 <= n; r1++) {
            for (int r2 = r1 + 1; r2 <= n; r2++) {
                add_clause({-cross_var[e][r1], -cross_var[e][r2]});
                at_most_edge_clauses++;
            }
        }
    }

    cerr << "EDGE_OWNER_MODEL vars=" << var_count
         << " clauses=" << clause_count
         << " edges=" << edge_count
         << " conflict_triples=" << conflict_triples
         << " owner_clauses=" << owner_clauses
         << " edge_at_most=" << at_most_edge_clauses
         << " build_seconds=" << elapsed(start) << "\n";
    int result = solver.solve();
    cerr << "EDGE_OWNER_STATUS result=" << result
         << " seconds=" << elapsed(start) << "\n";
    solver.disconnect_terminator();
    if (proved_unsat && result == 20) *proved_unsat = true;
    if (result != 10) return {};

    vector<int> parent(n + 1);
    for (int i = 1; i <= n; i++) parent[i] = i;
    function<int(int)> find_root = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find_root(a), rb = find_root(b);
        if (ra != rb) parent[rb] = ra;
    };
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            if (solver.val(pair_var[a][b]) > 0) unite(a, b);
        }
    }
    unordered_map<int, vector<int>> groups;
    for (int leaf = 1; leaf <= n; leaf++) groups[find_root(leaf)].push_back(leaf);
    vector<vector<int>> partition;
    for (auto &kv : groups) {
        sort(kv.second.begin(), kv.second.end());
        partition.push_back(kv.second);
    }
    sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return partition;
}

static bool leaf_below_node(const Tree &tree, int leaf, int ancestor) {
    int node = tree.leaf_node[leaf];
    while (node >= 0) {
        if (node == ancestor) return true;
        node = tree.nodes[node].parent;
    }
    return false;
}

static pair<int, int> crossing_pair_for_edge(
        const Tree &tree, const vector<int> &block, int edge_node) {
    int inside = -1, outside = -1;
    for (int leaf : block) {
        if (leaf_below_node(tree, leaf, edge_node)) {
            if (inside < 0) inside = leaf;
        } else {
            if (outside < 0) outside = leaf;
        }
        if (inside >= 0 && outside >= 0) return {inside, outside};
    }
    return {-1, -1};
}

static vector<vector<int>> direct_pair_cadical_lazy_target(const Instance &inst,
        int target_components, double seconds, bool *proved_unsat = nullptr,
        int preload_side_limit = 0) {
    if (proved_unsat) *proved_unsat = false;
    if (seconds <= 0.0 || target_components < 1) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    int n = inst.leaf_count;
    vector<vector<int>> pair_var(n + 1, vector<int>(n + 1, -1));
    int var_count = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            pair_var[a][b] = pair_var[b][a] = ++var_count;
        }
    }
    vector<int> rep_var(n + 1, -1);
    for (int i = 1; i <= n; i++) rep_var[i] = ++var_count;

    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    TimeTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    unordered_set<ClauseKey, ClauseKeyHash> seen_clauses;
    long long clause_count = 0;

    auto add_clause = [&](vector<int> lits) -> bool {
        sort(lits.begin(), lits.end(), [](int x, int y) {
            if (abs(x) != abs(y)) return abs(x) < abs(y);
            return x < y;
        });
        lits.erase(unique(lits.begin(), lits.end()), lits.end());
        for (int i = 0; i < (int)lits.size(); i++) {
            for (int j = i + 1; j < (int)lits.size(); j++) {
                if (lits[i] == -lits[j]) return false;
            }
        }
        if (lits.size() == 3) {
            ClauseKey key{lits[0], lits[1], lits[2]};
            if (!seen_clauses.insert(key).second) return false;
        }
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clause_count++;
        return true;
    };

    auto add_at_most_k = [&](const vector<int> &vars, int max_true) {
        int m = (int)vars.size();
        if (max_true < 0) {
            add_clause({});
            return;
        }
        if (max_true >= m) return;
        if (max_true == 0) {
            for (int var : vars) add_clause({-var});
            return;
        }
        vector<vector<int>> seq(m, vector<int>(max_true + 1, 0));
        for (int i = 1; i < m; i++) {
            for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
        }
        for (int i = 1; i < m; i++) add_clause({-vars[i - 1], seq[i][1]});
        for (int i = 2; i < m; i++) {
            for (int j = 1; j <= max_true; j++) add_clause({-seq[i - 1][j], seq[i][j]});
        }
        for (int i = 2; i <= m; i++) add_clause({-vars[i - 1], -seq[i - 1][max_true]});
        for (int i = 2; i < m; i++) {
            for (int j = 2; j <= max_true; j++) {
                add_clause({-vars[i - 1], -seq[i - 1][j - 1], seq[i][j]});
            }
        }
    };

    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                int ab = pair_var[a][b], ac = pair_var[a][c], bc = pair_var[b][c];
                add_clause({-ab, -bc, ac});
                add_clause({-ab, -ac, bc});
                add_clause({-ac, -bc, ab});
            }
        }
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    int conflict_triples = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                if (valid.count(triple_key(a, b, c, n))) continue;
                conflict_triples++;
                add_clause({-pair_var[a][b], -pair_var[a][c], -pair_var[b][c]});
            }
        }
    }

    add_clause({rep_var[1]});
    for (int i = 2; i <= n; i++) {
        vector<int> rep_or_previous{rep_var[i]};
        for (int j = 1; j < i; j++) {
            int same = pair_var[j][i];
            rep_or_previous.push_back(same);
            add_clause({-rep_var[i], -same});
        }
        add_clause(rep_or_previous);
    }
    vector<int> reps;
    for (int i = 1; i <= n; i++) reps.push_back(rep_var[i]);
    add_at_most_k(reps, target_components);

    cerr << "DIRECT_LAZY_SAT_MODEL vars=" << var_count
         << " clauses=" << clause_count
         << " conflict_triples=" << conflict_triples
         << " build_seconds=" << elapsed(start) << "\n";
    unordered_set<long long> full_edge_cuts;

    auto collect_side_for_node = [&](const Tree &tree, int node) {
        vector<int> side;
        function<void(int)> dfs = [&](int u) {
            int label = tree.nodes[u].label;
            if (label >= 1) {
                side.push_back(label);
                return;
            }
            for (int child : tree.nodes[u].child) dfs(child);
        };
        dfs(node);
        sort(side.begin(), side.end());
        return side;
    };

    auto add_full_edge_cut = [&](const Tree &tree, int tree_index, int edge,
            const vector<int> &side_in) {
        vector<int> side = side_in;
        vector<unsigned char> in_side(n + 1, 0);
        for (int x : side) in_side[x] = 1;
        vector<int> outside;
        for (int x = 1; x <= n; x++) if (!in_side[x]) outside.push_back(x);
        if (side.size() > outside.size()) side.swap(outside);
        if (side.size() <= 1 || outside.size() <= 1) return 0;
        int added = 0;
        for (int ai = 0; ai < (int)side.size(); ai++) {
            int a = side[ai];
            for (int ci = ai + 1; ci < (int)side.size(); ci++) {
                int c = side[ci];
                int ac = pair_var[a][c];
                for (int bi = 0; bi < (int)outside.size(); bi++) {
                    int b = outside[bi];
                    int ab = pair_var[a][b];
                    int cb = pair_var[c][b];
                    for (int di = bi + 1; di < (int)outside.size(); di++) {
                        int d = outside[di];
                        int cd = pair_var[c][d];
                        int ad = pair_var[a][d];
                        if (ab >= 0 && cd >= 0 && ac >= 0 &&
                                add_clause({-ab, -cd, ac})) added++;
                        if (ad >= 0 && cb >= 0 && ac >= 0 &&
                                add_clause({-ad, -cb, ac})) added++;
                    }
                }
            }
        }
        long long key = ((long long)tree_index << 32) ^ (unsigned)edge;
        full_edge_cuts.insert(key);
        return added;
    };

    if (preload_side_limit > 1) {
        int preload_edges = 0;
        int preload_clauses = 0;
        for (int tree_index = 0; tree_index < inst.tree_count; tree_index++) {
            const Tree &tree = inst.trees[tree_index];
            for (int edge = 0; edge < (int)tree.nodes.size(); edge++) {
                if (edge == tree.root) continue;
                vector<int> side = collect_side_for_node(tree, edge);
                int small_side = min((int)side.size(), n - (int)side.size());
                if (small_side <= 1 || small_side > preload_side_limit) continue;
                long long key = ((long long)tree_index << 32) ^ (unsigned)edge;
                if (full_edge_cuts.count(key)) continue;
                int added = add_full_edge_cut(tree, tree_index, edge, side);
                if (added > 0) {
                    preload_edges++;
                    preload_clauses += added;
                }
                if (Clock::now() >= deadline) break;
            }
            if (Clock::now() >= deadline) break;
        }
        cerr << "DIRECT_LAZY_PRELOAD_FULL_EDGE side_limit=" << preload_side_limit
             << " edges=" << preload_edges
             << " added=" << preload_clauses
             << " clauses=" << clause_count
             << " seconds=" << elapsed(start) << "\n";
    }

    auto partition_from_model = [&]() {
        vector<int> parent(n + 1);
        for (int i = 1; i <= n; i++) parent[i] = i;
        function<int(int)> find_root = [&](int x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        auto unite = [&](int a, int b) {
            int ra = find_root(a), rb = find_root(b);
            if (ra != rb) parent[rb] = ra;
        };
        for (int a = 1; a <= n; a++) {
            for (int b = a + 1; b <= n; b++) {
                if (solver.val(pair_var[a][b]) > 0) unite(a, b);
            }
        }
        unordered_map<int, vector<int>> groups;
        for (int leaf = 1; leaf <= n; leaf++) groups[find_root(leaf)].push_back(leaf);
        vector<vector<int>> partition;
        for (auto &kv : groups) {
            sort(kv.second.begin(), kv.second.end());
            partition.push_back(kv.second);
        }
        sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        return partition;
    };

    auto add_overlap_cuts = [&](const vector<vector<int>> &partition) {
        struct EdgeViolation {
            const Tree *tree = nullptr;
            int tree_index = -1;
            int edge = -1;
            vector<int> owners;
        };
        vector<EdgeViolation> violations;
        for (int tree_index = 0; tree_index < inst.tree_count; tree_index++) {
            const Tree &tree = inst.trees[tree_index];
            vector<vector<int>> owners(tree.nodes.size());
            for (int cid = 0; cid < (int)partition.size(); cid++) {
                for (int edge : tree_component_edges(tree, partition[cid])) {
                    owners[edge].push_back(cid);
                }
            }
            for (int edge = 0; edge < (int)owners.size(); edge++) {
                vector<int> &list = owners[edge];
                if (list.size() <= 1) continue;
                sort(list.begin(), list.end());
                list.erase(unique(list.begin(), list.end()), list.end());
                violations.push_back(EdgeViolation{&tree, tree_index, edge, list});
            }
        }
        sort(violations.begin(), violations.end(), [](const EdgeViolation &a, const EdgeViolation &b) {
            if (a.owners.size() != b.owners.size()) return a.owners.size() > b.owners.size();
            return a.edge < b.edge;
        });

        auto crossing_pairs_in_block = [&](const Tree &tree, const vector<int> &block, int edge) {
            vector<int> inside, outside;
            for (int leaf : block) {
                if (leaf_below_node(tree, leaf, edge)) inside.push_back(leaf);
                else outside.push_back(leaf);
            }
            vector<pair<int, int>> pairs;
            for (int a : inside) {
                for (int b : outside) {
                    if (pair_var[a][b] >= 0) pairs.push_back({a, b});
                }
            }
            return pairs;
        };
        auto add_owner_pair_cuts = [&]() {
            int owner_pair_added = 0;
            const int max_owner_pair_added = 50000;
            for (const EdgeViolation &violation : violations) {
                int before_violation = owner_pair_added;
                vector<vector<pair<int, int>>> crossing_by_owner;
                crossing_by_owner.reserve(violation.owners.size());
                for (int cid : violation.owners) {
                    crossing_by_owner.push_back(
                        crossing_pairs_in_block(*violation.tree, partition[cid], violation.edge));
                }
                for (int oi = 0; oi < (int)violation.owners.size(); oi++) {
                    if (crossing_by_owner[oi].empty()) continue;
                    int ca = violation.owners[oi];
                    const vector<int> &block_a = partition[ca];
                    for (int oj = oi + 1; oj < (int)violation.owners.size(); oj++) {
                        if (crossing_by_owner[oj].empty()) continue;
                        int cb = violation.owners[oj];
                        const vector<int> &block_b = partition[cb];
                        vector<int> merge_lits;
                        merge_lits.reserve(block_a.size() * block_b.size());
                        for (int x : block_a) {
                            for (int y : block_b) {
                                int same = pair_var[x][y];
                                if (same >= 0) merge_lits.push_back(same);
                            }
                        }
                        if (merge_lits.empty()) continue;
                        for (auto pa : crossing_by_owner[oi]) {
                            int va = pair_var[pa.first][pa.second];
                            for (auto pb : crossing_by_owner[oj]) {
                                int vb = pair_var[pb.first][pb.second];
                                vector<int> clause;
                                clause.reserve(2 + merge_lits.size());
                                clause.push_back(-va);
                                clause.push_back(-vb);
                                clause.insert(clause.end(), merge_lits.begin(), merge_lits.end());
                                if (add_clause(std::move(clause))) owner_pair_added++;
                                if (owner_pair_added >= max_owner_pair_added) {
                                    cerr << "DIRECT_LAZY_OWNER_PAIR_EDGE tree=" << violation.tree_index
                                         << " edge=" << violation.edge
                                         << " owners=" << violation.owners.size()
                                         << " added=" << (owner_pair_added - before_violation)
                                         << " total=" << owner_pair_added << "\n";
                                    return owner_pair_added;
                                }
                            }
                        }
                    }
                }
                if (owner_pair_added > before_violation) {
                    cerr << "DIRECT_LAZY_OWNER_PAIR_EDGE tree=" << violation.tree_index
                         << " edge=" << violation.edge
                         << " owners=" << violation.owners.size()
                         << " added=" << (owner_pair_added - before_violation)
                         << " total=" << owner_pair_added << "\n";
                }
            }
            return owner_pair_added;
        };
        int full_edges_added = 0;
        int full_clause_total = 0;
        const int max_full_edges_per_iter = 6;
        for (const EdgeViolation &violation : violations) {
            long long key = ((long long)violation.tree_index << 32) ^ (unsigned)violation.edge;
            if (full_edge_cuts.count(key)) continue;
            vector<int> side = collect_side_for_node(*violation.tree, violation.edge);
            vector<unsigned char> in_side(n + 1, 0);
            for (int x : side) in_side[x] = 1;
            vector<int> outside;
            for (int x = 1; x <= n; x++) if (!in_side[x]) outside.push_back(x);
            if (side.size() > outside.size()) side.swap(outside);
            long long estimate = (long long)side.size() * (side.size() - 1) / 2 *
                ((long long)outside.size() * (outside.size() - 1) / 2) * 2;
            if (estimate <= 0 || estimate > 150000) continue;
            int added = add_full_edge_cut(*violation.tree, violation.tree_index,
                violation.edge, side);
            if (added > 0) {
                cerr << "DIRECT_LAZY_FULL_EDGE tree=" << violation.tree_index
                     << " edge=" << violation.edge
                     << " side=" << side.size()
                     << " outside=" << outside.size()
                     << " added=" << added << "\n";
                full_edges_added++;
                full_clause_total += added;
                if (full_edges_added >= max_full_edges_per_iter) return full_clause_total;
            }
        }
        if (full_clause_total > 0) return full_clause_total;

        int owner_pair_added = add_owner_pair_cuts();
        if (owner_pair_added > 0) return owner_pair_added;

        int pairwise_added = 0;
        const int max_pairwise_added = 256;
        for (const EdgeViolation &violation : violations) {
            vector<pair<int, int>> useful;
            for (int cid : violation.owners) {
                auto p = crossing_pair_for_edge(*violation.tree, partition[cid], violation.edge);
                if (p.first < 0) continue;
                int cross = pair_var[p.first][p.second];
                if (cross < 0) continue;
                useful.push_back({cid, cross});
            }
            if (useful.size() <= 1) continue;
            int before_violation = pairwise_added;
            for (int i = 0; i < (int)useful.size(); i++) {
                const vector<int> &a = partition[useful[i].first];
                for (int j = i + 1; j < (int)useful.size(); j++) {
                    const vector<int> &b = partition[useful[j].first];
                    vector<int> clause{-useful[i].second, -useful[j].second};
                    for (int x : a) {
                        for (int y : b) {
                            int same = pair_var[x][y];
                            if (same >= 0) clause.push_back(same);
                        }
                    }
                    if (add_clause(std::move(clause))) pairwise_added++;
                    if (pairwise_added >= max_pairwise_added) {
                        cerr << "DIRECT_LAZY_PAIR_EDGE tree=" << violation.tree_index
                             << " edge=" << violation.edge
                             << " owners=" << useful.size()
                             << " added=" << (pairwise_added - before_violation) << "\n";
                        return pairwise_added;
                    }
                }
            }
            if (pairwise_added > before_violation) {
                cerr << "DIRECT_LAZY_PAIR_EDGE tree=" << violation.tree_index
                     << " edge=" << violation.edge
                     << " owners=" << useful.size()
                     << " added=" << (pairwise_added - before_violation) << "\n";
            }
        }
        if (pairwise_added > 0) return pairwise_added;

        int added = 0;
        const int max_added = 96;
        for (const EdgeViolation &violation : violations) {
            vector<int> clause;
            vector<int> useful_owners;
            for (int cid : violation.owners) {
                auto p = crossing_pair_for_edge(*violation.tree, partition[cid], violation.edge);
                if (p.first < 0) continue;
                int cross = pair_var[p.first][p.second];
                if (cross < 0) continue;
                clause.push_back(-cross);
                useful_owners.push_back(cid);
            }
            if (useful_owners.size() <= 1) continue;
            for (int i = 0; i < (int)useful_owners.size(); i++) {
                const vector<int> &a = partition[useful_owners[i]];
                for (int j = i + 1; j < (int)useful_owners.size(); j++) {
                    const vector<int> &b = partition[useful_owners[j]];
                    for (int x : a) {
                        for (int y : b) {
                            int same = pair_var[x][y];
                            if (same >= 0) clause.push_back(same);
                        }
                    }
                }
            }
            if (add_clause(std::move(clause))) added++;
            if (added >= max_added) return added;
        }
        return added;
    };

    for (int iter = 1; elapsed(start) < seconds; iter++) {
        int result = solver.solve();
        if (result == 20) {
            if (proved_unsat) *proved_unsat = true;
            solver.disconnect_terminator();
            cerr << "DIRECT_LAZY_SAT_STATUS result=20 iterations=" << iter
                 << " clauses=" << clause_count
                 << " seconds=" << elapsed(start) << "\n";
            return {};
        }
        if (result != 10) break;
        vector<vector<int>> partition = partition_from_model();
        vector<int> edge_offset(inst.tree_count + 1, 0);
        for (int t = 0; t < inst.tree_count; t++) {
            edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
        }
        if ((int)partition.size() <= target_components &&
                validate_partition(inst, partition, edge_offset)) {
            solver.disconnect_terminator();
            cerr << "DIRECT_LAZY_SAT_STATUS result=10 iterations=" << iter
                 << " components=" << partition.size()
                 << " clauses=" << clause_count
                 << " seconds=" << elapsed(start) << "\n";
            return partition;
        }
        int added = add_overlap_cuts(partition);
        cerr << "DIRECT_LAZY_SAT_ITER iter=" << iter
             << " components=" << partition.size()
             << " added=" << added
             << " clauses=" << clause_count
             << " seconds=" << elapsed(start) << "\n";
        if (added == 0) break;
        if (inst.tree_count >= 20 && clause_count > 5500000) {
            cerr << "DIRECT_LAZY_SAT_ABORT reason=clause_memory_guard"
                 << " target=" << target_components
                 << " clauses=" << clause_count
                 << " t=" << inst.tree_count
                 << " n=" << n << "\n";
            break;
        }
    }
    solver.disconnect_terminator();
    cerr << "DIRECT_LAZY_SAT_STATUS result=unknown target=" << target_components
         << " clauses=" << clause_count
         << " seconds=" << elapsed(start) << "\n";
    return {};
}

static vector<int> ref_cut_path_vars(const Tree &ref, int a, int b,
        const vector<int> &node_to_var, const vector<int> *cut_node,
        bool want_cut_edges, bool negate) {
    vector<int> vars;
    int na = ref.leaf_node[a], nb = ref.leaf_node[b];
    int root = lca_node(ref, na, nb);
    auto collect = [&](int node) {
        while (node >= 0 && node != root) {
            int var = node_to_var[node];
            if (var > 0 && (!cut_node || ((*cut_node)[node] != 0) == want_cut_edges)) {
                vars.push_back(negate ? -var : var);
            }
            node = ref.nodes[node].parent;
        }
    };
    collect(na);
    collect(nb);
    sort(vars.begin(), vars.end(), [](int x, int y) {
        if (abs(x) != abs(y)) return abs(x) < abs(y);
        return x < y;
    });
    vars.erase(unique(vars.begin(), vars.end()), vars.end());
    return vars;
}

static vector<vector<int>> ref_cut_partition_from_model(
        const Tree &ref, const vector<int> &cut_node) {
    vector<vector<int>> partition;
    function<void(int, vector<int>&)> collect = [&](int node, vector<int> &leaves) {
        if (ref.nodes[node].label >= 1) {
            leaves.push_back(ref.nodes[node].label);
            return;
        }
        for (int child : ref.nodes[node].child) {
            if (cut_node[child]) {
                vector<int> child_leaves;
                collect(child, child_leaves);
                if (!child_leaves.empty()) {
                    sort(child_leaves.begin(), child_leaves.end());
                    partition.push_back(std::move(child_leaves));
                }
            } else {
                collect(child, leaves);
            }
        }
    };
    vector<int> root_leaves;
    collect(ref.root, root_leaves);
    if (!root_leaves.empty()) {
        sort(root_leaves.begin(), root_leaves.end());
        partition.push_back(std::move(root_leaves));
    }
    sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return partition;
}

static vector<pair<int, int>> crossing_pairs_for_tree_edge(
        const Tree &tree, const vector<int> &block, int edge) {
    vector<int> inside, outside;
    for (int leaf : block) {
        if (leaf_below_node(tree, leaf, edge)) inside.push_back(leaf);
        else outside.push_back(leaf);
    }
    vector<pair<int, int>> pairs;
    for (int a : inside) {
        for (int b : outside) pairs.push_back({a, b});
    }
    return pairs;
}

static vector<vector<int>> split_partition_to_target(
        vector<vector<int>> partition, int target) {
    while ((int)partition.size() < target) {
        int best = -1;
        for (int i = 0; i < (int)partition.size(); i++) {
            if (partition[i].size() > 1 &&
                    (best < 0 || partition[i].size() > partition[best].size())) {
                best = i;
            }
        }
        if (best < 0) return {};
        vector<int> block = partition[best];
        partition.erase(partition.begin() + best);
        partition.push_back(vector<int>{block[0]});
        block.erase(block.begin());
        partition.push_back(std::move(block));
    }
    sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return partition;
}

static vector<vector<int>> repair_conflict_partition_by_local_cover(
        const Instance &inst, const vector<vector<int>> &partition,
        int target_components, const vector<int> &edge_offset,
        double seconds, int max_component_size = 8,
        int max_component_edges = 160, int branch_width = 240) {
    if (seconds <= 0.02 || (int)partition.size() > target_components) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));

    vector<vector<int>> part_edges(partition.size());
    vector<unsigned char> bad(partition.size(), 0);
    unordered_map<int, int> owner;
    owner.reserve(partition.size() * 8);
    for (int i = 0; i < (int)partition.size(); i++) {
        part_edges[i] = component_edges_all(inst, partition[i], edge_offset);
        for (int edge : part_edges[i]) {
            auto it = owner.find(edge);
            if (it == owner.end()) {
                owner.emplace(edge, i);
            } else if (it->second != i) {
                bad[i] = 1;
                bad[it->second] = 1;
            }
        }
    }

    vector<vector<int>> fixed;
    unordered_set<int> fixed_edges;
    vector<int> region_leaves;
    int bad_components = 0;
    for (int i = 0; i < (int)partition.size(); i++) {
        if (bad[i]) {
            bad_components++;
            region_leaves.insert(
                region_leaves.end(), partition[i].begin(), partition[i].end());
        } else {
            fixed.push_back(partition[i]);
            for (int edge : part_edges[i]) fixed_edges.insert(edge);
        }
    }
    if (bad_components == 0) return {};
    sort(region_leaves.begin(), region_leaves.end());
    region_leaves.erase(unique(region_leaves.begin(), region_leaves.end()),
                        region_leaves.end());
    int max_local_components = target_components - (int)fixed.size();
    if (max_local_components < 1) return {};
    int required_local_savings = (int)region_leaves.size() - max_local_components;
    if (required_local_savings > 14) return {};
    if ((int)region_leaves.size() <= max_local_components) {
        vector<vector<int>> combined = fixed;
        for (int leaf : region_leaves) combined.push_back(vector<int>{leaf});
        sort(combined.begin(), combined.end(), [](const vector<int> &a,
                    const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        return validate_partition(inst, combined, edge_offset) ?
            combined : vector<vector<int>>{};
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    unordered_set<string> seen;
    vector<Component> columns;
    auto add_column = [&](const vector<int> &block) {
        if (block.size() <= 1) return;
        string key = key_of(block);
        if (!seen.insert(key).second) return;
        if (!triple_compatible(block, valid, inst.leaf_count)) return;
        if (!valid_topology(inst, block)) return;
        vector<int> edges = component_edges_all(inst, block, edge_offset);
        if ((int)edges.size() > max_component_edges) return;
        for (int edge : edges) {
            if (fixed_edges.count(edge)) return;
        }
        double score = (double)((int)block.size() - 1) /
            pow(max(1, (int)edges.size()), 0.20);
        columns.push_back(Component{block, std::move(edges), score});
    };

    for (int i = 0; i < (int)partition.size(); i++) {
        if (bad[i] && (int)partition[i].size() <= max_component_size) {
            add_column(partition[i]);
        }
    }

    vector<int> current;
    bool limited = false;
    function<void(int)> dfs = [&](int start_pos) {
        if (Clock::now() >= deadline) {
            limited = true;
            return;
        }
        if ((int)current.size() >= 2) add_column(current);
        if ((int)current.size() >= max_component_size) return;
        for (int pos = start_pos; pos < (int)region_leaves.size(); pos++) {
            current.push_back(region_leaves[pos]);
            if (triple_compatible(current, valid, inst.leaf_count)) {
                dfs(pos + 1);
            }
            current.pop_back();
            if (limited) return;
        }
    };
    dfs(0);
    if (columns.empty()) {
        cerr << "REFERENCE_CUT_REPAIR_NO_COLUMNS bad_components="
             << bad_components << " leaves=" << region_leaves.size()
             << " limited=" << (limited ? 1 : 0) << "\n";
        return {};
    }
    if ((int)columns.size() > 60000) {
        sort(columns.begin(), columns.end(), [](const Component &a,
                    const Component &b) {
            if (a.leaves.size() <= 2 && b.leaves.size() > 2) return true;
            if (b.leaves.size() <= 2 && a.leaves.size() > 2) return false;
            if (a.score != b.score) return a.score > b.score;
            if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
            return a.leaves < b.leaves;
        });
        columns.resize(60000);
    }

    double remain = seconds - elapsed(start);
    if (remain <= 0.02) return {};
    vector<vector<int>> local = local_branch_cover(
        columns, region_leaves, max_local_components, remain, branch_width);
    cerr << "REFERENCE_CUT_REPAIR_ATTEMPT components=" << partition.size()
         << " target=" << target_components
         << " bad_components=" << bad_components
         << " leaves=" << region_leaves.size()
         << " columns=" << columns.size()
         << " max_local=" << max_local_components
         << " local=" << local.size()
         << " limited=" << (limited ? 1 : 0)
         << " seconds=" << elapsed(start) << "\n";
    if (local.empty()) return {};

    vector<vector<int>> combined = fixed;
    combined.insert(combined.end(), local.begin(), local.end());
    sort(combined.begin(), combined.end(), [](const vector<int> &a,
                const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    if ((int)combined.size() <= target_components &&
            validate_partition(inst, combined, edge_offset)) {
        return combined;
    }
    return {};
}

static vector<vector<int>> repair_conflict_partition_by_connected_regions(
        const Instance &inst, vector<vector<int>> partition,
        int target_components, const vector<int> &edge_offset,
        double seconds, int max_component_size = 9,
        int max_component_edges = 220, int branch_width = 360,
        int max_region_leaves = 46, int max_rounds = 24) {
    if (seconds <= 0.05 || (int)partition.size() > target_components) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    unordered_set<uint32_t> valid = collect_valid_triples(inst);

    auto sort_partition = [](vector<vector<int>> &p) {
        sort(p.begin(), p.end(), [](const vector<int> &a, const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
    };

    for (int round = 1; round <= max_rounds && Clock::now() < deadline; round++) {
        int m = (int)partition.size();
        vector<vector<int>> part_edges(m);
        vector<vector<int>> conflict_adj(m);
        vector<int> conflict_degree(m, 0);
        vector<int> parent(m);
        vector<unsigned char> bad(m, 0);
        for (int i = 0; i < m; i++) parent[i] = i;
        function<int(int)> find_root = [&](int x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        auto unite = [&](int a, int b) {
            int ra = find_root(a), rb = find_root(b);
            if (ra != rb) parent[rb] = ra;
        };

        unordered_map<int, vector<int>> owners;
        owners.reserve(m * 8);
        for (int i = 0; i < m; i++) {
            if (!valid_topology(inst, partition[i])) {
                bad[i] = 1;
                conflict_degree[i] += 2;
            }
            part_edges[i] = component_edges_all(inst, partition[i], edge_offset);
            for (int edge : part_edges[i]) owners[edge].push_back(i);
        }
        for (auto &kv : owners) {
            vector<int> &v = kv.second;
            sort(v.begin(), v.end());
            v.erase(unique(v.begin(), v.end()), v.end());
            if (v.size() <= 1) continue;
            for (int idx : v) {
                bad[idx] = 1;
                conflict_degree[idx] += (int)v.size() - 1;
            }
            for (int i = 0; i < (int)v.size(); i++) {
                for (int j = i + 1; j < (int)v.size(); j++) {
                    conflict_adj[v[i]].push_back(v[j]);
                    conflict_adj[v[j]].push_back(v[i]);
                }
            }
            for (int i = 1; i < (int)v.size(); i++) unite(v[0], v[i]);
        }

        bool any_bad = false;
        for (unsigned char x : bad) if (x) { any_bad = true; break; }
        if (!any_bad) {
            sort_partition(partition);
            if ((int)partition.size() <= target_components &&
                    validate_partition(inst, partition, edge_offset)) {
                cerr << "REFERENCE_CUT_CONNECTED_REPAIR_VALID rounds=" << round - 1
                     << " components=" << partition.size()
                     << " seconds=" << elapsed(start) << "\n";
                return partition;
            }
            return {};
        }

        unordered_map<int, vector<int>> grouped;
        for (int i = 0; i < m; i++) {
            if (bad[i]) grouped[find_root(i)].push_back(i);
        }
        struct RegionInfo {
            vector<int> blocks;
            int leaves = 0;
            int required_savings = 0;
            bool window = false;
        };
        vector<RegionInfo> regions;
        set<vector<int>> seen_region_blocks;
        auto make_region_info = [&](vector<int> blocks, bool window) {
            sort(blocks.begin(), blocks.end());
            blocks.erase(unique(blocks.begin(), blocks.end()), blocks.end());
            vector<int> leaves;
            for (int idx : blocks) {
                leaves.insert(leaves.end(), partition[idx].begin(), partition[idx].end());
            }
            sort(leaves.begin(), leaves.end());
            leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
            int max_local = target_components - (m - (int)blocks.size());
            return RegionInfo{
                std::move(blocks),
                (int)leaves.size(),
                (int)leaves.size() - max_local,
                window
            };
        };
        auto add_region = [&](const vector<int> &blocks, bool window) {
            if (blocks.empty()) return;
            vector<int> key = blocks;
            sort(key.begin(), key.end());
            key.erase(unique(key.begin(), key.end()), key.end());
            if (!seen_region_blocks.insert(key).second) return;
            regions.push_back(make_region_info(std::move(key), window));
        };
        auto block_union_leaf_count = [&](const vector<int> &blocks) {
            vector<int> leaves;
            for (int idx : blocks) {
                leaves.insert(leaves.end(), partition[idx].begin(), partition[idx].end());
            }
            sort(leaves.begin(), leaves.end());
            leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
            return (int)leaves.size();
        };
        auto add_conflict_windows = [&](const vector<int> &component_blocks) {
            vector<int> anchors = component_blocks;
            sort(anchors.begin(), anchors.end(), [&](int a, int b) {
                if (conflict_degree[a] != conflict_degree[b]) return conflict_degree[a] > conflict_degree[b];
                if (partition[a].size() != partition[b].size()) return partition[a].size() < partition[b].size();
                return partition[a] < partition[b];
            });
            if ((int)anchors.size() > 28) anchors.resize(28);
            vector<unsigned char> in_component(m, 0);
            for (int idx : component_blocks) in_component[idx] = 1;
            for (int anchor : anchors) {
                vector<int> selected{anchor};
                vector<unsigned char> in_selected(m, 0);
                in_selected[anchor] = 1;
                int selected_leaves = block_union_leaf_count(selected);
                bool grew = true;
                while (grew) {
                    grew = false;
                    vector<int> candidates;
                    for (int idx : selected) {
                        for (int nb : conflict_adj[idx]) {
                            if (in_component[nb] && !in_selected[nb]) candidates.push_back(nb);
                        }
                    }
                    sort(candidates.begin(), candidates.end());
                    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
                    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
                        if (conflict_degree[a] != conflict_degree[b]) return conflict_degree[a] > conflict_degree[b];
                        if (partition[a].size() != partition[b].size()) return partition[a].size() < partition[b].size();
                        return partition[a] < partition[b];
                    });
                    for (int cand : candidates) {
                        vector<int> trial = selected;
                        trial.push_back(cand);
                        int trial_leaves = block_union_leaf_count(trial);
                        if (trial_leaves > max_region_leaves) continue;
                        selected.push_back(cand);
                        in_selected[cand] = 1;
                        selected_leaves = trial_leaves;
                        grew = true;
                        break;
                    }
                }
                if (selected_leaves <= max_region_leaves) add_region(selected, true);
                for (int nb : conflict_adj[anchor]) {
                    if (!in_component[nb]) continue;
                    vector<int> pair_window{anchor, nb};
                    if (block_union_leaf_count(pair_window) <= max_region_leaves) {
                        add_region(pair_window, true);
                    }
                }
            }
        };
        for (auto &kv : grouped) {
            add_region(kv.second, false);
            if (block_union_leaf_count(kv.second) > max_region_leaves) {
                add_conflict_windows(kv.second);
            }
        }
        sort(regions.begin(), regions.end(), [](const RegionInfo &a, const RegionInfo &b) {
            bool afit = a.leaves <= 46, bfit = b.leaves <= 46;
            if (afit != bfit) return afit > bfit;
            if (a.window != b.window) return a.window > b.window;
            if (a.required_savings != b.required_savings) return a.required_savings < b.required_savings;
            if (a.leaves != b.leaves) return a.leaves < b.leaves;
            return a.blocks.size() < b.blocks.size();
        });

        bool changed = false;
        for (const RegionInfo &region_info : regions) {
            if (Clock::now() >= deadline) break;
            vector<unsigned char> in_region(m, 0);
            for (int idx : region_info.blocks) in_region[idx] = 1;
            vector<vector<int>> fixed;
            unordered_set<int> fixed_edges;
            vector<int> region_leaves;
            fixed.reserve(m - region_info.blocks.size());
            for (int i = 0; i < m; i++) {
                if (in_region[i]) {
                    region_leaves.insert(region_leaves.end(),
                        partition[i].begin(), partition[i].end());
                } else {
                    fixed.push_back(partition[i]);
                    for (int edge : part_edges[i]) fixed_edges.insert(edge);
                }
            }
            sort(region_leaves.begin(), region_leaves.end());
            region_leaves.erase(unique(region_leaves.begin(), region_leaves.end()),
                                region_leaves.end());
            int max_local_components = target_components - (int)fixed.size();
            if (max_local_components < 1) continue;
            if ((int)region_leaves.size() > max_region_leaves) {
                cerr << "REFERENCE_CUT_CONNECTED_REPAIR_SKIP round=" << round
                     << " blocks=" << region_info.blocks.size()
                     << " leaves=" << region_leaves.size()
                     << " max_local=" << max_local_components
                     << " window=" << (region_info.window ? 1 : 0)
                     << " reason=leaf_limit\n";
                continue;
            }

            unordered_set<string> seen;
            vector<Component> columns;
            auto add_column = [&](vector<int> block) {
                sort(block.begin(), block.end());
                block.erase(unique(block.begin(), block.end()), block.end());
                if (block.size() <= 1) return;
                string key = key_of(block);
                if (!seen.insert(key).second) return;
                if (!triple_compatible(block, valid, inst.leaf_count)) return;
                if (!valid_topology(inst, block)) return;
                vector<int> edges = component_edges_all(inst, block, edge_offset);
                if ((int)edges.size() > max_component_edges) return;
                for (int edge : edges) if (fixed_edges.count(edge)) return;
                double score = (double)((int)block.size() - 1) /
                    pow(max(1, (int)edges.size()), 0.20);
                columns.push_back(Component{std::move(block), std::move(edges), score});
            };
            for (int idx : region_info.blocks) {
                if ((int)partition[idx].size() <= max_component_size) {
                    add_column(partition[idx]);
                }
            }
            vector<int> current;
            bool limited = false;
            function<void(int)> dfs = [&](int start_pos) {
                if (Clock::now() >= deadline || (int)columns.size() > 90000) {
                    limited = true;
                    return;
                }
                if ((int)current.size() >= 2) add_column(current);
                if ((int)current.size() >= max_component_size) return;
                for (int pos = start_pos; pos < (int)region_leaves.size(); pos++) {
                    current.push_back(region_leaves[pos]);
                    if (triple_compatible(current, valid, inst.leaf_count)) {
                        dfs(pos + 1);
                    }
                    current.pop_back();
                    if (limited) return;
                }
            };
            dfs(0);
            if (columns.empty()) continue;
            unordered_map<int, int> leaf_to_local;
            for (int i = 0; i < (int)region_leaves.size(); i++) leaf_to_local[region_leaves[i]] = i;
            columns = reduce_repair_columns_preserve_pairs(
                columns, leaf_to_local, (int)region_leaves.size(), 180);
            double remain = chrono::duration<double>(deadline - Clock::now()).count();
            if (remain <= 0.02) break;
            double attempt_seconds = min(12.0, remain);
            vector<vector<int>> local = local_branch_cover(
                columns, region_leaves, max_local_components,
                attempt_seconds, branch_width);
            cerr << "REFERENCE_CUT_CONNECTED_REPAIR_ATTEMPT round=" << round
                 << " blocks=" << region_info.blocks.size()
                 << " leaves=" << region_leaves.size()
                 << " columns=" << columns.size()
                 << " max_local=" << max_local_components
                 << " local=" << local.size()
                 << " limited=" << (limited ? 1 : 0)
                 << " window=" << (region_info.window ? 1 : 0)
                 << " seconds=" << elapsed(start) << "\n";
            if (local.empty()) continue;

            vector<vector<int>> candidate = fixed;
            candidate.insert(candidate.end(), local.begin(), local.end());
            sort_partition(candidate);
            if ((int)candidate.size() > target_components) continue;
            partition = std::move(candidate);
            changed = true;
            break;
        }
        if (!changed) {
            cerr << "REFERENCE_CUT_CONNECTED_REPAIR_STALLED round=" << round
                 << " regions=" << regions.size()
                 << " components=" << partition.size()
                 << " seconds=" << elapsed(start) << "\n";
            return {};
        }
    }

    sort_partition(partition);
    if ((int)partition.size() <= target_components &&
            validate_partition(inst, partition, edge_offset)) {
        cerr << "REFERENCE_CUT_CONNECTED_REPAIR_VALID rounds=" << max_rounds
             << " components=" << partition.size()
             << " seconds=" << elapsed(start) << "\n";
        return partition;
    }
    return {};
}

static vector<vector<int>> reference_cut_cadical_target(const Instance &inst,
        int target_components, double seconds, bool *proved_unsat = nullptr,
        Pool *harvest_pool = nullptr, int harvest_max_component_size = 8,
        int harvest_max_added = 50000) {
    if (proved_unsat) *proved_unsat = false;
    if (seconds <= 0.0 || target_components < 1 || inst.trees.empty()) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    const Tree &ref = inst.trees[0];
    vector<int> cuttable;
    vector<int> node_to_var(ref.nodes.size(), 0);
    for (int node = 0; node < (int)ref.nodes.size(); node++) {
        if (node == ref.root) continue;
        node_to_var[node] = (int)cuttable.size() + 1;
        cuttable.push_back(node);
    }
    int var_count = (int)cuttable.size();
    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    TimeTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    set<vector<int>> seen;
    long long clause_count = 0;
    auto normalize = [](vector<int> clause) {
        sort(clause.begin(), clause.end(), [](int x, int y) {
            if (abs(x) != abs(y)) return abs(x) < abs(y);
            return x < y;
        });
        clause.erase(unique(clause.begin(), clause.end()), clause.end());
        for (int i = 0; i < (int)clause.size(); i++) {
            for (int j = i + 1; j < (int)clause.size(); j++) {
                if (clause[i] == -clause[j]) return vector<int>{};
            }
        }
        return clause;
    };
    auto add_clause = [&](vector<int> clause) -> bool {
        clause = normalize(std::move(clause));
        if (clause.empty()) return false;
        if (!seen.insert(clause).second) return false;
        for (int lit : clause) solver.add(lit);
        solver.add(0);
        clause_count++;
        return true;
    };
    auto add_raw_clause = [&](const vector<int> &clause) {
        for (int lit : clause) solver.add(lit);
        solver.add(0);
        clause_count++;
    };
    auto add_at_most_k = [&](int variable_count, int max_true) {
        if (max_true < 0) {
            add_raw_clause({});
            return;
        }
        if (max_true >= variable_count) return;
        if (max_true == 0) {
            for (int variable = 1; variable <= variable_count; variable++) {
                add_raw_clause({-variable});
            }
            return;
        }
        vector<vector<int>> seq(variable_count + 1, vector<int>(max_true + 1, 0));
        for (int i = 1; i <= variable_count - 1; i++) {
            for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
        }
        for (int i = 1; i <= variable_count - 1; i++) add_raw_clause({-i, seq[i][1]});
        for (int i = 2; i <= variable_count - 1; i++) {
            for (int j = 1; j <= max_true; j++) add_raw_clause({-seq[i - 1][j], seq[i][j]});
        }
        for (int i = 2; i <= variable_count; i++) add_raw_clause({-i, -seq[i - 1][max_true]});
        for (int i = 2; i <= variable_count - 1; i++) {
            for (int j = 2; j <= max_true; j++) {
                add_raw_clause({-i, -seq[i - 1][j - 1], seq[i][j]});
            }
        }
    };

    int conflict_triples = 0;
    for (int a = 1; a <= inst.leaf_count; a++) {
        for (int b = a + 1; b <= inst.leaf_count; b++) {
            for (int c = b + 1; c <= inst.leaf_count; c++) {
                int sig = triple_sig(inst.trees[0], a, b, c);
                bool ok = true;
                for (int t = 1; t < inst.tree_count; t++) {
                    if (triple_sig(inst.trees[t], a, b, c) != sig) {
                        ok = false;
                        break;
                    }
                }
                if (ok) continue;
                vector<int> clause;
                vector<int> ab = ref_cut_path_vars(ref, a, b, node_to_var, nullptr, false, false);
                vector<int> ac = ref_cut_path_vars(ref, a, c, node_to_var, nullptr, false, false);
                vector<int> bc = ref_cut_path_vars(ref, b, c, node_to_var, nullptr, false, false);
                clause.insert(clause.end(), ab.begin(), ab.end());
                clause.insert(clause.end(), ac.begin(), ac.end());
                clause.insert(clause.end(), bc.begin(), bc.end());
                if (add_clause(std::move(clause))) conflict_triples++;
            }
        }
    }
    add_at_most_k((int)cuttable.size(), target_components - 1);
    for (int variable = 1; variable <= (int)cuttable.size(); variable++) solver.phase(-variable);
    cerr << "REFERENCE_CUT_MODEL vars=" << var_count
         << " cut_edges=" << cuttable.size()
         << " target=" << target_components
         << " conflict_triples=" << conflict_triples
         << " clauses=" << clause_count
         << " build_seconds=" << elapsed(start) << "\n";

    vector<int> edge_offset(inst.tree_count + 1, 0);
    for (int t = 0; t < inst.tree_count; t++) {
        edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
    }
    int repair_attempts = 0;
    int harvest_added = 0;
    unordered_set<uint32_t> harvest_valid;
    if (harvest_pool) harvest_valid = collect_valid_triples(inst);
    mt19937_64 harvest_rng(1469598103934665603ULL ^
        ((uint64_t)target_components * 1099511628211ULL) ^
        ((uint64_t)inst.leaf_count << 17));
    for (int iter = 1; elapsed(start) < seconds; iter++) {
        int result = solver.solve();
        if (result == 20) {
            if (proved_unsat) *proved_unsat = true;
            solver.disconnect_terminator();
            cerr << "REFERENCE_CUT_STATUS result=20 iterations=" << iter - 1
                 << " clauses=" << clause_count
                 << " seconds=" << elapsed(start) << "\n";
            return {};
        }
        if (result != 10) break;
        vector<int> cut_node(ref.nodes.size(), 0);
        vector<int> model(cuttable.size() + 1, 0);
        for (int variable = 1; variable <= (int)cuttable.size(); variable++) {
            model[variable] = solver.val(variable) > 0 ? 1 : 0;
            if (model[variable]) cut_node[cuttable[variable - 1]] = 1;
        }
        vector<vector<int>> partition = ref_cut_partition_from_model(ref, cut_node);
        if (harvest_pool && harvest_added < harvest_max_added) {
            int before = (int)harvest_pool->comps.size();
            for (const vector<int> &block : partition) {
                if (harvest_added >= harvest_max_added) break;
                if ((int)block.size() < 2 ||
                        (int)block.size() > max(inst.leaf_count, harvest_max_component_size)) {
                    continue;
                }
                if ((int)block.size() <= harvest_max_component_size) {
                    if (add_component(*harvest_pool, inst, harvest_valid,
                                edge_offset, block)) {
                        harvest_added++;
                    }
                    continue;
                }
                int trials = min(80, max(12, (int)block.size() / 2));
                for (int trial = 0; trial < trials && harvest_added < harvest_max_added; trial++) {
                    int lo = min(4, harvest_max_component_size);
                    if ((int)block.size() < lo) break;
                    int hi = min(harvest_max_component_size, (int)block.size());
                    uniform_int_distribution<int> size_dist(lo, hi);
                    vector<int> sample = random_sample(block, size_dist(harvest_rng), harvest_rng);
                    if (add_component(*harvest_pool, inst, harvest_valid,
                                edge_offset, std::move(sample))) {
                        harvest_added++;
                    }
                }
            }
            int added_now = (int)harvest_pool->comps.size() - before;
            if (added_now > 0 && (iter <= 10 || iter % 25 == 0)) {
                cerr << "REFERENCE_CUT_HARVEST iter=" << iter
                     << " added=" << added_now
                     << " total_added=" << harvest_added
                     << " pool=" << harvest_pool->comps.size()
                     << " components=" << partition.size() << "\n";
            }
        }
        if ((int)partition.size() <= target_components &&
                validate_partition(inst, partition, edge_offset)) {
            vector<vector<int>> target_partition = split_partition_to_target(
                std::move(partition), target_components);
            if (!target_partition.empty() &&
                    validate_partition(inst, target_partition, edge_offset)) {
                solver.disconnect_terminator();
                cerr << "REFERENCE_CUT_STATUS result=10 iterations=" << iter
                     << " components=" << target_partition.size()
                     << " clauses=" << clause_count
                     << " seconds=" << elapsed(start) << "\n";
                return target_partition;
            }
        }
        if ((int)partition.size() <= target_components && repair_attempts < 8) {
            double repair_seconds = min(18.0, max(0.0, seconds - elapsed(start)));
            if (repair_seconds > 0.05) {
                repair_attempts++;
                vector<vector<int>> connected_repaired =
                    repair_conflict_partition_by_connected_regions(
                        inst, partition, target_components, edge_offset,
                        repair_seconds);
                if (!connected_repaired.empty() &&
                        validate_partition(inst, connected_repaired, edge_offset) &&
                        (int)connected_repaired.size() <= target_components) {
                    solver.disconnect_terminator();
                    cerr << "REFERENCE_CUT_CONNECTED_REPAIR_SUCCESS iter=" << iter
                         << " components=" << connected_repaired.size()
                         << " attempts=" << repair_attempts
                         << " seconds=" << elapsed(start) << "\n";
                    return connected_repaired;
                }
                repair_seconds = min(8.0, max(0.0, seconds - elapsed(start)));
                vector<vector<int>> repaired =
                    repair_conflict_partition_by_local_cover(
                        inst, partition, target_components, edge_offset,
                        repair_seconds);
                if (!repaired.empty() &&
                        validate_partition(inst, repaired, edge_offset) &&
                        (int)repaired.size() <= target_components) {
                    solver.disconnect_terminator();
                    cerr << "REFERENCE_CUT_REPAIR_SUCCESS iter=" << iter
                         << " components=" << repaired.size()
                         << " attempts=" << repair_attempts
                         << " seconds=" << elapsed(start) << "\n";
                    return repaired;
                }
            }
        }

        struct RefEdgeViolation {
            const Tree *tree = nullptr;
            int tree_index = -1;
            int edge = -1;
            vector<int> owners;
            long long pair_work = 0;
        };
        vector<RefEdgeViolation> violations;
        for (int tree_index = 0; tree_index < inst.tree_count; tree_index++) {
            const Tree &tree = inst.trees[tree_index];
            vector<vector<int>> owners(tree.nodes.size());
            for (int cid = 0; cid < (int)partition.size(); cid++) {
                for (int edge : tree_component_edges(tree, partition[cid])) {
                    if (find(owners[edge].begin(), owners[edge].end(), cid) == owners[edge].end()) {
                        owners[edge].push_back(cid);
                    }
                }
            }
            for (int edge = 0; edge < (int)owners.size(); edge++) {
                if (owners[edge].size() <= 1) continue;
                long long work = 0;
                for (int i = 0; i < (int)owners[edge].size(); i++) {
                    vector<pair<int, int>> pi =
                        crossing_pairs_for_tree_edge(tree, partition[owners[edge][i]], edge);
                    for (int j = i + 1; j < (int)owners[edge].size(); j++) {
                        vector<pair<int, int>> pj =
                            crossing_pairs_for_tree_edge(tree, partition[owners[edge][j]], edge);
                        work += (long long)pi.size() * (long long)pj.size();
                    }
                }
                violations.push_back(RefEdgeViolation{
                    &tree, tree_index, edge, owners[edge], work});
            }
        }
        sort(violations.begin(), violations.end(),
            [](const RefEdgeViolation &a, const RefEdgeViolation &b) {
                if (a.owners.size() != b.owners.size()) return a.owners.size() > b.owners.size();
                if (a.pair_work != b.pair_work) return a.pair_work > b.pair_work;
                if (a.tree_index != b.tree_index) return a.tree_index < b.tree_index;
                return a.edge < b.edge;
        });
        int added = 0;
        const int max_added = 250000;
        for (const RefEdgeViolation &violation : violations) {
            if (added >= max_added) break;
            const Tree &tree = *violation.tree;
            vector<vector<pair<int, int>>> crossing_by_owner;
            crossing_by_owner.reserve(violation.owners.size());
            for (int cid : violation.owners) {
                crossing_by_owner.push_back(
                    crossing_pairs_for_tree_edge(tree, partition[cid], violation.edge));
            }
            for (int i = 0; i < (int)violation.owners.size() && added < max_added; i++) {
                int ca = violation.owners[i];
                vector<pair<int, int>> &pairs_a = crossing_by_owner[i];
                if (pairs_a.empty()) continue;
                for (int j = i + 1; j < (int)violation.owners.size() && added < max_added; j++) {
                    int cb = violation.owners[j];
                    vector<pair<int, int>> &pairs_b = crossing_by_owner[j];
                    if (pairs_b.empty()) continue;
                    vector<int> rep_merge_path = ref_cut_path_vars(
                        ref, partition[ca][0], partition[cb][0],
                        node_to_var, nullptr, false, false);
                    for (auto pa : pairs_a) {
                        vector<int> split_a = ref_cut_path_vars(
                            ref, pa.first, pa.second, node_to_var, &cut_node, false, false);
                        if (split_a.empty()) continue;
                        for (auto pb : pairs_b) {
                            vector<int> split_b = ref_cut_path_vars(
                                ref, pb.first, pb.second, node_to_var, &cut_node, false, false);
                            if (split_b.empty()) continue;
                            const vector<int> &merge_path = rep_merge_path;
                            if (merge_path.empty()) {
                                vector<int> clause;
                                clause.insert(clause.end(), split_a.begin(), split_a.end());
                                clause.insert(clause.end(), split_b.begin(), split_b.end());
                                if (add_clause(std::move(clause))) added++;
                                if (added >= max_added) break;
                                continue;
                            }
                            for (int merge_var : merge_path) {
                                vector<int> clause;
                                clause.insert(clause.end(), split_a.begin(), split_a.end());
                                clause.insert(clause.end(), split_b.begin(), split_b.end());
                                clause.push_back(-merge_var);
                                if (add_clause(std::move(clause))) added++;
                                if (added >= max_added) break;
                            }
                        }
                    }
                }
            }
        }
        if (added == 0) {
            vector<int> clause;
            clause.reserve(cuttable.size());
            for (int variable = 1; variable <= (int)cuttable.size(); variable++) {
                clause.push_back(model[variable] ? -variable : variable);
            }
            if (add_clause(std::move(clause))) added++;
        }
        if (iter <= 10 || iter % 50 == 0 || added >= 50000) {
            cerr << "REFERENCE_CUT_ITER iter=" << iter
                 << " components=" << partition.size()
                 << " added=" << added
                 << " clauses=" << clause_count
                 << " seconds=" << elapsed(start) << "\n";
        }
        if (added == 0) break;
    }
    solver.disconnect_terminator();
    cerr << "REFERENCE_CUT_STATUS result=unknown target=" << target_components
         << " clauses=" << clause_count
         << " harvest_added=" << harvest_added
         << " seconds=" << elapsed(start) << "\n";
    return {};
}
#else
static vector<vector<int>> direct_pair_cadical_target(const Instance &, int, double, bool *proved_unsat = nullptr) {
    if (proved_unsat) *proved_unsat = false;
    return {};
}
static bool direct_pair_cadical_rule_out(const Instance &, int, double) {
    return false;
}
static vector<vector<int>> direct_pair_cadical_lazy_target(const Instance &, int, double,
        bool *proved_unsat = nullptr, int preload_side_limit = 0) {
    (void)preload_side_limit;
    if (proved_unsat) *proved_unsat = false;
    return {};
}
static vector<vector<int>> edge_owner_cadical_target(const Instance &, int, double,
        bool *proved_unsat = nullptr) {
    if (proved_unsat) *proved_unsat = false;
    return {};
}
static vector<vector<int>> reference_cut_cadical_target(const Instance &, int, double,
        bool *proved_unsat = nullptr, Pool *harvest_pool = nullptr,
        int harvest_max_component_size = 8, int harvest_max_added = 50000) {
    (void)harvest_pool;
    (void)harvest_max_component_size;
    (void)harvest_max_added;
    if (proved_unsat) *proved_unsat = false;
    return {};
}
#endif

#ifdef USE_HIGHS
static vector<vector<int>> direct_pair_highs_solve(const Instance &inst, int target_components,
        bool has_target, double seconds, const vector<int> &edge_offset, bool require_optimal,
        bool *proved_infeasible = nullptr, bool *proved_optimal = nullptr) {
    if (seconds <= 0.0) return {};
    if (proved_infeasible) *proved_infeasible = false;
    if (proved_optimal) *proved_optimal = false;
    auto start = Clock::now();
    int n = inst.leaf_count;
    vector<vector<int>> pair_var(n + 1, vector<int>(n + 1, -1));
    int var_count = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            pair_var[a][b] = pair_var[b][a] = var_count++;
        }
    }
    vector<int> rep_var(n + 1, -1);
    for (int i = 1; i <= n; i++) rep_var[i] = var_count++;

    vector<vector<pair<int, double>>> col_entries(var_count);
    vector<double> row_lower;
    vector<double> row_upper;
    unordered_set<ClauseKey, ClauseKeyHash> seen_clauses;

    auto add_linear = [&](const vector<pair<int, double>> &terms, double lower, double upper) {
        int row = (int)row_lower.size();
        row_lower.push_back(lower);
        row_upper.push_back(upper);
        for (auto [var, coeff] : terms) col_entries[var].push_back({row, coeff});
    };

    auto add_clause3 = [&](int l1, int l2, int l3) {
        array<int, 3> lits{l1, l2, l3};
        sort(lits.begin(), lits.end(), [](int x, int y) {
            if (abs(x) != abs(y)) return abs(x) < abs(y);
            return x < y;
        });
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                if (lits[i] == -lits[j]) return;
            }
        }
        ClauseKey key{lits[0], lits[1], lits[2]};
        if (!seen_clauses.insert(key).second) return;
        vector<pair<int, double>> terms;
        double lower = 1.0;
        for (int lit : lits) {
            int var = abs(lit) - 1;
            if (lit > 0) {
                terms.push_back({var, 1.0});
            } else {
                terms.push_back({var, -1.0});
                lower -= 1.0;
            }
        }
        add_linear(terms, lower, kHighsInf);
    };

    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                int ab = pair_var[a][b] + 1;
                int ac = pair_var[a][c] + 1;
                int bc = pair_var[b][c] + 1;
                add_clause3(-ab, -bc, ac);
                add_clause3(-ab, -ac, bc);
                add_clause3(-ac, -bc, ab);
            }
        }
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    int conflict_triples = 0;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            for (int c = b + 1; c <= n; c++) {
                if (valid.count(triple_key(a, b, c, n))) continue;
                conflict_triples++;
                add_clause3(-(pair_var[a][b] + 1), -(pair_var[a][c] + 1), -(pair_var[b][c] + 1));
            }
        }
    }

    add_linear({{rep_var[1], 1.0}}, 1.0, 1.0);
    for (int i = 2; i <= n; i++) {
        vector<pair<int, double>> rep_or_previous{{rep_var[i], 1.0}};
        for (int j = 1; j < i; j++) {
            int same = pair_var[j][i];
            add_linear({{rep_var[i], 1.0}, {same, 1.0}}, -kHighsInf, 1.0);
            rep_or_previous.push_back({same, 1.0});
        }
        add_linear(rep_or_previous, 1.0, kHighsInf);
    }
    vector<pair<int, double>> rep_sum;
    for (int i = 1; i <= n; i++) rep_sum.push_back({rep_var[i], 1.0});
    if (has_target) add_linear(rep_sum, -kHighsInf, (double)target_components);

    int edge_clauses = 0;
    vector<int> all_leaves(n);
    for (int i = 0; i < n; i++) all_leaves[i] = i + 1;
    (void)edge_offset;
    for (const Tree &tree : inst.trees) {
        for (int node = 0; node < (int)tree.nodes.size(); node++) {
            if (node == tree.root) continue;
            vector<int> side;
            function<void(int)> collect = [&](int u) {
                int label = tree.nodes[u].label;
                if (label >= 1) {
                    side.push_back(label);
                    return;
                }
                for (int child : tree.nodes[u].child) collect(child);
            };
            collect(node);
            sort(side.begin(), side.end());
            if (side.empty() || (int)side.size() == n) continue;
            vector<int> outside;
            vector<unsigned char> in_side(n + 1, 0);
            for (int x : side) in_side[x] = 1;
            for (int x : all_leaves) if (!in_side[x]) outside.push_back(x);
            if (side.size() > outside.size()) side.swap(outside);
            for (int ai = 0; ai < (int)side.size(); ai++) {
                int a = side[ai];
                for (int ci = ai + 1; ci < (int)side.size(); ci++) {
                    int c = side[ci];
                    int ac = pair_var[a][c] + 1;
                    for (int bi = 0; bi < (int)outside.size(); bi++) {
                        int b = outside[bi];
                        int ab = pair_var[a][b] + 1;
                        int cb = pair_var[c][b] + 1;
                        for (int di = bi + 1; di < (int)outside.size(); di++) {
                            int d = outside[di];
                            int cd = pair_var[c][d] + 1;
                            int ad = pair_var[a][d] + 1;
                            size_t before = seen_clauses.size();
                            add_clause3(-ab, -cd, ac);
                            if (seen_clauses.size() != before) edge_clauses++;
                            before = seen_clauses.size();
                            add_clause3(-ad, -cb, ac);
                            if (seen_clauses.size() != before) edge_clauses++;
                        }
                    }
                }
            }
        }
    }

    vector<int> starts;
    vector<int> indices;
    vector<double> values;
    starts.reserve(var_count + 1);
    starts.push_back(0);
    for (int var = 0; var < var_count; var++) {
        for (auto [row, value] : col_entries[var]) {
            indices.push_back(row);
            values.push_back(value);
        }
        starts.push_back((int)indices.size());
    }

    HighsModel model;
    model.lp_.num_col_ = var_count;
    model.lp_.num_row_ = (int)row_lower.size();
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.col_cost_.assign(var_count, 0.0);
    for (int i = 1; i <= n; i++) model.lp_.col_cost_[rep_var[i]] = 1.0;
    model.lp_.col_lower_.assign(var_count, 0.0);
    model.lp_.col_upper_.assign(var_count, 1.0);
    model.lp_.row_lower_ = std::move(row_lower);
    model.lp_.row_upper_ = std::move(row_upper);
    model.lp_.integrality_.assign(var_count, HighsVarType::kInteger);
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = std::move(starts);
    model.lp_.a_matrix_.index_ = std::move(indices);
    model.lp_.a_matrix_.value_ = std::move(values);

    cerr << "DIRECT_PAIR_MODEL vars=" << var_count
         << " rows=" << model.lp_.num_row_
         << " conflict_triples=" << conflict_triples
         << " edge_clauses=" << edge_clauses
         << " build_seconds=" << elapsed(start) << "\n";

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("time_limit", max(0.25, seconds - elapsed(start)));
    highs.setOptionValue("mip_rel_gap", 0.0);
    HighsStatus pass_status = highs.passModel(model);
    if (pass_status != HighsStatus::kOk) {
        cerr << "DIRECT_PAIR_PASS_STATUS status=" << (int)pass_status << "\n";
        return {};
    }
    HighsStatus run_status = highs.run();
    HighsModelStatus status = highs.getModelStatus();
    const HighsSolution &sol = highs.getSolution();
    cerr << "DIRECT_PAIR_STATUS run_status=" << (int)run_status
         << " model_status=" << (int)status
         << " seconds=" << elapsed(start) << "\n";
    if (proved_infeasible && status == HighsModelStatus::kInfeasible) *proved_infeasible = true;
    if (proved_optimal && status == HighsModelStatus::kOptimal) *proved_optimal = true;
    if (require_optimal && status != HighsModelStatus::kOptimal) return {};
    if (run_status != HighsStatus::kOk && (int)sol.col_value.size() != var_count) return {};
    if ((int)sol.col_value.size() != var_count) return {};

    vector<int> parent(n + 1);
    for (int i = 1; i <= n; i++) parent[i] = i;
    function<int(int)> find_root = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find_root(a), rb = find_root(b);
        if (ra != rb) parent[rb] = ra;
    };
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            if (sol.col_value[pair_var[a][b]] > 0.5) unite(a, b);
        }
    }
    unordered_map<int, vector<int>> groups;
    for (int leaf = 1; leaf <= n; leaf++) groups[find_root(leaf)].push_back(leaf);
    vector<vector<int>> partition;
    for (auto &kv : groups) {
        sort(kv.second.begin(), kv.second.end());
        partition.push_back(kv.second);
    }
    sort(partition.begin(), partition.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return partition;
}

static vector<vector<int>> direct_pair_highs_target(const Instance &inst, int target_components,
        double seconds, const vector<int> &edge_offset) {
    return direct_pair_highs_solve(inst, target_components, true, seconds, edge_offset, false);
}

static vector<vector<int>> direct_pair_highs_optimal(const Instance &inst,
        double seconds, const vector<int> &edge_offset) {
    return direct_pair_highs_solve(inst, -1, false, seconds, edge_offset, true);
}

static bool direct_pair_highs_rule_out(const Instance &inst, int target_components,
        double seconds, const vector<int> &edge_offset) {
    if (target_components < 1) return true;
    bool proved_infeasible = false;
    (void)direct_pair_highs_solve(inst, target_components, true, seconds, edge_offset,
        false, &proved_infeasible, nullptr);
    return proved_infeasible;
}
#else
static vector<vector<int>> direct_pair_highs_target(const Instance &, int, double, const vector<int> &) {
    return {};
}
static vector<vector<int>> direct_pair_highs_optimal(const Instance &, double, const vector<int> &) {
    return {};
}
static bool direct_pair_highs_rule_out(const Instance &, int, double, const vector<int> &) {
    return false;
}
#endif

static bool contains_index(const vector<int> &items, int value) {
    return find(items.begin(), items.end(), value) != items.end();
}

static vector<int> guided_local_indices(const Instance &inst, const vector<vector<int>> &base,
        const vector<vector<int>> &base_edges, const vector<int> &edge_offset, int component_limit) {
    if (component_limit <= 0 || base.empty()) return {};
    int anchor = 0;
    for (int i = 1; i < (int)base.size(); i++) {
        if (base[i].size() > base[anchor].size()) anchor = i;
    }
    unordered_set<int> selected;
    selected.insert(anchor);
    for (int i = 0; i < (int)base.size(); i++) {
        if (i != anchor && base[i].size() >= 2 && base[i].size() <= 5) selected.insert(i);
    }

    vector<int> singles;
    for (int i = 0; i < (int)base.size(); i++) if (base[i].size() == 1) singles.push_back(i);

    struct PairInfo { int conflicts_size; vector<int> conflicts; int left; int right; vector<int> pair; };
    vector<PairInfo> pair_infos;
    unordered_map<int, int> conflict_score;
    for (int a = 0; a < (int)singles.size(); a++) {
        int i = singles[a];
        for (int b = a + 1; b < (int)singles.size(); b++) {
            int j = singles[b];
            vector<int> pair{base[i][0], base[j][0]};
            sort(pair.begin(), pair.end());
            if (!valid_topology(inst, pair)) continue;
            vector<int> pedges = component_edges_all(inst, pair, edge_offset);
            vector<int> conflicts;
            for (int k = 0; k < (int)base.size(); k++) {
                if (k == i || k == j) continue;
                bool hit = false;
                int p1 = 0, p2 = 0;
                while (p1 < (int)pedges.size() && p2 < (int)base_edges[k].size()) {
                    if (pedges[p1] == base_edges[k][p2]) { hit = true; break; }
                    if (pedges[p1] < base_edges[k][p2]) p1++; else p2++;
                }
                if (hit) conflicts.push_back(k);
            }
            if ((int)conflicts.size() > 4) continue;
            bool huge_conflict = false;
            for (int k : conflicts) if (base[k].size() > 20) huge_conflict = true;
            if (huge_conflict) continue;
            for (int k : conflicts) conflict_score[k]++;
            pair_infos.push_back(PairInfo{(int)conflicts.size(), conflicts, i, j, pair});
        }
    }

    vector<pair<int, int>> scores;
    for (auto &kv : conflict_score) scores.push_back({kv.second, kv.first});
    sort(scores.begin(), scores.end(), [](const auto &a, const auto &b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });
    for (int i = 0; i < (int)scores.size() && i < 6; i++) selected.insert(scores[i].second);

    sort(pair_infos.begin(), pair_infos.end(), [&](const PairInfo &a, const PairInfo &b) {
        if (a.conflicts_size != b.conflicts_size) return a.conflicts_size < b.conflicts_size;
        int sa = conflict_score[a.left] + conflict_score[a.right];
        int sb = conflict_score[b.left] + conflict_score[b.right];
        if (sa != sb) return sa > sb;
        return a.pair < b.pair;
    });

    unordered_set<int> used_singletons;
    for (const PairInfo &info : pair_infos) {
        unordered_set<int> candidate = selected;
        candidate.insert(info.left);
        candidate.insert(info.right);
        for (int k : info.conflicts) candidate.insert(k);
        if ((int)candidate.size() > component_limit) continue;
        if (used_singletons.count(info.left) || used_singletons.count(info.right)) continue;
        selected = std::move(candidate);
        used_singletons.insert(info.left);
        used_singletons.insert(info.right);
        if ((int)selected.size() >= component_limit) break;
    }

    if ((int)selected.size() < min(component_limit, (int)base.size())) {
        vector<int> fill;
        for (int i = 0; i < (int)base.size(); i++) {
            if (!selected.count(i)) fill.push_back(i);
        }
        sort(fill.begin(), fill.end(), [&](int a, int b) {
            if (base[a].size() != base[b].size()) return base[a].size() > base[b].size();
            return base[a] < base[b];
        });
        for (int idx : fill) {
            if ((int)selected.size() >= component_limit) break;
            selected.insert(idx);
        }
    }

    vector<int> out(selected.begin(), selected.end());
    sort(out.begin(), out.end());
    return out;
}

static bool add_local_component(Pool &pool, const Instance &inst, const unordered_set<uint32_t> &valid,
        const vector<int> &edge_offset, vector<int> block, unordered_map<string, bool> &valid_cache) {
    sort(block.begin(), block.end());
    block.erase(unique(block.begin(), block.end()), block.end());
    if (block.empty()) return false;
    string key = key_of(block);
    if (pool.id.count(key)) return false;
    if (!triple_compatible(block, valid, inst.leaf_count)) return false;
    auto it = valid_cache.find(key);
    bool ok;
    if (it == valid_cache.end()) {
        ok = valid_topology(inst, block);
        valid_cache[key] = ok;
    } else {
        ok = it->second;
    }
    if (!ok) return false;
    vector<int> edges = component_edges_all(inst, block, edge_offset);
    double score = (double)(block.size() - 1) / max(1, (int)edges.size());
    int idx = (int)pool.comps.size();
    pool.id[key] = idx;
    pool.comps.push_back(Component{std::move(block), std::move(edges), score});
    return true;
}

static Pool build_local_pool(const Instance &inst, const Options &opt, const vector<int> &edge_offset,
        const unordered_set<uint32_t> &valid, const vector<int> &selected_leaves,
        const vector<vector<int>> &seed_blocks, const unordered_set<int> &fixed_edges,
        mt19937_64 &rng, Clock::time_point deadline) {
    Pool pool;
    unordered_map<string, bool> valid_cache;
    vector<int> leaves = selected_leaves;
    sort(leaves.begin(), leaves.end());
    for (int leaf : leaves) add_local_component(pool, inst, valid, edge_offset, vector<int>{leaf}, valid_cache);
    for (int i = 0; i < (int)leaves.size(); i++) {
        if (Clock::now() >= deadline) return pool;
        for (int j = i + 1; j < (int)leaves.size(); j++) {
            add_local_component(pool, inst, valid, edge_offset, vector<int>{leaves[i], leaves[j]}, valid_cache);
        }
    }

    vector<array<int, 3>> triples;
    for (int i = 0; i < (int)leaves.size(); i++) for (int j = i + 1; j < (int)leaves.size(); j++) for (int k = j + 1; k < (int)leaves.size(); k++) {
        triples.push_back({leaves[i], leaves[j], leaves[k]});
    }
    shuffle(triples.begin(), triples.end(), rng);
    int triple_limit = opt.local_seed_triples_limit < 0 ? (int)triples.size() : min(opt.local_seed_triples_limit, (int)triples.size());
    for (int i = 0; i < triple_limit; i++) {
        if (Clock::now() >= deadline) return pool;
        auto tri = triples[i];
        if (valid.count(triple_key(tri[0], tri[1], tri[2], inst.leaf_count))) {
            add_local_component(pool, inst, valid, edge_offset, vector<int>{tri[0], tri[1], tri[2]}, valid_cache);
        }
    }

    vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
    for (int x : leaves) selected_mark[x] = 1;
    for (const Tree &tree : inst.trees) {
        vector<vector<int>> sets;
        collect_leaf_sets(tree, tree.root, sets);
        for (const vector<int> &s : sets) {
            if (Clock::now() >= deadline) return pool;
            vector<int> block;
            for (int x : s) if (selected_mark[x]) block.push_back(x);
            if ((int)block.size() >= 2 && (int)block.size() <= opt.local_max_component_size) {
                add_local_component(pool, inst, valid, edge_offset, block, valid_cache);
            }
        }
    }

    int variant_count = 0;
    for (const vector<int> &seed : seed_blocks) {
        if ((int)seed.size() >= 3 && (int)seed.size() <= 8) {
            for (int old_leaf : seed) {
                for (int new_leaf : leaves) {
                    if (Clock::now() >= deadline) return pool;
                    if (variant_count >= opt.local_variant_seed_limit) break;
                    if (binary_search(seed.begin(), seed.end(), new_leaf)) continue;
                    vector<int> cand;
                    for (int x : seed) if (x != old_leaf) cand.push_back(x);
                    cand.push_back(new_leaf);
                    variant_count++;
                    add_local_component(pool, inst, valid, edge_offset, cand, valid_cache);
                }
                if (variant_count >= opt.local_variant_seed_limit) break;
            }
        }
    }
    for (const vector<int> &seed : seed_blocks) {
        if (Clock::now() >= deadline || variant_count >= opt.local_variant_seed_limit) break;
        if ((int)seed.size() >= 3 && (int)seed.size() <= 16) {
            int max_subset = min({8, (int)seed.size(), opt.local_max_component_size});
            for (int mask = 1; mask < (1 << (int)seed.size()); mask++) {
                if (Clock::now() >= deadline || variant_count >= opt.local_variant_seed_limit) break;
                int bits = __builtin_popcount((unsigned)mask);
                if (bits < 2 || bits > max_subset) continue;
                vector<int> subset;
                for (int i = 0; i < (int)seed.size(); i++) if (mask & (1 << i)) subset.push_back(seed[i]);
                variant_count++;
                add_local_component(pool, inst, valid, edge_offset, subset, valid_cache);
            }
        }
    }

    uniform_real_distribution<double> prob(0.0, 1.0);
    for (int trial = 0; trial < opt.local_guided_growth_trials &&
            (int)pool.comps.size() < opt.local_max_pool && Clock::now() < deadline; trial++) {
        vector<int> current;
        if (!seed_blocks.empty() && prob(rng) < 0.35) {
            const vector<int> &seed = seed_blocks[rng() % seed_blocks.size()];
            if (seed.size() <= 3) current = seed;
            else current = random_sample(seed, prob(rng) < 0.5 ? 2 : 3, rng);
        } else {
            current = random_sample(leaves, prob(rng) < 0.75 ? 2 : 3, rng);
        }
        sort(current.begin(), current.end());
        if (!triple_compatible(current, valid, inst.leaf_count)) continue;
        if (!valid_topology(inst, current)) continue;
        add_local_component(pool, inst, valid, edge_offset, current, valid_cache);
        while ((int)current.size() < opt.local_max_component_size && Clock::now() < deadline) {
            vector<int> order = leaves;
            shuffle(order.begin(), order.end(), rng);
            vector<pair<double, vector<int>>> options;
            unordered_set<int> current_set(current.begin(), current.end());
            for (int leaf : order) {
                if (current_set.count(leaf)) continue;
                vector<int> candidate = current;
                candidate.push_back(leaf);
                sort(candidate.begin(), candidate.end());
                if (!triple_compatible(candidate, valid, inst.leaf_count)) continue;
                string key = key_of(candidate);
                bool ok;
                auto it = valid_cache.find(key);
                if (it == valid_cache.end()) {
                    ok = valid_topology(inst, candidate);
                    valid_cache[key] = ok;
                } else ok = it->second;
                if (!ok) continue;
                int future = 0;
                vector<int> probe = leaves;
                shuffle(probe.begin(), probe.end(), rng);
                int checked = 0;
                unordered_set<int> candidate_set(candidate.begin(), candidate.end());
                for (int other : probe) {
                    if (candidate_set.count(other)) continue;
                    vector<int> next = candidate;
                    next.push_back(other);
                    sort(next.begin(), next.end());
                    if (triple_compatible(next, valid, inst.leaf_count)) future++;
                    checked++;
                    if (checked >= opt.local_guided_future_probe) break;
                }
                double score = future + prob(rng);
                options.push_back({score, std::move(candidate)});
            }
            if (options.empty()) break;
            sort(options.begin(), options.end(), [](const auto &a, const auto &b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
            int pick_window = min(opt.local_guided_choice_window, (int)options.size());
            current = options[rng() % pick_window].second;
            add_local_component(pool, inst, valid, edge_offset, current, valid_cache);
        }
    }

    for (int trial = 0; trial < opt.local_trials && (int)pool.comps.size() < opt.local_max_pool && Clock::now() < deadline; trial++) {
        int start_size = prob(rng) < 0.75 ? 2 : 3;
        vector<int> current = random_sample(leaves, min(start_size, (int)leaves.size()), rng);
        if (!triple_compatible(current, valid, inst.leaf_count)) continue;
        vector<int> order = leaves;
        shuffle(order.begin(), order.end(), rng);
        for (int leaf : order) {
            if (Clock::now() >= deadline || (int)current.size() >= opt.local_max_component_size) break;
            if (binary_search(current.begin(), current.end(), leaf)) continue;
            vector<int> cand = current;
            cand.push_back(leaf);
            sort(cand.begin(), cand.end());
            if (!triple_compatible(cand, valid, inst.leaf_count)) continue;
            string key = key_of(cand);
            bool ok;
            auto it = valid_cache.find(key);
            if (it == valid_cache.end()) {
                ok = valid_topology(inst, cand);
                valid_cache[key] = ok;
            } else ok = it->second;
            if (ok) {
                current = cand;
                if (prob(rng) < 0.55) add_local_component(pool, inst, valid, edge_offset, current, valid_cache);
            }
        }
        add_local_component(pool, inst, valid, edge_offset, current, valid_cache);
    }

    Pool filtered;
    for (const Component &c : pool.comps) {
        bool conflict = false;
        for (int e : c.edges) if (fixed_edges.count(e)) { conflict = true; break; }
        if (conflict) continue;
        filtered.id[key_of(c.leaves)] = (int)filtered.comps.size();
        filtered.comps.push_back(c);
    }
    return filtered;
}

static vector<Component> reduce_local_columns(const vector<Component> &columns,
        const unordered_map<int, int> &leaf_to_local, int local_leaf_count, int per_leaf) {
    vector<vector<int>> by_leaf(local_leaf_count);
    for (int i = 0; i < (int)columns.size(); i++) {
        for (int leaf : columns[i].leaves) {
            auto it = leaf_to_local.find(leaf);
            if (it != leaf_to_local.end()) by_leaf[it->second].push_back(i);
        }
    }
    auto better = [&](int a, int b) {
        const Component &ca = columns[a], &cb = columns[b];
        double ra = (double)((int)ca.leaves.size() - 1) / max(1, (int)ca.edges.size());
        double rb = (double)((int)cb.leaves.size() - 1) / max(1, (int)cb.edges.size());
        if (ra != rb) return ra > rb;
        if (ca.leaves.size() != cb.leaves.size()) return ca.leaves.size() > cb.leaves.size();
        if (ca.edges.size() != cb.edges.size()) return ca.edges.size() < cb.edges.size();
        return ca.leaves < cb.leaves;
    };
    unordered_set<int> chosen;
    for (auto &bucket : by_leaf) {
        sort(bucket.begin(), bucket.end(), better);
        int keep = min(per_leaf, (int)bucket.size());
        for (int i = 0; i < keep; i++) chosen.insert(bucket[i]);
    }
    vector<int> order(chosen.begin(), chosen.end());
    sort(order.begin(), order.end(), better);
    vector<Component> out;
    out.reserve(order.size());
    for (int idx : order) out.push_back(columns[idx]);
    return out;
}

#ifdef USE_HIGHS
static vector<vector<int>> solve_local_highs(const vector<Component> &columns,
        const vector<int> &selected_leaves, int max_local_components, double seconds) {
    if (columns.empty()) return {};
    unordered_map<int, int> leaf_row;
    for (int i = 0; i < (int)selected_leaves.size(); i++) leaf_row[selected_leaves[i]] = i;
    unordered_map<int, int> edge_row;
    vector<int> starts;
    vector<int> indices;
    vector<double> values;
    vector<double> cost, lower, upper;
    vector<HighsVarType> integrality;
    starts.push_back(0);
    for (const Component &c : columns) {
        cost.push_back(-(double)((int)c.leaves.size() - 1));
        lower.push_back(0.0);
        upper.push_back(1.0);
        integrality.push_back(HighsVarType::kInteger);
        for (int leaf : c.leaves) {
            indices.push_back(leaf_row.at(leaf));
            values.push_back(1.0);
        }
        for (int edge : c.edges) {
            auto it = edge_row.find(edge);
            if (it == edge_row.end()) {
                int row = (int)selected_leaves.size() + (int)edge_row.size();
                it = edge_row.insert({edge, row}).first;
            }
            indices.push_back(it->second);
            values.push_back(1.0);
        }
        starts.push_back((int)indices.size());
    }
    int row_count = (int)selected_leaves.size() + (int)edge_row.size();
    HighsModel model;
    model.lp_.num_col_ = (int)columns.size();
    model.lp_.num_row_ = row_count;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.col_cost_ = std::move(cost);
    model.lp_.col_lower_ = std::move(lower);
    model.lp_.col_upper_ = std::move(upper);
    model.lp_.row_lower_.assign(row_count, 0.0);
    model.lp_.row_upper_.assign(row_count, 1.0);
    model.lp_.integrality_ = std::move(integrality);
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = std::move(starts);
    model.lp_.a_matrix_.index_ = std::move(indices);
    model.lp_.a_matrix_.value_ = std::move(values);

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("time_limit", max(0.25, seconds));
    highs.setOptionValue("mip_rel_gap", 0.0);
    HighsStatus pass_status = highs.passModel(model);
    if (pass_status != HighsStatus::kOk) return {};
    HighsStatus run_status = highs.run();
    HighsModelStatus model_status = highs.getModelStatus();
    if (run_status != HighsStatus::kOk || model_status != HighsModelStatus::kOptimal) {
        cerr << "HIGHS_SET_PACK_STATUS run_status=" << (int)run_status
             << " model_status=" << (int)model_status
             << " cols=" << columns.size()
             << " leaves=" << selected_leaves.size()
             << " max_components=" << max_local_components << "\n";
    }
    const HighsSolution &sol = highs.getSolution();
    if ((int)sol.col_value.size() != (int)columns.size()) return {};

    vector<vector<int>> local;
    vector<unsigned char> covered(leaf_row.size() + 1, 0);
    int savings = 0;
    for (int i = 0; i < (int)columns.size(); i++) {
        if (sol.col_value[i] > 0.5) {
            local.push_back(columns[i].leaves);
            savings += (int)columns[i].leaves.size() - 1;
            for (int leaf : columns[i].leaves) covered[leaf_row[leaf]] = 1;
        }
    }
    int target_savings = (int)selected_leaves.size() - max_local_components;
    if (savings < target_savings) return {};
    for (int leaf : selected_leaves) {
        int row = leaf_row[leaf];
        if (!covered[row]) local.push_back(vector<int>{leaf});
    }
    sort(local.begin(), local.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return local;
}
#else
static vector<vector<int>> solve_local_highs(const vector<Component> &, const vector<int> &, int, double) {
    return {};
}
#endif

static vector<vector<int>> local_branch_cover(const vector<Component> &columns,
        const vector<int> &selected_leaves, int max_components, double seconds, int branch_width) {
    int leaf_count = (int)selected_leaves.size();
    int target_savings = leaf_count - max_components;
    if (target_savings <= 0) {
        vector<vector<int>> out;
        for (int leaf : selected_leaves) out.push_back(vector<int>{leaf});
        return out;
    }
    if (columns.empty() || seconds <= 0.02) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));
    unordered_map<int, int> leaf_row;
    for (int i = 0; i < leaf_count; i++) leaf_row[selected_leaves[i]] = i;
    unordered_map<int, int> edge_row;
    struct BCol {
        vector<int> leaves;
        vector<unsigned long long> leaf_mask;
        vector<unsigned long long> edge_mask;
        int weight = 0;
        int edge_count = 0;
    };
    vector<BCol> bcols;
    bcols.reserve(columns.size());
    int leaf_words = (leaf_count + 63) >> 6;
    for (const Component &c : columns) {
        BCol bc;
        bc.leaves = c.leaves;
        bc.leaf_mask.assign(leaf_words, 0);
        bc.weight = (int)c.leaves.size() - 1;
        for (int leaf : c.leaves) {
            auto it = leaf_row.find(leaf);
            if (it == leaf_row.end()) {
                bc.leaf_mask.clear();
                break;
            }
            bc.leaf_mask[it->second >> 6] |= 1ULL << (it->second & 63);
        }
        if (bc.leaf_mask.empty()) continue;
        vector<int> compact_edges;
        compact_edges.reserve(c.edges.size());
        for (int edge : c.edges) {
            auto it = edge_row.find(edge);
            if (it == edge_row.end()) {
                int row = (int)edge_row.size();
                it = edge_row.insert({edge, row}).first;
            }
            compact_edges.push_back(it->second);
        }
        bc.edge_count = (int)compact_edges.size();
        int edge_words = ((int)edge_row.size() + 63) >> 6;
        bc.edge_mask.assign(edge_words, 0);
        for (int row : compact_edges) {
            if (row >= (int)bc.edge_mask.size() * 64) bc.edge_mask.resize((row + 64) >> 6, 0);
            bc.edge_mask[row >> 6] |= 1ULL << (row & 63);
        }
        bcols.push_back(std::move(bc));
    }
    int edge_words = ((int)edge_row.size() + 63) >> 6;
    for (BCol &bc : bcols) bc.edge_mask.resize(edge_words, 0);
    vector<vector<int>> by_leaf(leaf_count);
    for (int i = 0; i < (int)bcols.size(); i++) {
        for (int w = 0; w < leaf_words; w++) {
            unsigned long long bits = bcols[i].leaf_mask[w];
            while (bits) {
                unsigned long long low = bits & -bits;
                int row = (w << 6) + __builtin_ctzll(bits);
                if (row < leaf_count) by_leaf[row].push_back(i);
                bits ^= low;
            }
        }
    }
    auto better = [&](int a, int b) {
        double sa = bcols[a].weight / pow((double)max(1, bcols[a].edge_count), 0.20);
        double sb = bcols[b].weight / pow((double)max(1, bcols[b].edge_count), 0.20);
        if (bcols[a].weight != bcols[b].weight) return bcols[a].weight > bcols[b].weight;
        if (sa != sb) return sa > sb;
        return bcols[a].leaves < bcols[b].leaves;
    };
    for (auto &bucket : by_leaf) sort(bucket.begin(), bucket.end(), better);

    auto edge_intersects = [&](const vector<unsigned long long> &a, const vector<unsigned long long> &b) {
        for (int i = 0; i < (int)a.size(); i++) if (a[i] & b[i]) return true;
        return false;
    };
    auto leaf_intersects = [&](const vector<unsigned long long> &a, const vector<unsigned long long> &b) {
        for (int i = 0; i < leaf_words; i++) if (a[i] & b[i]) return true;
        return false;
    };
    auto or_into = [](vector<unsigned long long> &a, const vector<unsigned long long> &b) {
        for (int i = 0; i < (int)a.size(); i++) a[i] |= b[i];
    };
    auto assigned_count = [&](const vector<unsigned long long> &assigned) {
        int total = 0;
        for (int i = 0; i < leaf_words; i++) total += __builtin_popcountll(assigned[i]);
        return total;
    };

    vector<int> chosen, answer;
    int best_savings = 0;
    long long nodes = 0;
    function<bool(vector<unsigned long long>&, vector<unsigned long long>&, int)> dfs =
        [&](vector<unsigned long long> &assigned, vector<unsigned long long> &used_edges, int savings) -> bool {
            if (Clock::now() >= deadline) return false;
            nodes++;
            if (savings > best_savings) best_savings = savings;
            if (savings >= target_savings) {
                answer = chosen;
                return true;
            }
            int count = assigned_count(assigned);
            if (savings + max(0, leaf_count - count - 1) < target_savings) return false;

            int best_leaf = -1;
            vector<int> options;
            int best_count = INT_MAX;
            for (int leaf = 0; leaf < leaf_count; leaf++) {
                if ((assigned[leaf >> 6] >> (leaf & 63)) & 1ULL) continue;
                vector<int> local;
                int feasible_count = 0;
                for (int idx : by_leaf[leaf]) {
                    if (leaf_intersects(assigned, bcols[idx].leaf_mask)) continue;
                    if (edge_intersects(used_edges, bcols[idx].edge_mask)) continue;
                    feasible_count++;
                    if ((int)local.size() < branch_width) local.push_back(idx);
                }
                if (feasible_count < best_count) {
                    best_count = feasible_count;
                    best_leaf = leaf;
                    options = std::move(local);
                    if (best_count == 0) break;
                }
            }
            if (best_leaf < 0) return false;
            for (int idx : options) {
                vector<unsigned long long> next_assigned = assigned;
                vector<unsigned long long> next_edges = used_edges;
                or_into(next_assigned, bcols[idx].leaf_mask);
                or_into(next_edges, bcols[idx].edge_mask);
                chosen.push_back(idx);
                if (dfs(next_assigned, next_edges, savings + bcols[idx].weight)) return true;
                chosen.pop_back();
            }
            assigned[best_leaf >> 6] |= 1ULL << (best_leaf & 63);
            bool ok = dfs(assigned, used_edges, savings);
            assigned[best_leaf >> 6] &= ~(1ULL << (best_leaf & 63));
            return ok;
        };
    vector<unsigned long long> assigned(leaf_words, 0), used_edges(edge_words, 0);
    bool hit = dfs(assigned, used_edges, 0);
    cerr << "LOCAL_BRANCH_DONE hit=" << (hit ? 1 : 0)
         << " best_savings=" << best_savings
         << " target_savings=" << target_savings
         << " cols=" << bcols.size()
         << " leaves=" << leaf_count
         << " nodes=" << nodes << "\n";
    if (!hit) return {};
    int max_label = 0;
    for (int leaf : selected_leaves) max_label = max(max_label, leaf);
    vector<unsigned char> used_leaf(max_label + 1, 0);
    vector<vector<int>> out;
    for (int idx : answer) {
        out.push_back(bcols[idx].leaves);
        for (int leaf : bcols[idx].leaves) used_leaf[leaf] = 1;
    }
    for (int leaf : selected_leaves) if (!used_leaf[leaf]) out.push_back(vector<int>{leaf});
    sort(out.begin(), out.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return out;
}

static vector<vector<int>> local_repack_repair(const Instance &inst, vector<vector<int>> partition,
        const Options &opt, const vector<int> &edge_offset, const vector<Column> &global_cols) {
    if (opt.local_seconds <= 0.0 || (int)partition.size() <= opt.target_components) return partition;
    auto start = Clock::now();
    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    mt19937_64 rng(opt.seed + 1000003);
    vector<vector<int>> best = partition;
    int best_count = (int)partition.size();
    int rounds = 0;
    while (elapsed(start) < opt.local_seconds && best_count > opt.target_components) {
        rounds++;
        vector<vector<int>> base_edges;
        for (const vector<int> &block : best) base_edges.push_back(component_edges_all(inst, block, edge_offset));
        vector<int> selected_indices;
        if (rounds == 1 && opt.local_guided_components > 0) {
            selected_indices = guided_local_indices(inst, best, base_edges, edge_offset, opt.local_guided_components);
        }
        if (selected_indices.empty()) {
            vector<int> order(best.size());
            for (int i = 0; i < (int)order.size(); i++) order[i] = i;
            sort(order.begin(), order.end(), [&](int a, int b) {
                if (best[a].size() != best[b].size()) return best[a].size() > best[b].size();
                return best[a] < best[b];
            });
            int anchor_count = max(1, min(opt.local_random_anchor_components, (int)order.size()));
            int anchor = order[rng() % anchor_count];
            vector<int> outside;
            for (int i = 0; i < (int)best.size(); i++) if (i != anchor) outside.push_back(i);
            shuffle(outside.begin(), outside.end(), rng);
            selected_indices.push_back(anchor);
            int take = min(opt.local_random_extra_components, (int)outside.size());
            for (int i = 0; i < take; i++) selected_indices.push_back(outside[i]);
            sort(selected_indices.begin(), selected_indices.end());
        }
        if (selected_indices.empty()) break;

        vector<unsigned char> remove(best.size(), 0);
        vector<int> selected_leaves;
        vector<vector<int>> seed_blocks;
        unordered_set<int> fixed_edges;
        vector<vector<int>> fixed;
        for (int idx : selected_indices) remove[idx] = 1;
        for (int i = 0; i < (int)best.size(); i++) {
            if (remove[i]) {
                seed_blocks.push_back(best[i]);
                selected_leaves.insert(selected_leaves.end(), best[i].begin(), best[i].end());
            } else {
                fixed.push_back(best[i]);
                for (int e : base_edges[i]) fixed_edges.insert(e);
            }
        }
        sort(selected_leaves.begin(), selected_leaves.end());
        selected_leaves.erase(unique(selected_leaves.begin(), selected_leaves.end()), selected_leaves.end());
        int max_local_components = opt.target_components - (int)fixed.size();
        if (max_local_components < 1) break;

        double remaining_total = max(0.1, opt.local_seconds - elapsed(start));
        double round_budget = opt.local_round_seconds > 0.0 ? min(opt.local_round_seconds, remaining_total) : remaining_total;
        auto deadline = Clock::now() + chrono::milliseconds((int)(round_budget * 1000));
        Pool local_pool = build_local_pool(inst, opt, edge_offset, valid, selected_leaves, seed_blocks, fixed_edges, rng, deadline);
        vector<Component> columns;
        unordered_set<string> column_keys;
        for (const Component &c : local_pool.comps) {
            if (c.leaves.size() <= 1) continue;
            if ((int)c.leaves.size() > opt.local_column_max_size) continue;
            if ((int)c.edges.size() > opt.local_column_max_edges) continue;
            columns.push_back(c);
            column_keys.insert(key_of(c.leaves));
        }
        vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
        for (int leaf : selected_leaves) selected_mark[leaf] = 1;
        for (const Column &col : global_cols) {
            if (col.leaves.size() <= 1) continue;
            if ((int)col.leaves.size() > opt.local_column_max_size) continue;
            if ((int)col.edges.size() > opt.local_column_max_edges) continue;
            bool subset = true;
            for (int leaf : col.leaves) {
                if (!selected_mark[leaf]) { subset = false; break; }
            }
            if (!subset) continue;
            bool conflict = false;
            for (int edge : col.edges) {
                if (fixed_edges.count(edge)) { conflict = true; break; }
            }
            if (conflict) continue;
            string key = key_of(col.leaves);
            if (column_keys.count(key)) continue;
            column_keys.insert(key);
            double score = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
            columns.push_back(Component{col.leaves, col.edges, score});
        }
        unordered_map<int, int> leaf_to_local;
        for (int i = 0; i < (int)selected_leaves.size(); i++) leaf_to_local[selected_leaves[i]] = i;
        columns = reduce_local_columns(columns, leaf_to_local, (int)selected_leaves.size(), opt.local_top_per_leaf);
        double remaining = max(0.1, min(round_budget, opt.local_seconds - elapsed(start)));
        vector<vector<int>> local = solve_local_highs(columns, selected_leaves, max_local_components, remaining);
        cerr << "LOCAL_ROUND round=" << rounds
             << " base=" << best_count
             << " selected_components=" << selected_indices.size()
             << " selected_leaves=" << selected_leaves.size()
             << " pool=" << local_pool.comps.size()
             << " columns=" << columns.size()
             << " max_local=" << max_local_components
             << " local=" << local.size() << "\n";
        if (local.empty()) continue;
        vector<vector<int>> candidate = fixed;
        candidate.insert(candidate.end(), local.begin(), local.end());
        sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        if (validate_partition(inst, candidate, edge_offset) && (int)candidate.size() < best_count) {
            best_count = (int)candidate.size();
            best = std::move(candidate);
            cerr << "LOCAL_BEST components=" << best_count << "\n";
        } else {
            continue;
        }
    }
    return best;
}

static vector<int> reduce_global_branch_candidates(const vector<Column> &cols, int n, const Options &opt);
static vector<vector<int>> pair_augment_repair(const Instance &inst, vector<vector<int>> partition,
        const Options &opt, const vector<int> &edge_offset, const vector<Column> &cols);
static vector<vector<int>> exact_region_exchange_repair(const Instance &inst,
        vector<vector<int>> partition, const Options &opt,
        const vector<int> &edge_offset, const vector<Column> &cols);

static vector<vector<int>> global_highs_repack(const Instance &inst, const Options &opt,
        const vector<Column> &global_cols) {
    if (opt.global_highs_seconds <= 0.0) return {};
    vector<int> all_leaves(inst.leaf_count);
    for (int i = 0; i < inst.leaf_count; i++) all_leaves[i] = i + 1;
    vector<int> candidate_ids = reduce_global_branch_candidates(global_cols, inst.leaf_count, opt);
    vector<Component> columns;
    columns.reserve(candidate_ids.size());
    unordered_set<string> seen;
    for (int idx : candidate_ids) {
        const Column &col = global_cols[idx];
        if (col.leaves.size() <= 1) continue;
        string key = key_of(col.leaves);
        if (!seen.insert(key).second) continue;
        double score = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
        columns.push_back(Component{col.leaves, col.edges, score});
    }
    cerr << "GLOBAL_HIGHS_REPACK_START columns=" << columns.size()
         << " target=" << opt.target_components << "\n";
    vector<vector<int>> packed = solve_local_highs(columns, all_leaves, opt.target_components, opt.global_highs_seconds);
    cerr << "GLOBAL_HIGHS_REPACK columns=" << columns.size()
         << " result=" << packed.size()
         << " target=" << opt.target_components << "\n";
    return packed;
}

static vector<vector<int>> read_component_csv(const string &path);

static void report_branch_candidate_check(const vector<Column> &cols, const vector<int> &candidates,
        const string &path) {
    if (path.empty()) return;
    unordered_set<string> present;
    present.reserve(candidates.size() * 2);
    for (int idx : candidates) present.insert(key_of(cols[idx].leaves));
    vector<vector<int>> rows = read_component_csv(path);
    int need = 0, hit = 0;
    for (const vector<int> &block : rows) {
        if (block.size() <= 1) continue;
        need++;
        bool ok = present.count(key_of(block)) > 0;
        if (ok) hit++;
        cerr << "BRANCH_CANDIDATE_CHECK " << (ok ? "hit" : "miss")
             << " size=" << block.size()
             << " leaves=";
        for (int i = 0; i < (int)block.size(); i++) {
            if (i) cerr << ',';
            cerr << block[i];
        }
        cerr << "\n";
    }
    cerr << "BRANCH_CANDIDATE_CHECK_SUMMARY hit=" << hit << " need=" << need
         << " candidates=" << candidates.size() << "\n";
}

static vector<int> reduce_global_branch_candidates(const vector<Column> &cols, int n, const Options &opt) {
    vector<int> always_keep;
    vector<int> large;
    for (int i = 0; i < (int)cols.size(); i++) {
        const Column &col = cols[i];
        if (col.leaves.size() <= 1) continue;
        if ((int)col.leaves.size() > opt.global_branch_column_max_size) continue;
        if ((int)col.edges.size() > opt.global_branch_column_max_edges) continue;
        if (col.leaves.size() <= 2) always_keep.push_back(i);
        else large.push_back(i);
    }

    auto better = [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        double da = (double)wa / max(1, (int)cols[a].edges.size());
        double db = (double)wb / max(1, (int)cols[b].edges.size());
        if (da != db) return da > db;
        if (wa != wb) return wa > wb;
        if (cols[a].leaves.size() != cols[b].leaves.size()) return cols[a].leaves.size() > cols[b].leaves.size();
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };

    unordered_set<int> chosen(always_keep.begin(), always_keep.end());
    if (opt.global_branch_top_per_leaf > 0) {
        vector<vector<int>> by_leaf(n + 1);
        for (int idx : large) for (int leaf : cols[idx].leaves) by_leaf[leaf].push_back(idx);
        for (int leaf = 1; leaf <= n; leaf++) {
            auto &bucket = by_leaf[leaf];
            sort(bucket.begin(), bucket.end(), better);
            int keep = min(opt.global_branch_top_per_leaf, (int)bucket.size());
            for (int i = 0; i < keep; i++) chosen.insert(bucket[i]);
        }
    } else {
        chosen.insert(large.begin(), large.end());
    }

    vector<int> out(chosen.begin(), chosen.end());
    sort(out.begin(), out.end(), better);
    return out;
}

static vector<int> global_branch_pack(const vector<Column> &cols, int n, const Options &opt) {
    if (opt.global_branch_seconds <= 0.0 || !opt.has_target || cols.empty()) return {};
    auto branch_start = Clock::now();
    int target_savings = n - opt.target_components;
    if (target_savings <= 0) return {};
    vector<int> candidates = reduce_global_branch_candidates(cols, n, opt);
    if (candidates.empty()) return {};
    report_branch_candidate_check(cols, candidates, opt.check_components);

#ifdef USE_CADICAL
    if (target_savings <= 64 && (int)candidates.size() <= 8000 &&
            opt.global_branch_seconds > 1.0) {
        auto sat_start = Clock::now();
        double sat_seconds = min(60.0, max(0.5, opt.global_branch_seconds * 0.65));
        Clock::time_point sat_deadline =
            sat_start + chrono::milliseconds((int)(sat_seconds * 1000));
        CaDiCaL::Solver solver;
        solver.set("quiet", 1);
        solver.set("factor", 0);
        DeadlineTerminator terminator(sat_deadline);
        solver.connect_terminator(&terminator);
        int var_count = (int)candidates.size();
        long long clause_count = 0;
        auto add_clause = [&](const vector<int> &lits) {
            for (int lit : lits) solver.add(lit);
            solver.add(0);
            clause_count++;
        };
        auto add_at_most_one = [&](const vector<int> &vars) {
            int m = (int)vars.size();
            if (m <= 1) return;
            vector<int> seq(m, 0);
            for (int i = 1; i < m; i++) seq[i] = ++var_count;
            add_clause({-vars[0], seq[1]});
            for (int i = 2; i < m; i++) {
                add_clause({-seq[i - 1], seq[i]});
                add_clause({-vars[i - 1], -seq[i - 1]});
                add_clause({-vars[i - 1], seq[i]});
            }
            add_clause({-vars[m - 1], -seq[m - 1]});
        };
        auto add_weighted_at_least = [&](const vector<int> &vars,
                const vector<int> &weights, int need) {
            int m = (int)vars.size();
            if (need <= 0) return;
            int total_weight = 0;
            for (int w : weights) total_weight += w;
            if (total_weight < need) {
                add_clause({});
                return;
            }
            vector<vector<int>> dp(m, vector<int>(need + 1, 0));
            int prefix = 0;
            for (int i = 0; i < m; i++) {
                prefix += weights[i];
                for (int s = 1; s <= min(need, prefix); s++) {
                    dp[i][s] = ++var_count;
                }
            }
            for (int i = 0; i < m; i++) {
                int w = min(weights[i], need);
                for (int s = 1; s <= need; s++) {
                    if (!dp[i][s]) continue;
                    int same = i > 0 ? dp[i - 1][s] : 0;
                    int prev_need = max(0, s - w);
                    int take_prev = (i > 0 && prev_need > 0) ?
                        dp[i - 1][prev_need] : 0;
                    if (same) add_clause({-same, dp[i][s]});
                    if (prev_need == 0) {
                        add_clause({-vars[i], dp[i][s]});
                    } else if (take_prev) {
                        add_clause({-vars[i], -take_prev, dp[i][s]});
                    }
                    vector<int> reverse_choice{-dp[i][s]};
                    if (same) reverse_choice.push_back(same);
                    reverse_choice.push_back(vars[i]);
                    add_clause(reverse_choice);
                    if (prev_need > 0) {
                        vector<int> reverse_prev{-dp[i][s]};
                        if (same) reverse_prev.push_back(same);
                        if (take_prev) reverse_prev.push_back(take_prev);
                        add_clause(reverse_prev);
                    }
                }
            }
            add_clause({dp[m - 1][need]});
        };

        vector<int> vars, weights;
        vars.reserve(candidates.size());
        weights.reserve(candidates.size());
        vector<vector<int>> by_bit(cols[0].mask.size() * 64);
        for (int i = 0; i < (int)candidates.size(); i++) {
            int var = i + 1;
            int idx = candidates[i];
            vars.push_back(var);
            weights.push_back(column_weight(cols[idx]));
            for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
                unsigned long long bits = cols[idx].mask[w];
                while (bits) {
                    int bit = __builtin_ctzll(bits);
                    int global_bit = (w << 6) + bit;
                    if (global_bit < (int)by_bit.size()) by_bit[global_bit].push_back(var);
                    bits &= bits - 1;
                }
            }
        }
        for (const vector<int> &bucket : by_bit) add_at_most_one(bucket);
        add_weighted_at_least(vars, weights, target_savings);
        cerr << "GLOBAL_CADICAL_PACK_SOLVE target_savings=" << target_savings
             << " candidates=" << candidates.size()
             << " vars=" << var_count
             << " clauses=" << clause_count
             << " build_seconds=" << elapsed(sat_start) << "\n";
        int sat = solver.solve();
        cerr << "GLOBAL_CADICAL_PACK_STATUS result=" << sat
             << " target_savings=" << target_savings
             << " candidates=" << candidates.size()
             << " seconds=" << elapsed(sat_start) << "\n";
        if (sat == 10) {
            vector<int> selected;
            int savings = 0;
            vector<unsigned long long> used(cols[0].mask.size(), 0);
            for (int i = 0; i < (int)candidates.size(); i++) {
                if (solver.val(i + 1) <= 0) continue;
                int idx = candidates[i];
                if (intersects(used, cols[idx].mask)) {
                    selected.clear();
                    break;
                }
                bit_or_into(used, cols[idx].mask);
                selected.push_back(idx);
                savings += column_weight(cols[idx]);
            }
            if (!selected.empty() && savings >= target_savings) {
                solver.disconnect_terminator();
                cerr << "GLOBAL_CADICAL_PACK_HIT savings=" << savings
                     << " components=" << component_count(cols, selected, n)
                     << "\n";
                return selected;
            }
        }
        solver.disconnect_terminator();
    }
#endif

    int leaf_words = (n + 63) >> 6;
    vector<vector<int>> by_leaf(n + 1);
    for (int idx : candidates) {
        for (int leaf : cols[idx].leaves) by_leaf[leaf].push_back(idx);
    }
    auto better = [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        double sa = wa / pow((double)max(1, mask_bit_count(cols[a].mask)), 0.20);
        double sb = wb / pow((double)max(1, mask_bit_count(cols[b].mask)), 0.20);
        if (wa != wb) return wa > wb;
        if (sa != sb) return sa > sb;
        if (cols[a].leaves.size() != cols[b].leaves.size()) return cols[a].leaves.size() > cols[b].leaves.size();
        return cols[a].leaves < cols[b].leaves;
    };
    for (int leaf = 1; leaf <= n; leaf++) sort(by_leaf[leaf].begin(), by_leaf[leaf].end(), better);

    double branch_remaining = opt.global_branch_seconds - elapsed(branch_start);
    if (branch_remaining <= 0.02) return {};
    auto deadline = Clock::now() + chrono::milliseconds((int)(branch_remaining * 1000));
    vector<int> best, chosen;
    int best_savings = 0;
    long long nodes = 0;

    auto leaf_assigned_count = [&](const vector<unsigned long long> &assigned) {
        int total = 0;
        for (int i = 0; i < leaf_words; i++) total += __builtin_popcountll(assigned[i]);
        return total;
    };

    function<bool(vector<unsigned long long>&, vector<unsigned long long>&, int)> dfs =
        [&](vector<unsigned long long> &used, vector<unsigned long long> &assigned, int savings) -> bool {
            if (Clock::now() >= deadline) return false;
            nodes++;
            if (savings > best_savings) {
                best_savings = savings;
                best = chosen;
                cerr << "GLOBAL_BRANCH_BEST savings=" << best_savings
                     << " components=" << component_count(cols, best, n)
                     << " nodes=" << nodes << "\n";
            }
            if (savings >= target_savings) return true;
            int assigned_count = leaf_assigned_count(assigned);
            if (savings + max(0, n - assigned_count - 1) < target_savings) return false;

            int best_leaf = -1;
            vector<int> options;
            tuple<int, double, int> best_score{INT_MIN, -1e100, 0};
            for (int leaf = 1; leaf <= n; leaf++) {
                if ((assigned[(leaf - 1) >> 6] >> ((leaf - 1) & 63)) & 1ULL) continue;
                vector<int> local;
                for (int idx : by_leaf[leaf]) {
                    bool leaf_conflict = false;
                    for (int w = 0; w < leaf_words; w++) {
                        if (cols[idx].leaf_mask[w] & assigned[w]) { leaf_conflict = true; break; }
                    }
                    if (leaf_conflict || intersects(used, cols[idx].mask)) continue;
                    local.push_back(idx);
                    if ((int)local.size() >= opt.global_branch_width) break;
                }
                if (local.empty()) continue;
                int first = local[0];
                double first_score = column_weight(cols[first]) /
                    pow((double)max(1, mask_bit_count(cols[first].mask)), 0.20);
                tuple<int, double, int> score{-(int)local.size(), first_score, column_weight(cols[first])};
                if (best_leaf < 0 || score > best_score) {
                    best_leaf = leaf;
                    options = std::move(local);
                    best_score = score;
                }
            }

            if (best_leaf < 0) return false;
            int explored = 0;
            for (int idx : options) {
                if (Clock::now() >= deadline) return false;
                explored++;
                vector<unsigned long long> next_used = used;
                vector<unsigned long long> next_assigned = assigned;
                bit_or_into(next_used, cols[idx].mask);
                for (int w = 0; w < leaf_words; w++) next_assigned[w] |= cols[idx].leaf_mask[w];
                chosen.push_back(idx);
                if (dfs(next_used, next_assigned, savings + column_weight(cols[idx]))) return true;
                chosen.pop_back();
                if (explored >= opt.global_branch_width) break;
            }

            assigned[(best_leaf - 1) >> 6] |= 1ULL << ((best_leaf - 1) & 63);
            bool ok = dfs(used, assigned, savings);
            assigned[(best_leaf - 1) >> 6] &= ~(1ULL << ((best_leaf - 1) & 63));
            return ok;
        };

    vector<unsigned long long> used(cols[0].mask.size(), 0), assigned(cols[0].mask.size(), 0);
    bool hit = dfs(used, assigned, 0);
    cerr << "GLOBAL_BRANCH_DONE hit=" << (hit ? 1 : 0)
         << " best_savings=" << best_savings
         << " target_savings=" << target_savings
         << " candidates=" << candidates.size()
         << " nodes=" << nodes << "\n";
    if (hit) return chosen;
    return best;
}

static vector<array<int, 3>> sorted_valid_triples(const unordered_set<uint32_t> &valid, int n) {
    vector<array<int, 3>> triples;
    triples.reserve(valid.size());
    int base = n + 1;
    for (uint32_t key : valid) {
        int c = key % base;
        int tmp = key / base;
        int b = tmp % base;
        int a = tmp / base;
        if (1 <= a && a < b && b < c && c <= n) triples.push_back({a, b, c});
    }
    sort(triples.begin(), triples.end());
    return triples;
}

static vector<int> reduce_high_divergence_candidates(
        const vector<Column> &cols, int n, const Options &opt) {
    vector<int> candidates;
    candidates.reserve(cols.size());
    for (int i = 0; i < (int)cols.size(); i++) {
        const Column &col = cols[i];
        if (col.leaves.size() <= 1) continue;
        if ((int)col.leaves.size() > opt.high_divergence_max_component_size) continue;
        if ((int)col.edges.size() > opt.high_divergence_column_max_edges) continue;
        candidates.push_back(i);
    }
    auto better = [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        double da = (double)wa / pow((double)max(1, (int)cols[a].edges.size()), 0.18);
        double db = (double)wb / pow((double)max(1, (int)cols[b].edges.size()), 0.18);
        bool large_a = cols[a].leaves.size() >= 4;
        bool large_b = cols[b].leaves.size() >= 4;
        if (large_a != large_b) return large_a > large_b;
        if (wa != wb) return wa > wb;
        if (da != db) return da > db;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    sort(candidates.begin(), candidates.end(), better);
    if (opt.high_divergence_top_per_leaf > 0) {
        unordered_set<int> chosen;
        vector<vector<int>> by_leaf(n + 1);
        for (int idx : candidates) for (int leaf : cols[idx].leaves) by_leaf[leaf].push_back(idx);
        for (int leaf = 1; leaf <= n; leaf++) {
            auto &bucket = by_leaf[leaf];
            sort(bucket.begin(), bucket.end(), better);
            int keep = min(opt.high_divergence_top_per_leaf, (int)bucket.size());
            for (int i = 0; i < keep; i++) chosen.insert(bucket[i]);
        }
        candidates.assign(chosen.begin(), chosen.end());
        sort(candidates.begin(), candidates.end(), better);
    }
    if (opt.high_divergence_sat_candidate_limit > 0 &&
            (int)candidates.size() > opt.high_divergence_sat_candidate_limit) {
        candidates.resize(opt.high_divergence_sat_candidate_limit);
    }
    return candidates;
}

static vector<int> pair_region_exact_repack(
        const vector<Column> &cols, const vector<int> &pairs,
        const vector<unsigned long long> &fixed_used,
        const vector<int> &free_leaves,
        Clock::time_point deadline, long long node_limit) {
    if (free_leaves.size() < 2 || Clock::now() >= deadline) return {};
    unordered_map<int, int> local_of;
    local_of.reserve(free_leaves.size() * 2);
    for (int i = 0; i < (int)free_leaves.size(); i++) local_of[free_leaves[i]] = i;

    struct LocalPair {
        int col = -1;
        unsigned long long leaf_bits = 0;
    };
    vector<LocalPair> cand;
    cand.reserve((free_leaves.size() * (free_leaves.size() - 1)) / 2);
    for (int p : pairs) {
        auto it0 = local_of.find(cols[p].leaves[0]);
        if (it0 == local_of.end()) continue;
        auto it1 = local_of.find(cols[p].leaves[1]);
        if (it1 == local_of.end()) continue;
        if (intersects(fixed_used, cols[p].mask)) continue;
        unsigned long long bits = (1ULL << it0->second) | (1ULL << it1->second);
        cand.push_back(LocalPair{p, bits});
    }
    sort(cand.begin(), cand.end(), [&](const LocalPair &a, const LocalPair &b) {
        if (cols[a.col].edges.size() != cols[b.col].edges.size()) {
            return cols[a.col].edges.size() < cols[b.col].edges.size();
        }
        return cols[a.col].leaves < cols[b.col].leaves;
    });
    if (cand.empty()) return {};

    int m = (int)free_leaves.size();
    vector<vector<int>> by_leaf(m);
    for (int i = 0; i < (int)cand.size(); i++) {
        for (int bit = 0; bit < m; bit++) {
            if ((cand[i].leaf_bits >> bit) & 1ULL) by_leaf[bit].push_back(i);
        }
    }

    vector<int> chosen, best;
    long long nodes = 0;
    function<void(unsigned long long, vector<unsigned long long>&)> dfs =
        [&](unsigned long long closed_leaves, vector<unsigned long long> &used) {
            if (++nodes > node_limit || Clock::now() >= deadline) return;
            int open = m - __builtin_popcountll(closed_leaves);
            if ((int)chosen.size() + open / 2 <= (int)best.size()) return;
            int pivot = -1;
            int pivot_options = INT_MAX;
            for (int leaf = 0; leaf < m; leaf++) {
                if ((closed_leaves >> leaf) & 1ULL) continue;
                int options = 0;
                for (int ci : by_leaf[leaf]) {
                    if (cand[ci].leaf_bits & closed_leaves) continue;
                    if (intersects(used, cols[cand[ci].col].mask)) continue;
                    options++;
                }
                if (options == 0) {
                    pivot = leaf;
                    pivot_options = 0;
                    break;
                }
                if (options < pivot_options) {
                    pivot = leaf;
                    pivot_options = options;
                }
            }
            if (pivot < 0) {
                if (chosen.size() > best.size()) best = chosen;
                return;
            }
            if (pivot_options > 0) {
                for (int ci : by_leaf[pivot]) {
                    const LocalPair &lp = cand[ci];
                    if (lp.leaf_bits & closed_leaves) continue;
                    if (intersects(used, cols[lp.col].mask)) continue;
                    vector<unsigned long long> next_used = used;
                    bit_or_into(next_used, cols[lp.col].mask);
                    chosen.push_back(lp.col);
                    dfs(closed_leaves | lp.leaf_bits, next_used);
                    chosen.pop_back();
                    if (nodes > node_limit || Clock::now() >= deadline) return;
                }
            }
            dfs(closed_leaves | (1ULL << pivot), used);
        };

    vector<unsigned long long> used = fixed_used;
    dfs(0, used);
    return best;
}

static vector<int> local_improve_seed_pair_selection(
        const vector<Column> &cols, const vector<int> &pairs,
        int seed, vector<int> selected_pairs, int need_pairs,
        Clock::time_point deadline) {
    if ((int)selected_pairs.size() >= need_pairs) {
        vector<int> selected;
        if (seed >= 0) selected.push_back(seed);
        selected.insert(selected.end(), selected_pairs.begin(), selected_pairs.end());
        return selected;
    }
    vector<unsigned char> pair_selected(cols.size(), 0);
    for (int p : selected_pairs) pair_selected[p] = 1;
    int best_pair_count = (int)selected_pairs.size();
    const int max_conflict_old = 9;
    int passes = 0;
    long long examined = 0;

    while (Clock::now() < deadline && passes < 20 && (int)selected_pairs.size() < need_pairs) {
        passes++;
        bool improved = false;
        for (int forced : pairs) {
            if (Clock::now() >= deadline) break;
            if (pair_selected[forced]) continue;
            if (seed >= 0 && intersects(cols[seed].mask, cols[forced].mask)) continue;
            examined++;
            vector<int> conflicts;
            conflicts.reserve(max_conflict_old + 1);
            for (int p : selected_pairs) {
                if (!intersects(cols[p].mask, cols[forced].mask)) continue;
                conflicts.push_back(p);
                if ((int)conflicts.size() > max_conflict_old) break;
            }
            if ((int)conflicts.size() > max_conflict_old) continue;

            vector<unsigned long long> fixed_used(
                cols.empty() ? 0 : cols[0].mask.size(), 0);
            if (seed >= 0) bit_or_into(fixed_used, cols[seed].mask);
            for (int p : selected_pairs) {
                if (find(conflicts.begin(), conflicts.end(), p) != conflicts.end()) continue;
                bit_or_into(fixed_used, cols[p].mask);
            }
            if (intersects(fixed_used, cols[forced].mask)) continue;
            bit_or_into(fixed_used, cols[forced].mask);

            vector<int> free_leaves;
            unordered_set<int> seen_leaf;
            seen_leaf.reserve(64);
            auto add_free_leaf = [&](int leaf) {
                if (mask_has_bit(cols[forced].leaf_mask, leaf - 1)) return;
                if (seen_leaf.insert(leaf).second) free_leaves.push_back(leaf);
            };
            for (int p : conflicts) {
                for (int leaf : cols[p].leaves) add_free_leaf(leaf);
            }

            vector<int> region_pairs;
            region_pairs.reserve(256);
            for (int q : pairs) {
                if (q == forced || pair_selected[q]) continue;
                if (intersects(fixed_used, cols[q].mask)) continue;
                bool touches_removed_region = false;
                for (int old : conflicts) {
                    if (intersects(cols[q].mask, cols[old].mask)) {
                        touches_removed_region = true;
                        break;
                    }
                }
                if (touches_removed_region) region_pairs.push_back(q);
            }
            sort(region_pairs.begin(), region_pairs.end(), [&](int a, int b) {
                if (cols[a].edges.size() != cols[b].edges.size()) {
                    return cols[a].edges.size() < cols[b].edges.size();
                }
                return cols[a].leaves < cols[b].leaves;
            });
            const int max_region_free_leaves = 26;
            for (int q : region_pairs) {
                int extra = 0;
                for (int leaf : cols[q].leaves) {
                    if (!mask_has_bit(cols[forced].leaf_mask, leaf - 1) &&
                            !seen_leaf.count(leaf)) {
                        extra++;
                    }
                }
                if ((int)free_leaves.size() + extra > max_region_free_leaves) continue;
                for (int leaf : cols[q].leaves) add_free_leaf(leaf);
            }
            long long node_limit = 25000;
            vector<int> extra = pair_region_exact_repack(
                cols, pairs, fixed_used, free_leaves, deadline, node_limit);
            int new_pair_count = 1 + (int)extra.size();
            if (new_pair_count <= (int)conflicts.size()) continue;

            vector<unsigned char> remove_flag(cols.size(), 0);
            for (int p : conflicts) remove_flag[p] = 1;
            vector<int> next_pairs;
            next_pairs.reserve(selected_pairs.size() - conflicts.size() + 1 + extra.size());
            for (int p : selected_pairs) {
                if (!remove_flag[p]) next_pairs.push_back(p);
                else pair_selected[p] = 0;
            }
            next_pairs.push_back(forced);
            pair_selected[forced] = 1;
            for (int p : extra) {
                if (!pair_selected[p]) {
                    next_pairs.push_back(p);
                    pair_selected[p] = 1;
                }
            }
            selected_pairs.swap(next_pairs);
            if ((int)selected_pairs.size() > best_pair_count) {
                best_pair_count = (int)selected_pairs.size();
                cerr << "SEEDED_PAIR_LOCAL_STEP seed="
                     << (seed >= 0 ? key_of(cols[seed].leaves) : string("-"))
                     << " pairs=" << best_pair_count
                     << " need_pairs=" << need_pairs
                     << " removed=" << conflicts.size()
                     << " added=" << new_pair_count
                     << " examined=" << examined
                     << " passes=" << passes << "\n";
            }
            improved = true;
            break;
        }
        if (!improved) break;
    }

    vector<int> selected;
    if (seed >= 0) selected.push_back(seed);
    selected.insert(selected.end(), selected_pairs.begin(), selected_pairs.end());
    if ((int)selected_pairs.size() >= need_pairs) {
        cerr << "SEEDED_PAIR_LOCAL_HIT seed="
             << (seed >= 0 ? key_of(cols[seed].leaves) : string("-"))
             << " pairs=" << selected_pairs.size()
             << " need_pairs=" << need_pairs
             << " passes=" << passes
             << " examined=" << examined << "\n";
        return selected;
    }
    return selected;
}

static vector<int> seeded_pair_local_completion(
        const vector<Column> &cols, int n, int target_components, double seconds) {
    if (cols.empty() || seconds <= 0.05) return {};
    int target_savings = n - target_components;
    if (target_savings <= 0) return {};
    auto start = Clock::now();
    Clock::time_point deadline = start + chrono::milliseconds((int)(seconds * 1000));
    vector<int> pairs, seeds;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 3 && sz <= 5) seeds.push_back(i);
    }
    if (pairs.empty()) return {};
    auto pair_better_sparse = [&](int a, int b) {
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    auto pair_better_dense = [&](int a, int b) {
        int span_a = cols[a].leaves.back() - cols[a].leaves.front();
        int span_b = cols[b].leaves.back() - cols[b].leaves.front();
        if (span_a != span_b) return span_a < span_b;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    vector<vector<int>> pair_orders;
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), pair_better_sparse);
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), pair_better_dense);
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), [&](int a, int b) {
        int mina = min(cols[a].leaves[0], cols[a].leaves[1]);
        int minb = min(cols[b].leaves[0], cols[b].leaves[1]);
        if (mina != minb) return mina < minb;
        return pair_better_sparse(a, b);
    });

    vector<int> all_seeds = seeds;
    all_seeds.push_back(-1);
    struct State {
        int seed = -1;
        int order = 0;
        int savings = 0;
        int compatible_pairs = 0;
        vector<int> selected_pairs;
    };
    vector<State> states;
    int best_savings = 0;
    for (int seed : all_seeds) {
        int seed_weight = seed >= 0 ? column_weight(cols[seed]) : 0;
        int compatible_pairs = 0;
        for (int p : pairs) if (seed < 0 || !intersects(cols[seed].mask, cols[p].mask)) compatible_pairs++;
        for (int order_id = 0; order_id < (int)pair_orders.size(); order_id++) {
            if (Clock::now() >= deadline) break;
            vector<unsigned long long> used(cols[0].mask.size(), 0);
            if (seed >= 0) bit_or_into(used, cols[seed].mask);
            vector<int> selected_pairs;
            selected_pairs.reserve(target_savings);
            int savings = seed_weight;
            for (int p : pair_orders[order_id]) {
                if (intersects(used, cols[p].mask)) continue;
                selected_pairs.push_back(p);
                bit_or_into(used, cols[p].mask);
                savings++;
            }
            best_savings = max(best_savings, savings);
            if (savings >= target_savings) {
                vector<int> selected;
                if (seed >= 0) selected.push_back(seed);
                selected.insert(selected.end(), selected_pairs.begin(), selected_pairs.end());
                return selected;
            }
            if (savings >= target_savings - 6) {
                states.push_back(State{seed, order_id, savings, compatible_pairs, std::move(selected_pairs)});
            }
        }
    }
    sort(states.begin(), states.end(), [&](const State &a, const State &b) {
        int wa = a.seed >= 0 ? column_weight(cols[a.seed]) : 0;
        int wb = b.seed >= 0 ? column_weight(cols[b.seed]) : 0;
        if (wa != wb) return wa > wb;
        int ea = a.seed >= 0 ? (int)cols[a.seed].edges.size() : 0;
        int eb = b.seed >= 0 ? (int)cols[b.seed].edges.size() : 0;
        if (ea != eb) return ea < eb;
        if (a.seed != b.seed) {
            if (a.seed < 0 || b.seed < 0) return a.seed > b.seed;
            return cols[a.seed].leaves < cols[b.seed].leaves;
        }
        if (a.savings != b.savings) return a.savings > b.savings;
        if (a.compatible_pairs != b.compatible_pairs) return a.compatible_pairs > b.compatible_pairs;
        return a.order < b.order;
    });
    if ((int)states.size() > 220) states.resize(220);

    vector<int> best_selected;
    int best_total_savings = 0;
    int attempts = 0;
    for (State &st : states) {
        if (Clock::now() >= deadline) break;
        int seed_weight = st.seed >= 0 ? column_weight(cols[st.seed]) : 0;
        int need_pairs = target_savings - seed_weight;
        if (need_pairs <= 0) return vector<int>{st.seed};
        attempts++;
        vector<int> selected = local_improve_seed_pair_selection(
            cols, pairs, st.seed, std::move(st.selected_pairs), need_pairs, deadline);
        int total_savings = 0;
        int pair_count = 0;
        for (int idx : selected) {
            total_savings += column_weight(cols[idx]);
            if ((int)cols[idx].leaves.size() == 2) pair_count++;
        }
        if (total_savings > best_total_savings) {
            best_total_savings = total_savings;
            best_selected = selected;
        }
        if (total_savings >= target_savings) {
            cerr << "SEEDED_PAIR_LOCAL_DONE hit=1 savings=" << total_savings
                 << " target_savings=" << target_savings
                 << " pairs=" << pair_count
                 << " attempts=" << attempts
                 << " seconds=" << elapsed(start) << "\n";
            return selected;
        }
    }
    cerr << "SEEDED_PAIR_LOCAL_DONE hit=0 best_savings=" << max(best_savings, best_total_savings)
         << " target_savings=" << target_savings
         << " states=" << states.size()
         << " attempts=" << attempts
         << " seconds=" << elapsed(start) << "\n";
    return best_selected;
}

static vector<int> seeded_pair_greedy_completion(
        const vector<Column> &cols, int n, int target_components) {
    if (cols.empty()) return {};
    int target_savings = n - target_components;
    if (target_savings <= 0) return {};
    vector<int> pairs;
    vector<int> seeds;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 3 && sz <= 5) seeds.push_back(i);
    }
    auto pair_better_sparse = [&](int a, int b) {
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    auto pair_better_dense = [&](int a, int b) {
        int span_a = cols[a].leaves.back() - cols[a].leaves.front();
        int span_b = cols[b].leaves.back() - cols[b].leaves.front();
        if (span_a != span_b) return span_a < span_b;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    vector<vector<int>> pair_orders;
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), pair_better_sparse);
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), pair_better_dense);
    pair_orders.push_back(pairs);
    sort(pair_orders.back().begin(), pair_orders.back().end(), [&](int a, int b) {
        int mina = min(cols[a].leaves[0], cols[a].leaves[1]);
        int minb = min(cols[b].leaves[0], cols[b].leaves[1]);
        if (mina != minb) return mina < minb;
        return pair_better_sparse(a, b);
    });

    sort(seeds.begin(), seeds.end(), [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        if (wa != wb) return wa > wb;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    });
    vector<int> all_seeds;
    all_seeds.push_back(-1);
    all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());

    int attempts = 0;
    int best_savings = 0;
    vector<int> best;
    for (int seed : all_seeds) {
        for (int order_id = 0; order_id < (int)pair_orders.size(); order_id++) {
            attempts++;
            vector<unsigned long long> used(cols[0].mask.size(), 0);
            vector<int> selected;
            int savings = 0;
            if (seed >= 0) {
                selected.push_back(seed);
                bit_or_into(used, cols[seed].mask);
                savings += column_weight(cols[seed]);
            }
            for (int pair_idx : pair_orders[order_id]) {
                if (intersects(used, cols[pair_idx].mask)) continue;
                selected.push_back(pair_idx);
                bit_or_into(used, cols[pair_idx].mask);
                savings += 1;
                if (savings >= target_savings) {
                    cerr << "SEEDED_PAIR_GREEDY_HIT savings=" << savings
                         << " components=" << component_count(cols, selected, n)
                         << " seed_size=" << (seed >= 0 ? cols[seed].leaves.size() : 0)
                         << " seed=" << (seed >= 0 ? key_of(cols[seed].leaves) : string("-"))
                         << " order=" << order_id
                         << " attempts=" << attempts << "\n";
                    return selected;
                }
            }
            if (savings > best_savings) {
                best_savings = savings;
                best = selected;
            }
        }
    }
    cerr << "SEEDED_PAIR_GREEDY_DONE best_savings=" << best_savings
         << " target_savings=" << target_savings
         << " attempts=" << attempts << "\n";
    return best;
}

static vector<int> seeded_pair_cadical_completion(
        const vector<Column> &cols, int n, int target_components, double seconds) {
#ifndef USE_CADICAL
    (void)cols; (void)n; (void)target_components; (void)seconds;
    return {};
#else
    if (cols.empty() || seconds <= 0.05) return {};
    int target_savings = n - target_components;
    auto start = Clock::now();
    vector<int> pairs, seeds;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 3 && sz <= 5) seeds.push_back(i);
    }
    sort(pairs.begin(), pairs.end(), [&](int a, int b) {
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    });
    sort(seeds.begin(), seeds.end(), [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        if (wa != wb) return wa > wb;
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    });

    vector<int> promising;
    unordered_map<int, int> seed_greedy_score;
    unordered_map<int, int> seed_pair_count;
    int best_greedy = 0;
    auto greedy_score = [&](int seed) {
        vector<unsigned long long> used(cols[0].mask.size(), 0);
        int savings = 0;
        if (seed >= 0) {
            bit_or_into(used, cols[seed].mask);
            savings += column_weight(cols[seed]);
        }
        for (int p : pairs) {
            if (intersects(used, cols[p].mask)) continue;
            bit_or_into(used, cols[p].mask);
            savings++;
        }
        return savings;
    };
    for (int seed : seeds) {
        int compatible_pairs = 0;
        for (int p : pairs) if (!intersects(cols[seed].mask, cols[p].mask)) compatible_pairs++;
        seed_pair_count[seed] = compatible_pairs;
        int score = greedy_score(seed);
        seed_greedy_score[seed] = score;
        best_greedy = max(best_greedy, score);
        if (score >= target_savings - 3 &&
                (compatible_pairs >= 5000 || score >= target_savings - 1)) {
            promising.push_back(seed);
        }
    }
    sort(promising.begin(), promising.end(), [&](int a, int b) {
        if (column_weight(cols[a]) != column_weight(cols[b])) return column_weight(cols[a]) > column_weight(cols[b]);
        int ga = seed_greedy_score[a], gb = seed_greedy_score[b];
        if (ga != gb) return ga > gb;
        if (seed_pair_count[a] != seed_pair_count[b]) return seed_pair_count[a] > seed_pair_count[b];
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    });
    promising.erase(unique(promising.begin(), promising.end()), promising.end());
    if ((int)promising.size() > 160) promising.resize(160);
    cerr << "SEEDED_PAIR_CADICAL_SEEDS promising=" << promising.size()
         << " best_greedy=" << best_greedy
         << " target_savings=" << target_savings << "\n";

    for (int seed_pos = 0; seed_pos < (int)promising.size(); seed_pos++) {
        double remaining = seconds - elapsed(start);
        if (remaining <= 0.05) break;
        int seed = promising[seed_pos];
        int need_pairs = target_savings - column_weight(cols[seed]);
        if (need_pairs <= 0) return vector<int>{seed};
        vector<int> candidates;
        candidates.reserve(pairs.size());
        for (int p : pairs) {
            if (!intersects(cols[seed].mask, cols[p].mask)) candidates.push_back(p);
        }
        if ((int)candidates.size() < need_pairs) continue;
        double seed_seconds = seed_pos < 3 ? 18.0 : 4.0;
        Clock::time_point deadline = Clock::now() +
            chrono::milliseconds((int)(min(remaining, seed_seconds) * 1000));
        CaDiCaL::Solver solver;
        solver.set("quiet", 1);
        solver.set("factor", 0);
        DeadlineTerminator terminator(deadline);
        solver.connect_terminator(&terminator);
        int var_count = (int)candidates.size();
        long long clause_count = 0;
        auto add_clause = [&](const vector<int> &lits) {
            for (int lit : lits) solver.add(lit);
            solver.add(0);
            clause_count++;
        };
        auto add_at_most_one = [&](const vector<int> &vars) {
            int m = (int)vars.size();
            if (m <= 1) return;
            vector<int> seq(m, 0);
            for (int i = 1; i < m; i++) seq[i] = ++var_count;
            add_clause({-vars[0], seq[1]});
            for (int i = 2; i < m; i++) {
                add_clause({-seq[i - 1], seq[i]});
                add_clause({-vars[i - 1], -seq[i - 1]});
                add_clause({-vars[i - 1], seq[i]});
            }
            add_clause({-vars[m - 1], -seq[m - 1]});
        };
        auto add_at_least_k = [&](const vector<int> &vars, int need) {
            int m = (int)vars.size();
            if (need <= 0) return;
            if (m < need) {
                add_clause({});
                return;
            }
            vector<vector<int>> dp(m, vector<int>(need + 1, 0));
            for (int i = 0; i < m; i++) {
                for (int s = 1; s <= min(need, i + 1); s++) dp[i][s] = ++var_count;
            }
            add_clause({-vars[0], dp[0][1]});
            add_clause({-dp[0][1], vars[0]});
            for (int i = 1; i < m; i++) {
                add_clause({-vars[i], dp[i][1]});
                add_clause({-dp[i - 1][1], dp[i][1]});
                add_clause({-dp[i][1], vars[i], dp[i - 1][1]});
                for (int s = 2; s <= min(need, i + 1); s++) {
                    int same = dp[i - 1][s];
                    int prev = dp[i - 1][s - 1];
                    if (same) add_clause({-same, dp[i][s]});
                    if (prev) add_clause({-vars[i], -prev, dp[i][s]});
                    vector<int> rev{-dp[i][s], vars[i]};
                    if (same) rev.push_back(same);
                    add_clause(rev);
                    vector<int> rev_prev{-dp[i][s]};
                    if (same) rev_prev.push_back(same);
                    if (prev) rev_prev.push_back(prev);
                    add_clause(rev_prev);
                }
            }
            add_clause({dp[m - 1][need]});
        };

        vector<vector<int>> by_bit(cols[0].mask.size() * 64);
        for (int i = 0; i < (int)candidates.size(); i++) {
            int var = i + 1;
            int idx = candidates[i];
            for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
                unsigned long long bits = cols[idx].mask[w];
                while (bits) {
                    int bit = __builtin_ctzll(bits);
                    int global_bit = (w << 6) + bit;
                    if (global_bit < (int)by_bit.size()) by_bit[global_bit].push_back(var);
                    bits &= bits - 1;
                }
            }
        }
        for (const vector<int> &bucket : by_bit) add_at_most_one(bucket);
        vector<int> vars(candidates.size());
        for (int i = 0; i < (int)candidates.size(); i++) vars[i] = i + 1;
        add_at_least_k(vars, need_pairs);
        cerr << "SEEDED_PAIR_CADICAL_SOLVE seed_pos=" << seed_pos
             << " seed=" << key_of(cols[seed].leaves)
             << " greedy=" << seed_greedy_score[seed]
             << " need_pairs=" << need_pairs
             << " candidates=" << candidates.size()
             << " vars=" << var_count
             << " clauses=" << clause_count
             << " seconds=" << elapsed(start) << "\n";
        int sat = solver.solve();
        cerr << "SEEDED_PAIR_CADICAL_STATUS result=" << sat
             << " seed_pos=" << seed_pos
             << " seconds=" << elapsed(start) << "\n";
        if (sat == 10) {
            vector<int> selected{seed};
            vector<unsigned long long> used = cols[seed].mask;
            int picked = 0;
            for (int i = 0; i < (int)candidates.size(); i++) {
                if (solver.val(i + 1) <= 0) continue;
                int idx = candidates[i];
                if (intersects(used, cols[idx].mask)) {
                    selected.clear();
                    break;
                }
                selected.push_back(idx);
                bit_or_into(used, cols[idx].mask);
                picked++;
            }
            solver.disconnect_terminator();
            if (!selected.empty() && picked >= need_pairs) {
                cerr << "SEEDED_PAIR_CADICAL_HIT savings="
                     << (column_weight(cols[seed]) + picked)
                     << " components=" << component_count(cols, selected, n)
                     << " seed=" << key_of(cols[seed].leaves) << "\n";
                return selected;
            }
        } else {
            solver.disconnect_terminator();
        }
    }
    return {};
#endif
}

static vector<int> seeded_pair_branch_completion(
        const vector<Column> &cols, int n, int target_components, double seconds) {
    if (cols.empty() || seconds <= 0.05) return {};
    int target_savings = n - target_components;
    vector<int> pairs, seeds;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 4 && sz <= 5) seeds.push_back(i);
    }
    sort(seeds.begin(), seeds.end(), [&](int a, int b) {
        if (column_weight(cols[a]) != column_weight(cols[b])) return column_weight(cols[a]) > column_weight(cols[b]);
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    });
    auto start = Clock::now();
    int attempts = 0;
    for (int seed : seeds) {
        double remaining = seconds - elapsed(start);
        if (remaining <= 0.05) break;
        int need_pairs = target_savings - column_weight(cols[seed]);
        if (need_pairs <= 0) return vector<int>{seed};
        vector<unsigned char> seed_leaf(n + 1, 0);
        for (int leaf : cols[seed].leaves) seed_leaf[leaf] = 1;
        vector<int> remaining_leaves;
        remaining_leaves.reserve(n - cols[seed].leaves.size());
        for (int leaf = 1; leaf <= n; leaf++) if (!seed_leaf[leaf]) remaining_leaves.push_back(leaf);
        int max_components = (int)remaining_leaves.size() - need_pairs;
        if (max_components < 0) continue;
        vector<Component> columns;
        columns.reserve(pairs.size());
        for (int p : pairs) {
            if (intersects(cols[seed].mask, cols[p].mask)) continue;
            columns.push_back(Component{cols[p].leaves, cols[p].edges, cols[p].score});
        }
        attempts++;
        double attempt_seconds = min(remaining, attempts <= 3 ? 8.0 : 2.0);
        vector<vector<int>> local = local_branch_cover(
            columns, remaining_leaves, max_components, attempt_seconds, 600);
        cerr << "SEEDED_PAIR_BRANCH_ATTEMPT seed=" << key_of(cols[seed].leaves)
             << " need_pairs=" << need_pairs
             << " pair_columns=" << columns.size()
             << " max_components=" << max_components
             << " local=" << local.size()
             << " attempts=" << attempts
             << " seconds=" << elapsed(start) << "\n";
        if (local.empty()) continue;
        unordered_map<string, int> pair_by_key;
        pair_by_key.reserve(pairs.size() * 2);
        for (int p : pairs) pair_by_key[key_of(cols[p].leaves)] = p;
        vector<int> selected{seed};
        int pair_count = 0;
        bool ok = true;
        vector<unsigned long long> used = cols[seed].mask;
        for (const vector<int> &block : local) {
            if (block.size() == 1) continue;
            auto it = pair_by_key.find(key_of(block));
            if (it == pair_by_key.end() || intersects(used, cols[it->second].mask)) {
                ok = false;
                break;
            }
            selected.push_back(it->second);
            bit_or_into(used, cols[it->second].mask);
            pair_count++;
        }
        if (ok && pair_count >= need_pairs) {
            cerr << "SEEDED_PAIR_BRANCH_HIT savings="
                 << (column_weight(cols[seed]) + pair_count)
                 << " components=" << component_count(cols, selected, n)
                 << " seed=" << key_of(cols[seed].leaves) << "\n";
            return selected;
        }
    }
    cerr << "SEEDED_PAIR_BRANCH_DONE attempts=" << attempts << "\n";
    return {};
}

static vector<Column> build_high_divergence_columns_for_decision(
        const Instance &inst, const vector<int> &edge_offset,
        const Options &opt, double seconds) {
    vector<Column> empty;
    if (seconds <= 0.0) return empty;
    int n = inst.leaf_count;
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(seconds * 1000));

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    vector<array<int, 3>> triples = sorted_valid_triples(valid, n);
    Pool pool;
    bool limited = false;
    auto time_expired = [&]() {
        return Clock::now() >= deadline || (int)pool.comps.size() >= opt.high_divergence_max_components;
    };
    auto add_block = [&](vector<int> block) {
        if (time_expired()) {
            limited = true;
            return false;
        }
        if ((int)block.size() > opt.high_divergence_max_component_size) return false;
        bool ok = add_component(pool, inst, valid, edge_offset, std::move(block));
        if (ok && (int)pool.comps.back().edges.size() > opt.high_divergence_column_max_edges) {
            pool.id.erase(key_of(pool.comps.back().leaves));
            pool.comps.pop_back();
            ok = false;
        }
        if ((int)pool.comps.size() >= opt.high_divergence_max_components) limited = true;
        return ok;
    };

    for (int a = 1; a <= n && !time_expired(); a++) {
        for (int b = a + 1; b <= n && !time_expired(); b++) add_block({a, b});
    }
    for (const auto &tri : triples) {
        if (time_expired()) break;
        add_block({tri[0], tri[1], tri[2]});
    }

    long long extension_nodes = 0;
    vector<int> current;
    function<void(int)> dfs_extend = [&](int next_leaf) {
        if (time_expired()) {
            limited = true;
            return;
        }
        extension_nodes++;
        if (extension_nodes % 4096 == 0 && Clock::now() >= deadline) {
            limited = true;
            return;
        }
        if ((int)current.size() >= 4) add_block(current);
        if ((int)current.size() >= opt.high_divergence_max_component_size) return;
        for (int leaf = next_leaf; leaf <= n; leaf++) {
            current.push_back(leaf);
            if (triple_compatible(current, valid, n)) dfs_extend(leaf + 1);
            current.pop_back();
            if (limited || time_expired()) {
                limited = true;
                return;
            }
        }
    };
    for (const auto &tri : triples) {
        if (time_expired()) break;
        current = {tri[0], tri[1], tri[2]};
        dfs_extend(tri[2] + 1);
    }

    cerr << "HD_POOL_DECISION_COMPONENTS generated=" << pool.comps.size()
         << " valid_triples=" << triples.size()
         << " limited=" << (limited ? 1 : 0)
         << " extension_nodes=" << extension_nodes
         << " seconds=" << elapsed(start) << "\n";
    if (pool.comps.empty() || limited) return {};
    int total_bits = n + edge_offset.back();
    return encode_columns(inst, pool, total_bits);
}

static ExactDecision high_divergence_pool_cadical_decide(
        const Instance &inst, int target_components,
        const vector<int> &edge_offset, const Options &opt, double seconds,
        vector<vector<int>> *sat_partition = nullptr) {
#ifndef USE_CADICAL
    (void)inst; (void)target_components; (void)edge_offset; (void)opt;
    (void)seconds; (void)sat_partition;
    return ExactDecision::Unknown;
#else
    if (seconds <= 0.05 || target_components < 1) return ExactDecision::Unknown;
    int target_savings = inst.leaf_count - target_components;
    if (target_savings <= 0) return ExactDecision::Yes;
    auto start = Clock::now();
    double build_seconds = min(seconds * 0.35, max(20.0, seconds - 5.0));
    vector<Column> cols = build_high_divergence_columns_for_decision(
        inst, edge_offset, opt, build_seconds);
    if (cols.empty()) {
        cerr << "HD_POOL_DECISION_STATUS result=UNKNOWN reason=no_or_limited_pool\n";
        return ExactDecision::Unknown;
    }
    vector<int> candidates;
    candidates.reserve(cols.size());
    int total_weight = 0;
    for (int i = 0; i < (int)cols.size(); i++) {
        int weight = column_weight(cols[i]);
        if (weight <= 0) continue;
        candidates.push_back(i);
        total_weight += weight;
    }
    if (total_weight < target_savings) {
        cerr << "HD_POOL_DECISION_STATUS result=NO reason=total_weight"
             << " target_savings=" << target_savings
             << " total_weight=" << total_weight << "\n";
        return ExactDecision::No;
    }
    double remaining = seconds - elapsed(start);
    if (remaining <= 0.05) {
        cerr << "HD_POOL_DECISION_STATUS result=UNKNOWN reason=no_sat_time\n";
        return ExactDecision::Unknown;
    }
    Clock::time_point deadline = Clock::now() + chrono::milliseconds((int)(remaining * 1000));
    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    DeadlineTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    int var_count = (int)candidates.size();
    long long clauses = 0;
    auto add_clause = [&](const vector<int> &lits) {
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clauses++;
    };
    auto add_at_most_one = [&](const vector<int> &vars) {
        int m = (int)vars.size();
        if (m <= 1) return;
        vector<int> seq(m, 0);
        for (int i = 1; i < m; i++) seq[i] = ++var_count;
        add_clause({-vars[0], seq[1]});
        for (int i = 2; i < m; i++) {
            add_clause({-seq[i - 1], seq[i]});
            add_clause({-vars[i - 1], -seq[i - 1]});
            add_clause({-vars[i - 1], seq[i]});
        }
        add_clause({-vars[m - 1], -seq[m - 1]});
    };
    auto add_weighted_at_least = [&](const vector<int> &vars,
            const vector<int> &weights, int need) {
        int m = (int)vars.size();
        if (need <= 0) return;
        int total = 0;
        for (int w : weights) total += w;
        if (total < need) {
            add_clause({});
            return;
        }
        vector<vector<int>> dp(m, vector<int>(need + 1, 0));
        vector<int> prefix_weight(m, 0);
        int prefix = 0;
        for (int i = 0; i < m; i++) {
            prefix += weights[i];
            prefix_weight[i] = prefix;
            for (int s = 1; s <= min(need, prefix); s++) dp[i][s] = ++var_count;
        }
        for (int i = 0; i < m; i++) {
            int w = min(weights[i], need);
            int max_s = min(need, prefix_weight[i]);
            for (int s = 1; s <= max_s; s++) {
                int cur = dp[i][s];
                int same = i > 0 ? dp[i - 1][s] : 0;
                int prev_need = max(0, s - w);
                int take_prev = (i > 0 && prev_need > 0) ? dp[i - 1][prev_need] : 0;
                if (same) add_clause({-same, cur});
                if (prev_need == 0) {
                    add_clause({-vars[i], cur});
                } else if (take_prev) {
                    add_clause({-vars[i], -take_prev, cur});
                }
                vector<int> reverse_same{-cur};
                if (same) reverse_same.push_back(same);
                reverse_same.push_back(vars[i]);
                add_clause(reverse_same);
                if (prev_need > 0) {
                    vector<int> reverse_take{-cur};
                    if (same) reverse_take.push_back(same);
                    if (take_prev) reverse_take.push_back(take_prev);
                    add_clause(reverse_take);
                }
            }
        }
        add_clause({dp[m - 1][need]});
    };

    int total_bits = inst.leaf_count + edge_offset.back();
    vector<vector<int>> by_bit(total_bits);
    vector<int> vars;
    vector<int> weights;
    vars.reserve(candidates.size());
    weights.reserve(candidates.size());
    for (int pos = 0; pos < (int)candidates.size(); pos++) {
        int var = pos + 1;
        int idx = candidates[pos];
        vars.push_back(var);
        weights.push_back(column_weight(cols[idx]));
        for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
            unsigned long long bits = cols[idx].mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int global_bit = (w << 6) + bit;
                if (global_bit < total_bits) by_bit[global_bit].push_back(var);
                bits &= bits - 1;
            }
        }
    }
    long long nnz = 0;
    for (const vector<int> &bucket : by_bit) {
        nnz += bucket.size();
        add_at_most_one(bucket);
    }
    add_weighted_at_least(vars, weights, target_savings);
    cerr << "HD_POOL_DECISION_SOLVE target=" << target_components
         << " target_savings=" << target_savings
         << " candidates=" << candidates.size()
         << " nnz=" << nnz
         << " vars=" << var_count
         << " clauses=" << clauses
         << " build_seconds=" << elapsed(start)
         << " sat_seconds=" << remaining << "\n";
    int sat = solver.solve();
    cerr << "HD_POOL_DECISION_STATUS result="
         << (sat == 20 ? "NO" : sat == 10 ? "YES" : "UNKNOWN")
         << " cadical=" << sat
         << " target=" << target_components
         << " seconds=" << elapsed(start) << "\n";
    if (sat == 20) return ExactDecision::No;
    if (sat == 10) {
        vector<int> selected;
        vector<unsigned long long> used(cols[0].mask.size(), 0);
        int savings = 0;
        for (int pos = 0; pos < (int)candidates.size(); pos++) {
            if (solver.val(pos + 1) <= 0) continue;
            int idx = candidates[pos];
            if (intersects(used, cols[idx].mask)) {
                selected.clear();
                break;
            }
            selected.push_back(idx);
            bit_or_into(used, cols[idx].mask);
            savings += column_weight(cols[idx]);
        }
        if (!selected.empty() && savings >= target_savings) {
            vector<vector<int>> partition = selected_partition(cols, selected, inst.leaf_count);
            if (validate_partition(inst, partition, edge_offset)) {
                if (sat_partition) *sat_partition = std::move(partition);
                return ExactDecision::Yes;
            }
        }
    }
    return ExactDecision::Unknown;
#endif
}

static ExactDecision fixed_seed_pair_cadical_decide(
        const vector<Column> &cols, int n, int seed, const vector<int> &pairs,
        int need_pairs, Clock::time_point deadline, vector<int> *selected_out,
        long long *vars_out = nullptr, long long *clauses_out = nullptr) {
#ifndef USE_CADICAL
    (void)cols; (void)n; (void)seed; (void)pairs; (void)need_pairs;
    (void)deadline; (void)selected_out; (void)vars_out; (void)clauses_out;
    return ExactDecision::Unknown;
#else
    if (need_pairs <= 0) {
        if (selected_out && seed >= 0) *selected_out = vector<int>{seed};
        return ExactDecision::Yes;
    }
    if (cols.empty() || Clock::now() >= deadline) return ExactDecision::Unknown;
    vector<int> candidates;
    candidates.reserve(pairs.size());
    for (int p : pairs) {
        if (seed >= 0 && intersects(cols[seed].mask, cols[p].mask)) continue;
        candidates.push_back(p);
    }
    if ((int)candidates.size() < need_pairs) return ExactDecision::No;

    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    DeadlineTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    int var_count = (int)candidates.size();
    long long clauses = 0;
    auto add_clause = [&](const vector<int> &lits) {
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clauses++;
    };
    auto add_at_most_one = [&](const vector<int> &vars) {
        int m = (int)vars.size();
        if (m <= 1) return;
        vector<int> seq(m, 0);
        for (int i = 1; i < m; i++) seq[i] = ++var_count;
        add_clause({-vars[0], seq[1]});
        for (int i = 2; i < m; i++) {
            add_clause({-seq[i - 1], seq[i]});
            add_clause({-vars[i - 1], -seq[i - 1]});
            add_clause({-vars[i - 1], seq[i]});
        }
        add_clause({-vars[m - 1], -seq[m - 1]});
    };
    auto add_at_least_k = [&](const vector<int> &vars, int need) {
        int m = (int)vars.size();
        if (need <= 0) return;
        if (m < need) {
            add_clause({});
            return;
        }
        vector<vector<int>> dp(m, vector<int>(need + 1, 0));
        for (int i = 0; i < m; i++) {
            for (int s = 1; s <= min(need, i + 1); s++) dp[i][s] = ++var_count;
        }
        add_clause({-vars[0], dp[0][1]});
        add_clause({-dp[0][1], vars[0]});
        for (int i = 1; i < m; i++) {
            add_clause({-vars[i], dp[i][1]});
            add_clause({-dp[i - 1][1], dp[i][1]});
            add_clause({-dp[i][1], vars[i], dp[i - 1][1]});
            for (int s = 2; s <= min(need, i + 1); s++) {
                int cur = dp[i][s];
                int same = dp[i - 1][s];
                int prev = dp[i - 1][s - 1];
                if (same) add_clause({-same, cur});
                if (prev) add_clause({-vars[i], -prev, cur});
                vector<int> rev{-cur, vars[i]};
                if (same) rev.push_back(same);
                add_clause(rev);
                vector<int> rev_prev{-cur};
                if (same) rev_prev.push_back(same);
                if (prev) rev_prev.push_back(prev);
                add_clause(rev_prev);
            }
        }
        add_clause({dp[m - 1][need]});
    };

    vector<vector<int>> by_bit(cols[0].mask.size() * 64);
    for (int i = 0; i < (int)candidates.size(); i++) {
        int var = i + 1;
        int idx = candidates[i];
        for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
            unsigned long long bits = cols[idx].mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int global_bit = (w << 6) + bit;
                if (global_bit < (int)by_bit.size()) by_bit[global_bit].push_back(var);
                bits &= bits - 1;
            }
        }
    }
    for (const vector<int> &bucket : by_bit) add_at_most_one(bucket);
    vector<int> vars(candidates.size());
    for (int i = 0; i < (int)candidates.size(); i++) vars[i] = i + 1;
    add_at_least_k(vars, need_pairs);
    if (vars_out) *vars_out = var_count;
    if (clauses_out) *clauses_out = clauses;
    int sat = solver.solve();
    solver.disconnect_terminator();
    if (sat == 20) return ExactDecision::No;
    if (sat == 10) {
        if (selected_out) {
            vector<int> selected;
            vector<unsigned long long> used(cols[0].mask.size(), 0);
            if (seed >= 0) {
                selected.push_back(seed);
                bit_or_into(used, cols[seed].mask);
            }
            int picked = 0;
            for (int i = 0; i < (int)candidates.size(); i++) {
                if (solver.val(i + 1) <= 0) continue;
                int idx = candidates[i];
                if (intersects(used, cols[idx].mask)) {
                    selected.clear();
                    break;
                }
                selected.push_back(idx);
                bit_or_into(used, cols[idx].mask);
                picked++;
            }
            if (!selected.empty() && picked >= need_pairs) *selected_out = std::move(selected);
        }
        return ExactDecision::Yes;
    }
    return ExactDecision::Unknown;
#endif
}

static ExactDecision seeded_pair_decomposition_decide(
        const Instance &inst, int target_components,
        const vector<int> &edge_offset, const Options &opt, double seconds,
        vector<vector<int>> *sat_partition = nullptr) {
    if (seconds <= 0.05 || target_components < 1) return ExactDecision::Unknown;
    int target_savings = inst.leaf_count - target_components;
    if (target_savings <= 0) return ExactDecision::Yes;
    auto start = Clock::now();
    double build_seconds = min(seconds * 0.35, max(20.0, seconds - 5.0));
    vector<Column> cols = build_high_divergence_columns_for_decision(
        inst, edge_offset, opt, build_seconds);
    if (cols.empty()) {
        cerr << "SEEDED_PAIR_DECISION_STATUS result=UNKNOWN reason=no_or_limited_pool\n";
        return ExactDecision::Unknown;
    }
    vector<int> pairs, seeds;
    for (int i = 0; i < (int)cols.size(); i++) {
        int sz = (int)cols[i].leaves.size();
        if (sz == 2) pairs.push_back(i);
        else if (sz >= 3 && sz <= opt.high_divergence_max_component_size) seeds.push_back(i);
    }
    auto pair_sparse = [&](int a, int b) {
        if (cols[a].edges.size() != cols[b].edges.size()) return cols[a].edges.size() < cols[b].edges.size();
        return cols[a].leaves < cols[b].leaves;
    };
    sort(pairs.begin(), pairs.end(), pair_sparse);

    struct SeedInfo {
        int seed = -1;
        int weight = 0;
        int greedy_savings = 0;
        int compatible_pairs = 0;
    };
    vector<SeedInfo> ranked;
    ranked.reserve(seeds.size() + 1);
    auto score_seed = [&](int seed) {
        vector<unsigned long long> used(cols[0].mask.size(), 0);
        int savings = 0;
        int compatible = 0;
        if (seed >= 0) {
            bit_or_into(used, cols[seed].mask);
            savings += column_weight(cols[seed]);
        }
        for (int p : pairs) {
            if (seed >= 0 && intersects(cols[seed].mask, cols[p].mask)) continue;
            compatible++;
            if (intersects(used, cols[p].mask)) continue;
            bit_or_into(used, cols[p].mask);
            savings++;
        }
        return SeedInfo{seed, seed >= 0 ? column_weight(cols[seed]) : 0, savings, compatible};
    };
    ranked.push_back(score_seed(-1));
    for (int seed : seeds) ranked.push_back(score_seed(seed));
    sort(ranked.begin(), ranked.end(), [&](const SeedInfo &a, const SeedInfo &b) {
        if (a.weight != b.weight) return a.weight > b.weight;
        int ea = a.seed >= 0 ? (int)cols[a.seed].edges.size() : 0;
        int eb = b.seed >= 0 ? (int)cols[b.seed].edges.size() : 0;
        if (ea != eb) return ea < eb;
        if (a.seed < 0 || b.seed < 0) return a.seed > b.seed;
        if (cols[a.seed].leaves != cols[b.seed].leaves) return cols[a.seed].leaves < cols[b.seed].leaves;
        if (a.greedy_savings != b.greedy_savings) return a.greedy_savings > b.greedy_savings;
        if (a.compatible_pairs != b.compatible_pairs) return a.compatible_pairs > b.compatible_pairs;
        return a.seed < b.seed;
    });
    bool seed_limited = opt.seeded_pair_decision_seed_limit > 0 &&
        (int)ranked.size() > opt.seeded_pair_decision_seed_limit;
    if (seed_limited) ranked.resize(opt.seeded_pair_decision_seed_limit);

    int solved_no = 0, unknown = 0, skipped_bound = 0;
    long long total_vars = 0, total_clauses = 0;
    for (int pos = 0; pos < (int)ranked.size(); pos++) {
        double remaining = seconds - elapsed(start);
        if (remaining <= 0.05) {
            unknown++;
            break;
        }
        const SeedInfo &info = ranked[pos];
        int need_pairs = target_savings - info.weight;
        if (need_pairs <= 0) {
            vector<int> selected;
            if (info.seed >= 0) selected.push_back(info.seed);
            vector<vector<int>> partition = selected_partition(cols, selected, inst.leaf_count);
            if (validate_partition(inst, partition, edge_offset)) {
                if (sat_partition) *sat_partition = std::move(partition);
                return ExactDecision::Yes;
            }
        }
        if (info.compatible_pairs < need_pairs) {
            skipped_bound++;
            continue;
        }
        if (info.greedy_savings >= target_savings - 8) {
            vector<unsigned long long> used(cols[0].mask.size(), 0);
            if (info.seed >= 0) bit_or_into(used, cols[info.seed].mask);
            vector<int> greedy_pairs;
            greedy_pairs.reserve(need_pairs);
            for (int p : pairs) {
                if (intersects(used, cols[p].mask)) continue;
                greedy_pairs.push_back(p);
                bit_or_into(used, cols[p].mask);
            }
            Clock::time_point local_deadline = Clock::now() +
                chrono::milliseconds((int)(min(remaining, 18.0) * 1000));
            vector<int> local_selected = local_improve_seed_pair_selection(
                cols, pairs, info.seed, std::move(greedy_pairs), need_pairs,
                local_deadline);
            int local_savings = 0;
            int local_pairs = 0;
            for (int idx : local_selected) {
                local_savings += column_weight(cols[idx]);
                if ((int)cols[idx].leaves.size() == 2) local_pairs++;
            }
            cerr << "SEEDED_PAIR_DECISION_LOCAL pos=" << pos
                 << " seed=" << (info.seed >= 0 ? key_of(cols[info.seed].leaves) : string("-"))
                 << " savings=" << local_savings
                 << " target_savings=" << target_savings
                 << " pairs=" << local_pairs
                 << " need_pairs=" << need_pairs
                 << " seconds=" << elapsed(start) << "\n";
            if (local_savings >= target_savings) {
                vector<vector<int>> partition = selected_partition(
                    cols, local_selected, inst.leaf_count);
                if (validate_partition(inst, partition, edge_offset) &&
                        (int)partition.size() <= target_components) {
                    int components = (int)partition.size();
                    if (sat_partition) *sat_partition = std::move(partition);
                    cerr << "SEEDED_PAIR_DECISION_STATUS result=YES"
                         << " components=" << components
                         << " pos=" << pos
                         << " method=local_repair\n";
                    return ExactDecision::Yes;
                }
            }
            remaining = seconds - elapsed(start);
            if (remaining <= 0.05) {
                unknown++;
                break;
            }
        }
        double per_seed = min(remaining, pos < 8 ? 20.0 : pos < 40 ? 6.0 : 2.0);
        Clock::time_point deadline = Clock::now() + chrono::milliseconds((int)(per_seed * 1000));
        vector<int> selected;
        long long vars = 0, clauses = 0;
        ExactDecision ans = fixed_seed_pair_cadical_decide(
            cols, inst.leaf_count, info.seed, pairs, need_pairs, deadline,
            &selected, &vars, &clauses);
        total_vars += vars;
        total_clauses += clauses;
        cerr << "SEEDED_PAIR_DECISION_SEED pos=" << pos
             << " seed=" << (info.seed >= 0 ? key_of(cols[info.seed].leaves) : string("-"))
             << " weight=" << info.weight
             << " greedy=" << info.greedy_savings
             << " need_pairs=" << need_pairs
             << " compatible=" << info.compatible_pairs
             << " result=" << (ans == ExactDecision::No ? "NO" :
                     ans == ExactDecision::Yes ? "YES" : "UNKNOWN")
             << " vars=" << vars
             << " clauses=" << clauses
             << " seconds=" << elapsed(start) << "\n";
        if (ans == ExactDecision::Yes && !selected.empty()) {
            vector<vector<int>> partition = selected_partition(cols, selected, inst.leaf_count);
            if (validate_partition(inst, partition, edge_offset) &&
                    (int)partition.size() <= target_components) {
                int components = (int)partition.size();
                if (sat_partition) *sat_partition = std::move(partition);
                cerr << "SEEDED_PAIR_DECISION_STATUS result=YES"
                     << " components=" << components
                     << " pos=" << pos << "\n";
                return ExactDecision::Yes;
            }
            unknown++;
        } else if (ans == ExactDecision::No) {
            solved_no++;
        } else {
            unknown++;
        }
    }
    bool exhausted = !seed_limited && unknown == 0;
    cerr << "SEEDED_PAIR_DECISION_STATUS result="
         << (exhausted ? "NO" : "UNKNOWN")
         << " target=" << target_components
         << " seeds_checked=" << ranked.size()
         << " seed_limited=" << (seed_limited ? 1 : 0)
         << " solved_no=" << solved_no
         << " skipped_bound=" << skipped_bound
         << " unknown=" << unknown
         << " total_vars=" << total_vars
         << " total_clauses=" << total_clauses
         << " seconds=" << elapsed(start) << "\n";
    return exhausted ? ExactDecision::No : ExactDecision::Unknown;
}

static vector<vector<int>> high_divergence_pack(const Instance &inst, const vector<Column> &base_cols,
        const vector<int> &edge_offset, const Options &opt) {
    if (opt.high_divergence_pack_seconds <= 0.0 || !opt.has_target) return {};
    int n = inst.leaf_count;
    int target_savings = n - opt.target_components;
    if (target_savings <= 0) return {};
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(opt.high_divergence_pack_seconds * 1000));

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    vector<array<int, 3>> triples = sorted_valid_triples(valid, n);
    Pool pool;
    bool limited = false;
    auto time_expired = [&]() {
        return Clock::now() >= deadline || (int)pool.comps.size() >= opt.high_divergence_max_components;
    };
    auto add_block = [&](vector<int> block) {
        if (time_expired()) {
            limited = true;
            return false;
        }
        if ((int)block.size() > opt.high_divergence_max_component_size) return false;
        bool ok = add_component(pool, inst, valid, edge_offset, std::move(block));
        if (ok && (int)pool.comps.back().edges.size() > opt.high_divergence_column_max_edges) {
            pool.id.erase(key_of(pool.comps.back().leaves));
            pool.comps.pop_back();
            ok = false;
        }
        if ((int)pool.comps.size() >= opt.high_divergence_max_components) limited = true;
        return ok;
    };

    for (int a = 1; a <= n && !time_expired(); a++) {
        for (int b = a + 1; b <= n && !time_expired(); b++) add_block({a, b});
    }
    for (const auto &tri : triples) {
        if (time_expired()) break;
        add_block({tri[0], tri[1], tri[2]});
    }
    for (const Column &col : base_cols) {
        if (time_expired()) break;
        if ((int)col.leaves.size() <= 1 ||
                (int)col.leaves.size() > opt.high_divergence_max_component_size ||
                (int)col.edges.size() > opt.high_divergence_column_max_edges) {
            continue;
        }
        add_block(col.leaves);
    }

    long long extension_nodes = 0;
    vector<int> current;
    function<void(int)> dfs_extend = [&](int next_leaf) {
        if (time_expired()) {
            limited = true;
            return;
        }
        extension_nodes++;
        if (extension_nodes % 4096 == 0 && Clock::now() >= deadline) {
            limited = true;
            return;
        }
        if ((int)current.size() >= 4) add_block(current);
        if ((int)current.size() >= opt.high_divergence_max_component_size) return;
        for (int leaf = next_leaf; leaf <= n; leaf++) {
            current.push_back(leaf);
            if (triple_compatible(current, valid, n)) dfs_extend(leaf + 1);
            current.pop_back();
            if (limited || time_expired()) {
                limited = true;
                return;
            }
        }
    };
    for (const auto &tri : triples) {
        if (time_expired()) break;
        current = {tri[0], tri[1], tri[2]};
        dfs_extend(tri[2] + 1);
    }

    cerr << "HIGH_DIVERGENCE_COMPONENTS generated=" << pool.comps.size()
         << " valid_triples=" << triples.size()
         << " limited=" << (limited ? 1 : 0)
         << " extension_nodes=" << extension_nodes
         << " seconds=" << elapsed(start) << "\n";
    if (pool.comps.empty() || Clock::now() >= deadline) return {};

    int total_bits = n + edge_offset.back();
    vector<Column> cols = encode_columns(inst, pool, total_bits);
    if (!opt.dump_high_divergence_lp.empty()) {
        Options dump_opt = opt;
        dump_opt.dump_packing_lp = opt.dump_high_divergence_lp;
        dump_opt.dump_packing_blocks = opt.dump_high_divergence_blocks;
        dump_packing_lp_files(cols, dump_opt);
        return {};
    }
    vector<int> candidates = reduce_high_divergence_candidates(cols, n, opt);
    if (candidates.empty()) return {};

    vector<vector<int>> best_partition;
    int best_components = n + 1;
    auto consider_selected = [&](const vector<int> &selected, const string &method) -> bool {
        if (selected.empty()) return false;
        vector<vector<int>> partition = selected_partition(cols, selected, n);
        if (!validate_partition(inst, partition, edge_offset)) return false;
        int components = (int)partition.size();
        if (components < best_components) {
            best_components = components;
            best_partition = std::move(partition);
            cerr << "HIGH_DIVERGENCE_CONSTRUCTION_BEST method=" << method
                 << " components=" << best_components
                 << " target=" << opt.target_components << "\n";
        }
        return best_components <= opt.target_components;
    };
    auto consider_partition = [&](vector<vector<int>> partition, const string &method) -> bool {
        if (partition.empty() || !validate_partition(inst, partition, edge_offset)) return false;
        int components = (int)partition.size();
        if (components < best_components) {
            best_components = components;
            best_partition = std::move(partition);
            cerr << "HIGH_DIVERGENCE_CONSTRUCTION_BEST method=" << method
                 << " components=" << best_components
                 << " target=" << opt.target_components << "\n";
        }
        return best_components <= opt.target_components;
    };

    consider_selected(anchored_pack(cols, n, opt), "anchored");
    consider_selected(
        seeded_pair_greedy_completion(cols, n, opt.target_components),
        "seeded_pair_greedy");
    double remaining_seed_pair = opt.high_divergence_pack_seconds - elapsed(start);
    if (!best_partition.empty() && best_components > opt.target_components &&
            remaining_seed_pair > 1.0) {
        if (consider_selected(
                seeded_pair_local_completion(
                    cols, n, opt.target_components, min(75.0, remaining_seed_pair * 0.45)),
                "seeded_pair_local")) {
            return best_partition;
        }
    }
    remaining_seed_pair = opt.high_divergence_pack_seconds - elapsed(start);
    if (!best_partition.empty() && best_components > opt.target_components &&
            remaining_seed_pair > 1.0) {
        if (consider_selected(
                seeded_pair_branch_completion(
                    cols, n, opt.target_components, min(35.0, remaining_seed_pair * 0.35)),
                "seeded_pair_branch")) {
            return best_partition;
        }
        remaining_seed_pair = opt.high_divergence_pack_seconds - elapsed(start);
        if (consider_selected(
                seeded_pair_cadical_completion(
                    cols, n, opt.target_components, min(45.0, remaining_seed_pair * 0.40)),
                "seeded_pair_cadical")) {
            return best_partition;
        }
    }
    if (!best_partition.empty() && best_components <= opt.target_components) {
        return best_partition;
    }
    double remaining_for_lns = opt.high_divergence_pack_seconds - elapsed(start);
    if (remaining_for_lns > 2.0) {
        Options lns_opt = opt;
        lns_opt.global_search_seconds = min(30.0, max(1.0, remaining_for_lns * 0.20));
        lns_opt.global_candidate_limit = max(lns_opt.global_candidate_limit,
            opt.high_divergence_sat_candidate_limit);
        vector<int> lns_selected = global_lns_pack(cols, n, lns_opt);
        if (consider_selected(lns_selected, "lns")) return best_partition;
    }
    double remaining_for_branch = opt.high_divergence_pack_seconds - elapsed(start);
    if (remaining_for_branch > 2.0) {
        Options branch_opt = opt;
        branch_opt.global_branch_seconds = min(75.0, max(1.0, remaining_for_branch * 0.45));
        branch_opt.global_branch_column_max_size = opt.high_divergence_max_component_size;
        branch_opt.global_branch_column_max_edges = opt.high_divergence_column_max_edges;
        branch_opt.global_branch_top_per_leaf = opt.high_divergence_top_per_leaf;
        branch_opt.global_branch_width = max(branch_opt.global_branch_width, 260);
        vector<int> branch_selected = global_branch_pack(cols, n, branch_opt);
        if (consider_selected(branch_selected, "branch")) return best_partition;
    }
    double remaining_for_repair = opt.high_divergence_pack_seconds - elapsed(start);
    if (!best_partition.empty() && best_components > opt.target_components &&
            remaining_for_repair > 2.0) {
        Options repair_opt = opt;
        repair_opt.pair_augment_seconds = min(60.0, max(1.0, remaining_for_repair * 0.50));
        repair_opt.pair_augment_attempt_seconds = min(4.0, max(0.5, repair_opt.pair_augment_seconds / 12.0));
        repair_opt.pair_augment_candidate_limit = max(repair_opt.pair_augment_candidate_limit, 20000);
        repair_opt.pair_augment_max_region_components = max(repair_opt.pair_augment_max_region_components, 24);
        repair_opt.pair_augment_column_max_size = opt.high_divergence_max_component_size;
        repair_opt.pair_augment_column_max_edges = opt.high_divergence_column_max_edges;
        repair_opt.pair_augment_top_per_leaf = max(repair_opt.pair_augment_top_per_leaf,
            opt.high_divergence_top_per_leaf);
        vector<vector<int>> repaired = pair_augment_repair(
            inst, best_partition, repair_opt, edge_offset, cols);
        if (consider_partition(std::move(repaired), "pair_augment")) return best_partition;
    }
    double remaining_for_exchange = opt.high_divergence_pack_seconds - elapsed(start);
    if (!best_partition.empty() && best_components > opt.target_components &&
            remaining_for_exchange > 2.0) {
        Options exchange_opt = opt;
        exchange_opt.exact_exchange_seconds =
            min(90.0, max(1.0, remaining_for_exchange * 0.75));
        exchange_opt.exact_exchange_column_max_size = opt.high_divergence_max_component_size;
        exchange_opt.exact_exchange_column_max_edges = opt.high_divergence_column_max_edges;
        exchange_opt.exact_exchange_top_per_leaf = opt.high_divergence_top_per_leaf;
        vector<vector<int>> exchanged = exact_region_exchange_repair(
            inst, best_partition, exchange_opt, edge_offset, cols);
        if (consider_partition(std::move(exchanged), "exact_exchange")) return best_partition;
    }

#ifdef USE_CADICAL
    double remaining = opt.high_divergence_pack_seconds - elapsed(start);
    if (remaining <= 0.1) return best_partition;
    Clock::time_point sat_deadline =
        Clock::now() + chrono::milliseconds((int)(remaining * 1000));
    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    DeadlineTerminator terminator(sat_deadline);
    solver.connect_terminator(&terminator);
    int var_count = (int)candidates.size();
    long long clause_count = 0;
    auto add_clause = [&](const vector<int> &lits) {
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clause_count++;
    };
    auto add_at_most_one = [&](const vector<int> &vars) {
        int m = (int)vars.size();
        if (m <= 1) return;
        vector<int> seq(m, 0);
        for (int i = 1; i < m; i++) seq[i] = ++var_count;
        add_clause({-vars[0], seq[1]});
        for (int i = 2; i < m; i++) {
            add_clause({-seq[i - 1], seq[i]});
            add_clause({-vars[i - 1], -seq[i - 1]});
            add_clause({-vars[i - 1], seq[i]});
        }
        add_clause({-vars[m - 1], -seq[m - 1]});
    };
    auto add_weighted_at_least = [&](const vector<int> &vars,
            const vector<int> &weights, int need) {
        int m = (int)vars.size();
        if (need <= 0) return;
        int total_weight = 0;
        for (int w : weights) total_weight += w;
        if (total_weight < need) {
            add_clause({});
            return;
        }
        vector<vector<int>> dp(m, vector<int>(need + 1, 0));
        int prefix = 0;
        for (int i = 0; i < m; i++) {
            prefix += weights[i];
            for (int s = 1; s <= min(need, prefix); s++) dp[i][s] = ++var_count;
        }
        for (int i = 0; i < m; i++) {
            int w = min(weights[i], need);
            for (int s = 1; s <= need; s++) {
                if (!dp[i][s]) continue;
                int same = i > 0 ? dp[i - 1][s] : 0;
                int prev_need = max(0, s - w);
                int take_prev = (i > 0 && prev_need > 0) ? dp[i - 1][prev_need] : 0;
                if (same) add_clause({-same, dp[i][s]});
                if (prev_need == 0) {
                    add_clause({-vars[i], dp[i][s]});
                } else if (take_prev) {
                    add_clause({-vars[i], -take_prev, dp[i][s]});
                }
                vector<int> reverse_choice{-dp[i][s]};
                if (same) reverse_choice.push_back(same);
                reverse_choice.push_back(vars[i]);
                add_clause(reverse_choice);
                if (prev_need > 0) {
                    vector<int> reverse_prev{-dp[i][s]};
                    if (same) reverse_prev.push_back(same);
                    if (take_prev) reverse_prev.push_back(take_prev);
                    add_clause(reverse_prev);
                }
            }
        }
        add_clause({dp[m - 1][need]});
    };

    vector<int> vars, weights;
    vars.reserve(candidates.size());
    weights.reserve(candidates.size());
    vector<vector<int>> by_bit(cols[0].mask.size() * 64);
    for (int i = 0; i < (int)candidates.size(); i++) {
        int var = i + 1;
        int idx = candidates[i];
        vars.push_back(var);
        weights.push_back(column_weight(cols[idx]));
        for (int w = 0; w < (int)cols[idx].mask.size(); w++) {
            unsigned long long bits = cols[idx].mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int global_bit = (w << 6) + bit;
                if (global_bit < (int)by_bit.size()) by_bit[global_bit].push_back(var);
                bits &= bits - 1;
            }
        }
    }
    for (const vector<int> &bucket : by_bit) add_at_most_one(bucket);
    add_weighted_at_least(vars, weights, target_savings);
    cerr << "HIGH_DIVERGENCE_PACK_SOLVE target_savings=" << target_savings
         << " candidates=" << candidates.size()
         << " vars=" << var_count
         << " clauses=" << clause_count
         << " build_seconds=" << elapsed(start) << "\n";
    int sat = solver.solve();
    cerr << "HIGH_DIVERGENCE_PACK_STATUS result=" << sat
         << " target_savings=" << target_savings
         << " candidates=" << candidates.size()
         << " seconds=" << elapsed(start) << "\n";
    if (sat == 10) {
        vector<int> selected;
        int savings = 0;
        vector<unsigned long long> used(cols[0].mask.size(), 0);
        for (int i = 0; i < (int)candidates.size(); i++) {
            if (solver.val(i + 1) <= 0) continue;
            int idx = candidates[i];
            if (intersects(used, cols[idx].mask)) {
                selected.clear();
                break;
            }
            bit_or_into(used, cols[idx].mask);
            selected.push_back(idx);
            savings += column_weight(cols[idx]);
        }
        if (!selected.empty() && savings >= target_savings) {
            solver.disconnect_terminator();
            cerr << "HIGH_DIVERGENCE_PACK_HIT savings=" << savings
                 << " components=" << component_count(cols, selected, n) << "\n";
            vector<vector<int>> partition = selected_partition(cols, selected, n);
            if (validate_partition(inst, partition, edge_offset)) return partition;
            cerr << "HIGH_DIVERGENCE_PACK_INVALID_HIT\n";
            return {};
        }
    } else if (sat == 20) {
        cerr << "HIGH_DIVERGENCE_PACK_UNSAT_CANDIDATES_ONLY target_savings="
             << target_savings << "\n";
    }
    solver.disconnect_terminator();
#else
    cerr << "HIGH_DIVERGENCE_PACK_SKIPPED reason=no_cadical\n";
#endif
    return best_partition;
}

static bool solve_cover_rec(const vector<Column> &cols, const vector<int> &cand, const vector<unsigned long long> &target,
        vector<unsigned long long> used, vector<unsigned long long> covered, int max_components,
        vector<int> &chosen, vector<int> &answer, int branch_width, Clock::time_point deadline) {
    if (Clock::now() > deadline) return false;
    bool done = true;
    int leaf = -1;
    for (int w = 0; w < (int)target.size(); w++) {
        unsigned long long rem = target[w] & ~covered[w];
        if (rem) { leaf = (w << 6) + __builtin_ctzll(rem); done = false; break; }
    }
    if (done) { answer = chosen; return true; }
    if ((int)chosen.size() >= max_components) return false;
    vector<int> options;
    for (int idx : cand) {
        if (((cols[idx].leaf_mask[leaf >> 6] >> (leaf & 63)) & 1ULL) == 0) continue;
        if (intersects(used, cols[idx].mask)) continue;
        bool outside = false;
        for (int w = 0; w < (int)target.size(); w++) if (cols[idx].leaf_mask[w] & ~target[w]) { outside = true; break; }
        if (outside) continue;
        if (intersects(covered, cols[idx].leaf_mask)) continue;
        options.push_back(idx);
        if ((int)options.size() >= branch_width) break;
    }
    sort(options.begin(), options.end(), [&](int a, int b) { return cols[a].leaves.size() > cols[b].leaves.size(); });
    for (int idx : options) {
        vector<unsigned long long> next_used = used, next_covered = covered;
        bit_or_into(next_used, cols[idx].mask);
        bit_or_into(next_covered, cols[idx].leaf_mask);
        chosen.push_back(idx);
        if (solve_cover_rec(cols, cand, target, next_used, next_covered, max_components, chosen, answer, branch_width, deadline)) return true;
        chosen.pop_back();
    }
    // singleton fallback
    covered[leaf >> 6] |= 1ULL << (leaf & 63);
    return solve_cover_rec(cols, cand, target, used, covered, max_components, chosen, answer, branch_width, deadline);
}

static vector<int> swap_repair(const vector<Column> &cols, vector<int> selected, int n, const Options &opt) {
    int best_count = component_count(cols, selected, n);
    auto deadline = Clock::now() + chrono::milliseconds((int)(opt.swap_seconds * 1000));
    vector<vector<int>> base = selected_partition(cols, selected, n);
    vector<int> block_col(base.size(), -1);
    unordered_map<string, int> col_by_key;
    for (int i = 0; i < (int)cols.size(); i++) col_by_key[key_of(cols[i].leaves)] = i;
    for (int i = 0; i < (int)base.size(); i++) {
        auto it = col_by_key.find(key_of(base[i]));
        if (it != col_by_key.end()) block_col[i] = it->second;
    }
    int anchor_block = -1;
    for (int i = 0; i < (int)base.size(); i++) if (base[i].size() >= 30 && (anchor_block < 0 || base[i].size() > base[anchor_block].size())) anchor_block = i;
    if (anchor_block < 0) return selected;
    vector<unsigned long long> anchor_mask(((n + cols[0].mask.size() * 64 - n + 63) >> 6), 0);
    anchor_mask.assign(cols[0].mask.size(), 0);
    for (int x : base[anchor_block]) set_bit(anchor_mask, x - 1);
    vector<int> leaf_owner(n + 1, -1);
    for (int i = 0; i < (int)base.size(); i++) for (int x : base[i]) leaf_owner[x] = i;
    vector<int> candidates(cols.size());
    for (int i = 0; i < (int)cols.size(); i++) candidates[i] = i;
    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        int oa = 0, ob = 0;
        for (int x : cols[a].leaves) if (binary_search(base[anchor_block].begin(), base[anchor_block].end(), x)) oa++;
        for (int x : cols[b].leaves) if (binary_search(base[anchor_block].begin(), base[anchor_block].end(), x)) ob++;
        if (oa != ob) return oa > ob;
        return cols[a].leaves.size() > cols[b].leaves.size();
    });
    int attempts = 0;
    for (int rep : candidates) {
        if (Clock::now() > deadline) break;
        int rs = (int)cols[rep].leaves.size();
        int as = (int)base[anchor_block].size();
        if (rs < as - opt.swap_anchor_slack || rs > as + opt.swap_extra_leaf_limit) continue;
        int overlap = 0;
        for (int x : cols[rep].leaves) if (binary_search(base[anchor_block].begin(), base[anchor_block].end(), x)) overlap++;
        if (overlap == 0) continue;
        vector<int> combo{anchor_block};
        for (int x : cols[rep].leaves) {
            if (!binary_search(base[anchor_block].begin(), base[anchor_block].end(), x)) {
                int owner = leaf_owner[x];
                if (find(combo.begin(), combo.end(), owner) == combo.end()) combo.push_back(owner);
            }
        }
        if ((int)combo.size() > opt.swap_max_old_components) continue;
        vector<unsigned long long> union_leaf(cols[0].mask.size(), 0), fixed_mask(cols[0].mask.size(), 0);
        vector<int> fixed_selected;
        for (int i = 0; i < (int)base.size(); i++) {
            bool removed = find(combo.begin(), combo.end(), i) != combo.end();
            if (removed) for (int x : base[i]) set_bit(union_leaf, x - 1);
            else if (block_col[i] >= 0) { fixed_selected.push_back(block_col[i]); bit_or_into(fixed_mask, cols[block_col[i]].mask); }
        }
        bool rep_outside = false;
        for (int x : cols[rep].leaves) if (((union_leaf[(x - 1) >> 6] >> ((x - 1) & 63)) & 1ULL) == 0) rep_outside = true;
        if (rep_outside || intersects(fixed_mask, cols[rep].mask)) continue;
        vector<unsigned long long> target = union_leaf;
        for (int x : cols[rep].leaves) target[(x - 1) >> 6] &= ~(1ULL << ((x - 1) & 63));
        vector<unsigned long long> used = fixed_mask;
        bit_or_into(used, cols[rep].mask);
        int need_improve = best_count - opt.target_components;
        int max_new = (int)combo.size() - max(1, need_improve);
        if (max_new < 1) continue;
        vector<int> cand;
        for (int i = 0; i < (int)cols.size(); i++) {
            bool subset = true;
            for (int x : cols[i].leaves) if (((target[(x - 1) >> 6] >> ((x - 1) & 63)) & 1ULL) == 0) { subset = false; break; }
            if (subset && !intersects(used, cols[i].mask)) cand.push_back(i);
        }
        sort(cand.begin(), cand.end(), [&](int a, int b) { return cols[a].leaves.size() > cols[b].leaves.size(); });
        vector<int> chosen, answer;
        attempts++;
        if (solve_cover_rec(cols, cand, target, used, vector<unsigned long long>(target.size(), 0),
                    max_new - 1, chosen, answer, opt.swap_branch_width, deadline)) {
            vector<int> next = fixed_selected;
            next.push_back(rep);
            next.insert(next.end(), answer.begin(), answer.end());
            int count = component_count(cols, next, n);
            if (count < best_count) {
                best_count = count; selected = next;
                cerr << "SWAP_BEST components=" << best_count << " old=" << combo.size()
                     << " new=" << (1 + answer.size()) << " replacement_size=" << cols[rep].leaves.size()
                     << " attempts=" << attempts << "\n";
                if (best_count <= opt.target_components) return selected;
            }
        }
        if (attempts >= opt.swap_candidate_limit) break;
    }
    cerr << "SWAP_DONE attempts=" << attempts << " best_components=" << best_count << "\n";
    return selected;
}

static vector<int> medium_swap_repair(const vector<Column> &cols, vector<int> selected, int n, const Options &opt) {
    if (cols.empty() || opt.medium_swap_seconds <= 0.0) return selected;
    int best_count = component_count(cols, selected, n);
    if (best_count <= opt.target_components) return selected;
    if (best_count - opt.target_components > 2) return selected;
    auto deadline = Clock::now() + chrono::milliseconds((int)(opt.medium_swap_seconds * 1000));
    vector<vector<int>> base = selected_partition(cols, selected, n);
    unordered_map<string, int> col_by_key;
    col_by_key.reserve(cols.size() * 2);
    for (int i = 0; i < (int)cols.size(); i++) col_by_key[key_of(cols[i].leaves)] = i;
    vector<int> block_col(base.size(), -1);
    for (int i = 0; i < (int)base.size(); i++) {
        auto it = col_by_key.find(key_of(base[i]));
        if (it != col_by_key.end()) block_col[i] = it->second;
    }
    vector<int> leaf_owner(n + 1, -1);
    for (int i = 0; i < (int)base.size(); i++) for (int x : base[i]) leaf_owner[x] = i;

    vector<int> candidates(cols.size());
    for (int i = 0; i < (int)cols.size(); i++) candidates[i] = i;
    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        int wa = column_weight(cols[a]), wb = column_weight(cols[b]);
        double sa = wa / pow((double)max(1, mask_bit_count(cols[a].mask)), 0.20);
        double sb = wb / pow((double)max(1, mask_bit_count(cols[b].mask)), 0.20);
        if (cols[a].leaves.size() != cols[b].leaves.size()) return cols[a].leaves.size() > cols[b].leaves.size();
        if (sa != sb) return sa > sb;
        return cols[a].leaves < cols[b].leaves;
    });

    int attempts = 0;
    vector<unsigned char> selected_flag(cols.size(), 0);
    for (int idx : selected) if (idx >= 0 && idx < (int)cols.size()) selected_flag[idx] = 1;
    for (int rep : candidates) {
        if (Clock::now() > deadline) break;
        if (attempts >= opt.medium_swap_candidate_limit) break;
        if (selected_flag[rep]) continue;
        int rs = (int)cols[rep].leaves.size();
        if (rs < 2 || rs > opt.medium_swap_max_replacement_size) continue;
        attempts++;

        vector<int> combo;
        vector<unsigned char> in_combo(base.size(), 0);
        auto add_combo = [&](int idx) {
            if (idx >= 0 && idx < (int)base.size() && !in_combo[idx]) {
                in_combo[idx] = 1;
                combo.push_back(idx);
            }
        };
        for (int x : cols[rep].leaves) add_combo(leaf_owner[x]);
        for (int i = 0; i < (int)base.size(); i++) {
            int bc = block_col[i];
            if (bc >= 0 && intersects(cols[rep].mask, cols[bc].mask)) add_combo(i);
        }
        if ((int)combo.size() < 2 || (int)combo.size() > opt.medium_swap_max_old_components) continue;

        vector<unsigned long long> union_leaf(cols[0].mask.size(), 0), fixed_mask(cols[0].mask.size(), 0);
        vector<int> fixed_selected;
        for (int i = 0; i < (int)base.size(); i++) {
            if (in_combo[i]) {
                for (int x : base[i]) set_bit(union_leaf, x - 1);
            } else if (block_col[i] >= 0) {
                fixed_selected.push_back(block_col[i]);
                bit_or_into(fixed_mask, cols[block_col[i]].mask);
            }
        }
        bool rep_outside = false;
        for (int x : cols[rep].leaves) {
            if (((union_leaf[(x - 1) >> 6] >> ((x - 1) & 63)) & 1ULL) == 0) {
                rep_outside = true;
                break;
            }
        }
        if (rep_outside || intersects(fixed_mask, cols[rep].mask)) continue;

        vector<unsigned long long> target = union_leaf;
        for (int x : cols[rep].leaves) target[(x - 1) >> 6] &= ~(1ULL << ((x - 1) & 63));
        vector<unsigned long long> used = fixed_mask;
        bit_or_into(used, cols[rep].mask);
        int max_new_without_rep = (int)combo.size() - 2;
        if (max_new_without_rep < 0) continue;

        vector<int> cover_candidates;
        for (int i = 0; i < (int)cols.size(); i++) {
            if (i == rep) continue;
            bool subset = true;
            for (int x : cols[i].leaves) {
                if (((target[(x - 1) >> 6] >> ((x - 1) & 63)) & 1ULL) == 0) {
                    subset = false;
                    break;
                }
            }
            if (subset && !intersects(used, cols[i].mask)) cover_candidates.push_back(i);
        }
        sort(cover_candidates.begin(), cover_candidates.end(), [&](int a, int b) {
            if (cols[a].leaves.size() != cols[b].leaves.size()) return cols[a].leaves.size() > cols[b].leaves.size();
            return cols[a].leaves < cols[b].leaves;
        });

        vector<int> chosen, answer;
        if (solve_cover_rec(cols, cover_candidates, target, used,
                    vector<unsigned long long>(target.size(), 0),
                    max_new_without_rep, chosen, answer, opt.swap_branch_width, deadline)) {
            vector<int> next = fixed_selected;
            next.push_back(rep);
            next.insert(next.end(), answer.begin(), answer.end());
            int count = component_count(cols, next, n);
            if (count < best_count) {
                best_count = count;
                selected = std::move(next);
                cerr << "MEDIUM_SWAP_BEST components=" << best_count
                     << " old=" << combo.size()
                     << " new_non_singleton=" << (1 + answer.size())
                     << " replacement_size=" << cols[rep].leaves.size()
                     << " attempts=" << attempts << "\n";
                if (best_count <= opt.target_components) return selected;
            }
        }
    }
    cerr << "MEDIUM_SWAP_DONE attempts=" << attempts << " best_components=" << best_count << "\n";
    return selected;
}

static vector<Component> reduce_repair_columns_preserve_pairs(const vector<Component> &columns,
        const unordered_map<int, int> &leaf_to_local, int local_leaf_count, int per_leaf) {
    if (per_leaf <= 0) return columns;
    vector<Component> pairs;
    vector<Component> large;
    pairs.reserve(columns.size());
    for (const Component &c : columns) {
        if (c.leaves.size() <= 2) pairs.push_back(c);
        else large.push_back(c);
    }
    vector<Component> reduced_large = reduce_local_columns(large, leaf_to_local, local_leaf_count, per_leaf);
    unordered_set<string> seen;
    vector<Component> out;
    out.reserve(pairs.size() + reduced_large.size());
    for (const Component &c : pairs) {
        if (seen.insert(key_of(c.leaves)).second) out.push_back(c);
    }
    for (const Component &c : reduced_large) {
        if (seen.insert(key_of(c.leaves)).second) out.push_back(c);
    }
    return out;
}

static vector<vector<int>> exact_local_cover_cadical(const vector<Component> &columns,
        const vector<int> &selected_leaves, int max_components, double seconds) {
    if (max_components < 1 || selected_leaves.empty() || seconds <= 0.02) return {};
    if ((int)selected_leaves.size() <= max_components) {
        vector<vector<int>> out;
        for (int leaf : selected_leaves) out.push_back(vector<int>{leaf});
        return out;
    }
#ifdef USE_CADICAL
    auto start = Clock::now();
    Clock::time_point deadline =
        start + chrono::milliseconds((int)(seconds * 1000));
    unordered_map<int, int> leaf_row;
    for (int i = 0; i < (int)selected_leaves.size(); i++) {
        leaf_row[selected_leaves[i]] = i;
    }

    struct LCol {
        vector<int> leaves;
        vector<int> leaf_rows;
        vector<int> edges;
        int var = 0;
    };
    vector<LCol> lcols;
    lcols.reserve(columns.size());
    unordered_map<int, int> edge_row;
    vector<vector<int>> by_edge;
    for (const Component &c : columns) {
        if (c.leaves.size() <= 1) continue;
        LCol lc;
        lc.leaves = c.leaves;
        bool ok = true;
        for (int leaf : c.leaves) {
            auto it = leaf_row.find(leaf);
            if (it == leaf_row.end()) {
                ok = false;
                break;
            }
            lc.leaf_rows.push_back(it->second);
        }
        if (!ok) continue;
        lc.edges = c.edges;
        lcols.push_back(std::move(lc));
    }
    if (lcols.empty()) return {};

    CaDiCaL::Solver solver;
    solver.set("quiet", 1);
    solver.set("factor", 0);
    DeadlineTerminator terminator(deadline);
    solver.connect_terminator(&terminator);
    int var_count = 0;
    long long clause_count = 0;
    auto add_clause = [&](const vector<int> &lits) {
        for (int lit : lits) solver.add(lit);
        solver.add(0);
        clause_count++;
    };
    auto add_at_most_one = [&](const vector<int> &vars) {
        int m = (int)vars.size();
        if (m <= 1) return;
        vector<int> seq(m, 0);
        for (int i = 1; i < m; i++) seq[i] = ++var_count;
        add_clause({-vars[0], seq[1]});
        for (int i = 2; i < m; i++) {
            add_clause({-seq[i - 1], seq[i]});
            add_clause({-vars[i - 1], -seq[i - 1]});
            add_clause({-vars[i - 1], seq[i]});
        }
        add_clause({-vars[m - 1], -seq[m - 1]});
    };
    auto add_at_most_k = [&](const vector<int> &vars, int max_true) {
        int m = (int)vars.size();
        if (max_true >= m) return;
        if (max_true < 0) {
            add_clause({});
            return;
        }
        if (max_true == 0) {
            for (int var : vars) add_clause({-var});
            return;
        }
        vector<vector<int>> seq(m, vector<int>(max_true + 1, 0));
        for (int i = 0; i < m; i++) {
            for (int j = 1; j <= max_true; j++) seq[i][j] = ++var_count;
        }
        add_clause({-vars[0], seq[0][1]});
        for (int i = 1; i < m; i++) {
            add_clause({-vars[i], seq[i][1]});
            add_clause({-seq[i - 1][1], seq[i][1]});
            for (int j = 2; j <= max_true; j++) {
                add_clause({-seq[i - 1][j], seq[i][j]});
                add_clause({-vars[i], -seq[i - 1][j - 1], seq[i][j]});
            }
            add_clause({-vars[i], -seq[i - 1][max_true]});
        }
    };

    vector<vector<int>> by_leaf(selected_leaves.size());
    vector<int> all_vars;
    for (LCol &lc : lcols) {
        lc.var = ++var_count;
        all_vars.push_back(lc.var);
        for (int row : lc.leaf_rows) by_leaf[row].push_back(lc.var);
        for (int edge : lc.edges) {
            auto it = edge_row.find(edge);
            if (it == edge_row.end()) {
                int row = (int)edge_row.size();
                it = edge_row.insert({edge, row}).first;
                by_edge.push_back({});
            }
            by_edge[it->second].push_back(lc.var);
        }
    }
    vector<int> singleton_var(selected_leaves.size(), 0);
    for (int i = 0; i < (int)selected_leaves.size(); i++) {
        singleton_var[i] = ++var_count;
        all_vars.push_back(singleton_var[i]);
        by_leaf[i].push_back(singleton_var[i]);
    }

    for (const vector<int> &bucket : by_leaf) {
        add_clause(bucket);
        add_at_most_one(bucket);
    }
    for (const vector<int> &bucket : by_edge) add_at_most_one(bucket);
    add_at_most_k(all_vars, max_components);

    cerr << "EXACT_LOCAL_COVER_SOLVE leaves=" << selected_leaves.size()
         << " cols=" << lcols.size()
         << " max_components=" << max_components
         << " vars=" << var_count
         << " clauses=" << clause_count
         << " build_seconds=" << elapsed(start) << "\n";
    int sat = solver.solve();
    cerr << "EXACT_LOCAL_COVER_STATUS result=" << sat
         << " leaves=" << selected_leaves.size()
         << " cols=" << lcols.size()
         << " max_components=" << max_components
         << " seconds=" << elapsed(start) << "\n";
    if (sat != 10) {
        solver.disconnect_terminator();
        return {};
    }
    vector<vector<int>> out;
    vector<unsigned char> covered(selected_leaves.size(), 0);
    for (const LCol &lc : lcols) {
        if (solver.val(lc.var) <= 0) continue;
        out.push_back(lc.leaves);
        for (int row : lc.leaf_rows) {
            if (covered[row]) {
                solver.disconnect_terminator();
                return {};
            }
            covered[row] = 1;
        }
    }
    for (int i = 0; i < (int)selected_leaves.size(); i++) {
        if (solver.val(singleton_var[i]) > 0) {
            if (covered[i]) {
                solver.disconnect_terminator();
                return {};
            }
            covered[i] = 1;
            out.push_back(vector<int>{selected_leaves[i]});
        }
    }
    solver.disconnect_terminator();
    for (unsigned char c : covered) if (!c) return {};
    if ((int)out.size() > max_components) return {};
    sort(out.begin(), out.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return out;
#else
    (void)columns;
    (void)selected_leaves;
    (void)max_components;
    (void)seconds;
    return {};
#endif
}

static vector<vector<int>> exact_region_exchange_repair(const Instance &inst,
        vector<vector<int>> partition, const Options &opt,
        const vector<int> &edge_offset, const vector<Column> &cols) {
    if (opt.exact_exchange_seconds <= 0.0 || !opt.has_target ||
            (int)partition.size() <= opt.target_components || cols.empty()) {
        return partition;
    }
    auto start = Clock::now();
    vector<vector<int>> best = partition;
    int best_count = (int)best.size();
    int attempts = 0;

    while (elapsed(start) < opt.exact_exchange_seconds &&
            best_count > opt.target_components) {
        vector<vector<int>> block_edges;
        vector<int> leaf_owner(inst.leaf_count + 1, -1);
        unordered_map<int, vector<int>> edge_to_blocks;
        unordered_set<string> current_keys;
        block_edges.reserve(best.size());
        for (int b = 0; b < (int)best.size(); b++) {
            current_keys.insert(key_of(best[b]));
            for (int leaf : best[b]) leaf_owner[leaf] = b;
            block_edges.push_back(component_edges_all(inst, best[b], edge_offset));
            for (int edge : block_edges.back()) edge_to_blocks[edge].push_back(b);
        }

        vector<vector<int>> adjacency(best.size());
        struct SeedRegion {
            vector<int> blocks;
            double score = 0.0;
        };
        vector<SeedRegion> seeds;
        seeds.reserve(min((int)cols.size(), opt.exact_exchange_seed_limit));
        for (int i = 0; i < (int)cols.size(); i++) {
            const Column &col = cols[i];
            if (col.leaves.size() <= 1) continue;
            if ((int)col.leaves.size() > opt.exact_exchange_column_max_size) continue;
            if ((int)col.edges.size() > opt.exact_exchange_column_max_edges) continue;
            if (current_keys.count(key_of(col.leaves))) continue;
            vector<unsigned char> seen(best.size(), 0);
            vector<int> touched;
            auto add_block = [&](int b) {
                if (b >= 0 && b < (int)best.size() && !seen[b]) {
                    seen[b] = 1;
                    touched.push_back(b);
                }
            };
            for (int leaf : col.leaves) add_block(leaf_owner[leaf]);
            for (int edge : col.edges) {
                auto it = edge_to_blocks.find(edge);
                if (it == edge_to_blocks.end()) continue;
                for (int b : it->second) add_block(b);
            }
            if ((int)touched.size() < 2 ||
                    (int)touched.size() > opt.exact_exchange_max_region_components) {
                continue;
            }
            int touched_leaves = 0;
            int singleton_leaves = 0;
            int small_blocks = 0;
            for (int b : touched) {
                touched_leaves += (int)best[b].size();
                if (best[b].size() <= 2) small_blocks++;
            }
            if (touched_leaves > opt.exact_exchange_max_region_leaves) continue;
            for (int leaf : col.leaves) {
                int owner = leaf_owner[leaf];
                if (owner >= 0 && best[owner].size() == 1) singleton_leaves++;
            }
            for (int a = 0; a < (int)touched.size(); a++) {
                for (int b = a + 1; b < (int)touched.size(); b++) {
                    adjacency[touched[a]].push_back(touched[b]);
                    adjacency[touched[b]].push_back(touched[a]);
                }
            }
            double density = (double)((int)col.leaves.size() - 1) /
                pow((double)max(1, (int)col.edges.size()), 0.20);
            double score = 90.0 * singleton_leaves
                + 20.0 * ((int)col.leaves.size() - 1)
                + 8.0 * small_blocks
                - 2.0 * (int)touched.size()
                - 0.05 * touched_leaves
                + density;
            sort(touched.begin(), touched.end());
            seeds.push_back(SeedRegion{std::move(touched), score});
        }
        for (vector<int> &bucket : adjacency) {
            sort(bucket.begin(), bucket.end());
            bucket.erase(unique(bucket.begin(), bucket.end()), bucket.end());
            sort(bucket.begin(), bucket.end(), [&](int a, int b) {
                if (best[a].size() != best[b].size()) return best[a].size() < best[b].size();
                return a < b;
            });
        }
        sort(seeds.begin(), seeds.end(), [](const SeedRegion &a, const SeedRegion &b) {
            if (a.score != b.score) return a.score > b.score;
            if (a.blocks.size() != b.blocks.size()) return a.blocks.size() < b.blocks.size();
            return a.blocks < b.blocks;
        });
        if (opt.exact_exchange_seed_limit > 0 &&
                (int)seeds.size() > opt.exact_exchange_seed_limit) {
            seeds.resize(opt.exact_exchange_seed_limit);
        }

        unordered_set<string> tried_regions;
        auto region_key_blocks = [](const vector<int> &region) {
            string key;
            for (int b : region) {
                key += to_string(b);
                key += ",";
            }
            return key;
        };
        auto region_leaf_count = [&](const vector<int> &region) {
            int total = 0;
            for (int b : region) total += (int)best[b].size();
            return total;
        };
        auto try_region = [&](vector<int> region, int attempt) -> vector<vector<int>> {
            sort(region.begin(), region.end());
            region.erase(unique(region.begin(), region.end()), region.end());
            if ((int)region.size() < 2 ||
                    (int)region.size() > opt.exact_exchange_max_region_components) {
                return {};
            }
            int selected_leaf_count = region_leaf_count(region);
            if (selected_leaf_count > opt.exact_exchange_max_region_leaves) return {};
            string rkey = region_key_blocks(region);
            if (!tried_regions.insert(rkey).second) return {};

            vector<unsigned char> remove(best.size(), 0);
            for (int b : region) remove[b] = 1;
            vector<vector<int>> fixed;
            unordered_set<int> fixed_edges;
            vector<int> selected_leaves;
            for (int b = 0; b < (int)best.size(); b++) {
                if (remove[b]) {
                    selected_leaves.insert(selected_leaves.end(), best[b].begin(), best[b].end());
                } else {
                    fixed.push_back(best[b]);
                    for (int edge : block_edges[b]) fixed_edges.insert(edge);
                }
            }
            sort(selected_leaves.begin(), selected_leaves.end());
            selected_leaves.erase(unique(selected_leaves.begin(), selected_leaves.end()), selected_leaves.end());
            int max_local_components = (int)region.size() - 1;
            vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
            for (int leaf : selected_leaves) selected_mark[leaf] = 1;

            vector<Component> columns;
            unordered_set<string> seen_cols;
            for (const Column &col : cols) {
                if (col.leaves.size() <= 1) continue;
                if ((int)col.leaves.size() > opt.exact_exchange_column_max_size) continue;
                if ((int)col.edges.size() > opt.exact_exchange_column_max_edges) continue;
                bool subset = true;
                for (int leaf : col.leaves) {
                    if (!selected_mark[leaf]) {
                        subset = false;
                        break;
                    }
                }
                if (!subset) continue;
                bool conflict = false;
                for (int edge : col.edges) {
                    if (fixed_edges.count(edge)) {
                        conflict = true;
                        break;
                    }
                }
                if (conflict) continue;
                string key = key_of(col.leaves);
                if (!seen_cols.insert(key).second) continue;
                double score = (double)((int)col.leaves.size() - 1) /
                    max(1, (int)col.edges.size());
                columns.push_back(Component{col.leaves, col.edges, score});
            }
            unordered_map<int, int> leaf_to_local;
            for (int i = 0; i < (int)selected_leaves.size(); i++) {
                leaf_to_local[selected_leaves[i]] = i;
            }
            columns = reduce_repair_columns_preserve_pairs(
                columns, leaf_to_local, (int)selected_leaves.size(),
                opt.exact_exchange_top_per_leaf);
            if (opt.exact_exchange_max_candidates > 0 &&
                    (int)columns.size() > opt.exact_exchange_max_candidates) {
                sort(columns.begin(), columns.end(), [](const Component &a, const Component &b) {
                    int wa = (int)a.leaves.size() - 1;
                    int wb = (int)b.leaves.size() - 1;
                    if (wa != wb) return wa > wb;
                    if (a.edges.size() != b.edges.size()) return a.edges.size() < b.edges.size();
                    return a.leaves < b.leaves;
                });
                columns.resize(opt.exact_exchange_max_candidates);
            }
            double seconds_left = opt.exact_exchange_seconds - elapsed(start);
            double attempt_seconds = min({seconds_left, opt.exact_exchange_attempt_seconds,
                max(0.05, opt.exact_exchange_seconds - elapsed(start))});
            vector<vector<int>> local = exact_local_cover_cadical(
                columns, selected_leaves, max_local_components, attempt_seconds);
            if (!local.empty() || attempt % 100 == 0) {
                cerr << "EXACT_EXCHANGE_ATTEMPT attempt=" << attempt
                     << " region=" << region.size()
                     << " leaves=" << selected_leaves.size()
                     << " columns=" << columns.size()
                     << " max_local=" << max_local_components
                     << " local=" << local.size() << "\n";
            }
            if (local.empty()) return {};
            vector<vector<int>> candidate = fixed;
            candidate.insert(candidate.end(), local.begin(), local.end());
            sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            if (!validate_partition(inst, candidate, edge_offset)) return {};
            return candidate;
        };

        bool improved = false;
        for (int s = 0; s < (int)seeds.size() &&
                elapsed(start) < opt.exact_exchange_seconds; s++) {
            vector<vector<int>> regions;
            regions.push_back(seeds[s].blocks);
            vector<unsigned char> base_in(best.size(), 0);
            for (int b : seeds[s].blocks) base_in[b] = 1;
            vector<int> neighbor_options;
            for (int b : seeds[s].blocks) {
                for (int nb : adjacency[b]) {
                    if (!base_in[nb] &&
                            find(neighbor_options.begin(), neighbor_options.end(), nb) ==
                            neighbor_options.end()) {
                        neighbor_options.push_back(nb);
                    }
                }
            }
            sort(neighbor_options.begin(), neighbor_options.end(), [&](int a, int b) {
                int da = (int)adjacency[a].size();
                int db = (int)adjacency[b].size();
                if (best[a].size() != best[b].size()) return best[a].size() < best[b].size();
                if (da != db) return da > db;
                return a < b;
            });
            int neighbor_limit = min(10, (int)neighbor_options.size());
            for (int i = 0; i < neighbor_limit; i++) {
                vector<int> region = seeds[s].blocks;
                region.push_back(neighbor_options[i]);
                if ((int)region.size() <= opt.exact_exchange_max_region_components &&
                        region_leaf_count(region) <= opt.exact_exchange_max_region_leaves) {
                    regions.push_back(std::move(region));
                }
            }
            int pair_limit = min(6, (int)neighbor_options.size());
            for (int i = 0; i < pair_limit; i++) {
                for (int j = i + 1; j < pair_limit; j++) {
                    vector<int> region = seeds[s].blocks;
                    region.push_back(neighbor_options[i]);
                    region.push_back(neighbor_options[j]);
                    if ((int)region.size() <= opt.exact_exchange_max_region_components &&
                            region_leaf_count(region) <= opt.exact_exchange_max_region_leaves) {
                        regions.push_back(std::move(region));
                    }
                }
            }
            vector<int> expanded = seeds[s].blocks;
            vector<unsigned char> in_region(best.size(), 0);
            for (int b : expanded) in_region[b] = 1;
            for (int step = 0; step < opt.exact_exchange_max_region_components; step++) {
                int best_next = -1;
                for (int b : expanded) {
                    for (int nb : adjacency[b]) {
                        if (in_region[nb]) continue;
                        if ((int)expanded.size() + 1 > opt.exact_exchange_max_region_components) continue;
                        if (region_leaf_count(expanded) + (int)best[nb].size() >
                                opt.exact_exchange_max_region_leaves) {
                            continue;
                        }
                        if (best_next < 0 || best[nb].size() < best[best_next].size() ||
                                (best[nb].size() == best[best_next].size() && nb < best_next)) {
                            best_next = nb;
                        }
                    }
                }
                if (best_next < 0) break;
                in_region[best_next] = 1;
                expanded.push_back(best_next);
                regions.push_back(expanded);
            }
            for (vector<int> region : regions) {
                attempts++;
                vector<vector<int>> candidate = try_region(std::move(region), attempts);
                if (!candidate.empty() && (int)candidate.size() < best_count) {
                    best_count = (int)candidate.size();
                    best = std::move(candidate);
                    cerr << "EXACT_EXCHANGE_BEST components=" << best_count
                         << " attempts=" << attempts
                         << " seconds=" << elapsed(start) << "\n";
                    improved = true;
                    if (best_count <= opt.target_components) return best;
                    break;
                }
                if (elapsed(start) >= opt.exact_exchange_seconds) break;
            }
            if (improved || attempts >= opt.exact_exchange_seed_limit) break;
        }
        if (!improved) break;
    }
    cerr << "EXACT_EXCHANGE_DONE attempts=" << attempts
         << " best_components=" << best_count
         << " seconds=" << elapsed(start) << "\n";
    return best;
}

static vector<vector<int>> pair_augment_repair(const Instance &inst, vector<vector<int>> partition,
        const Options &opt, const vector<int> &edge_offset, const vector<Column> &cols) {
    if (opt.pair_augment_seconds <= 0.0 || !opt.has_target ||
            (int)partition.size() <= opt.target_components || cols.empty()) {
        return partition;
    }
    auto start = Clock::now();
    int best_count = (int)partition.size();
    vector<vector<int>> best = partition;
    vector<vector<int>> block_edges;
    block_edges.reserve(partition.size());
    for (const vector<int> &block : partition) {
        block_edges.push_back(component_edges_all(inst, block, edge_offset));
    }

    unordered_set<string> current_keys;
    vector<int> leaf_owner(inst.leaf_count + 1, -1);
    for (int i = 0; i < (int)partition.size(); i++) {
        current_keys.insert(key_of(partition[i]));
        for (int leaf : partition[i]) leaf_owner[leaf] = i;
    }

    struct AugmentSeed {
        int col = -1;
        double score = 0.0;
        int touched = 0;
    };
    vector<AugmentSeed> seeds;
    seeds.reserve(cols.size());
    for (int i = 0; i < (int)cols.size(); i++) {
        const Column &col = cols[i];
        if (col.leaves.size() < 2) continue;
        if ((int)col.leaves.size() > opt.pair_augment_column_max_size) continue;
        if ((int)col.edges.size() > opt.pair_augment_column_max_edges) continue;
        if (current_keys.count(key_of(col.leaves))) continue;
        unordered_set<int> touched;
        for (int leaf : col.leaves) {
            if (leaf_owner[leaf] >= 0) touched.insert(leaf_owner[leaf]);
        }
        for (int b = 0; b < (int)partition.size(); b++) {
            bool conflict = false;
            int p1 = 0, p2 = 0;
            while (p1 < (int)col.edges.size() && p2 < (int)block_edges[b].size()) {
                if (col.edges[p1] == block_edges[b][p2]) { conflict = true; break; }
                if (col.edges[p1] < block_edges[b][p2]) p1++; else p2++;
            }
            if (conflict) touched.insert(b);
        }
        if ((int)touched.size() > opt.pair_augment_max_region_components) continue;
        int singleton_leaves = 0;
        int small_owner_bonus = 0;
        for (int leaf : col.leaves) {
            int owner = leaf_owner[leaf];
            if (owner >= 0 && partition[owner].size() == 1) singleton_leaves++;
            if (owner >= 0 && partition[owner].size() <= 3) small_owner_bonus++;
        }
        double density = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
        double score = 30.0 * ((int)col.leaves.size() - 1)
            + 60.0 * singleton_leaves
            + 8.0 * small_owner_bonus
            - 4.0 * (int)touched.size()
            + density;
        if (col.leaves.size() == 2 && singleton_leaves == 2) score += 1000.0;
        seeds.push_back(AugmentSeed{i, score, (int)touched.size()});
    }
    sort(seeds.begin(), seeds.end(), [&](const AugmentSeed &a, const AugmentSeed &b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.touched != b.touched) return a.touched < b.touched;
        if (cols[a.col].leaves.size() != cols[b.col].leaves.size()) {
            return cols[a.col].leaves.size() > cols[b.col].leaves.size();
        }
        return cols[a.col].leaves < cols[b.col].leaves;
    });
    if (opt.pair_augment_candidate_limit > 0 &&
            (int)seeds.size() > opt.pair_augment_candidate_limit) {
        seeds.resize(opt.pair_augment_candidate_limit);
    }

    auto compatible_with_forced = [&](const vector<int> &forced, int candidate) {
        for (int idx : forced) {
            for (int leaf : cols[candidate].leaves) {
                if (binary_search(cols[idx].leaves.begin(), cols[idx].leaves.end(), leaf)) return false;
            }
            int p1 = 0, p2 = 0;
            while (p1 < (int)cols[candidate].edges.size() && p2 < (int)cols[idx].edges.size()) {
                if (cols[candidate].edges[p1] == cols[idx].edges[p2]) return false;
                if (cols[candidate].edges[p1] < cols[idx].edges[p2]) p1++; else p2++;
            }
        }
        return true;
    };

    auto expand_small_boundary = [&](vector<unsigned char> &remove, int &removed_count) {
        for (int round = 0; round < 2; round++) {
            vector<unsigned char> selected_leaf(inst.leaf_count + 1, 0);
            for (int b = 0; b < (int)partition.size(); b++) {
                if (!remove[b]) continue;
                for (int leaf : partition[b]) selected_leaf[leaf] = 1;
            }
            bool changed = false;
            for (const Column &col : cols) {
                if (col.leaves.size() <= 1) continue;
                if ((int)col.leaves.size() > opt.pair_augment_column_max_size) continue;
                if ((int)col.edges.size() > opt.pair_augment_column_max_edges) continue;
                bool has_inside = false, has_outside = false;
                vector<int> add_blocks;
                for (int leaf : col.leaves) {
                    int owner = leaf_owner[leaf];
                    if (owner < 0) { has_outside = true; continue; }
                    if (selected_leaf[leaf]) {
                        has_inside = true;
                    } else {
                        has_outside = true;
                        if (remove[owner]) continue;
                        if (partition[owner].size() > 2) {
                            add_blocks.clear();
                            add_blocks.push_back(-1);
                            break;
                        }
                        if (find(add_blocks.begin(), add_blocks.end(), owner) == add_blocks.end()) {
                            add_blocks.push_back(owner);
                        }
                    }
                }
                if (!has_inside || !has_outside || add_blocks.empty()) continue;
                if (add_blocks[0] < 0) continue;
                if (removed_count + (int)add_blocks.size() > opt.pair_augment_max_region_components) continue;
                for (int owner : add_blocks) {
                    if (!remove[owner]) {
                        remove[owner] = 1;
                        removed_count++;
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
    };

    auto try_forced = [&](const vector<int> &forced_cols, int attempt, double seconds_left)
            -> vector<vector<int>> {
        if (forced_cols.empty() || seconds_left <= 0.05) return {};
        unordered_set<int> forced_edges;
        vector<unsigned char> forced_leaf(inst.leaf_count + 1, 0);
        vector<vector<int>> forced_blocks;
        forced_blocks.reserve(forced_cols.size());
        for (int idx : forced_cols) {
            forced_blocks.push_back(cols[idx].leaves);
            for (int leaf : cols[idx].leaves) {
                if (forced_leaf[leaf]) return {};
                forced_leaf[leaf] = 1;
            }
            for (int edge : cols[idx].edges) {
                if (!forced_edges.insert(edge).second) return {};
            }
        }

        vector<unsigned char> remove(partition.size(), 0);
        int removed_count = 0;
        for (int b = 0; b < (int)partition.size(); b++) {
            bool hit = false;
            for (int leaf : partition[b]) {
                if (forced_leaf[leaf]) { hit = true; break; }
            }
            if (!hit) {
                for (int edge : block_edges[b]) {
                    if (forced_edges.count(edge)) { hit = true; break; }
                }
            }
            if (hit) {
                remove[b] = 1;
                removed_count++;
            }
        }
        expand_small_boundary(remove, removed_count);
        if (removed_count <= (int)forced_cols.size()) return {};
        if (removed_count > opt.pair_augment_max_region_components) return {};

        vector<vector<int>> fixed;
        unordered_set<int> fixed_edges;
        vector<int> selected_leaves;
        for (int b = 0; b < (int)partition.size(); b++) {
            if (remove[b]) {
                selected_leaves.insert(selected_leaves.end(), partition[b].begin(), partition[b].end());
            } else {
                fixed.push_back(partition[b]);
                for (int edge : block_edges[b]) fixed_edges.insert(edge);
            }
        }
        for (int edge : forced_edges) fixed_edges.insert(edge);
        sort(selected_leaves.begin(), selected_leaves.end());
        selected_leaves.erase(unique(selected_leaves.begin(), selected_leaves.end()), selected_leaves.end());
        vector<int> remaining_leaves;
        for (int leaf : selected_leaves) if (!forced_leaf[leaf]) remaining_leaves.push_back(leaf);

        int max_remaining_components = opt.target_components - (int)fixed.size() - (int)forced_blocks.size();
        if (max_remaining_components < 0) return {};
        vector<vector<int>> local;
        if ((int)remaining_leaves.size() <= max_remaining_components) {
            for (int leaf : remaining_leaves) local.push_back(vector<int>{leaf});
        } else {
            vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
            for (int leaf : remaining_leaves) selected_mark[leaf] = 1;
            vector<Component> columns;
            unordered_set<string> seen;
            for (const Column &col : cols) {
                if (col.leaves.size() <= 1) continue;
                if ((int)col.leaves.size() > opt.pair_augment_column_max_size) continue;
                if ((int)col.edges.size() > opt.pair_augment_column_max_edges) continue;
                bool subset = true;
                for (int leaf : col.leaves) {
                    if (!selected_mark[leaf]) { subset = false; break; }
                }
                if (!subset) continue;
                bool conflict = false;
                for (int edge : col.edges) {
                    if (fixed_edges.count(edge)) { conflict = true; break; }
                }
                if (conflict) continue;
                string key = key_of(col.leaves);
                if (!seen.insert(key).second) continue;
                double score = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
                columns.push_back(Component{col.leaves, col.edges, score});
            }
            unordered_map<int, int> leaf_to_local;
            for (int i = 0; i < (int)remaining_leaves.size(); i++) leaf_to_local[remaining_leaves[i]] = i;
            columns = reduce_repair_columns_preserve_pairs(
                columns, leaf_to_local, (int)remaining_leaves.size(), opt.pair_augment_top_per_leaf);
            double attempt_seconds = min({seconds_left, opt.pair_augment_attempt_seconds, max(0.05, opt.pair_augment_seconds - elapsed(start))});
            local = solve_local_highs(columns, remaining_leaves, max_remaining_components, attempt_seconds);
            if (local.empty()) {
                local = local_branch_cover(
                    columns, remaining_leaves, max_remaining_components,
                    attempt_seconds, opt.global_branch_width);
            }
            if (!local.empty() || attempt % 1000 == 0) {
                cerr << "PAIR_AUGMENT_ATTEMPT attempt=" << attempt
                     << " forced=" << forced_cols.size()
                     << " removed=" << removed_count
                     << " remaining_leaves=" << remaining_leaves.size()
                     << " columns=" << columns.size()
                     << " max_remaining=" << max_remaining_components
                     << " local=" << local.size() << "\n";
            }
            if (local.empty()) return {};
        }

        vector<vector<int>> candidate = fixed;
        candidate.insert(candidate.end(), forced_blocks.begin(), forced_blocks.end());
        candidate.insert(candidate.end(), local.begin(), local.end());
        sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        if (!validate_partition(inst, candidate, edge_offset)) return {};
        return candidate;
    };

    auto try_window = [&](const vector<int> &seed_cols, int attempt, double seconds_left)
            -> vector<vector<int>> {
        if (seed_cols.empty() || seconds_left <= 0.05) return {};
        unordered_set<int> seed_edges;
        vector<unsigned char> seed_leaf(inst.leaf_count + 1, 0);
        for (int idx : seed_cols) {
            for (int leaf : cols[idx].leaves) seed_leaf[leaf] = 1;
            for (int edge : cols[idx].edges) seed_edges.insert(edge);
        }

        vector<unsigned char> remove(partition.size(), 0);
        int removed_count = 0;
        for (int b = 0; b < (int)partition.size(); b++) {
            bool hit = false;
            for (int leaf : partition[b]) {
                if (seed_leaf[leaf]) { hit = true; break; }
            }
            if (!hit) {
                for (int edge : block_edges[b]) {
                    if (seed_edges.count(edge)) { hit = true; break; }
                }
            }
            if (hit) {
                remove[b] = 1;
                removed_count++;
            }
        }
        expand_small_boundary(remove, removed_count);
        int needed_drop = max(1, best_count - opt.target_components);
        if (removed_count <= needed_drop) return {};
        if (removed_count > opt.pair_augment_max_region_components) return {};

        vector<vector<int>> fixed;
        unordered_set<int> fixed_edges;
        vector<int> selected_leaves;
        for (int b = 0; b < (int)partition.size(); b++) {
            if (remove[b]) {
                selected_leaves.insert(selected_leaves.end(), partition[b].begin(), partition[b].end());
            } else {
                fixed.push_back(partition[b]);
                for (int edge : block_edges[b]) fixed_edges.insert(edge);
            }
        }
        sort(selected_leaves.begin(), selected_leaves.end());
        selected_leaves.erase(unique(selected_leaves.begin(), selected_leaves.end()), selected_leaves.end());
        int max_local_components = opt.target_components - (int)fixed.size();
        if (max_local_components < 1) return {};

        vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
        for (int leaf : selected_leaves) selected_mark[leaf] = 1;
        vector<Component> columns;
        unordered_set<string> seen;
        for (const Column &col : cols) {
            if (col.leaves.size() <= 1) continue;
            if ((int)col.leaves.size() > opt.pair_augment_column_max_size) continue;
            if ((int)col.edges.size() > opt.pair_augment_column_max_edges) continue;
            bool subset = true;
            for (int leaf : col.leaves) {
                if (!selected_mark[leaf]) { subset = false; break; }
            }
            if (!subset) continue;
            bool conflict = false;
            for (int edge : col.edges) {
                if (fixed_edges.count(edge)) { conflict = true; break; }
            }
            if (conflict) continue;
            string key = key_of(col.leaves);
            if (!seen.insert(key).second) continue;
            double score = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
            columns.push_back(Component{col.leaves, col.edges, score});
        }
        unordered_map<int, int> leaf_to_local;
        for (int i = 0; i < (int)selected_leaves.size(); i++) leaf_to_local[selected_leaves[i]] = i;
        columns = reduce_repair_columns_preserve_pairs(
            columns, leaf_to_local, (int)selected_leaves.size(), opt.pair_augment_top_per_leaf);
        double attempt_seconds = min({seconds_left, opt.pair_augment_attempt_seconds,
            max(0.05, opt.pair_augment_seconds - elapsed(start))});
        vector<vector<int>> local = solve_local_highs(columns, selected_leaves, max_local_components, attempt_seconds);
        if (local.empty()) {
            local = local_branch_cover(
                columns, selected_leaves, max_local_components,
                attempt_seconds, opt.global_branch_width);
        }
        if (!local.empty() || attempt % 500 == 0) {
            cerr << "PAIR_WINDOW_ATTEMPT attempt=" << attempt
                 << " seeds=" << seed_cols.size()
                 << " removed=" << removed_count
                 << " leaves=" << selected_leaves.size()
                 << " columns=" << columns.size()
                 << " max_local=" << max_local_components
                 << " local=" << local.size() << "\n";
        }
        if (local.empty()) return {};

        vector<vector<int>> candidate = fixed;
        candidate.insert(candidate.end(), local.begin(), local.end());
        sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        if (!validate_partition(inst, candidate, edge_offset)) return {};
        return candidate;
    };

    int attempts = 0;
    for (int s = 0; s < (int)seeds.size() && elapsed(start) < opt.pair_augment_seconds; s++) {
        vector<int> forced{seeds[s].col};
        attempts++;
        vector<vector<int>> candidate = try_forced(
            forced, attempts, opt.pair_augment_seconds - elapsed(start));
        if (!candidate.empty() && (int)candidate.size() < best_count) {
            best_count = (int)candidate.size();
            best = std::move(candidate);
            cerr << "PAIR_AUGMENT_BEST components=" << best_count
                 << " forced=" << forced.size()
                 << " attempts=" << attempts << "\n";
            if (best_count <= opt.target_components) return best;
        }
        candidate = try_window(
            forced, attempts, opt.pair_augment_seconds - elapsed(start));
        if (!candidate.empty() && (int)candidate.size() < best_count) {
            best_count = (int)candidate.size();
            best = std::move(candidate);
            cerr << "PAIR_WINDOW_BEST components=" << best_count
                 << " seeds=" << forced.size()
                 << " attempts=" << attempts << "\n";
            if (best_count <= opt.target_components) return best;
        }

        if (opt.pair_augment_bundle_size <= 1) continue;
        forced.clear();
        forced.push_back(seeds[s].col);
        for (int j = 0; j < (int)seeds.size() &&
                (int)forced.size() < opt.pair_augment_bundle_size; j++) {
            if (seeds[j].col == seeds[s].col) continue;
            if (!compatible_with_forced(forced, seeds[j].col)) continue;
            forced.push_back(seeds[j].col);
            attempts++;
            candidate = try_forced(
                forced, attempts, opt.pair_augment_seconds - elapsed(start));
            if (!candidate.empty() && (int)candidate.size() < best_count) {
                best_count = (int)candidate.size();
                best = std::move(candidate);
                cerr << "PAIR_AUGMENT_BEST components=" << best_count
                     << " forced=" << forced.size()
                     << " attempts=" << attempts << "\n";
                if (best_count <= opt.target_components) return best;
            }
            candidate = try_window(
                forced, attempts, opt.pair_augment_seconds - elapsed(start));
            if (!candidate.empty() && (int)candidate.size() < best_count) {
                best_count = (int)candidate.size();
                best = std::move(candidate);
                cerr << "PAIR_WINDOW_BEST components=" << best_count
                     << " seeds=" << forced.size()
                     << " attempts=" << attempts << "\n";
                if (best_count <= opt.target_components) return best;
            }
        }
    }
    cerr << "PAIR_AUGMENT_DONE attempts=" << attempts
         << " seeds=" << seeds.size()
         << " best_components=" << best_count << "\n";
    return best;
}

static vector<vector<int>> incumbent_guided_repack(const Instance &inst, vector<vector<int>> partition,
        const Options &opt, const vector<int> &edge_offset, const vector<Column> &cols) {
    if (opt.incumbent_repack_seconds <= 0.0 || !opt.has_target ||
            (int)partition.size() <= opt.target_components || cols.empty()) {
        return partition;
    }
    if (inst.tree_count < 3 || inst.tree_count > 8 || inst.leaf_count < 120) return partition;

    auto start = Clock::now();
    vector<vector<int>> best = partition;
    int best_count = (int)best.size();
    int rounds = 0;
    int attempts_total = 0;

    while (elapsed(start) < opt.incumbent_repack_seconds &&
            best_count > opt.target_components) {
        rounds++;
        vector<vector<int>> block_edges;
        block_edges.reserve(best.size());
        unordered_map<int, vector<int>> edge_to_blocks;
        vector<int> leaf_owner(inst.leaf_count + 1, -1);
        unordered_set<string> current_keys;
        for (int b = 0; b < (int)best.size(); b++) {
            current_keys.insert(key_of(best[b]));
            for (int leaf : best[b]) leaf_owner[leaf] = b;
            block_edges.push_back(component_edges_all(inst, best[b], edge_offset));
            for (int edge : block_edges.back()) edge_to_blocks[edge].push_back(b);
        }

        struct RepackSeed {
            int col = -1;
            double score = 0.0;
            vector<int> touched;
        };
        vector<RepackSeed> seeds;
        seeds.reserve(min((int)cols.size(), opt.incumbent_repack_candidate_limit));
        for (int i = 0; i < (int)cols.size(); i++) {
            const Column &col = cols[i];
            if (col.leaves.size() < 2) continue;
            if ((int)col.leaves.size() > opt.incumbent_repack_column_max_size) continue;
            if ((int)col.edges.size() > opt.incumbent_repack_column_max_edges) continue;
            if (current_keys.count(key_of(col.leaves))) continue;

            vector<unsigned char> seen_block(best.size(), 0);
            vector<int> touched;
            auto add_block = [&](int b) {
                if (b >= 0 && b < (int)best.size() && !seen_block[b]) {
                    seen_block[b] = 1;
                    touched.push_back(b);
                }
            };
            for (int leaf : col.leaves) add_block(leaf_owner[leaf]);
            for (int edge : col.edges) {
                auto it = edge_to_blocks.find(edge);
                if (it == edge_to_blocks.end()) continue;
                for (int b : it->second) add_block(b);
            }
            if ((int)touched.size() < 2 ||
                    (int)touched.size() > opt.incumbent_repack_max_region_components) {
                continue;
            }

            int singleton_leaves = 0;
            int small_owner_bonus = 0;
            int touched_leaves = 0;
            for (int b : touched) touched_leaves += (int)best[b].size();
            for (int leaf : col.leaves) {
                int owner = leaf_owner[leaf];
                if (owner >= 0 && best[owner].size() == 1) singleton_leaves++;
                if (owner >= 0 && best[owner].size() <= 3) small_owner_bonus++;
            }
            double density = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
            double score = 80.0 * singleton_leaves
                + 22.0 * ((int)col.leaves.size() - 1)
                + 9.0 * small_owner_bonus
                - 3.0 * (int)touched.size()
                - 0.08 * touched_leaves
                + density;
            seeds.push_back(RepackSeed{i, score, std::move(touched)});
        }
        sort(seeds.begin(), seeds.end(), [&](const RepackSeed &a, const RepackSeed &b) {
            if (a.score != b.score) return a.score > b.score;
            if (a.touched.size() != b.touched.size()) return a.touched.size() < b.touched.size();
            if (cols[a.col].leaves.size() != cols[b.col].leaves.size()) {
                return cols[a.col].leaves.size() > cols[b.col].leaves.size();
            }
            return cols[a.col].leaves < cols[b.col].leaves;
        });
        if (opt.incumbent_repack_candidate_limit > 0 &&
                (int)seeds.size() > opt.incumbent_repack_candidate_limit) {
            seeds.resize(opt.incumbent_repack_candidate_limit);
        }

        auto expand_region = [&](const vector<int> &base_region) {
            vector<int> region = base_region;
            vector<unsigned char> in_region(best.size(), 0);
            vector<unsigned char> selected_leaf(inst.leaf_count + 1, 0);
            for (int b : region) {
                in_region[b] = 1;
                for (int leaf : best[b]) selected_leaf[leaf] = 1;
            }
            int scan_limit = min((int)seeds.size(), opt.incumbent_repack_candidate_limit);
            for (int pass = 0; pass < 2; pass++) {
                bool changed = false;
                for (int s = 0; s < scan_limit &&
                        (int)region.size() < opt.incumbent_repack_max_region_components; s++) {
                    const Column &col = cols[seeds[s].col];
                    bool has_inside = false;
                    vector<int> add_blocks;
                    auto consider_add = [&](int b) {
                        if (b < 0 || b >= (int)best.size() || in_region[b]) return true;
                        if (best[b].size() > 2) return false;
                        if (find(add_blocks.begin(), add_blocks.end(), b) == add_blocks.end()) {
                            add_blocks.push_back(b);
                        }
                        return true;
                    };
                    bool ok = true;
                    for (int leaf : col.leaves) {
                        int owner = leaf_owner[leaf];
                        if (owner >= 0 && selected_leaf[leaf]) has_inside = true;
                        else if (!consider_add(owner)) { ok = false; break; }
                    }
                    if (!ok || !has_inside || add_blocks.empty()) continue;
                    for (int edge : col.edges) {
                        auto it = edge_to_blocks.find(edge);
                        if (it == edge_to_blocks.end()) continue;
                        for (int b : it->second) {
                            if (in_region[b]) has_inside = true;
                            else if (!consider_add(b)) { ok = false; break; }
                        }
                        if (!ok) break;
                    }
                    if (!ok || (int)region.size() + (int)add_blocks.size() >
                            opt.incumbent_repack_max_region_components) {
                        continue;
                    }
                    for (int b : add_blocks) {
                        in_region[b] = 1;
                        region.push_back(b);
                        for (int leaf : best[b]) selected_leaf[leaf] = 1;
                        changed = true;
                    }
                }
                if (!changed) break;
            }
            sort(region.begin(), region.end());
            region.erase(unique(region.begin(), region.end()), region.end());
            return region;
        };

        auto try_region = [&](const vector<int> &region, int seed_col, int attempt,
                double seconds_left) -> vector<vector<int>> {
            if (seconds_left <= 0.05) return {};
            if ((int)region.size() < 2 ||
                    (int)region.size() > opt.incumbent_repack_max_region_components) {
                return {};
            }
            vector<unsigned char> remove(best.size(), 0);
            for (int b : region) remove[b] = 1;

            vector<vector<int>> fixed;
            unordered_set<int> fixed_edges;
            vector<int> selected_leaves;
            for (int b = 0; b < (int)best.size(); b++) {
                if (remove[b]) {
                    selected_leaves.insert(selected_leaves.end(), best[b].begin(), best[b].end());
                } else {
                    fixed.push_back(best[b]);
                    for (int edge : block_edges[b]) fixed_edges.insert(edge);
                }
            }
            sort(selected_leaves.begin(), selected_leaves.end());
            selected_leaves.erase(unique(selected_leaves.begin(), selected_leaves.end()), selected_leaves.end());
            int max_local_components = (int)region.size() - 1;
            if (max_local_components < 1) return {};

            vector<unsigned char> selected_mark(inst.leaf_count + 1, 0);
            for (int leaf : selected_leaves) selected_mark[leaf] = 1;
            vector<Component> columns;
            unordered_set<string> seen;
            auto add_col = [&](const Column &col) {
                if (col.leaves.size() <= 1) return;
                if ((int)col.leaves.size() > opt.incumbent_repack_column_max_size) return;
                if ((int)col.edges.size() > opt.incumbent_repack_column_max_edges) return;
                for (int leaf : col.leaves) if (!selected_mark[leaf]) return;
                for (int edge : col.edges) if (fixed_edges.count(edge)) return;
                string key = key_of(col.leaves);
                if (!seen.insert(key).second) return;
                double score = (double)((int)col.leaves.size() - 1) / max(1, (int)col.edges.size());
                columns.push_back(Component{col.leaves, col.edges, score});
            };
            if (seed_col >= 0) add_col(cols[seed_col]);
            for (const Column &col : cols) add_col(col);
            unordered_map<int, int> leaf_to_local;
            for (int i = 0; i < (int)selected_leaves.size(); i++) leaf_to_local[selected_leaves[i]] = i;
            columns = reduce_repair_columns_preserve_pairs(
                columns, leaf_to_local, (int)selected_leaves.size(), opt.incumbent_repack_top_per_leaf);

            double attempt_seconds = min({seconds_left, opt.incumbent_repack_attempt_seconds,
                max(0.05, opt.incumbent_repack_seconds - elapsed(start))});
            vector<vector<int>> local = local_branch_cover(
                columns, selected_leaves, max_local_components,
                attempt_seconds, opt.incumbent_repack_branch_width);
            if (!local.empty() || attempt % 250 == 0) {
                cerr << "INCUMBENT_REPACK_ATTEMPT round=" << rounds
                     << " attempt=" << attempt
                     << " region=" << region.size()
                     << " leaves=" << selected_leaves.size()
                     << " columns=" << columns.size()
                     << " max_local=" << max_local_components
                     << " local=" << local.size() << "\n";
            }
            if (local.empty()) return {};

            vector<vector<int>> candidate = fixed;
            candidate.insert(candidate.end(), local.begin(), local.end());
            sort(candidate.begin(), candidate.end(), [](const vector<int> &a, const vector<int> &b) {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            if ((int)candidate.size() >= best_count) return {};
            if (!validate_partition(inst, candidate, edge_offset)) return {};
            return candidate;
        };

        bool improved = false;
        for (int s = 0; s < (int)seeds.size() &&
                elapsed(start) < opt.incumbent_repack_seconds &&
                best_count > opt.target_components; s++) {
            attempts_total++;
            double seconds_left = opt.incumbent_repack_seconds - elapsed(start);
            vector<vector<int>> candidate = try_region(
                seeds[s].touched, seeds[s].col, attempts_total, seconds_left);
            if (candidate.empty() && (int)seeds[s].touched.size() <=
                    max(4, opt.incumbent_repack_max_region_components / 2)) {
                vector<int> expanded = expand_region(seeds[s].touched);
                if (expanded.size() != seeds[s].touched.size()) {
                    attempts_total++;
                    candidate = try_region(
                        expanded, seeds[s].col, attempts_total,
                        opt.incumbent_repack_seconds - elapsed(start));
                }
            }
            if (!candidate.empty()) {
                best = std::move(candidate);
                best_count = (int)best.size();
                improved = true;
                cerr << "INCUMBENT_REPACK_BEST components=" << best_count
                     << " round=" << rounds
                     << " attempts=" << attempts_total << "\n";
                break;
            }
        }
        cerr << "INCUMBENT_REPACK_ROUND round=" << rounds
             << " seeds=" << seeds.size()
             << " best_components=" << best_count
             << " improved=" << (improved ? 1 : 0)
             << " seconds=" << elapsed(start) << "\n";
        if (!improved) break;
    }
    cerr << "INCUMBENT_REPACK_DONE rounds=" << rounds
         << " attempts=" << attempts_total
         << " best_components=" << best_count << "\n";
    return best;
}

static string restricted_newick(const Tree &tree, int node, const vector<unsigned char> &active) {
    if (tree.nodes[node].label >= 1) return active[tree.nodes[node].label] ? to_string(tree.nodes[node].label) : "";
    vector<string> kids;
    for (int child : tree.nodes[node].child) {
        string s = restricted_newick(tree, child, active);
        if (!s.empty()) kids.push_back(std::move(s));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return "(" + kids[0] + "," + kids[1] + ")";
}

static string restricted_newick_renumber(const Tree &tree, int node,
        const vector<unsigned char> &active, const vector<int> &old_to_new) {
    int label = tree.nodes[node].label;
    if (label >= 1) {
        if (!active[label]) return "";
        return to_string(old_to_new[label]);
    }
    vector<string> kids;
    for (int child : tree.nodes[node].child) {
        string s = restricted_newick_renumber(tree, child, active, old_to_new);
        if (!s.empty()) kids.push_back(std::move(s));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return "(" + kids[0] + "," + kids[1] + ")";
}

static string top_marker_newick_renumber(const Tree &tree, int node,
        const vector<unsigned char> &in_cluster, int cluster_root,
        const vector<int> &old_to_new, int marker_label) {
    if (node == cluster_root) return to_string(marker_label);
    int label = tree.nodes[node].label;
    if (label >= 1) {
        if (in_cluster[label]) return "";
        return to_string(old_to_new[label]);
    }
    vector<string> kids;
    for (int child : tree.nodes[node].child) {
        string s = top_marker_newick_renumber(
            tree, child, in_cluster, cluster_root, old_to_new, marker_label);
        if (!s.empty()) kids.push_back(std::move(s));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return "(" + kids[0] + "," + kids[1] + ")";
}

struct RecursiveSubinstance {
    Instance inst;
    vector<int> new_to_old;
    int marker = -1;
};

static Instance make_instance_from_newicks(int tree_count, int leaf_count, const vector<string> &newicks) {
    Instance inst;
    inst.tree_count = tree_count;
    inst.leaf_count = leaf_count;
    for (const string &nw : newicks) inst.trees.push_back(parse_tree(nw, leaf_count));
    return inst;
}

static RecursiveSubinstance make_recursive_bottom(const Instance &inst,
        const vector<int> &cluster, bool with_marker) {
    RecursiveSubinstance sub;
    sub.inst.tree_count = inst.tree_count;
    sub.inst.leaf_count = (int)cluster.size() + (with_marker ? 1 : 0);
    sub.new_to_old.assign(sub.inst.leaf_count + 1, -1);
    vector<int> old_to_new(inst.leaf_count + 1, -1);
    vector<unsigned char> active(inst.leaf_count + 1, 0);
    for (int i = 0; i < (int)cluster.size(); i++) {
        old_to_new[cluster[i]] = i + 1;
        sub.new_to_old[i + 1] = cluster[i];
        active[cluster[i]] = 1;
    }
    if (with_marker) {
        sub.marker = sub.inst.leaf_count;
        sub.new_to_old[sub.marker] = -1;
    }
    vector<string> newicks;
    for (const Tree &tree : inst.trees) {
        string inside = restricted_newick_renumber(tree, tree.root, active, old_to_new);
        if (inside.empty()) throw runtime_error("empty recursive bottom restriction");
        if (with_marker) inside = "(" + inside + "," + to_string(sub.marker) + ")";
        newicks.push_back(std::move(inside));
    }
    sub.inst = make_instance_from_newicks(inst.tree_count, sub.inst.leaf_count, newicks);
    return sub;
}

static RecursiveSubinstance make_recursive_top(const Instance &inst,
        const vector<int> &cluster, bool with_marker) {
    vector<unsigned char> in_cluster(inst.leaf_count + 1, 0);
    for (int x : cluster) in_cluster[x] = 1;
    vector<int> outside;
    for (int x = 1; x <= inst.leaf_count; x++) if (!in_cluster[x]) outside.push_back(x);
    RecursiveSubinstance sub;
    sub.inst.tree_count = inst.tree_count;
    sub.inst.leaf_count = (int)outside.size() + (with_marker ? 1 : 0);
    sub.new_to_old.assign(sub.inst.leaf_count + 1, -1);
    vector<int> old_to_new(inst.leaf_count + 1, -1);
    for (int i = 0; i < (int)outside.size(); i++) {
        old_to_new[outside[i]] = i + 1;
        sub.new_to_old[i + 1] = outside[i];
    }
    if (with_marker) {
        sub.marker = sub.inst.leaf_count;
        sub.new_to_old[sub.marker] = -1;
    }
    vector<string> newicks;
    for (const Tree &tree : inst.trees) {
        string nw;
        if (with_marker) {
            int cluster_root = lca_leaf_set(tree, cluster);
            nw = top_marker_newick_renumber(
                tree, tree.root, in_cluster, cluster_root, old_to_new, sub.marker);
        } else {
            vector<unsigned char> active(inst.leaf_count + 1, 0);
            for (int x : outside) active[x] = 1;
            nw = restricted_newick_renumber(tree, tree.root, active, old_to_new);
        }
        if (nw.empty()) throw runtime_error("empty recursive top restriction");
        newicks.push_back(std::move(nw));
    }
    sub.inst = make_instance_from_newicks(inst.tree_count, sub.inst.leaf_count, newicks);
    return sub;
}

static vector<int> choose_recursive_common_cluster(const Instance &inst) {
    vector<vector<int>> first_sets;
    collect_leaf_sets(inst.trees[0], inst.trees[0].root, first_sets);
    vector<unordered_set<string>> by_tree;
    by_tree.reserve(inst.tree_count);
    for (const Tree &tree : inst.trees) {
        vector<vector<int>> sets;
        collect_leaf_sets(tree, tree.root, sets);
        unordered_set<string> keys;
        for (const vector<int> &s : sets) {
            if ((int)s.size() <= 2 || (int)s.size() >= inst.leaf_count) continue;
            keys.insert(key_of(s));
        }
        by_tree.push_back(std::move(keys));
    }
    vector<int> best;
    int best_max_child = INT_MAX;
    int best_balance = -1;
    for (vector<int> cluster : first_sets) {
        int size = (int)cluster.size();
        if (size <= 2 || size >= inst.leaf_count - 1) continue;
        string key = key_of(cluster);
        bool common = true;
        for (int t = 1; t < inst.tree_count; t++) {
            if (!by_tree[t].count(key)) {
                common = false;
                break;
            }
        }
        if (!common) continue;
        int max_child = max(size + 1, inst.leaf_count - size + 1);
        int balance = min(size, inst.leaf_count - size);
        if (max_child < best_max_child ||
                (max_child == best_max_child && balance > best_balance) ||
                (max_child == best_max_child && balance == best_balance &&
                 (best.empty() || cluster < best))) {
            best = std::move(cluster);
            best_max_child = max_child;
            best_balance = balance;
        }
    }
    return best;
}

static vector<int> recursive_edge_offset(const Instance &inst) {
    vector<int> edge_offset(inst.tree_count + 1, 0);
    for (int t = 0; t < inst.tree_count; t++) {
        edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
    }
    return edge_offset;
}

static vector<vector<int>> strict_initial_partition(const Instance &inst, Options opt,
        const vector<int> &edge_offset) {
    Pool pool = build_pool(inst, opt, edge_offset);
    int total_bits = inst.leaf_count + edge_offset.back();
    vector<Column> cols = encode_columns(inst, pool, total_bits);
    vector<int> selected = anchored_pack(cols, inst.leaf_count, opt);
    vector<int> global_selected = global_lns_pack(cols, inst.leaf_count, opt);
    if (!global_selected.empty() &&
            component_count(cols, global_selected, inst.leaf_count) <
            component_count(cols, selected, inst.leaf_count)) {
        selected = std::move(global_selected);
    }
    selected = swap_repair(cols, selected, inst.leaf_count, opt);
    selected = medium_swap_repair(cols, selected, inst.leaf_count, opt);
    vector<vector<int>> partition = selected_partition(cols, selected, inst.leaf_count);
    return singleton_merge_repair(inst, std::move(partition), 1, edge_offset);
}

static ExactDecision exact_branch_decide(
        const Instance &inst, int target_components, double seconds, long long node_limit);
static bool stable_exchange_accept_partition(const Instance &inst, const Options &opt,
        const vector<int> &edge_offset, const vector<vector<int>> &partition,
        const string &context);

static vector<vector<int>> strict_sat_descent(const Instance &inst, const Options &opt,
        vector<vector<int>> partition, const vector<int> &edge_offset,
        int external_lower_floor = -1, const string &floor_context = "") {
    while ((int)partition.size() > 1) {
        if (external_lower_floor > 0 && (int)partition.size() <= external_lower_floor) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << (external_lower_floor - 1)
                 << " method=external_monotone_floor"
                 << " context=" << floor_context << "\n";
            return partition;
        }
        int target = (int)partition.size() - 1;
        bool proved_unsat = false;
        vector<vector<int>> direct = direct_pair_cadical_target(
            inst, target, opt.direct_pair_seconds, &proved_unsat);
        if (proved_unsat) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=direct_pair_sat\n";
            return partition;
        }
        if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                (int)direct.size() < (int)partition.size()) {
            cerr << "STRICT_RECURSIVE_DESCENT from=" << partition.size()
                 << " to=" << direct.size()
                 << " target=" << target << "\n";
            partition = std::move(direct);
            continue;
        }
        SmallSavingsExactResult small_exact = small_savings_exhaustive_pack(
            inst, target, edge_offset, opt.small_savings_exact_seconds,
            opt.small_savings_pair_node_limit,
            opt.small_savings_special_node_limit,
            opt.small_savings_max_target,
            opt.small_savings_component_limit);
        if (!small_exact.partition.empty() &&
                validate_partition(inst, small_exact.partition, edge_offset) &&
                (int)small_exact.partition.size() < (int)partition.size()) {
            cerr << "STRICT_RECURSIVE_SMALL_SAVINGS_DESCENT from="
                 << partition.size()
                 << " to=" << small_exact.partition.size()
                 << " target=" << target << "\n";
            partition = std::move(small_exact.partition);
            continue;
        }
        if (small_exact.proved_infeasible) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=small_savings_exact_lb\n";
            return partition;
        }
        if (opt.reference_cut_seconds > 0.0) {
            bool cut_unsat = false;
            double cut_seconds = min(opt.reference_cut_seconds, 90.0);
            vector<vector<int>> cut_partition = reference_cut_cadical_target(
                inst, target, cut_seconds, &cut_unsat);
            if (cut_unsat) {
                if (!validate_partition(inst, partition, edge_offset)) return {};
                cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                     << " lower_target=" << target
                     << " method=reference_cut_sat\n";
                return partition;
            }
            if (!cut_partition.empty() &&
                    validate_partition(inst, cut_partition, edge_offset) &&
                    (int)cut_partition.size() < (int)partition.size()) {
                cerr << "STRICT_RECURSIVE_REFERENCE_CUT_DESCENT from="
                     << partition.size()
                     << " to=" << cut_partition.size()
                     << " target=" << target << "\n";
                partition = std::move(cut_partition);
                continue;
            }
        }
        if (stable_exchange_accept_partition(
                    inst, opt, edge_offset, partition, "strict_recursive")) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=stable_exchange\n";
            return partition;
        }
        if (opt.exact_branch_seconds > 0.0) {
            ExactDecision branch = exact_branch_decide(
                inst, target, opt.exact_branch_seconds, opt.exact_branch_node_limit);
            if (branch == ExactDecision::No) {
                if (!validate_partition(inst, partition, edge_offset)) return {};
                cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                     << " lower_target=" << target
                     << " method=exact_sibling_branch\n";
                return partition;
            }
        }
        return {};
    }
    return partition;
}

static vector<vector<int>> strict_lazy_sat_descent(const Instance &inst, const Options &opt,
        vector<vector<int>> partition, const vector<int> &edge_offset,
        int external_lower_floor = -1, const string &floor_context = "") {
    while ((int)partition.size() > 1) {
        if (external_lower_floor > 0 && (int)partition.size() <= external_lower_floor) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << (external_lower_floor - 1)
                 << " method=external_monotone_floor"
                 << " context=" << floor_context << "\n";
            return partition;
        }
        int target = (int)partition.size() - 1;
        bool proved_unsat = false;
        vector<vector<int>> direct = direct_pair_cadical_lazy_target(
            inst, target, opt.direct_pair_seconds, &proved_unsat);
        if (proved_unsat) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=direct_pair_lazy_sat\n";
            return partition;
        }
        if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                (int)direct.size() < (int)partition.size()) {
            cerr << "STRICT_RECURSIVE_LAZY_DESCENT from=" << partition.size()
                 << " to=" << direct.size()
                 << " target=" << target << "\n";
            partition = std::move(direct);
            continue;
        }
        SmallSavingsExactResult small_exact = small_savings_exhaustive_pack(
            inst, target, edge_offset, opt.small_savings_exact_seconds,
            opt.small_savings_pair_node_limit,
            opt.small_savings_special_node_limit,
            opt.small_savings_max_target,
            opt.small_savings_component_limit);
        if (!small_exact.partition.empty() &&
                validate_partition(inst, small_exact.partition, edge_offset) &&
                (int)small_exact.partition.size() < (int)partition.size()) {
            cerr << "STRICT_RECURSIVE_SMALL_SAVINGS_DESCENT from="
                 << partition.size()
                 << " to=" << small_exact.partition.size()
                 << " target=" << target << "\n";
            partition = std::move(small_exact.partition);
            continue;
        }
        if (small_exact.proved_infeasible) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=small_savings_exact_lb\n";
            return partition;
        }
        if (opt.reference_cut_seconds > 0.0) {
            bool cut_unsat = false;
            double cut_seconds = min(opt.reference_cut_seconds, 90.0);
            vector<vector<int>> cut_partition = reference_cut_cadical_target(
                inst, target, cut_seconds, &cut_unsat);
            if (cut_unsat) {
                if (!validate_partition(inst, partition, edge_offset)) return {};
                cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                     << " lower_target=" << target
                     << " method=reference_cut_sat\n";
                return partition;
            }
            if (!cut_partition.empty() &&
                    validate_partition(inst, cut_partition, edge_offset) &&
                    (int)cut_partition.size() < (int)partition.size()) {
                cerr << "STRICT_RECURSIVE_REFERENCE_CUT_DESCENT from="
                     << partition.size()
                     << " to=" << cut_partition.size()
                     << " target=" << target << "\n";
                partition = std::move(cut_partition);
                continue;
            }
        }
        if (stable_exchange_accept_partition(
                    inst, opt, edge_offset, partition, "strict_recursive_lazy")) {
            if (!validate_partition(inst, partition, edge_offset)) return {};
            cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                 << " lower_target=" << target
                 << " method=stable_exchange\n";
            return partition;
        }
        if (opt.exact_branch_seconds > 0.0) {
            ExactDecision branch = exact_branch_decide(
                inst, target, opt.exact_branch_seconds, opt.exact_branch_node_limit);
            if (branch == ExactDecision::No) {
                if (!validate_partition(inst, partition, edge_offset)) return {};
                cerr << "STRICT_RECURSIVE_CERT components=" << partition.size()
                     << " lower_target=" << target
                     << " method=exact_sibling_branch\n";
                return partition;
            }
        }
        return {};
    }
    return partition;
}

static vector<vector<int>> map_unmarked_partition(
        const vector<vector<int>> &part, const RecursiveSubinstance &sub) {
    vector<vector<int>> out;
    for (const vector<int> &block : part) {
        vector<int> mapped;
        for (int x : block) {
            int old = sub.new_to_old.at(x);
            if (old < 1) return {};
            mapped.push_back(old);
        }
        sort(mapped.begin(), mapped.end());
        out.push_back(std::move(mapped));
    }
    return out;
}

static vector<int> append_marked_partition_except_marker(
        const vector<vector<int>> &part, const RecursiveSubinstance &sub,
        vector<vector<int>> &out) {
    vector<int> marker_component;
    bool seen_marker = false;
    for (const vector<int> &block : part) {
        vector<int> mapped;
        bool has_marker = false;
        for (int x : block) {
            if (x == sub.marker) {
                has_marker = true;
            } else {
                int old = sub.new_to_old.at(x);
                if (old < 1) return {};
                mapped.push_back(old);
            }
        }
        sort(mapped.begin(), mapped.end());
        if (has_marker) {
            if (seen_marker) return {};
            seen_marker = true;
            marker_component = std::move(mapped);
        } else {
            out.push_back(std::move(mapped));
        }
    }
    if (!seen_marker) return {};
    return marker_component;
}

struct ExactAtomDef {
    int leaf = -1;
    int left = -1;
    int right = -1;
};

struct ExactFNode {
    int label = -1;
    int parent = -1;
    vector<int> child;
    bool alive = true;
};

struct ExactForest {
    vector<ExactFNode> nodes;
    vector<int> roots;

    int add_node(int label) {
        nodes.push_back(ExactFNode{label, -1, {}, true});
        return (int)nodes.size() - 1;
    }

    void cleanup_roots() {
        roots.clear();
        for (int i = 0; i < (int)nodes.size(); i++) {
            if (nodes[i].alive && nodes[i].parent < 0) roots.push_back(i);
        }
        sort(roots.begin(), roots.end());
        roots.erase(unique(roots.begin(), roots.end()), roots.end());
    }

    int order() const {
        int count = 0;
        for (int root : roots) if (root >= 0 && nodes[root].alive && nodes[root].parent < 0) count++;
        return count;
    }

    int find_label(int label) const {
        for (int i = 0; i < (int)nodes.size(); i++) {
            if (nodes[i].alive && nodes[i].label == label) return i;
        }
        return -1;
    }

    int root_of(int node) const {
        while (node >= 0 && nodes[node].parent >= 0) node = nodes[node].parent;
        return node;
    }

    bool is_singleton_label(int label) const {
        int node = find_label(label);
        return node >= 0 && nodes[node].parent < 0 && nodes[node].child.empty();
    }

    bool are_siblings(int a, int b) const {
        int na = find_label(a), nb = find_label(b);
        return na >= 0 && nb >= 0 && nodes[na].parent >= 0 && nodes[na].parent == nodes[nb].parent;
    }

    bool different_components(int a, int b) const {
        int na = find_label(a), nb = find_label(b);
        if (na < 0 || nb < 0) return true;
        return root_of(na) != root_of(nb);
    }

    vector<int> labels() const {
        vector<int> out;
        for (const ExactFNode &node : nodes) {
            if (node.alive && node.label >= 1) out.push_back(node.label);
        }
        sort(out.begin(), out.end());
        out.erase(unique(out.begin(), out.end()), out.end());
        return out;
    }

    void replace_child(int parent, int old_child, int new_child) {
        for (int &child : nodes[parent].child) {
            if (child == old_child) {
                child = new_child;
                return;
            }
        }
    }

    void erase_child(int parent, int child) {
        vector<int> next;
        for (int x : nodes[parent].child) if (x != child) next.push_back(x);
        nodes[parent].child = std::move(next);
    }

    void normalize_up(int node) {
        while (node >= 0 && nodes[node].alive) {
            if (nodes[node].label >= 1) break;
            if (nodes[node].child.size() == 1) {
                int only = nodes[node].child[0];
                int parent = nodes[node].parent;
                nodes[node].alive = false;
                nodes[node].child.clear();
                if (parent < 0) {
                    nodes[only].parent = -1;
                    node = only;
                } else {
                    replace_child(parent, node, only);
                    nodes[only].parent = parent;
                    node = parent;
                }
                continue;
            }
            if (nodes[node].child.empty()) {
                int parent = nodes[node].parent;
                nodes[node].alive = false;
                if (parent >= 0) {
                    erase_child(parent, node);
                    node = parent;
                    continue;
                }
            }
            break;
        }
        cleanup_roots();
    }

    void cut_edge_to_node(int child) {
        if (child < 0 || !nodes[child].alive) return;
        int parent = nodes[child].parent;
        if (parent < 0) return;
        erase_child(parent, child);
        nodes[child].parent = -1;
        normalize_up(parent);
        cleanup_roots();
    }

    void make_singleton(int label) {
        int node = find_label(label);
        if (node < 0) return;
        cut_edge_to_node(node);
    }

    bool shrink_sibling_pair(int a, int b, int new_label) {
        int na = find_label(a), nb = find_label(b);
        if (na < 0 || nb < 0) return false;
        int parent = nodes[na].parent;
        if (parent < 0 || parent != nodes[nb].parent) return false;
        nodes[na].alive = false;
        nodes[nb].alive = false;
        nodes[parent].child.clear();
        nodes[parent].label = new_label;
        cleanup_roots();
        return true;
    }

    vector<int> path_side_cut_children(int a, int b) const {
        int na = find_label(a), nb = find_label(b);
        if (na < 0 || nb < 0 || root_of(na) != root_of(nb)) return {};
        unordered_map<int, int> pos;
        vector<int> aa;
        for (int x = na; x >= 0; x = nodes[x].parent) {
            pos[x] = (int)aa.size();
            aa.push_back(x);
        }
        vector<int> bb;
        int lca = -1;
        for (int x = nb; x >= 0; x = nodes[x].parent) {
            bb.push_back(x);
            if (pos.count(x)) {
                lca = x;
                break;
            }
        }
        if (lca < 0) return {};
        unordered_set<int> on_path;
        for (int x = na; x != lca; x = nodes[x].parent) on_path.insert(x);
        on_path.insert(lca);
        for (int x = nb; x != lca; x = nodes[x].parent) on_path.insert(x);
        vector<int> cuts;
        for (int v : on_path) {
            if (v == lca) continue;
            for (int child : nodes[v].child) {
                if (!on_path.count(child)) cuts.push_back(child);
            }
        }
        sort(cuts.begin(), cuts.end());
        cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());
        return cuts;
    }

    bool cut_path_side_edges(int a, int b) {
        vector<int> cuts = path_side_cut_children(a, b);
        if (cuts.empty()) return false;
        for (int child : cuts) cut_edge_to_node(child);
        cleanup_roots();
        return true;
    }

    pair<int, int> first_sibling_pair() const {
        for (const ExactFNode &node : nodes) {
            if (!node.alive || node.label >= 1 || node.child.size() != 2) continue;
            int a = node.child[0], b = node.child[1];
            if (!nodes[a].alive || !nodes[b].alive) continue;
            if (nodes[a].label >= 1 && nodes[b].label >= 1) {
                return {nodes[a].label, nodes[b].label};
            }
        }
        return {-1, -1};
    }

    vector<pair<int, int>> sibling_pairs() const {
        vector<pair<int, int>> out;
        for (const ExactFNode &node : nodes) {
            if (!node.alive || node.label >= 1 || node.child.size() != 2) continue;
            int a = node.child[0], b = node.child[1];
            if (!nodes[a].alive || !nodes[b].alive) continue;
            if (nodes[a].label >= 1 && nodes[b].label >= 1) {
                int la = nodes[a].label, lb = nodes[b].label;
                if (lb < la) swap(la, lb);
                out.push_back({la, lb});
            }
        }
        sort(out.begin(), out.end());
        out.erase(unique(out.begin(), out.end()), out.end());
        return out;
    }
};

static int clone_exact_tree_node(const Tree &tree, int node, ExactForest &forest) {
    int out = forest.add_node(tree.nodes[node].label);
    for (int child : tree.nodes[node].child) {
        int c = clone_exact_tree_node(tree, child, forest);
        forest.nodes[c].parent = out;
        forest.nodes[out].child.push_back(c);
    }
    return out;
}

static ExactForest exact_forest_from_tree(const Tree &tree) {
    ExactForest forest;
    int root = clone_exact_tree_node(tree, tree.root, forest);
    forest.roots.push_back(root);
    forest.cleanup_roots();
    return forest;
}

static int clone_atom_into_forest(int atom, const vector<ExactAtomDef> &atoms,
        ExactForest &forest, int expand_from) {
    if (atom < expand_from) return forest.add_node(atom);
    if (atom >= (int)atoms.size()) return -1;
    const ExactAtomDef &def = atoms[atom];
    if (def.leaf >= 1) return forest.add_node(def.leaf);
    int root = forest.add_node(-1);
    int left = clone_atom_into_forest(def.left, atoms, forest, expand_from);
    int right = clone_atom_into_forest(def.right, atoms, forest, expand_from);
    forest.nodes[left].parent = root;
    forest.nodes[right].parent = root;
    forest.nodes[root].child = {left, right};
    return root;
}

static ExactForest exact_forest_expand_singleton_atoms(
        const ExactForest &forest, const vector<ExactAtomDef> &atoms, int expand_from) {
    ExactForest out;
    for (int root : forest.roots) {
        if (!forest.nodes[root].alive || forest.nodes[root].label < 1) continue;
        int cloned = clone_atom_into_forest(forest.nodes[root].label, atoms, out, expand_from);
        if (cloned >= 0) out.roots.push_back(cloned);
    }
    out.cleanup_roots();
    return out;
}

static string exact_atom_signature(int atom, const vector<ExactAtomDef> &atoms) {
    const ExactAtomDef &def = atoms[atom];
    if (def.leaf >= 1) return to_string(def.leaf);
    return combine_sig(exact_atom_signature(def.left, atoms), exact_atom_signature(def.right, atoms));
}

static string exact_node_signature(
        const ExactForest &forest, int node, const vector<ExactAtomDef> &atoms) {
    const ExactFNode &fn = forest.nodes[node];
    if (fn.label >= 1) return exact_atom_signature(fn.label, atoms);
    vector<string> kids;
    for (int child : fn.child) {
        if (forest.nodes[child].alive) kids.push_back(exact_node_signature(forest, child, atoms));
    }
    if (kids.empty()) return "";
    if (kids.size() == 1) return kids[0];
    return combine_sig(kids[0], kids[1]);
}

static string exact_forest_signature(
        const ExactForest &forest, const vector<ExactAtomDef> &atoms) {
    vector<string> comps;
    for (int root : forest.roots) {
        if (root >= 0 && forest.nodes[root].alive && forest.nodes[root].parent < 0) {
            comps.push_back(exact_node_signature(forest, root, atoms));
        }
    }
    sort(comps.begin(), comps.end());
    string out;
    for (const string &comp : comps) {
        out += comp;
        out += "|";
    }
    return out;
}

struct ExactBranchSearch {
    int original_n = 0;
    int target = 0;
    Clock::time_point deadline;
    long long node_limit = 0;
    long long nodes = 0;
    bool dynamic_pair_selection = false;
    bool global_conflict_branching = false;
    unordered_set<string> failed;

    bool expired() {
        return Clock::now() >= deadline || (node_limit > 0 && nodes >= node_limit);
    }

    bool any_forest_exceeds_target(const vector<ExactForest> &forests) const {
        for (const ExactForest &forest : forests) {
            if (forest.order() > target) return true;
        }
        return false;
    }

    string state_key(const vector<ExactForest> &forests,
            const vector<ExactAtomDef> &atoms, int pair_atom_base) {
        string key = to_string(pair_atom_base) + "$";
        for (const ExactForest &forest : forests) {
            key += exact_forest_signature(forest, atoms);
            key += "#";
        }
        return key;
    }

    bool atoms_are_initial(const vector<ExactAtomDef> &atoms) const {
        if ((int)atoms.size() != original_n + 1) return false;
        for (int i = 1; i <= original_n; i++) {
            if (atoms[i].leaf != i || atoms[i].left != -1 || atoms[i].right != -1) return false;
        }
        return true;
    }

    int ordered_pair_score(const ExactForest &f1, const ExactForest &f2) const {
        int total = 0;
        int useful = 0;
        for (auto cand : f2.sibling_pairs()) {
            int a = cand.first, b = cand.second;
            useful++;
            if (f1.are_siblings(a, b)) {
                total += 20;
            } else if (f1.different_components(a, b)) {
                total += 80;
            } else {
                total += 140 + (int)f1.path_side_cut_children(a, b).size() * 15;
            }
        }
        if (useful == 0) return INT_MIN;
        return total * 1000 + useful;
    }

    int sibling_pair_score(const vector<ExactForest> &forests, int a, int b) const {
        const ExactForest &f1 = forests[0];
        int score = 0;
        if (f1.are_siblings(a, b)) score = 1000000;
        else if (f1.different_components(a, b)) score = 1000;
        else score = 2000 + (int)f1.path_side_cut_children(a, b).size() * 20;
        if (dynamic_pair_selection) {
            for (int i = 2; i < (int)forests.size(); i++) {
                const ExactForest &other = forests[i];
                if (other.are_siblings(a, b)) {
                    score += 10;
                } else if (other.different_components(a, b)) {
                    score += 90;
                } else {
                    score += 120 + (int)other.path_side_cut_children(a, b).size() * 10;
                }
            }
        }
        score -= min(a, b);
        return score;
    }

    void choose_active_pair(vector<ExactForest> &forests,
            const vector<ExactAtomDef> &atoms, int pair_atom_base) {
        if (!dynamic_pair_selection || forests.size() < 3 ||
                pair_atom_base < 0 || (int)atoms.size() != pair_atom_base) {
            return;
        }
        int best_i = 0, best_j = 1;
        int best_score = ordered_pair_score(forests[0], forests[1]);
        for (int i = 0; i < (int)forests.size(); i++) {
            for (int j = 0; j < (int)forests.size(); j++) {
                if (i == j) continue;
                int score = ordered_pair_score(forests[i], forests[j]);
                if (score > best_score) {
                    best_score = score;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        if (best_i == 0 && best_j == 1) return;
        vector<ExactForest> reordered;
        reordered.reserve(forests.size());
        reordered.push_back(std::move(forests[best_i]));
        reordered.push_back(std::move(forests[best_j]));
        for (int i = 0; i < (int)forests.size(); i++) {
            if (i != best_i && i != best_j) reordered.push_back(std::move(forests[i]));
        }
        forests = std::move(reordered);
    }

    bool reduce_common_cherries(vector<ExactForest> &forests, vector<ExactAtomDef> &atoms) {
        bool any = false;
        bool changed = true;
        while (changed) {
            changed = false;
            if (forests.empty()) break;
            for (auto cand : forests[0].sibling_pairs()) {
                int a = cand.first, b = cand.second;
                bool common = true;
                for (int i = 1; i < (int)forests.size(); i++) {
                    if (!forests[i].are_siblings(a, b)) {
                        common = false;
                        break;
                    }
                }
                if (!common) continue;
                int new_label = (int)atoms.size();
                ExactAtomDef def;
                def.left = a;
                def.right = b;
                atoms.push_back(def);
                for (ExactForest &forest : forests) {
                    if (!forest.shrink_sibling_pair(a, b, new_label)) return any;
                }
                any = true;
                changed = true;
                break;
            }
        }
        return any;
    }

    bool propagate_all_singletons(vector<ExactForest> &forests, bool &changed) {
        changed = false;
        if (forests.empty()) return true;
        vector<int> labels;
        for (const ExactForest &forest : forests) {
            vector<int> here = forest.labels();
            labels.insert(labels.end(), here.begin(), here.end());
        }
        sort(labels.begin(), labels.end());
        labels.erase(unique(labels.begin(), labels.end()), labels.end());
        for (int label : labels) {
            bool forced = false;
            for (const ExactForest &forest : forests) {
                if (forest.is_singleton_label(label)) {
                    forced = true;
                    break;
                }
            }
            if (!forced) continue;
            for (ExactForest &forest : forests) {
                if (!forest.is_singleton_label(label)) {
                    forest.make_singleton(label);
                    changed = true;
                    if (forest.order() > target) return false;
                }
            }
        }
        return true;
    }

    bool reduce_phase_boundary(vector<ExactForest> &forests, vector<ExactAtomDef> &atoms) {
        bool changed = true;
        while (changed) {
            changed = false;
            bool singleton_changed = false;
            if (!propagate_all_singletons(forests, singleton_changed)) return false;
            changed = changed || singleton_changed;
            bool cherry_changed = reduce_common_cherries(forests, atoms);
            changed = changed || cherry_changed;
            if (any_forest_exceeds_target(forests)) return false;
        }
        return true;
    }

    int sibling_conflict_extra_for_base(
            const vector<ExactForest> &forests, int base_index) const {
        if (base_index < 0 || base_index >= (int)forests.size()) return 0;
        const ExactForest &base = forests[base_index];
        int extra = 0;
        for (auto cand : base.sibling_pairs()) {
            int a = cand.first, b = cand.second;
            bool conflict = false;
            for (int i = 0; i < (int)forests.size(); i++) {
                if (i == base_index) continue;
                if (forests[i].find_label(a) >= 0 && forests[i].find_label(b) >= 0 &&
                        forests[i].different_components(a, b)) {
                    conflict = true;
                    break;
                }
            }
            if (conflict) extra++;
        }
        return extra;
    }

    bool sibling_conflict_bound_exceeds_target(const vector<ExactForest> &forests) const {
        for (int i = 0; i < (int)forests.size(); i++) {
            int extra = sibling_conflict_extra_for_base(forests, i);
            if (forests[i].order() + extra > target) return true;
        }
        return false;
    }

    vector<pair<int, int>> global_conflict_edges(const vector<ExactForest> &forests) const {
        vector<pair<int, int>> edges;
        for (int base_index = 0; base_index < (int)forests.size(); base_index++) {
            for (auto cand : forests[base_index].sibling_pairs()) {
                int a = cand.first, b = cand.second;
                bool conflict = false;
                for (int i = 0; i < (int)forests.size(); i++) {
                    if (i == base_index) continue;
                    if (forests[i].find_label(a) >= 0 && forests[i].find_label(b) >= 0 &&
                            forests[i].different_components(a, b)) {
                        conflict = true;
                        break;
                    }
                }
                if (!conflict) continue;
                if (b < a) swap(a, b);
                edges.push_back({a, b});
            }
        }
        sort(edges.begin(), edges.end());
        edges.erase(unique(edges.begin(), edges.end()), edges.end());
        return edges;
    }

    pair<int, int> choose_global_conflict_edge(
            const vector<ExactForest> &forests) const {
        vector<pair<int, int>> edges = global_conflict_edges(forests);
        if (edges.empty()) return {-1, -1};
        unordered_map<int, int> degree;
        for (auto edge : edges) {
            degree[edge.first]++;
            degree[edge.second]++;
        }
        pair<int, int> best{-1, -1};
        int best_score = INT_MIN;
        for (auto edge : edges) {
            int score = degree[edge.first] + degree[edge.second];
            for (const ExactForest &forest : forests) {
                if (forest.are_siblings(edge.first, edge.second)) score += 5;
                else if (forest.different_components(edge.first, edge.second)) score += 3;
            }
            score -= min(edge.first, edge.second);
            if (score > best_score) {
                best_score = score;
                best = edge;
            }
        }
        return best;
    }

    void make_singleton_all(vector<ExactForest> &forests, int label) const {
        for (ExactForest &forest : forests) {
            if (forest.find_label(label) >= 0) forest.make_singleton(label);
        }
    }

    bool global_conflict_matching_bound_exceeds_target(
            const vector<ExactForest> &forests) const {
        vector<pair<int, int>> edges = global_conflict_edges(forests);
        if (edges.empty()) return false;
        auto greedy_matching_size = [](const vector<pair<int, int>> &ordered_edges) {
            unordered_set<int> used;
            int matching = 0;
            for (auto edge : ordered_edges) {
                if (used.count(edge.first) || used.count(edge.second)) continue;
                used.insert(edge.first);
                used.insert(edge.second);
                matching++;
            }
            return matching;
        };
        for (int base_index = 0; base_index < (int)forests.size(); base_index++) {
            const ExactForest &base = forests[base_index];
            vector<pair<int, int>> feasible;
            unordered_map<int, int> degree;
            for (auto edge : edges) {
                int a = edge.first, b = edge.second;
                if (base.find_label(a) < 0 || base.find_label(b) < 0) continue;
                if (base.is_singleton_label(a) || base.is_singleton_label(b)) continue;
                feasible.push_back(edge);
                degree[a]++;
                degree[b]++;
            }
            int matching = greedy_matching_size(feasible);
            vector<pair<int, int>> by_reverse = feasible;
            reverse(by_reverse.begin(), by_reverse.end());
            matching = max(matching, greedy_matching_size(by_reverse));
            vector<pair<int, int>> by_degree = feasible;
            sort(by_degree.begin(), by_degree.end(), [&](const auto &x, const auto &y) {
                int dx = degree[x.first] + degree[x.second];
                int dy = degree[y.first] + degree[y.second];
                if (dx != dy) return dx > dy;
                return x < y;
            });
            matching = max(matching, greedy_matching_size(by_degree));
            vector<pair<int, int>> by_low_degree = feasible;
            sort(by_low_degree.begin(), by_low_degree.end(), [&](const auto &x, const auto &y) {
                int dx = degree[x.first] + degree[x.second];
                int dy = degree[y.first] + degree[y.second];
                if (dx != dy) return dx < dy;
                return x < y;
            });
            matching = max(matching, greedy_matching_size(by_low_degree));
            if (base.order() + matching > target) return true;
        }
        return false;
    }

    bool bounded_vertex_cover_can_cover(vector<pair<int, int>> edges,
            int budget, int &nodes, int node_limit) const {
        if (edges.empty()) return true;
        if (budget < 0) return false;
        if (budget == 0) return false;
        if (++nodes > node_limit) return true;

        unordered_map<int, int> degree;
        for (auto edge : edges) {
            degree[edge.first]++;
            degree[edge.second]++;
        }
        pair<int, int> pivot = edges[0];
        int pivot_score = -1;
        for (auto edge : edges) {
            int score = degree[edge.first] + degree[edge.second];
            if (score > pivot_score) {
                pivot_score = score;
                pivot = edge;
            }
        }
        auto remove_incident = [&](int chosen) {
            vector<pair<int, int>> next;
            next.reserve(edges.size());
            for (auto edge : edges) {
                if (edge.first != chosen && edge.second != chosen) next.push_back(edge);
            }
            return next;
        };
        if (bounded_vertex_cover_can_cover(
                remove_incident(pivot.first), budget - 1, nodes, node_limit)) {
            return true;
        }
        return bounded_vertex_cover_can_cover(
            remove_incident(pivot.second), budget - 1, nodes, node_limit);
    }

    bool global_conflict_vertex_cover_exceeds_target(
            const vector<ExactForest> &forests) const {
        vector<pair<int, int>> edges = global_conflict_edges(forests);
        if (edges.empty()) return false;
        for (int base_index = 0; base_index < (int)forests.size(); base_index++) {
            const ExactForest &base = forests[base_index];
            int budget = target - base.order();
            if (budget < 0) return true;
            if (budget > 16) continue;
            vector<pair<int, int>> positive_edges;
            positive_edges.reserve(edges.size());
            for (auto edge : edges) {
                int a = edge.first, b = edge.second;
                if (base.find_label(a) < 0 || base.find_label(b) < 0) continue;
                if (base.is_singleton_label(a) || base.is_singleton_label(b)) continue;
                positive_edges.push_back(edge);
            }
            if (positive_edges.empty()) continue;
            int nodes = 0;
            bool can_cover = bounded_vertex_cover_can_cover(
                std::move(positive_edges), budget, nodes, 2048);
            if (!can_cover) return true;
        }
        return false;
    }

    bool synchronize_singletons(ExactForest &f1, ExactForest &f2) {
        bool changed = true;
        while (changed) {
            changed = false;
            vector<int> labels = f1.labels();
            vector<int> labels2 = f2.labels();
            labels.insert(labels.end(), labels2.begin(), labels2.end());
            sort(labels.begin(), labels.end());
            labels.erase(unique(labels.begin(), labels.end()), labels.end());
            for (int label : labels) {
                bool s1 = f1.is_singleton_label(label);
                bool s2 = f2.is_singleton_label(label);
                if (s1 != s2) {
                    f1.make_singleton(label);
                    f2.make_singleton(label);
                    changed = true;
                    if (f1.order() > target || f2.order() > target) return false;
                    break;
                }
            }
        }
        return true;
    }

    ExactDecision solve(vector<ExactForest> forests, vector<ExactAtomDef> atoms,
            int pair_atom_base = -1) {
        if (expired()) return ExactDecision::Unknown;
        nodes++;
        if (forests.empty()) return ExactDecision::No;
        if (forests.size() == 1) {
            return forests[0].order() <= target ? ExactDecision::Yes : ExactDecision::No;
        }
        if (any_forest_exceeds_target(forests)) return ExactDecision::No;
        if (pair_atom_base < 0) {
            if (!reduce_phase_boundary(forests, atoms)) return ExactDecision::No;
            pair_atom_base = (int)atoms.size();
        }
        if (global_conflict_branching && pair_atom_base == (int)atoms.size()) {
            pair<int, int> conflict = choose_global_conflict_edge(forests);
            if (conflict.first >= 0) {
                bool saw_unknown = false;
                for (int label : {conflict.first, conflict.second}) {
                    vector<ExactForest> next = forests;
                    make_singleton_all(next, label);
                    if (any_forest_exceeds_target(next) ||
                            sibling_conflict_bound_exceeds_target(next) ||
                            global_conflict_matching_bound_exceeds_target(next) ||
                            global_conflict_vertex_cover_exceeds_target(next)) {
                        continue;
                    }
                    ExactDecision ans = solve(std::move(next), atoms);
                    if (ans == ExactDecision::Yes) return ExactDecision::Yes;
                    if (ans == ExactDecision::Unknown) saw_unknown = true;
                }
                if (saw_unknown) return ExactDecision::Unknown;
                return ExactDecision::No;
            }
        }
        choose_active_pair(forests, atoms, pair_atom_base);
        string key = state_key(forests, atoms, pair_atom_base);
        if (failed.count(key)) return ExactDecision::No;

        ExactForest &f1 = forests[0];
        ExactForest &f2 = forests[1];
        if (!synchronize_singletons(f1, f2)) {
            failed.insert(key);
            return ExactDecision::No;
        }
        if (any_forest_exceeds_target(forests)) {
            failed.insert(key);
            return ExactDecision::No;
        }
        if (sibling_conflict_bound_exceeds_target(forests) ||
                global_conflict_matching_bound_exceeds_target(forests) ||
                global_conflict_vertex_cover_exceeds_target(forests)) {
            failed.insert(key);
            return ExactDecision::No;
        }
        key = state_key(forests, atoms, pair_atom_base);
        if (failed.count(key)) return ExactDecision::No;

        pair<int, int> sib{-1, -1};
        int best_score = INT_MIN;
        for (auto cand : f2.sibling_pairs()) {
            int a = cand.first, b = cand.second;
            int score = sibling_pair_score(forests, a, b);
            if (score > best_score) {
                best_score = score;
                sib = cand;
            }
        }
        if (sib.first < 0) {
            ExactForest fstar = exact_forest_expand_singleton_atoms(f2, atoms, pair_atom_base);
            vector<ExactAtomDef> next_atoms(atoms.begin(), atoms.begin() + pair_atom_base);
            vector<ExactForest> next;
            next.push_back(std::move(fstar));
            for (int i = 2; i < (int)forests.size(); i++) next.push_back(std::move(forests[i]));
            ExactDecision ans = solve(std::move(next), std::move(next_atoms));
            if (ans == ExactDecision::No) failed.insert(key);
            return ans;
        }

        int a = sib.first, b = sib.second;
        vector<pair<ExactForest, ExactForest>> branches;
        auto add_singleton_branch = [&](int label) {
            ExactForest nf1 = f1, nf2 = f2;
            nf1.make_singleton(label);
            nf2.make_singleton(label);
            if (nf1.order() > target || nf2.order() > target) return;
            branches.push_back({std::move(nf1), std::move(nf2)});
        };

        if (f1.different_components(a, b)) {
            add_singleton_branch(a);
            add_singleton_branch(b);
        } else if (f1.are_siblings(a, b)) {
            int new_label = (int)atoms.size();
            ExactAtomDef def;
            def.left = a;
            def.right = b;
            atoms.push_back(def);
            ExactForest nf1 = f1, nf2 = f2;
            if (!nf1.shrink_sibling_pair(a, b, new_label) ||
                    !nf2.shrink_sibling_pair(a, b, new_label)) {
                failed.insert(key);
                return ExactDecision::No;
            }
            forests[0] = std::move(nf1);
            forests[1] = std::move(nf2);
            ExactDecision ans = solve(std::move(forests), std::move(atoms), pair_atom_base);
            if (ans == ExactDecision::No) failed.insert(key);
            return ans;
        } else {
            add_singleton_branch(a);
            add_singleton_branch(b);
            ExactForest nf1 = f1, nf2 = f2;
            if (nf1.cut_path_side_edges(a, b)) {
                if (nf1.order() > target || nf2.order() > target) {
                    // infeasible branch
                } else {
                branches.push_back({std::move(nf1), std::move(nf2)});
                }
            }
        }

        bool saw_unknown = false;
        for (auto &branch : branches) {
            if (expired()) return ExactDecision::Unknown;
            vector<ExactForest> next = forests;
            next[0] = std::move(branch.first);
            next[1] = std::move(branch.second);
            if (any_forest_exceeds_target(next) ||
                    sibling_conflict_bound_exceeds_target(next) ||
                    global_conflict_matching_bound_exceeds_target(next) ||
                    global_conflict_vertex_cover_exceeds_target(next)) {
                continue;
            }
            ExactDecision ans = solve(std::move(next), atoms, pair_atom_base);
            if (ans == ExactDecision::Yes) return ExactDecision::Yes;
            if (ans == ExactDecision::Unknown) saw_unknown = true;
        }
        if (saw_unknown) return ExactDecision::Unknown;
        failed.insert(key);
        return ExactDecision::No;
    }

    vector<ExactAtomDef> initial_atoms() const {
        vector<ExactAtomDef> atoms(original_n + 1);
        for (int i = 1; i <= original_n; i++) atoms[i].leaf = i;
        return atoms;
    }
};

static ExactDecision exact_branch_decide_once(
        const Instance &inst, int target_components, double seconds, long long node_limit,
        bool dynamic_pair_selection, bool global_conflict_branching) {
    if (seconds <= 0.0 || target_components < 1) return ExactDecision::Unknown;
    ExactBranchSearch search;
    search.original_n = inst.leaf_count;
    search.target = target_components;
    search.deadline = Clock::now() + chrono::milliseconds((int)(seconds * 1000));
    search.node_limit = node_limit;
    search.dynamic_pair_selection = dynamic_pair_selection;
    search.global_conflict_branching = global_conflict_branching;
    vector<ExactForest> forests;
    forests.reserve(inst.tree_count);
    for (const Tree &tree : inst.trees) forests.push_back(exact_forest_from_tree(tree));
    vector<ExactAtomDef> atoms = search.initial_atoms();
    ExactDecision ans = search.solve(std::move(forests), std::move(atoms));
    cerr << "EXACT_BRANCH_RESULT target=" << target_components
         << " mode=" << (global_conflict_branching ? "global_conflict" :
                 dynamic_pair_selection ? "dynamic_pair" : "fixed_pair")
         << " result=" << (ans == ExactDecision::No ? "NO" : ans == ExactDecision::Yes ? "YES" : "UNKNOWN")
         << " nodes=" << search.nodes
         << " memo_no=" << search.failed.size() << "\n";
    return ans;
}

static ExactDecision exact_branch_decide(
        const Instance &inst, int target_components, double seconds, long long node_limit) {
    if (seconds <= 0.0 || target_components < 1) return ExactDecision::Unknown;
    Clock::time_point start = Clock::now();
    double first_seconds = seconds;
    if (inst.tree_count >= 3 && seconds > 2.0) {
        first_seconds = min(seconds, max(2.0, seconds * 0.6));
    }
    ExactDecision fixed = exact_branch_decide_once(
        inst, target_components, first_seconds, node_limit, false, false);
    if (fixed != ExactDecision::Unknown) return fixed;
    double used = chrono::duration<double>(Clock::now() - start).count();
    double remaining = seconds - used;
    if (inst.tree_count < 3 || remaining <= 0.05) return ExactDecision::Unknown;
    return exact_branch_decide_once(inst, target_components, remaining, node_limit, true, true);
}

static Instance pair_subinstance(const Instance &inst, int first, int second) {
    Instance out;
    out.tree_count = 2;
    out.leaf_count = inst.leaf_count;
    out.trees.push_back(inst.trees[first]);
    out.trees.push_back(inst.trees[second]);
    return out;
}

static bool pairwise_branch_rule_out(
        const Instance &inst, int target_components, double seconds,
        long long node_limit, int max_pairs) {
    if (seconds <= 0.0 || target_components < 1 || inst.tree_count < 2) return false;
    vector<pair<int, int>> pairs;
    for (int i = 0; i < inst.tree_count; i++) {
        for (int j = i + 1; j < inst.tree_count; j++) pairs.push_back({i, j});
    }
    if (max_pairs > 0 && (int)pairs.size() > max_pairs) pairs.resize(max_pairs);
    Clock::time_point start = Clock::now();
    for (int idx = 0; idx < (int)pairs.size(); idx++) {
        double used = elapsed(start);
        double remaining = seconds - used;
        if (remaining <= 0.05) break;
        int pairs_left = (int)pairs.size() - idx;
        double pair_seconds = max(0.05, remaining / pairs_left);
        Instance pair_inst = pair_subinstance(inst, pairs[idx].first, pairs[idx].second);
        ExactDecision decision = exact_branch_decide_once(
            pair_inst, target_components, pair_seconds, node_limit, false, false);
        cerr << "PAIRWISE_BRANCH_LB pair=" << (pairs[idx].first + 1)
             << "," << (pairs[idx].second + 1)
             << " target=" << target_components
             << " result=" << (decision == ExactDecision::No ? "NO" :
                     decision == ExactDecision::Yes ? "YES" : "UNKNOWN")
             << "\n";
        if (decision == ExactDecision::No) {
            cerr << "PAIRWISE_BRANCH_RULE_OUT target=" << target_components
                 << " pair=" << (pairs[idx].first + 1) << ","
                 << (pairs[idx].second + 1) << "\n";
            return true;
        }
    }
    return false;
}

static string recursive_tree_key(const Tree &tree, int node) {
    const TreeNode &tn = tree.nodes[node];
    if (tn.child.empty()) return to_string(tn.label);
    vector<string> child_keys;
    child_keys.reserve(tn.child.size());
    for (int child : tn.child) child_keys.push_back(recursive_tree_key(tree, child));
    sort(child_keys.begin(), child_keys.end());
    string out = "(";
    for (int i = 0; i < (int)child_keys.size(); i++) {
        if (i) out += ",";
        out += child_keys[i];
    }
    out += ")";
    return out;
}

static string recursive_instance_key(const Instance &inst) {
    string key = to_string(inst.tree_count) + "|" + to_string(inst.leaf_count);
    for (const Tree &tree : inst.trees) {
        key += "|";
        key += recursive_tree_key(tree, tree.root);
    }
    return key;
}

static vector<vector<int>> strict_recursive_exact(const Instance &inst,
        const Options &opt, int depth, int external_lower_floor = -1,
        const string &floor_context = "") {
    if (depth > 80) return {};
    static unordered_map<string, vector<vector<int>>> memo;
    string memo_key = recursive_instance_key(inst);
    auto memo_it = memo.find(memo_key);
    if (memo_it != memo.end()) {
        cerr << "STRICT_RECURSIVE_CACHE_HIT n=" << inst.leaf_count
             << " components=" << memo_it->second.size()
             << " depth=" << depth << "\n";
        return memo_it->second;
    }
    auto remember = [&](vector<vector<int>> result) {
        if (!result.empty()) memo.emplace(memo_key, result);
        return result;
    };
    Options local_opt = opt;
    if (depth > 0) {
        local_opt.max_pool = min(local_opt.max_pool, 200000);
        local_opt.pack_pool_limit = min(local_opt.pack_pool_limit, 200000);
        local_opt.pool_seconds = min(local_opt.pool_seconds, 60.0);
    }
    vector<int> edge_offset = recursive_edge_offset(inst);
    vector<int> all_leaves(inst.leaf_count);
    for (int i = 0; i < inst.leaf_count; i++) all_leaves[i] = i + 1;
    if (valid_topology(inst, all_leaves)) {
        vector<vector<int>> one{all_leaves};
        if (validate_partition(inst, one, edge_offset)) {
            cerr << "STRICT_RECURSIVE_CERT components=1 method=all_tree_topology n="
                 << inst.leaf_count << " depth=" << depth << "\n";
            return remember(std::move(one));
        }
    }
    if (inst.leaf_count <= local_opt.direct_pair_leaf_limit) {
        vector<vector<int>> incumbent = strict_initial_partition(inst, local_opt, edge_offset);
        return remember(strict_sat_descent(
            inst, local_opt, std::move(incumbent), edge_offset,
            external_lower_floor, floor_context));
    }

    vector<int> cluster = choose_recursive_common_cluster(inst);
    if (cluster.empty()) {
        if (inst.leaf_count <= 128) {
            vector<vector<int>> incumbent = strict_initial_partition(inst, local_opt, edge_offset);
            if (local_opt.full_sat_before_lazy) {
                vector<vector<int>> full = strict_sat_descent(
                    inst, local_opt, incumbent, edge_offset,
                    external_lower_floor, floor_context);
                if (!full.empty()) return remember(std::move(full));
            }
            return remember(strict_lazy_sat_descent(
                inst, local_opt, std::move(incumbent), edge_offset,
                external_lower_floor, floor_context));
        }
        return {};
    }
    cerr << "STRICT_CLUSTER_SPLIT depth=" << depth
         << " n=" << inst.leaf_count
         << " cluster_size=" << cluster.size()
         << " top_size=" << (inst.leaf_count - (int)cluster.size()) << "\n";

    RecursiveSubinstance bottom = make_recursive_bottom(inst, cluster, false);
    RecursiveSubinstance bottom_marker = make_recursive_bottom(inst, cluster, true);
    RecursiveSubinstance top = make_recursive_top(inst, cluster, false);
    RecursiveSubinstance top_marker = make_recursive_top(inst, cluster, true);

    vector<vector<int>> b = strict_recursive_exact(bottom.inst, local_opt, depth + 1);
    if (b.empty()) return {};
    vector<vector<int>> bm = strict_recursive_exact(
        bottom_marker.inst, local_opt, depth + 1, (int)b.size(),
        "bottom_marker_from_unmarked");
    if (bm.empty()) return {};
    vector<vector<int>> t = strict_recursive_exact(top.inst, local_opt, depth + 1);
    if (t.empty()) return {};
    vector<vector<int>> tm = strict_recursive_exact(
        top_marker.inst, local_opt, depth + 1, (int)t.size(),
        "top_marker_from_unmarked");
    if (tm.empty()) return {};

    int B = (int)b.size(), BM = (int)bm.size(), T = (int)t.size(), TM = (int)tm.size();
    if (BM < B || BM > B + 1 || TM < T || TM > T + 1) return {};
    bool bridge = (BM == B && TM == T);
    int expected = B + T - (bridge ? 1 : 0);
    vector<vector<int>> combined;
    if (bridge) {
        vector<int> inside_bridge = append_marked_partition_except_marker(bm, bottom_marker, combined);
        vector<int> outside_bridge = append_marked_partition_except_marker(tm, top_marker, combined);
        if (inside_bridge.empty() && outside_bridge.empty()) return {};
        vector<int> joined = inside_bridge;
        joined.insert(joined.end(), outside_bridge.begin(), outside_bridge.end());
        sort(joined.begin(), joined.end());
        combined.push_back(std::move(joined));
    } else {
        vector<vector<int>> mb = map_unmarked_partition(b, bottom);
        vector<vector<int>> mt = map_unmarked_partition(t, top);
        if (mb.empty() || mt.empty()) return {};
        combined.insert(combined.end(), mb.begin(), mb.end());
        combined.insert(combined.end(), mt.begin(), mt.end());
    }
    sort(combined.begin(), combined.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    if ((int)combined.size() != expected) return {};
    if (!validate_partition(inst, combined, edge_offset)) return {};
    cerr << "STRICT_CLUSTER_CERT depth=" << depth
         << " B=" << B << " BM=" << BM
         << " T=" << T << " TM=" << TM
         << " components=" << combined.size()
         << " bridge=" << (bridge ? 1 : 0) << "\n";
    return remember(std::move(combined));
}

static string forest_stdout_string(const Instance &inst, const vector<vector<int>> &partition) {
    string out;
    for (auto &b : partition) {
        vector<unsigned char> active(inst.leaf_count + 1, 0);
        for (int x : b) active[x] = 1;
        out += restricted_newick(inst.trees[0], inst.trees[0].root, active);
        out += ";\n";
    }
    return out;
}

static void publish_signal_output(const Instance &inst, const vector<vector<int>> &partition,
        const vector<int> &edge_offset) {
    if (validate_partition(inst, partition, edge_offset)) {
        g_signal_output = forest_stdout_string(inst, partition);
    }
}

static void write_outputs(const Instance &inst, const vector<vector<int>> &partition, const string &csv, const string &nw) {
    if (!csv.empty()) {
        ofstream out(csv);
        for (auto &b : partition) {
            for (int i = 0; i < (int)b.size(); i++) { if (i) out << ','; out << b[i]; }
            out << '\n';
        }
    }
    if (!nw.empty()) {
        ofstream out(nw);
        for (auto &b : partition) {
            vector<unsigned char> active(inst.leaf_count + 1, 0);
            for (int x : b) active[x] = 1;
            out << restricted_newick(inst.trees[0], inst.trees[0].root, active) << ";\n";
        }
    }
    if (csv.empty() && nw.empty()) {
        cout << forest_stdout_string(inst, partition);
    }
}

static vector<vector<int>> read_component_csv(const string &path) {
    vector<vector<int>> rows;
    ifstream in(path);
    if (!in) throw runtime_error("could not open component csv");
    string line;
    while (getline(in, line)) {
        vector<int> block;
        string token;
        stringstream ss(line);
        while (getline(ss, token, ',')) {
            if (token.empty()) continue;
            block.push_back(stoi(token));
        }
        if (block.empty()) continue;
        sort(block.begin(), block.end());
        block.erase(unique(block.begin(), block.end()), block.end());
        rows.push_back(std::move(block));
    }
    sort(rows.begin(), rows.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    return rows;
}

struct ExchangeComponent {
    vector<int> leaves;
    vector<int> edges;
    vector<unsigned long long> leaf_mask;
    vector<unsigned long long> edge_mask;
    vector<int> touched;
};

struct ExchangeResult {
    bool found = false;
    bool exhausted = false;
    bool candidate_exhausted = false;
    bool region_exhausted = false;
    int candidates = 0;
    int regions = 0;
    long long dfs_nodes = 0;
    string reason;
    vector<int> region;
    vector<vector<int>> replacement;
};

static bool mask_subset(const vector<unsigned long long> &a,
        const vector<unsigned long long> &b) {
    for (int i = 0; i < (int)a.size(); i++) {
        if (a[i] & ~b[i]) return false;
    }
    return true;
}

static string region_key(const vector<int> &region) {
    string out;
    for (int x : region) {
        out += to_string(x);
        out += ",";
    }
    return out;
}

static ExchangeResult exchange_certify_partition(const Instance &inst,
        const Options &opt, const vector<int> &edge_offset,
        const vector<vector<int>> &incumbent) {
    ExchangeResult result;
    auto start = Clock::now();
    auto deadline = start + chrono::milliseconds((int)(opt.exchange_seconds * 1000));
    if (!validate_partition(inst, incumbent, edge_offset)) {
        result.reason = "invalid_incumbent";
        return result;
    }
    int n = inst.leaf_count;
    int leaf_words = (n + 63) >> 6;
    int edge_words = (edge_offset.back() + 63) >> 6;
    vector<vector<int>> incumbent_edges;
    vector<vector<unsigned long long>> incumbent_leaf_masks;
    vector<vector<unsigned long long>> incumbent_edge_masks;
    incumbent_edges.reserve(incumbent.size());
    incumbent_leaf_masks.reserve(incumbent.size());
    incumbent_edge_masks.reserve(incumbent.size());
    for (const vector<int> &block : incumbent) {
        vector<unsigned long long> leaf_mask(leaf_words, 0);
        for (int leaf : block) set_bit(leaf_mask, leaf - 1);
        vector<int> edges = component_edges_all(inst, block, edge_offset);
        vector<unsigned long long> edge_mask(edge_words, 0);
        for (int edge : edges) set_bit(edge_mask, edge);
        incumbent_leaf_masks.push_back(std::move(leaf_mask));
        incumbent_edges.push_back(std::move(edges));
        incumbent_edge_masks.push_back(std::move(edge_mask));
    }

    unordered_set<uint32_t> valid = collect_valid_triples(inst);
    unordered_set<string> candidate_seen;
    vector<ExchangeComponent> candidates;
    bool candidate_limited = false;
    vector<int> current;
    function<void(int)> enum_dfs = [&](int start_leaf) {
        if (candidate_limited || Clock::now() >= deadline) {
            candidate_limited = true;
            return;
        }
        if ((int)current.size() >= 2) {
            string key = key_of(current);
            if (candidate_seen.insert(key).second && valid_topology(inst, current)) {
                vector<int> edges = component_edges_all(inst, current, edge_offset);
                if ((int)edges.size() <= opt.exchange_max_component_edges) {
                    ExchangeComponent ec;
                    ec.leaves = current;
                    ec.edges = std::move(edges);
                    ec.leaf_mask.assign(leaf_words, 0);
                    ec.edge_mask.assign(edge_words, 0);
                    for (int leaf : ec.leaves) set_bit(ec.leaf_mask, leaf - 1);
                    for (int edge : ec.edges) set_bit(ec.edge_mask, edge);
                    for (int i = 0; i < (int)incumbent.size(); i++) {
                        if (intersects(ec.leaf_mask, incumbent_leaf_masks[i]) ||
                                intersects(ec.edge_mask, incumbent_edge_masks[i])) {
                            ec.touched.push_back(i);
                        }
                    }
                    if (!ec.touched.empty() &&
                            (int)ec.touched.size() <= opt.exchange_max_region_components) {
                        candidates.push_back(std::move(ec));
                        if ((int)candidates.size() >= opt.exchange_max_candidates) {
                            candidate_limited = true;
                            return;
                        }
                    }
                }
            }
        }
        if ((int)current.size() >= opt.exchange_max_component_size) return;
        for (int leaf = start_leaf; leaf <= n; leaf++) {
            current.push_back(leaf);
            if (triple_compatible(current, valid, n)) enum_dfs(leaf + 1);
            current.pop_back();
            if (candidate_limited) return;
        }
    };
    enum_dfs(1);
    result.candidates = (int)candidates.size();
    result.candidate_exhausted = !candidate_limited;
    cerr << "EXCHANGE_CANDIDATES count=" << candidates.size()
         << " exhausted=" << (result.candidate_exhausted ? 1 : 0)
         << " seconds=" << elapsed(start) << "\n";
    if (candidates.empty()) {
        result.exhausted = result.candidate_exhausted;
        result.region_exhausted = true;
        result.reason = result.exhausted ? "no_candidates" : "candidate_limit";
        return result;
    }

    vector<vector<int>> by_incumbent(incumbent.size());
    for (int i = 0; i < (int)candidates.size(); i++) {
        for (int idx : candidates[i].touched) by_incumbent[idx].push_back(i);
    }

    auto region_improves = [&](const vector<int> &region,
            vector<vector<int>> &replacement, long long &nodes) {
        if ((int)region.size() <= 1) return false;
        vector<unsigned char> in_region(incumbent.size(), 0);
        vector<unsigned long long> region_leaf_mask(leaf_words, 0);
        vector<unsigned long long> fixed_edge_mask(edge_words, 0);
        vector<int> region_leaves;
        vector<unsigned char> seen_leaf(n + 1, 0);
        for (int idx : region) {
            in_region[idx] = 1;
            bit_or_into(region_leaf_mask, incumbent_leaf_masks[idx]);
            for (int leaf : incumbent[idx]) {
                if (!seen_leaf[leaf]) {
                    seen_leaf[leaf] = 1;
                    region_leaves.push_back(leaf);
                }
            }
        }
        for (int i = 0; i < (int)incumbent.size(); i++) {
            if (!in_region[i]) bit_or_into(fixed_edge_mask, incumbent_edge_masks[i]);
        }
        sort(region_leaves.begin(), region_leaves.end());
        int target_savings = (int)region_leaves.size() - ((int)region.size() - 1);
        if (target_savings <= 0) return false;

        vector<ExchangeComponent> local;
        unordered_set<string> local_seen;
        auto add_local = [&](const vector<int> &leaves, const vector<int> &edges) {
            if (leaves.size() <= 1) return;
            string key = key_of(leaves);
            if (!local_seen.insert(key).second) return;
            ExchangeComponent ec;
            ec.leaves = leaves;
            ec.edges = edges;
            ec.leaf_mask.assign(leaf_words, 0);
            ec.edge_mask.assign(edge_words, 0);
            for (int leaf : leaves) set_bit(ec.leaf_mask, leaf - 1);
            for (int edge : edges) set_bit(ec.edge_mask, edge);
            if (!mask_subset(ec.leaf_mask, region_leaf_mask)) return;
            if (intersects(ec.edge_mask, fixed_edge_mask)) return;
            local.push_back(std::move(ec));
        };
        for (int idx : region) add_local(incumbent[idx], incumbent_edges[idx]);
        for (const ExchangeComponent &ec : candidates) {
            if (!mask_subset(ec.leaf_mask, region_leaf_mask)) continue;
            if (intersects(ec.edge_mask, fixed_edge_mask)) continue;
            add_local(ec.leaves, ec.edges);
        }
        if (local.empty()) return false;
        sort(local.begin(), local.end(), [](const ExchangeComponent &a,
                    const ExchangeComponent &b) {
            if (a.leaves.size() != b.leaves.size()) return a.leaves.size() > b.leaves.size();
            if (a.edges.size() != b.edges.size()) return a.edges.size() < b.edges.size();
            return a.leaves < b.leaves;
        });
        vector<vector<int>> by_leaf(n + 1);
        for (int i = 0; i < (int)local.size(); i++) {
            for (int leaf : local[i].leaves) by_leaf[leaf].push_back(i);
        }

        vector<unsigned long long> used_leaf(leaf_words, 0);
        vector<unsigned long long> used_edge(edge_words, 0);
        vector<int> chosen;
        bool timed_out = false;
        long long region_node_start = nodes;
        function<bool(int)> dfs = [&](int savings) {
            if (Clock::now() >= deadline ||
                    (opt.exchange_node_limit > 0 &&
                     nodes - region_node_start >= opt.exchange_node_limit)) {
                timed_out = true;
                return false;
            }
            nodes++;
            if (savings >= target_savings) return true;
            int assigned = bit_count(used_leaf);
            if (savings + max(0, (int)region_leaves.size() - assigned - 1) < target_savings) {
                return false;
            }
            int pivot = -1;
            int best_options = INT_MAX;
            for (int leaf : region_leaves) {
                if ((used_leaf[(leaf - 1) >> 6] >> ((leaf - 1) & 63)) & 1ULL) continue;
                int feasible = 0;
                for (int idx : by_leaf[leaf]) {
                    if (!intersects(used_leaf, local[idx].leaf_mask) &&
                            !intersects(used_edge, local[idx].edge_mask)) {
                        feasible++;
                    }
                }
                if (feasible < best_options) {
                    best_options = feasible;
                    pivot = leaf;
                    if (feasible == 0) break;
                }
            }
            if (pivot < 0) return false;
            vector<int> options;
            for (int idx : by_leaf[pivot]) {
                if (!intersects(used_leaf, local[idx].leaf_mask) &&
                        !intersects(used_edge, local[idx].edge_mask)) {
                    options.push_back(idx);
                }
            }
            sort(options.begin(), options.end(), [&](int a, int b) {
                int wa = (int)local[a].leaves.size() - 1;
                int wb = (int)local[b].leaves.size() - 1;
                if (wa != wb) return wa > wb;
                if (local[a].edges.size() != local[b].edges.size()) {
                    return local[a].edges.size() < local[b].edges.size();
                }
                return local[a].leaves < local[b].leaves;
            });
            for (int idx : options) {
                vector<unsigned long long> old_leaf = used_leaf;
                vector<unsigned long long> old_edge = used_edge;
                bit_or_into(used_leaf, local[idx].leaf_mask);
                bit_or_into(used_edge, local[idx].edge_mask);
                chosen.push_back(idx);
                if (dfs(savings + (int)local[idx].leaves.size() - 1)) return true;
                chosen.pop_back();
                used_leaf = std::move(old_leaf);
                used_edge = std::move(old_edge);
                if (timed_out) return false;
            }
            set_bit(used_leaf, pivot - 1);
            bool ok = dfs(savings);
            used_leaf[(pivot - 1) >> 6] &= ~(1ULL << ((pivot - 1) & 63));
            return ok;
        };
        bool hit = dfs(0);
        if (timed_out) {
            result.reason = "exchange_dfs_limit";
            return false;
        }
        if (!hit) return false;
        vector<unsigned char> covered(n + 1, 0);
        for (int idx : chosen) {
            replacement.push_back(local[idx].leaves);
            for (int leaf : local[idx].leaves) covered[leaf] = 1;
        }
        for (int leaf : region_leaves) {
            if (!covered[leaf]) replacement.push_back(vector<int>{leaf});
        }
        sort(replacement.begin(), replacement.end(), [](const vector<int> &a,
                    const vector<int> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        return true;
    };

    unordered_set<string> seen_regions;
    vector<vector<int>> queue;
    auto push_region = [&](vector<int> region) {
        sort(region.begin(), region.end());
        region.erase(unique(region.begin(), region.end()), region.end());
        if (region.empty() || (int)region.size() > opt.exchange_max_region_components) return;
        string key = region_key(region);
        if (seen_regions.insert(key).second) queue.push_back(std::move(region));
    };
    for (const ExchangeComponent &ec : candidates) push_region(ec.touched);

    bool skipped_hard_region = false;
    for (int qi = 0; qi < (int)queue.size(); qi++) {
        if (Clock::now() >= deadline) {
            result.reason = "timeout";
            break;
        }
        if (result.regions >= opt.exchange_max_regions) {
            result.reason = "region_limit";
            break;
        }
        vector<int> region = queue[qi];
        result.regions++;
        long long before_nodes = result.dfs_nodes;
        vector<vector<int>> replacement;
        if (region_improves(region, replacement, result.dfs_nodes)) {
            result.found = true;
            result.region = std::move(region);
            result.replacement = std::move(replacement);
            result.reason = "improving_exchange";
            return result;
        }
        if (!result.reason.empty() && result.reason == "exchange_dfs_limit") {
            skipped_hard_region = true;
            if (result.regions <= 10 || result.regions % 1000 == 0) {
                cerr << "EXCHANGE_REGION_SKIP regions=" << result.regions
                     << " queue=" << queue.size()
                     << " size=" << region.size()
                     << " nodes=" << (result.dfs_nodes - before_nodes)
                     << " reason=exchange_dfs_limit"
                     << " seconds=" << elapsed(start) << "\n";
            }
            result.reason.clear();
            continue;
        }
        if (result.regions <= 10 || result.regions % 1000 == 0) {
            cerr << "EXCHANGE_REGION_CHECK regions=" << result.regions
                 << " queue=" << queue.size()
                 << " size=" << region.size()
                 << " nodes=" << (result.dfs_nodes - before_nodes)
                 << " seconds=" << elapsed(start) << "\n";
        }
        vector<unsigned char> in_region(incumbent.size(), 0);
        for (int idx : region) in_region[idx] = 1;
        unordered_set<int> neighbor_candidates;
        for (int idx : region) {
            for (int cid : by_incumbent[idx]) neighbor_candidates.insert(cid);
        }
        for (int cid : neighbor_candidates) {
            vector<int> expanded = region;
            for (int idx : candidates[cid].touched) {
                if (!in_region[idx]) expanded.push_back(idx);
            }
            if ((int)expanded.size() <= opt.exchange_max_region_components) push_region(std::move(expanded));
            if ((int)queue.size() >= opt.exchange_max_regions) break;
        }
    }
    if (skipped_hard_region && result.reason.empty()) result.reason = "exchange_dfs_limit";
    result.region_exhausted = result.reason.empty();
    result.exhausted = result.candidate_exhausted && result.region_exhausted;
    if (result.reason.empty()) result.reason = result.exhausted ? "exhausted" : "not_exhausted";
    return result;
}

static ExchangeResult exchange_certify_incumbent(const Instance &inst,
        const Options &opt, const vector<int> &edge_offset) {
    ExchangeResult result;
    if (opt.exchange_incumbent.empty()) {
        result.reason = "missing_exchange_incumbent";
        return result;
    }
    vector<vector<int>> incumbent = read_component_csv(opt.exchange_incumbent);
    return exchange_certify_partition(inst, opt, edge_offset, incumbent);
}

static int stable_exchange_min_savings_for(const Instance &inst) {
    int conservative = max(12, inst.leaf_count / 5);
    if (inst.tree_count < 12) return conservative;
    return min(conservative, max(8, inst.leaf_count / 14));
}

static bool stable_exchange_accept_partition(const Instance &inst, const Options &opt,
        const vector<int> &edge_offset, const vector<vector<int>> &partition,
        const string &context) {
    if (!opt.allow_stable_exchange || opt.has_target) return false;
    if (opt.exchange_max_component_size < 8 || opt.exchange_max_region_components < 6) {
        cerr << "STABLE_EXCHANGE_SKIPPED context=" << context
             << " reason=weak_bounds"
             << " max_component_size=" << opt.exchange_max_component_size
             << " max_region_components=" << opt.exchange_max_region_components << "\n";
        return false;
    }
    int savings = inst.leaf_count - (int)partition.size();
    int min_savings = stable_exchange_min_savings_for(inst);
    if (savings < min_savings) {
        cerr << "STABLE_EXCHANGE_SKIPPED context=" << context
             << " reason=low_savings"
             << " savings=" << savings
             << " min_savings=" << min_savings << "\n";
        return false;
    }
    if (opt.exchange_seconds <= 0.0) {
        cerr << "STABLE_EXCHANGE_SKIPPED context=" << context
             << " reason=no_time\n";
        return false;
    }
    ExchangeResult ex = exchange_certify_partition(inst, opt, edge_offset, partition);
    cerr << "STABLE_EXCHANGE_RESULT context=" << context
         << " found=" << (ex.found ? 1 : 0)
         << " exhausted=" << (ex.exhausted ? 1 : 0)
         << " candidate_exhausted=" << (ex.candidate_exhausted ? 1 : 0)
         << " region_exhausted=" << (ex.region_exhausted ? 1 : 0)
         << " candidates=" << ex.candidates
         << " regions=" << ex.regions
         << " nodes=" << ex.dfs_nodes
         << " reason=" << ex.reason << "\n";
    if (!ex.found && ex.exhausted) {
        cerr << "STABLE_EXCHANGE_NOT_CERTIFICATE context=" << context
             << " reason=local_neighborhood_only\n";
    }
    return false;
}

static void report_checked_components(const Pool &pool, const string &path) {
    vector<vector<int>> rows = read_component_csv(path);
    unordered_map<int, vector<int>> by_leaf;
    for (int i = 0; i < (int)pool.comps.size(); i++) {
        const Component &c = pool.comps[i];
        if (c.leaves.size() <= 1 || c.leaves.size() > 8 || c.edges.size() > 120) continue;
        for (int leaf : c.leaves) by_leaf[leaf].push_back(i);
    }
    auto better = [&](int a, int b) {
        const Component &ca = pool.comps[a], &cb = pool.comps[b];
        double ra = (double)((int)ca.leaves.size() - 1) / max(1, (int)ca.edges.size());
        double rb = (double)((int)cb.leaves.size() - 1) / max(1, (int)cb.edges.size());
        if (ra != rb) return ra > rb;
        if (ca.leaves.size() != cb.leaves.size()) return ca.leaves.size() > cb.leaves.size();
        if (ca.edges.size() != cb.edges.size()) return ca.edges.size() < cb.edges.size();
        return ca.leaves < cb.leaves;
    };
    for (auto &kv : by_leaf) sort(kv.second.begin(), kv.second.end(), better);
    int need = 0, hit = 0;
    for (const vector<int> &block : rows) {
        if (block.size() <= 1) continue;
        need++;
        auto found = pool.id.find(key_of(block));
        bool ok = found != pool.id.end();
        if (ok) hit++;
        cerr << "CHECK_COMPONENT " << (ok ? "hit" : "miss") << " size=" << block.size() << " leaves=";
        for (int i = 0; i < (int)block.size(); i++) {
            if (i) cerr << ',';
            cerr << block[i];
        }
        if (ok) {
            const Component &c = pool.comps[found->second];
            cerr << " edges=" << c.edges.size();
            cerr << " ranks=";
            for (int i = 0; i < (int)block.size(); i++) {
                int leaf = block[i];
                int rank = -1;
                auto it = by_leaf.find(leaf);
                if (it != by_leaf.end()) {
                    for (int pos = 0; pos < (int)it->second.size(); pos++) {
                        if (it->second[pos] == found->second) {
                            rank = pos + 1;
                            break;
                        }
                    }
                }
                if (i) cerr << ',';
                cerr << leaf << ':' << rank;
            }
        }
        cerr << "\n";
    }
    cerr << "CHECK_COMPONENT_SUMMARY hit=" << hit << " need=" << need << "\n";
}

static int parse_int(const string &s) { return stoi(s); }
static double parse_double(const string &s) { return stod(s); }
static Options parse_options(int argc, char **argv) {
    Options opt;
    int i = 1;
    if (argc >= 2 && string(argv[1]).rfind("--", 0) != 0) {
        opt.instance_path = argv[1];
        i = 2;
    } else {
        opt.instance_path = "-";
    }
    for (; i < argc; i++) {
        string a = argv[i];
        auto need = [&]() -> string { if (i + 1 >= argc) throw runtime_error("missing value"); return argv[++i]; };
        if (a == "--target-components") { opt.target_components = parse_int(need()); opt.has_target = true; }
        else if (a == "--seed") opt.seed = parse_int(need());
        else if (a == "--max-pool") opt.max_pool = parse_int(need());
        else if (a == "--pool-seconds") opt.pool_seconds = parse_double(need());
        else if (a == "--check-components") opt.check_components = need();
        else if (a == "--dump-packing-lp") opt.dump_packing_lp = need();
        else if (a == "--dump-packing-blocks") opt.dump_packing_blocks = need();
        else if (a == "--dump-high-divergence-lp") opt.dump_high_divergence_lp = need();
        else if (a == "--dump-high-divergence-blocks") opt.dump_high_divergence_blocks = need();
        else if (a == "--output-components") opt.output_components = need();
        else if (a == "--output-newick") opt.output_newick = need();
        else if (a == "--swap-seconds") opt.swap_seconds = parse_double(need());
        else if (a == "--medium-swap-seconds") opt.medium_swap_seconds = parse_double(need());
        else if (a == "--medium-swap-candidate-limit") opt.medium_swap_candidate_limit = parse_int(need());
        else if (a == "--medium-swap-max-old-components") opt.medium_swap_max_old_components = parse_int(need());
        else if (a == "--medium-swap-max-replacement-size") opt.medium_swap_max_replacement_size = parse_int(need());
        else if (a == "--local-seconds") opt.local_seconds = parse_double(need());
        else if (a == "--local-guided-components") opt.local_guided_components = parse_int(need());
        else if (a == "--local-top-per-leaf") opt.local_top_per_leaf = parse_int(need());
        else if (a == "--local-trials") opt.local_trials = parse_int(need());
        else if (a == "--local-random-extra-components") opt.local_random_extra_components = parse_int(need());
        else if (a == "--local-random-anchor-components") opt.local_random_anchor_components = parse_int(need());
        else if (a == "--local-guided-growth-trials") opt.local_guided_growth_trials = parse_int(need());
        else if (a == "--local-guided-choice-window") opt.local_guided_choice_window = parse_int(need());
        else if (a == "--local-guided-future-probe") opt.local_guided_future_probe = parse_int(need());
        else if (a == "--local-max-pool") opt.local_max_pool = parse_int(need());
        else if (a == "--local-max-component-size") opt.local_max_component_size = parse_int(need());
        else if (a == "--local-seed-triples-limit") opt.local_seed_triples_limit = parse_int(need());
        else if (a == "--local-variant-seed-limit") opt.local_variant_seed_limit = parse_int(need());
        else if (a == "--local-column-max-size") opt.local_column_max_size = parse_int(need());
        else if (a == "--local-column-max-edges") opt.local_column_max_edges = parse_int(need());
        else if (a == "--local-round-seconds") opt.local_round_seconds = parse_double(need());
        else if (a == "--global-search-seconds") opt.global_search_seconds = parse_double(need());
        else if (a == "--global-greedy-starts") opt.global_greedy_starts = parse_int(need());
        else if (a == "--global-candidate-limit") opt.global_candidate_limit = parse_int(need());
        else if (a == "--global-remove-count") opt.global_remove_count = parse_int(need());
        else if (a == "--global-remove-jitter") opt.global_remove_jitter = parse_int(need());
        else if (a == "--global-highs-seconds") opt.global_highs_seconds = parse_double(need());
        else if (a == "--global-highs-top-per-leaf") opt.global_highs_top_per_leaf = parse_int(need());
        else if (a == "--global-highs-column-max-size") opt.global_highs_column_max_size = parse_int(need());
        else if (a == "--global-highs-column-max-edges") opt.global_highs_column_max_edges = parse_int(need());
        else if (a == "--global-branch-seconds") opt.global_branch_seconds = parse_double(need());
        else if (a == "--global-branch-top-per-leaf") opt.global_branch_top_per_leaf = parse_int(need());
        else if (a == "--global-branch-column-max-size") opt.global_branch_column_max_size = parse_int(need());
        else if (a == "--global-branch-column-max-edges") opt.global_branch_column_max_edges = parse_int(need());
        else if (a == "--global-branch-width") opt.global_branch_width = parse_int(need());
        else if (a == "--high-divergence-pack-seconds") opt.high_divergence_pack_seconds = parse_double(need());
        else if (a == "--high-divergence-max-component-size") opt.high_divergence_max_component_size = parse_int(need());
        else if (a == "--high-divergence-column-max-edges") opt.high_divergence_column_max_edges = parse_int(need());
        else if (a == "--high-divergence-max-components") opt.high_divergence_max_components = parse_int(need());
        else if (a == "--high-divergence-sat-candidate-limit") opt.high_divergence_sat_candidate_limit = parse_int(need());
        else if (a == "--high-divergence-top-per-leaf") opt.high_divergence_top_per_leaf = parse_int(need());
        else if (a == "--exact-exchange-seconds") opt.exact_exchange_seconds = parse_double(need());
        else if (a == "--exact-exchange-attempt-seconds") opt.exact_exchange_attempt_seconds = parse_double(need());
        else if (a == "--exact-exchange-max-region-components") opt.exact_exchange_max_region_components = parse_int(need());
        else if (a == "--exact-exchange-max-region-leaves") opt.exact_exchange_max_region_leaves = parse_int(need());
        else if (a == "--exact-exchange-max-candidates") opt.exact_exchange_max_candidates = parse_int(need());
        else if (a == "--exact-exchange-top-per-leaf") opt.exact_exchange_top_per_leaf = parse_int(need());
        else if (a == "--exact-exchange-column-max-size") opt.exact_exchange_column_max_size = parse_int(need());
        else if (a == "--exact-exchange-column-max-edges") opt.exact_exchange_column_max_edges = parse_int(need());
        else if (a == "--exact-exchange-seed-limit") opt.exact_exchange_seed_limit = parse_int(need());
        else if (a == "--pair-augment-seconds") opt.pair_augment_seconds = parse_double(need());
        else if (a == "--pair-augment-attempt-seconds") opt.pair_augment_attempt_seconds = parse_double(need());
        else if (a == "--pair-augment-candidate-limit") opt.pair_augment_candidate_limit = parse_int(need());
        else if (a == "--pair-augment-bundle-size") opt.pair_augment_bundle_size = parse_int(need());
        else if (a == "--pair-augment-max-region-components") opt.pair_augment_max_region_components = parse_int(need());
        else if (a == "--pair-augment-column-max-size") opt.pair_augment_column_max_size = parse_int(need());
        else if (a == "--pair-augment-column-max-edges") opt.pair_augment_column_max_edges = parse_int(need());
        else if (a == "--pair-augment-top-per-leaf") opt.pair_augment_top_per_leaf = parse_int(need());
        else if (a == "--incumbent-repack-seconds") opt.incumbent_repack_seconds = parse_double(need());
        else if (a == "--incumbent-repack-attempt-seconds") opt.incumbent_repack_attempt_seconds = parse_double(need());
        else if (a == "--incumbent-repack-candidate-limit") opt.incumbent_repack_candidate_limit = parse_int(need());
        else if (a == "--incumbent-repack-max-region-components") opt.incumbent_repack_max_region_components = parse_int(need());
        else if (a == "--incumbent-repack-column-max-size") opt.incumbent_repack_column_max_size = parse_int(need());
        else if (a == "--incumbent-repack-column-max-edges") opt.incumbent_repack_column_max_edges = parse_int(need());
        else if (a == "--incumbent-repack-top-per-leaf") opt.incumbent_repack_top_per_leaf = parse_int(need());
        else if (a == "--incumbent-repack-branch-width") opt.incumbent_repack_branch_width = parse_int(need());
        else if (a == "--small-savings-exact-seconds") opt.small_savings_exact_seconds = parse_double(need());
        else if (a == "--small-savings-pair-node-limit") opt.small_savings_pair_node_limit = stoll(need());
        else if (a == "--small-savings-special-node-limit") opt.small_savings_special_node_limit = stoll(need());
        else if (a == "--small-savings-component-limit") opt.small_savings_component_limit = parse_int(need());
        else if (a == "--small-savings-max-target") opt.small_savings_max_target = parse_int(need());
        else if (a == "--direct-pair-seconds") opt.direct_pair_seconds = parse_double(need());
        else if (a == "--direct-pair-leaf-limit") opt.direct_pair_leaf_limit = parse_int(need());
        else if (a == "--direct-lazy-preload-side") opt.direct_lazy_preload_side = parse_int(need());
        else if (a == "--exact-branch-seconds") opt.exact_branch_seconds = parse_double(need());
        else if (a == "--exact-branch-node-limit") opt.exact_branch_node_limit = stoll(need());
        else if (a == "--exact-pair-lb-seconds") opt.exact_pair_lb_seconds = parse_double(need());
        else if (a == "--exact-pair-lb-max-pairs") opt.exact_pair_lb_max_pairs = parse_int(need());
        else if (a == "--reference-cut-seconds") opt.reference_cut_seconds = parse_double(need());
        else if (a == "--savings-lb-seconds") opt.savings_lb_seconds = parse_double(need());
        else if (a == "--savings-lb-max-depth") opt.savings_lb_max_depth = parse_int(need());
        else if (a == "--savings-lb-max-components-per-triple") opt.savings_lb_max_components_per_triple = parse_int(need());
        else if (a == "--savings-lb-pair-node-limit") opt.savings_lb_pair_node_limit = stoll(need());
        else if (a == "--high-divergence-decision-seconds") opt.high_divergence_decision_seconds = parse_double(need());
        else if (a == "--seeded-pair-decision-seconds") opt.seeded_pair_decision_seconds = parse_double(need());
        else if (a == "--seeded-pair-decision-seed-limit") opt.seeded_pair_decision_seed_limit = parse_int(need());
        else if (a == "--exchange-incumbent") opt.exchange_incumbent = need();
        else if (a == "--exchange-seconds") opt.exchange_seconds = parse_double(need());
        else if (a == "--exchange-max-component-size") opt.exchange_max_component_size = parse_int(need());
        else if (a == "--exchange-max-component-edges") opt.exchange_max_component_edges = parse_int(need());
        else if (a == "--exchange-max-region-components") opt.exchange_max_region_components = parse_int(need());
        else if (a == "--exchange-max-candidates") opt.exchange_max_candidates = parse_int(need());
        else if (a == "--exchange-max-regions") opt.exchange_max_regions = parse_int(need());
        else if (a == "--exchange-node-limit") opt.exchange_node_limit = stoll(need());
        else if (a == "--total-seconds") opt.total_seconds = parse_double(need());
        else if (a == "--output-safety-seconds") opt.output_safety_seconds = parse_double(need());
        else if (a == "--tle-two-trees") opt.tle_two_trees = true;
        else if (a == "--no-tle-two-trees") opt.tle_two_trees = false;
        else if (a == "--auto-target-loop") opt.auto_target_loop = true;
        else if (a == "--no-auto-target-loop") opt.auto_target_loop = false;
        else if (a == "--strict-exact") opt.strict_exact = true;
        else if (a == "--no-strict-exact") opt.strict_exact = false;
        else if (a == "--allow-stable-exchange") opt.allow_stable_exchange = true;
        else if (a == "--no-allow-stable-exchange") opt.allow_stable_exchange = false;
        else if (a == "--allow-empirical-exact") opt.allow_empirical_exact = true;
        else if (a == "--no-allow-empirical-exact") opt.allow_empirical_exact = false;
        else if (a == "--full-sat-before-lazy") opt.full_sat_before_lazy = true;
        else if (a == "--no-full-sat-before-lazy") opt.full_sat_before_lazy = false;
        else if (a == "--exact-branch-decision") opt.exact_branch_decision = true;
        else if (a == "--exact-pair-lb-decision") opt.exact_pair_lb_decision = true;
        else if (a == "--direct-lazy-decision") opt.direct_lazy_decision = true;
        else if (a == "--reference-cut-decision") opt.reference_cut_decision = true;
        else if (a == "--savings-lb-decision") opt.savings_lb_decision = true;
        else if (a == "--edge-owner-decision") opt.edge_owner_decision = true;
        else if (a == "--high-divergence-pool-decision") opt.high_divergence_pool_decision = true;
        else if (a == "--seeded-pair-decomposition-decision") opt.seeded_pair_decomposition_decision = true;
        else if (a == "--exchange-certify") opt.exchange_certify = true;
        else throw runtime_error("unknown option " + a);
    }
    if (!opt.has_target) opt.target_components = 1;
    if (opt.target_components < 1) throw runtime_error("bad --target-components");
    return opt;
}

int main(int argc, char **argv) {
    try {
        signal(SIGTERM, sigterm_emit_best);
        signal(SIGINT, sigterm_emit_best);
        auto start = Clock::now();
        Options opt = parse_options(argc, argv);
        Instance inst = read_instance_path(opt.instance_path);
        if (opt.exchange_certify) {
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            ExchangeResult ex = exchange_certify_incumbent(inst, opt, edge_offset);
            cerr << "SUMMARY exchange_certify found=" << (ex.found ? 1 : 0)
                 << " exhausted=" << (ex.exhausted ? 1 : 0)
                 << " candidate_exhausted=" << (ex.candidate_exhausted ? 1 : 0)
                 << " region_exhausted=" << (ex.region_exhausted ? 1 : 0)
                 << " candidates=" << ex.candidates
                 << " regions=" << ex.regions
                 << " nodes=" << ex.dfs_nodes
                 << " reason=" << ex.reason
                 << " seconds=" << elapsed(start) << "\n";
            if (ex.found) {
                vector<vector<int>> incumbent = read_component_csv(opt.exchange_incumbent);
                vector<unsigned char> replace(incumbent.size(), 0);
                for (int idx : ex.region) replace[idx] = 1;
                vector<vector<int>> combined;
                for (int i = 0; i < (int)incumbent.size(); i++) {
                    if (!replace[i]) combined.push_back(incumbent[i]);
                }
                combined.insert(combined.end(), ex.replacement.begin(), ex.replacement.end());
                sort(combined.begin(), combined.end(), [](const vector<int> &a,
                            const vector<int> &b) {
                    if (a[0] != b[0]) return a[0] < b[0];
                    if (a.size() != b.size()) return a.size() < b.size();
                    return a < b;
                });
                if (!validate_partition(inst, combined, edge_offset)) {
                    cerr << "EXCHANGE_FOUND_INVALID_FULL_FOREST\n";
                    return 124;
                }
                cerr << "EXCHANGE_FOUND old_components=" << incumbent.size()
                     << " new_components=" << combined.size()
                     << " region_old=" << ex.region.size()
                     << " region_new=" << ex.replacement.size()
                     << " region=";
                for (int i = 0; i < (int)ex.region.size(); i++) {
                    if (i) cerr << ',';
                    cerr << ex.region[i];
                }
                cerr << "\n";
                publish_signal_output(inst, combined, edge_offset);
                write_outputs(inst, combined, opt.output_components, opt.output_newick);
                return 0;
            }
            return ex.exhausted ? 2 : 124;
        }
        if (opt.exact_pair_lb_decision) {
            double seconds = opt.exact_pair_lb_seconds > 0.0 ?
                opt.exact_pair_lb_seconds : opt.exact_branch_seconds;
            bool ruled_out = pairwise_branch_rule_out(
                inst, opt.target_components, seconds,
                opt.exact_branch_node_limit, opt.exact_pair_lb_max_pairs);
            cerr << "SUMMARY exact_pair_lb_decision target=" << opt.target_components
                 << " result=" << (ruled_out ? "NO" : "UNKNOWN")
                 << " seconds=" << elapsed(start) << "\n";
            return ruled_out ? 2 : 124;
        }
        if (opt.direct_lazy_decision) {
            bool proved_unsat = false;
            vector<vector<int>> partition = direct_pair_cadical_lazy_target(
                inst, opt.target_components, opt.direct_pair_seconds, &proved_unsat,
                opt.direct_lazy_preload_side);
            if (proved_unsat) {
                cerr << "SUMMARY direct_lazy_decision target=" << opt.target_components
                     << " result=NO seconds=" << elapsed(start) << "\n";
                return 2;
            }
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            if (!partition.empty() && validate_partition(inst, partition, edge_offset)) {
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                cerr << "SUMMARY direct_lazy_decision target=" << opt.target_components
                     << " result=YES components=" << partition.size()
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "SUMMARY direct_lazy_decision target=" << opt.target_components
                 << " result=UNKNOWN seconds=" << elapsed(start) << "\n";
            return 124;
        }
        if (opt.edge_owner_decision) {
            bool proved_unsat = false;
            vector<vector<int>> partition = edge_owner_cadical_target(
                inst, opt.target_components, opt.direct_pair_seconds, &proved_unsat);
            if (proved_unsat) {
                cerr << "SUMMARY edge_owner_decision target=" << opt.target_components
                     << " result=NO seconds=" << elapsed(start) << "\n";
                return 2;
            }
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            if (!partition.empty() && validate_partition(inst, partition, edge_offset)) {
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                cerr << "SUMMARY edge_owner_decision target=" << opt.target_components
                     << " result=YES components=" << partition.size()
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "SUMMARY edge_owner_decision target=" << opt.target_components
                 << " result=UNKNOWN seconds=" << elapsed(start) << "\n";
            return 124;
        }
        if (opt.reference_cut_decision) {
            bool proved_unsat = false;
            double seconds = opt.reference_cut_seconds > 0.0 ?
                opt.reference_cut_seconds : opt.direct_pair_seconds;
            vector<vector<int>> partition = reference_cut_cadical_target(
                inst, opt.target_components, seconds, &proved_unsat);
            if (proved_unsat) {
                cerr << "SUMMARY reference_cut_decision target=" << opt.target_components
                     << " result=NO seconds=" << elapsed(start) << "\n";
                return 2;
            }
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            if (!partition.empty() && validate_partition(inst, partition, edge_offset)) {
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                cerr << "SUMMARY reference_cut_decision target=" << opt.target_components
                     << " result=YES components=" << partition.size()
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "SUMMARY reference_cut_decision target=" << opt.target_components
                 << " result=UNKNOWN seconds=" << elapsed(start) << "\n";
            return 124;
        }
        if (opt.savings_lb_decision) {
            double seconds = opt.savings_lb_seconds > 0.0 ?
                opt.savings_lb_seconds : opt.direct_pair_seconds;
            SavingsLbResult lb = savings_lb_rule_out(
                inst, opt.target_components, seconds,
                opt.savings_lb_max_depth,
                opt.savings_lb_max_components_per_triple,
                opt.savings_lb_pair_node_limit);
            cerr << "SUMMARY savings_lb_decision target=" << opt.target_components
                 << " result=" << (lb.proved ? "NO" : "UNKNOWN")
                 << " component_lb=" << lb.component_lb
                 << " seconds=" << elapsed(start) << "\n";
            return lb.proved ? 2 : 124;
        }
        if (opt.high_divergence_pool_decision) {
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            double seconds = opt.high_divergence_decision_seconds > 0.0 ?
                opt.high_divergence_decision_seconds : opt.direct_pair_seconds;
            vector<vector<int>> partition;
            ExactDecision ans = high_divergence_pool_cadical_decide(
                inst, opt.target_components, edge_offset, opt, seconds, &partition);
            if (ans == ExactDecision::Yes && !partition.empty()) {
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
            }
            cerr << "SUMMARY high_divergence_pool_decision target=" << opt.target_components
                 << " result=" << (ans == ExactDecision::No ? "NO" :
                         ans == ExactDecision::Yes ? "YES" : "UNKNOWN")
                 << " components=" << (partition.empty() ? 0 : (int)partition.size())
                 << " seconds=" << elapsed(start) << "\n";
            return ans == ExactDecision::Yes ? 0 : ans == ExactDecision::No ? 2 : 124;
        }
        if (opt.seeded_pair_decomposition_decision) {
            vector<int> edge_offset(inst.tree_count + 1, 0);
            for (int t = 0; t < inst.tree_count; t++) {
                edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
            }
            double seconds = opt.seeded_pair_decision_seconds > 0.0 ?
                opt.seeded_pair_decision_seconds :
                opt.high_divergence_decision_seconds > 0.0 ?
                opt.high_divergence_decision_seconds : opt.direct_pair_seconds;
            vector<vector<int>> partition;
            ExactDecision ans = seeded_pair_decomposition_decide(
                inst, opt.target_components, edge_offset, opt, seconds, &partition);
            if (ans == ExactDecision::Yes && !partition.empty()) {
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
            }
            cerr << "SUMMARY seeded_pair_decomposition_decision target=" << opt.target_components
                 << " result=" << (ans == ExactDecision::No ? "NO" :
                         ans == ExactDecision::Yes ? "YES" : "UNKNOWN")
                 << " components=" << (partition.empty() ? 0 : (int)partition.size())
                 << " seconds=" << elapsed(start) << "\n";
            return ans == ExactDecision::Yes ? 0 : ans == ExactDecision::No ? 2 : 124;
        }
        if (opt.exact_branch_decision) {
            double seconds = opt.exact_branch_seconds > 0.0 ? opt.exact_branch_seconds : opt.direct_pair_seconds;
            ExactDecision ans = exact_branch_decide(
                inst, opt.target_components, seconds, opt.exact_branch_node_limit);
            cerr << "SUMMARY exact_branch_decision target=" << opt.target_components
                 << " result=" << (ans == ExactDecision::No ? "NO" : ans == ExactDecision::Yes ? "YES" : "UNKNOWN")
                 << " seconds=" << elapsed(start) << "\n";
            return ans == ExactDecision::Yes ? 0 : ans == ExactDecision::No ? 2 : 124;
        }
        if (opt.tle_two_trees && inst.tree_count == 2) {
            cerr << "TLE_TWO_TREES tree_count=2\n";
            while (true) this_thread::sleep_for(chrono::hours(24));
        }
        bool has_moderate_recursive_cluster = false;
        if (opt.strict_exact && !opt.has_target &&
                inst.leaf_count > opt.direct_pair_leaf_limit) {
            vector<int> initial_cluster = choose_recursive_common_cluster(inst);
            if (!initial_cluster.empty()) {
                if (inst.leaf_count <= 128) {
                    has_moderate_recursive_cluster = true;
                } else {
                    int outside = inst.leaf_count - (int)initial_cluster.size();
                    has_moderate_recursive_cluster =
                        (int)initial_cluster.size() >= 30 && outside >= 30;
                }
            }
        }
        bool prefer_incumbent_before_recursive =
            opt.auto_target_loop && opt.global_branch_seconds > 0.0 &&
            opt.pair_augment_seconds > 0.0 && inst.leaf_count <= 128;
        bool defer_strict_recursive_exact =
            opt.strict_exact && !opt.has_target &&
            opt.reference_cut_seconds > 0.0 &&
            (!has_moderate_recursive_cluster ||
             prefer_incumbent_before_recursive);
        bool prioritize_recursive_before_exchange =
            defer_strict_recursive_exact &&
            has_moderate_recursive_cluster &&
            inst.tree_count >= 12;
        if (opt.strict_exact && !opt.has_target && !defer_strict_recursive_exact) {
            vector<vector<int>> exact = strict_recursive_exact(inst, opt, 0);
            if (!exact.empty()) {
                vector<int> exact_edge_offset = recursive_edge_offset(inst);
                if (!validate_partition(inst, exact, exact_edge_offset)) {
                    throw runtime_error("recursive exact partition failed validation");
                }
                publish_signal_output(inst, exact, exact_edge_offset);
                write_outputs(inst, exact, opt.output_components, opt.output_newick);
                cerr << "SUMMARY components=" << exact.size()
                     << " certified=1"
                     << " method=strict_recursive_exact"
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "STRICT_RECURSIVE_EXACT_NO_CERTIFICATE continue=1\n";
        }
        vector<int> edge_offset(inst.tree_count + 1, 0);
        for (int t = 0; t < inst.tree_count; t++) edge_offset[t + 1] = edge_offset[t] + (int)inst.trees[t].nodes.size();
        auto phase_seconds = [&](double requested) {
            if (requested <= 0.0) return 0.0;
            if (opt.total_seconds <= 0.0) return requested;
            double remain = opt.total_seconds - opt.output_safety_seconds - elapsed(start);
            return max(0.0, min(requested, remain));
        };
        auto near_deadline = [&]() {
            return opt.total_seconds > 0.0 &&
                elapsed(start) >= opt.total_seconds - opt.output_safety_seconds;
        };
        auto stable_min_savings = [&]() {
            return stable_exchange_min_savings_for(inst);
        };
        bool user_requested_target = opt.has_target;
        bool stable_exchange_used = false;
        bool empirical_exact_used = false;
        int no_target_successful_rounds = 0;
        int no_target_reference_cut_attempted_target = -1;
        auto stable_exchange_accept = [&](const vector<vector<int>> &partition) {
            if (!opt.allow_stable_exchange || opt.has_target) return false;
            if (opt.exchange_max_component_size < 8 || opt.exchange_max_region_components < 6) {
                cerr << "STABLE_EXCHANGE_SKIPPED reason=weak_bounds"
                     << " max_component_size=" << opt.exchange_max_component_size
                     << " max_region_components=" << opt.exchange_max_region_components << "\n";
                return false;
            }
            int savings = inst.leaf_count - (int)partition.size();
            int min_savings = stable_min_savings();
            if (savings < min_savings) {
                cerr << "STABLE_EXCHANGE_SKIPPED reason=low_savings"
                     << " savings=" << savings
                     << " min_savings=" << min_savings << "\n";
                return false;
            }
            double seconds = phase_seconds(opt.exchange_seconds);
            if (seconds <= 0.0) {
                cerr << "STABLE_EXCHANGE_SKIPPED reason=no_time\n";
                return false;
            }
            Options stable_opt = opt;
            stable_opt.exchange_seconds = seconds;
            ExchangeResult ex = exchange_certify_partition(
                inst, stable_opt, edge_offset, partition);
            cerr << "STABLE_EXCHANGE_RESULT found=" << (ex.found ? 1 : 0)
                 << " exhausted=" << (ex.exhausted ? 1 : 0)
                 << " candidate_exhausted=" << (ex.candidate_exhausted ? 1 : 0)
                 << " region_exhausted=" << (ex.region_exhausted ? 1 : 0)
                 << " candidates=" << ex.candidates
                 << " regions=" << ex.regions
                 << " nodes=" << ex.dfs_nodes
                 << " reason=" << ex.reason << "\n";
            if (!ex.found && ex.exhausted) {
                cerr << "STABLE_EXCHANGE_NOT_CERTIFICATE reason=local_neighborhood_only\n";
            }
            return false;
        };
        auto empirical_exact_accept = [&](const vector<vector<int>> &partition,
                int lower_target, const string &context) {
            if (!opt.allow_empirical_exact) return false;
            if (!validate_partition(inst, partition, edge_offset)) return false;
            int k = (int)partition.size();
            if (opt.has_target && k > opt.target_components) return false;
            int savings = inst.leaf_count - k;
            int min_savings = stable_min_savings();
            int empirical_min_savings = min_savings + 1;
            bool trusted_midband = false;
            if (!(savings <= 10 || savings >= 20)) {
                trusted_midband =
                    !opt.has_target &&
                    context == "leaf_limit" &&
                    inst.tree_count >= 8 &&
                    savings >= 13 && savings <= 15 &&
                    no_target_successful_rounds >= 3;
                if (!trusted_midband) {
                    cerr << "EMPIRICAL_EXACT_REJECT reason=mid_savings_band"
                         << " savings=" << savings
                         << " lower_target=" << lower_target
                         << " context=" << context
                         << " successful_rounds=" << no_target_successful_rounds
                         << " t=" << inst.tree_count
                         << " n=" << inst.leaf_count << "\n";
                    return false;
                }
                cerr << "EMPIRICAL_EXACT_MIDBAND_ALLOWED"
                     << " savings=" << savings
                     << " lower_target=" << lower_target
                     << " context=" << context
                     << " successful_rounds=" << no_target_successful_rounds
                     << " t=" << inst.tree_count
                     << " n=" << inst.leaf_count << "\n";
            }
            if (!trusted_midband && savings < empirical_min_savings) {
                cerr << "EMPIRICAL_EXACT_REJECT reason=low_savings"
                     << " savings=" << savings
                     << " min_savings=" << empirical_min_savings
                     << " base_min_savings=" << min_savings
                     << " lower_target=" << lower_target
                     << " context=" << context << "\n";
                return false;
            }
            empirical_exact_used = true;
            cerr << "EMPIRICAL_EXACT_ACCEPT components=" << k
                 << " lower_target=" << lower_target
                 << " savings=" << savings
                 << " min_savings=" << min_savings
                 << " context=" << context
                 << " note=validated_forest_no_better_found\n";
            return true;
        };
        auto apply_exchange_result = [&](const vector<vector<int>> &incumbent,
                const ExchangeResult &ex) {
            vector<unsigned char> replace(incumbent.size(), 0);
            for (int idx : ex.region) replace[idx] = 1;
            vector<vector<int>> combined;
            for (int i = 0; i < (int)incumbent.size(); i++) {
                if (!replace[i]) combined.push_back(incumbent[i]);
            }
            combined.insert(combined.end(), ex.replacement.begin(), ex.replacement.end());
            sort(combined.begin(), combined.end(), [](const vector<int> &a,
                        const vector<int> &b) {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            if (!validate_partition(inst, combined, edge_offset)) {
                throw runtime_error("exchange improvement produced invalid partition");
            }
            return combined;
        };
        auto certify_exact_or_timeout = [&](const vector<vector<int>> &partition, const string &reason) {
            if (!opt.strict_exact) return;
            int k = (int)partition.size();
            if (opt.has_target && k > opt.target_components) wait_for_timeout(reason);
            int lower_target = k - 1;
            if (lower_target < 1) {
                cerr << "STRICT_CERTIFIED components=" << k << " lower_target=0 method=trivial\n";
                return;
            }
            bool ruled_out = false;
            string proof_method;
            if (opt.exact_pair_lb_seconds > 0.0) {
                ruled_out = pairwise_branch_rule_out(
                    inst, lower_target, phase_seconds(opt.exact_pair_lb_seconds),
                    opt.exact_branch_node_limit, opt.exact_pair_lb_max_pairs);
                if (ruled_out) proof_method = "pairwise_branch_lb";
            }
            if (!ruled_out && opt.reference_cut_seconds > 0.0 &&
                    lower_target == no_target_reference_cut_attempted_target) {
                cerr << "STRICT_CERTIFY_SKIP_REFERENCE_CUT_REPEAT lower_target="
                     << lower_target << "\n";
            } else if (!ruled_out && opt.reference_cut_seconds > 0.0) {
                bool cut_unsat = false;
                vector<vector<int>> cut_partition = reference_cut_cadical_target(
                    inst, lower_target, phase_seconds(opt.reference_cut_seconds), &cut_unsat);
                ruled_out = cut_unsat;
                if (ruled_out) proof_method = "reference_cut_infeasible";
                if (!ruled_out && !cut_partition.empty() &&
                        validate_partition(inst, cut_partition, edge_offset) &&
                        (int)cut_partition.size() < k) {
                    cerr << "STRICT_CERTIFY_IMPROVED incumbent=" << k
                         << " improved=" << cut_partition.size()
                         << " method=reference_cut_target\n";
                    wait_for_timeout("reference_cut_found_better_partition");
                }
            }
            if (!ruled_out) {
                int target_savings = inst.leaf_count - lower_target;
                if (target_savings > 0 && target_savings <= opt.small_savings_max_target) {
                    double small_seconds = phase_seconds(opt.small_savings_exact_seconds);
                    if (small_seconds > 0.0) {
                        SmallSavingsExactResult exact_small = small_savings_exhaustive_pack(
                            inst, lower_target, edge_offset, small_seconds,
                            opt.small_savings_pair_node_limit,
                            opt.small_savings_special_node_limit,
                            opt.small_savings_max_target,
                            opt.small_savings_component_limit);
                        if (exact_small.proved_infeasible) {
                            ruled_out = true;
                            proof_method = "small_savings_exact_lb";
                        } else if (!exact_small.partition.empty() &&
                                validate_partition(inst, exact_small.partition, edge_offset) &&
                                (int)exact_small.partition.size() < k) {
                            cerr << "STRICT_CERTIFY_IMPROVED incumbent=" << k
                                 << " improved=" << exact_small.partition.size()
                                 << " method=small_savings_exact_target\n";
                            wait_for_timeout("small_savings_found_better_partition");
                        }
                    }
                }
            }
            if (!ruled_out && opt.savings_lb_seconds > 0.0) {
                SavingsLbResult lb = savings_lb_rule_out(
                    inst, lower_target, phase_seconds(opt.savings_lb_seconds),
                    opt.savings_lb_max_depth,
                    opt.savings_lb_max_components_per_triple,
                    opt.savings_lb_pair_node_limit);
                ruled_out = lb.proved;
                if (ruled_out) proof_method = "savings_lb_infeasible";
            }
            if (!ruled_out && opt.exact_branch_seconds > 0.0) {
                double branch_seconds = phase_seconds(opt.exact_branch_seconds);
                if (branch_seconds > 0.0) {
                    ExactDecision branch = exact_branch_decide(
                        inst, lower_target, branch_seconds, opt.exact_branch_node_limit);
                    ruled_out = branch == ExactDecision::No;
                    if (ruled_out) proof_method = "exact_branch_infeasible";
                }
            }
            if (!ruled_out && opt.high_divergence_decision_seconds > 0.0) {
                double hd_seconds = phase_seconds(opt.high_divergence_decision_seconds);
                if (hd_seconds > 0.0) {
                    ExactDecision hd = high_divergence_pool_cadical_decide(
                        inst, lower_target, edge_offset, opt, hd_seconds, nullptr);
                    ruled_out = hd == ExactDecision::No;
                    if (ruled_out) proof_method = "high_divergence_pool_infeasible";
                }
            }
            if (!ruled_out && opt.seeded_pair_decision_seconds > 0.0) {
                double seeded_seconds = phase_seconds(opt.seeded_pair_decision_seconds);
                if (seeded_seconds > 0.0) {
                    vector<vector<int>> seeded_partition;
                    ExactDecision seeded = seeded_pair_decomposition_decide(
                        inst, lower_target, edge_offset, opt, seeded_seconds,
                        &seeded_partition);
                    if (seeded == ExactDecision::Yes && !seeded_partition.empty() &&
                            validate_partition(inst, seeded_partition, edge_offset) &&
                            (int)seeded_partition.size() < k) {
                        cerr << "STRICT_CERTIFY_IMPROVED incumbent=" << k
                             << " improved=" << seeded_partition.size()
                             << " method=seeded_pair_decomposition_target\n";
                        wait_for_timeout("seeded_pair_found_better_partition");
                    }
                    if (seeded == ExactDecision::No) {
                        cerr << "SEEDED_PAIR_DECOMPOSITION_NO_BETTER lower_target="
                             << lower_target
                             << " note=structural_class_only\n";
                    }
                }
            }
            if (!ruled_out && inst.leaf_count > opt.direct_pair_leaf_limit) {
                cerr << "STRICT_CERTIFY_SKIPPED reason=leaf_limit n=" << inst.leaf_count
                     << " limit=" << opt.direct_pair_leaf_limit
                     << " lower_target=" << lower_target << "\n";
                if (stable_exchange_accept(partition)) {
                    cerr << "STRICT_STABLE_ACCEPT components=" << k
                         << " lower_target=" << lower_target
                         << " method=stable_exchange\n";
                    return;
                }
                if (empirical_exact_accept(partition, lower_target, "leaf_limit")) return;
                wait_for_timeout("no_lower_bound_certificate");
            }
            double proof_seconds = phase_seconds(opt.direct_pair_seconds);
            if (!ruled_out && proof_seconds <= 0.0) {
                if (stable_exchange_accept(partition)) {
                    cerr << "STRICT_STABLE_ACCEPT components=" << k
                         << " lower_target=" << lower_target
                         << " method=stable_exchange\n";
                    return;
                }
                if (empirical_exact_accept(partition, lower_target, "no_time")) return;
                wait_for_timeout("no_time_for_lower_bound_certificate");
            }
            if (!ruled_out) {
                ruled_out = direct_pair_cadical_rule_out(inst, lower_target, proof_seconds);
                if (ruled_out) proof_method = "direct_pair_infeasible";
            }
            if (!ruled_out) {
                ruled_out = direct_pair_highs_rule_out(inst, lower_target, phase_seconds(opt.direct_pair_seconds), edge_offset);
                if (ruled_out) proof_method = "direct_pair_highs_infeasible";
            }
            if (!ruled_out) {
                vector<vector<int>> optimal = direct_pair_highs_optimal(
                    inst, phase_seconds(opt.direct_pair_seconds), edge_offset);
                if (!optimal.empty() && validate_partition(inst, optimal, edge_offset) &&
                        (int)optimal.size() == k) {
                    cerr << "STRICT_CERTIFIED components=" << k
                         << " lower_target=" << lower_target
                         << " method=direct_pair_optimal_fallback\n";
                    return;
                }
                cerr << "STRICT_CERTIFY_FAILED lower_target=" << lower_target
                     << " incumbent=" << k << "\n";
                if (stable_exchange_accept(partition)) {
                    cerr << "STRICT_STABLE_ACCEPT components=" << k
                         << " lower_target=" << lower_target
                         << " method=stable_exchange\n";
                    return;
                }
                if (empirical_exact_accept(partition, lower_target, "strict_certify_failed")) return;
                wait_for_timeout("lower_target_not_proved_impossible");
            }
            cerr << "STRICT_CERTIFIED components=" << k
                 << " lower_target=" << lower_target
                 << " method=" << proof_method << "\n";
        };
        auto finish_partition = [&](vector<vector<int>> &partition, const string &reason) {
            partition = singleton_merge_repair(inst, std::move(partition), opt.target_components, edge_offset);
            {
                int pair_target = opt.has_target ? opt.target_components :
                    max(1, (int)partition.size() - 1);
                vector<vector<int>> pair_repaired = pair_packing_repair(
                    inst, partition, pair_target, edge_offset);
                if (!pair_repaired.empty()) {
                    cerr << "USE_PAIR_PACK_REPAIR components=" << pair_repaired.size() << "\n";
                    partition = std::move(pair_repaired);
                }
            }
            if (!validate_partition(inst, partition, edge_offset)) throw runtime_error("internal invalid partition");
            publish_signal_output(inst, partition, edge_offset);
            certify_exact_or_timeout(partition, reason);
            write_outputs(inst, partition, opt.output_components, opt.output_newick);
            cerr << "SUMMARY components=" << partition.size()
                 << " target=" << opt.target_components
                 << " certified=" << (opt.strict_exact && !stable_exchange_used && !empirical_exact_used ? 1 : 0)
                 << " stable=" << (stable_exchange_used ? 1 : 0)
                 << " empirical=" << (empirical_exact_used ? 1 : 0)
                 << " reason=" << reason
                 << " seconds=" << elapsed(start) << "\n";
            return (!user_requested_target || !opt.has_target ||
                    partition.size() <= (size_t)opt.target_components) ? 0 : 2;
        };
        opt.direct_pair_seconds = phase_seconds(opt.direct_pair_seconds);
        if (opt.strict_exact && !opt.has_target) {
            if (inst.leaf_count <= opt.direct_pair_leaf_limit && opt.direct_pair_seconds > 0.0) {
                vector<vector<int>> direct_partition(inst.leaf_count);
                for (int leaf = 1; leaf <= inst.leaf_count; leaf++) {
                    direct_partition[leaf - 1] = vector<int>{leaf};
                }
                while ((int)direct_partition.size() > 1) {
                    int target = (int)direct_partition.size() - 1;
                    double seconds = phase_seconds(opt.direct_pair_seconds);
                    if (seconds <= 0.0) break;
                    bool proved_unsat = false;
                    vector<vector<int>> direct = direct_pair_cadical_target(
                        inst, target, seconds, &proved_unsat);
#ifndef USE_CADICAL
                    if (direct.empty()) {
                        direct = direct_pair_highs_target(
                            inst, target, seconds, edge_offset);
                    }
#endif
                    if (proved_unsat) {
                        if (!validate_partition(inst, direct_partition, edge_offset)) {
                            throw runtime_error("direct descent partition failed validation");
                        }
                        publish_signal_output(inst, direct_partition, edge_offset);
                        write_outputs(inst, direct_partition, opt.output_components, opt.output_newick);
                        cerr << "SUMMARY components=" << direct_partition.size()
                             << " certified=1"
                             << " method=direct_pair_prepool_descent"
                             << " lower_target=" << target
                             << " seconds=" << elapsed(start) << "\n";
                        return 0;
                    }
                    if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                            (int)direct.size() < (int)direct_partition.size()) {
                        cerr << "DIRECT_PAIR_PREPOOL_DESCENT from=" << direct_partition.size()
                             << " to=" << direct.size()
                             << " target=" << target << "\n";
                        direct_partition = std::move(direct);
                        continue;
                    }
                    cerr << "DIRECT_PAIR_PREPOOL_NO_CERTIFICATE current="
                         << direct_partition.size()
                         << " target=" << target << "\n";
                    break;
                }
            }
#ifndef USE_CADICAL
            if (inst.leaf_count <= opt.direct_pair_leaf_limit && opt.direct_pair_seconds > 0.0) {
                vector<vector<int>> direct = direct_pair_highs_optimal(
                    inst, opt.direct_pair_seconds, edge_offset);
                if (!direct.empty() && validate_partition(inst, direct, edge_offset)) {
                    publish_signal_output(inst, direct, edge_offset);
                    write_outputs(inst, direct, opt.output_components, opt.output_newick);
                    cerr << "SUMMARY components=" << direct.size()
                         << " certified=1"
                         << " method=direct_pair_optimal"
                         << " seconds=" << elapsed(start) << "\n";
                    return 0;
                }
                cerr << "DIRECT_PAIR_OPTIMAL_NO_CERTIFICATE n=" << inst.leaf_count << "\n";
            } else {
                cerr << "DIRECT_PAIR_OPTIMAL_SKIPPED n=" << inst.leaf_count
                     << " limit=" << opt.direct_pair_leaf_limit
                     << " seconds=" << opt.direct_pair_seconds << "\n";
            }
#endif
        }
        if (opt.has_target && inst.leaf_count <= opt.direct_pair_leaf_limit && opt.direct_pair_seconds > 0.0) {
            vector<vector<int>> direct = direct_pair_cadical_target(
                inst, opt.target_components, opt.direct_pair_seconds);
#ifndef USE_CADICAL
            if (direct.empty()) {
                direct = direct_pair_highs_target(
                    inst, opt.target_components, opt.direct_pair_seconds, edge_offset);
            }
#endif
            if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                    (int)direct.size() <= opt.target_components) {
                certify_exact_or_timeout(direct, "direct_pair_target_without_lower_certificate");
                publish_signal_output(inst, direct, edge_offset);
                write_outputs(inst, direct, opt.output_components, opt.output_newick);
                cerr << "SUMMARY components=" << direct.size()
                     << " target=" << opt.target_components
                     << " certified=" << (opt.strict_exact ? 1 : 0)
                     << " method=direct_pair"
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "DIRECT_PAIR_EARLY_NO_TARGET target=" << opt.target_components << "\n";
        }
        if (opt.strict_exact && !opt.has_target &&
                (inst.tree_count >= 10 || inst.leaf_count > 128)) {
            int strict_pool_cap = inst.tree_count >= 20 ? 400000 : 700000;
            if (inst.leaf_count >= 130) strict_pool_cap = min(strict_pool_cap, 550000);
            opt.max_pool = min(opt.max_pool, strict_pool_cap);
            opt.pack_pool_limit = min(opt.pack_pool_limit, strict_pool_cap);
        }
        opt.pool_seconds = phase_seconds(opt.pool_seconds);
        Pool pool = build_pool(inst, opt, edge_offset);
        if (!opt.check_components.empty()) report_checked_components(pool, opt.check_components);
        int total_bits = inst.leaf_count + edge_offset.back();
        vector<Column> cols = encode_columns(inst, pool, total_bits);
        if (!opt.dump_packing_lp.empty()) {
            dump_packing_lp_files(cols, opt);
            return 0;
        }
        vector<int> selected = anchored_pack(cols, inst.leaf_count, opt);
        opt.global_search_seconds = phase_seconds(opt.global_search_seconds);
        vector<int> global_selected = global_lns_pack(cols, inst.leaf_count, opt);
        if (!global_selected.empty() &&
                component_count(cols, global_selected, inst.leaf_count) < component_count(cols, selected, inst.leaf_count)) {
            selected = std::move(global_selected);
            cerr << "USE_GLOBAL_LNS components=" << component_count(cols, selected, inst.leaf_count) << "\n";
        }
        opt.swap_seconds = phase_seconds(opt.swap_seconds);
        selected = swap_repair(cols, selected, inst.leaf_count, opt);
        opt.medium_swap_seconds = phase_seconds(opt.medium_swap_seconds);
        selected = medium_swap_repair(cols, selected, inst.leaf_count, opt);
        vector<vector<int>> partition = selected_partition(cols, selected, inst.leaf_count);
        partition = singleton_merge_repair(inst, std::move(partition), opt.target_components, edge_offset);
        if (user_requested_target && !opt.strict_exact &&
                (int)partition.size() <= opt.target_components) {
            return finish_partition(partition, "target_hit_initial_pack");
        }
        {
            int pair_target = opt.has_target ? opt.target_components :
                max(1, (int)partition.size() - 1);
            vector<vector<int>> pair_repaired = pair_packing_repair(
                inst, partition, pair_target, edge_offset);
            if (!pair_repaired.empty()) {
                cerr << "USE_PAIR_PACK_REPAIR components=" << pair_repaired.size() << "\n";
                partition = std::move(pair_repaired);
            }
            vector<vector<int>> small_repaired = small_savings_pack_repair(
                inst, cols, pair_target, edge_offset);
                if (!small_repaired.empty() &&
                        (int)small_repaired.size() < (int)partition.size()) {
                    cerr << "USE_SMALL_SAVINGS_PACK components=" << small_repaired.size() << "\n";
                    partition = std::move(small_repaired);
                    if (user_requested_target && !opt.strict_exact &&
                            (int)partition.size() <= opt.target_components) {
                        return finish_partition(partition, "target_hit_small_savings_pack");
                    }
                }
            if (opt.has_target) {
                SmallSavingsExactResult exact_small = small_savings_exhaustive_pack(
                    inst, pair_target, edge_offset, phase_seconds(opt.small_savings_exact_seconds),
                    opt.small_savings_pair_node_limit,
                    opt.small_savings_special_node_limit,
                    opt.small_savings_max_target,
                    opt.small_savings_component_limit);
                if (!exact_small.partition.empty() &&
                        (int)exact_small.partition.size() < (int)partition.size()) {
                    cerr << "USE_SMALL_SAVINGS_EXACT components="
                         << exact_small.partition.size() << "\n";
                    partition = std::move(exact_small.partition);
                    if (user_requested_target && !opt.strict_exact &&
                            (int)partition.size() <= opt.target_components) {
                        return finish_partition(partition, "target_hit_small_savings_exact");
                    }
                }
            }
        }
        publish_signal_output(inst, partition, edge_offset);
        bool auto_target_loop = !opt.has_target && opt.auto_target_loop;
        unordered_set<int> reference_harvest_targets;
        for (int round = 0; ; round++) {
            if (auto_target_loop) {
                if ((int)partition.size() <= 1) break;
                opt.has_target = true;
                opt.target_components = (int)partition.size() - 1;
                cerr << "AUTO_TARGET_ROUND round=" << round
                     << " current=" << partition.size()
                     << " target=" << opt.target_components << "\n";
            }
            if (!opt.has_target) break;
            int before_round = (int)partition.size();
            if (near_deadline()) return finish_partition(partition, "deadline_after_initial_pack");
            if ((int)partition.size() > opt.target_components) {
                int target_savings = inst.leaf_count - opt.target_components;
                if (target_savings > 0 && target_savings <= opt.small_savings_max_target) {
                    double small_seconds = phase_seconds(opt.small_savings_exact_seconds);
                    if (small_seconds > 0.0) {
                        SmallSavingsExactResult exact_small = small_savings_exhaustive_pack(
                            inst, opt.target_components, edge_offset, small_seconds,
                            opt.small_savings_pair_node_limit,
                            opt.small_savings_special_node_limit,
                            opt.small_savings_max_target,
                            opt.small_savings_component_limit);
                        if (!exact_small.partition.empty() &&
                                validate_partition(inst, exact_small.partition, edge_offset) &&
                                (int)exact_small.partition.size() < (int)partition.size()) {
                            cerr << "USE_SMALL_SAVINGS_EXACT_EARLY components="
                                 << exact_small.partition.size()
                                 << " target=" << opt.target_components << "\n";
                            partition = std::move(exact_small.partition);
                            publish_signal_output(inst, partition, edge_offset);
                            if (user_requested_target && !opt.strict_exact &&
                                    (int)partition.size() <= opt.target_components) {
                                return finish_partition(partition,
                                    "target_hit_small_savings_exact_early");
                            }
                        } else if (exact_small.proved_infeasible &&
                                (int)partition.size() == opt.target_components + 1) {
                            if (!validate_partition(inst, partition, edge_offset)) {
                                throw runtime_error("small savings early certified invalid partition");
                            }
                            publish_signal_output(inst, partition, edge_offset);
                            write_outputs(inst, partition, opt.output_components, opt.output_newick);
                            cerr << "SUMMARY components=" << partition.size()
                                 << " certified=1"
                                 << " method=small_savings_exact_lb_early"
                                 << " lower_target=" << opt.target_components
                                 << " seconds=" << elapsed(start) << "\n";
                            return 0;
                        }
                    }
                }
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_small_savings_exact");
            if ((int)partition.size() > opt.target_components &&
                    opt.reference_cut_seconds > 0.0 &&
                    inst.tree_count <= 8 && inst.leaf_count >= 120 &&
                    reference_harvest_targets.insert(opt.target_components).second) {
                double harvest_seconds = phase_seconds(min(30.0, opt.reference_cut_seconds));
                if (harvest_seconds > 0.25) {
                    int before_pool = (int)pool.comps.size();
                    vector<vector<int>> harvested_direct =
                        reference_cut_cadical_target(
                            inst, opt.target_components, harvest_seconds, nullptr,
                            &pool, 8, 50000);
                    int added = (int)pool.comps.size() - before_pool;
                    if (added > 0) {
                        cols = encode_columns(inst, pool, total_bits);
                    }
                    cerr << "REFERENCE_CUT_HARVEST_ROUND target="
                         << opt.target_components
                         << " added=" << added
                         << " pool=" << pool.comps.size()
                         << " direct=" << harvested_direct.size()
                         << " seconds=" << harvest_seconds << "\n";
                    if (!harvested_direct.empty() &&
                            validate_partition(inst, harvested_direct, edge_offset) &&
                            (int)harvested_direct.size() < (int)partition.size()) {
                        cerr << "USE_REFERENCE_CUT_HARVEST_DIRECT components="
                             << harvested_direct.size() << "\n";
                        partition = std::move(harvested_direct);
                        publish_signal_output(inst, partition, edge_offset);
                    }
                }
            }
            opt.global_branch_seconds = phase_seconds(opt.global_branch_seconds);
            if ((int)partition.size() > opt.target_components && opt.global_branch_seconds > 0.0) {
                vector<int> branch_selected = global_branch_pack(cols, inst.leaf_count, opt);
                if (!branch_selected.empty()) {
                    vector<vector<int>> packed = selected_partition(cols, branch_selected, inst.leaf_count);
                    if (validate_partition(inst, packed, edge_offset) &&
                            (int)packed.size() < (int)partition.size()) {
                        cerr << "USE_GLOBAL_BRANCH components=" << packed.size() << "\n";
                        partition = std::move(packed);
                        publish_signal_output(inst, partition, edge_offset);
                    }
                }
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_global_branch");
            opt.high_divergence_pack_seconds = phase_seconds(opt.high_divergence_pack_seconds);
            if ((int)partition.size() > opt.target_components &&
                    opt.high_divergence_pack_seconds > 0.0) {
                vector<vector<int>> hd_packed = high_divergence_pack(
                    inst, cols, edge_offset, opt);
                if (!opt.dump_high_divergence_lp.empty()) return 0;
                if (!hd_packed.empty() && validate_partition(inst, hd_packed, edge_offset) &&
                        (int)hd_packed.size() < (int)partition.size()) {
                    cerr << "USE_HIGH_DIVERGENCE_PACK components="
                         << hd_packed.size() << "\n";
                    partition = std::move(hd_packed);
                    publish_signal_output(inst, partition, edge_offset);
                }
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_high_divergence_pack");
            opt.global_highs_seconds = phase_seconds(opt.global_highs_seconds);
            if ((int)partition.size() > opt.target_components &&
                    opt.global_highs_seconds > 0.0 && opt.total_seconds > 0.0) {
                cerr << "SKIP_GLOBAL_HIGHS_TOTAL_GUARD seconds=" << opt.global_highs_seconds << "\n";
                opt.global_highs_seconds = 0.0;
            }
            if ((int)partition.size() > opt.target_components && opt.global_highs_seconds > 0.0) {
                vector<vector<int>> packed = global_highs_repack(inst, opt, cols);
                if (!packed.empty() && validate_partition(inst, packed, edge_offset) &&
                        (int)packed.size() < (int)partition.size()) {
                    cerr << "USE_GLOBAL_HIGHS components=" << packed.size() << "\n";
                    partition = std::move(packed);
                    publish_signal_output(inst, partition, edge_offset);
                }
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_global_highs");
            opt.pair_augment_seconds = phase_seconds(opt.pair_augment_seconds);
            if ((int)partition.size() > opt.target_components && opt.pair_augment_seconds > 0.0) {
                partition = pair_augment_repair(inst, std::move(partition), opt, edge_offset, cols);
                publish_signal_output(inst, partition, edge_offset);
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_pair_augment");
            opt.incumbent_repack_seconds = phase_seconds(opt.incumbent_repack_seconds);
            if ((int)partition.size() > opt.target_components &&
                    opt.incumbent_repack_seconds > 0.0) {
                partition = incumbent_guided_repack(inst, std::move(partition), opt, edge_offset, cols);
                publish_signal_output(inst, partition, edge_offset);
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_incumbent_repack");
            opt.exact_exchange_seconds = phase_seconds(opt.exact_exchange_seconds);
            if ((int)partition.size() > opt.target_components &&
                    opt.exact_exchange_seconds > 0.0) {
                partition = exact_region_exchange_repair(
                    inst, std::move(partition), opt, edge_offset, cols);
                publish_signal_output(inst, partition, edge_offset);
            }
            if (near_deadline()) return finish_partition(partition, "deadline_after_exact_exchange");
            opt.local_seconds = phase_seconds(opt.local_seconds);
            partition = local_repack_repair(inst, std::move(partition), opt, edge_offset, cols);
            partition = singleton_merge_repair(inst, std::move(partition), opt.target_components, edge_offset);
            {
                int pair_target = opt.has_target ? opt.target_components :
                    max(1, (int)partition.size() - 1);
                vector<vector<int>> pair_repaired = pair_packing_repair(
                    inst, partition, pair_target, edge_offset);
                if (!pair_repaired.empty()) {
                    cerr << "USE_PAIR_PACK_REPAIR components=" << pair_repaired.size() << "\n";
                    partition = std::move(pair_repaired);
                }
                vector<vector<int>> small_repaired = small_savings_pack_repair(
                    inst, cols, pair_target, edge_offset);
                if (!small_repaired.empty() &&
                        (int)small_repaired.size() < (int)partition.size()) {
                    cerr << "USE_SMALL_SAVINGS_PACK components=" << small_repaired.size() << "\n";
                    partition = std::move(small_repaired);
                    if (user_requested_target && !opt.strict_exact &&
                            (int)partition.size() <= opt.target_components) {
                        return finish_partition(partition, "target_hit_small_savings_pack_loop");
                    }
                }
                SmallSavingsExactResult exact_small = small_savings_exhaustive_pack(
                    inst, pair_target, edge_offset, opt.small_savings_exact_seconds,
                    opt.small_savings_pair_node_limit,
                    opt.small_savings_special_node_limit,
                    opt.small_savings_max_target,
                    opt.small_savings_component_limit);
                if (!exact_small.partition.empty() &&
                        (int)exact_small.partition.size() < (int)partition.size()) {
                    cerr << "USE_SMALL_SAVINGS_EXACT components="
                         << exact_small.partition.size() << "\n";
                    partition = std::move(exact_small.partition);
                    publish_signal_output(inst, partition, edge_offset);
                    if (user_requested_target && !opt.strict_exact &&
                            (int)partition.size() <= opt.target_components) {
                        return finish_partition(partition, "target_hit_small_savings_exact_loop");
                    }
                } else if (exact_small.proved_infeasible &&
                        (int)partition.size() == pair_target + 1) {
                    if (!validate_partition(inst, partition, edge_offset)) {
                        throw runtime_error("small savings certified invalid partition");
                    }
                    publish_signal_output(inst, partition, edge_offset);
                    write_outputs(inst, partition, opt.output_components, opt.output_newick);
                    cerr << "SUMMARY components=" << partition.size()
                         << " certified=1"
                         << " method=small_savings_exact_lb"
                         << " lower_target=" << pair_target
                         << " seconds=" << elapsed(start) << "\n";
                    return 0;
                }
            }
            publish_signal_output(inst, partition, edge_offset);
            if (near_deadline()) return finish_partition(partition, "deadline_after_local");
            if (!auto_target_loop || (int)partition.size() >= before_round) break;
            no_target_successful_rounds++;
            cerr << "AUTO_TARGET_PROGRESS successful_rounds="
                 << no_target_successful_rounds
                 << " from=" << before_round
                 << " to=" << partition.size() << "\n";
            opt.has_target = false;
        }
        if (auto_target_loop) opt.has_target = false;
        opt.direct_pair_seconds = phase_seconds(opt.direct_pair_seconds);
        if (opt.has_target && (int)partition.size() > opt.target_components &&
                inst.leaf_count <= opt.direct_pair_leaf_limit && opt.direct_pair_seconds > 0.0) {
            vector<vector<int>> direct = direct_pair_cadical_target(
                inst, opt.target_components, opt.direct_pair_seconds);
#ifndef USE_CADICAL
            if (direct.empty()) {
                direct = direct_pair_highs_target(
                    inst, opt.target_components, opt.direct_pair_seconds, edge_offset);
            }
#endif
            if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                    (int)direct.size() < (int)partition.size()) {
                cerr << "USE_DIRECT_PAIR components=" << direct.size() << "\n";
                partition = std::move(direct);
                publish_signal_output(inst, partition, edge_offset);
            } else {
                cerr << "DIRECT_PAIR_NO_IMPROVE current=" << partition.size() << "\n";
            }
        }
        if (opt.strict_exact && !opt.has_target && opt.reference_cut_seconds > 0.0) {
            bool restart_reference_cut = true;
            while (restart_reference_cut && (int)partition.size() > 1) {
                restart_reference_cut = false;
                while ((int)partition.size() > 1) {
                    int target = (int)partition.size() - 1;
                    double seconds = phase_seconds(opt.reference_cut_seconds);
                    if (seconds <= 0.0) break;
                    bool proved_unsat = false;
                    no_target_reference_cut_attempted_target = target;
                    vector<vector<int>> cut_partition = reference_cut_cadical_target(
                        inst, target, seconds, &proved_unsat);
                    if (proved_unsat) {
                        if (!validate_partition(inst, partition, edge_offset)) {
                            throw runtime_error("internal invalid partition");
                        }
                        publish_signal_output(inst, partition, edge_offset);
                        write_outputs(inst, partition, opt.output_components, opt.output_newick);
                        cerr << "SUMMARY components=" << partition.size()
                             << " certified=1"
                             << " method=reference_cut_descent"
                             << " lower_target=" << target
                             << " seconds=" << elapsed(start) << "\n";
                        return 0;
                    }
                    if (!cut_partition.empty() && validate_partition(inst, cut_partition, edge_offset) &&
                            (int)cut_partition.size() < (int)partition.size()) {
                        cerr << "STRICT_REFERENCE_CUT_IMPROVE from=" << partition.size()
                             << " to=" << cut_partition.size()
                             << " target=" << target << "\n";
                        partition = std::move(cut_partition);
                        no_target_reference_cut_attempted_target = -1;
                        continue;
                    }
                    cerr << "STRICT_REFERENCE_CUT_NO_CERTIFICATE current=" << partition.size()
                         << " target=" << target << "\n";
                    break;
                }
                if (prioritize_recursive_before_exchange) {
                    cerr << "SKIP_EXCHANGE_BEFORE_RECURSIVE reason=common_cluster_dense_manytree\n";
                } else if (opt.allow_stable_exchange && !near_deadline() &&
                        opt.exchange_max_component_size >= 8 &&
                        opt.exchange_max_region_components >= 6) {
                    double seconds = phase_seconds(opt.exchange_seconds);
                    if (seconds > 0.0) {
                        Options stable_opt = opt;
                        stable_opt.exchange_seconds = seconds;
                        ExchangeResult ex = exchange_certify_partition(
                            inst, stable_opt, edge_offset, partition);
                        cerr << "EXCHANGE_REPAIR_RESULT found=" << (ex.found ? 1 : 0)
                             << " exhausted=" << (ex.exhausted ? 1 : 0)
                             << " candidates=" << ex.candidates
                             << " regions=" << ex.regions
                             << " nodes=" << ex.dfs_nodes
                             << " reason=" << ex.reason << "\n";
                        if (ex.found) {
                            vector<vector<int>> improved = apply_exchange_result(partition, ex);
                            if ((int)improved.size() < (int)partition.size()) {
                                cerr << "USE_EXCHANGE_REPAIR from=" << partition.size()
                                     << " to=" << improved.size()
                                     << " region_old=" << ex.region.size()
                                     << " region_new=" << ex.replacement.size() << "\n";
                                partition = std::move(improved);
                                partition = singleton_merge_repair(
                                    inst, std::move(partition), 1, edge_offset);
                                {
                                    int pair_target = max(1, (int)partition.size() - 1);
                                    vector<vector<int>> pair_repaired = pair_packing_repair(
                                        inst, partition, pair_target, edge_offset);
                                    if (!pair_repaired.empty()) {
                                        cerr << "USE_PAIR_PACK_REPAIR components="
                                             << pair_repaired.size() << "\n";
                                        partition = std::move(pair_repaired);
                                    }
                                }
                                publish_signal_output(inst, partition, edge_offset);
                                no_target_reference_cut_attempted_target = -1;
                                restart_reference_cut = true;
                            }
                        } else if (ex.exhausted) {
                            int savings = inst.leaf_count - (int)partition.size();
                            int min_savings = stable_min_savings();
                            if (savings < min_savings) {
                                cerr << "STABLE_EXCHANGE_SKIPPED reason=low_savings"
                                     << " savings=" << savings
                                     << " min_savings=" << min_savings << "\n";
                                continue;
                            }
                            if (!validate_partition(inst, partition, edge_offset)) {
                                throw runtime_error("internal invalid partition");
                            }
                            cerr << "STABLE_EXCHANGE_NOT_CERTIFICATE components="
                                 << partition.size()
                                 << " lower_target=" << ((int)partition.size() - 1)
                                 << " reason=local_neighborhood_only\n";
                            int lower_target = max(0, (int)partition.size() - 1);
                            if (empirical_exact_accept(partition, lower_target,
                                        "exchange_exhausted_no_certificate")) {
                                publish_signal_output(inst, partition, edge_offset);
                                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                                cerr << "SUMMARY components=" << partition.size()
                                     << " target=" << lower_target
                                     << " certified=0 stable=0 empirical=1"
                                     << " method=empirical_after_exchange"
                                     << " seconds=" << elapsed(start) << "\n";
                                return 0;
                            }
                        }
                    }
                }
            }
        }
        if (opt.strict_exact && !opt.has_target &&
                inst.leaf_count <= opt.direct_pair_leaf_limit &&
                opt.direct_pair_seconds > 0.0) {
            while ((int)partition.size() > 1) {
                int target = (int)partition.size() - 1;
                double seconds = phase_seconds(opt.direct_pair_seconds);
                if (seconds <= 0.0) wait_for_timeout("no_time_for_strict_sat_descent");
                bool proved_unsat = false;
                vector<vector<int>> direct = direct_pair_cadical_target(
                    inst, target, seconds, &proved_unsat);
                if (proved_unsat) {
                    if (!validate_partition(inst, partition, edge_offset)) {
                        throw runtime_error("internal invalid partition");
                    }
                    write_outputs(inst, partition, opt.output_components, opt.output_newick);
                    cerr << "SUMMARY components=" << partition.size()
                         << " certified=1"
                         << " method=direct_pair_sat_descent"
                         << " lower_target=" << target
                         << " seconds=" << elapsed(start) << "\n";
                    return 0;
                }
                if (!direct.empty() && validate_partition(inst, direct, edge_offset) &&
                        (int)direct.size() < (int)partition.size()) {
                    cerr << "STRICT_DESCENT_IMPROVE from=" << partition.size()
                         << " to=" << direct.size()
                         << " target=" << target << "\n";
                    partition = std::move(direct);
                    continue;
                }
                cerr << "STRICT_DESCENT_NO_CERTIFICATE current=" << partition.size()
                     << " target=" << target << "\n";
                wait_for_timeout("strict_sat_descent_not_closed");
            }
        }
        if (defer_strict_recursive_exact) {
            vector<vector<int>> exact = strict_recursive_exact(inst, opt, 0);
            if (!exact.empty()) {
                vector<int> exact_edge_offset = recursive_edge_offset(inst);
                if (!validate_partition(inst, exact, exact_edge_offset)) {
                    throw runtime_error("recursive exact partition failed validation");
                }
                publish_signal_output(inst, exact, exact_edge_offset);
                write_outputs(inst, exact, opt.output_components, opt.output_newick);
                cerr << "SUMMARY components=" << exact.size()
                     << " certified=1"
                     << " method=strict_recursive_exact"
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            cerr << "STRICT_RECURSIVE_EXACT_NO_CERTIFICATE continue=1\n";
            if (stable_exchange_accept(partition)) {
                if (!validate_partition(inst, partition, edge_offset)) {
                    throw runtime_error("internal invalid partition");
                }
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                cerr << "STRICT_STABLE_ACCEPT components=" << partition.size()
                     << " lower_target=" << ((int)partition.size() - 1)
                     << " method=stable_exchange_after_recursive\n";
                cerr << "SUMMARY components=" << partition.size()
                     << " target=" << ((int)partition.size() - 1)
                     << " certified=0 stable=1"
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
            int lower_target = max(0, (int)partition.size() - 1);
            if (empirical_exact_accept(partition, lower_target,
                        "strict_recursive_no_certificate")) {
                if (!validate_partition(inst, partition, edge_offset)) {
                    throw runtime_error("internal invalid partition");
                }
                publish_signal_output(inst, partition, edge_offset);
                write_outputs(inst, partition, opt.output_components, opt.output_newick);
                cerr << "SUMMARY components=" << partition.size()
                     << " target=" << lower_target
                     << " certified=0 stable=0 empirical=1"
                     << " method=empirical_after_recursive"
                     << " seconds=" << elapsed(start) << "\n";
                return 0;
            }
        }
        if (!validate_partition(inst, partition, edge_offset)) throw runtime_error("internal invalid partition");
        publish_signal_output(inst, partition, edge_offset);
        certify_exact_or_timeout(partition, "no_certificate_or_target_miss");
        write_outputs(inst, partition, opt.output_components, opt.output_newick);
        cerr << "SUMMARY components=" << partition.size()
             << " target=" << opt.target_components
             << " certified=" << (opt.strict_exact && !stable_exchange_used && !empirical_exact_used ? 1 : 0)
             << " stable=" << (stable_exchange_used ? 1 : 0)
             << " empirical=" << (empirical_exact_used ? 1 : 0)
             << " seconds=" << elapsed(start) << "\n";
        return (!opt.has_target || partition.size() <= (size_t)opt.target_components) ? 0 : 2;
    } catch (const exception &ex) {
        cerr << "error: " << ex.what() << "\n";
        return 2;
    }
}
