#include "whidden/whidden_bridge.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

static uint32_t g_portable_rand_state = 1;

void portable_srand(unsigned seed) {
    constexpr uint32_t modulus = 2147483647u;
    g_portable_rand_state = seed % modulus;
    if (g_portable_rand_state == 0) g_portable_rand_state = 1;
}

int portable_rand() {
    constexpr uint64_t multiplier = 16807u;
    constexpr uint64_t modulus = 2147483647u;
    g_portable_rand_state =
        static_cast<uint32_t>((multiplier * g_portable_rand_state) % modulus);
    return static_cast<int>(g_portable_rand_state);
}

}  // namespace

#define srand portable_srand
#define rand portable_rand

#define Forest WhiddenForest
#define Node WhiddenNode
#define ClusterForest WhiddenClusterForest
#define ClusterInstance WhiddenClusterInstance
#define LCA WhiddenLCA
#define SiblingPair WhiddenSiblingPair
#define UndoMachine WhiddenUndoMachine

class WhiddenForest;
class WhiddenNode;
class WhiddenClusterForest;

#ifndef H25_GREEDY_SELECTOR_WEIGHT
#define H25_GREEDY_SELECTOR_WEIGHT 2
#endif

static int H25_GLOBAL_GREEDY_DISTANCE = -1;
static int h25_local_greedy_distance(
    WhiddenForest* first, WhiddenForest* second,
    WhiddenForest** out_first = nullptr,
    WhiddenForest** out_second = nullptr);
static int h25_greedy_subtree_score(
    WhiddenNode* candidate, WhiddenForest* first, WhiddenForest* second);

#include "../../third_party/whidden/rspr.h"

#undef Forest
#undef Node
#undef ClusterForest
#undef ClusterInstance
#undef LCA
#undef SiblingPair
#undef UndoMachine
#undef srand
#undef rand

#include "forest/forest.h"
#include "heuristic/greedy_merge.h"
#include "io/newick.h"

static std::string h25_restricted_local_newick(
    WhiddenNode* root,
    const std::unordered_map<WhiddenNode*, int>& local_label) {
    if (root == nullptr) return "";
    if (root->get_children().empty()) {
        auto found = local_label.find(root);
        return found == local_label.end() ? "" : std::to_string(found->second);
    }

    std::vector<std::string> children;
    children.reserve(root->get_children().size());
    for (WhiddenNode* child : root->get_children()) {
        std::string restricted =
            h25_restricted_local_newick(child, local_label);
        if (!restricted.empty()) children.push_back(std::move(restricted));
    }
    if (children.empty()) return "";
    if (children.size() == 1) return children.front();

    std::string result = "(";
    for (size_t i = 0; i < children.size(); ++i) {
        if (i > 0) result.push_back(',');
        result += children[i];
    }
    result.push_back(')');
    return result;
}

static int h25_local_greedy_distance(
    WhiddenForest* first, WhiddenForest* second,
    WhiddenForest** out_first,
    WhiddenForest** out_second) {
    (void)first;
    (void)second;
    if (out_first != nullptr) *out_first = nullptr;
    if (out_second != nullptr) *out_second = nullptr;
    return -1;
}

static int h25_greedy_subtree_score(
    WhiddenNode* candidate, WhiddenForest* first, WhiddenForest* second) {
    if (candidate == nullptr || first == nullptr || second == nullptr) {
        return -1;
    }

    std::vector<WhiddenNode*> first_leaves = candidate->find_leaves();
    if (first_leaves.size() <= 1) return 0;

    std::vector<WhiddenNode*> all_first_leaves;
    std::vector<WhiddenNode*> all_second_leaves;
    for (int cid = 0; cid < first->num_components(); ++cid) {
        std::vector<WhiddenNode*> leaves =
            first->get_component(cid)->find_leaves();
        all_first_leaves.insert(
            all_first_leaves.end(), leaves.begin(), leaves.end());
    }
    for (int cid = 0; cid < second->num_components(); ++cid) {
        std::vector<WhiddenNode*> leaves =
            second->get_component(cid)->find_leaves();
        all_second_leaves.insert(
            all_second_leaves.end(), leaves.begin(), leaves.end());
    }

    std::unordered_map<std::string, int> first_name_count;
    std::unordered_map<std::string, std::vector<WhiddenNode*>> second_by_name;
    std::unordered_map<WhiddenNode*, WhiddenNode*> second_by_original_first;
    first_name_count.reserve(all_first_leaves.size() * 2 + 1);
    second_by_name.reserve(all_second_leaves.size() * 2 + 1);
    second_by_original_first.reserve(all_second_leaves.size() * 2 + 1);
    for (WhiddenNode* leaf : all_first_leaves) {
        first_name_count[leaf->get_name()]++;
    }
    for (WhiddenNode* leaf : all_second_leaves) {
        second_by_name[leaf->get_name()].push_back(leaf);
        WhiddenNode* original_first = leaf->get_twin();
        if (original_first != nullptr &&
            !second_by_original_first.count(original_first)) {
            second_by_original_first[original_first] = leaf;
        }
    }

    std::unordered_map<WhiddenNode*, int> first_label;
    std::unordered_map<WhiddenNode*, int> second_label;
    first_label.reserve(first_leaves.size() * 2 + 1);
    second_label.reserve(first_leaves.size() * 2 + 1);
    for (int i = 0; i < static_cast<int>(first_leaves.size()); ++i) {
        WhiddenNode* matched_second = nullptr;
        const std::string name = first_leaves[i]->get_name();
        auto same_name = second_by_name.find(name);
        if (first_name_count[name] == 1 &&
            same_name != second_by_name.end() &&
            same_name->second.size() == 1) {
            matched_second = same_name->second.front();
        } else {
            WhiddenNode* original_second = first_leaves[i]->get_twin();
            WhiddenNode* original_first =
                original_second == nullptr ? nullptr : original_second->get_twin();
            auto found = second_by_original_first.find(original_first);
            if (original_first != nullptr &&
                found != second_by_original_first.end()) {
                matched_second = found->second;
            }
        }
        if (matched_second == nullptr || second_label.count(matched_second)) {
            return -1;
        }
        first_label[first_leaves[i]] = i + 1;
        second_label[matched_second] = i + 1;
    }

    std::string first_newick =
        h25_restricted_local_newick(candidate, first_label);
    std::vector<std::string> second_parts;
    second_parts.reserve(second->num_components());
    for (int cid = 0; cid < second->num_components(); ++cid) {
        std::string part = h25_restricted_local_newick(
            second->get_component(cid), second_label);
        if (!part.empty()) second_parts.push_back(std::move(part));
    }
    if (first_newick.empty() || second_parts.empty()) return -1;

    std::string second_newick;
    if (second_parts.size() == 1) {
        second_newick = std::move(second_parts.front());
    } else {
        second_newick.push_back('(');
        for (size_t i = 0; i < second_parts.size(); ++i) {
            if (i > 0) second_newick.push_back(',');
            second_newick += second_parts[i];
        }
        second_newick.push_back(')');
    }

    const int leaf_count = static_cast<int>(first_leaves.size());
    Tree t1 = parse_newick_tree_h45_order(first_newick, leaf_count);
    Tree t2 = parse_newick_tree_h45_order(second_newick, leaf_count);
    Forest singleton = make_singleton_forest(leaf_count);
    if (!rebuild_forest_data(singleton, t1, t2)) return -1;
    Forest greedy =
        run_h45_greedy_merge(std::move(singleton), t1, t2, 0, false, false);
    return std::max(0, greedy.active_count - 1);
}

namespace {

void configure_whidden_bridge(int cluster_tune, int split_threshold,
                              int random_seed) {
    CUT_ALL_B = true;
    CUT_ONE_B = true;
    REVERSE_CUT_ONE_B = true;
    REVERSE_CUT_ONE_B_3 = true;
    CUT_TWO_B = true;
    CUT_AC_SEPARATE_COMPONENTS = true;
    EDGE_PROTECTION = true;
    EDGE_PROTECTION_TWO_B = true;
    NEAR_PREORDER_SIBLING_PAIRS = true;
    LEAF_REDUCTION = true;
    LEAF_REDUCTION2 = true;
    PREFER_NONBRANCHING = true;
    APPROX_CUT_ONE_B = true;
    APPROX_CUT_TWO_B = true;
    APPROX_REVERSE_CUT_ONE_B = true;
    DEEPEST_PROTECTED_ORDER = true;
    DEEPEST_ORDER = true;
    PREORDER_SIBLING_PAIRS = true;
    BB = true;
    PREFER_RHO = true;
    CLUSTER_TUNE = cluster_tune;
    SPLIT_APPROX = true;
    SPLIT_APPROX_THRESHOLD = split_threshold;
    portable_srand(static_cast<unsigned>(random_seed));
}

bool parse_leaf_label(const std::string& text, int leaf_count, int& value) {
    if (text.empty()) return false;
    long long parsed = 0;
    for (unsigned char ch : text) {
        if (!std::isdigit(ch)) return false;
        parsed = parsed * 10 + (ch - '0');
        if (parsed > leaf_count) return false;
    }
    if (parsed < 1) return false;
    value = static_cast<int>(parsed);
    return true;
}

bool component_leaves(WhiddenNode* component, int leaf_count,
                      std::vector<int>& leaves) {
    leaves.clear();
    if (component == nullptr) return false;
    std::vector<WhiddenNode*> found = component->find_leaves();
    leaves.reserve(found.size());
    for (WhiddenNode* leaf : found) {
        int value = 0;
        if (leaf == nullptr ||
            !parse_leaf_label(leaf->get_name(), leaf_count, value)) {
            return false;
        }
        leaves.push_back(value);
    }
    std::sort(leaves.begin(), leaves.end());
    return std::adjacent_find(leaves.begin(), leaves.end()) == leaves.end();
}

bool convert_whidden_forest(WhiddenForest* source, int leaf_count,
                            std::map<int, std::string>* reverse_label_map,
                            std::vector<std::vector<int>>& components) {
    if (source == nullptr) return false;
    source->numbers_to_labels(reverse_label_map);
    expand_contracted_nodes(source);

    components.clear();
    std::vector<unsigned char> seen(leaf_count + 1, 0);
    for (int i = 0; i < source->num_components(); ++i) {
        std::vector<int> leaves;
        if (!component_leaves(source->get_component(i), leaf_count, leaves)) {
            return false;
        }
        if (leaves.empty()) continue;
        for (int leaf : leaves) {
            if (leaf < 1 || leaf > leaf_count || seen[leaf]) return false;
            seen[leaf] = 1;
        }
        components.push_back(std::move(leaves));
    }
    for (int leaf = 1; leaf <= leaf_count; ++leaf) {
        if (!seen[leaf]) return false;
    }
    return true;
}

}  // namespace

bool whidden_seed_components(const std::string& first_newick,
                             const std::string& second_newick,
                             int leaf_count,
                             int greedy_upper_bound,
                             int cluster_tune,
                             int split_threshold,
                             int random_seed,
                             std::vector<std::vector<int>>& components) {
    configure_whidden_bridge(cluster_tune, split_threshold, random_seed);
    H25_GLOBAL_GREEDY_DISTANCE = greedy_upper_bound;

    WhiddenNode* first = build_tree(first_newick);
    WhiddenNode* second = build_tree(second_newick);
    if (first == nullptr || second == nullptr) {
        H25_GLOBAL_GREEDY_DISTANCE = -1;
        return false;
    }

    std::map<std::string, int> label_map;
    std::map<int, std::string> reverse_label_map;
    first->labels_to_numbers(&label_map, &reverse_label_map);
    second->labels_to_numbers(&label_map, &reverse_label_map);

    WhiddenForest* out_first = nullptr;
    WhiddenForest* out_second = nullptr;
    int whidden_k = rSPR_branch_and_bound_simple_clustering(
        first, second, false, &label_map, &reverse_label_map,
        -1, -1, &out_first, &out_second);
    (void)whidden_k;
    H25_GLOBAL_GREEDY_DISTANCE = -1;

    first->delete_tree();
    second->delete_tree();
    if (out_first == nullptr) {
        delete out_second;
        return false;
    }

    bool ok = convert_whidden_forest(out_first, leaf_count,
                                     &reverse_label_map, components);
    delete out_first;
    delete out_second;
    return ok;
}
