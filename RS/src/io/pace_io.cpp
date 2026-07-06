#include "io/pace_io.h"

#include <cctype>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/newick.h"

namespace {

std::string trim(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

int max_leaf_label(const std::string& a, const std::string& b) {
    int best = 0;
    for (const std::string* text : {&a, &b}) {
        for (size_t i = 0; i < text->size();) {
            if (!std::isdigit(static_cast<unsigned char>((*text)[i]))) {
                ++i;
                continue;
            }
            int value = 0;
            while (i < text->size() &&
                   std::isdigit(static_cast<unsigned char>((*text)[i]))) {
                value = value * 10 + ((*text)[i] - '0');
                ++i;
            }
            best = std::max(best, value);
        }
    }
    return best;
}

}  // namespace

Instance read_instance(std::istream& in) {
    std::vector<std::string> trees;
    std::string line;
    int n = 0;
    bool has_lower_bound_rule = false;
    double lower_bound_factor = 0.0;
    int lower_bound_slack = 0;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.rfind("#p", 0) == 0) {
                std::stringstream ss(line);
                std::string hash;
                int tree_count = 0;
                ss >> hash >> tree_count >> n;
                if (hash != "#p" || tree_count != 2 || n <= 0) {
                    throw std::runtime_error("expected header: #p 2 n");
                }
            } else if (line.rfind("#a", 0) == 0) {
                std::stringstream ss(line);
                std::string hash;
                ss >> hash >> lower_bound_factor >> lower_bound_slack;
                if (hash == "#a" && !ss.fail() && lower_bound_factor >= 1.0 &&
                    lower_bound_slack >= 0) {
                    has_lower_bound_rule = true;
                }
            }
            continue;
        }
        trees.push_back(line);
    }
    if (trees.empty()) throw std::runtime_error("empty input");
    if (trees.size() < 2) throw std::runtime_error("expected two trees");

    if (n == 0) n = max_leaf_label(trees[0], trees[1]);
    if (n <= 0) throw std::runtime_error("could not infer leaf count");

    Instance instance;
    instance.n = n;
    instance.has_lower_bound_rule = has_lower_bound_rule;
    instance.lower_bound_factor = lower_bound_factor;
    instance.lower_bound_slack = lower_bound_slack;
    instance.t1 = parse_newick_tree_h45_order(trees[0], n);
    instance.t2 = parse_newick_tree_h45_order(trees[1], n);
    return instance;
}
