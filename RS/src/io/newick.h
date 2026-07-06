#pragma once

#include <string>

#include "tree/tree.h"

Tree parse_newick_tree(const std::string& text, int leaf_count);
Tree parse_newick_tree_h45_order(const std::string& text, int leaf_count);
