#pragma once

#include <istream>

#include "tree/tree.h"

struct Instance {
    int n = 0;
    bool has_lower_bound_rule = false;
    double lower_bound_factor = 0.0;
    int lower_bound_slack = 0;
    Tree t1;
    Tree t2;
};

Instance read_instance(std::istream& in);
