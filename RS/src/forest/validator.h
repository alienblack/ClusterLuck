#pragma once

#include <string>

#include "forest/forest.h"

struct ValidationResult {
    bool ok = false;
    std::string message;
};

ValidationResult validate_forest(const Forest& forest, const Tree& t1, const Tree& t2);

