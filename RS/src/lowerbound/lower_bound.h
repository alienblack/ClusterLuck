#pragma once

#include "io/pace_io.h"

struct LowerBoundCert {
    bool enabled = false;
    int lower_bound = 1;
    int target_components = 0;
};

LowerBoundCert compute_lower_bound_certificate(const Instance& instance);
bool lower_bound_target_met(int components, const LowerBoundCert& cert);
