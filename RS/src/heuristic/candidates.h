#pragma once

#include <vector>

#include "forest/forest.h"

struct MergeCandidate {
    CompId a = -1;
    CompId b = -1;
    int score = 0;
};

enum class CandidateMode {
    SiblingAllPairs,
    WindowScored,
    WindowTight,
    WindowWide,
    DistanceScored,
    LocalDistance,
};

const char* candidate_mode_name(CandidateMode mode);

std::vector<MergeCandidate> generate_candidates(const Forest& forest, const Tree& t1,
                                                const Tree& t2, int all_pair_cap,
                                                CandidateMode mode, int jitter_seed = 0,
                                                int jitter_span = 0);
