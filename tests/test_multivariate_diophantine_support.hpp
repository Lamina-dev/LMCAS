#pragma once

#include "test_multivariate_support.hpp"
#include <rapidcheck.h>
using namespace LMCAS;

inline bool verify_diophantine_solution(
    const std::vector<MultiPoly>& factors,
    const std::vector<MultiPoly>& solution,
    const MultiPoly& target)
{
    if (solution.size() != factors.size()) { return false; }

    const auto& vars = factors[0].variables();
    MultiPoly sum(Rational(0), vars);
    for (size_t i = 0; i < factors.size(); ++i) {
        sum = sum + solution[i] * factors[i];
    }

    return sum == target;
}
