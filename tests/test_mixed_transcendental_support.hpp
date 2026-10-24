#include "expr.hpp"
#pragma once

#include "test_common.hpp"
#include "solve_mixed_transcendental.hpp"
#include "solve_strategies.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace LMCAS;

inline std::optional<double> eval_at(const std::shared_ptr<SymbolicExpr>& expr,
                                     const std::string& var, double val) {
    try {
        auto substituted = LMCAS::substitute(expr, var, SymbolicExpr::number(val)).value();
        double result = substituted->to_numeric();
        if (std::isnan(result) || std::isinf(result)) return std::nullopt;
        return result;
    } catch (...) {
        return std::nullopt;
    }
}

inline std::shared_ptr<SymbolicExpr> compute_derivative(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !LMCAS::detail::node(expr)) return nullptr;
    try {
        auto derivative = LMCAS::differentiate(expr, var);
        return derivative ? std::move(derivative.value()) : nullptr;
    } catch (...) {
        return nullptr;
    }
}
