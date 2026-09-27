#pragma once

#include "test_common.hpp"
#include "solve_mixed_transcendental.hpp"
#include "solve_strategies.hpp"
#include "internal/visitors/differentiation_visitor.hpp"

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
        auto substituted = expr->substitute(var, SymbolicExpr::number(val));
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
        DifferentiationVisitor dv(var);
        LMCAS::detail::node(expr)->accept(dv);
        auto result_node = dv.get_result();
        if (!result_node) return nullptr;
        return LMCAS::detail::make_expression_ptr(result_node);
    } catch (...) {
        return nullptr;
    }
}
