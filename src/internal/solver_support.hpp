#pragma once

#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::solver_detail {

inline lmmc_real_t finite_midpoint(lmmc_real_t lower, lmmc_real_t upper) {
    if (lower < 0.0 && upper > 0.0) {
        return lower * 0.5 + upper * 0.5;
    }
    return lower + (upper - lower) * 0.5;
}

std::shared_ptr<SymbolicExpr> to_ptr(const SymbolicExpr& expression);
std::shared_ptr<SymbolicExpr> multiply_no_expand(
    const std::shared_ptr<const SymbolicNode>& term,
    const std::vector<std::shared_ptr<const SymbolicNode>>& denominator_factors);
bool is_polynomial_node(const std::shared_ptr<const SymbolicNode>& node);
std::shared_ptr<SymbolicExpr> multiply_factors(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors);
bool collect_denominator_factors(
    const std::shared_ptr<const SymbolicNode>& node,
    std::vector<std::shared_ptr<const SymbolicNode>>& denominator_factors,
    std::vector<std::shared_ptr<SymbolicExpr>>& denominator_constraints);

bool contains_imaginary(const std::shared_ptr<SymbolicExpr>& expression);
bool solver_number_value(const NumberNode& number, double& value);

} // namespace LMCAS::solver_detail
