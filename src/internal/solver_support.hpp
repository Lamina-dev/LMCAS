#pragma once

#include "symbolic.hpp"
#include "symbolic_ast.hpp"

namespace lamina::solver_detail {

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

} // namespace lamina::solver_detail
