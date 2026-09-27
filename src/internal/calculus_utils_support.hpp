#pragma once

#include "calculus_utils.hpp"

namespace LMCAS::calculus_utils_detail {

Result<void> calculus_utils_validate_expr(const std::shared_ptr<SymbolicExpr>& expr,
                                          const std::string& var,
                                          ComputationContext& context,
                                          const std::string& operation);

Result<void> calculus_utils_validate_two_exprs(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation);

Result<void> calculus_utils_validate_expr_target(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& target,
    ComputationContext& context,
    const std::string& operation);

Result<void> calculus_utils_validate_expr_bounds(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation);

bool calculus_utils_is_infinity(const std::shared_ptr<SymbolicExpr>& expr);
bool calculus_utils_expr_equal(const std::shared_ptr<SymbolicExpr>& a,
                               const std::shared_ptr<SymbolicExpr>& b);
std::shared_ptr<SymbolicExpr> calculus_utils_make_abs(
    const std::shared_ptr<SymbolicExpr>& expr);

}
