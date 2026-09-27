#pragma once

#include "transform_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include <utility>
#include <optional>

namespace LMCAS::transform_detail {

bool te_depends_on(const std::shared_ptr<SymbolicExpr>& expression,
                   const std::string& variable);
Result<void> te_validate_expr_vars(const std::shared_ptr<SymbolicExpr>& expression,
                                   const std::string& input_var,
                                   const std::string& output_var,
                                   ComputationContext& context,
                                   const std::string& operation);
bool te_contains_transform(const std::shared_ptr<const SymbolicNode>& node);
TransformEngineResult te_wrap_transform_result(std::shared_ptr<SymbolicExpr> expression,
                                               const std::string& operation);
std::shared_ptr<SymbolicExpr> te_gt_condition(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs);
bool te_is_zero(const std::shared_ptr<SymbolicExpr>& expression);
std::optional<double> te_number_value(
    const std::shared_ptr<const SymbolicNode>& node);
BigInt te_factorial(const BigInt& n, ComputationContext& context);
std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>> te_split_coeff(
    const std::shared_ptr<SymbolicExpr>& expression, const std::string& variable);
std::shared_ptr<SymbolicExpr> te_linear_coefficient(
    const std::shared_ptr<const SymbolicNode>& argument, const std::string& variable);
bool te_is_trig(const std::shared_ptr<SymbolicExpr>& expression,
                const std::string& variable, bool& is_sin,
                std::shared_ptr<SymbolicExpr>& frequency);

struct GaussianMatch {
    std::shared_ptr<SymbolicExpr> coefficient;
    double rate;
};

std::optional<double> te_match_pure_gaussian(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);
std::optional<GaussianMatch> te_match_scaled_gaussian(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);
std::shared_ptr<SymbolicExpr> laplace_transform_core(
    const std::shared_ptr<SymbolicExpr>& expression, const std::string& input_var,
    const std::string& output_var, ComputationContext& context);

}
