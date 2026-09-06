#pragma once

#include "symbolic_ode_engine.hpp"

namespace LMCAS {

double try_eval_double(const std::shared_ptr<SymbolicExpr>& expression);

Result<void> validate_ode_expr_var_pair(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation);

Result<void> validate_ode_pair_var_pair(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation);

ODESolutionResult wrap_ode_solution(
    ODESolution solution,
    ODEType expected_method,
    const std::string& operation);

Result<void> validate_ode_variables(
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation);

Result<void> validate_numeric_ode_coefficients(
    const std::vector<double>& coefficients,
    std::size_t minimum_size,
    std::size_t maximum_size,
    const std::string& operation);

Result<void> validate_ode_three_expr_one_var(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::shared_ptr<SymbolicExpr>& third,
    const std::string& x,
    ComputationContext& context,
    const std::string& operation);

Result<void> validate_ode_two_expr_point(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    ComputationContext& context,
    const std::string& operation);

} // namespace LMCAS
