#pragma once

#include "vector_calculus.hpp"
#include "symbolic_ast.hpp"

namespace LMCAS::vector_calculus_detail {

std::vector<std::string> vector_calculus_coord_vars(std::size_t dimension);

std::shared_ptr<SymbolicExpr> vector_calculus_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context);

VectorCalculusExprResult curve_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& function,
    const VectorField& parametrization,
    const std::string& parameter,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult curve_integral_vector_strict(
    const VectorField& field,
    const VectorField& parametrization,
    const std::string& parameter,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult surface_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& function,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult surface_integral_vector_strict(
    const VectorField& field,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult greens_theorem_strict(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::vector<std::string>& variables,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult greens_theorem_area_strict(
    const VectorField& parametrization,
    const std::string& parameter,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult divergence_theorem_strict(
    const VectorField& field,
    const std::vector<std::string>& variables,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context,
    const std::string& operation);

VectorCalculusExprResult stokes_theorem_strict(
    const VectorField& field,
    const std::vector<std::string>& variables,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>,
                    std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context,
    const std::string& operation);

bool vector_calculus_checked_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    double& value,
    ComputationContext* context = nullptr);

bool vector_calculus_contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node,
    std::size_t depth = 0);

Result<void> vector_calculus_validate_expr_vars(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_field_vars(
    const VectorField& field,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_functions_vars(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_curve_parametrization(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_curve_scalar_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_curve_vector_inputs(
    const VectorField& F,
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_surface_parametrization(
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_surface_scalar_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_surface_vector_inputs(
    const VectorField& F,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation);

Result<void> vector_calculus_validate_bound_pair(
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& bounds,
    ComputationContext& context,
    const std::string& operation,
    const std::string& label);

Result<void> vector_calculus_validate_distinct_vars(
    const std::vector<std::string>& vars,
    std::size_t expected,
    ComputationContext& context,
    const std::string& operation);

bool vector_calculus_expr_zero_after_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& point);

bool vector_calculus_point_has_vars(
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& point,
    const std::vector<std::string>& vars);

VectorCalculusExprResult vector_calculus_wrap_expr(
    std::shared_ptr<SymbolicExpr> expr,
    const std::string& operation);

VectorCalculusFieldResult vector_calculus_wrap_field(
    VectorField field,
    const std::string& operation);

} // namespace LMCAS::vector_calculus_detail

namespace LMCAS {

VectorCalculusExprResult vector_calculus_simplify_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message);

VectorCalculusExprResult vector_calculus_differentiate_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable,
    const std::string& operation);

} // namespace LMCAS
