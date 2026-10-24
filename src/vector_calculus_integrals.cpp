#include "internal/expression_transform.hpp"
#include "internal/vector_calculus_support.hpp"
#include "integration.hpp"
#include "internal/integration_support.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "internal/symbolic_ast.hpp"
#include "residual_verification.hpp"
#include "internal/expression_analysis.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace LMCAS {

using namespace vector_calculus_detail;

std::vector<std::string> vector_calculus_detail::vector_calculus_coord_vars(
    size_t dim)
{
    if (dim == 2) return {"x", "y"};
    return {"x", "y", "z"};
}

static VectorCalculusExprResult vector_calculus_try_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    detail::RewriteBudget& transform_budget)
{
    Integrator integrator;
    auto integrated = integrator.integrate_def_checked(
        *integrand, var, *a, *b, context);
    if (!integrated) {
        return VectorCalculusExprResult::failure(integrated.error());
    }
    SymbolicExpr result = std::move(integrated.value());

    if (contains_unevaluated_integral(LMCAS::detail::node(result))) {
        return std::shared_ptr<SymbolicExpr>{};
    }
    auto res = LMCAS::detail::make_expression_ptr(result);
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(res), transform_budget));
    if (simplified &&
        contains_unevaluated_integral(LMCAS::detail::node(simplified))) {
        return std::shared_ptr<SymbolicExpr>{};
    }
    return simplified ? simplified : res;
}

static VectorCalculusExprResult vector_calculus_numerical_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    detail::RewriteBudget& transform_budget)
{
    double a_val = 0.0;
    double b_val = 0.0;
    if (!vector_calculus_checked_finite_numeric(a, a_val, context) ||
        !vector_calculus_checked_finite_numeric(b, b_val, context)) {
        return std::shared_ptr<SymbolicExpr>{};
    }

    int n = 1000;
    double h = (b_val - a_val) / n;
    if (!std::isfinite(h)) return std::shared_ptr<SymbolicExpr>{};
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        auto step_result =
            context.consume_steps(1, "vector_calculus.numeric_integral");
        if (!step_result) {
            return VectorCalculusExprResult::failure(step_result.error());
        }
        double xi = a_val + i * h;
        if (!std::isfinite(xi)) return std::shared_ptr<SymbolicExpr>{};
        auto xi_expr = SymbolicExpr::number(xi);
        auto fi = detail::make_expression_ptr(detail::substitute_node(detail::node(integrand), var, detail::node(xi_expr), transform_budget));
        double fi_val = 0.0;
        if (!vector_calculus_checked_finite_numeric(fi, fi_val, context)) {
            return std::shared_ptr<SymbolicExpr>{};
        }

        double weight = 1.0;
        if (i == 0 || i == n) {
            weight = 1.0;
        } else if (i % 2 == 1) {
            weight = 4.0;
        } else {
            weight = 2.0;
        }
        sum += weight * fi_val;
        if (!std::isfinite(sum)) return std::shared_ptr<SymbolicExpr>{};
    }
    sum *= h / 3.0;
    if (!std::isfinite(sum)) return std::shared_ptr<SymbolicExpr>{};

    return SymbolicExpr::number(sum);
}

std::shared_ptr<SymbolicExpr>
vector_calculus_detail::vector_calculus_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    detail::RewriteBudget& transform_budget)
{
    auto symbolic_result = vector_calculus_try_definite(integrand, var, a, b, context, transform_budget);
    if (symbolic_result && symbolic_result.value()) {
        return std::move(symbolic_result.value());
    }
    if (!symbolic_result) throw symbolic_result.error();
    auto numerical_result = vector_calculus_numerical_definite(integrand, var, a, b, context, transform_budget);
    if (!numerical_result) throw numerical_result.error();
    return std::move(numerical_result.value());
}

static VectorCalculusExprResult vector_calculus_inconclusive(
    const std::string& operation,
    const std::string& message)
{
    return VectorCalculusExprResult::failure(CasErrc::Inconclusive,
                                             message, operation);
}

VectorCalculusExprResult vector_calculus_detail::vector_calculus_simplify_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message,
    detail::RewriteBudget& transform_budget)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return vector_calculus_inconclusive(operation, message);
    }
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(expr), transform_budget));
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return vector_calculus_inconclusive(operation, message);
    }
    return VectorCalculusExprResult::success(std::move(simplified));
}

VectorCalculusExprResult vector_calculus_detail::vector_calculus_differentiate_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
    try {
        auto derivative = detail::make_expression_ptr(detail::differentiate_node(detail::node(expr), var, transform_budget));
        return vector_calculus_simplify_strict(derivative, operation,
        "expression derivative is outside the supported domain", transform_budget);
    } catch (const CasError& error) {
        if (error.code != CasErrc::UnsupportedExpression) {
            return VectorCalculusExprResult::failure(error);
        }
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
}

static VectorCalculusExprResult vector_calculus_substitute_coords_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& coord_vars,
    const VectorField& values,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    if (!expr || !LMCAS::detail::node(expr) || coord_vars.size() != values.size()) {
        return vector_calculus_inconclusive(
            operation, "coordinate substitution is outside the supported domain");
    }
    auto substituted = expr;
    for (size_t i = 0; i < coord_vars.size(); ++i) {
        if (!values[i] || !LMCAS::detail::node(values[i])) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
        substituted = detail::make_expression_ptr(detail::substitute_node(detail::node(substituted), coord_vars[i], detail::node(values[i]), transform_budget));
        if (!substituted || !LMCAS::detail::node(substituted)) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
    }
    return vector_calculus_simplify_strict(substituted, operation,
    "coordinate substitution is outside the supported domain", transform_budget);
}

static VectorCalculusExprResult vector_calculus_definite_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    auto result = vector_calculus_try_definite(integrand, var, a, b, context, transform_budget);
    if (!result) return result;
    if (!result.value() || !LMCAS::detail::node(result.value())) {
        return vector_calculus_inconclusive(
            operation,
            "integral could not be evaluated exactly in the supported domain");
    }
    return result;
}

static VectorCalculusExprResult vector_calculus_multiple_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::vector<IntegrationStep>& steps,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    if (!integrand || !LMCAS::detail::node(integrand)) {
        return vector_calculus_inconclusive(
            operation, "integrand construction is outside the supported domain");
    }
    Integrator integrator;
    auto integrated = integrate_multiple_checked(*integrand, steps, integrator, context);
    if (!integrated) return VectorCalculusExprResult::failure(integrated.error());
    auto result = LMCAS::detail::make_expression_ptr(integrated.value());
    return vector_calculus_simplify_strict(result, operation,
    "integral result simplification is outside the supported domain", transform_budget);
}

static Result<VectorField> vector_calculus_cross_product_partials_strict(
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    VectorField r_u;
    VectorField r_v;
    r_u.reserve(3);
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto du = vector_calculus_differentiate_strict(parametrization[i], u, operation, transform_budget);
        if (!du) return Result<VectorField>::failure(du.error());
        auto dv = vector_calculus_differentiate_strict(parametrization[i], v, operation, transform_budget);
        if (!dv) return Result<VectorField>::failure(dv.error());
        r_u.push_back(std::move(du.value()));
        r_v.push_back(std::move(dv.value()));
    }

    VectorField cross;
    cross.reserve(3);
    for (std::size_t index = 0; index < 3; ++index) {
        auto component = vector_calculus_cross_component(
            r_u, r_v, index);
        auto checked = vector_calculus_simplify_strict(component, operation,
        "surface normal construction is outside the supported domain", transform_budget);
        if (!checked) {
            return Result<VectorField>::failure(checked.error());
        }
        cross.push_back(std::move(checked.value()));
    }
    return Result<VectorField>::success(std::move(cross));
}

VectorCalculusExprResult
vector_calculus_detail::curve_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    auto f_composed = vector_calculus_substitute_coords_strict(f, coord_vars, parametrization, operation, transform_budget);
    if (!f_composed) {
        return f_composed;
    }

    std::shared_ptr<SymbolicExpr> speed_sq;
    for (const auto& component : parametrization) {
        auto derivative = vector_calculus_differentiate_strict(component, t, operation, transform_budget);
        if (!derivative) {
            return derivative;
        }
        auto sq = SymbolicExpr::power(derivative.value(), SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(sq, operation, "curve speed construction is outside the supported domain", transform_budget);
        if (!sq_checked) {
            return sq_checked;
        }
        speed_sq = speed_sq ? SymbolicExpr::add(speed_sq, sq_checked.value())
                            : sq_checked.value();
        auto speed_sq_checked = vector_calculus_simplify_strict(speed_sq, operation, "curve speed construction is outside the supported domain", transform_budget);
        if (!speed_sq_checked) {
            return speed_sq_checked;
        }
        speed_sq = std::move(speed_sq_checked.value());
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    auto speed_checked = vector_calculus_simplify_strict(speed, operation, "curve speed construction is outside the supported domain", transform_budget);
    if (!speed_checked) {
        return speed_checked;
    }

    auto integrand = SymbolicExpr::multiply(f_composed.value(), speed_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(integrand, operation, "curve integrand construction is outside the supported domain", transform_budget);
    if (!integrand_checked) {
        return integrand_checked;
    }

    return vector_calculus_definite_integral_strict(integrand_checked.value(), t, a, b, context, operation, transform_budget);
}

VectorCalculusExprResult
vector_calculus_detail::curve_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < F.size(); ++i) {
        auto Fi_composed = vector_calculus_substitute_coords_strict(F[i], coord_vars, parametrization, operation, transform_budget);
        if (!Fi_composed) return Fi_composed;

        auto ri_prime = vector_calculus_differentiate_strict(parametrization[i], t, operation, transform_budget);
        if (!ri_prime) return ri_prime;

        auto term = SymbolicExpr::multiply(Fi_composed.value(), ri_prime.value());
        auto term_checked = vector_calculus_simplify_strict(term, operation, "curve integrand construction is outside the supported domain", transform_budget);
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(dot_product, operation, "curve integrand construction is outside the supported domain", transform_budget);
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    return vector_calculus_definite_integral_strict(dot_product, t, a, b, context, operation, transform_budget);
}

VectorCalculusExprResult
vector_calculus_detail::surface_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto f_composed = vector_calculus_substitute_coords_strict(f, coord_vars, parametrization, operation, transform_budget);
    if (!f_composed) return f_composed;

    auto cross = vector_calculus_cross_product_partials_strict(parametrization, u, v, operation, transform_budget);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    std::shared_ptr<SymbolicExpr> mag_sq;
    for (const auto& component : cross.value()) {
        auto sq = SymbolicExpr::power(component, SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(sq, operation, "surface magnitude construction is outside the supported domain", transform_budget);
        if (!sq_checked) return sq_checked;
        mag_sq = mag_sq ? SymbolicExpr::add(mag_sq, sq_checked.value())
                        : sq_checked.value();
        auto mag_sq_checked = vector_calculus_simplify_strict(mag_sq, operation, "surface magnitude construction is outside the supported domain", transform_budget);
        if (!mag_sq_checked) return mag_sq_checked;
        mag_sq = std::move(mag_sq_checked.value());
    }

    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    auto magnitude_checked = vector_calculus_simplify_strict(magnitude, operation, "surface magnitude construction is outside the supported domain", transform_budget);
    if (!magnitude_checked) return magnitude_checked;

    auto integrand = SymbolicExpr::multiply(f_composed.value(), magnitude_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(integrand, operation, "surface integrand construction is outside the supported domain", transform_budget);
    if (!integrand_checked) return integrand_checked;

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(integrand_checked.value(), steps, context, operation, transform_budget);
}

static VectorCalculusExprResult vector_calculus_surface_dot_strict(
    const VectorField& field,
    const std::vector<std::string>& coordinates,
    const VectorField& parametrization,
    const VectorField& normal,
    const std::string& message,
    const std::string& operation,
    detail::RewriteBudget& transform_budget) {
    std::shared_ptr<SymbolicExpr> dot_product;
    for (std::size_t i = 0; i < 3; ++i) {
        auto composed = vector_calculus_substitute_coords_strict(field[i], coordinates, parametrization, operation, transform_budget);
        if (!composed) {
            return composed;
        }
        auto term = SymbolicExpr::multiply(
            composed.value(), normal[i]);
        auto term_checked = vector_calculus_simplify_strict(term, operation, message, transform_budget);
        if (!term_checked) {
            return term_checked;
        }
        dot_product = dot_product
            ? SymbolicExpr::add(dot_product, term_checked.value())
            : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(dot_product, operation, message, transform_budget);
        if (!dot_checked) {
            return dot_checked;
        }
        dot_product = std::move(dot_checked.value());
    }
    return VectorCalculusExprResult::success(std::move(dot_product));
}

VectorCalculusExprResult
vector_calculus_detail::surface_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto cross = vector_calculus_cross_product_partials_strict(parametrization, u, v, operation, transform_budget);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    auto dot_product = vector_calculus_surface_dot_strict(F, coord_vars, parametrization, cross.value(),
    "surface integrand construction is outside the supported domain",
    operation, transform_budget);
    if (!dot_product) return dot_product;

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(dot_product.value(), steps, context, operation, transform_budget);
}

VectorCalculusExprResult vector_calculus_detail::greens_theorem_strict(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    auto dQ_dx = vector_calculus_differentiate_strict(Q, vars[0], operation, transform_budget);
    if (!dQ_dx) return dQ_dx;
    auto dP_dy = vector_calculus_differentiate_strict(P, vars[1], operation, transform_budget);
    if (!dP_dy) return dP_dy;

    auto integrand = SymbolicExpr::add(
        dQ_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dP_dy.value()));
    auto integrand_checked = vector_calculus_simplify_strict(integrand, operation, "Green's theorem integrand is outside the supported domain", transform_budget);
    if (!integrand_checked) return integrand_checked;

    std::vector<IntegrationStep> steps;
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(integrand_checked.value(), steps, context, operation, transform_budget);
}

static VectorCalculusExprResult vector_calculus_constant_integrand(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& variable,
    ComputationContext& context,
    detail::RewriteBudget& transform_budget) {
    auto candidate = detail::make_expression_ptr(detail::substitute_node(detail::node(integrand), variable, detail::node(SymbolicExpr::number(0)), transform_budget));
    candidate = candidate ? detail::make_expression_ptr(detail::normalize_node(detail::node(candidate), transform_budget)) : nullptr;
    if (!candidate || expression_depends_on_variable(detail::node(candidate), variable)) {
        return VectorCalculusExprResult::success(integrand);
    }

    EqvOptions trig_options;
    trig_options.profile = EqvProfile::TrigBasic;
    auto identity = check_equivalent(
        integrand, candidate, context, trig_options);
    if (!identity) {
        return VectorCalculusExprResult::failure(identity.error());
    }
    if (!std::holds_alternative<ProvedZeroResidual>(
            identity.value())) {
        return VectorCalculusExprResult::success(integrand);
    }
    return VectorCalculusExprResult::success(std::move(candidate));
}

VectorCalculusExprResult vector_calculus_detail::greens_theorem_area_strict(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    auto dx_dt = vector_calculus_differentiate_strict(parametrization[0], t, operation, transform_budget);
    if (!dx_dt) return dx_dt;
    auto dy_dt = vector_calculus_differentiate_strict(parametrization[1], t, operation, transform_budget);
    if (!dy_dt) return dy_dt;

    auto term1 = SymbolicExpr::multiply(parametrization[0], dy_dt.value());
    auto term2 = SymbolicExpr::multiply(parametrization[1], dx_dt.value());
    auto integrand = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    auto integrand_checked = vector_calculus_simplify_strict(integrand, operation, "Green's area integrand is outside the supported domain", transform_budget);
    if (!integrand_checked) return integrand_checked;
    integrand_checked = vector_calculus_constant_integrand(integrand_checked.value(), t, context, transform_budget);
    if (!integrand_checked) return integrand_checked;

    auto integral = vector_calculus_definite_integral_strict(integrand_checked.value(), t, a, b, context, operation, transform_budget);
    if (!integral) return integral;

    auto half = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral.value());
    return vector_calculus_simplify_strict(area, operation, "Green's area result is outside the supported domain", transform_budget);
}

VectorCalculusExprResult vector_calculus_detail::divergence_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    std::shared_ptr<SymbolicExpr> div_F;
    for (size_t i = 0; i < 3; ++i) {
        auto partial = vector_calculus_differentiate_strict(F[i], vars[i], operation, transform_budget);
        if (!partial) return partial;
        div_F = div_F ? SymbolicExpr::add(div_F, partial.value()) : partial.value();
        auto div_checked = vector_calculus_simplify_strict(div_F, operation, "divergence theorem integrand is outside the supported domain", transform_budget);
        if (!div_checked) return div_checked;
        div_F = std::move(div_checked.value());
    }

    std::vector<IntegrationStep> steps;
    steps.push_back({vars[2], z_bounds.first, z_bounds.second});
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(div_F, steps, context, operation, transform_budget);
}

static Result<VectorField> vector_calculus_curl_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    auto dFz_dy = vector_calculus_differentiate_strict(F[2], vars[1], operation, transform_budget);
    if (!dFz_dy) {
        return Result<VectorField>::failure(dFz_dy.error());
    }
    auto dFy_dz = vector_calculus_differentiate_strict(F[1], vars[2], operation, transform_budget);
    if (!dFy_dz) {
        return Result<VectorField>::failure(dFy_dz.error());
    }
    auto dFx_dz = vector_calculus_differentiate_strict(F[0], vars[2], operation, transform_budget);
    if (!dFx_dz) {
        return Result<VectorField>::failure(dFx_dz.error());
    }
    auto dFz_dx = vector_calculus_differentiate_strict(F[2], vars[0], operation, transform_budget);
    if (!dFz_dx) {
        return Result<VectorField>::failure(dFz_dx.error());
    }
    auto dFy_dx = vector_calculus_differentiate_strict(F[1], vars[0], operation, transform_budget);
    if (!dFy_dx) {
        return Result<VectorField>::failure(dFy_dx.error());
    }
    auto dFx_dy = vector_calculus_differentiate_strict(F[0], vars[1], operation, transform_budget);
    if (!dFx_dy) {
        return Result<VectorField>::failure(dFx_dy.error());
    }

    auto cx = SymbolicExpr::add(
        dFz_dy.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFy_dz.value()));
    auto cx_checked = vector_calculus_simplify_strict(cx, operation, "curl construction is outside the supported domain", transform_budget);
    if (!cx_checked) {
        return Result<VectorField>::failure(cx_checked.error());
    }

    auto cy = SymbolicExpr::add(
        dFx_dz.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFz_dx.value()));
    auto cy_checked = vector_calculus_simplify_strict(cy, operation, "curl construction is outside the supported domain", transform_budget);
    if (!cy_checked) {
        return Result<VectorField>::failure(cy_checked.error());
    }

    auto cz = SymbolicExpr::add(
        dFy_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFx_dy.value()));
    auto cz_checked = vector_calculus_simplify_strict(cz, operation, "curl construction is outside the supported domain", transform_budget);
    if (!cz_checked) {
        return Result<VectorField>::failure(cz_checked.error());
    }

    return Result<VectorField>::success(
        VectorField{std::move(cx_checked.value()),
                    std::move(cy_checked.value()),
                    std::move(cz_checked.value())});
}

VectorCalculusExprResult vector_calculus_detail::stokes_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context,
    const std::string& operation,
    detail::RewriteBudget& transform_budget)
{
    auto curl_F = vector_calculus_curl_strict(F, vars, operation, transform_budget);
    if (!curl_F) return VectorCalculusExprResult::failure(curl_F.error());

    auto cross = vector_calculus_cross_product_partials_strict(parametrization, u, v, operation, transform_budget);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    auto dot_product = vector_calculus_surface_dot_strict(curl_F.value(), vars, parametrization, cross.value(),
    "Stokes theorem integrand is outside the supported domain",
    operation, transform_budget);
    if (!dot_product) return dot_product;

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_bounds.first, v_bounds.second});
    steps.push_back({u, u_bounds.first, u_bounds.second});
    return vector_calculus_multiple_integral_strict(dot_product.value(), steps, context, operation, transform_budget);
}


}
