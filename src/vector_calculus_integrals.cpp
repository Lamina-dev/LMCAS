#include "internal/vector_calculus_support.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "symbolic_ast.hpp"
#include "residual_verification.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace lamina {

using namespace vector_calculus_detail;

static std::vector<std::string> vector_calculus_coord_vars(size_t dim)
{
    if (dim == 2) return {"x", "y"};
    return {"x", "y", "z"};
}

/**
 * @internal
 * @brief 尝试符号定积分，若结果仍含未求值积分节点则返回 nullptr。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_try_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    Integrator integrator;
    SymbolicExpr result = detail::propagate_result(
        integrator.integrate_def(*integrand, var, *a, *b));

    if (vector_calculus_contains_unevaluated_integral(lamina::detail::node(result))) {
        return nullptr;
    }
    auto res = lamina::detail::make_expression_ptr(result);
    auto simplified = res->simplify();
    if (simplified && vector_calculus_contains_unevaluated_integral(lamina::detail::node(simplified))) {
        return nullptr;
    }
    return simplified ? simplified : res;
}

/**
 * @internal
 * @brief 数值定积分回退（复合 Simpson 法）。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_numerical_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    double a_val = 0.0;
    double b_val = 0.0;
    if (!vector_calculus_checked_finite_numeric(a, a_val, &context) ||
        !vector_calculus_checked_finite_numeric(b, b_val, &context)) {
        return nullptr;
    }

    int n = 1000;
    double h = (b_val - a_val) / n;
    if (!std::isfinite(h)) return nullptr;
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        double xi = a_val + i * h;
        if (!std::isfinite(xi)) return nullptr;
        auto xi_expr = SymbolicExpr::number(xi);
        auto fi = integrand->substitute(var, xi_expr);
        double fi_val = 0.0;
        if (!vector_calculus_checked_finite_numeric(fi, fi_val, &context)) {
            return nullptr;
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
        if (!std::isfinite(sum)) return nullptr;
    }
    sum *= h / 3.0;
    if (!std::isfinite(sum)) return nullptr;

    return SymbolicExpr::number(sum);
}

/**
 * @internal
 * @brief 符号积分优先，失败时回退到数值积分。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    auto symbolic_result = vector_calculus_try_definite(integrand, var, a, b);
    if (symbolic_result) {
        return symbolic_result;
    }
    return vector_calculus_numerical_definite(integrand, var, a, b);
}

static VectorCalculusExprResult vector_calculus_inconclusive(
    const std::string& operation,
    const std::string& message)
{
    return VectorCalculusExprResult::failure(CasErrc::Inconclusive,
                                             message, operation);
}

VectorCalculusExprResult vector_calculus_simplify_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message)
{
    if (!expr || !lamina::detail::node(expr)) {
        return vector_calculus_inconclusive(operation, message);
    }
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return vector_calculus_inconclusive(operation, message);
    }
    return VectorCalculusExprResult::success(std::move(simplified));
}

VectorCalculusExprResult vector_calculus_differentiate_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr)) {
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
    try {
        auto derivative = expr->differentiate(var);
        return vector_calculus_simplify_strict(
            derivative, operation,
            "expression derivative is outside the supported domain");
    } catch (const std::exception&) {
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
}

static VectorCalculusExprResult vector_calculus_substitute_coords_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& coord_vars,
    const VectorField& values,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr) || coord_vars.size() != values.size()) {
        return vector_calculus_inconclusive(
            operation, "coordinate substitution is outside the supported domain");
    }
    auto substituted = expr;
    for (size_t i = 0; i < coord_vars.size(); ++i) {
        if (!values[i] || !lamina::detail::node(values[i])) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
        substituted = substituted->substitute(coord_vars[i], values[i]);
        if (!substituted || !lamina::detail::node(substituted)) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
    }
    return vector_calculus_simplify_strict(
        substituted, operation,
        "coordinate substitution is outside the supported domain");
}

static VectorCalculusExprResult vector_calculus_definite_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    auto result = vector_calculus_try_definite(integrand, var, a, b);
    if (!result || !lamina::detail::node(result)) {
        return vector_calculus_inconclusive(
            operation, "integral could not be evaluated exactly in the supported domain");
    }
    return VectorCalculusExprResult::success(std::move(result));
}

static VectorCalculusExprResult vector_calculus_multiple_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::vector<IntegrationStep>& steps,
    ComputationContext& context,
    const std::string& operation)
{
    if (!integrand || !lamina::detail::node(integrand)) {
        return vector_calculus_inconclusive(
            operation, "integrand construction is outside the supported domain");
    }
    Integrator integrator;
    auto integrated = integrate_multiple_checked(*integrand, steps, integrator, context);
    if (!integrated) return VectorCalculusExprResult::failure(integrated.error());
    auto result = lamina::detail::make_expression_ptr(integrated.value());
    return vector_calculus_simplify_strict(
        result, operation,
        "integral result simplification is outside the supported domain");
}

static Result<VectorField> vector_calculus_cross_product_partials_strict(
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::string& operation)
{
    VectorField r_u;
    VectorField r_v;
    r_u.reserve(3);
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto du = vector_calculus_differentiate_strict(
            parametrization[i], u, operation);
        if (!du) return Result<VectorField>::failure(du.error());
        auto dv = vector_calculus_differentiate_strict(
            parametrization[i], v, operation);
        if (!dv) return Result<VectorField>::failure(dv.error());
        r_u.push_back(std::move(du.value()));
        r_v.push_back(std::move(dv.value()));
    }

    auto cross_x = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[1], r_v[2]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[2], r_v[1])));
    auto cross_x_checked = vector_calculus_simplify_strict(
        cross_x, operation, "surface normal construction is outside the supported domain");
    if (!cross_x_checked) return Result<VectorField>::failure(cross_x_checked.error());

    auto cross_y = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[2], r_v[0]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[0], r_v[2])));
    auto cross_y_checked = vector_calculus_simplify_strict(
        cross_y, operation, "surface normal construction is outside the supported domain");
    if (!cross_y_checked) return Result<VectorField>::failure(cross_y_checked.error());

    auto cross_z = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[0], r_v[1]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[1], r_v[0])));
    auto cross_z_checked = vector_calculus_simplify_strict(
        cross_z, operation, "surface normal construction is outside the supported domain");
    if (!cross_z_checked) return Result<VectorField>::failure(cross_z_checked.error());

    return Result<VectorField>::success(
        VectorField{std::move(cross_x_checked.value()),
                    std::move(cross_y_checked.value()),
                    std::move(cross_z_checked.value())});
}

static VectorCalculusExprResult curve_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    auto f_composed = vector_calculus_substitute_coords_strict(
        f, coord_vars, parametrization, operation);
    if (!f_composed) return f_composed;

    std::shared_ptr<SymbolicExpr> speed_sq;
    for (const auto& component : parametrization) {
        auto derivative = vector_calculus_differentiate_strict(
            component, t, operation);
        if (!derivative) return derivative;
        auto sq = SymbolicExpr::power(derivative.value(), SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(
            sq, operation, "curve speed construction is outside the supported domain");
        if (!sq_checked) return sq_checked;
        speed_sq = speed_sq ? SymbolicExpr::add(speed_sq, sq_checked.value())
                            : sq_checked.value();
        auto speed_sq_checked = vector_calculus_simplify_strict(
            speed_sq, operation, "curve speed construction is outside the supported domain");
        if (!speed_sq_checked) return speed_sq_checked;
        speed_sq = std::move(speed_sq_checked.value());
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    auto speed_checked = vector_calculus_simplify_strict(
        speed, operation, "curve speed construction is outside the supported domain");
    if (!speed_checked) return speed_checked;

    auto integrand = SymbolicExpr::multiply(f_composed.value(), speed_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "curve integrand construction is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    return vector_calculus_definite_integral_strict(
        integrand_checked.value(), t, a, b, operation);
}

static VectorCalculusExprResult curve_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < F.size(); ++i) {
        auto Fi_composed = vector_calculus_substitute_coords_strict(
            F[i], coord_vars, parametrization, operation);
        if (!Fi_composed) return Fi_composed;

        auto ri_prime = vector_calculus_differentiate_strict(
            parametrization[i], t, operation);
        if (!ri_prime) return ri_prime;

        auto term = SymbolicExpr::multiply(Fi_composed.value(), ri_prime.value());
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "curve integrand construction is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "curve integrand construction is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    return vector_calculus_definite_integral_strict(
        dot_product, t, a, b, operation);
}

static VectorCalculusExprResult surface_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto f_composed = vector_calculus_substitute_coords_strict(
        f, coord_vars, parametrization, operation);
    if (!f_composed) return f_composed;

    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    std::shared_ptr<SymbolicExpr> mag_sq;
    for (const auto& component : cross.value()) {
        auto sq = SymbolicExpr::power(component, SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(
            sq, operation, "surface magnitude construction is outside the supported domain");
        if (!sq_checked) return sq_checked;
        mag_sq = mag_sq ? SymbolicExpr::add(mag_sq, sq_checked.value())
                        : sq_checked.value();
        auto mag_sq_checked = vector_calculus_simplify_strict(
            mag_sq, operation, "surface magnitude construction is outside the supported domain");
        if (!mag_sq_checked) return mag_sq_checked;
        mag_sq = std::move(mag_sq_checked.value());
    }

    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    auto magnitude_checked = vector_calculus_simplify_strict(
        magnitude, operation, "surface magnitude construction is outside the supported domain");
    if (!magnitude_checked) return magnitude_checked;

    auto integrand = SymbolicExpr::multiply(f_composed.value(), magnitude_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "surface integrand construction is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(
        integrand_checked.value(), steps, context, operation);
}

static VectorCalculusExprResult surface_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < 3; ++i) {
        auto Fi_composed = vector_calculus_substitute_coords_strict(
            F[i], coord_vars, parametrization, operation);
        if (!Fi_composed) return Fi_composed;

        auto term = SymbolicExpr::multiply(Fi_composed.value(), cross.value()[i]);
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "surface integrand construction is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "surface integrand construction is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(dot_product, steps, context, operation);
}

static VectorCalculusExprResult greens_theorem_strict(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context,
    const std::string& operation)
{
    auto dQ_dx = vector_calculus_differentiate_strict(Q, vars[0], operation);
    if (!dQ_dx) return dQ_dx;
    auto dP_dy = vector_calculus_differentiate_strict(P, vars[1], operation);
    if (!dP_dy) return dP_dy;

    auto integrand = SymbolicExpr::add(
        dQ_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dP_dy.value()));
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "Green's theorem integrand is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    std::vector<IntegrationStep> steps;
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(
        integrand_checked.value(), steps, context, operation);
}

static VectorCalculusExprResult greens_theorem_area_strict(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    auto dx_dt = vector_calculus_differentiate_strict(
        parametrization[0], t, operation);
    if (!dx_dt) return dx_dt;
    auto dy_dt = vector_calculus_differentiate_strict(
        parametrization[1], t, operation);
    if (!dy_dt) return dy_dt;

    auto term1 = SymbolicExpr::multiply(parametrization[0], dy_dt.value());
    auto term2 = SymbolicExpr::multiply(parametrization[1], dx_dt.value());
    auto integrand = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "Green's area integrand is outside the supported domain");
    if (!integrand_checked) return integrand_checked;
    lsr::EqvOptions trig_options;
    trig_options.profile = lsr::EqvProfile::TrigBasic;
    ComputationContext identity_context;
    auto unit_identity = check_equivalent(
        integrand_checked.value(), SymbolicExpr::number(1),
        identity_context, trig_options);
    if (unit_identity &&
        std::holds_alternative<ProvedZeroResidual>(
            unit_identity.value())) {
        integrand_checked = VectorCalculusExprResult::success(
            SymbolicExpr::number(1));
    }

    auto integral = vector_calculus_definite_integral_strict(
        integrand_checked.value(), t, a, b, operation);
    if (!integral) return integral;

    auto half = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral.value());
    return vector_calculus_simplify_strict(
        area, operation, "Green's area result is outside the supported domain");
}

static VectorCalculusExprResult divergence_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context,
    const std::string& operation)
{
    std::shared_ptr<SymbolicExpr> div_F;
    for (size_t i = 0; i < 3; ++i) {
        auto partial = vector_calculus_differentiate_strict(F[i], vars[i], operation);
        if (!partial) return partial;
        div_F = div_F ? SymbolicExpr::add(div_F, partial.value()) : partial.value();
        auto div_checked = vector_calculus_simplify_strict(
            div_F, operation, "divergence theorem integrand is outside the supported domain");
        if (!div_checked) return div_checked;
        div_F = std::move(div_checked.value());
    }

    std::vector<IntegrationStep> steps;
    steps.push_back({vars[2], z_bounds.first, z_bounds.second});
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(div_F, steps, context, operation);
}

static Result<VectorField> vector_calculus_curl_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::string& operation)
{
    auto dFz_dy = vector_calculus_differentiate_strict(F[2], vars[1], operation);
    if (!dFz_dy) return Result<VectorField>::failure(dFz_dy.error());
    auto dFy_dz = vector_calculus_differentiate_strict(F[1], vars[2], operation);
    if (!dFy_dz) return Result<VectorField>::failure(dFy_dz.error());
    auto dFx_dz = vector_calculus_differentiate_strict(F[0], vars[2], operation);
    if (!dFx_dz) return Result<VectorField>::failure(dFx_dz.error());
    auto dFz_dx = vector_calculus_differentiate_strict(F[2], vars[0], operation);
    if (!dFz_dx) return Result<VectorField>::failure(dFz_dx.error());
    auto dFy_dx = vector_calculus_differentiate_strict(F[1], vars[0], operation);
    if (!dFy_dx) return Result<VectorField>::failure(dFy_dx.error());
    auto dFx_dy = vector_calculus_differentiate_strict(F[0], vars[1], operation);
    if (!dFx_dy) return Result<VectorField>::failure(dFx_dy.error());

    auto cx = SymbolicExpr::add(
        dFz_dy.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFy_dz.value()));
    auto cx_checked = vector_calculus_simplify_strict(
        cx, operation, "curl construction is outside the supported domain");
    if (!cx_checked) return Result<VectorField>::failure(cx_checked.error());

    auto cy = SymbolicExpr::add(
        dFx_dz.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFz_dx.value()));
    auto cy_checked = vector_calculus_simplify_strict(
        cy, operation, "curl construction is outside the supported domain");
    if (!cy_checked) return Result<VectorField>::failure(cy_checked.error());

    auto cz = SymbolicExpr::add(
        dFy_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFx_dy.value()));
    auto cz_checked = vector_calculus_simplify_strict(
        cz, operation, "curl construction is outside the supported domain");
    if (!cz_checked) return Result<VectorField>::failure(cz_checked.error());

    return Result<VectorField>::success(
        VectorField{std::move(cx_checked.value()),
                    std::move(cy_checked.value()),
                    std::move(cz_checked.value())});
}

static VectorCalculusExprResult stokes_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context,
    const std::string& operation)
{
    auto curl_F = vector_calculus_curl_strict(F, vars, operation);
    if (!curl_F) return VectorCalculusExprResult::failure(curl_F.error());

    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < 3; ++i) {
        auto curl_i_composed = vector_calculus_substitute_coords_strict(
            curl_F.value()[i], coord_vars, parametrization, operation);
        if (!curl_i_composed) return curl_i_composed;

        auto term = SymbolicExpr::multiply(curl_i_composed.value(), cross.value()[i]);
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "Stokes theorem integrand is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "Stokes theorem integrand is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    std::vector<IntegrationStep> steps;
    steps.push_back({v, v_bounds.first, v_bounds.second});
    steps.push_back({u, u_bounds.first, u_bounds.second});
    return vector_calculus_multiple_integral_strict(dot_product, steps, context, operation);
}


VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "curve_integral_scalar";
    auto valid = vector_calculus_validate_curve_scalar_inputs(
        f, parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(parametrization.size() * 4 + 4, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return curve_integral_scalar_strict(
            f, parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return curve_integral_scalar_checked(f, parametrization, t, a, b, context);
}

std::shared_ptr<SymbolicExpr> curve_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (!f) {
        throw std::invalid_argument("curve_integral_scalar: f must not be null");
    }
    if (parametrization.empty()) {
        throw std::invalid_argument("curve_integral_scalar: parametrization must not be empty");
    }
    if (!a || !b) {
        throw std::invalid_argument("curve_integral_scalar: bounds must not be null");
    }

    size_t dim = parametrization.size();
    auto coord_vars = vector_calculus_coord_vars(dim);

    /// 将 f 中的坐标变量替换为参数化表达式: f(r(t))
    auto f_composed = f;
    for (size_t i = 0; i < dim; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    /// 计算 |r'(t)| = √(∑(r_i'(t))²)
    std::shared_ptr<SymbolicExpr> speed_sq = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        auto deriv = parametrization[i]->differentiate(t);
        if (!deriv) continue;
        deriv = deriv->simplify();
        auto sq = SymbolicExpr::power(deriv, SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!speed_sq) {
            speed_sq = sq;
        } else {
            speed_sq = SymbolicExpr::add(speed_sq, sq);
            speed_sq = speed_sq->simplify();
        }
    }

    if (!speed_sq) {
        return SymbolicExpr::number(0);
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    speed = speed->simplify();

    /// 被积函数: f(r(t)) · |r'(t)|
    auto integrand = SymbolicExpr::multiply(f_composed, speed);
    integrand = integrand->simplify();

    /// 计算定积分 ∫ₐᵇ f(r(t))·|r'(t)| dt
    return vector_calculus_integrate_with_fallback(integrand, t, a, b);
}


VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "curve_integral_vector";
    auto valid = vector_calculus_validate_curve_vector_inputs(
        F, parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 3 + 4,
                                      operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return curve_integral_vector_strict(
            F, parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return curve_integral_vector_checked(F, parametrization, t, a, b, context);
}

std::shared_ptr<SymbolicExpr> curve_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (F.empty()) {
        throw std::invalid_argument("curve_integral_vector: F must not be empty");
    }
    if (parametrization.empty()) {
        throw std::invalid_argument("curve_integral_vector: parametrization must not be empty");
    }
    if (F.size() != parametrization.size()) {
        throw std::invalid_argument(
            "curve_integral_vector: F and parametrization must have the same dimension");
    }
    if (!a || !b) {
        throw std::invalid_argument("curve_integral_vector: bounds must not be null");
    }

    size_t dim = parametrization.size();
    auto coord_vars = vector_calculus_coord_vars(dim);

    /// 计算 F(r(t))·r'(t) = ∑ Fᵢ(r(t)) · rᵢ'(t)
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        if (!F[i]) continue;

        /// 将 Fᵢ 中的坐标变量替换为参数化表达式
        auto Fi_composed = F[i];
        for (size_t j = 0; j < dim; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

        /// rᵢ'(t)
        auto ri_prime = parametrization[i]->differentiate(t);
        if (!ri_prime) continue;
        ri_prime = ri_prime->simplify();

        /// Fᵢ(r(t)) · rᵢ'(t)
        auto term = SymbolicExpr::multiply(Fi_composed, ri_prime);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    /// 计算定积分 ∫ₐᵇ F(r(t))·r'(t) dt
    return vector_calculus_integrate_with_fallback(dot_product, t, a, b);
}


/**
 * @internal
 * @brief 计算参数化曲面的法向量 r_u × r_v。
 *
 * 参数化曲面 r(u,v) = [x(u,v), y(u,v), z(u,v)]，
 * 返回叉积 r_u × r_v 的三个分量。
 */
static VectorField vector_calculus_cross_product_partials(
    const VectorField& parametrization,
    const std::string& u, const std::string& v)
{
    /// r_u = [∂x/∂u, ∂y/∂u, ∂z/∂u]
    VectorField r_u;
    r_u.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(u);
        r_u.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    /// r_v = [∂x/∂v, ∂y/∂v, ∂z/∂v]
    VectorField r_v;
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(v);
        r_v.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    /// 叉积: r_u × r_v = (r_u[1]*r_v[2] - r_u[2]*r_v[1],
    ///                     r_u[2]*r_v[0] - r_u[0]*r_v[2],
    ///                     r_u[0]*r_v[1] - r_u[1]*r_v[0])
    auto cross_x = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[1], r_v[2]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[2], r_v[1])));
    cross_x = cross_x->simplify();

    auto cross_y = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[2], r_v[0]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[0], r_v[2])));
    cross_y = cross_y->simplify();

    auto cross_z = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[0], r_v[1]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[1], r_v[0])));
    cross_z = cross_z->simplify();

    return VectorField{cross_x, cross_y, cross_z};
}


VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context)
{
    const std::string operation = "surface_integral_scalar";
    auto valid = vector_calculus_validate_surface_scalar_inputs(
        f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(parametrization.size() * 6 + 8, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return surface_integral_scalar_strict(
            f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
            context, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "surface integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    ComputationContext context;
    return surface_integral_scalar_checked(
        f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper, context);
}

std::shared_ptr<SymbolicExpr> surface_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    if (!f) {
        throw std::invalid_argument("surface_integral_scalar: f must not be null");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_scalar: parametrization must have 3 components");
    }
    if (!u_lower || !u_upper || !v_lower || !v_upper) {
        throw std::invalid_argument("surface_integral_scalar: bounds must not be null");
    }

    std::vector<std::string> coord_vars = {"x", "y", "z"};

    /// 将 f 中的坐标变量替换为参数化表达式: f(r(u,v))
    auto f_composed = f;
    for (size_t i = 0; i < 3; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    /// 计算 r_u × r_v
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// |r_u × r_v| = √(cross_x² + cross_y² + cross_z²)
    std::shared_ptr<SymbolicExpr> mag_sq = nullptr;
    for (size_t i = 0; i < 3; ++i) {
        auto sq = SymbolicExpr::power(cross[i], SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!mag_sq) {
            mag_sq = sq;
        } else {
            mag_sq = SymbolicExpr::add(mag_sq, sq);
            mag_sq = mag_sq->simplify();
        }
    }

    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    magnitude = magnitude->simplify();

    /// 被积函数: f(r(u,v)) · |r_u × r_v|
    auto integrand = SymbolicExpr::multiply(f_composed, magnitude);
    integrand = integrand->simplify();

    Integrator integrator;
    ComputationContext context;
    std::vector<IntegrationStep> steps = {
        {v, v_lower, v_upper},
        {u, u_lower, u_upper},
    };
    auto result = integrate_multiple_checked(*integrand, steps, integrator, context);
    if (!result) throw std::runtime_error(result.error().message);
    return lamina::detail::make_expression_ptr(result.value());
}


VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context)
{
    const std::string operation = "surface_integral_vector";
    auto valid = vector_calculus_validate_surface_vector_inputs(
        F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 6 + 8,
                                      operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return surface_integral_vector_strict(
            F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
            context, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "surface integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    ComputationContext context;
    return surface_integral_vector_checked(
        F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper, context);
}

std::shared_ptr<SymbolicExpr> surface_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    if (F.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_vector: F must have 3 components");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_vector: parametrization must have 3 components");
    }
    if (!u_lower || !u_upper || !v_lower || !v_upper) {
        throw std::invalid_argument("surface_integral_vector: bounds must not be null");
    }

    std::vector<std::string> coord_vars = {"x", "y", "z"};

    /// 计算 r_u × r_v
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// 将 F 中的坐标变量替换为参数化表达式，然后计算 F·(r_u × r_v)
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < 3; ++i) {
        if (!F[i]) continue;

        /// 将 Fᵢ 中的坐标变量替换为参数化表达式
        auto Fi_composed = F[i];
        for (size_t j = 0; j < 3; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

        /// Fᵢ(r(u,v)) · (r_u × r_v)ᵢ
        auto term = SymbolicExpr::multiply(Fi_composed, cross[i]);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    Integrator integrator;
    ComputationContext context;
    std::vector<IntegrationStep> steps = {
        {v, v_lower, v_upper},
        {u, u_lower, u_upper},
    };
    auto result = integrate_multiple_checked(*dot_product, steps, integrator, context);
    if (!result) throw std::runtime_error(result.error().message);
    return lamina::detail::make_expression_ptr(result.value());
}


VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context)
{
    const std::string operation = "greens_theorem";
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 2, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    if (!P || !lamina::detail::node(P) || !Q || !lamina::detail::node(Q)) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's theorem vector-field components cannot be null",
            operation);
    }
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) return VectorCalculusExprResult::failure(x_valid.error());
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) return VectorCalculusExprResult::failure(y_valid.error());
    auto step = context.consume_steps(10, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return greens_theorem_strict(P, Q, vars, x_bounds, y_bounds, context, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Green's theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds)
{
    ComputationContext context;
    return greens_theorem_checked(P, Q, vars, x_bounds, y_bounds, context);
}

std::shared_ptr<SymbolicExpr> greens_theorem(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds)
{
    if (!P || !Q) {
        throw std::invalid_argument("greens_theorem: P and Q must not be null");
    }
    if (vars.size() != 2) {
        throw std::invalid_argument("greens_theorem: vars must have exactly 2 elements");
    }
    if (!x_bounds.first || !x_bounds.second ||
        !y_bounds.first || !y_bounds.second) {
        throw std::invalid_argument("greens_theorem: bounds must not be null");
    }

    const std::string& x_var = vars[0];
    const std::string& y_var = vars[1];

    /// 计算被积函数: ∂Q/∂x - ∂P/∂y
    auto dQ_dx = Q->differentiate(x_var);
    auto dP_dy = P->differentiate(y_var);

    if (!dQ_dx) dQ_dx = SymbolicExpr::number(0);
    if (!dP_dy) dP_dy = SymbolicExpr::number(0);

    dQ_dx = dQ_dx->simplify();
    dP_dy = dP_dy->simplify();

    auto integrand = SymbolicExpr::add(
        dQ_dx,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dP_dy));
    integrand = integrand->simplify();

    Integrator integrator;
    ComputationContext context;
    std::vector<IntegrationStep> steps = {
        {y_var, y_bounds.first, y_bounds.second},
        {x_var, x_bounds.first, x_bounds.second},
    };
    auto result = integrate_multiple_checked(*integrand, steps, integrator, context);
    if (!result) throw std::runtime_error(result.error().message);
    return lamina::detail::make_expression_ptr(result.value());
}


VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "greens_theorem_area";
    auto valid = vector_calculus_validate_curve_parametrization(
        parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    if (parametrization.size() != 2) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's area parametrization must be two-dimensional",
            operation);
    }
    auto step = context.consume_steps(10, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return greens_theorem_area_strict(parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Green's area allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return greens_theorem_area_checked(parametrization, t, a, b, context);
}

std::shared_ptr<SymbolicExpr> greens_theorem_area(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (parametrization.size() != 2) {
        throw std::invalid_argument(
            "greens_theorem_area: parametrization must have 2 components (x(t), y(t))");
    }
    if (!a || !b) {
        throw std::invalid_argument("greens_theorem_area: bounds must not be null");
    }

    auto x_t = parametrization[0];
    auto y_t = parametrization[1];

    if (!x_t || !y_t) {
        throw std::invalid_argument("greens_theorem_area: parametrization components must not be null");
    }

    /// A = (1/2) ∮ (x dy - y dx)
    // = (1/2) ∫ₐᵇ (x(t)·y'(t) - y(t)·x'(t)) dt
    auto dx_dt = x_t->differentiate(t);
    auto dy_dt = y_t->differentiate(t);

    if (!dx_dt) dx_dt = SymbolicExpr::number(0);
    if (!dy_dt) dy_dt = SymbolicExpr::number(0);

    dx_dt = dx_dt->simplify();
    dy_dt = dy_dt->simplify();

    /// x(t)·y'(t) - y(t)·x'(t)
    auto term1 = SymbolicExpr::multiply(x_t, dy_dt);
    auto term2 = SymbolicExpr::multiply(y_t, dx_dt);
    auto integrand = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    integrand = integrand->simplify();

    /// 计算定积分
    auto integral_result = vector_calculus_integrate_with_fallback(integrand, t, a, b);
    if (!integral_result) {
        return nullptr;
    }

    /// 乘以 1/2
    auto half = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral_result);
    area = area->simplify();

    return area;
}


VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context)
{
    const std::string operation = "divergence_theorem";
    auto field_valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!field_valid) return VectorCalculusExprResult::failure(field_valid.error());
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) return VectorCalculusExprResult::failure(x_valid.error());
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) return VectorCalculusExprResult::failure(y_valid.error());
    auto z_valid = vector_calculus_validate_bound_pair(z_bounds, context, operation, "z");
    if (!z_valid) return VectorCalculusExprResult::failure(z_valid.error());
    auto step = context.consume_steps(14, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return divergence_theorem_strict(
            F, vars, x_bounds, y_bounds, z_bounds, context, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "divergence theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds)
{
    ComputationContext context;
    return divergence_theorem_checked(F, vars, x_bounds, y_bounds, z_bounds, context);
}

std::shared_ptr<SymbolicExpr> divergence_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds)
{
    if (F.size() != 3) {
        throw std::invalid_argument("divergence_theorem: F must have 3 components");
    }
    if (vars.size() != 3) {
        throw std::invalid_argument("divergence_theorem: vars must have exactly 3 elements");
    }
    if (!x_bounds.first || !x_bounds.second ||
        !y_bounds.first || !y_bounds.second ||
        !z_bounds.first || !z_bounds.second) {
        throw std::invalid_argument("divergence_theorem: bounds must not be null");
    }

    /// 计算散度 ∇·F = ∂F₁/∂x + ∂F₂/∂y + ∂F₃/∂z
    auto div_F = divergence(F, vars);

    Integrator integrator;
    ComputationContext context;
    std::vector<IntegrationStep> steps = {
        {vars[2], z_bounds.first, z_bounds.second},
        {vars[1], y_bounds.first, y_bounds.second},
        {vars[0], x_bounds.first, x_bounds.second},
    };
    auto result = integrate_multiple_checked(*div_F, steps, integrator, context);
    if (!result) throw std::runtime_error(result.error().message);
    return lamina::detail::make_expression_ptr(result.value());
}


VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context)
{
    const std::string operation = "stokes_theorem";
    auto field_valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!field_valid) return VectorCalculusExprResult::failure(field_valid.error());
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    auto param_valid = vector_calculus_validate_surface_parametrization(
        parametrization, u, v, u_bounds.first, u_bounds.second,
        v_bounds.first, v_bounds.second, context, operation);
    if (!param_valid) return VectorCalculusExprResult::failure(param_valid.error());
    auto step = context.consume_steps(18, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return stokes_theorem_strict(
            F, vars, parametrization, u, v, u_bounds, v_bounds, context, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Stokes theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds)
{
    ComputationContext context;
    return stokes_theorem_checked(F, vars, parametrization, u, v,
                                  u_bounds, v_bounds, context);
}

std::shared_ptr<SymbolicExpr> stokes_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds)
{
    if (F.size() != 3) {
        throw std::invalid_argument("stokes_theorem: F must have 3 components");
    }
    if (vars.size() != 3) {
        throw std::invalid_argument("stokes_theorem: vars must have exactly 3 elements");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "stokes_theorem: parametrization must have 3 components");
    }
    if (!u_bounds.first || !u_bounds.second ||
        !v_bounds.first || !v_bounds.second) {
        throw std::invalid_argument("stokes_theorem: bounds must not be null");
    }

    /// 计算旋度 ∇×F
    auto curl_F = curl(F, vars);

    /// 计算 r_u × r_v（曲面法向量）
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// 将 (∇×F) 中的坐标变量替换为参数化表达式，然后计算 (∇×F)·(r_u × r_v)
    std::vector<std::string> coord_vars = {"x", "y", "z"};
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;

    for (size_t i = 0; i < 3; ++i) {
        if (!curl_F[i]) continue;

        /// 将 (∇×F)ᵢ 中的坐标变量替换为参数化表达式
        auto curl_i_composed = curl_F[i];
        for (size_t j = 0; j < 3; ++j) {
            curl_i_composed = curl_i_composed->substitute(coord_vars[j], parametrization[j]);
        }
        curl_i_composed = curl_i_composed->simplify();

        /// (∇×F)ᵢ(r(u,v)) · (r_u × r_v)ᵢ
        auto term = SymbolicExpr::multiply(curl_i_composed, cross[i]);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    Integrator integrator;
    ComputationContext context;
    std::vector<IntegrationStep> steps = {
        {v, v_bounds.first, v_bounds.second},
        {u, u_bounds.first, u_bounds.second},
    };
    auto result = integrate_multiple_checked(*dot_product, steps, integrator, context);
    if (!result) throw std::runtime_error(result.error().message);
    return lamina::detail::make_expression_ptr(result.value());
}
} // namespace lamina
