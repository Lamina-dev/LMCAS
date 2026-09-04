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
            f, parametrization, t, a, b, context, operation);
    } catch (const detail::ResultPropagation& propagation) {
        return VectorCalculusExprResult::failure(propagation.error());
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
    ComputationContext context;
    return vector_calculus_integrate_with_fallback(
        integrand, t, a, b, context);
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
            F, parametrization, t, a, b, context, operation);
    } catch (const detail::ResultPropagation& propagation) {
        return VectorCalculusExprResult::failure(propagation.error());
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
    ComputationContext context;
    return vector_calculus_integrate_with_fallback(
        dot_product, t, a, b, context);
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
        return greens_theorem_area_strict(
            parametrization, t, a, b, context, operation);
    } catch (const detail::ResultPropagation& propagation) {
        return VectorCalculusExprResult::failure(propagation.error());
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
    ComputationContext context;
    auto integral_result = vector_calculus_integrate_with_fallback(
        integrand, t, a, b, context);
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
