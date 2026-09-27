#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_theorems.hpp"
#include "vector_calculus_differential.hpp"
#include "integration.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

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
    if (!vars_valid) {
        return VectorCalculusExprResult::failure(vars_valid.error());
    }
    if (!P || !LMCAS::detail::node(P) || !Q || !LMCAS::detail::node(Q)) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's theorem vector-field components cannot be null",
            operation);
    }
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) {
        return VectorCalculusExprResult::failure(x_valid.error());
    }
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) {
        return VectorCalculusExprResult::failure(y_valid.error());
    }
    auto step = context.consume_steps(10, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

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
    return LMCAS::detail::make_expression_ptr(result.value());
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    if (parametrization.size() != 2) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's area parametrization must be two-dimensional",
            operation);
    }
    auto step = context.consume_steps(10, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

    try {
        return greens_theorem_area_strict(
            parametrization, t, a, b, context, operation);
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

    ComputationContext context;
    auto area = greens_theorem_area_checked(
        parametrization, t, a, b, context);
    return area ? std::move(area.value()) : nullptr;
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
    if (!field_valid) {
        return VectorCalculusExprResult::failure(field_valid.error());
    }
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) {
        return VectorCalculusExprResult::failure(vars_valid.error());
    }
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) {
        return VectorCalculusExprResult::failure(x_valid.error());
    }
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) {
        return VectorCalculusExprResult::failure(y_valid.error());
    }
    auto z_valid = vector_calculus_validate_bound_pair(z_bounds, context, operation, "z");
    if (!z_valid) {
        return VectorCalculusExprResult::failure(z_valid.error());
    }
    auto step = context.consume_steps(14, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

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
    return LMCAS::detail::make_expression_ptr(result.value());
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
    if (!field_valid) {
        return VectorCalculusExprResult::failure(field_valid.error());
    }
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) {
        return VectorCalculusExprResult::failure(vars_valid.error());
    }
    auto param_valid = vector_calculus_validate_surface_parametrization(
        parametrization, u, v, u_bounds.first, u_bounds.second,
        v_bounds.first, v_bounds.second, context, operation);
    if (!param_valid) {
        return VectorCalculusExprResult::failure(param_valid.error());
    }
    auto step = context.consume_steps(18, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

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
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;

    for (size_t i = 0; i < 3; ++i) {
        if (!curl_F[i]) {
            continue;
        }

        /// 将 (∇×F)ᵢ 中的坐标变量替换为参数化表达式
        auto curl_i_composed = curl_F[i];
        for (size_t j = 0; j < 3; ++j) {
            curl_i_composed = curl_i_composed->substitute(vars[j], parametrization[j]);
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
    return LMCAS::detail::make_expression_ptr(result.value());
}

}
