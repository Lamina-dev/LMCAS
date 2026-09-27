#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_surfaces.hpp"
#include "integration.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorField vector_calculus_detail::vector_calculus_cross_product_partials(
    const VectorField& parametrization,
    const std::string& u, const std::string& v)
{
    VectorField r_u;
    r_u.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(u);
        r_u.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    VectorField r_v;
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(v);
        r_v.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    VectorField result;
    result.reserve(3);
    for (std::size_t index = 0; index < 3; ++index) {
        result.push_back(
            vector_calculus_cross_component(r_u, r_v, index)->simplify());
    }
    return result;
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }

    auto step = context.consume_steps(parametrization.size() * 6 + 8, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

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

    auto f_composed = f;
    for (size_t i = 0; i < 3; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

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
    return LMCAS::detail::make_expression_ptr(result.value());
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 6 + 8,
                                      operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

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

    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < 3; ++i) {
        if (!F[i]) {
            continue;
        }

        auto Fi_composed = F[i];
        for (size_t j = 0; j < 3; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

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
    return LMCAS::detail::make_expression_ptr(result.value());
}

}
