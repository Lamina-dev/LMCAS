#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_curves.hpp"
#include "integration.hpp"
#include <stdexcept>

namespace LMCAS {

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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }

    auto step = context.consume_steps(parametrization.size() * 4 + 4, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

    try {
        return curve_integral_scalar_strict(
            f, parametrization, t, a, b, context, operation);
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

    auto f_composed = f;
    for (size_t i = 0; i < dim; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    std::shared_ptr<SymbolicExpr> speed_sq = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        auto deriv = parametrization[i]->differentiate(t);
        if (!deriv) {
            continue;
        }
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

    auto integrand = SymbolicExpr::multiply(f_composed, speed);
    integrand = integrand->simplify();

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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 3 + 4,
                                      operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }

    try {
        return curve_integral_vector_strict(
            F, parametrization, t, a, b, context, operation);
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

    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        if (!F[i]) {
            continue;
        }

        auto Fi_composed = F[i];
        for (size_t j = 0; j < dim; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

        auto ri_prime = parametrization[i]->differentiate(t);
        if (!ri_prime) {
            continue;
        }
        ri_prime = ri_prime->simplify();

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

    ComputationContext context;
    return vector_calculus_integrate_with_fallback(
        dot_product, t, a, b, context);
}

}
