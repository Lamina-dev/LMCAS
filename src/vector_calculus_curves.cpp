#include "internal/expression_transform.hpp"
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
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "curve_integral_scalar");
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
        return curve_integral_scalar_strict(f, parametrization, t, a, b, context, operation, transform_budget);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
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
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
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
        f_composed = detail::make_expression_ptr(detail::substitute_node(detail::node(f_composed), coord_vars[i], detail::node(parametrization[i]), transform_budget));
    }
    f_composed = detail::make_expression_ptr(detail::normalize_node(detail::node(f_composed), transform_budget));

    std::shared_ptr<SymbolicExpr> speed_sq = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        auto deriv = detail::make_expression_ptr(detail::differentiate_node(detail::node(parametrization[i]), t, transform_budget));
        if (!deriv) {
            continue;
        }
        deriv = detail::make_expression_ptr(detail::normalize_node(detail::node(deriv), transform_budget));
        auto sq = SymbolicExpr::power(deriv, SymbolicExpr::number(2));
        sq = detail::make_expression_ptr(detail::normalize_node(detail::node(sq), transform_budget));
        if (!speed_sq) {
            speed_sq = sq;
        } else {
            speed_sq = SymbolicExpr::add(speed_sq, sq);
            speed_sq = detail::make_expression_ptr(detail::normalize_node(detail::node(speed_sq), transform_budget));
        }
    }

    if (!speed_sq) {
        return SymbolicExpr::number(0);
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    speed = detail::make_expression_ptr(detail::normalize_node(detail::node(speed), transform_budget));

    auto integrand = SymbolicExpr::multiply(f_composed, speed);
    integrand = detail::make_expression_ptr(detail::normalize_node(detail::node(integrand), transform_budget));


    return vector_calculus_integrate_with_fallback(integrand, t, a, b, context, transform_budget);
}


VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "curve_integral_vector");
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
        return curve_integral_vector_strict(F, parametrization, t, a, b, context, operation, transform_budget);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
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
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
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
            Fi_composed = detail::make_expression_ptr(detail::substitute_node(detail::node(Fi_composed), coord_vars[j], detail::node(parametrization[j]), transform_budget));
        }
        Fi_composed = detail::make_expression_ptr(detail::normalize_node(detail::node(Fi_composed), transform_budget));

        auto ri_prime = detail::make_expression_ptr(detail::differentiate_node(detail::node(parametrization[i]), t, transform_budget));
        if (!ri_prime) {
            continue;
        }
        ri_prime = detail::make_expression_ptr(detail::normalize_node(detail::node(ri_prime), transform_budget));

        auto term = SymbolicExpr::multiply(Fi_composed, ri_prime);
        term = detail::make_expression_ptr(detail::normalize_node(detail::node(term), transform_budget));

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = detail::make_expression_ptr(detail::normalize_node(detail::node(dot_product), transform_budget));
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }


    return vector_calculus_integrate_with_fallback(dot_product, t, a, b, context, transform_budget);
}

}
