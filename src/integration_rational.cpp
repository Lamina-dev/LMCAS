#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "internal/integration_rational_support.hpp"

namespace LMCAS {

namespace {
void rd_reduce_reciprocal(Polynomial<Rational>& P, Polynomial<Rational>& Q) {
    if (!P.is_zero() && P.degree() > 0 && P.degree() < Q.degree()) {
        auto [reduced_denominator, remainder] = Q.div_mod(P);
        if (remainder.is_zero()) {
            P = Polynomial<Rational>(
                std::vector<Rational>{Rational(1)}, P.variable_name);
            Q = std::move(reduced_denominator);
        }
    }
}

std::shared_ptr<SymbolicExpr> rd_integrate_polynomial(
    const Polynomial<Rational>& quot, const std::string& var) {
    auto result = SymbolicExpr::number(0);
        if (!quot.is_zero()) {
            for (size_t k = 0; k < quot.coeffs.size(); ++k) {
                if (quot.coeffs[k] == Rational(0)) { continue; }
                const BigInt integrated_power = BigInt(k) + BigInt(1);
                Rational coeff = quot.coeffs[k] / Rational(integrated_power);
                std::shared_ptr<SymbolicExpr> term;
                auto v = SymbolicExpr::variable(var);
                auto pw = SymbolicExpr::power(v, SymbolicExpr::number(integrated_power));
                if (coeff == Rational(1)) {
                    term = pw;
                } else {
                    term = SymbolicExpr::multiply(SymbolicExpr::number(coeff), pw);
                }
                result = SymbolicExpr::add(result, term);
            }
        }
    return result;
}

detail::IntegralParts rd_add_fraction_primitives(
    const std::string& var,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    const std::vector<Polynomial<Rational>>& numerators,
    detail::IntegralParts result, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                            context.limits().max_ast_nodes, "transform");
    size_t idx = 0;
    for (const auto& [factor, mult] : factors) {
        for (int power = 1; power <= mult; ++power, ++idx) {
            auto term = detail::integrate_rational_term(numerators[idx], factor, power, var);
            result.primitive = *SymbolicExpr::add(
                detail::make_expression_ptr(result.primitive),
                detail::make_expression_ptr(term.primitive));
            result.remainder = *SymbolicExpr::add(
                detail::make_expression_ptr(result.remainder),
                detail::make_expression_ptr(term.remainder));
        }
    }
    result.primitive = detail::expression_from_node(
        detail::normalize_node(detail::node(result.primitive), transform_budget));
    result.remainder = detail::expression_from_node(
        detail::normalize_node(detail::node(result.remainder), transform_budget));
    return result;
}

Result<detail::IntegralParts> rd_integrate_proper_fraction(
    const std::string& var,
    const Polynomial<Rational>& P, const Polynomial<Rational>& Q,
    const SymbolicExpr& zero, ComputationContext& computation) {
    Polynomial<Rational> quot(Q.variable_name), rem = P;
    if (P.degree() >= Q.degree()) {
        auto divided = P.div_mod(Q);
        quot = std::move(divided.first);
        rem = std::move(divided.second);
    }

    detail::IntegralParts result{*rd_integrate_polynomial(quot, var), zero};
    if (rem.is_zero()) { return result; }

    std::vector<std::pair<Polynomial<Rational>, int>> factors;
    if (!detail::factor_denominator(Q, factors)) {
        result.remainder = *SymbolicExpr::divide(poly_to_symbolic(rem), poly_to_symbolic(Q));
        return result;
    }

    std::vector<Polynomial<Rational>> numerators;
    auto coefficients = detail::solve_coefficients(rem, Q, factors, numerators, computation);
    if (!coefficients) {
        return Result<detail::IntegralParts>::failure(coefficients.error());
    }
    if (!coefficients.value()) {
        result.remainder = *SymbolicExpr::divide(poly_to_symbolic(rem), poly_to_symbolic(Q));
        return result;
    }
    return rd_add_fraction_primitives(var, factors, numerators, std::move(result), computation);
}
}

Result<std::optional<detail::IntegralParts>> Integrator::try_rational_decomposition(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int depth) try {
    detail::RewriteBudget transform_budget(computation, computation.limits().max_recursion_depth,
                                            computation.limits().max_ast_nodes, "transform");
    (void)depth;

    Polynomial<Rational> P, Q;
    auto extracted = detail::extract_rational(expr, var, P, Q);
    if (!extracted) {
        return Result<std::optional<detail::IntegralParts>>::failure(extracted.error());
    }
    if (!extracted.value()) { return std::optional<detail::IntegralParts>{}; }
    rd_reduce_reciprocal(P, Q);
    if (Q.is_zero() || Q.degree() < 1) { return std::optional<detail::IntegralParts>{}; }

    if (Q.degree() == 1 && P.degree() <= 0 &&
        !P.coeffs.empty() && Q.coeffs.size() > 1 &&
        Q.coeffs[1] != Rational(0)) {
        auto coefficient = SymbolicExpr::number(P.coeffs[0] / Q.coeffs[1]);
        auto logarithm = SymbolicExpr::ln(poly_to_symbolic(Q));
        auto primitive = detail::expression_from_node(detail::normalize_node(
            detail::node(SymbolicExpr::multiply(coefficient, logarithm)), transform_budget));
        return std::optional<detail::IntegralParts>{{std::move(primitive), zero_}};
    }
    auto integrated = rd_integrate_proper_fraction(var, P, Q, zero_, computation);
    if (!integrated) {
        return Result<std::optional<detail::IntegralParts>>::failure(integrated.error());
    }
    return std::optional<detail::IntegralParts>{std::move(integrated.value())};
} catch (const CasError& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::optional<detail::IntegralParts>>::failure(
        CasErrc::ResourceLimit, "allocation failed", "integrate.rational");
} catch (const std::length_error& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(
        CasErrc::ResourceLimit, error.what(), "integrate.rational");
} catch (const std::exception& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(
        CasErrc::InternalInvariant, error.what(), "integrate.rational");
}

}
