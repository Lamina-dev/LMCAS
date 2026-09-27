#include "internal/integration_rational_support.hpp"

namespace LMCAS {

namespace {
void rd_reduce_reciprocal(RationalDecompositionStrategy& strategy,
    Polynomial<Rational>& P, Polynomial<Rational>& Q) {
    if (!P.is_zero() && P.degree() > 0 && P.degree() < Q.degree()) {
        Polynomial<Rational> reduced_denominator;
        Polynomial<Rational> remainder;
        strategy.poly_divide(Q, P, reduced_denominator, remainder);
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

std::shared_ptr<SymbolicExpr> rd_add_fraction_primitives(
    RationalDecompositionStrategy& strategy, const SymbolicExpr& expr, const std::string& var,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    const std::vector<Polynomial<Rational>>& numerators,
    std::shared_ptr<SymbolicExpr> result) {
        if (!numerators.empty()) {
            size_t idx = 0;
            for (size_t i = 0; i < factors.size(); ++i) {
                const auto& [fpoly, mult] = factors[i];
                for (int l = 1; l <= mult; ++l, ++idx) {
                    if (idx >= numerators.size()) { break; }
                    auto term = strategy.integrate_term(numerators[idx], fpoly, l, var);
                    if (!term) {
                        return LMCAS::detail::make_expression_ptr(
                            LMCAS::detail::make_node<IntegralNode>(
                                LMCAS::detail::node(expr), var));
                    }
                    result = SymbolicExpr::add(result, term);
                }
            }
        }
        auto simplified = result->simplify();
        if (simplified) {
            return simplified;
        }
        return result;
}

Result<std::shared_ptr<SymbolicExpr>> rd_integrate_proper_fraction(
    RationalDecompositionStrategy& strategy, const SymbolicExpr& expr, const std::string& var,
    const Polynomial<Rational>& P, const Polynomial<Rational>& Q,
    ComputationContext& computation) {
        Polynomial<Rational> quot, rem;
        if (P.degree() >= Q.degree()) {
            strategy.poly_divide(P, Q, quot, rem);
        } else {
            quot = Polynomial<Rational>(Q.variable_name);
            rem = P;
        }

        std::vector<std::pair<Polynomial<Rational>, int>> factors;
        if (!strategy.factor_denominator(Q, factors)) {
            /// Q 上因式分解未决时保留未求值积分节点.
            return Integrator::depends_on(expr, var)
                ? LMCAS::detail::make_expression_ptr(
                      LMCAS::detail::make_node<IntegralNode>(
                          LMCAS::detail::node(expr), var))
                : nullptr;
        }
        if (factors.empty() && quot.is_zero()) {
            return SymbolicExpr::number(0);
        }

        std::vector<Polynomial<Rational>> numerators;
        if (!rem.is_zero()) {
            auto coefficients = strategy.solve_coefficients(
                rem, Q, factors, numerators, computation);
            if (!coefficients) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    coefficients.error());
            }
            if (!coefficients.value()) {
                return LMCAS::detail::make_expression_ptr(
                    LMCAS::detail::make_node<IntegralNode>(
                        LMCAS::detail::node(expr), var));
            }
        }
    auto result = rd_integrate_polynomial(quot, var);
    return rd_add_fraction_primitives(strategy, expr, var, factors, numerators, result);
}
}

Result<std::shared_ptr<SymbolicExpr>> RationalDecompositionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {
    (void)ctx;
    (void)depth;

    Polynomial<Rational> P, Q;
    auto extracted = extract_rational(expr, var, P, Q);
    if (!extracted) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(extracted.error());
    }
    if (!extracted.value()) { return nullptr; }
    rd_reduce_reciprocal(*this, P, Q);
    if (Q.is_zero() || Q.degree() < 1) { return nullptr; }

    /**
     * @brief 常数分子与一次分母直接积分，其余情形分解后逐项积分。
     * @note 不可约二次因子的原函数可含 arctan 与 ln。
     */
    if (Q.degree() == 1 && P.degree() <= 0 &&
        !P.coeffs.empty() && Q.coeffs.size() > 1 &&
        Q.coeffs[1] != Rational(0)) {
        auto coefficient = SymbolicExpr::number(
            P.coeffs[0] / Q.coeffs[1]);
        auto logarithm = SymbolicExpr::ln(poly_to_symbolic(Q));
        return SymbolicExpr::multiply(
            coefficient, logarithm)->simplify();
    }
    if (Q.degree() < 2) { return nullptr; }
    return rd_integrate_proper_fraction(*this, expr, var, P, Q, computation);
}

}
