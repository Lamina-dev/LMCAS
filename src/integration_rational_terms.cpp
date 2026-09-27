#include "internal/integration_rational_support.hpp"

namespace LMCAS {

namespace {
std::shared_ptr<SymbolicExpr> rd_integrate_linear_term(
    const Polynomial<Rational>& numerator, const Polynomial<Rational>& factor,
    int power, const std::string& var) {
        Rational r = Rational(0) - factor.coeffs[0];
        Rational A = (numerator.coeffs.size() > 0) ? numerator.coeffs[0] : Rational(0);
        if (A == Rational(0)) { return SymbolicExpr::number(0); }

        if (power == 1) {
            auto inner = rd_var_minus(var, r);
            return SymbolicExpr::multiply(SymbolicExpr::number(A), SymbolicExpr::ln(inner));
        }
        Rational coeff = Rational(0) - A / Rational(BigInt(power - 1));
        auto base = rd_var_minus(var, r);
        auto exp = SymbolicExpr::number(BigInt(power - 1));
        auto pw = SymbolicExpr::power(base, exp);
        auto inv_pw = SymbolicExpr::power(pw, SymbolicExpr::number(BigInt(-1)));
        return SymbolicExpr::multiply(SymbolicExpr::number(coeff), inv_pw);
}

std::shared_ptr<SymbolicExpr> rd_unevaluated_quadratic_term(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    int power, const std::string& var) {
    auto numerator_expr = poly_to_symbolic(numerator);
    auto denominator = poly_to_symbolic(factor);
    auto denominator_power = SymbolicExpr::power(
        denominator, SymbolicExpr::number(BigInt(power)));
    auto inverse = SymbolicExpr::power(
        denominator_power, SymbolicExpr::number(BigInt(-1)));
    auto integrand = SymbolicExpr::multiply(numerator_expr, inverse);
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(integrand), var));
}

std::shared_ptr<SymbolicExpr> rd_integrate_repeated_quadratic(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    const Rational& p, const Rational& B,
    const Rational& residual, const std::string& var) {
    Polynomial<Rational> linear(
        std::vector<Rational>{p / Rational(2), Rational(1)},
        factor.variable_name);
    Polynomial<Rational> linear_numerator(
        std::vector<Rational>{B}, numerator.variable_name);
    Polynomial<Rational> repeated_numerator(
        std::vector<Rational>{residual}, numerator.variable_name);
    auto logarithmic = rd_integrate_linear_term(
        linear_numerator, linear, 1, var);
    auto repeated = rd_integrate_linear_term(
        repeated_numerator, linear, 2, var);
    return SymbolicExpr::add(logarithmic, repeated);
}

std::shared_ptr<SymbolicExpr> rd_quadratic_logarithm(
    const Polynomial<Rational>& factor, const Rational& B) {
    const Rational half = B / Rational(2);
    auto result = SymbolicExpr::number(0);
    if (half == Rational(0)) {
        return result;
    }
    auto denominator = poly_to_symbolic(factor);
    auto logarithm = SymbolicExpr::multiply(
        SymbolicExpr::number(half), SymbolicExpr::ln(denominator));
    return SymbolicExpr::add(result, logarithm);
}

std::shared_ptr<SymbolicExpr> rd_shifted_double_variable(
    const Rational& p, const std::string& var) {
    auto two_x = SymbolicExpr::multiply(
        SymbolicExpr::number(BigInt(2)), SymbolicExpr::variable(var));
    return p == Rational(0)
        ? two_x
        : SymbolicExpr::add(two_x, SymbolicExpr::number(p));
}

std::shared_ptr<SymbolicExpr> rd_integrate_negative_discriminant(
    const std::shared_ptr<SymbolicExpr>& result,
    const std::shared_ptr<SymbolicExpr>& shifted,
    const Rational& residual, const Rational& discriminant) {
    auto square_root = SymbolicExpr::sqrt(
        SymbolicExpr::number(Rational(0) - discriminant));
    auto inverse_root = SymbolicExpr::power(
        square_root, SymbolicExpr::number(BigInt(-1)));
    auto arctangent = make_arctan(
        SymbolicExpr::multiply(shifted, inverse_root));
    auto coefficient = SymbolicExpr::multiply(
        SymbolicExpr::multiply(
            SymbolicExpr::number(residual),
            SymbolicExpr::number(BigInt(2))),
        inverse_root);
    return SymbolicExpr::add(
        result, SymbolicExpr::multiply(coefficient, arctangent));
}

std::shared_ptr<SymbolicExpr> rd_integrate_positive_discriminant(
    const std::shared_ptr<SymbolicExpr>& result,
    const std::shared_ptr<SymbolicExpr>& shifted,
    const Rational& residual, const Rational& discriminant) {
    auto square_root =
        SymbolicExpr::sqrt(SymbolicExpr::number(discriminant));
    auto negative_root = SymbolicExpr::multiply(
        SymbolicExpr::number(BigInt(-1)), square_root);
    auto ratio = SymbolicExpr::divide(
        SymbolicExpr::add(shifted, negative_root),
        SymbolicExpr::add(shifted, square_root));
    auto inverse_root = SymbolicExpr::power(
        square_root, SymbolicExpr::number(BigInt(-1)));
    auto coefficient = SymbolicExpr::multiply(
        SymbolicExpr::number(residual), inverse_root);
    return SymbolicExpr::add(
        result,
        SymbolicExpr::multiply(coefficient, SymbolicExpr::ln(ratio)));
}

std::shared_ptr<SymbolicExpr> rd_integrate_quadratic_term(
    const Polynomial<Rational>& numerator, const Polynomial<Rational>& factor,
    int power, const std::string& var) {
    if (power >= 2) {
        return rd_unevaluated_quadratic_term(
            numerator, factor, power, var);
    }

    const Rational q = factor.coeffs.empty()
        ? Rational(0) : factor.coeffs[0];
    const Rational p = factor.coeffs.size() > 1
        ? factor.coeffs[1] : Rational(0);
    const Rational discriminant = p * p - Rational(4) * q;
    const Rational B = numerator.coeffs.size() > 1
        ? numerator.coeffs[1] : Rational(0);
    const Rational C = numerator.coeffs.empty()
        ? Rational(0) : numerator.coeffs[0];
    const Rational residual = C - (B * p) / Rational(2);

    if (discriminant == Rational(0)) {
        return rd_integrate_repeated_quadratic(
            numerator, factor, p, B, residual, var);
    }

    auto result = rd_quadratic_logarithm(factor, B);
    if (residual == Rational(0)) {
        return result;
    }
    auto shifted = rd_shifted_double_variable(p, var);
    if (discriminant.get_numerator().is_negative()) {
        return rd_integrate_negative_discriminant(
            result, shifted, residual, discriminant);
    }
    return rd_integrate_positive_discriminant(
        result, shifted, residual, discriminant);
}
}


std::shared_ptr<SymbolicExpr> RationalDecompositionStrategy::integrate_term(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    int power, const std::string& var) {

    if (numerator.is_zero()) { return SymbolicExpr::number(0); }

    if (factor.degree() == 1) { return rd_integrate_linear_term(numerator, factor, power, var); }
    if (factor.degree() == 2) { return rd_integrate_quadratic_term(numerator, factor, power, var); }
    return SymbolicExpr::number(0);
}

}
