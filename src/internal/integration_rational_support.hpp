#pragma once
#include "internal/integration_support.hpp"

namespace LMCAS::detail {

Result<bool> extract_rational(
    const SymbolicExpr& expr, const std::string& var,
    Polynomial<Rational>& numerator, Polynomial<Rational>& denominator);

bool factor_denominator(
    const Polynomial<Rational>& denominator,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors);

Result<bool> solve_coefficients(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& denominator,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    std::vector<Polynomial<Rational>>& numerators,
    ComputationContext& context);

IntegralParts integrate_rational_term(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    int power, const std::string& var);

}
