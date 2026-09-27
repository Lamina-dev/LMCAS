/**
 * @file poly_utils.hpp
 * @brief 多项式转换与精确符号多项式代数。
 */
#pragma once

#include "polynomial_conversion.hpp"

namespace LMCAS {

using SymbolicGcdResult = Result<std::shared_ptr<SymbolicExpr>>;

/**
 * Compute the monic GCD of two exact rational multivariate polynomials.
 *
 * The expressions may contain any number of symbolic variables, exact integer
 * or rational coefficients, addition, multiplication, and nonnegative integer
 * powers. Approximate numbers and non-polynomial nodes are rejected. A
 * successful result is verified to divide both inputs exactly.
 */
LMCAS_API SymbolicGcdResult symbolic_polynomial_gcd(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    ComputationContext& context);

/** Return the positive coefficient content of an exact rational polynomial. */
LMCAS_API Result<Rational> symbolic_polynomial_content(
    const SymbolicExpr& expression,
    ComputationContext& context);

} // namespace LMCAS
