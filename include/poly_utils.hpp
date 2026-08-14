/**
 * @file poly_utils.hpp
 * @brief Symbolic expression and univariate polynomial conversion utilities.
 */
#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "symbolic.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace lamina {

/**
 * @brief Convert a symbolic expression to a univariate polynomial.
 *
 * The stable API supports BigInt, Rational, and SymbolicPolyCoeff
 * coefficients. Expressions outside the supported conversion domain produce the
 * zero polynomial; callers that require proof must use
 * recognize_rational_polynomial().
 */
template <typename T>
LAMINA_API Polynomial<T> symbolic_to_poly(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);

/**
 * @brief Convert a supported univariate polynomial to a symbolic expression.
 *
 * The stable API supports BigInt, Rational, and SymbolicPolyCoeff
 * coefficients.
 */
template <typename T>
LAMINA_API std::shared_ptr<SymbolicExpr> poly_to_symbolic(
    const Polynomial<T>& polynomial);

using OptionalRationalPolynomial = std::optional<Polynomial<Rational>>;

using SymbolicGcdResult = Result<std::shared_ptr<SymbolicExpr>>;

/**
 * Compute the monic GCD of two exact rational multivariate polynomials.
 *
 * The expressions may contain any number of symbolic variables, exact integer
 * or rational coefficients, addition, multiplication, and nonnegative integer
 * powers. Approximate numbers and non-polynomial nodes are rejected. A
 * successful result is verified to divide both inputs exactly.
 */
LAMINA_API SymbolicGcdResult symbolic_polynomial_gcd(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    ComputationContext& context);

/** Return the positive coefficient content of an exact rational polynomial. */
LAMINA_API Result<Rational> symbolic_polynomial_content(
    const SymbolicExpr& expression,
    ComputationContext& context);

/**
 * @brief Prove that an expression is an exact univariate rational polynomial.
 *
 * A successful empty optional means that the expression is outside this
 * structural support domain. Approximate numbers are never converted to exact
 * coefficients. Traversal and expansion consume the supplied context budget.
 */
LAMINA_API Result<OptionalRationalPolynomial> recognize_rational_polynomial(
    const SymbolicExpr& expression,
    const std::string& variable,
    ComputationContext& context);

/** @brief Symbolic expression coefficient for symbolic polynomial algorithms. */
struct SymbolicPolyCoeff {
    std::shared_ptr<SymbolicExpr> val;

    SymbolicPolyCoeff() : val(SymbolicExpr::number(0)) {}
    explicit SymbolicPolyCoeff(int value) : val(SymbolicExpr::number(value)) {}
    SymbolicPolyCoeff(std::shared_ptr<SymbolicExpr> value) : val(std::move(value)) {}

    bool operator==(const SymbolicPolyCoeff& other) const {
        if (!val || !other.val) return false;
        if (val == other.val) return true;
        auto difference = SymbolicExpr::add(
            val,
            SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)));
        return difference->simplify()->is_zero();
    }

    bool operator!=(const SymbolicPolyCoeff& other) const {
        return !(*this == other);
    }

    SymbolicPolyCoeff operator+(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(val, other.val));
    }

    SymbolicPolyCoeff operator-(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(
            val,
            SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1))));
    }

    SymbolicPolyCoeff operator*(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, other.val));
    }

    SymbolicPolyCoeff operator/(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::divide(val, other.val));
    }

    SymbolicPolyCoeff operator-() const {
        return SymbolicPolyCoeff(
            SymbolicExpr::multiply(val, SymbolicExpr::number(-1)));
    }

    std::string ToString() const {
        return val ? val->to_string() : "0";
    }

    friend SymbolicPolyCoeff abs(const SymbolicPolyCoeff& value) {
        return value;
    }

    bool operator<(const SymbolicPolyCoeff&) const {
        return false;
    }
};

/** @brief Extract a supported coefficient from a symbolic expression. */
template <typename T>
T extract_coeff_value(const std::shared_ptr<SymbolicExpr>& coefficient);

template <>
inline SymbolicPolyCoeff extract_coeff_value<SymbolicPolyCoeff>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    return SymbolicPolyCoeff(coefficient);
}

template <>
LAMINA_API BigInt extract_coeff_value<BigInt>(
    const std::shared_ptr<SymbolicExpr>& coefficient);

template <>
LAMINA_API Rational extract_coeff_value<Rational>(
    const std::shared_ptr<SymbolicExpr>& coefficient);

/** @brief Return whether an expression contains a free occurrence of variable. */
LAMINA_API bool contains(
    const SymbolicExpr& expression,
    const std::string& variable);

/**
 * @brief Perform Gaussian elimination on an augmented symbolic matrix.
 */
LAMINA_API void gaussian_eliminate(
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& matrix,
    std::size_t rows,
    std::size_t columns,
    std::vector<std::size_t>& pivot_column_for_row,
    int& sign);

} // namespace lamina
