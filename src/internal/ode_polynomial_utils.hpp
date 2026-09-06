#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace LMCAS::ode_root_detail {

using RealPolynomial = std::vector<double>;

struct CharRoot {
    double real_part;
    double imag_part;
    int multiplicity;
    bool is_complex;
};

struct Complex {
    double re;
    double im;

    Complex(double real = 0.0, double imaginary = 0.0)
        : re(real), im(imaginary) {}
    Complex operator*(const Complex& other) const {
        return {re * other.re - im * other.im,
                re * other.im + im * other.re};
    }
    Complex operator+(const Complex& other) const {
        return {re + other.re, im + other.im};
    }
    Complex operator-(const Complex& other) const {
        return {re - other.re, im - other.im};
    }
    Complex operator/(const Complex& other) const {
        const double denominator =
            other.re * other.re + other.im * other.im;
        if (denominator < 1e-300)
            return {INFINITY, INFINITY};
        return {(re * other.re + im * other.im) / denominator,
                (im * other.re - re * other.im) / denominator};
    }
    double magnitude() const {
        return std::hypot(re, im);
    }
};

inline double polynomial_root_backward_error(
    const RealPolynomial& coefficients, const Complex& root) {
    Complex value{coefficients.front(), 0.0};
    double scale = std::abs(coefficients.front());
    const double root_magnitude = root.magnitude();
    for (size_t i = 1; i < coefficients.size(); ++i) {
        value = value * root + Complex{coefficients[i], 0.0};
        scale = scale * root_magnitude + std::abs(coefficients[i]);
    }
    if (!std::isfinite(value.re) || !std::isfinite(value.im) ||
        !std::isfinite(scale))
        return INFINITY;
    if (scale == 0.0) {
        // At an exact zero root with zero constant term, both the residual
        // and its coefficientwise bound vanish: no perturbation is needed.
        // A nonzero root can instead underflow the bound, which is not proof.
        return root.re == 0.0 && root.im == 0.0 &&
                       coefficients.back() == 0.0
                   ? 0.0
                   : INFINITY;
    }
    return value.magnitude() / scale;
}

inline double polynomial_scale(const RealPolynomial& polynomial) {
    double scale = 0.0;
    for (double coefficient : polynomial) {
        scale = std::max(scale, std::abs(coefficient));
    }
    return scale;
}

inline void trim_polynomial(RealPolynomial& polynomial) {
    const double scale = polynomial_scale(polynomial);
    const double tolerance = std::max(1.0, scale) * 1e-10;
    auto first = polynomial.begin();
    while (first + 1 != polynomial.end() && std::abs(*first) <= tolerance) {
        ++first;
    }
    polynomial.erase(polynomial.begin(), first);
    if (polynomial.empty()) polynomial.push_back(0.0);
}

inline std::pair<RealPolynomial, RealPolynomial> divide_polynomials(
    const RealPolynomial& dividend, const RealPolynomial& divisor) {
    if (divisor.empty() || divisor.front() == 0.0 ||
        dividend.size() < divisor.size()) {
        return {{0.0}, dividend};
    }

    RealPolynomial remainder = dividend;
    RealPolynomial quotient(
        dividend.size() - divisor.size() + 1, 0.0);
    for (size_t i = 0; i < quotient.size(); ++i) {
        const double factor = remainder[i] / divisor.front();
        quotient[i] = factor;
        for (size_t j = 0; j < divisor.size(); ++j) {
            remainder[i + j] -= factor * divisor[j];
        }
    }
    remainder.erase(
        remainder.begin(), remainder.begin() + quotient.size());
    if (remainder.empty()) remainder.push_back(0.0);
    trim_polynomial(quotient);
    trim_polynomial(remainder);
    return {std::move(quotient), std::move(remainder)};
}

inline RealPolynomial monic_polynomial(RealPolynomial polynomial) {
    trim_polynomial(polynomial);
    if (polynomial.empty() || polynomial.front() == 0.0) return {0.0};
    const double leading = polynomial.front();
    for (double& coefficient : polynomial) coefficient /= leading;
    return polynomial;
}

inline RealPolynomial polynomial_derivative(
    const RealPolynomial& polynomial) {
    if (polynomial.size() <= 1) return {0.0};
    RealPolynomial derivative;
    derivative.reserve(polynomial.size() - 1);
    const size_t degree = polynomial.size() - 1;
    for (size_t i = 0; i < degree; ++i) {
        derivative.push_back(
            polynomial[i] * static_cast<double>(degree - i));
    }
    return derivative;
}

inline RealPolynomial approximate_polynomial_gcd(
    RealPolynomial lhs, RealPolynomial rhs) {
    lhs = monic_polynomial(std::move(lhs));
    rhs = monic_polynomial(std::move(rhs));
    for (int iteration = 0; iteration < 8 && rhs.size() > 1; ++iteration) {
        auto division = divide_polynomials(lhs, rhs);
        RealPolynomial remainder = std::move(division.second);
        if (remainder.size() == 1 && remainder.front() == 0.0) {
            return rhs;
        }
        lhs = std::move(rhs);
        rhs = monic_polynomial(std::move(remainder));
    }
    return {1.0};
}

inline bool polynomial_remainder_is_small(
    const RealPolynomial& remainder, const RealPolynomial& dividend) {
    return polynomial_scale(remainder) <=
           std::max(1.0, polynomial_scale(dividend)) * 1e-7;
}

} // namespace LMCAS::ode_root_detail
