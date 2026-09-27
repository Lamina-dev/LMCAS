#pragma once

#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include <cmath>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace LMCAS;

inline std::shared_ptr<SymbolicExpr> linear(int a, int b) {
    auto x = SymbolicExpr::variable("x");
    auto ax = SymbolicExpr::multiply(SymbolicExpr::number(a), x);
    return SymbolicExpr::add(ax, SymbolicExpr::number(b));
}

inline std::shared_ptr<SymbolicExpr> negate(std::shared_ptr<SymbolicExpr> e) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), e);
}

inline constexpr InequalityType inequality_types[] = {
    InequalityType::GreaterThan,
    InequalityType::GreaterEqual,
    InequalityType::LessThan,
    InequalityType::LessEqual
};

inline std::shared_ptr<SymbolicExpr> build_poly(const std::vector<int>& coeffs) {
    auto x = SymbolicExpr::variable("x");
    auto result = SymbolicExpr::number(coeffs[0]);
    for (size_t i = 1; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) continue;
        auto term = SymbolicExpr::multiply(
            SymbolicExpr::number(coeffs[i]),
            SymbolicExpr::power(x, SymbolicExpr::number(static_cast<int>(i))));
        result = SymbolicExpr::add(result, term);
    }
    return result;
}

inline double eval_poly(const std::shared_ptr<SymbolicExpr>& poly, double point) {
    auto substituted = poly->substitute("x", SymbolicExpr::number(point));
    return substituted->to_numeric();
}

inline bool satisfies(double value, InequalityType type) {
    switch (type) {
        case InequalityType::GreaterThan: return value > 0;
        case InequalityType::GreaterEqual: return value >= 0;
        case InequalityType::LessThan: return value < 0;
        case InequalityType::LessEqual: return value <= 0;
    }
    return false;
}

inline std::vector<int> generate_coefficients(std::mt19937& rng) {
    std::uniform_int_distribution<int> degree_dist(2, 2);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);
    int degree = degree_dist(rng);
    std::vector<int> coeffs(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        coeffs[i] = coeff_dist(rng);
    }
    while (coeffs[degree] == 0) {
        coeffs[degree] = coeff_dist(rng);
    }
    return coeffs;
}

inline std::string format_coefficients(const std::vector<int>& coeffs) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (i > 0) out << ",";
        out << coeffs[i];
    }
    out << "]";
    return out.str();
}
