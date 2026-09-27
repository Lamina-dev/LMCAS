#pragma once
#include "test_common.hpp"
#include "newton_raphson.hpp"
#include "solve_polynomial.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <set>
#include <sstream>

using namespace LMCAS;

inline std::shared_ptr<SymbolicExpr> num_expr(int n) { return SymbolicExpr::number(n); }

inline LMCAS::Polynomial<Rational> poly_from_roots(const std::vector<int>& roots) {
    LMCAS::Polynomial<Rational> result({Rational(1)}, "x");
    for (int r : roots) {

        LMCAS::Polynomial<Rational> factor({Rational(-r), Rational(1)}, "x");
        result = result * factor;
    }
    return result;
}

inline double eval_poly_at_double(const LMCAS::Polynomial<Rational>& poly, double x) {
    double result = 0.0;
    double x_pow = 1.0;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        result += poly.coeffs[i].to_double() * x_pow;
        x_pow *= x;
    }
    return result;
}
