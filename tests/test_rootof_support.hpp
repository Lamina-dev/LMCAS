#pragma once
#include "test_common.hpp"
#include "root_of_utils.hpp"
#include "solve_polynomial.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace LMCAS;

using LMCAS::Polynomial;

inline std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

inline double eval_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) { return std::nan(""); }
    return expr->to_numeric();
}
