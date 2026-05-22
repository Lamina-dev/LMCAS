#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include "solve_strategies.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <utility>

namespace lamina {

// 数值求解单变量方程
LAMINA_API std::vector<NumericRoot> solve_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

// 对多项式进行实根隔离 (Sturm 序列)
LAMINA_API std::vector<std::pair<Rational, Rational>> isolate_real_roots(
    const Polynomial<Rational>& poly);

// 单次 Newton-Raphson 求根（给定初始点）
LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts);

// Newton-Raphson with bracket for bisection fallback
// When the derivative is near zero, falls back to bisection on [lo, hi]
LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts);

// Bisection method on a bracket [lo, hi]
// Used as fallback when Newton-Raphson encounters near-zero derivative
LAMINA_API std::optional<NumericRoot> bisection(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts);

} // namespace lamina
