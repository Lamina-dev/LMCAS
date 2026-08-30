/**
 * @file newton_raphson.hpp
 * @brief 数值求根：Newton-Raphson 迭代、二分法、Sturm 实根隔离。
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include "solve_strategies.hpp"
#include "computation_context.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <utility>

namespace lamina {

using NumericRootResult = Result<std::optional<NumericRoot>>;
using NumericRootsResult = Result<std::vector<NumericRoot>>;

/**
 * @brief Numerically solve in the supported univariate domain.
 *
 * Polynomial inputs use Sturm isolation followed by checked refinement.
 * Other inputs require a Newton initial value; zero is used when none is
 * supplied. Evaluation, cancellation, and resource failures are
 * returned to the caller.
 */
LAMINA_API NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts = {});

/**
 * @brief Numerically solves using an isolated default computation context.
 */
LAMINA_API NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

/**
 * @brief Isolates all distinct real roots of an exact rational polynomial.
 */
LAMINA_API Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(
    const Polynomial<Rational>& poly,
    ComputationContext& context);

LAMINA_API Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(const Polynomial<Rational>& poly);

/**
 * @brief 带错误传播和资源预算的 Newton-Raphson 迭代。
 */
LAMINA_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    ComputationContext& context,
    const SolveOptions& opts);

LAMINA_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts);

/**
 * @brief 带区间、错误传播和资源预算的 Newton-Raphson 迭代。
 */
LAMINA_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    ComputationContext& context,
    const SolveOptions& opts);

LAMINA_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts);

/**
 * @brief 带错误传播和资源预算的二分法求根。
 */
LAMINA_API NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    ComputationContext& context,
    const SolveOptions& opts);

LAMINA_API NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts);

}
