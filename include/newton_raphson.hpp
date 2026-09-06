/**
 * @file newton_raphson.hpp
 * @brief 数值求根:Newton-Raphson 迭代,二分法,Sturm 实根隔离.
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

namespace LMCAS {

using NumericRootResult = Result<std::optional<NumericRoot>>;
using NumericRootsResult = Result<std::vector<NumericRoot>>;

/**
 * @brief 在受支持的一元表达式域内执行数值求根.
 *
 * 多项式输入先由 Sturm 序列隔离实根,再执行受检精化；隔离区间通过
 * 溢出安全的中点生成 Newton 初值，包括最大有限浮点端点。
 * 其他输入使用 Newton 初值,缺省初值为零；求值、取消与资源诊断
 * 均通过 Result 通道返回调用方。
 *
 * @see Jacques Charles François Sturm,
 *      "Mémoire sur la résolution des équations numériques," 1829.
 */
LMCAS_API NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts = {});

/**
 * @brief Numerically solves using an isolated default computation context.
 */
LMCAS_API NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

/**
 * @brief Isolates all distinct real roots of an exact rational polynomial.
 */
LMCAS_API Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(
    const Polynomial<Rational>& poly,
    ComputationContext& context);

LMCAS_API Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(const Polynomial<Rational>& poly);

/**
 * @brief 带错误传播和资源预算的 Newton-Raphson 迭代.
 */
LMCAS_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    ComputationContext& context,
    const SolveOptions& opts);

LMCAS_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts);

/**
 * @brief 带区间,错误传播和资源预算的 Newton-Raphson 迭代.
 */
LMCAS_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    ComputationContext& context,
    const SolveOptions& opts);

LMCAS_API NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts);

/**
 * @brief 带错误传播和资源预算的二分法求根.
 */
LMCAS_API NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    ComputationContext& context,
    const SolveOptions& opts);

LMCAS_API NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts);

}
