/**
 * @file solve_transcendental.hpp
 * @brief 超越方程求解：三角、指数、对数方程的符号反演。
 */
#pragma once

#include "symbolic.hpp"
#include "solve_strategies.hpp"
#include "result.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace LMCAS {

/** @brief 换元结果，记录换元表达式、换元后的多项式及换元变量名。 */
struct SubstitutionResult {
    std::shared_ptr<SymbolicExpr> u_expr;      ///< 换元表达式（u = f(x) 中的 f(x)）
    std::shared_ptr<SymbolicExpr> poly_in_u;   ///< 换元后关于 u 的多项式
    std::string u_var;                         ///< 换元变量名
};

/**
 * @brief 求解三角、指数、对数类超越方程。
 * @param expr 待求解的表达式（视为等于零）。
 * @param var 求解变量名。
 * @return 完整实数解集。三角函数解族使用精确符号 pi 和新引入的整数参数；
 * 未求出的原像保留原始谓词，超出处理范围的变换返回 Inconclusive 而非部分根。
 */
LMCAS_API SolveResult solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var,
    ComputationContext& context, const SolveOptions& options = {});

LMCAS_API SolveResult solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var,
    const SolveOptions& options = {});

/**
 * @brief 检测表达式中可用的换元模式。
 * @param expr 待分析的表达式
 * @param var 目标变量名
 * @return 成功时返回有效换元或空 optional；资源等计算错误保留原错误分类与来源。
 */
LMCAS_API Result<std::optional<SubstitutionResult>> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

}
