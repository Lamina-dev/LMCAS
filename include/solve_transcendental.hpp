/**
 * @file solve_transcendental.hpp
 * @brief 超越方程求解：三角、指数、对数方程的符号反演。
 */
#pragma once

#include "symbolic.hpp"
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
 * @brief 求解超越方程（三角、指数、对数类型）。
 * @param expr 待求解的表达式（视为等于零）
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

/**
 * @brief 检测表达式中可用的换元模式。
 * @param expr 待分析的表达式
 * @param var 目标变量名
 * @return 若检测到有效换元则返回 SubstitutionResult，否则返回 nullopt
 */
LMCAS_API std::optional<SubstitutionResult> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

}
