/**
 * @file numerical_integration.hpp
 * @brief 数值积分算法：辛普森法则和高斯求积（符号-数值混合桥接）。
 */
#pragma once
#include "symbolic_ast.hpp"
#include <memory>
#include <string>

class SymbolicExpr;


namespace lamina {

/**
 * @brief 辛普森法则数值积分（返回符号常数或未求值表达式）
 * @param f 被积函数
 * @param var 积分变量
 * @param a 积分下限
 * @param b 积分上限
 * @param n 区间等分数（必须为偶数）
 * @return 积分近似值的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> quadrature_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

/**
 * @brief 高斯-勒让德求积法则数值积分
 * @param f 被积函数
 * @param var 积分变量
 * @param a 积分下限
 * @param b 积分上限
 * @param n 积分点数（通常为 5, 7, 或 10）
 * @return 积分近似值的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> quadrature_gaussian(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 5
);

/**
 * @brief 自适应辛普森积分：递归细分直至误差估计低于容差。
 * @param f 被积函数
 * @param var 积分变量
 * @param a 下限（数值）
 * @param b 上限（数值）
 * @param tol 误差容差
 * @return 积分近似值（NumberNode）
 */
LAMINA_API std::shared_ptr<SymbolicExpr> adaptive_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    double tol = 1e-10
);

/**
 * @brief 通用数值积分入口（默认复合辛普森）。
 * @param f 被积函数
 * @param var 积分变量
 * @param a 下限
 * @param b 上限
 * @param n 子区间数（偶数）
 * @return 积分近似值的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> numerical_integrate(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

} // namespace lamina
