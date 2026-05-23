/**
 * @file newton_raphson.hpp
 * @brief 数值求根：Newton-Raphson 迭代、二分法、Sturm 实根隔离。
 */
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

/**
 * @brief 对表达式进行数值求根，综合使用 Newton-Raphson 和根隔离方法。
 * @param expr 待求根的表达式（视为等于零）
 * @param var 求解变量名
 * @param opts 求解选项（容差、最大迭代次数等）
 * @return 数值根列表
 */
LAMINA_API std::vector<NumericRoot> solve_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

/**
 * @brief 使用 Sturm 序列隔离有理系数多项式的所有实根。
 * @param poly 有理系数多项式
 * @return 隔离区间列表，每个区间包含恰好一个实根
 */
LAMINA_API std::vector<std::pair<Rational, Rational>> isolate_real_roots(
    const Polynomial<Rational>& poly);

/**
 * @brief Newton-Raphson 迭代求根（无区间约束）。
 * @param f 目标函数表达式
 * @param df 目标函数的导数表达式
 * @param var 求解变量名
 * @param x0 迭代初始值
 * @param opts 求解选项
 * @return 收敛时返回数值根，否则返回 nullopt
 */
LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts);

/**
 * @brief Newton-Raphson 迭代求根（带区间约束，迭代超出区间时回退到二分法）。
 * @param f 目标函数表达式
 * @param df 目标函数的导数表达式
 * @param var 求解变量名
 * @param x0 迭代初始值
 * @param bracket_lo 搜索区间下界
 * @param bracket_hi 搜索区间上界
 * @param opts 求解选项
 * @return 收敛时返回数值根，否则返回 nullopt
 */
LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts);

/**
 * @brief 二分法求根。
 * @param f 目标函数表达式
 * @param var 求解变量名
 * @param lo 区间下界（要求 f(lo) 与 f(hi) 异号）
 * @param hi 区间上界
 * @param opts 求解选项
 * @return 收敛时返回数值根，否则返回 nullopt
 */
LAMINA_API std::optional<NumericRoot> bisection(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts);

}
