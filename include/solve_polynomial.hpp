/**
 * @file solve_polynomial.hpp
 * @brief 多项式方程闭式求解：三次、四次公式，有理根检验，无平方因子分解。
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "rational.hpp"
#include <vector>
#include <utility>
#include <memory>
#include <string>

namespace LMCAS {

/**
 * @brief 使用三次公式求解三次方程 ax³ + bx² + cx + d = 0。
 *
 * 当高阶系数为零时依次退化为二次或一次方程；非零常数方程返回空根集。
 * @throws std::invalid_argument 方程恒等于零或退化后不再依赖 @p var，
 *         因而无法用有限根列表表示。
 * @param a 三次项系数
 * @param b 二次项系数
 * @param c 一次项系数
 * @param d 常数项
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string& var);

/**
 * @brief 使用四次公式求解四次方程 ax⁴ + bx³ + cx² + dx + e = 0。
 *
 * 纯数值双二次分支使用避免相消的二次求根公式，以同时保留尺度悬殊的根。
 * @param a 四次项系数
 * @param b 三次项系数
 * @param c 二次项系数
 * @param d 一次项系数
 * @param e 常数项
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var);

/**
 * @brief 求解双二次方程 ax⁴ + bx² + c = 0。
 *
 * 纯数值系数使用避免相消的二次求根公式求解 @f$u=x^2@f$；
 * 判别式或 @f$u@f$ 为负时保留对应的复根，而不按固定绝对容差截断。
 * 当 @p a 为零时退化为 @f$bx^2+c=0@f$；非零常数方程返回空根集。
 * @throws std::invalid_argument 方程恒等于零或退化后不再依赖 @p var，
 *         因而无法用有限根列表表示。
 * @param a 四次项系数
 * @param b 二次项系数
 * @param c 常数项
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string& var);

/**
 * @brief 查找有理系数多项式的所有有理根。
 * @param poly 有理系数多项式
 * @return 有理根列表，保留重数；零多项式和常数多项式返回空列表。
 * 公共非零有理系数先通过首一化消去；一次余式直接精确求解。
 */
LMCAS_API std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly);

/**
 * @brief 对有理系数多项式进行无平方因子分解。
 * @param poly 有理系数多项式
 * @return 因子与重数的列表
 */
LMCAS_API std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly);

/**
 * @brief 通过因式分解求解多项式方程。
 * @param poly 符号系数多项式
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
