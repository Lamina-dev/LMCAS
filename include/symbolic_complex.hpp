/**
 * @file symbolic_complex.hpp
 * @brief 符号复数运算：加减乘除、共轭、模、辐角、极坐标形式。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lamina {

/** @brief 符号复数，由实部和虚部组成 */
struct ComplexSymbolic {
    std::shared_ptr<SymbolicExpr> real;  ///< 实部
    std::shared_ptr<SymbolicExpr> imag;  ///< 虚部
};

using ComplexSymbolicResult = Result<ComplexSymbolic>;
using ComplexRootsResult = Result<std::vector<ComplexSymbolic>>;

/**
 * @brief 构造符号复数
 * @param real 实部表达式
 * @param imag 虚部表达式
 * @return 符号复数对象
 */
inline ComplexSymbolic make_complex(
    std::shared_ptr<SymbolicExpr> real,
    std::shared_ptr<SymbolicExpr> imag
) {
    return ComplexSymbolic{real, imag};
}

/**
 * @brief 复数加法
 * @param a 加数
 * @param b 加数
 * @return 和
 */

LAMINA_API ComplexSymbolicResult complex_add_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_add_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b);

/**
 * @brief 复数减法
 * @param a 被减数
 * @param b 减数
 * @return 差
 */

LAMINA_API ComplexSymbolicResult complex_sub_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_sub_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b);

/**
 * @brief 复数乘法
 * @param a 乘数
 * @param b 乘数
 * @return 积
 */

LAMINA_API ComplexSymbolicResult complex_mul_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_mul_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b);

/**
 * @brief 复数除法
 * @param a 被除数
 * @param b 除数
 * @return 商
 */

LAMINA_API ComplexSymbolicResult complex_div_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_div_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b);

/**
 * @brief 计算复数的共轭
 * @param z 输入复数
 * @return 共轭复数
 */

LAMINA_API ComplexSymbolicResult complex_conj_checked(
    const ComplexSymbolic& z,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_conj_checked(const ComplexSymbolic& z);

/**
 * @brief 计算复数的模
 * @param z 输入复数
 * @return 模的符号表达式
 */

LAMINA_API ExpressionResult complex_abs_checked(
    const ComplexSymbolic& z,
    ComputationContext& context);

LAMINA_API ExpressionResult complex_abs_checked(const ComplexSymbolic& z);

/**
 * @brief 计算复数的辐角
 * @param z 输入复数
 * @return 辐角的符号表达式
 */

LAMINA_API ExpressionResult complex_arg_checked(
    const ComplexSymbolic& z,
    ComputationContext& context);

LAMINA_API ExpressionResult complex_arg_checked(const ComplexSymbolic& z);

/**
 * @brief 将极坐标形式转换为指数形式复数 r*e^(i*theta)
 * @param r 模
 * @param theta 辐角
 * @return 对应的符号复数
 */

LAMINA_API ComplexSymbolicResult complex_exp_form_checked(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_exp_form_checked(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta);

/**
 * @brief 将极坐标形式转换为三角形式复数 r*(cos(theta) + i*sin(theta))
 * @param r 模
 * @param theta 辐角
 * @return 对应的符号复数
 */

LAMINA_API ComplexSymbolicResult complex_trig_form_checked(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta,
    ComputationContext& context);

LAMINA_API ComplexSymbolicResult complex_trig_form_checked(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta);

/**
 * @brief 求复数 c 的 n 次方根
 * @param c 被开方的复数表达式
 * @param n 根的次数
 * @return n 个复数根的列表
 */

LAMINA_API ComplexRootsResult solve_complex_nth_root_checked(
    std::shared_ptr<SymbolicExpr> c,
    int n,
    ComputationContext& context);

LAMINA_API ComplexRootsResult solve_complex_nth_root_checked(
    std::shared_ptr<SymbolicExpr> c,
    int n);

/**
 * @brief 求复系数一元二次方程 az^2 + bz + c = 0 的根.
 */
LAMINA_API ComplexRootsResult solve_complex_quadratic_checked(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c,
    ComputationContext& context);

LAMINA_API ComplexRootsResult solve_complex_quadratic_checked(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c);

LAMINA_API ExpressionResult complex_locus_circle_checked(
    const ComplexSymbolic& a,
    std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var,
    ComputationContext& context);

LAMINA_API ExpressionResult complex_locus_circle_checked(
    const ComplexSymbolic& a,
    std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var = "z");

LAMINA_API ExpressionResult complex_locus_perpendicular_bisector_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    const std::string& z_var,
    ComputationContext& context);

LAMINA_API ExpressionResult complex_locus_perpendicular_bisector_checked(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    const std::string& z_var = "z");

}
