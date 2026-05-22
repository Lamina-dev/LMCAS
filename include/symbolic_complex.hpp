/**
 * @file symbolic_complex.hpp
 * @brief 符号复数运算：加减乘除、共轭、模、辐角、极坐标形式。
 */
#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

/** @brief 符号复数，由实部和虚部组成 */
struct ComplexSymbolic {
    std::shared_ptr<SymbolicExpr> real;  ///< 实部
    std::shared_ptr<SymbolicExpr> imag;  ///< 虚部
};

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
ComplexSymbolic complex_add(const ComplexSymbolic& a, const ComplexSymbolic& b);

/**
 * @brief 复数减法
 * @param a 被减数
 * @param b 减数
 * @return 差
 */
ComplexSymbolic complex_sub(const ComplexSymbolic& a, const ComplexSymbolic& b);

/**
 * @brief 复数乘法
 * @param a 乘数
 * @param b 乘数
 * @return 积
 */
ComplexSymbolic complex_mul(const ComplexSymbolic& a, const ComplexSymbolic& b);

/**
 * @brief 复数除法
 * @param a 被除数
 * @param b 除数
 * @return 商
 */
ComplexSymbolic complex_div(const ComplexSymbolic& a, const ComplexSymbolic& b);

/**
 * @brief 计算复数的共轭
 * @param z 输入复数
 * @return 共轭复数
 */
ComplexSymbolic complex_conj(const ComplexSymbolic& z);

/**
 * @brief 计算复数的模
 * @param z 输入复数
 * @return 模的符号表达式
 */
std::shared_ptr<SymbolicExpr> complex_abs(const ComplexSymbolic& z);

/**
 * @brief 计算复数的辐角
 * @param z 输入复数
 * @return 辐角的符号表达式
 */
std::shared_ptr<SymbolicExpr> complex_arg(const ComplexSymbolic& z);

/**
 * @brief 将极坐标形式转换为指数形式复数 r*e^(i*theta)
 * @param r 模
 * @param theta 辐角
 * @return 对应的符号复数
 */
ComplexSymbolic complex_exp_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

/**
 * @brief 将极坐标形式转换为三角形式复数 r*(cos(theta) + i*sin(theta))
 * @param r 模
 * @param theta 辐角
 * @return 对应的符号复数
 */
ComplexSymbolic complex_trig_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

/**
 * @brief 求复数 c 的 n 次方根
 * @param c 被开方的复数表达式
 * @param n 根的次数
 * @return n 个复数根的列表
 */
std::vector<ComplexSymbolic> solve_complex_nth_root(
    std::shared_ptr<SymbolicExpr> c,
    int n
);

/**
 * @brief 求复系数一元二次方程 az^2 + bz + c = 0 的根
 * @param a 二次项系数
 * @param b 一次项系数
 * @param c 常数项
 * @return 两个复数根
 */
std::vector<ComplexSymbolic> solve_complex_quadratic(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c
);

/**
 * @brief 生成以 a 为圆心、r 为半径的复数轨迹方程 |z - a| = r
 * @param a 圆心复数
 * @param r 半径
 * @param z_var 复变量名，默认为 "z"
 * @return 轨迹方程的符号表达式
 */
std::shared_ptr<SymbolicExpr> complex_locus_circle(
    const ComplexSymbolic& a,
    std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var = "z"
);

/**
 * @brief 生成复数 a 和 b 的中垂线轨迹方程 |z - a| = |z - b|
 * @param a 第一个复数点
 * @param b 第二个复数点
 * @param z_var 复变量名，默认为 "z"
 * @return 轨迹方程的符号表达式
 */
std::shared_ptr<SymbolicExpr> complex_locus_perpendicular_bisector(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    const std::string& z_var = "z"
);

}
