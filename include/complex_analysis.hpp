/**
 * @file complex_analysis.hpp
 * @brief 复变函数分析：留数计算、柯西积分公式、解析延拓。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

using ComplexExprResult = Result<std::shared_ptr<SymbolicExpr>>;
using ComplexBoolResult = Result<bool>;

/**
 * @brief 计算复变函数在某极点处的留数
 * @param f 被积函数
 * @param z 复变量名
 * @param z0 极点
 * @param order 极点的阶数（默认为1，简单极点）
 * @return 留数的符号表达式
 */
LAMINA_API ComplexExprResult calculate_residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文计算留数，并显式报告无效输入和未覆盖域。
 */
LAMINA_API ComplexExprResult calculate_residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order = 1
);

LAMINA_API std::shared_ptr<SymbolicExpr> calculate_residue(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order = 1
);

/**
 * @brief 使用柯西积分公式计算闭合路径积分 ∮ f(z)/(z-z0)^n dz
 * @param f 解析函数部分
 * @param z 复变量名
 * @param z0 奇点
 * @param n 阶数
 * @return 积分结果的符号表达式
 */
LAMINA_API ComplexExprResult cauchy_integral_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文应用柯西积分公式，并显式报告无效输入和未覆盖域。
 */
LAMINA_API ComplexExprResult cauchy_integral_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n = 1
);

LAMINA_API std::shared_ptr<SymbolicExpr> cauchy_integral(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n = 1
);

/**
 * @brief 尝试对函数进行解析延拓
 * @param f 输入函数
 * @param z 复变量名
 * @return 延拓后的函数（或返回自身若无法显式延拓）
 */
LAMINA_API std::shared_ptr<SymbolicExpr> analytic_continuation(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z
);

/**
 * @brief 提取表达式的实部 Re(expr)。
 *
 * 将复变量 z = x + iy（或显式 i 单位）展开，分离实部。
 * 对实值表达式直接返回自身。
 *
 * @param expr 输入表达式
 * @return 实部表达式
 */
LAMINA_API ComplexExprResult real_part_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文提取表达式实部，并显式报告无效输入。
 */
LAMINA_API ComplexExprResult real_part_checked(
    const std::shared_ptr<SymbolicExpr>& expr
);

LAMINA_API std::shared_ptr<SymbolicExpr> real_part(
    const std::shared_ptr<SymbolicExpr>& expr
);

/**
 * @brief 提取表达式的虚部 Im(expr)。
 * @param expr 输入表达式
 * @return 虚部表达式
 */
LAMINA_API ComplexExprResult imag_part_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文提取表达式虚部，并显式报告无效输入。
 */
LAMINA_API ComplexExprResult imag_part_checked(
    const std::shared_ptr<SymbolicExpr>& expr
);

LAMINA_API std::shared_ptr<SymbolicExpr> imag_part(
    const std::shared_ptr<SymbolicExpr>& expr
);

/**
 * @brief 计算表达式的复共轭 conj(expr)，将 i 替换为 -i。
 * @param expr 输入表达式
 * @return 共轭表达式
 */
LAMINA_API ComplexExprResult conjugate_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文计算共轭，并显式报告无效输入。
 */
LAMINA_API ComplexExprResult conjugate_checked(
    const std::shared_ptr<SymbolicExpr>& expr
);

LAMINA_API std::shared_ptr<SymbolicExpr> conjugate(
    const std::shared_ptr<SymbolicExpr>& expr
);

/**
 * @brief 检查函数 f(z) 在变量 z 上是否解析（满足 Cauchy-Riemann 方程）。
 *
 * 设 f = u(x,y) + i·v(x,y)，检查 ∂u/∂x = ∂v/∂y 且 ∂u/∂y = -∂v/∂x。
 *
 * @param f 复变函数表达式
 * @param z 复变量名（实部记为 z_re，虚部记为 z_im 进行分解）
 * @return 解析返回 true
 */
LAMINA_API ComplexBoolResult is_analytic_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    ComputationContext& context
);

/**
 * @brief 使用默认计算上下文检查解析性，并显式报告无效输入。
 */
LAMINA_API ComplexBoolResult is_analytic_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z
);

LAMINA_API bool is_analytic(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z
);

/**
 * @brief calculate_residue 的别名，匹配规范命名 residue。
 */
LAMINA_API ComplexExprResult residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context
);

LAMINA_API ComplexExprResult residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order = 1
);

LAMINA_API std::shared_ptr<SymbolicExpr> residue(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order = 1
);

} // namespace lamina
