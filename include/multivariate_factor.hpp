/**
 * @file multivariate_factor.hpp
 * @brief 多元多项式因式分解器公共 API。
 *
 * 采用 Wang-EEZ（Enhanced Extended Zassenhaus）算法，将多元多项式分解为
 * 不可约因子的乘积。核心流程：预处理（容度/本原/无平方）→ 求值同态降维 →
 * 一元分解 → 多元 Hensel 提升 → 试除验证。
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." Chapter 15.
 */
#pragma once

#include "multivariate_poly.hpp"

#include <optional>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#ifdef LAMINA_CORE_EXPORTS
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif

namespace lamina {

/**
 * @brief 多元因式分解结果
 *
 * 满足不变量：constant * ∏(factors[i] ^ multiplicities[i]) == 原多项式。
 * 每个 factors[i] 为本原不可约多项式，首项系数为正。
 */
struct MultiFactorResult {
    Rational constant;               ///< 数值常数因子
    std::vector<MultiPoly> factors;  ///< 不可约因子列表（各因子本原）
    std::vector<int> multiplicities; ///< 各因子的重数
};

/**
 * @brief 多元多项式因式分解主入口
 *
 * 采用 Wang-EEZ 算法，将多元多项式分解为不可约因子的乘积。
 * 返回 constant * factors[0]^mult[0] * factors[1]^mult[1] * ... = poly。
 *
 * 对于零多项式返回 {0, [], []}；对于常数多项式返回 {constant, [], []}。
 * 内部错误时返回原多项式作为不可约因子，不抛异常。
 *
 * @param[in] poly 待分解的多元多项式
 * @return 分解结果
 */
LAMINA_API MultiFactorResult factor_multivariate(const MultiPoly& poly);

/// ---------------------------------------------------------------------------
/// 内部算法组件（供测试使用）
/// ---------------------------------------------------------------------------

/**
 * @brief 计算多元容度
 *
 * 将多项式视为主变量的一元多项式，计算所有系数多项式（辅助变量的多项式）的 GCD。
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 容度多项式（辅助变量的多项式）
 */
LAMINA_API MultiPoly multivariate_content(const MultiPoly& poly,
                                          const std::string& main_var);

/**
 * @brief 计算多元本原部分
 *
 * 将多项式除以其关于主变量的容度，得到本原部分。
 * 满足：multivariate_content(poly, var) * multivariate_primitive_part(poly, var) == poly（至多差符号）。
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 本原部分多项式
 */
LAMINA_API MultiPoly multivariate_primitive_part(const MultiPoly& poly,
                                                 const std::string& main_var);

/**
 * @brief 计算多元多项式 GCD
 *
 * 采用子结式 PRS（Polynomial Remainder Sequence）算法计算两个多元多项式的最大公因式。
 *
 * @param[in] a 第一个多项式
 * @param[in] b 第二个多项式
 * @return gcd(a, b)
 */
LAMINA_API MultiPoly multivariate_gcd(const MultiPoly& a, const MultiPoly& b);

/**
 * @brief 无平方因子分解结果
 *
 * components[i] 的重数为 i+1。
 * 满足：f = components[0] * components[1]² * components[2]³ * ...（至多差常数）。
 */
struct SquareFreeDecomp {
    std::vector<MultiPoly> components; ///< 各无平方分量，components[i] 重数为 i+1
};

/**
 * @brief 多元无平方因子分解
 *
 * 通过计算 gcd(f, ∂f/∂x_main) 检测重因子，将多项式分解为
 * f = f₁ · f₂² · f₃³ · ... 的形式，其中每个 fᵢ 无平方。
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 无平方因子分解结果
 */
LAMINA_API SquareFreeDecomp square_free_decompose(const MultiPoly& poly,
                                                  const std::string& main_var);

/**
 * @brief 多元 Hensel 提升（单变量步）
 *
 * 将一元因子逐变量提升，从 f(x₁, a₂, ..., aₙ) 的因子恢复到包含 lift_var 的因子。
 * 每次引入一个辅助变量，计算修正项直到达到次数上界。
 *
 * @param[in] poly              原始多元多项式
 * @param[in] univariate_factors 一元分解得到的因子列表
 * @param[in] lift_var          当前提升的变量名
 * @param[in] eval_point        该变量的求值点
 * @param[in] degree_bound      提升次数上界
 * @return 提升后的多元因子列表
 */
LAMINA_API std::vector<MultiPoly> multivariate_hensel_lift(
    const MultiPoly& poly,
    const std::vector<Polynomial<Rational>>& univariate_factors,
    const std::string& lift_var,
    const Rational& eval_point,
    int degree_bound);

/**
 * @brief 多元丢番图方程求解器
 *
 * 求解 s₁·f₁ + s₂·f₂ + ... + sᵣ·fᵣ ≡ target (mod ideal)，
 * 其中 ideal 由 (var - eval_point)^degree_bound 生成。
 * 用于 Hensel 提升过程中计算修正项。
 *
 * @param[in] factors      互素因子列表
 * @param[in] target       目标多项式
 * @param[in] var          变量名
 * @param[in] eval_point   求值点
 * @param[in] degree_bound 次数上界
 * @return 解多项式列表 [s₁, s₂, ..., sᵣ]
 */
LAMINA_API std::vector<MultiPoly> multivariate_diophantine(
    const std::vector<MultiPoly>& factors,
    const MultiPoly& target,
    const std::string& var,
    const Rational& eval_point,
    int degree_bound);

/// ---------------------------------------------------------------------------
/// 快速路径检测（detail 命名空间）
/// ---------------------------------------------------------------------------

namespace detail {

/**
 * @brief 检测是否为线性多项式
 *
 * 判断多项式关于主变量的次数是否 ≤ 1，线性多项式直接视为不可约。
 *
 * @param[in] poly     输入多项式
 * @param[in] main_var 主变量名
 * @return 若关于 main_var 的次数 ≤ 1 则返回 true
 */
LAMINA_API bool is_linear(const MultiPoly& poly, const std::string& main_var);

/**
 * @brief 检测差平方模式 a² - b²
 *
 * 若多项式可表示为两个多项式平方之差，返回 (a, b)，
 * 使得 poly = a² - b² = (a+b)(a-b)。
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (a, b)，否则返回 nullopt
 */
LAMINA_API std::optional<std::pair<MultiPoly, MultiPoly>>
detect_difference_of_squares(const MultiPoly& poly);

/**
 * @brief 检测二项式幂模式 xⁿ - yⁿ
 *
 * 若多项式为两个变量的纯幂之差（如 x³ - y³），返回变量名和幂次，
 * 以便应用分圆多项式分解公式。
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (变量1, 变量2, 幂次)，否则返回 nullopt
 */
LAMINA_API std::optional<std::tuple<std::string, std::string, int>>
detect_binomial_power(const MultiPoly& poly);

/**
 * @brief 提取公因子单项式
 *
 * 计算所有项的 GCD 单项式并提取，返回公因子单项式和剩余多项式。
 *
 * @param[in] poly 输入多项式
 * @return pair(公因子单项式, 提取后的多项式)
 */
LAMINA_API std::pair<Monomial, MultiPoly>
extract_common_monomial(const MultiPoly& poly);

/**
 * @brief 齐次二元多项式分解
 *
 * 对齐次二元多项式，利用 f(x,1) 的一元分解结果推导二元因子。
 * 仅当多项式恰好含两个变量且为齐次时适用。
 *
 * @param[in] poly 输入多项式（须为齐次二元）
 * @return 若适用则返回分解结果，否则返回 nullopt
 */
LAMINA_API std::optional<MultiFactorResult>
factor_homogeneous_bivariate(const MultiPoly& poly);

} // namespace detail

} // namespace lamina
