/**
 * @file transcendental_factor.hpp
 * @brief 混合超越方程不可约因式分解器公共 API.
 *
 * 将 Berlekamp/Zassenhaus 风格的因式分解算法从纯多项式域推广到含超越函数的
 * 表达式空间.核心流程:换元 -> 多项式构造 -> 模分解 -> Hensel 提升 -> 因子组合 -> 逆换元.
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "modular_arithmetic.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "result.hpp"
#include "computation_context.hpp"

#include <memory>
#include <string>
#include <vector>
#include <cstdint>

namespace LMCAS {

/// 换元映射条目:超越子表达式 -> 代数不定元名
struct TransSubstitution {
    std::shared_ptr<SymbolicExpr> trans_expr;  ///< 原始超越子表达式 (e.g., sin(x))
    std::string indeterminate;                  ///< 代数不定元名 (e.g., "u0")
};

/// 换元结果
struct TransSubstitutionResult {
    std::vector<TransSubstitution> mappings;    ///< 所有换元映射
    std::shared_ptr<SymbolicExpr> poly_expr;   ///< 换元后的多项式表达式
    std::vector<std::shared_ptr<SymbolicExpr>> constraints; ///< 代数约束 (e.g., u0^2+u1^2-1=0)
};

/// Berlekamp 分解结果
struct BerlekampResult {
    int64_t prime;                              ///< 使用的素数
    std::vector<Polynomial<ModInt>> factors;   ///< 模 p 下的不可约因子
    int null_space_dim;                         ///< 零空间维度(= 不可约因子数)
    std::vector<std::vector<int64_t>> null_space_basis; ///< 零空间基向量

    BerlekampResult() : prime(-1), null_space_dim(0) {}
};

/**
 * @brief 混合超越表达式不可约因式分解主入口.
 *
 * 对含超越函数(sin,cos,exp,ln 等)的表达式执行因式分解,返回不可约因子列表.
 * 若输入为纯多项式,则委托给现有多项式分解;若分解失败或超时,返回原表达式本身.
 *
 * @param[in] expr 待分解的符号表达式
 * @param[in] var  目标变量名
 * @return 不可约因子的列表(乘积等于原表达式,可能含常数因子)
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

/**
 * @brief 检测表达式中的超越函数换元模式.
 *
 * 遍历表达式 AST,识别所有依赖目标变量的超越子表达式,为每个分配代数不定元,
 * 并记录不定元之间的代数约束(如三角恒等式 u_sin^2 + u_cos^2 = 1).
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 换元结果,包含映射表,换元后多项式表达式及约束列表
 */
LMCAS_API TransSubstitutionResult detect_trans_substitutions(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

/**
 * @brief Berlekamp 模素数分解(内部接口).
 *
 * 在有限域 F_p 上对有理系数多项式执行 Berlekamp 算法,
 * 返回模 p 下的不可约因子列表.
 *
 * @param[in] poly  有理系数多项式
 * @param[in] prime 选定的素数 p
 * @return Berlekamp 分解结果
 *
 * @internal
 */
LMCAS_API BerlekampResult berlekamp_factor(
    const Polynomial<Rational>& poly,
    int64_t prime);

using HenselLiftResult = Result<std::vector<Polynomial<BigInt>>>;

/**
 * @brief 执行带显式前置条件诊断的 Hensel 提升.
 *
 * 空因子列表,非法素数,零或常数输入,多项式零因子,以及模 p 因子乘积
 * 与原多项式不一致等输入均通过 CasError 返回具体诊断.
 */
LMCAS_API HenselLiftResult hensel_lift_checked(
    const Polynomial<BigInt>& poly,
    const std::vector<Polynomial<ModInt>>& mod_factors,
    int64_t prime,
    int lift_bound);

using ZassenhausResult =
    Result<MathResult<std::vector<Polynomial<Rational>>>>;

/**
 * @brief 受检且有界的精确 Zassenhaus 因子组合.
 *
 * 要求重构模数为正,且提升因子乘积与 @p poly 的整数归一化结果同余.
 * 候选子集按字典序枚举,不受机器字宽限制.计算预算耗尽时,精确部分分解
 * 以 `Inconclusive` 返回.
 */
LMCAS_API ZassenhausResult zassenhaus_combine_checked(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& reconstruction_modulus,
    ComputationContext& context);

LMCAS_API ZassenhausResult zassenhaus_combine_checked(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& reconstruction_modulus);

/**
 * @brief 二因子二次 Hensel 提升的状态结构.
 *
 * 存储一对因子 g, h 及其 Bezout 系数 s, t,满足:
 * - f == g * h (mod modulus)
 * - s * g + t * h == 1 (mod modulus)
 *
 * 系数按升幂存储:coeffs[i] 对应 x^i 的系数.
 */
struct HenselLiftPair {
    std::vector<BigInt> g;      ///< 第一个因子系数向量
    std::vector<BigInt> h;      ///< 第二个因子系数向量
    std::vector<BigInt> s;      ///< g 的 Bezout 系数向量
    std::vector<BigInt> t;      ///< h 的 Bezout 系数向量
    BigInt modulus;              ///< 当前模数
};

/**
 * @brief 执行一步二次 Hensel 提升:mod m -> mod m^2.
 *
 * 给定 f == g*h (mod m) 且 s*g + t*h == 1 (mod m),
 * 计算 g', h', s', t' 使得 f == g'*h' (mod m^2) 且 s'*g' + t'*h' == 1 (mod m^2).
 *
 * @param[in] f       原始多项式系数向量(升幂排列)
 * @param[in] current 当前提升状态(g, h, s, t, modulus=m)
 * @return 提升后的状态(g', h', s', t', modulus=m^2)
 *
 * @pre f == g*h (mod m)
 * @pre s*g + t*h == 1 (mod m)
 * @pre deg(s) < deg(h), deg(t) < deg(g)
 *
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 */
LMCAS_API HenselLiftPair hl_two_factor_lift(
    const std::vector<BigInt>& f,
    const HenselLiftPair& current);

/// 无平方因子预处理结果
struct TfSquareFreeResult {
    Polynomial<Rational> square_free;          ///< 无平方因子部分
    Polynomial<Rational> repeated_factor;      ///< 重复因子 gcd(f, f')
    bool had_repeated_factors;                 ///< 是否存在重复因子

    TfSquareFreeResult() : square_free("x"), repeated_factor("x"), had_repeated_factors(false) {}
};

/**
 * @brief 计算多项式的无平方因子部分.
 *
 * 通过计算 gcd(f, f') 检测并去除重复因子,返回无平方因子多项式.
 * 用于模分解前的预处理,确保 Berlekamp 算法的输入满足 square-free 条件.
 *
 * @param[in] poly 输入的有理系数多项式
 * @return 无平方因子预处理结果,包含 square-free 部分和重复因子信息
 *
 * @internal
 */
LMCAS_API TfSquareFreeResult tf_square_free(const Polynomial<Rational>& poly);

/// 多项式构造结果
struct TfPolyBuildResult {
    bool success;                              ///< 转换是否成功
    Polynomial<Rational> poly;                 ///< 主变量多项式(单变量或主变量策略)
    std::string main_variable;                 ///< 选定的主变量名
    std::vector<std::string> param_variables;  ///< 参数变量列表(非主变量的不定元)

    TfPolyBuildResult() : success(false), poly("x") {}
};

/**
 * @brief 将换元后的表达式构造为有理系数多项式.
 *
 * 对换元后的多项式表达式,选择次数最高的不定元作为主变量,
 * 通过 symbolic_to_poly 转换为 Polynomial<Rational>.
 * 对于多元情形,非主变量将被视为参数(其系数需为有理数).
 *
 * @param[in] poly_expr       换元后的符号表达式
 * @param[in] indeterminates  不定元名称列表(如 {"u0", "u1"})
 * @param[in] original_var    原始目标变量名(如 "x")
 * @return 多项式构造结果
 *
 * @internal
 */
LMCAS_API TfPolyBuildResult tf_build_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& indeterminates,
    const std::string& original_var);

/**
 * @brief 对逆换元后的因子列表执行化简与常数乘子提取.
 *
 * 对每个因子调用 simplify() 规范化,提取数值前导系数,
 * 将所有常数乘子合并为单一数值因子.
 *
 * @param[in,out] factors 因子列表
 * @return 化简后的因子列表(可能含首位常数因子)
 *
 * @internal
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors(
    std::vector<std::shared_ptr<SymbolicExpr>>& factors);

} // namespace LMCAS
