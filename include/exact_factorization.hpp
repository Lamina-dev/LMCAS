#pragma once

#include "polynomial.hpp"
#include "modular_arithmetic.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "result.hpp"
#include "computation_context.hpp"
#include <cstdint>
#include <vector>

namespace LMCAS {

/** @brief Berlekamp 分解结果。 */
struct BerlekampResult {
    int64_t prime;                              /**< 使用的素数。 */
    std::vector<Polynomial<ModInt>> factors;   /**< 模 p 下的不可约因子。 */
    int null_space_dim;                         /**< 零空间维度，等于不可约因子数。 */
    std::vector<std::vector<int64_t>> null_space_basis; /**< 零空间基向量。 */

    BerlekampResult() : prime(-1), null_space_dim(0) {}
};

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
 * @brief 执行 Hensel 提升并诊断前置条件。
 * 空因子列表、非法素数、零或常数输入、零因子及模 p 因子乘积
 * 与原多项式不符时，通过 CasError 返回具体诊断。
 */
LMCAS_API HenselLiftResult hensel_lift_checked(
    const Polynomial<BigInt>& poly,
    const std::vector<Polynomial<ModInt>>& mod_factors,
    int64_t prime,
    int lift_bound);

using ZassenhausResult =
    Result<MathResult<std::vector<Polynomial<Rational>>>>;

/**
 * @brief 在计算预算内执行精确 Zassenhaus 因子组合。
 * 要求重构模数为正，且提升因子乘积与 @p poly 的整数归一化结果同余。
 * 候选子集按字典序枚举，支持任意因子数；预算耗尽时以 Inconclusive 返回精确部分分解。
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
    std::vector<BigInt> g;      /**< 第一个因子的系数。 */
    std::vector<BigInt> h;      /**< 第二个因子的系数。 */
    std::vector<BigInt> s;      /**< g 的 Bezout 系数。 */
    std::vector<BigInt> t;      /**< h 的 Bezout 系数。 */
    BigInt modulus;              /**< 当前模数。 */
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

/** @brief 无平方因子预处理结果。 */
struct TfSquareFreeResult {
    Polynomial<Rational> square_free;          /**< 无平方因子部分。 */
    Polynomial<Rational> repeated_factor;      /**< 重复因子 gcd(f, f')。 */
    bool had_repeated_factors;                 /**< 是否存在重复因子。 */

    TfSquareFreeResult() : square_free("x"), repeated_factor("x"), had_repeated_factors(false) {}
};

/**
 * @brief 计算多项式的无平方因子部分.
 *
 * 以 gcd(f, f') 分离重复因子，为 Berlekamp 分解提供无平方因子部分。
 *
 * @param[in] poly 输入的有理系数多项式
 * @return 无平方因子预处理结果,包含 square-free 部分和重复因子信息
 *
 * @internal
 */
LMCAS_API TfSquareFreeResult tf_square_free(const Polynomial<Rational>& poly);

}
