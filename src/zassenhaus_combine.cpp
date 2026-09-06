/**
 * @file transcendental_factor.cpp
 * @brief 混合超越方程不可约因式分解:换元检测与主入口实现.
 *
 * 本文件实现 Phase 1(换元检测)的核心逻辑:遍历表达式 AST,
 * 收集依赖目标变量的超越子表达式,去重后分配代数不定元.
 */

#include "transcendental_factor.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>


namespace LMCAS {


/**
 * @brief 将 BigInt 系数归约到对称表示 [-m/2, m/2).
 *
 * @param[in] c 待归约的系数
 * @param[in] m 模数(正整数)
 * @return 对称表示下的归约值
 * @internal
 */
static BigInt zc_symmetric_mod(const BigInt& c, const BigInt& m) {
    if (m.is_zero()) return c;

    BigInt r = c % m;
    if (r.IsNegative()) {
        r = r + m;
    }

    BigInt half_m = m / BigInt(2);
    if (r > half_m) {
        r = r - m;
    }
    return r;
}

/**
 * 在模 m 下乘以按字典序枚举的提升因子子集.
 */
static Result<std::vector<BigInt>> zc_subset_product(
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const std::vector<size_t>& indices,
    const BigInt& mod,
    ComputationContext& context) {
    std::vector<BigInt> product = {BigInt(1)};

    for (size_t index : indices) {
        auto step = context.consume_steps(1, "zassenhaus_combine");
        if (!step) {
            return Result<std::vector<BigInt>>::failure(step.error());
        }
        if (index >= lifted_factors.size() ||
            lifted_factors[index].coeffs.empty()) {
            return Result<std::vector<BigInt>>::failure(
                CasErrc::InvalidArgument, "invalid lifted-factor index or zero factor",
                "zassenhaus_combine");
        }
        const auto& factor_coeffs = lifted_factors[index].coeffs;
        std::vector<BigInt> next(product.size() + factor_coeffs.size() - 1,
                                 BigInt(0));
        for (size_t i = 0; i < product.size(); ++i) {
            if (product[i].is_zero()) continue;
            for (size_t j = 0; j < factor_coeffs.size(); ++j) {
                if (factor_coeffs[j].is_zero()) continue;
                next[i + j] = next[i + j] + product[i] * factor_coeffs[j];
            }
        }
        for (auto& coefficient : next) {
            coefficient = zc_symmetric_mod(coefficient, mod);
        }
        while (!next.empty() && next.back().is_zero()) next.pop_back();
        product = std::move(next);
    }
    return Result<std::vector<BigInt>>::success(std::move(product));
}

/**
 * @brief 将 BigInt 系数向量转换为有理系数多项式.
 *
 * 每个 BigInt 系数直接转换为 Rational(分母为 1).
 *
 * @param[in] coeffs BigInt 系数向量
 * @param[in] var    变量名
 * @return 有理系数多项式
 * @internal
 */
static Polynomial<Rational> zc_bigint_to_rational_poly(
    const std::vector<BigInt>& coeffs,
    const std::string& var) {

    std::vector<Rational> rat_coeffs;
    rat_coeffs.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        rat_coeffs.emplace_back(c);
    }
    return Polynomial<Rational>(rat_coeffs, var);
}

/**
 * @brief 对 BigInt 系数执行有理数重构.
 *
 * 给定整数 a 和模数 m,使用扩展欧几里得算法寻找有理数 p/q 满足:
 * - a == p/q (mod m)
 * - |p| <= floor(sqrt(m/2)), |q| <= floor(sqrt(m/2))
 *
 * 当模数超出 int64_t 范围时使用 BigInt 算术;否则委托给
 * modular_arithmetic.hpp 中的 int64_t 版本以获得更好性能.
 *
 * @param[in] a 待重构的整数系数(已归约到对称表示)
 * @param[in] m 模数(正整数,通常为 p^k)
 * @param[out] num 输出分子
 * @param[out] den 输出分母
 * @return 满足界约束且 gcd=1 时返回 true;false 表示重构未决
 * @internal
 */
static bool zc_rational_reconstruction(
    const BigInt& a,
    const BigInt& m,
    BigInt& num,
    BigInt& den) {

    if (m.is_zero() || m.IsNegative()) return false;

    /// 对于可精确表示的小模数,委托给 int64_t 版本
    const auto a_small = a.try_to_int64();
    const auto m_small = m.try_to_int64();
    if (a_small && m_small) {
        const int64_t a_val = *a_small;
        const int64_t m_val = *m_small;

        auto reconstructed =
            LMCAS::rational_reconstruction_checked(a_val, m_val);
        if (!reconstructed) return false;
        const auto [p, q] = reconstructed.value();
        num = BigInt(static_cast<long long>(p));
        den = BigInt(static_cast<long long>(q));
        return true;
    }

    /// BigInt 版本的有理重构
    /// 将 a 归约到 [0, m)
    BigInt a_mod = a % m;
    if (a_mod.IsNegative()) {
        a_mod = a_mod + m;
    }

    /// 计算界 bound = floor(sqrt(m/2))
    BigInt half_m = m / BigInt(2);
    BigInt bound = half_m.sqrt();
    if (bound.is_zero()) bound = BigInt(1);

    /// 扩展欧几里得算法
    BigInt r0 = m, r1 = a_mod;
    BigInt s0 = BigInt(0), s1 = BigInt(1);

    while (r1 > bound) {
        BigInt q = r0 / r1;
        BigInt r_new = r0 - q * r1;
        BigInt s_new = s0 - q * s1;

        r0 = r1;
        r1 = r_new;
        s0 = s1;
        s1 = s_new;
    }

    BigInt p = r1;
    BigInt q_val = s1;

    /// 确保分母为正
    if (q_val.IsNegative()) {
        p = -p;
        q_val = -q_val;
    }

    /// 验证界约束
    if (p.Abs() > bound || q_val.is_zero() || q_val > bound) {
        return false;
    }

    /// 验证 gcd(|p|, q) == 1
    BigInt g = BigInt::gcd(p.Abs(), q_val);
    if (g != BigInt(1)) {
        return false;
    }

    num = p;
    den = q_val;
    return true;
}

/**
 * @brief 对子集乘积的所有系数执行有理重构,构造有理系数多项式.
 *
 * 对 product_coeffs 中的每个 BigInt 系数调用 zc_rational_reconstruction,
 * 若所有系数均成功重构,则返回对应的有理系数多项式.
 * 若任一系数重构失败,则回退到整数系数直接转换.
 *
 * @param[in] product_coeffs 子集乘积的 BigInt 系数向量(对称表示)
 * @param[in] mod            模数 p^k
 * @param[in] var            变量名
 * @return 有理系数多项式(通过有理重构或整数直接转换)
 * @internal
 */
static Polynomial<Rational> zc_reconstruct_candidate(
    const std::vector<BigInt>& product_coeffs,
    const BigInt& mod,
    const std::string& var) {

    std::vector<Rational> rat_coeffs;
    rat_coeffs.reserve(product_coeffs.size());

    /// 尝试对每个系数执行有理重构
    for (const auto& c : product_coeffs) {
        BigInt num, den;
        if (zc_rational_reconstruction(c, mod, num, den)) {
            rat_coeffs.emplace_back(num, den);
        } else {
            /// 任一系数重构失败,回退到整数系数
            return zc_bigint_to_rational_poly(product_coeffs, var);
        }
    }

    return Polynomial<Rational>(rat_coeffs, var);
}

/**
 * @brief 检验候选因子是否整除原多项式(精确有理除法).
 *
 * 执行多项式带余除法 f / candidate,若余数为零则整除.
 *
 * @param[in] f         原多项式
 * @param[in] candidate 候选因子
 * @return 整除返回 true
 * @internal
 */
static bool zc_divides_exactly(
    const Polynomial<Rational>& f,
    const Polynomial<Rational>& candidate) {

    if (candidate.is_zero()) return false;
    if (candidate.degree() > f.degree()) return false;
    if (candidate.degree() == 0) return true;

    auto [quotient, remainder] = f.div_mod(candidate);
    return remainder.is_zero();
}

/**
 * @brief 使候选多项式首一化(首项系数归一).
 *
 * @param[in] poly 输入多项式
 * @return 首一多项式
 * @internal
 */
static Polynomial<Rational> zc_make_primitive(const Polynomial<Rational>& poly) {
    if (poly.is_zero()) return poly;
    return poly.make_monic();
}

static bool zc_next_combination(std::vector<size_t>& positions, size_t count) {
    const size_t width = positions.size();
    for (size_t cursor = width; cursor > 0; --cursor) {
        const size_t i = cursor - 1;
        if (positions[i] < count - width + i) {
            ++positions[i];
            for (size_t j = i + 1; j < width; ++j) {
                positions[j] = positions[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

static std::vector<size_t> zc_select_indices(
    const std::vector<size_t>& active,
    const std::vector<size_t>& positions) {
    std::vector<size_t> indices;
    indices.reserve(positions.size());
    for (size_t position : positions) indices.push_back(active[position]);
    return indices;
}

static Result<void> zc_validate_lifted_product(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& modulus,
    ComputationContext& context) {
    if (lifted_factors.empty()) return Result<void>::success();

    std::vector<size_t> all_indices(lifted_factors.size());
    for (size_t i = 0; i < all_indices.size(); ++i) all_indices[i] = i;
    auto product_result =
        zc_subset_product(lifted_factors, all_indices, modulus, context);
    if (!product_result) return Result<void>::failure(product_result.error());
    auto product = std::move(product_result.value());

    BigInt denominator_lcm(1);
    for (const Rational& coefficient : poly.coeffs) {
        denominator_lcm =
            BigInt::lcm(denominator_lcm, coefficient.get_denominator());
    }
    std::vector<BigInt> expected;
    expected.reserve(poly.coeffs.size());
    for (const Rational& coefficient : poly.coeffs) {
        BigInt integer_coefficient =
            coefficient.get_numerator() *
            (denominator_lcm / coefficient.get_denominator());
        expected.push_back(zc_symmetric_mod(integer_coefficient, modulus));
    }
    while (!expected.empty() && expected.back().is_zero()) expected.pop_back();
    while (!product.empty() && product.back().is_zero()) product.pop_back();
    if (!(product == expected)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "lifted factors are not congruent to the normalized polynomial",
            "zassenhaus_combine");
    }
    return Result<void>::success();
}

static ZassenhausResult zc_finalize(
    const Polynomial<Rational>& original,
    std::vector<Polynomial<Rational>> factors,
    Polynomial<Rational> remaining,
    Completeness completeness,
    std::string reason) {
    if (!remaining.is_zero() && remaining.degree() > 0) {
        factors.push_back(remaining.make_monic());
    }

    if (!original.is_zero()) {
        Polynomial<Rational> product({Rational(1)}, original.variable_name);
        for (const auto& factor : factors) product = product * factor;
        if (!(product == original.make_monic())) {
            return ZassenhausResult::failure(
                CasErrc::InternalInvariant,
                "accepted Zassenhaus factors do not reconstruct the polynomial",
                "zassenhaus_combine");
        }
    } else if (!factors.empty()) {
        return ZassenhausResult::failure(
            CasErrc::InternalInvariant,
            "zero polynomial cannot have a nonempty monic factorization",
            "zassenhaus_combine");
    }

    return ZassenhausResult::success(
        MathResult<std::vector<Polynomial<Rational>>>{
            std::move(factors), completeness, std::move(reason)});
}

static ZassenhausResult zassenhaus_combine_impl(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& reconstruction_modulus,
    ComputationContext& context) {
    if (reconstruction_modulus.is_zero() ||
        reconstruction_modulus.IsNegative()) {
        return ZassenhausResult::failure(
            CasErrc::InvalidArgument,
            "reconstruction modulus must be positive",
            "zassenhaus_combine");
    }
    if (poly.is_zero()) {
        if (!lifted_factors.empty()) {
            return ZassenhausResult::failure(
                CasErrc::InvalidArgument,
                "zero polynomial requires an empty lifted-factor list",
                "zassenhaus_combine");
        }
        return ZassenhausResult::success(
            MathResult<std::vector<Polynomial<Rational>>>{{},
                Completeness::Complete, {}});
    }
    if (lifted_factors.empty()) {
        return zc_finalize(poly, {}, poly.make_monic(),
                           Completeness::Inconclusive,
                           "no lifted factors were provided");
    }

    auto validation = zc_validate_lifted_product(
        poly, lifted_factors, reconstruction_modulus, context);
    if (!validation) {
        if (validation.error().code == CasErrc::ResourceLimit) {
            return zc_finalize(poly, {}, poly.make_monic(),
                               Completeness::Inconclusive,
                               "budget exhausted while validating lifted factors");
        }
        return ZassenhausResult::failure(validation.error());
    }

    std::vector<Polynomial<Rational>> true_factors;
    std::vector<size_t> active(lifted_factors.size());
    for (size_t i = 0; i < active.size(); ++i) active[i] = i;
    Polynomial<Rational> remaining = poly.make_monic();

    for (;;) {
        bool found = false;
        if (active.size() <= 1 || remaining.degree() <= 1) break;

        for (size_t subset_size = 1;
             subset_size <= active.size() / 2 && !found;
             ++subset_size) {
            std::vector<size_t> positions(subset_size);
            for (size_t i = 0; i < subset_size; ++i) positions[i] = i;
            do {
                auto candidate_step =
                    context.consume_steps(1, "zassenhaus_combine");
                if (!candidate_step) {
                    if (candidate_step.error().code == CasErrc::ResourceLimit) {
                        return zc_finalize(
                            poly, std::move(true_factors), std::move(remaining),
                            Completeness::Inconclusive,
                            "budget exhausted during bounded combination enumeration");
                    }
                    return ZassenhausResult::failure(candidate_step.error());
                }

                std::vector<size_t> indices =
                    zc_select_indices(active, positions);
                auto product_result = zc_subset_product(
                    lifted_factors, indices, reconstruction_modulus, context);
                if (!product_result) {
                    if (product_result.error().code == CasErrc::ResourceLimit) {
                        return zc_finalize(
                            poly, std::move(true_factors), std::move(remaining),
                            Completeness::Inconclusive,
                            "budget exhausted during candidate multiplication");
                    }
                    return ZassenhausResult::failure(product_result.error());
                }
                if (product_result.value().empty()) continue;

                Polynomial<Rational> candidate = zc_make_primitive(
                    zc_reconstruct_candidate(product_result.value(),
                                             reconstruction_modulus,
                                             poly.variable_name));
                if (candidate.is_zero() || candidate.degree() <= 0 ||
                    candidate.degree() >= remaining.degree() ||
                    !zc_divides_exactly(remaining, candidate)) {
                    continue;
                }

                auto [quotient, remainder] = remaining.div_mod(candidate);
                if (!remainder.is_zero()) {
                    return ZassenhausResult::failure(
                        CasErrc::InternalInvariant,
                        "exact divisor produced a nonzero remainder",
                        "zassenhaus_combine");
                }
                true_factors.push_back(candidate);
                remaining = std::move(quotient);
                for (size_t cursor = positions.size(); cursor > 0; --cursor) {
                    active.erase(active.begin() +
                                 static_cast<std::ptrdiff_t>(positions[cursor - 1]));
                }
                found = true;
                break;
            } while (zc_next_combination(positions, active.size()));
        }
        if (!found) break;
    }

    return zc_finalize(poly, std::move(true_factors), std::move(remaining),
                       Completeness::Complete, {});
}

ZassenhausResult zassenhaus_combine_checked(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& reconstruction_modulus,
    ComputationContext& context) {
    try {
        return zassenhaus_combine_impl(
            poly, lifted_factors, reconstruction_modulus, context);
    } catch (const std::bad_alloc&) {
        return ZassenhausResult::failure(
            CasErrc::ResourceLimit, "allocation failed", "zassenhaus_combine");
    } catch (const std::exception& error) {
        return ZassenhausResult::failure(
            CasErrc::InternalInvariant, error.what(), "zassenhaus_combine");
    }
}

ZassenhausResult zassenhaus_combine_checked(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& reconstruction_modulus) {
    ComputationContext context;
    return zassenhaus_combine_checked(
        poly, lifted_factors, reconstruction_modulus, context);
}

} // namespace LMCAS
