/**
 * @file transcendental_factor.cpp
 * @brief 混合超越方程不可约因式分解：换元检测与主入口实现。
 *
 * 本文件实现 Phase 1（换元检测）的核心逻辑：遍历表达式 AST，
 * 收集依赖目标变量的超越子表达式，去重后分配代数不定元。
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


namespace lamina {

static int zc_popcount(uint64_t mask) {
    int count = 0;
    while (mask) {
        count += static_cast<int>(mask & 1);
        mask >>= 1;
    }
    return count;
}

/**
 * @brief 将 BigInt 系数归约到对称表示 [-m/2, m/2)。
 *
 * @param[in] c 待归约的系数
 * @param[in] m 模数（正整数）
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
 * @brief 计算一组 Hensel 提升因子的乘积，系数模 prime_power 归约。
 *
 * 对给定的因子子集（由位掩码指定），逐个相乘并对每个系数取模归约到对称表示。
 *
 * @param[in] lifted_factors 所有提升后的因子
 * @param[in] mask           子集位掩码（第 i 位为 1 表示选取第 i 个因子）
 * @param[in] mod            模数 p^k
 * @return 子集因子乘积的系数向量（对称表示）
 * @internal
 */
static std::vector<BigInt> zc_subset_product(
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    uint64_t mask,
    const BigInt& mod) {

    std::vector<BigInt> product = {BigInt(1)};

    for (size_t i = 0; i < lifted_factors.size(); ++i) {
        if (!((mask >> i) & 1)) continue;

        const auto& factor_coeffs = lifted_factors[i].coeffs;
        if (factor_coeffs.empty()) return {};

        /// 多项式乘法
        size_t new_size = product.size() + factor_coeffs.size() - 1;
        std::vector<BigInt> new_product(new_size, BigInt(0));

        for (size_t a = 0; a < product.size(); ++a) {
            if (product[a].is_zero()) continue;
            for (size_t b = 0; b < factor_coeffs.size(); ++b) {
                if (factor_coeffs[b].is_zero()) continue;
                new_product[a + b] = new_product[a + b] + product[a] * factor_coeffs[b];
            }
        }

        /// 对称归约
        for (auto& c : new_product) {
            c = zc_symmetric_mod(c, mod);
        }

        /// 去除高次零系数
        while (!new_product.empty() && new_product.back().is_zero()) {
            new_product.pop_back();
        }

        product = std::move(new_product);
    }

    return product;
}

/**
 * @brief 将 BigInt 系数向量转换为有理系数多项式。
 *
 * 每个 BigInt 系数直接转换为 Rational（分母为 1）。
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
 * @brief 对 BigInt 系数执行有理数重构。
 *
 * 给定整数 a 和模数 m，使用扩展欧几里得算法寻找有理数 p/q 满足：
 * - a ≡ p/q (mod m)
 * - |p| ≤ floor(sqrt(m/2)), |q| ≤ floor(sqrt(m/2))
 *
 * 当模数超出 int64_t 范围时使用 BigInt 算术；否则委托给
 * modular_arithmetic.hpp 中的 int64_t 版本以获得更好性能。
 *
 * @param[in] a 待重构的整数系数（已归约到对称表示）
 * @param[in] m 模数（正整数，通常为 p^k）
 * @param[out] num 输出分子
 * @param[out] den 输出分母
 * @return 满足界约束且 gcd=1 时返回 true；false 表示重构未决
 * @internal
 */
static bool zc_rational_reconstruction(
    const BigInt& a,
    const BigInt& m,
    BigInt& num,
    BigInt& den) {

    if (m.is_zero() || m.IsNegative()) return false;

    /// 对于可精确表示的小模数，委托给 int64_t 版本
    const auto a_small = a.try_to_int64();
    const auto m_small = m.try_to_int64();
    if (a_small && m_small) {
        const int64_t a_val = *a_small;
        const int64_t m_val = *m_small;

        auto reconstructed =
            lamina::rational_reconstruction_checked(a_val, m_val);
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
 * @brief 对子集乘积的所有系数执行有理重构，构造有理系数多项式。
 *
 * 对 product_coeffs 中的每个 BigInt 系数调用 zc_rational_reconstruction，
 * 若所有系数均成功重构，则返回对应的有理系数多项式。
 * 若任一系数重构失败，则回退到整数系数直接转换。
 *
 * @param[in] product_coeffs 子集乘积的 BigInt 系数向量（对称表示）
 * @param[in] mod            模数 p^k
 * @param[in] var            变量名
 * @return 有理系数多项式（通过有理重构或整数直接转换）
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
            /// 任一系数重构失败，回退到整数系数
            return zc_bigint_to_rational_poly(product_coeffs, var);
        }
    }

    return Polynomial<Rational>(rat_coeffs, var);
}

/**
 * @brief 检验候选因子是否整除原多项式（精确有理除法）。
 *
 * 执行多项式带余除法 f / candidate，若余数为零则整除。
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
 * @brief 使候选多项式首一化（首项系数归一）。
 *
 * @param[in] poly 输入多项式
 * @return 首一多项式
 * @internal
 */
static Polynomial<Rational> zc_make_primitive(const Polynomial<Rational>& poly) {
    if (poly.is_zero()) return poly;
    return poly.make_monic();
}

/**
 * @brief 计算多项式系数的 L1 范数（系数绝对值之和）。
 *
 * @param[in] coeffs BigInt 系数向量
 * @return L1 范数
 * @internal
 */
static BigInt zc_l1_norm(const std::vector<BigInt>& coeffs) {
    BigInt norm(0);
    for (const auto& c : coeffs) {
        norm = norm + c.Abs();
    }
    return norm;
}

/**
 * @brief 计算 Mignotte 界：用于判断候选因子系数的合理范围。
 *
 * Mignotte 界为 C(n, floor(n/2)) * ||f||_2 的近似上界。
 * 此处使用简化估计：2^n * max|coeff(f)| 作为保守上界。
 *
 * @param[in] poly 原多项式
 * @return Mignotte 界的近似值
 * @internal
 */
static BigInt zc_mignotte_bound(const Polynomial<Rational>& poly) {
    BigInt max_coeff(0);
    for (const auto& c : poly.coeffs) {
        BigInt abs_num = c.get_numerator().Abs();
        /// 取 |numerator| 作为系数大小的近似上界
        if (abs_num > max_coeff) {
            max_coeff = abs_num;
        }
    }

    int n = poly.degree();
    if (n <= 0) return max_coeff;

    /// 2^n * max_coeff 作为保守 Mignotte 界
    BigInt bound = max_coeff;
    for (int i = 0; i < n; ++i) {
        bound = bound * BigInt(2);
    }
    return bound;
}

/**
 * @brief 基于度数和范数剪枝的因子组合（用于因子数 > 15 的情形）。
 *
 * 当模因子数超过 15 时，本函数使用启发式剪枝控制子集搜索规模：
 *
 * 1. 度数剪枝：预计算每个提升因子的次数，仅枚举总度数为原多项式
 *    次数的真因数的子集；
 * 2. 范数剪枝：候选子集乘积的 L1 范数超过 2 倍 Mignotte 界时，
 *    直接进入下一个候选。
 *
 * @note 完整的 LLL 格基约化算法可进一步优化此步骤，将搜索复杂度
 *       从指数级降至多项式级。对于典型用例（15 < r ≤ 30），
 *       度数+范数剪枝已足够高效。
 *
 * @param[in] poly           有理系数原多项式
 * @param[in] lifted_factors Hensel 提升后的整系数因子
 * @param[in] mod            模数 p^k
 * @return 有理系数真因子列表
 *
 * @see van Hoeij, M. "Factoring polynomials and the knapsack problem."
 *      Journal of Number Theory, 95(2), 2002.
 * @internal
 */
static std::vector<Polynomial<Rational>> zc_lll_pruned_combine(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& mod) {

    std::vector<Polynomial<Rational>> true_factors;
    size_t r = lifted_factors.size();
    std::string var = poly.variable_name;

    /// 预计算每个提升因子的次数
    std::vector<int> factor_degrees(r);
    for (size_t i = 0; i < r; ++i) {
        factor_degrees[i] = static_cast<int>(lifted_factors[i].coeffs.size()) - 1;
        if (factor_degrees[i] < 0) factor_degrees[i] = 0;
    }

    /// 计算 Mignotte 界用于范数剪枝
    BigInt mignotte = zc_mignotte_bound(poly);
    BigInt norm_threshold = mignotte * BigInt(2);

    /// 工作副本
    Polynomial<Rational> remaining = poly;
    uint64_t active_mask = (1ULL << r) - 1;

    /// TODO: 完整的 LLL 格基约化实现可将此搜索从指数级优化为多项式级。
    /// 当前使用度数+范数启发式剪枝，对典型用例（15 < r ≤ 30）已足够高效。
    /// 参考：van Hoeij (2002) 的 LLL-based 因子组合算法。

    bool found_factor = true;
    while (found_factor) {
        found_factor = false;

        int active_count = zc_popcount(active_mask);
        if (active_count <= 1) break;
        if (remaining.degree() <= 1) break;

        int remaining_deg = remaining.degree();
        int max_subset_size = active_count / 2;

        for (int subset_size = 1; subset_size <= max_subset_size; ++subset_size) {
            uint64_t subset = (1ULL << subset_size) - 1;

            while (subset <= active_mask) {
                if ((subset & active_mask) == subset && zc_popcount(subset) == subset_size) {

                    /// --- 度数剪枝 ---
                    /// 计算子集因子的总度数
                    int subset_degree = 0;
                    for (size_t i = 0; i < r; ++i) {
                        if ((subset >> i) & 1) {
                            subset_degree += factor_degrees[i];
                        }
                    }

                    /// 跳过总度数超过剩余多项式度数一半的子集
                    if (subset_degree > remaining_deg / 2) {
                        goto next_pruned_subset;
                    }

                    /// 跳过总度数为 0 的子集（不可能产生有意义的因子）
                    if (subset_degree <= 0) {
                        goto next_pruned_subset;
                    }

                    {
                        /// 计算子集因子乘积 mod p^k
                        std::vector<BigInt> product_coeffs = zc_subset_product(
                            lifted_factors, subset, mod);

                        if (product_coeffs.empty()) {
                            goto next_pruned_subset;
                        }

                        /// --- 范数剪枝 ---
                        /// 若乘积的 L1 范数超过 2 倍 Mignotte 界，拒绝此候选
                        BigInt candidate_norm = zc_l1_norm(product_coeffs);
                        if (candidate_norm > norm_threshold) {
                            goto next_pruned_subset;
                        }

                        {
                            /// 有理重构
                            Polynomial<Rational> candidate = zc_reconstruct_candidate(
                                product_coeffs, mod, var);

                            candidate = zc_make_primitive(candidate);

                            /// 整除性检验
                            if (!candidate.is_zero() && candidate.degree() > 0 &&
                                candidate.degree() < remaining.degree() &&
                                zc_divides_exactly(remaining, candidate)) {

                                true_factors.push_back(candidate);

                                auto [quotient, rem] = remaining.div_mod(candidate);
                                remaining = quotient;

                                active_mask &= ~subset;
                                found_factor = true;
                                goto restart_pruned_enumeration;
                            }
                        }
                    }
                }

                next_pruned_subset:
                if (subset == 0) break;
                uint64_t c = subset & (-static_cast<int64_t>(subset));
                uint64_t rr = subset + c;
                subset = (((rr ^ subset) >> 2) / c) | rr;
                if (subset > active_mask || subset == 0) break;
            }
        }

        restart_pruned_enumeration:;
    }

    /// 剩余因子形成最后一个真因子
    if (!remaining.is_zero() && remaining.degree() > 0) {
        true_factors.push_back(remaining.make_monic());
    } else if (!remaining.is_zero() && remaining.degree() == 0) {
        if (remaining.coeffs[0] != Rational(1)) {
            true_factors.push_back(remaining);
        }
    }

    return true_factors;
}

/**
 * @brief Zassenhaus 因子组合：从 Hensel 提升后的因子中筛选真因子。
 *
 * 算法：
 * 1. 将原多项式转为首一形式
 * 2. 按子集大小从 1 到 floor(r/2) 枚举 lifted_factors 的子集
 * 3. 对每个子集，计算其乘积 mod p^k（对称表示）
 * 4. 对乘积系数执行有理重构（rational_reconstruction），恢复有理系数
 * 5. 首一化后检验候选因子是否整除当前剩余多项式
 * 6. 若整除，记录为真因子，从因子池中移除已用因子，更新剩余多项式
 * 7. 剩余因子形成最后一个真因子
 *
 * 当因子数超过 15 时，使用度数和范数启发式剪枝加速搜索。
 *
 * @param[in] poly           有理系数原多项式
 * @param[in] lifted_factors Hensel 提升后的整系数因子
 * @param[in] prime_power    素数幂 p^k
 * @return 有理系数真因子列表
 *
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 */
std::vector<Polynomial<Rational>> zassenhaus_combine(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    int64_t prime_power) {

    std::vector<Polynomial<Rational>> true_factors;

    /// 边界情形
    if (poly.is_zero() || lifted_factors.empty()) {
        if (!poly.is_zero()) true_factors.push_back(poly);
        return true_factors;
    }

    /// 单因子：原多项式本身不可约
    if (lifted_factors.size() == 1) {
        true_factors.push_back(poly);
        return true_factors;
    }

    /// 因子数上限
    size_t r = lifted_factors.size();
    if (r > 30) {
        /// 因子数过多（> 30），即使剪枝也不实际，返回原多项式
        true_factors.push_back(poly);
        return true_factors;
    }

    /// 当因子数 > 15 时，使用度数+范数剪枝的组合策略
    if (r > 15) {
        BigInt mod(static_cast<long long>(prime_power));
        return zc_lll_pruned_combine(poly, lifted_factors, mod);
    }

    BigInt mod(static_cast<long long>(prime_power));
    std::string var = poly.variable_name;

    /// 工作副本：当前剩余多项式和因子池
    Polynomial<Rational> remaining = poly;

    /// 活跃因子索引集合（用位掩码表示）
    uint64_t active_mask = (1ULL << r) - 1;  // 所有因子初始活跃

    /// 按子集大小从 1 到 floor(active_count/2) 枚举
    bool found_factor = true;
    while (found_factor) {
        found_factor = false;

        int active_count = zc_popcount(active_mask);

        /// 早期终止：仅剩 1 个活跃因子，剩余多项式本身即为不可约因子
        if (active_count <= 1) break;

        /// 早期终止：剩余多项式为线性或常数，必然不可约
        if (remaining.degree() <= 1) break;

        int max_subset_size = active_count / 2;

        for (int subset_size = 1; subset_size <= max_subset_size; ++subset_size) {
            /// 枚举所有大小为 subset_size 的活跃因子子集
            /// 使用位掩码枚举：遍历 active_mask 的所有子集中 popcount == subset_size 的
            uint64_t subset = 0;

            /// Gosper's hack 初始化：最小的 subset_size 位子集
            subset = (1ULL << subset_size) - 1;

            while (subset <= active_mask) {
                /// 检查 subset 是否为 active_mask 的子集
                if ((subset & active_mask) == subset && zc_popcount(subset) == subset_size) {
                    /// 计算子集因子的乘积 mod p^k
                    std::vector<BigInt> product_coeffs = zc_subset_product(
                        lifted_factors, subset, mod);

                    if (product_coeffs.empty()) {
                        /// 跳过空乘积
                        goto next_subset;
                    }

                    {
                        /// 通过有理重构构造候选因子
                        Polynomial<Rational> candidate = zc_reconstruct_candidate(
                            product_coeffs, mod, var);

                        /// 首一化
                        candidate = zc_make_primitive(candidate);

                        /// 检验整除性
                        if (!candidate.is_zero() && candidate.degree() > 0 &&
                            candidate.degree() < remaining.degree() &&
                            zc_divides_exactly(remaining, candidate)) {

                            /// 找到真因子
                            true_factors.push_back(candidate);

                            /// 更新剩余多项式
                            auto [quotient, rem] = remaining.div_mod(candidate);
                            remaining = quotient;

                            /// 从活跃集中移除已用因子
                            active_mask &= ~subset;

                            found_factor = true;
                            goto restart_enumeration;
                        }
                    }
                }

                next_subset:
                /// 下一个子集（Gosper's hack）
                if (subset == 0) break;
                uint64_t c = subset & (-static_cast<int64_t>(subset));
                uint64_t rr = subset + c;
                subset = (((rr ^ subset) >> 2) / c) | rr;

                /// active_mask 为子集枚举设置数值上界。
                if (subset > active_mask || subset == 0) break;
            }
        }

        restart_enumeration:;
    }

    /// 剩余因子形成最后一个真因子
    if (!remaining.is_zero() && remaining.degree() > 0) {
        true_factors.push_back(remaining.make_monic());
    } else if (!remaining.is_zero() && remaining.degree() == 0) {
        /// 常数因子：若非 1 则记录
        if (remaining.coeffs[0] != Rational(1)) {
            true_factors.push_back(remaining);
        }
    }

    return true_factors;
}

} // namespace lamina
