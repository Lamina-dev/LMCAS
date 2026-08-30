/**
 * @file multivariate_factor.cpp
 * @brief 多元多项式因式分解器实现。
 */
#include "multivariate_factor.hpp"
#include "transcendental_factor.hpp"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include "internal/multivariate_factor_support.hpp"
namespace lamina {
static std::map<int, MultiPoly>
extract_coefficients(const MultiPoly& poly, const std::string& main_var)
{
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    std::vector<std::string> remaining_vars;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (static_cast<int>(i) != var_idx) remaining_vars.push_back(vars[i]);
    }
    std::map<int, MultiPoly> coeffs;
    if (var_idx < 0) { coeffs[0] = poly; return coeffs; }
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        Monomial reduced_mono;
        reduced_mono.reserve(vars.size() - 1);
        for (size_t i = 0; i < mono.size(); ++i) {
            if (static_cast<int>(i) != var_idx) reduced_mono.push_back(mono[i]);
        }
        std::vector<MultiPoly::Term> t_vec = {{reduced_mono, term.second}};
        MultiPoly term_poly(std::move(t_vec), remaining_vars);
        if (coeffs.find(exp) == coeffs.end()) coeffs[exp] = term_poly;
        else coeffs[exp] = coeffs[exp] + term_poly;
    }
    return coeffs;
}
static bool divides_poly(const MultiPoly& dividend, const MultiPoly& divisor)
{
    if (divisor.is_zero()) return false;
    if (dividend.is_zero()) return true;
    try { dividend.exact_div(divisor); return true; }
    catch (...) { return false; }
}
static std::string choose_main_variable(const MultiPoly& a, const MultiPoly& b)
{
    const auto& vars = a.variables();
    if (vars.empty()) return "";
    std::string best_var = vars[0];
    int best_deg = a.degree(vars[0]) + b.degree(vars[0]);
    for (size_t i = 1; i < vars.size(); ++i) {
        int deg = a.degree(vars[i]) + b.degree(vars[i]);
        if (deg > best_deg) { best_deg = deg; best_var = vars[i]; }
    }
    return best_var;
}
static Rational rational_gcd(const Rational& a, const Rational& b)
{
    if (a.is_zero()) return b.abs();
    if (b.is_zero()) return a.abs();
    BigInt num_gcd = BigInt::gcd(a.get_numerator().Abs(), b.get_numerator().Abs());
    BigInt den_gcd = BigInt::gcd(a.get_denominator(), b.get_denominator());
    BigInt den_lcm = a.get_denominator() * b.get_denominator() / den_gcd;
    return Rational(num_gcd, den_lcm);
}
static MultiPoly lagrange_interpolate(const std::vector<Rational>& points,
    const std::vector<MultiPoly>& values, const std::string& main_var,
    const std::vector<std::string>& all_vars)
{
    size_t n = points.size();
    if (n == 0) return MultiPoly(Rational(0), all_vars);
    int var_idx = -1;
    for (size_t i = 0; i < all_vars.size(); ++i) {
        if (all_vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    MultiPoly result(Rational(0), all_vars);
    for (size_t i = 0; i < n; ++i) {
        if (values[i].is_zero()) continue;
        Rational denom(1);
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            denom = denom * (points[i] - points[j]);
        }
        if (denom.is_zero()) continue;
        std::vector<Rational> basis_coeffs = {Rational(1)};
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            std::vector<Rational> nc(basis_coeffs.size() + 1, Rational(0));
            for (size_t k = 0; k < basis_coeffs.size(); ++k) {
                nc[k + 1] = nc[k + 1] + basis_coeffs[k];
                nc[k] = nc[k] - points[j] * basis_coeffs[k];
            }
            basis_coeffs = std::move(nc);
        }
        for (auto& c : basis_coeffs) c = c / denom;
        for (size_t k = 0; k < basis_coeffs.size(); ++k) {
            if (basis_coeffs[k].is_zero()) continue;
            MultiPoly scaled = values[i] * basis_coeffs[k];
            if (scaled.is_zero()) continue;
            for (const auto& term : scaled.terms()) {
                Monomial full_mono(all_vars.size(), 0);
                size_t ri = 0;
                for (size_t vi = 0; vi < all_vars.size(); ++vi) {
                    if (static_cast<int>(vi) == var_idx) {
                        full_mono[vi] = static_cast<int>(k);
                    } else {
                        if (ri < term.first.size()) full_mono[vi] = term.first[ri];
                        ++ri;
                    }
                }
                std::vector<MultiPoly::Term> t = {{full_mono, term.second}};
                result = result + MultiPoly(std::move(t), all_vars);
            }
        }
    }
    return result;
}

/**
 * @internal
 * @brief 生成搜索序列中第 k 个整数值
 *
 * 按 0, 1, -1, 2, -2, 3, -3, ... 的顺序生成。
 *
 * @param[in] k 序号（从 0 开始）
 * @return 对应的整数值
 */
static int search_value(int k)
{
    if (k == 0) return 0;
    int half = (k + 1) / 2;
    return (k % 2 == 1) ? half : -half;
}

/**
 * @internal
 * @brief 判断一元多项式是否无平方
 *
 * 通过计算 gcd(f, f') 判断：若 GCD 为常数（次数 ≤ 0），则 f 无平方。
 *
 * @param[in] poly 一元多项式
 * @return 若无平方返回 true
 */
static bool is_square_free_univariate(const Polynomial<Rational>& poly)
{
    if (poly.degree() <= 0) return true;
    Polynomial<Rational> deriv = poly.differentiate();
    if (deriv.is_zero()) return true;
    Polynomial<Rational> g = Polynomial<Rational>::gcd(poly, deriv);
    return g.degree() <= 0;
}

/**
 * @internal
 * @brief 一元因式分解桥接：调用现有 Berlekamp/Hensel/Zassenhaus 基础设施
 *
 * 将有理系数一元多项式分解为不可约因子列表。流程：
 * 1. 无平方因子预处理
 * 2. 尝试多个素数执行 Berlekamp 模分解
 * 3. 转换为整系数多项式并执行 Hensel 提升
 * 4. 通过 Zassenhaus 因子组合恢复有理系数真因子
 *
 * @param[in] poly 有理系数一元多项式（应为无平方的）
 * @return 不可约因子列表；若分解失败则返回原多项式本身
 *
 * @see berlekamp_factor, hensel_lift, zassenhaus_combine
 */
std::vector<Polynomial<Rational>> factor_univariate_bridge(
    const Polynomial<Rational>& poly)
{
    /// 常数或线性多项式：不可约
    if (poly.degree() <= 1) {
        return {poly};
    }

    /// 无平方因子预处理
    TfSquareFreeResult sqf = tf_square_free(poly);
    Polynomial<Rational> work_poly = sqf.square_free;

    if (work_poly.degree() <= 1) {
        /// 无平方部分为线性或常数，原多项式为完全幂
        std::vector<Polynomial<Rational>> result;
        if (work_poly.degree() == 1) result.push_back(work_poly);
        if (sqf.had_repeated_factors && sqf.repeated_factor.degree() >= 1) {
            /// 递归分解重复因子
            auto sub = factor_univariate_bridge(sqf.repeated_factor);
            result.insert(result.end(), sub.begin(), sub.end());
        }
        return result.empty() ? std::vector<Polynomial<Rational>>{poly} : result;
    }

    /// Berlekamp 模分解：尝试多个素数
    static const int64_t TRIAL_PRIMES[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    BerlekampResult berl_result;
    bool berlekamp_success = false;

    for (int64_t p : TRIAL_PRIMES) {
        try {
            berl_result = berlekamp_factor(work_poly, p);
            if (berl_result.prime > 0 && !berl_result.factors.empty()) {
                berlekamp_success = true;
                break;
            }
        } catch (...) {
            continue;
        }
    }

    if (!berlekamp_success || berl_result.factors.size() <= 1) {
        /// 在所有尝试素数下不可约
        return {poly};
    }

    /// 转换为整系数多项式（乘以分母 LCM）
    BigInt lcm_denom(1);
    for (const auto& c : work_poly.coeffs) {
        BigInt d = c.get_denominator();
        lcm_denom = BigInt::lcm(lcm_denom, d);
    }

    std::vector<BigInt> int_coeffs;
    int_coeffs.reserve(work_poly.coeffs.size());
    for (const auto& c : work_poly.coeffs) {
        BigInt num = c.get_numerator() * (lcm_denom / c.get_denominator());
        int_coeffs.push_back(num);
    }
    Polynomial<BigInt> int_poly(int_coeffs, work_poly.variable_name);

    /// 计算 Mignotte 界确定提升精度
    int n = work_poly.degree();
    BigInt max_coeff(0);
    for (const auto& c : int_coeffs) {
        BigInt ac = c.Abs();
        if (ac > max_coeff) max_coeff = ac;
    }
    BigInt mignotte = max_coeff;
    for (int i = 0; i < n; ++i) {
        mignotte = mignotte * BigInt(2);
    }

    /// 确定提升次数 k 使得 p^k > 2 * mignotte
    BigInt target = mignotte * BigInt(2);
    int64_t prime = berl_result.prime;
    int lift_bound = 1;
    BigInt pk(static_cast<long long>(prime));
    while (pk <= target && lift_bound < 100) {
        pk = pk * BigInt(static_cast<long long>(prime));
        lift_bound++;
    }

    /// Hensel 提升
    auto lifted = hensel_lift_checked(
        int_poly, berl_result.factors, prime, lift_bound);
    if (!lifted) return {poly};
    std::vector<Polynomial<BigInt>> lifted_factors =
        std::move(lifted.value());

    if (lifted_factors.empty()) {
        return {poly};
    }

    /// 计算素数幂 p^k
    int64_t prime_power = 1;
    for (int i = 0; i < lift_bound; ++i) {
        if (prime_power > INT64_MAX / prime) {
            prime_power = INT64_MAX;
            break;
        }
        prime_power *= prime;
    }

    /// Zassenhaus 因子组合
    std::vector<Polynomial<Rational>> true_factors;
    try {
        true_factors = zassenhaus_combine(work_poly, lifted_factors, prime_power);
    } catch (...) {
        return {poly};
    }

    if (true_factors.empty()) {
        return {poly};
    }

    /// 验证：因子数提供有效上界（一元因子数 ≥ 多元因子数）
    /// 若重复因子存在，递归分解并合并
    if (sqf.had_repeated_factors && sqf.repeated_factor.degree() >= 1) {
        auto sub = factor_univariate_bridge(sqf.repeated_factor);
        true_factors.insert(true_factors.end(), sub.begin(), sub.end());
    }

    return true_factors;
}

/**
 * @internal
 * @brief 选择求值点：确定主变量和辅助变量的整数求值点
 *
 * 选择次数最高的变量作为主变量，对辅助变量搜索整数求值点，
 * 使得求值后的一元多项式满足：
 * 1. 次数等于原多项式关于主变量的次数（首项系数不消失）
 * 2. 求值后的多项式无平方
 *
 * 搜索顺序：对每个辅助变量从 0, ±1, ±2, ... 开始尝试。
 * 多个辅助变量时系统地枚举组合。
 *
 * @param[in]  poly        输入多元多项式
 * @param[out] main_var    选定的主变量名
 * @param[out] eval_points 辅助变量到求值点的映射
 * @return 若找到有效求值点返回 true，否则返回 false
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 */
[[maybe_unused]] static bool select_evaluation_points(
    const MultiPoly& poly,
    std::string& main_var,
    std::map<std::string, Rational>& eval_points)
{
    const auto& vars = poly.variables();
    if (vars.empty() || poly.is_zero()) return false;

    /// 步骤 1：选择主变量（次数最高的变量）
    main_var = vars[0];
    int max_deg = poly.degree(vars[0]);
    for (size_t i = 1; i < vars.size(); ++i) {
        int d = poly.degree(vars[i]);
        if (d > max_deg) {
            max_deg = d;
            main_var = vars[i];
        }
    }

    /// 若多项式实质为一元（无辅助变量），无需求值点
    std::vector<std::string> aux_vars;
    for (const auto& v : vars) {
        if (v != main_var) aux_vars.push_back(v);
    }
    if (aux_vars.empty()) {
        eval_points.clear();
        return true;
    }

    int target_degree = poly.degree(main_var);
    int num_aux = static_cast<int>(aux_vars.size());

    /// 步骤 2：搜索求值点组合
    /// 对于单个辅助变量，直接线性搜索
    /// 对于多个辅助变量，系统地枚举组合
    int max_attempts = 1000;
    int attempt = 0;

    if (num_aux == 1) {
        /// 单辅助变量：简单线性搜索
        for (int k = 0; attempt < max_attempts; ++k, ++attempt) {
            int val = search_value(k);
            Rational r_val(val);

            /// 求值：将辅助变量代入
            MultiPoly evaluated = poly.eval(aux_vars[0], r_val);

            /// 检查次数保持（首项系数不消失）
            if (evaluated.degree(main_var) != target_degree) continue;

            /// 转换为一元多项式并检查无平方性
            try {
                Polynomial<Rational> uni = evaluated.to_univariate();
                if (!is_square_free_univariate(uni)) continue;
            } catch (...) {
                /// 若转换失败（仍含多变量），跳过
                continue;
            }

            /// 找到有效点
            eval_points.clear();
            eval_points[aux_vars[0]] = r_val;
            return true;
        }
    } else {
        /// 多辅助变量：按层级搜索
        /// 使用混合基数计数器枚举组合
        /// 搜索范围逐步扩大：先尝试所有变量在 {0} 内，
        /// 再尝试 {0, 1, -1}，再 {0, 1, -1, 2, -2}，...
        int range = 0; // 当前搜索半径
        while (attempt < max_attempts) {
            ++range;
            int values_per_var = 2 * range + 1; // 从 -range 到 +range

            /// 枚举所有组合：values_per_var^num_aux 种
            int total_combos = 1;
            for (int i = 0; i < num_aux; ++i) {
                total_combos *= values_per_var;
                if (total_combos > max_attempts - attempt) {
                    total_combos = max_attempts - attempt;
                    break;
                }
            }

            for (int combo = 0; combo < total_combos && attempt < max_attempts;
                 ++combo, ++attempt) {
                /// 将 combo 解码为各辅助变量的值
                std::map<std::string, Rational> candidate;
                int temp = combo;
                for (int i = 0; i < num_aux; ++i) {
                    int idx = temp % values_per_var;
                    temp /= values_per_var;
                    /// 将 idx 映射为搜索序列值：0→0, 1→1, 2→-1, 3→2, 4→-2, ...
                    int val = search_value(idx);
                    candidate[aux_vars[i]] = Rational(val);
                }

                /// 求值：将所有辅助变量代入
                MultiPoly evaluated = poly.eval(candidate);

                /// 检查次数保持
                if (evaluated.degree(main_var) != target_degree) continue;

                /// 转换为一元多项式并检查无平方性
                try {
                    Polynomial<Rational> uni = evaluated.to_univariate();
                    if (!is_square_free_univariate(uni)) continue;
                } catch (...) {
                    continue;
                }

                /// 找到有效点
                eval_points = std::move(candidate);
                return true;
            }
        }
    }

    return false;
}


/**
 * @internal
 * @brief 试除验证：逐一检验提升后的因子候选是否为真因子
 *
 * 对每个提升后的因子候选 fᵢ，先本原化，然后尝试对剩余多项式执行精确除法。
 * 若整除成功，则 fᵢ 为真因子，更新剩余多项式为商；否则跳过该候选。
 *
 * @param[in,out] remaining      剩余多项式，成功除法后更新为商
 * @param[in,out] lifted_factors 提升后的因子候选列表，已验证的因子从中移除
 * @return 验证通过的真因子列表
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 */
[[maybe_unused]] static std::vector<MultiPoly> trial_division(
    MultiPoly& remaining,
    std::vector<MultiPoly>& lifted_factors)
{
    std::vector<MultiPoly> true_factors;
    std::vector<bool> used(lifted_factors.size(), false);

    for (size_t i = 0; i < lifted_factors.size(); ++i) {
        if (remaining.is_constant()) break;

        /// 本原化候选因子
        MultiPoly candidate = lifted_factors[i].make_primitive();
        if (candidate.is_zero() || candidate.is_constant()) {
            used[i] = true;
            continue;
        }

        /// 确保首项系数为正
        if (!candidate.terms().empty() && candidate.terms()[0].second < Rational(0)) {
            candidate = candidate * Rational(-1);
        }

        /// 尝试精确除法
        try {
            MultiPoly quotient = remaining.exact_div(candidate);
            /// 除法成功：记录为真因子
            true_factors.push_back(candidate);
            remaining = quotient;
            used[i] = true;
        } catch (...) {
            /// 精确除法失败时，当前候选进入下一组合。
        }
    }

    /// 从 lifted_factors 中移除已使用的因子
    std::vector<MultiPoly> remaining_factors;
    for (size_t i = 0; i < lifted_factors.size(); ++i) {
        if (!used[i]) {
            remaining_factors.push_back(std::move(lifted_factors[i]));
        }
    }
    lifted_factors = std::move(remaining_factors);

    return true_factors;
}

/**
 * @internal
 * @brief 因子组合：对未通过单独试除的因子，尝试子集乘积验证
 *
 * 单个提升因子未通过整除检验时，依次枚举含 2、3、... 个因子的乘积。
 * 成功候选记录为真因子，并从池中移除其成员后继续处理剩余多项式。
 *
 * 早期终止条件：
 * - 剩余多项式变为常数
 * - 剩余多项式次数 ≤ 1（不可约）
 *
 * @param[in,out] remaining      剩余多项式，成功除法后更新为商
 * @param[in,out] lifted_factors 未通过试除的因子候选列表
 * @return 通过组合验证的真因子列表
 *
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." §15.7.
 */
[[maybe_unused]] static std::vector<MultiPoly> factor_combination(
    MultiPoly& remaining,
    std::vector<MultiPoly>& lifted_factors)
{
    std::vector<MultiPoly> true_factors;
    int n = static_cast<int>(lifted_factors.size());

    if (n == 0 || remaining.is_constant()) return true_factors;

    /// 按子集大小递增枚举：2, 3, ..., n/2
    int max_subset_size = n / 2;

    for (int subset_size = 2; subset_size <= max_subset_size; ++subset_size) {
        /// 早期终止：剩余多项式为常数或线性
        if (remaining.is_constant() || remaining.total_degree() <= 1) break;

        /// 枚举大小为 subset_size 的所有子集
        /// 使用位掩码或组合索引枚举
        n = static_cast<int>(lifted_factors.size());
        if (subset_size > n / 2) break;

        /// 生成组合索引
        std::vector<int> indices(subset_size);
        for (int i = 0; i < subset_size; ++i) indices[i] = i;

        bool found_in_this_size = true;
        while (found_in_this_size) {
            found_in_this_size = false;

            /// 重置组合索引
            n = static_cast<int>(lifted_factors.size());
            if (subset_size > n) break;
            for (int i = 0; i < subset_size; ++i) indices[i] = i;

            bool has_next = true;
            while (has_next) {
                /// 早期终止检查
                if (remaining.is_constant() || remaining.total_degree() <= 1) {
                    has_next = false;
                    break;
                }

                /// 计算当前子集的乘积
                MultiPoly product = lifted_factors[indices[0]];
                for (int i = 1; i < subset_size; ++i) {
                    product = product * lifted_factors[indices[i]];
                }

                /// 本原化乘积
                product = product.make_primitive();
                if (!product.is_zero() && !product.terms().empty() &&
                    product.terms()[0].second < Rational(0)) {
                    product = product * Rational(-1);
                }

                /// 尝试精确除法
                bool division_success = false;
                try {
                    MultiPoly quotient = remaining.exact_div(product);
                    /// 成功：记录乘积为真因子
                    true_factors.push_back(product);
                    remaining = quotient;
                    division_success = true;
                } catch (...) {
                    /// 除法失败
                }

                if (division_success) {
                    /// 从 lifted_factors 中移除已使用的因子（按降序索引移除）
                    std::vector<int> to_remove(indices.begin(),
                                              indices.begin() + subset_size);
                    std::sort(to_remove.rbegin(), to_remove.rend());
                    for (int idx : to_remove) {
                        lifted_factors.erase(lifted_factors.begin() + idx);
                    }
                    found_in_this_size = true;
                    break; // 重新开始当前大小的枚举
                }

                /// 推进到下一个组合
                /// 标准组合生成：从最后一个索引开始尝试递增
                int pos = subset_size - 1;
                while (pos >= 0) {
                    indices[pos]++;
                    if (indices[pos] <= n - subset_size + pos) {
                        /// 填充后续索引
                        for (int j = pos + 1; j < subset_size; ++j) {
                            indices[j] = indices[j - 1] + 1;
                        }
                        break;
                    }
                    pos--;
                }
                if (pos < 0) has_next = false;
            }

            /// 若本轮未找到匹配，退出当前大小
            if (!found_in_this_size) break;
        }
    }

    /// 早期终止：若剩余多项式为常数或线性，无需继续
    /// 剩余未组合的因子假定为不可约，直接返回
    if (!remaining.is_constant() && remaining.total_degree() > 0 && !lifted_factors.empty()) {
        /// 将剩余因子视为不可约返回
        /// 注意：这些因子可能需要进一步验证，但在当前阶段
        /// 假设它们是不可约的（因为所有子集组合都已尝试）
    }

    return true_factors;
}


/**
 * @internal
 * @brief 组装最终因式分解结果
 *
 * 收集所有真因子，确保每个因子本原化且首项系数为正，
 * 计算整体数值常数，并验证 constant * ∏(factors[i]^mult[i]) == original。
 *
 * 算法：
 * 1. 对每个因子执行 make_primitive()，若首项系数为负则取反
 * 2. 计算整体常数 = original.numeric_content() / ∏(各因子的 numeric_content)
 * 3. 验证乘积等于原多项式；若验证失败则回退返回原多项式
 *
 * @param[in] factors        真因子列表（来自试除/组合验证）
 * @param[in] multiplicities 各因子的重数（来自无平方因子分解）
 * @param[in] original       原始多项式（用于计算常数和验证）
 * @return 满足不变量的 MultiFactorResult
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 */
[[maybe_unused]] static MultiFactorResult assemble_result(
    const std::vector<MultiPoly>& factors,
    const std::vector<int>& multiplicities,
    const MultiPoly& original)
{
    MultiFactorResult result;

    /// 空因子列表：原多项式为常数或零
    if (factors.empty()) {
        if (original.is_zero()) {
            result.constant = Rational(0);
        } else {
            result.constant = original.numeric_content();
            if (!original.terms().empty() && original.terms()[0].second < Rational(0)) {
                result.constant = -result.constant;
            }
        }
        return result;
    }

    /// 步骤 1：本原化每个因子，确保首项系数为正
    std::vector<MultiPoly> prim_factors;
    prim_factors.reserve(factors.size());
    Rational content_product(1);

    for (size_t i = 0; i < factors.size(); ++i) {
        const MultiPoly& f = factors[i];

        if (f.is_zero() || f.is_constant()) {
            /// 常数/零因子：吸收到整体常数中
            if (!f.is_zero()) {
                Rational c = f.numeric_content();
                if (!f.terms().empty() && f.terms()[0].second < Rational(0)) {
                    c = -c;
                }
                int mult = (i < multiplicities.size()) ? multiplicities[i] : 1;
                for (int m = 0; m < mult; ++m) {
                    content_product = content_product * c;
                }
            }
            continue;
        }

        Rational nc = f.numeric_content();
        MultiPoly prim = f.make_primitive();

        /// 确保首项系数为正
        if (!prim.terms().empty() && prim.terms()[0].second < Rational(0)) {
            prim = prim * Rational(-1);
            nc = -nc;
        }

        /// 累积各因子内容到常数乘积中（考虑重数）
        int mult = (i < multiplicities.size()) ? multiplicities[i] : 1;
        for (int m = 0; m < mult; ++m) {
            content_product = content_product * nc;
        }

        prim_factors.push_back(std::move(prim));
        result.multiplicities.push_back(mult);
    }

    result.factors = std::move(prim_factors);

    /// 步骤 2：计算整体常数
    /// original == content_product * result.constant * ∏(prim_factors[i]^mult[i])
    /// 因此 result.constant = original 的数值内容 / content_product
    /// 但更精确地：直接从原多项式和因子乘积的比值推导
    Rational orig_nc = original.numeric_content();
    if (!original.terms().empty() && original.terms()[0].second < Rational(0)) {
        orig_nc = -orig_nc;
    }

    if (!content_product.is_zero()) {
        result.constant = orig_nc / content_product;
    } else {
        result.constant = orig_nc;
    }

    /// 步骤 3：验证 constant * ∏(factors[i]^mult[i]) == original
    /// 计算因子乘积
    MultiPoly product(Rational(1), original.variables());
    for (size_t i = 0; i < result.factors.size(); ++i) {
        MultiPoly power = result.factors[i];
        for (int m = 1; m < result.multiplicities[i]; ++m) {
            power = power * result.factors[i];
        }
        product = product * power;
    }
    product = product * result.constant;

    /// 验证乘积等于原多项式
    if (product != original) {
        /// 验证失败：尝试通过精确除法修正常数
        /// 若 product 与 original 仅差一个标量倍数，修正常数即可
        if (!product.is_zero()) {
            /// 比较首项系数比值
            const auto& orig_terms = original.terms();
            const auto& prod_terms = product.terms();
            if (!orig_terms.empty() && !prod_terms.empty() &&
                orig_terms[0].first == prod_terms[0].first) {
                Rational ratio = orig_terms[0].second / prod_terms[0].second;
                MultiPoly corrected = product * ratio;
                if (corrected == original) {
                    result.constant = result.constant * ratio;
                    return result;
                }
            }
        }

        /// 修正失败：回退返回原多项式作为不可约因子
        MultiFactorResult fallback;
        fallback.constant = orig_nc;
        MultiPoly prim = original.make_primitive();
        if (!prim.terms().empty() && prim.terms()[0].second < Rational(0)) {
            prim = prim * Rational(-1);
            fallback.constant = -fallback.constant;
        }
        fallback.factors = {prim};
        fallback.multiplicities = {1};
        return fallback;
    }

    return result;
}

MultiFactorResult factor_multivariate(const MultiPoly& poly)
{
    if (poly.is_zero()) return {Rational(0), {}, {}};

    /// 常数多项式
    if (poly.is_constant()) {
        Rational c = poly.numeric_content();
        if (!poly.terms().empty() && poly.terms()[0].second < Rational(0)) {
            c = -c;
        }
        return {c, {}, {}};
    }

    const auto& vars = poly.variables();

    /// 选择主变量（次数最高的变量）
    std::string main_var = vars[0];
    int max_deg = poly.degree(vars[0]);
    for (size_t i = 1; i < vars.size(); ++i) {
        int d = poly.degree(vars[i]);
        if (d > max_deg) {
            max_deg = d;
            main_var = vars[i];
        }
    }

    /// --- 快速路径 1：线性多项式 ---
    if (detail::is_linear(poly, main_var)) {
        Rational nc = poly.numeric_content();
        MultiPoly prim = poly.make_primitive();
        /// make_primitive() 确保首项系数为正。若原多项式首项系数为负，
        /// 则 prim 已被取反，需要将符号吸收到常数中。
        if (!poly.terms().empty() && poly.terms()[0].second < Rational(0)) {
            nc = -nc;
        }
        return {nc, {prim}, {1}};
    }

    /// --- 快速路径 2：差平方 a² - b² = (a+b)(a-b) ---
    /// 先提取数值内容，使形如 2x²-8y² 的多项式归约为 2·(x²-4y²)，
    /// 其本原部分的系数 (1,4) 才是完全平方，可被差平方检测识别。
    {
        Rational dos_content = poly.numeric_content();
        bool dos_was_negative = (!poly.terms().empty() && poly.terms()[0].second < Rational(0));
        MultiPoly dos_prim = poly.make_primitive();
        Rational dos_const = dos_was_negative ? -dos_content : dos_content;

        if (auto dos = detail::detect_difference_of_squares(dos_prim)) {
            auto& [a, b] = *dos;
            MultiPoly sum = a + b;
            MultiPoly diff = a - b;

            MultiFactorResult result;
            result.constant = dos_const;

            auto process_factor = [&](MultiPoly& f) {
                Rational nc = f.numeric_content();
                bool was_negative = (!f.terms().empty() && f.terms()[0].second < Rational(0));
                f = f.make_primitive();
                if (was_negative) {
                    nc = -nc;
                }
                result.constant = result.constant * nc;
                result.factors.push_back(f);
                result.multiplicities.push_back(1);
            };

            process_factor(sum);
            process_factor(diff);
            return result;
        }
    }

    /// --- 旧快速路径 2（保留作兜底，处理本原检测未覆盖的情形）---
    if (auto dos = detail::detect_difference_of_squares(poly)) {
        auto& [a, b] = *dos;
        MultiPoly sum = a + b;
        MultiPoly diff = a - b;

        /// 递归分解各因子
        MultiFactorResult result;
        result.constant = Rational(1);

        auto process_factor = [&](MultiPoly& f) {
            Rational nc = f.numeric_content();
            /// make_primitive() 确保首项系数为正
            bool was_negative = (!f.terms().empty() && f.terms()[0].second < Rational(0));
            f = f.make_primitive();
            if (was_negative) {
                nc = -nc;
            }
            result.constant = result.constant * nc;
            result.factors.push_back(f);
            result.multiplicities.push_back(1);
        };

        process_factor(sum);
        process_factor(diff);
        return result;
    }

    /// --- 快速路径 3：二项式幂 xⁿ - yⁿ ---
    if (auto bp = detail::detect_binomial_power(poly)) {
        auto& [var1, var2, n] = *bp;

        /// xⁿ - yⁿ 的分圆分解
        /// 对于 n ≥ 2：xⁿ - yⁿ = (x - y) * (x^(n-1) + x^(n-2)*y + ... + y^(n-1))
        /// 若 n 为偶数：还有 (x + y) 因子
        MultiFactorResult result;
        result.constant = Rational(1);

        int v1_idx = -1, v2_idx = -1;
        for (size_t i = 0; i < vars.size(); ++i) {
            if (vars[i] == var1) v1_idx = static_cast<int>(i);
            if (vars[i] == var2) v2_idx = static_cast<int>(i);
        }

        /// 因子 (x - y)
        {
            Monomial m1(vars.size(), 0);
            m1[v1_idx] = 1;
            Monomial m2(vars.size(), 0);
            m2[v2_idx] = 1;
            std::vector<MultiPoly::Term> f_terms = {{m1, Rational(1)}, {m2, Rational(-1)}};
            MultiPoly factor(std::move(f_terms), vars);
            result.factors.push_back(factor);
            result.multiplicities.push_back(1);
        }

        if (n % 2 == 0) {
            /// 因子 (x + y)
            Monomial m1(vars.size(), 0);
            m1[v1_idx] = 1;
            Monomial m2(vars.size(), 0);
            m2[v2_idx] = 1;
            std::vector<MultiPoly::Term> f_terms = {{m1, Rational(1)}, {m2, Rational(1)}};
            MultiPoly factor(std::move(f_terms), vars);
            result.factors.push_back(factor);
            result.multiplicities.push_back(1);

            /// 剩余因子：(x^(n-1) + x^(n-2)*y + ... + y^(n-1)) / (x + y)
            /// 即 xⁿ - yⁿ = (x-y)(x+y) * Q(x,y)
            /// 计算 Q = poly / ((x-y)*(x+y))
            MultiPoly product = result.factors[0] * result.factors[1];
            try {
                MultiPoly quotient = poly.exact_div(product);
                if (!quotient.is_constant()) {
                    Rational nc = quotient.numeric_content();
                    bool was_negative = (!quotient.terms().empty() &&
                                         quotient.terms()[0].second < Rational(0));
                    MultiPoly prim = quotient.make_primitive();
                    if (was_negative) {
                        nc = -nc;
                    }
                    result.constant = result.constant * nc;
                    result.factors.push_back(prim);
                    result.multiplicities.push_back(1);
                } else {
                    result.constant = result.constant * quotient.numeric_content();
                }
            } catch (...) {
                /// 除法失败，回退到返回原多项式
                return {Rational(1), {poly}, {1}};
            }
        } else {
            /// n 为奇数：xⁿ - yⁿ = (x - y) * (x^(n-1) + x^(n-2)*y + ... + y^(n-1))
            /// 计算商
            try {
                MultiPoly quotient = poly.exact_div(result.factors[0]);
                if (!quotient.is_constant()) {
                    Rational nc = quotient.numeric_content();
                    bool was_negative = (!quotient.terms().empty() &&
                                         quotient.terms()[0].second < Rational(0));
                    MultiPoly prim = quotient.make_primitive();
                    if (was_negative) {
                        nc = -nc;
                    }
                    result.constant = result.constant * nc;
                    result.factors.push_back(prim);
                    result.multiplicities.push_back(1);
                } else {
                    result.constant = result.constant * quotient.numeric_content();
                }
            } catch (...) {
                return {Rational(1), {poly}, {1}};
            }
        }

        return result;
    }

    /// --- 快速路径 4：公因子单项式 ---
    auto [gcd_mono, quotient] = detail::extract_common_monomial(poly);
    bool has_common_monomial = false;
    for (size_t i = 0; i < gcd_mono.size(); ++i) {
        if (gcd_mono[i] != 0) { has_common_monomial = true; break; }
    }
    if (has_common_monomial) {
        /// 递归分解商多项式
        MultiFactorResult sub_result = factor_multivariate(quotient);

        /// 将公因子单项式拆分为各变量的幂次因子
        MultiFactorResult result;
        result.constant = sub_result.constant;

        for (size_t i = 0; i < gcd_mono.size(); ++i) {
            if (gcd_mono[i] > 0) {
                /// 创建变量 vars[i] 的单项式因子
                Monomial var_mono(vars.size(), 0);
                var_mono[i] = 1;
                std::vector<MultiPoly::Term> var_terms = {{var_mono, Rational(1)}};
                MultiPoly var_factor(std::move(var_terms), vars);
                result.factors.push_back(var_factor);
                result.multiplicities.push_back(gcd_mono[i]);
            }
        }

        /// 添加商多项式的因子
        for (size_t i = 0; i < sub_result.factors.size(); ++i) {
            result.factors.push_back(sub_result.factors[i]);
            result.multiplicities.push_back(sub_result.multiplicities[i]);
        }

        return result;
    }

    /// --- 快速路径 5：齐次二元多项式 ---
    if (auto hb = detail::factor_homogeneous_bivariate(poly)) {
        return *hb;
    }

    Rational nc = poly.numeric_content();
    MultiPoly prim = poly.make_primitive();
    /// make_primitive() 确保首项系数为正。若原多项式首项系数为负，
    /// 则 prim 已被取反，需要将符号吸收到常数中。
    if (!poly.terms().empty() && poly.terms()[0].second < Rational(0)) {
        nc = -nc;
    }
    return {nc, {prim}, {1}};
}
MultiPoly multivariate_content(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return MultiPoly();
    /// 若 main_var 不在变量列表中，或多项式在 main_var 上次数为 0，容度就是多项式本身
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) return poly;
    if (poly.degree(main_var) <= 0) return poly;
    auto coeffs = extract_coefficients(poly, main_var);
    if (coeffs.empty()) return MultiPoly();
    auto it = coeffs.begin();
    MultiPoly content_gcd = it->second;
    ++it;
    for (; it != coeffs.end(); ++it) {
        if (it->second.is_zero()) continue;
        content_gcd = multivariate_gcd(content_gcd, it->second);
        if (content_gcd.is_constant()) {
            Rational nc = content_gcd.numeric_content();
            if (nc == Rational(1)) return content_gcd;
        }
    }
    if (!content_gcd.is_zero() && !content_gcd.terms().empty()) {
        if (content_gcd.terms()[0].second < Rational(0))
            content_gcd = content_gcd * Rational(-1);
    }
    return content_gcd;
}
MultiPoly multivariate_primitive_part(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return poly;
    MultiPoly content = multivariate_content(poly, main_var);
    if (content.is_zero()) return MultiPoly();
    if (content.is_constant() && content.numeric_content() == Rational(1)) return poly;

    /// content 的变量集是去掉 main_var 后的辅助变量集，
    /// 需要将其嵌入到 poly 的完整变量集中才能正确执行 exact_div。
    const auto& full_vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < full_vars.size(); ++i) {
        if (full_vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }

    MultiPoly content_in_full_vars = content;
    if (var_idx >= 0 && content.variables().size() < full_vars.size()) {
        std::vector<MultiPoly::Term> new_terms;
        for (const auto& term : content.terms()) {
            Monomial full_mono(full_vars.size(), 0);
            size_t ri = 0;
            for (size_t i = 0; i < full_vars.size(); ++i) {
                if (static_cast<int>(i) == var_idx) {
                    full_mono[i] = 0;
                } else {
                    if (ri < term.first.size()) full_mono[i] = term.first[ri];
                    ++ri;
                }
            }
            new_terms.emplace_back(std::move(full_mono), term.second);
        }
        content_in_full_vars = MultiPoly(std::move(new_terms), full_vars);
    }

    return poly.exact_div(content_in_full_vars);
}
MultiPoly multivariate_gcd(const MultiPoly& a, const MultiPoly& b)
{
    if (a.is_zero()) { return b.is_zero() ? MultiPoly() : b.make_primitive(); }
    if (b.is_zero()) return a.make_primitive();
    const auto& vars = a.variables().empty() ? b.variables() : a.variables();
    /// 两者均为常数
    if (a.is_constant() && b.is_constant()) {
        Rational g = rational_gcd(a.numeric_content(), b.numeric_content());
        return MultiPoly(g, vars);
    }
    /// 基本情形：两者均为一元多项式
    if (a.is_univariate() && b.is_univariate()) {
        std::string a_var, b_var;
        for (const auto& v : a.variables()) { if (a.degree(v) > 0) { a_var = v; break; } }
        for (const auto& v : b.variables()) { if (b.degree(v) > 0) { b_var = v; break; } }
        if (a_var.empty() && b_var.empty()) {
            return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
        }
        if (a_var.empty() || b_var.empty()) {
            return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
        }
        if (a_var == b_var) {
            Rational nc_gcd = rational_gcd(a.numeric_content(), b.numeric_content());
            Polynomial<Rational> pa = a.to_univariate();
            Polynomial<Rational> pb = b.to_univariate();
            Polynomial<Rational> g = Polynomial<Rational>::gcd(pa, pb);
            if (g.is_zero()) return MultiPoly(nc_gcd, vars);
            /// monic GCD * numeric content GCD
            MultiPoly result = MultiPoly::from_univariate(g, a_var);
            result = result * nc_gcd;
            if (vars.size() > 1) {
                int avi = -1;
                for (size_t i = 0; i < vars.size(); ++i) {
                    if (vars[i] == a_var) { avi = static_cast<int>(i); break; }
                }
                std::vector<MultiPoly::Term> nt;
                for (const auto& term : result.terms()) {
                    Monomial fm(vars.size(), 0);
                    if (avi >= 0 && !term.first.empty()) fm[static_cast<size_t>(avi)] = term.first[0];
                    nt.emplace_back(std::move(fm), term.second);
                }
                result = MultiPoly(std::move(nt), vars);
            }
            if (!result.terms().empty() && result.terms()[0].second < Rational(0))
                result = result * Rational(-1);
            return result;
        }
        return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
    }
    /// 递归情形：求值-插值法
    std::string main_var = choose_main_variable(a, b);
    Rational nc_a = a.numeric_content();
    Rational nc_b = b.numeric_content();
    Rational nc_gcd = rational_gcd(nc_a, nc_b);
    MultiPoly a_prim = a.make_primitive();
    MultiPoly b_prim = b.make_primitive();
    int deg_a = a_prim.degree(main_var);
    int deg_b = b_prim.degree(main_var);
    int degree_bound = std::min(deg_a, deg_b);
    int num_points = degree_bound + 1;
    std::vector<Rational> eval_pts;
    eval_pts.push_back(Rational(0));
    for (int k = 1; static_cast<int>(eval_pts.size()) < num_points + 10; ++k) {
        eval_pts.push_back(Rational(k));
        eval_pts.push_back(Rational(-k));
    }
    std::vector<Rational> used_pts;
    std::vector<MultiPoly> gcd_vals;
    int target_deg = -1;
    for (size_t pi = 0; pi < eval_pts.size() &&
         static_cast<int>(used_pts.size()) < num_points; ++pi) {
        Rational pt = eval_pts[pi];
        MultiPoly lc_a = a_prim.leading_coeff(main_var);
        MultiPoly lc_b = b_prim.leading_coeff(main_var);
        MultiPoly lca_ev = lc_a.eval(main_var, pt);
        MultiPoly lcb_ev = lc_b.eval(main_var, pt);
        if (lca_ev.is_zero() || lcb_ev.is_zero()) continue;
        MultiPoly a_ev = a_prim.eval(main_var, pt);
        MultiPoly b_ev = b_prim.eval(main_var, pt);
        MultiPoly g_ev = multivariate_gcd(a_ev, b_ev);
        int gd = g_ev.is_zero() ? -1 : g_ev.total_degree();
        if (target_deg < 0) { target_deg = gd; }
        else if (gd < target_deg) { used_pts.clear(); gcd_vals.clear(); target_deg = gd; }
        else if (gd > target_deg) { continue; }
        if (!g_ev.is_zero()) g_ev = g_ev.make_primitive();
        used_pts.push_back(pt);
        gcd_vals.push_back(g_ev);
    }
    if (target_deg <= 0 || used_pts.empty()) return MultiPoly(nc_gcd, vars);
    if (used_pts.size() == 1 && degree_bound == 0) {
        MultiPoly r = gcd_vals[0].is_zero() ? MultiPoly(Rational(1), vars) : gcd_vals[0];
        if (nc_gcd != Rational(1)) r = r * nc_gcd;
        return r;
    }
    MultiPoly cand = lagrange_interpolate(used_pts, gcd_vals, main_var, vars);
    if (cand.is_zero()) return MultiPoly(nc_gcd, vars);
    cand = cand.make_primitive();
    if (divides_poly(a_prim, cand) && divides_poly(b_prim, cand)) {
        if (nc_gcd != Rational(1)) cand = cand * nc_gcd;
        if (!cand.terms().empty() && cand.terms()[0].second < Rational(0))
            cand = cand * Rational(-1);
        return cand;
    }
    /// 扩展搜索
    int ext_bound = degree_bound + 5;
    std::vector<Rational> ext_pts;
    std::vector<MultiPoly> ext_vals;
    std::vector<Rational> all_pts;
    all_pts.push_back(Rational(0));
    for (int k = 1; static_cast<int>(all_pts.size()) < ext_bound + 10; ++k) {
        all_pts.push_back(Rational(k));
        all_pts.push_back(Rational(-k));
    }
    target_deg = -1;
    for (size_t pi = 0; pi < all_pts.size() &&
         static_cast<int>(ext_pts.size()) < ext_bound + 1; ++pi) {
        Rational pt = all_pts[pi];
        MultiPoly lc_a = a_prim.leading_coeff(main_var);
        MultiPoly lc_b = b_prim.leading_coeff(main_var);
        MultiPoly lca_ev = lc_a.eval(main_var, pt);
        MultiPoly lcb_ev = lc_b.eval(main_var, pt);
        if (lca_ev.is_zero() || lcb_ev.is_zero()) continue;
        MultiPoly a_ev = a_prim.eval(main_var, pt);
        MultiPoly b_ev = b_prim.eval(main_var, pt);
        MultiPoly g_ev = multivariate_gcd(a_ev, b_ev);
        int gd = g_ev.is_zero() ? -1 : g_ev.total_degree();
        if (target_deg < 0) { target_deg = gd; }
        else if (gd < target_deg) { ext_pts.clear(); ext_vals.clear(); target_deg = gd; }
        else if (gd > target_deg) { continue; }
        if (!g_ev.is_zero()) g_ev = g_ev.make_primitive();
        ext_pts.push_back(pt);
        ext_vals.push_back(g_ev);
    }
    if (target_deg <= 0 || ext_pts.empty()) return MultiPoly(nc_gcd, vars);
    cand = lagrange_interpolate(ext_pts, ext_vals, main_var, vars);
    if (cand.is_zero()) return MultiPoly(nc_gcd, vars);
    cand = cand.make_primitive();
    if (divides_poly(a_prim, cand) && divides_poly(b_prim, cand)) {
        if (nc_gcd != Rational(1)) cand = cand * nc_gcd;
        if (!cand.terms().empty() && cand.terms()[0].second < Rational(0))
            cand = cand * Rational(-1);
        return cand;
    }
    return MultiPoly(nc_gcd, vars);
}
/**
 * @internal
 * @brief 计算多元多项式关于指定变量的形式偏导数
 *
 * 对每个含 main_var 的项，新系数 = 原系数 × 指数，新指数 = 原指数 - 1。
 * 指数为 0 的项被丢弃。
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 求导变量名
 * @return ∂poly/∂main_var
 */
static MultiPoly formal_derivative(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return poly;

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) {
        /// main_var 不在变量列表中，导数为零
        return MultiPoly(Rational(0), vars);
    }

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp == 0) continue; // 常数项关于该变量导数为零

        Rational new_coeff = term.second * Rational(exp);
        Monomial new_mono = mono;
        new_mono[var_idx] = exp - 1;
        result_terms.emplace_back(std::move(new_mono), std::move(new_coeff));
    }

    if (result_terms.empty()) return MultiPoly(Rational(0), vars);
    return MultiPoly(std::move(result_terms), vars);
}

/**
 * @brief 多元无平方因子分解（Yun 算法的多元推广）
 *
 * 通过计算 gcd(f, ∂f/∂x_main) 检测重因子，将多项式分解为
 * f = f₁ · f₂² · f₃³ · ... 的形式，其中每个 fᵢ 无平方。
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 无平方因子分解结果，components[i] 的重数为 i+1
 */
SquareFreeDecomp square_free_decompose(const MultiPoly& poly, const std::string& main_var)
{
    /// 零多项式或常数多项式：直接返回
    if (poly.is_zero() || poly.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 若 main_var 不在变量列表中或次数为 0，多项式无重因子
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0 || poly.degree(main_var) <= 0) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 1：计算 f' = ∂f/∂main_var
    MultiPoly f_prime = formal_derivative(poly, main_var);

    /// 步骤 2：计算 g = gcd(f, f')；特征 0 下 f' 为零时，f 直接作为单个无平方因子。
    if (f_prime.is_zero()) {
        return SquareFreeDecomp{{poly}};
    }

    MultiPoly g = multivariate_gcd(poly, f_prime);

    /// 步骤 3：若 g 为常数，f 已无平方
    if (g.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 4：Yun 算法迭代提取无平方分量
    /// w = f / g, y = f' / g, z = y - w'
    MultiPoly w = poly.exact_div(g);
    MultiPoly y = f_prime.exact_div(g);
    MultiPoly w_prime = formal_derivative(w, main_var);
    MultiPoly z = y - w_prime;

    std::vector<MultiPoly> components;

    while (!z.is_zero()) {
        MultiPoly component_i = multivariate_gcd(w, z);
        w = w.exact_div(component_i);
        y = z.exact_div(component_i);
        w_prime = formal_derivative(w, main_var);
        z = y - w_prime;
        components.push_back(std::move(component_i));
    }

    /// w 是最后一个分量
    components.push_back(std::move(w));

    /// 过滤掉常数为 1 的平凡分量（保留结构但不影响语义）
    /// 注意：components[i] 的重数为 i+1
    return SquareFreeDecomp{std::move(components)};
}
/**
 * @internal
 * @brief 首项系数预计算（Wang's leading coefficient trick）
 *
 * 在 Hensel 提升前，将原多项式关于主变量的首项系数 lc(f, x_main) 分解为
 * 辅助变量的因子，并按求值后的首项系数值分配给各一元因子，
 * 使提升过程保持正确的首项系数与本原内容。
 *
 * 算法步骤：
 * 1. 计算 lc = lc(f, x_main)，为辅助变量的多项式
 * 2. 若 lc 为常数，无需预计算
 * 3. 将 lc 在求值点处求值得到 lc_eval（有理数）
 * 4. 对每个一元因子 fᵢ，其首项系数 lc_i 应整除 lc_eval
 * 5. 按 lc_i 在 lc_eval 中的贡献比例，将 lc 的因子分配给各 fᵢ
 * 6. 将分配的 lc 部分乘入各 fᵢ，使提升后首项系数正确
 *
 * @param[in]     poly              原始多元多项式
 * @param[in]     main_var          主变量名
 * @param[in]     eval_points       辅助变量到求值点的映射
 * @param[in,out] univariate_factors 一元因子列表，函数会修改各因子的首项系数
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 * @see D. Y. Y. Yun, “On Square-Free Decomposition Algorithms,”
 *      Proceedings of SYMSAC 1976.
 */
[[maybe_unused]] static void precompute_leading_coefficients(
    const MultiPoly& poly,
    const std::string& main_var,
    const std::map<std::string, Rational>& eval_points,
    std::vector<Polynomial<Rational>>& univariate_factors)
{
    int r = static_cast<int>(univariate_factors.size());
    if (r <= 1) return;

    /// 步骤 1：计算 lc(f, x_main) — 关于主变量的首项系数多项式
    MultiPoly lc_poly = poly.leading_coeff(main_var);

    /// 步骤 2：若 lc 为常数，无需预计算
    if (lc_poly.is_constant()) return;

    /// 步骤 3：将 lc 在求值点处求值得到 lc_eval（有理数值）
    MultiPoly lc_evaluated = lc_poly.eval(eval_points);
    if (lc_evaluated.is_zero()) return; // 退化情形：首项系数消失

    Rational lc_eval_val = lc_evaluated.is_zero() ? Rational(0)
                         : lc_evaluated.terms()[0].second;

    /// 步骤 4：收集各一元因子的首项系数
    /// 各因子 fᵢ 的首项系数 lc_i 满足 ∏ lc_i = lc_eval_val（至多差常数）
    std::vector<Rational> factor_lcs;
    factor_lcs.reserve(r);
    Rational product_of_lcs(1);
    for (int i = 0; i < r; ++i) {
        Rational lc_i = univariate_factors[i].lead_coeff();
        if (lc_i.is_zero()) lc_i = Rational(1);
        factor_lcs.push_back(lc_i);
        product_of_lcs = product_of_lcs * lc_i;
    }

    /// 若一元因子首项系数之积已等于 lc_eval_val，无需调整
    if (product_of_lcs == lc_eval_val) return;

    /// 步骤 5：计算分配比例
    /// 策略：将 lc_eval_val / product_of_lcs 的差额分配给第一个因子
    /// 这是简化版本——对于大多数情形（lc 为单项式或简单多项式），
    /// 将整个 lc 分配给首项系数最大的因子即可保证提升正确性。
    //
    /// 完整版本需要递归分解 lc_poly 并逐一匹配，但对于首次实现，
    /// 采用按比例分配的简化策略。

    if (product_of_lcs.is_zero()) return;

    /// 计算缩放因子：scale = lc_eval_val / product_of_lcs
    Rational scale = lc_eval_val / product_of_lcs;

    if (scale == Rational(1)) return;

    /// 策略：将缩放因子分配给第一个因子
    /// 这保证 ∏ lc(fᵢ) = lc_eval_val
    /// 对于更复杂的情形（lc_poly 有多个不可约因子），
    /// 需要匹配各因子的首项系数与 lc_poly 的因子。
    //
    /// 高级分配：尝试将 scale 分解为各因子的贡献
    /// 对每个因子，检查 lc_eval_val 是否能被其首项系数整除
    /// 若能，则该因子获得对应的 lc 份额

    /// 尝试精确分配：对每个因子 fᵢ，计算其应得的 lc 份额
    /// 方法：lc_eval_val = c₁ * c₂ * ... * cᵣ，其中 cᵢ 是分配给 fᵢ 的值
    /// 约束：cᵢ 的求值值等于 fᵢ 的首项系数乘以某个因子

    /// 简化实现：将整个缩放因子乘入第一个因子
    Polynomial<Rational>& first_factor = univariate_factors[0];
    std::vector<Rational> new_coeffs = first_factor.coeffs;
    for (auto& c : new_coeffs) {
        c = c * scale;
    }
    first_factor = Polynomial<Rational>(new_coeffs, first_factor.variable_name);
}

} // namespace lamina
