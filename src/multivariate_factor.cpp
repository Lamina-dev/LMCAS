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
namespace lamina {
// ============================================================
/// 前向声明
// ============================================================
static MultiPoly truncate_mod_var(const MultiPoly& poly, const std::string& var,
                                  int degree_bound);

// ============================================================
/// 内部辅助函数
// ============================================================
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
// ============================================================
/// 求值点选择
// ============================================================

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
static std::vector<Polynomial<Rational>> factor_univariate_bridge(
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
    std::vector<Polynomial<BigInt>> lifted_factors;
    try {
        lifted_factors = hensel_lift(int_poly, berl_result.factors, prime, lift_bound);
    } catch (...) {
        return {poly};
    }

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

// ============================================================
/// 试除验证与因子组合
// ============================================================

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
            /// 除法失败：该候选不是真因子，跳过
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
 * 当单个提升因子不能整除剩余多项式时，尝试 2 个、3 个、... 因子的乘积。
 * 按子集大小递增枚举，对每个子集计算乘积并本原化后尝试精确除法。
 * 成功时记录乘积为真因子，从池中移除已用因子，更新剩余多项式。
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

// ============================================================
/// 结果组装
// ============================================================

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

// ============================================================
/// 公共 API
// ============================================================
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

    /// --- 回退：返回原多项式作为不可约因子 ---
    /// （完整 Wang-EEZ 流程将在 task 10 中实现）
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

    /// 步骤 2：计算 g = gcd(f, f')
    /// 若 f' 为零（特征 0 下不应发生，但防御性处理），f 本身无平方
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
 * 辅助变量的因子，并按求值后的首项系数值分配给各一元因子。
 * 这确保提升过程产生正确的首项系数，避免引入虚假内容。
 *
 * 算法（参考 Wang, P.S. 1978）：
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

/**
 * @brief 多元 Hensel 提升（单变量步）
 *
 * 将一元因子逐变量提升，从 f(x₁, a₂, ..., aₙ) 的因子恢复到包含 lift_var 的因子。
 * 每次引入一个辅助变量，计算修正项直到达到次数上界。
 *
 * 算法（参考 Geddes, Czapor, Labahn §15.6）：
 * 1. 将一元因子嵌入到完整变量集中
 * 2. 对 k = 1, ..., degree_bound：
 *    a. 计算当前因子乘积（截断到 lift_var 次数 ≤ k）
 *    b. 计算误差 E = poly - product（提取 lift_var 次数恰为 k 的部分）
 *    c. 若误差为零则继续
 *    d. 用丢番图求解器分配误差到各因子
 *    e. 将修正项加到各因子上
 *    f. 验证：product of lifted factors ≡ poly mod (lift_var - eval_point)^(k+1)
 * 3. 返回提升后的因子
 *
 * @param[in] poly              原始多元多项式
 * @param[in] univariate_factors 一元分解得到的因子列表
 * @param[in] lift_var          当前提升的变量名
 * @param[in] eval_point        该变量的求值点
 * @param[in] degree_bound      提升次数上界
 * @return 提升后的多元因子列表
 *
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." §15.6.
 */
std::vector<MultiPoly> multivariate_hensel_lift(
    const MultiPoly& poly,
    const std::vector<Polynomial<Rational>>& univariate_factors,
    const std::string& lift_var,
    const Rational& eval_point,
    int degree_bound)
{
    int r = static_cast<int>(univariate_factors.size());
    if (r == 0) return {};
    if (r == 1) {
        return {poly};
    }
    if (degree_bound <= 0) {
        /// 无需提升：将一元因子嵌入完整变量集返回
        const auto& vars = poly.variables();
        std::vector<MultiPoly> result;
        result.reserve(r);
        for (const auto& uf : univariate_factors) {
            MultiPoly mp = MultiPoly::from_univariate(uf, uf.variable_name);
            if (mp.variables() != vars && !vars.empty()) {
                std::string uni_var = uf.variable_name;
                int uni_idx = -1;
                for (size_t i = 0; i < vars.size(); ++i) {
                    if (vars[i] == uni_var) { uni_idx = static_cast<int>(i); break; }
                }
                if (uni_idx >= 0) {
                    std::vector<MultiPoly::Term> new_terms;
                    for (const auto& term : mp.terms()) {
                        Monomial full_mono(vars.size(), 0);
                        if (!term.first.empty()) full_mono[uni_idx] = term.first[0];
                        new_terms.emplace_back(std::move(full_mono), term.second);
                    }
                    mp = MultiPoly(std::move(new_terms), vars);
                }
            }
            result.push_back(std::move(mp));
        }
        return result;
    }

    const auto& vars = poly.variables();

    /// 确定主变量（一元因子的变量）和提升变量在完整变量集中的位置
    std::string main_var = univariate_factors[0].variable_name;
    int main_var_idx = -1;
    int lift_var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) main_var_idx = static_cast<int>(i);
        if (vars[i] == lift_var) lift_var_idx = static_cast<int>(i);
    }

    /// 将一元因子嵌入到完整变量集中
    std::vector<MultiPoly> factors;
    factors.reserve(r);
    for (const auto& uf : univariate_factors) {
        MultiPoly mp = MultiPoly::from_univariate(uf, main_var);
        if (mp.variables().size() != vars.size() && !vars.empty()) {
            std::vector<MultiPoly::Term> new_terms;
            for (const auto& term : mp.terms()) {
                Monomial full_mono(vars.size(), 0);
                if (main_var_idx >= 0 && !term.first.empty()) {
                    full_mono[main_var_idx] = term.first[0];
                }
                new_terms.emplace_back(std::move(full_mono), term.second);
            }
            mp = MultiPoly(std::move(new_terms), vars);
        }
        factors.push_back(std::move(mp));
    }

    /// 若 lift_var 不在变量集中，无法提升，直接返回嵌入后的因子
    if (lift_var_idx < 0) {
        return factors;
    }

    /// 迭代提升：对 k = 1, ..., degree_bound
    /// 当 eval_point != 0 时，需要在 (lift_var - eval_point) 的幂次展开中工作
    /// 预计算 (lift_var - eval_point) 的幂次用于修正项构造
    std::vector<MultiPoly> shift_powers;
    shift_powers.reserve(degree_bound + 1);
    shift_powers.push_back(MultiPoly(Rational(1), vars)); // (lift_var - eval_point)^0 = 1
    if (lift_var_idx >= 0) {
        Monomial var_mono(vars.size(), 0);
        var_mono[lift_var_idx] = 1;
        std::vector<MultiPoly::Term> lin_terms;
        lin_terms.emplace_back(var_mono, Rational(1));
        if (!eval_point.is_zero()) {
            Monomial const_mono(vars.size(), 0);
            lin_terms.emplace_back(const_mono, -eval_point);
        }
        MultiPoly shift_lin(std::move(lin_terms), vars);
        shift_powers.push_back(shift_lin);
        for (int k = 2; k <= degree_bound; ++k) {
            shift_powers.push_back(shift_powers[k - 1] * shift_lin);
        }
    } else {
        for (int k = 1; k <= degree_bound; ++k) {
            shift_powers.push_back(MultiPoly(Rational(1), vars));
        }
    }

    for (int k = 1; k <= degree_bound; ++k) {
        /// 计算当前因子乘积，截断到 lift_var 次数 ≤ k
        MultiPoly product = factors[0];
        for (int i = 1; i < r; ++i) {
            product = product * factors[i];
            product = truncate_mod_var(product, lift_var, k + 1);
        }

        /// 计算误差 E = poly_trunc - product（仅比较 lift_var 次数 ≤ k 的部分）
        MultiPoly poly_trunc = truncate_mod_var(poly, lift_var, k + 1);
        MultiPoly error = poly_trunc - product;

        if (error.is_zero()) continue;

        /// 提取误差的 k 阶 Taylor 系数（关于 lift_var - eval_point）
        /// 方法：将 error 除以 (lift_var - eval_point)^k，然后在 lift_var = eval_point 处求值
        /// 对于 eval_point = 0，这等价于提取 lift_var^k 的系数
        MultiPoly error_k;
        if (eval_point.is_zero()) {
            /// 简化路径：直接提取 lift_var 指数恰为 k 的项
            std::vector<MultiPoly::Term> error_k_terms;
            for (const auto& term : error.terms()) {
                const Monomial& mono = term.first;
                int exp = (lift_var_idx >= 0 && static_cast<size_t>(lift_var_idx) < mono.size())
                          ? mono[lift_var_idx] : 0;
                if (exp == k) {
                    Monomial reduced = mono;
                    reduced[lift_var_idx] = 0;
                    error_k_terms.emplace_back(std::move(reduced), term.second);
                }
            }
            if (error_k_terms.empty()) continue;
            error_k = MultiPoly(std::move(error_k_terms), vars);
        } else {
            /// 非零求值点：计算 error / (lift_var - eval_point)^k 在 lift_var = eval_point 处的值
            /// 等价于对 error 关于 lift_var 求 k 次导数后除以 k! 再在 eval_point 处求值
            /// 实现：逐次除以 (lift_var - eval_point) 并求值
            MultiPoly remainder = error;
            for (int j = 0; j < k; ++j) {
                /// 除以 (lift_var - eval_point)：先在 eval_point 处求值确认余数为零
                /// 然后执行多项式除法
                MultiPoly eval_check = remainder.eval(lift_var, eval_point);
                if (!eval_check.is_zero()) {
                    /// 低阶误差未完全消除——不应发生，但防御性跳过
                    remainder = MultiPoly(Rational(0), vars);
                    break;
                }
                /// 执行除法：remainder / (lift_var - eval_point)
                try {
                    remainder = remainder.exact_div(shift_powers[1]);
                } catch (...) {
                    remainder = MultiPoly(Rational(0), vars);
                    break;
                }
            }
            /// 在 eval_point 处求值得到 error_k
            error_k = remainder.eval(lift_var, eval_point);
            if (error_k.is_zero()) continue;
        }

        /// 求解丢番图方程：找到修正项 δ₁, ..., δᵣ
        /// 使得 Σ δᵢ * cofactor_i ≡ error_k，其中 cofactor_i = ∏_{j≠i} fⱼ
        /// 因子在 lift_var = eval_point 处求值得到一元形式
        std::vector<MultiPoly> cofactors;
        cofactors.reserve(r);
        for (int i = 0; i < r; ++i) {
            MultiPoly cof(Rational(1), vars);
            for (int j = 0; j < r; ++j) {
                if (j == i) continue;
                MultiPoly fj_eval = factors[j].eval(lift_var, eval_point);
                cof = cof * fj_eval;
            }
            cofactors.push_back(cof);
        }

        std::vector<MultiPoly> corrections = multivariate_diophantine(
            cofactors, error_k, lift_var, eval_point, k + 1);

        if (corrections.size() != static_cast<size_t>(r)) continue;

        /// 更新因子：fᵢ = fᵢ + δᵢ * (lift_var - eval_point)^k
        for (int i = 0; i < r; ++i) {
            if (corrections[i].is_zero()) continue;

            /// 构造 δᵢ * (lift_var - eval_point)^k
            MultiPoly correction_poly = corrections[i] * shift_powers[k];
            factors[i] = factors[i] + correction_poly;
        }

        /// 验证：product of lifted factors ≡ poly mod (lift_var - eval_point)^(k+1)
        MultiPoly verify_product = factors[0];
        for (int i = 1; i < r; ++i) {
            verify_product = verify_product * factors[i];
            verify_product = truncate_mod_var(verify_product, lift_var, k + 1);
        }
        MultiPoly verify_error = poly_trunc - verify_product;
        if (!verify_error.is_zero()) {
            bool residual_ok = true;
            for (const auto& term : verify_error.terms()) {
                const Monomial& mono = term.first;
                int exp = (lift_var_idx >= 0 &&
                           static_cast<size_t>(lift_var_idx) < mono.size())
                          ? mono[lift_var_idx] : 0;
                if (exp <= k) {
                    residual_ok = false;
                    break;
                }
            }
            if (!residual_ok) {
                break;
            }
        }
    }

    return factors;
}
/**
 * @internal
 * @brief 一元多项式扩展 GCD（有理数域上）
 *
 * 给定 a, b ∈ Q[x]，计算 gcd(a, b) 及 Bezout 系数 s, t，
 * 使得 s*a + t*b = gcd(a, b)。
 *
 * @param[in]  a 第一个多项式
 * @param[in]  b 第二个多项式
 * @param[out] s Bezout 系数 s
 * @param[out] t Bezout 系数 t
 * @return gcd(a, b)
 */
static Polynomial<Rational> extended_gcd_poly(
    const Polynomial<Rational>& a,
    const Polynomial<Rational>& b,
    Polynomial<Rational>& s,
    Polynomial<Rational>& t)
{
    std::string var = a.variable_name;

    /// r0 = a, r1 = b
    /// s0*a + t0*b = r0
    /// s1*a + t1*b = r1
    Polynomial<Rational> r0 = a, r1 = b;
    Polynomial<Rational> s0({Rational(1)}, var), s1(var);  // s0=1, s1=0
    Polynomial<Rational> t0(var), t1({Rational(1)}, var);  // t0=0, t1=1

    while (!r1.is_zero()) {
        auto [q, r] = r0.div_mod(r1);

        Polynomial<Rational> r_new = r;
        Polynomial<Rational> s_new = s0 - q * s1;
        Polynomial<Rational> t_new = t0 - q * t1;

        r0 = r1; r1 = r_new;
        s0 = s1; s1 = s_new;
        t0 = t1; t1 = t_new;
    }

    /// 归一化使 gcd 为首一多项式
    if (!r0.is_zero()) {
        Rational lc = r0.lead_coeff();
        if (lc != Rational(1)) {
            Rational inv = Rational(1) / lc;
            for (auto& c : r0.coeffs) c = c * inv;
            for (auto& c : s0.coeffs) c = c * inv;
            for (auto& c : t0.coeffs) c = c * inv;
        }
    }

    s = s0;
    t = t0;
    return r0;
}

/**
 * @internal
 * @brief 将 MultiPoly 截断为关于指定变量次数 < degree_bound 的部分
 *
 * 保留所有项中指定变量指数严格小于 degree_bound 的项。
 *
 * @param[in] poly         输入多项式
 * @param[in] var          变量名
 * @param[in] degree_bound 次数上界（保留 < degree_bound 的项）
 * @return 截断后的多项式
 */
static MultiPoly truncate_mod_var(const MultiPoly& poly, const std::string& var,
                                  int degree_bound)
{
    if (poly.is_zero()) return poly;

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == var) { var_idx = static_cast<int>(i); break; }
    }
    /// 若变量不在列表中，多项式不含该变量，无需截断
    if (var_idx < 0) return poly;

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp < degree_bound) {
            result_terms.push_back(term);
        }
    }

    if (result_terms.empty()) return MultiPoly(Rational(0), vars);
    return MultiPoly(std::move(result_terms), vars);
}

/**
 * @brief 多元丢番图方程求解器
 *
 * 求解 s₁·f₁ + s₂·f₂ + ... + sᵣ·fᵣ ≡ target (mod (var - eval_point)^degree_bound)，
 * 其中各 fᵢ 两两互素。用于 Hensel 提升过程中计算修正项。
 *
 * 算法：
 * 1. 二因子情形（r=2）：使用扩展 GCD 求 Bezout 系数
 * 2. 一般情形（r>2）：递归归约——令 g = f₂·...·fᵣ，
 *    先解 s₁·f₁ + t·g ≡ target，再递归解 s₂·f₂ + ... + sᵣ·fᵣ ≡ t
 * 3. 对每个解截断为关于 var 的次数 < degree_bound
 *
 * @param[in] factors      互素因子列表 [f₁, ..., fᵣ]
 * @param[in] target       目标多项式 c
 * @param[in] var          变量名
 * @param[in] eval_point   求值点
 * @param[in] degree_bound 次数上界
 * @return 解多项式列表 [s₁, s₂, ..., sᵣ]
 *
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." §15.5.
 */
std::vector<MultiPoly> multivariate_diophantine(
    const std::vector<MultiPoly>& factors,
    const MultiPoly& target,
    const std::string& var,
    const Rational& eval_point,
    int degree_bound)
{
    int r = static_cast<int>(factors.size());
    if (r == 0) return {};
    if (r == 1) {
        /// 单因子情形：s₁ = target / f₁（精确除法后截断）
        try {
            MultiPoly s1 = target.exact_div(factors[0]);
            s1 = truncate_mod_var(s1, var, degree_bound);
            return {s1};
        } catch (...) {
            /// 若不整除，返回 target 本身（退化情形）
            return {truncate_mod_var(target, var, degree_bound)};
        }
    }

    const auto& vars = factors[0].variables();

    if (r == 2) {
        /// 二因子情形：使用扩展 GCD
        /// 将 MultiPoly 因子转换为一元 Polynomial<Rational> 进行 GCD 计算
        /// 因为在 Hensel 提升过程中，因子在求值点处本质为一元多项式

        /// 先对因子在 eval_point 处求值（去除 var 维度），得到一元多项式
        /// 但实际上 factors 可能已经是关于某个主变量的一元多项式
        /// 策略：尝试直接转换为一元多项式；若失败则在 var 处求值后转换

        Polynomial<Rational> f1_uni, f2_uni;
        bool converted = false;

        /// 尝试直接转换
        try {
            f1_uni = factors[0].to_univariate();
            f2_uni = factors[1].to_univariate();
            converted = true;
        } catch (...) {
            /// 若不是一元的，在 var 处求值后再转换
            try {
                MultiPoly f1_eval = factors[0].eval(var, eval_point);
                MultiPoly f2_eval = factors[1].eval(var, eval_point);
                f1_uni = f1_eval.to_univariate();
                f2_uni = f2_eval.to_univariate();
                converted = true;
            } catch (...) {
                converted = false;
            }
        }

        if (!converted) {
            /// 无法转换为一元多项式，返回零解作为退化处理
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        /// 计算扩展 GCD：s*f1 + t*f2 = gcd(f1, f2)
        Polynomial<Rational> s_coeff, t_coeff;
        Polynomial<Rational> g = extended_gcd_poly(f1_uni, f2_uni, s_coeff, t_coeff);

        /// 将 target 也转换为一元多项式
        Polynomial<Rational> target_uni;
        try {
            target_uni = target.to_univariate();
        } catch (...) {
            try {
                MultiPoly target_eval = target.eval(var, eval_point);
                target_uni = target_eval.to_univariate();
            } catch (...) {
                /// target 为常数
                if (target.is_constant()) {
                    Rational c = target.is_zero() ? Rational(0) : target.terms()[0].second;
                    target_uni = Polynomial<Rational>({c}, f1_uni.variable_name);
                } else {
                    return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
                }
            }
        }

        /// gcd 应为 1（因子互素），但处理一般情形
        /// s₁ = s * (target / gcd), s₂ = t * (target / gcd)
        Polynomial<Rational> scale;
        if (g.is_zero() || g.degree() < 0) {
            /// 退化情形
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        auto [quotient, remainder] = target_uni.div_mod(g);
        if (!remainder.is_zero()) {
            /// target 不能被 gcd 整除——不应发生于互素因子
            /// 退化处理：返回零
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        /// s₁ = s * quotient, s₂ = t * quotient
        Polynomial<Rational> s1_uni = s_coeff * quotient;
        Polynomial<Rational> s2_uni = t_coeff * quotient;

        /// 对 s₁ 取模 f₂，对 s₂ 取模 f₁，确保次数约束
        /// s₁*f₁ + s₂*f₂ = target，且 deg(s₁) < deg(f₂), deg(s₂) < deg(f₁)
        if (!f2_uni.is_zero() && s1_uni.degree() >= f2_uni.degree()) {
            auto [q1, r1] = s1_uni.div_mod(f2_uni);
            s1_uni = r1;
            /// 调整 s₂：s₂ = s₂ + q1 * f₁
            s2_uni = s2_uni + q1 * f1_uni;
        }

        /// 确定用于 from_univariate 的变量名
        std::string uni_var = f1_uni.variable_name;

        /// 转换回 MultiPoly
        MultiPoly s1_mp = MultiPoly::from_univariate(s1_uni, uni_var);
        MultiPoly s2_mp = MultiPoly::from_univariate(s2_uni, uni_var);

        /// 若变量集不匹配，嵌入到完整变量集
        if (s1_mp.variables() != vars && !vars.empty()) {
            /// 找到 uni_var 在 vars 中的位置
            int uni_idx = -1;
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i] == uni_var) { uni_idx = static_cast<int>(i); break; }
            }
            if (uni_idx >= 0) {
                std::vector<MultiPoly::Term> new_terms1, new_terms2;
                for (const auto& term : s1_mp.terms()) {
                    Monomial full_mono(vars.size(), 0);
                    if (!term.first.empty()) full_mono[uni_idx] = term.first[0];
                    new_terms1.emplace_back(std::move(full_mono), term.second);
                }
                for (const auto& term : s2_mp.terms()) {
                    Monomial full_mono(vars.size(), 0);
                    if (!term.first.empty()) full_mono[uni_idx] = term.first[0];
                    new_terms2.emplace_back(std::move(full_mono), term.second);
                }
                s1_mp = MultiPoly(std::move(new_terms1), vars);
                s2_mp = MultiPoly(std::move(new_terms2), vars);
            }
        }

        /// 截断为 degree < degree_bound
        s1_mp = truncate_mod_var(s1_mp, var, degree_bound);
        s2_mp = truncate_mod_var(s2_mp, var, degree_bound);

        return {s1_mp, s2_mp};
    }

    /// 一般情形（r > 2）：递归归约
    /// 令 g = f₂ * f₃ * ... * fᵣ
    MultiPoly g = factors[1];
    for (int i = 2; i < r; ++i) {
        g = g * factors[i];
    }

    /// 解二因子方程：s₁*f₁ + t*g ≡ target
    std::vector<MultiPoly> two_factors = {factors[0], g};
    std::vector<MultiPoly> two_solution = multivariate_diophantine(
        two_factors, target, var, eval_point, degree_bound);

    if (two_solution.size() != 2) {
        /// 退化情形
        std::vector<MultiPoly> result(r, MultiPoly(Rational(0), vars));
        return result;
    }

    /// s₁ 已确定
    MultiPoly s1 = two_solution[0];
    MultiPoly t_poly = two_solution[1];

    /// 递归解：s₂*f₂ + ... + sᵣ*fᵣ ≡ t_poly * g
    /// 因为原方程为 s₁*f₁ + t*g = target，
    /// 我们需要 s₂*f₂ + ... + sᵣ*fᵣ = t*g = target - s₁*f₁
    MultiPoly remaining_target = t_poly * g;
    remaining_target = truncate_mod_var(remaining_target, var, degree_bound);

    std::vector<MultiPoly> remaining_factors(factors.begin() + 1, factors.end());
    std::vector<MultiPoly> remaining_solution = multivariate_diophantine(
        remaining_factors, remaining_target, var, eval_point, degree_bound);

    /// 组装完整解
    std::vector<MultiPoly> result;
    result.push_back(std::move(s1));
    for (auto& sol : remaining_solution) {
        result.push_back(std::move(sol));
    }

    return result;
}
namespace detail {

/**
 * @brief 检测是否为线性多项式（主变量次数 ≤ 1）
 */
bool is_linear(const MultiPoly& poly, const std::string& main_var)
{
    return poly.degree(main_var) <= 1;
}

/**
 * @brief 检测差平方模式 a² - b²
 *
 * 判断多项式是否恰好有两项，一正一负，且两项的单项式各分量指数均为偶数。
 * 若匹配，返回 (a, b)，其中 a² 为正项，b² 为负项。
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (a, b)，否则返回 nullopt
 */
std::optional<std::pair<MultiPoly, MultiPoly>>
detect_difference_of_squares(const MultiPoly& poly)
{
    if (poly.num_terms() != 2) return std::nullopt;

    const auto& terms = poly.terms();
    const auto& vars = poly.variables();

    /// 确定哪项为正、哪项为负
    int pos_idx = -1, neg_idx = -1;
    if (terms[0].second > Rational(0) && terms[1].second < Rational(0)) {
        pos_idx = 0; neg_idx = 1;
    } else if (terms[1].second > Rational(0) && terms[0].second < Rational(0)) {
        pos_idx = 1; neg_idx = 0;
    } else {
        return std::nullopt;
    }

    /// 检查系数绝对值相等（都应为完全平方的系数）
    /// 对于差平方 a² - b²，正项系数和负项系数的绝对值不必相等
    /// 但两项的系数本身必须是完全平方数
    Rational pos_coeff = terms[pos_idx].second;
    Rational neg_coeff = -terms[neg_idx].second; // 取绝对值

    /// 检查系数是否为完全平方有理数
    /// 分子和分母都必须是完全平方数
    BigInt pos_num = pos_coeff.get_numerator();
    BigInt pos_den = pos_coeff.get_denominator();
    BigInt neg_num = neg_coeff.get_numerator();
    BigInt neg_den = neg_coeff.get_denominator();

    if (pos_num < BigInt(0) || neg_num < BigInt(0)) return std::nullopt;

    /// 检查是否为完全平方
    /// 注意：不能依赖 BigInt::sqrt()（该实现有 bug，对小整数返回错误值），
    /// 改用 double 估计取整后用精确乘法验证。
    auto is_perfect_square = [](const BigInt& n) -> std::pair<bool, BigInt> {
        if (n == BigInt(0)) return {true, BigInt(0)};
        if (n == BigInt(1)) return {true, BigInt(1)};
        if (n < BigInt(0)) return {false, BigInt(0)};
        double dn = n.to_double();
        long long est = (long long)std::llround(std::sqrt(dn));
        /// 在估计值附近搜索，抵消浮点误差
        for (long long cand = (est > 2 ? est - 2 : 0); cand <= est + 2; ++cand) {
            BigInt r(cand);
            if (r * r == n) return {true, r};
        }
        return {false, BigInt(0)};
    };

    auto [pos_num_sq, pos_num_root] = is_perfect_square(pos_num);
    if (!pos_num_sq) return std::nullopt;
    auto [pos_den_sq, pos_den_root] = is_perfect_square(pos_den);
    if (!pos_den_sq) return std::nullopt;
    auto [neg_num_sq, neg_num_root] = is_perfect_square(neg_num);
    if (!neg_num_sq) return std::nullopt;
    auto [neg_den_sq, neg_den_root] = is_perfect_square(neg_den);
    if (!neg_den_sq) return std::nullopt;

    /// 检查两项的单项式各分量指数均为偶数
    const Monomial& pos_mono = terms[pos_idx].first;
    const Monomial& neg_mono = terms[neg_idx].first;

    for (size_t i = 0; i < pos_mono.size(); ++i) {
        if (pos_mono[i] % 2 != 0) return std::nullopt;
    }
    for (size_t i = 0; i < neg_mono.size(); ++i) {
        if (neg_mono[i] % 2 != 0) return std::nullopt;
    }

    /// 构造 a：正项单项式指数减半，系数为 sqrt(pos_coeff)
    Monomial a_mono(vars.size(), 0);
    for (size_t i = 0; i < pos_mono.size(); ++i) {
        a_mono[i] = pos_mono[i] / 2;
    }
    Rational a_coeff(pos_num_root, pos_den_root);
    std::vector<MultiPoly::Term> a_terms = {{a_mono, a_coeff}};
    MultiPoly a(std::move(a_terms), vars);

    /// 构造 b：负项单项式指数减半，系数为 sqrt(neg_coeff)
    Monomial b_mono(vars.size(), 0);
    for (size_t i = 0; i < neg_mono.size(); ++i) {
        b_mono[i] = neg_mono[i] / 2;
    }
    Rational b_coeff(neg_num_root, neg_den_root);
    std::vector<MultiPoly::Term> b_terms = {{b_mono, b_coeff}};
    MultiPoly b(std::move(b_terms), vars);

    return std::make_pair(a, b);
}

/**
 * @brief 检测二项式幂模式 xⁿ - yⁿ
 *
 * 判断多项式是否为两个不同变量的纯幂之差，系数分别为 +1 和 -1，
 * 且幂次相同。若匹配，返回 (变量1, 变量2, 幂次)。
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (var1, var2, n)，否则返回 nullopt
 */
std::optional<std::tuple<std::string, std::string, int>>
detect_binomial_power(const MultiPoly& poly)
{
    if (poly.num_terms() != 2) return std::nullopt;

    const auto& terms = poly.terms();
    const auto& vars = poly.variables();

    /// 确定正项和负项
    int pos_idx = -1, neg_idx = -1;
    if (terms[0].second == Rational(1) && terms[1].second == Rational(-1)) {
        pos_idx = 0; neg_idx = 1;
    } else if (terms[1].second == Rational(1) && terms[0].second == Rational(-1)) {
        pos_idx = 1; neg_idx = 0;
    } else {
        return std::nullopt;
    }

    /// 检查每项的单项式是否为单变量的纯幂（恰好一个非零指数）
    auto get_single_var_power = [&](const Monomial& mono)
        -> std::optional<std::pair<int, int>> {
        int nonzero_count = 0;
        int var_index = -1;
        int exponent = 0;
        for (size_t i = 0; i < mono.size(); ++i) {
            if (mono[i] != 0) {
                nonzero_count++;
                var_index = static_cast<int>(i);
                exponent = mono[i];
            }
        }
        if (nonzero_count != 1 || exponent <= 1) return std::nullopt;
        return std::make_pair(var_index, exponent);
    };

    auto pos_info = get_single_var_power(terms[pos_idx].first);
    if (!pos_info) return std::nullopt;
    auto neg_info = get_single_var_power(terms[neg_idx].first);
    if (!neg_info) return std::nullopt;

    auto [pos_var_idx, pos_exp] = *pos_info;
    auto [neg_var_idx, neg_exp] = *neg_info;

    /// 幂次必须相同，变量必须不同
    if (pos_exp != neg_exp) return std::nullopt;
    if (pos_var_idx == neg_var_idx) return std::nullopt;

    return std::make_tuple(vars[pos_var_idx], vars[neg_var_idx], pos_exp);
}

/**
 * @brief 提取公因子单项式
 *
 * 计算所有项单项式的逐分量最小值（GCD 单项式），将其从每项中减去。
 *
 * @param[in] poly 输入多项式
 * @return pair(公因子单项式, 提取后的多项式)
 */
std::pair<Monomial, MultiPoly>
extract_common_monomial(const MultiPoly& poly)
{
    const auto& vars = poly.variables();
    Monomial identity(vars.size(), 0);

    if (poly.is_zero() || poly.num_terms() == 0) {
        return {identity, poly};
    }

    const auto& terms = poly.terms();

    /// 计算逐分量最小值
    Monomial gcd_mono = terms[0].first;
    for (size_t t = 1; t < terms.size(); ++t) {
        const Monomial& mono = terms[t].first;
        for (size_t i = 0; i < gcd_mono.size(); ++i) {
            int val = (i < mono.size()) ? mono[i] : 0;
            if (val < gcd_mono[i]) gcd_mono[i] = val;
        }
    }

    /// 检查是否为平凡（全零）
    bool is_trivial = true;
    for (size_t i = 0; i < gcd_mono.size(); ++i) {
        if (gcd_mono[i] != 0) { is_trivial = false; break; }
    }
    if (is_trivial) return {identity, poly};

    /// 从每项中减去 GCD 单项式
    std::vector<MultiPoly::Term> quotient_terms;
    quotient_terms.reserve(terms.size());
    for (const auto& term : terms) {
        Monomial new_mono(vars.size(), 0);
        for (size_t i = 0; i < vars.size(); ++i) {
            int orig = (i < term.first.size()) ? term.first[i] : 0;
            new_mono[i] = orig - gcd_mono[i];
        }
        quotient_terms.emplace_back(std::move(new_mono), term.second);
    }

    MultiPoly quotient(std::move(quotient_terms), vars);
    return {gcd_mono, quotient};
}

/**
 * @brief 齐次二元多项式分解
 *
 * 对齐次二元多项式 f(x, y)，令 y=1 得到一元多项式 f(x, 1)，
 * 对其进行一元因式分解，然后将每个因子重新齐次化。
 *
 * @param[in] poly 输入多项式（须为齐次二元）
 * @return 若适用则返回分解结果，否则返回 nullopt
 */
std::optional<MultiFactorResult>
factor_homogeneous_bivariate(const MultiPoly& poly)
{
    if (!poly.is_homogeneous() || poly.num_vars() != 2) return std::nullopt;

    const auto& vars = poly.variables();
    std::string var_x = vars[0];
    std::string var_y = vars[1];

    /// 去齐次化：令 y = 1，得到关于 x 的一元多项式
    MultiPoly dehomogenized = poly.eval(var_y, Rational(1));

    /// 转换为一元多项式
    Polynomial<Rational> uni_poly;
    try {
        uni_poly = dehomogenized.to_univariate();
    } catch (...) {
        return std::nullopt;
    }

    /// 一元因式分解
    std::vector<Polynomial<Rational>> uni_factors = factor_univariate_bridge(uni_poly);
    if (uni_factors.size() <= 1) return std::nullopt;

    /// 重新齐次化每个因子
    MultiFactorResult result;
    result.constant = poly.numeric_content();

    for (const auto& uf : uni_factors) {
        int d = uf.degree();
        /// 对一元因子的每项 c*x^k，齐次化为 c*x^k*y^(d-k)
        std::vector<MultiPoly::Term> homo_terms;
        for (int k = 0; k <= d; ++k) {
            Rational coeff = (k < static_cast<int>(uf.coeffs.size())) ? uf.coeffs[k] : Rational(0);
            if (coeff.is_zero()) continue;
            Monomial mono(2, 0);
            mono[0] = k;       // x 的指数
            mono[1] = d - k;   // y 的指数
            homo_terms.emplace_back(std::move(mono), coeff);
        }
        MultiPoly homo_factor(std::move(homo_terms), vars);
        homo_factor = homo_factor.make_primitive();
        result.factors.push_back(std::move(homo_factor));
        result.multiplicities.push_back(1);
    }

    /// 调整常数因子使乘积等于原多项式
    /// 计算所有因子的乘积
    MultiPoly product(Rational(1), vars);
    for (const auto& f : result.factors) {
        product = product * f;
    }
    if (!product.is_zero() && !poly.is_zero()) {
        /// constant = poly / product 的比例系数
        try {
            MultiPoly ratio = poly.exact_div(product);
            if (ratio.is_constant()) {
                result.constant = ratio.numeric_content();
                if (!ratio.terms().empty() && ratio.terms()[0].second < Rational(0)) {
                    result.constant = -result.constant;
                }
            }
        } catch (...) {
            /// 若除法失败，使用首项系数比
            if (!product.terms().empty() && !poly.terms().empty()) {
                result.constant = poly.terms()[0].second / product.terms()[0].second;
            }
        }
    }

    return result;
}

} // namespace detail
} // namespace lamina
