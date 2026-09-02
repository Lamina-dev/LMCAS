/**
 * @file multivariate_factor.cpp
 * @brief 多元多项式因式分解器实现.
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
static MultiPoly truncate_mod_var(const MultiPoly& poly, const std::string& var,
                                  int degree_bound);
/**
 * @brief 多元 Hensel 提升(单变量步)
 *
 * 将一元因子逐变量提升,从 f(x_1, a_2, ..., a_n) 的因子恢复到包含 lift_var 的因子.
 * 每次引入一个辅助变量,计算修正项直到达到次数上界.
 *
 * 算法(参考 Geddes, Czapor, Labahn Section15.6):
 * 1. 将一元因子嵌入到完整变量集中
 * 2. 对 k = 1, ..., degree_bound:
 *    a. 计算当前因子乘积(截断到 lift_var 次数 <= k)
 *    b. 计算误差 E = poly - product(提取 lift_var 次数恰为 k 的部分)
 *    c. 若误差为零则继续
 *    d. 用丢番图求解器分配误差到各因子
 *    e. 将修正项加到各因子上
 *    f. 验证:product of lifted factors == poly mod (lift_var - eval_point)^(k+1)
 * 3. 返回提升后的因子
 *
 * @param[in] poly              原始多元多项式
 * @param[in] univariate_factors 一元分解得到的因子列表
 * @param[in] lift_var          当前提升的变量名
 * @param[in] eval_point        该变量的求值点
 * @param[in] degree_bound      提升次数上界
 * @return 提升后的多元因子列表
 *
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." Section15.6.
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
        /// 无需提升:将一元因子嵌入完整变量集返回
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

    /// 确定主变量(一元因子的变量)和提升变量在完整变量集中的位置
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

    /// lift_var 位于变量域之外时,返回嵌入后的原因子.
    if (lift_var_idx < 0) {
        return factors;
    }

    /// 迭代提升:对 k = 1, ..., degree_bound
    /// 当 eval_point != 0 时,需要在 (lift_var - eval_point) 的幂次展开中工作
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
        /// 计算当前因子乘积,截断到 lift_var 次数 <= k
        MultiPoly product = factors[0];
        for (int i = 1; i < r; ++i) {
            product = product * factors[i];
            product = truncate_mod_var(product, lift_var, k + 1);
        }

        /// 计算误差 E = poly_trunc - product(仅比较 lift_var 次数 <= k 的部分)
        MultiPoly poly_trunc = truncate_mod_var(poly, lift_var, k + 1);
        MultiPoly error = poly_trunc - product;

        if (error.is_zero()) continue;

        /// 提取误差的 k 阶 Taylor 系数(关于 lift_var - eval_point)
        /// 方法:将 error 除以 (lift_var - eval_point)^k,然后在 lift_var = eval_point 处求值
        /// 对于 eval_point = 0,这等价于提取 lift_var^k 的系数
        MultiPoly error_k;
        if (eval_point.is_zero()) {
            /// 简化路径:直接提取 lift_var 指数恰为 k 的项
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
            /// 非零求值点:计算 error / (lift_var - eval_point)^k 在 lift_var = eval_point 处的值
            /// 等价于对 error 关于 lift_var 求 k 次导数后除以 k! 再在 eval_point 处求值
            /// 实现:逐次除以 (lift_var - eval_point) 并求值
            MultiPoly remainder = error;
            for (int j = 0; j < k; ++j) {
                /// 除以 (lift_var - eval_point):先在 eval_point 处求值确认余数为零
                /// 然后执行多项式除法
                MultiPoly eval_check = remainder.eval(lift_var, eval_point);
                if (!eval_check.is_zero()) {
                    /// 低阶误差残留时终止本轮提升,并将余式重置为零.
                    remainder = MultiPoly(Rational(0), vars);
                    break;
                }
                /// 执行除法:remainder / (lift_var - eval_point)
                try {
                    remainder = remainder.exact_div(shift_powers[1]);
                } catch (const std::runtime_error&) {
                    remainder = MultiPoly(Rational(0), vars);
                    break;
                }
            }
            /// 在 eval_point 处求值得到 error_k
            error_k = remainder.eval(lift_var, eval_point);
            if (error_k.is_zero()) continue;
        }

        /// 求解丢番图方程:找到修正项 delta_1, ..., deltaᵣ
        /// 使得 Σ deltaᵢ * cofactor_i == error_k,其中 cofactor_i = prod_{j!=i} fⱼ
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

        /// 更新因子:fᵢ = fᵢ + deltaᵢ * (lift_var - eval_point)^k
        for (int i = 0; i < r; ++i) {
            if (corrections[i].is_zero()) continue;

            /// 构造 deltaᵢ * (lift_var - eval_point)^k
            MultiPoly correction_poly = corrections[i] * shift_powers[k];
            factors[i] = factors[i] + correction_poly;
        }

        /// 验证:product of lifted factors == poly mod (lift_var - eval_point)^(k+1)
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
 * @brief 一元多项式扩展 GCD(有理数域上)
 *
 * 给定 a, b  in  Q[x],计算 gcd(a, b) 及 Bezout 系数 s, t,
 * 使得 s*a + t*b = gcd(a, b).
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
 * 保留所有项中指定变量指数严格小于 degree_bound 的项.
 *
 * @param[in] poly         输入多项式
 * @param[in] var          变量名
 * @param[in] degree_bound 次数上界(保留 < degree_bound 的项)
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
    /// 若变量不在列表中,多项式不含该变量,无需截断
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
 * 求解 s_1*f_1 + s_2*f_2 + ... + sᵣ*fᵣ == target (mod (var - eval_point)^degree_bound),
 * 其中各 fᵢ 两两互素.用于 Hensel 提升过程中计算修正项.
 *
 * 算法:
 * 1. 二因子情形(r=2):使用扩展 GCD 求 Bezout 系数
 * 2. 一般情形(r>2):递归归约--令 g = f_2*...*fᵣ,
 *    先解 s_1*f_1 + t*g == target,再递归解 s_2*f_2 + ... + sᵣ*fᵣ == t
 * 3. 对每个解截断为关于 var 的次数 < degree_bound
 *
 * @param[in] factors      互素因子列表 [f_1, ..., fᵣ]
 * @param[in] target       目标多项式 c
 * @param[in] var          变量名
 * @param[in] eval_point   求值点
 * @param[in] degree_bound 次数上界
 * @return 解多项式列表 [s_1, s_2, ..., sᵣ]
 *
 * @see Geddes, Czapor, Labahn. "Algorithms for Computer Algebra." Section15.5.
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
        /// 单因子情形:s_1 = target / f_1(精确除法后截断)
        try {
            MultiPoly s1 = target.exact_div(factors[0]);
            s1 = truncate_mod_var(s1, var, degree_bound);
            return {s1};
        } catch (const std::runtime_error&) {
            /// 若不整除,返回 target 本身(退化情形)
            return {truncate_mod_var(target, var, degree_bound)};
        }
    }

    const auto& vars = factors[0].variables();

    if (r == 2) {
        /// 二因子情形:使用扩展 GCD
        /// 将 MultiPoly 因子转换为一元 Polynomial<Rational> 进行 GCD 计算
        /// 因为在 Hensel 提升过程中,因子在求值点处本质为一元多项式

        /// 先对因子在 eval_point 处求值(去除 var 维度),得到一元多项式
        /// 但实际上 factors 可能已经是关于某个主变量的一元多项式
        /// 策略:尝试直接转换为一元多项式;若失败则在 var 处求值后转换

        Polynomial<Rational> f1_uni, f2_uni;
        bool converted = false;

        /// 尝试直接转换
        try {
            f1_uni = factors[0].to_univariate();
            f2_uni = factors[1].to_univariate();
            converted = true;
        } catch (const std::invalid_argument&) {
            /// 多元输入先在 var 处求值,再转换为一元多项式.
            try {
                MultiPoly f1_eval = factors[0].eval(var, eval_point);
                MultiPoly f2_eval = factors[1].eval(var, eval_point);
                f1_uni = f1_eval.to_univariate();
                f2_uni = f2_eval.to_univariate();
                converted = true;
            } catch (const std::invalid_argument&) {
                converted = false;
            }
        }

        if (!converted) {
            /// 一元转换未决时返回零解,表示退化提升结果.
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        /// 计算扩展 GCD:s*f1 + t*f2 = gcd(f1, f2)
        Polynomial<Rational> s_coeff, t_coeff;
        Polynomial<Rational> g = extended_gcd_poly(f1_uni, f2_uni, s_coeff, t_coeff);

        /// 将 target 也转换为一元多项式
        Polynomial<Rational> target_uni;
        try {
            target_uni = target.to_univariate();
            /// 直接转换只有在目标与因子位于同一一元变量域时有效。
            /// 若目标仅含提升变量，必须先在提升点求值。
            if (target_uni.degree() > 0 &&
                target_uni.variable_name != f1_uni.variable_name) {
                throw std::logic_error("target uses the lift variable");
            }
        } catch (const std::logic_error&) {
            try {
                MultiPoly target_eval = target.eval(var, eval_point);
                target_uni = target_eval.to_univariate();
            } catch (const std::invalid_argument&) {
                /// target 为常数
                if (target.is_constant()) {
                    Rational c = target.is_zero() ? Rational(0) : target.terms()[0].second;
                    target_uni = Polynomial<Rational>({c}, f1_uni.variable_name);
                } else {
                    return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
                }
            }
        }
        /// 常数多项式没有固有变量域；to_univariate() 会为其选择默认
        /// 变量名。后续一元除法必须使用因子所在的变量域。
        if (target_uni.degree() <= 0) {
            target_uni.variable_name = f1_uni.variable_name;
        } else if (target_uni.variable_name != f1_uni.variable_name) {
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }


        /// gcd 应为 1(因子互素),但处理一般情形
        /// s_1 = s * (target / gcd), s_2 = t * (target / gcd)
        Polynomial<Rational> scale;
        if (g.is_zero() || g.degree() < 0) {
            /// 退化情形
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        auto [quotient, remainder] = target_uni.div_mod(g);
        if (!remainder.is_zero()) {
            /// 互素前提失效时返回零解,表示 Bézout 提升退化.
            return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
        }

        /// s_1 = s * quotient, s_2 = t * quotient
        Polynomial<Rational> s1_uni = s_coeff * quotient;
        Polynomial<Rational> s2_uni = t_coeff * quotient;

        /// 对 s_1 取模 f_2,对 s_2 取模 f_1,确保次数约束
        /// s_1*f_1 + s_2*f_2 = target,且 deg(s_1) < deg(f_2), deg(s_2) < deg(f_1)
        if (!f2_uni.is_zero() && s1_uni.degree() >= f2_uni.degree()) {
            auto [q1, r1] = s1_uni.div_mod(f2_uni);
            s1_uni = r1;
            /// 调整 s_2:s_2 = s_2 + q1 * f_1
            s2_uni = s2_uni + q1 * f1_uni;
        }

        /// 确定用于 from_univariate 的变量名
        std::string uni_var = f1_uni.variable_name;

        /// 转换回 MultiPoly
        MultiPoly s1_mp = MultiPoly::from_univariate(s1_uni, uni_var);
        MultiPoly s2_mp = MultiPoly::from_univariate(s2_uni, uni_var);

        /// 若变量集不匹配,嵌入到完整变量集
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

    /// 一般情形(r > 2):递归归约
    /// 令 g = f_2 * f₃ * ... * fᵣ
    MultiPoly g = factors[1];
    for (int i = 2; i < r; ++i) {
        g = g * factors[i];
    }

    /// 解二因子方程:s_1*f_1 + t*g == target
    std::vector<MultiPoly> two_factors = {factors[0], g};
    std::vector<MultiPoly> two_solution = multivariate_diophantine(
        two_factors, target, var, eval_point, degree_bound);

    if (two_solution.size() != 2) {
        /// 退化情形
        std::vector<MultiPoly> result(r, MultiPoly(Rational(0), vars));
        return result;
    }

    /// s_1 已确定
    MultiPoly s1 = two_solution[0];
    MultiPoly t_poly = two_solution[1];

    /// 递归解:s_2*f_2 + ... + sᵣ*fᵣ == t_poly * g
    /// 因为原方程为 s_1*f_1 + t*g = target,
    /// 我们需要 s_2*f_2 + ... + sᵣ*fᵣ = t*g = target - s_1*f_1
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

} // namespace lamina
