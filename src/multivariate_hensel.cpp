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
namespace LMCAS {
static MultiPoly truncate_mod_var(const MultiPoly& poly, const std::string& var,
                                  int degree_bound);
static MultiPoly embed_hensel_terms(const MultiPoly& polynomial,
    const std::vector<std::string>& variables) {
    if (polynomial.variables() == variables) { return polynomial; }
    if (polynomial.is_constant()) {
        return MultiPoly(
            polynomial.is_zero() ? Rational(0) : polynomial.terms().front().second,
            variables);
    }
    std::vector<size_t> coordinates;
    coordinates.reserve(polynomial.variables().size());
    for (const auto& variable : polynomial.variables()) {
        auto position = std::find(variables.begin(), variables.end(), variable);
        if (position == variables.end()) {
            throw std::invalid_argument("Hensel embedding: variable is absent from target ring");
        }
        coordinates.push_back(static_cast<size_t>(position - variables.begin()));
    }
    std::vector<MultiPoly::Term> terms;
    terms.reserve(polynomial.terms().size());
    for (const auto& term : polynomial.terms()) {
        Monomial monomial(variables.size(), 0);
        for (size_t i = 0; i < term.first.size(); ++i) {
            monomial[coordinates[i]] = term.first[i];
        }
        terms.emplace_back(std::move(monomial), term.second);
    }
    return MultiPoly(std::move(terms), variables);
}


static std::vector<MultiPoly> embed_without_lifting(
    const std::vector<Polynomial<Rational>>& univariate_factors,
    const std::vector<std::string>& vars) {
    std::vector<MultiPoly> result;
    result.reserve(univariate_factors.size());
    for (const auto& uf : univariate_factors) {
        result.push_back(embed_hensel_terms(
            MultiPoly::from_univariate(uf, uf.variable_name), vars));
    }
    return result;
}


static std::vector<MultiPoly> lifting_shift_powers(
    const std::vector<std::string>& vars, int lift_var_idx,
    const Rational& eval_point, int degree_bound) {
    std::vector<MultiPoly> shift_powers;
    shift_powers.reserve(degree_bound + 1);
    shift_powers.push_back(MultiPoly(Rational(1), vars)); /**< 零次幂 (lift_var - eval_point)^0 = 1。 */
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
    return shift_powers;
}

static MultiPoly lifting_error_coefficient(const MultiPoly& error,
    const std::string& lift_var, int lift_var_idx, const Rational& eval_point,
    int k, const std::vector<MultiPoly>& shift_powers,
    const std::vector<std::string>& vars) {
    MultiPoly error_k;
    if (eval_point.is_zero()) {
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
        if (error_k_terms.empty()) { return MultiPoly(Rational(0), vars); }
        error_k = MultiPoly(std::move(error_k_terms), vars);
    } else {
        /**
         * @brief 逐次除以 (lift_var - eval_point)，再在 eval_point 处求值。
         * 得到 error / (lift_var - eval_point)^k 的值，
         * 等价于 error 的 k 阶导数在该点的值除以 k!。
         */
        MultiPoly remainder = error;
        for (int j = 0; j < k; ++j) {
            MultiPoly eval_check = remainder.eval(lift_var, eval_point);
            if (!eval_check.is_zero()) {
                remainder = MultiPoly(Rational(0), vars);
                break;
            }
            try {
                remainder = remainder.exact_div(shift_powers[1]);
            } catch (const std::runtime_error&) {
                remainder = MultiPoly(Rational(0), vars);
                break;
            }
        }
        error_k = embed_hensel_terms(remainder.eval(lift_var, eval_point), vars);
        if (error_k.is_zero()) { return error_k; }
    }
    return error_k;
}

static std::vector<MultiPoly> evaluated_lifting_cofactors(
    const std::vector<MultiPoly>& factors, const std::string& lift_var,
    const Rational& eval_point, const std::vector<std::string>& vars) {
    const int r = static_cast<int>(factors.size());
    std::vector<MultiPoly> cofactors;
    cofactors.reserve(r);
    for (int i = 0; i < r; ++i) {
        MultiPoly cof(Rational(1), vars);
        for (int j = 0; j < r; ++j) {
            if (j == i) { continue; }
            MultiPoly fj_eval = embed_hensel_terms(
                factors[j].eval(lift_var, eval_point), vars);
            cof = cof * fj_eval;
        }
        cofactors.push_back(cof);
    }
    return cofactors;
}

static bool lifting_residual_is_zero(const MultiPoly& verify_error,
                                     int lift_var_idx, int k) {
    if (verify_error.is_zero()) { return true; }
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
    return residual_ok;
}

static MultiPoly truncated_factor_product(const std::vector<MultiPoly>& factors,
    const std::string& variable, int degree_bound) {
    MultiPoly product = factors[0];
    for (size_t i = 1; i < factors.size(); ++i) {
        product = product * factors[i];
        product = truncate_mod_var(product, variable, degree_bound);
    }
    return product;
}

static bool apply_hensel_lift_steps(const MultiPoly& poly,
    std::vector<MultiPoly>& factors, const std::string& lift_var, int lift_var_idx,
    const Rational& eval_point, int degree_bound, int factor_count) {
    const auto& vars = poly.variables();
    auto shift_powers = lifting_shift_powers(vars, lift_var_idx,
                                             eval_point, degree_bound);
    for (int k = 1; k <= degree_bound; ++k) {
        MultiPoly product = truncated_factor_product(factors, lift_var, k + 1);
        MultiPoly poly_trunc = truncate_mod_var(poly, lift_var, k + 1);
        MultiPoly error = poly_trunc - product;
        if (error.is_zero()) { continue; }
        auto error_k = lifting_error_coefficient(error, lift_var, lift_var_idx,
                                                  eval_point, k, shift_powers, vars);
        if (error_k.is_zero()) { continue; }
        auto cofactors = evaluated_lifting_cofactors(factors, lift_var,
                                                       eval_point, vars);
        std::vector<MultiPoly> corrections = multivariate_diophantine(
            cofactors, error_k, lift_var, eval_point, k + 1);
        if (corrections.size() != static_cast<size_t>(factor_count)) { return false; }
        for (int i = 0; i < factor_count; ++i) {
            if (corrections[i].is_zero()) { continue; }
            MultiPoly correction_poly = corrections[i] * shift_powers[k];
            factors[i] = factors[i] + correction_poly;
        }
        MultiPoly verify_product = truncated_factor_product(factors, lift_var, k + 1);
        MultiPoly verify_error = poly_trunc - verify_product;
        if (!lifting_residual_is_zero(verify_error, lift_var_idx, k)) { break; }
    }
    return true;
}

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
    if (r == 0) { return {}; }
    if (r == 1) {
        return {poly};
    }
    if (degree_bound <= 0) {
        return embed_without_lifting(univariate_factors, poly.variables());
    }

    const auto& vars = poly.variables();

    /// 确定主变量(一元因子的变量)和提升变量在完整变量集中的位置
    std::string main_var = univariate_factors[0].variable_name;
    int lift_var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == lift_var) { lift_var_idx = static_cast<int>(i); }
    }

    for (const auto& factor : univariate_factors) {
        if (factor.variable_name != main_var) {
            throw std::invalid_argument("Hensel lifting: univariate factor ring mismatch");
        }
    }
    auto factors = embed_without_lifting(univariate_factors, vars);

    /// lift_var 位于变量域之外时,返回嵌入后的原因子.
    if (lift_var_idx < 0) {
        return factors;
    }

    MultiPoly base_product(Rational(1), vars);
    for (const auto& factor : factors) { base_product = base_product * factor; }
    if (base_product.eval(lift_var, eval_point) != poly.eval(lift_var, eval_point)) {
        return {};
    }

    if (!apply_hensel_lift_steps(
            poly, factors, lift_var, lift_var_idx, eval_point, degree_bound, r)) {
        return {};
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
    if (poly.is_zero()) { return poly; }

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == var) { var_idx = static_cast<int>(i); break; }
    }
    /// 若变量不在列表中,多项式不含该变量,无需截断
    if (var_idx < 0) { return poly; }

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp < degree_bound) {
            result_terms.push_back(term);
        }
    }

    if (result_terms.empty()) { return MultiPoly(Rational(0), vars); }
    return MultiPoly(std::move(result_terms), vars);
}

static bool diophantine_univariate_factors(
    const std::vector<MultiPoly>& factors, const std::string& var,
    const Rational& eval_point, Polynomial<Rational>& f1_uni,
    Polynomial<Rational>& f2_uni) {
    bool converted = false;

    try {
        f1_uni = factors[0].to_univariate();
        f2_uni = factors[1].to_univariate();
        converted = true;
    } catch (const std::invalid_argument&) {
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
    if (converted && f1_uni.variable_name != f2_uni.variable_name) {
        /**
         * @brief 常数转换的默认变量名不代表其系数环。
         * 仅常数可嵌入另一变量环；不同活跃变量须保持区分。
         */
        if (f1_uni.degree() <= 0) {
            f1_uni = Polynomial<Rational>(f1_uni.coeffs, f2_uni.variable_name);
        } else if (f2_uni.degree() <= 0) {
            f2_uni = Polynomial<Rational>(f2_uni.coeffs, f1_uni.variable_name);
        } else {
            return false;
        }
    }
    return converted;
}

static bool diophantine_target(const MultiPoly& target,
    const std::string& var, const Rational& eval_point,
    const std::string& main_var, Polynomial<Rational>& target_uni) {
    try {
        target_uni = target.to_univariate();
        /**
         * @brief 同一一元变量域内可直接转换。
         * 仅含提升变量的目标须先在提升点求值。
         */
        if (target_uni.degree() > 0 &&
            target_uni.variable_name != main_var) {
            throw std::logic_error("target uses the lift variable");
        }
    } catch (const std::logic_error&) {
        try {
            MultiPoly target_eval = target.eval(var, eval_point);
            target_uni = target_eval.to_univariate();
        } catch (const std::invalid_argument&) {
            if (target.is_constant()) {
                Rational c = target.is_zero() ? Rational(0) : target.terms()[0].second;
                target_uni = Polynomial<Rational>({c}, main_var);
            } else {
                return false;
            }
        }
    }
    if (target_uni.degree() <= 0) {
        target_uni.variable_name = main_var;
    } else if (target_uni.variable_name != main_var) {
        return false;
    }
    return true;
}

static std::vector<MultiPoly> solve_diophantine_pair(
    const std::vector<MultiPoly>& factors, const MultiPoly& target,
    const std::string& var, const Rational& eval_point, int degree_bound) {
    const auto& vars = factors[0].variables();
    Polynomial<Rational> f1_uni, f2_uni;
    if (!diophantine_univariate_factors(factors, var, eval_point,
                                       f1_uni, f2_uni)) {
        return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
    }
    Polynomial<Rational> s_coeff, t_coeff;
    Polynomial<Rational> g = extended_gcd_poly(f1_uni, f2_uni, s_coeff, t_coeff);

    Polynomial<Rational> target_uni;
    if (!diophantine_target(target, var, eval_point,
                            f1_uni.variable_name, target_uni)) {
        return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
    }


    /**
     * @brief 用 target/gcd 缩放 Bézout 系数，兼容 gcd 非 1 的情形。
     * 互素因子的 gcd 为 1；s_1 = s*(target/gcd)，s_2 = t*(target/gcd)。
     */
    Polynomial<Rational> scale;
    if (g.is_zero() || g.degree() < 0) {
        return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
    }

    auto [quotient, remainder] = target_uni.div_mod(g);
    if (!remainder.is_zero()) {
        return {MultiPoly(Rational(0), vars), MultiPoly(Rational(0), vars)};
    }

    Polynomial<Rational> s1_uni = s_coeff * quotient;
    Polynomial<Rational> s2_uni = t_coeff * quotient;

    /**
     * @brief 按 f_2、f_1 约化系数并保持 s_1*f_1 + s_2*f_2 = target。
     * 次数约束为 deg(s_1) < deg(f_2)、deg(s_2) < deg(f_1)。
     */
    if (!f2_uni.is_zero() && s1_uni.degree() >= f2_uni.degree()) {
        auto [q1, r1] = s1_uni.div_mod(f2_uni);
        s1_uni = r1;
        s2_uni = s2_uni + q1 * f1_uni;
    }
    if (!(s1_uni * f1_uni + s2_uni * f2_uni == target_uni)) { return {}; }

    std::string uni_var = f1_uni.variable_name;
    MultiPoly s1_mp = embed_hensel_terms(
        MultiPoly::from_univariate(s1_uni, uni_var), vars);
    MultiPoly s2_mp = embed_hensel_terms(
        MultiPoly::from_univariate(s2_uni, uni_var), vars);
    s1_mp = truncate_mod_var(s1_mp, var, degree_bound);
    s2_mp = truncate_mod_var(s2_mp, var, degree_bound);

    return {s1_mp, s2_mp};
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
    if (r == 0) { return {}; }
    if (r == 1) {
        /// 单因子情形:s_1 = target / f_1(精确除法后截断)
        try {
            MultiPoly s1 = target.exact_div(factors[0]);
            s1 = truncate_mod_var(s1, var, degree_bound);
            return {s1};
        } catch (const std::runtime_error&) {
            return {};
        }
    }

    const auto& vars = factors[0].variables();

    if (r == 2) {
        return solve_diophantine_pair(factors, target, var, eval_point,
                                      degree_bound);
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

} // namespace LMCAS
