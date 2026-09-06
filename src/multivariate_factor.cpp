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
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <new>
#include "internal/multivariate_factor_support.hpp"
namespace LMCAS {
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
    catch (const std::runtime_error&) { return false; }
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

MultiPoly multivariate_content(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return MultiPoly();
    /// 若 main_var 不在变量列表中,或多项式在 main_var 上次数为 0,容度就是多项式本身
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

    /// content 的变量集是去掉 main_var 后的辅助变量集,
    /// 需要将其嵌入到 poly 的完整变量集中才能正确执行 exact_div.
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
    /// 基本情形:两者均为一元多项式
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
    /// 递归情形:求值-插值法
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
 * 对每个含 main_var 的项,新系数 = 原系数 x 指数,新指数 = 原指数 - 1.
 * 指数为 0 的项被丢弃.
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 求导变量名
 * @return partialpoly/partialmain_var
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
        /// main_var 不在变量列表中,导数为零
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
 * @brief 多元无平方因子分解(Yun 算法的多元推广)
 *
 * 通过计算 gcd(f, partialf/partialx_main) 检测重因子,将多项式分解为
 * f = f_1 * f_2^2 * f₃^3 * ... 的形式,其中每个 fᵢ 无平方.
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 无平方因子分解结果,components[i] 的重数为 i+1
 */
SquareFreeDecomp square_free_decompose(const MultiPoly& poly, const std::string& main_var)
{
    /// 零多项式或常数多项式:直接返回
    if (poly.is_zero() || poly.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 若 main_var 不在变量列表中或次数为 0,多项式无重因子
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0 || poly.degree(main_var) <= 0) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 1:计算 f' = partialf/partialmain_var
    MultiPoly f_prime = formal_derivative(poly, main_var);

    /// 步骤 2:计算 g = gcd(f, f');特征 0 下 f' 为零时,f 直接作为单个无平方因子.
    if (f_prime.is_zero()) {
        return SquareFreeDecomp{{poly}};
    }

    MultiPoly g = multivariate_gcd(poly, f_prime);

    /// 步骤 3:若 g 为常数,f 已无平方
    if (g.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 4:Yun 算法迭代提取无平方分量
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

    /// 过滤掉常数为 1 的平凡分量(保留结构但不影响语义)
    /// 注意:components[i] 的重数为 i+1
    return SquareFreeDecomp{std::move(components)};
}
/**
 * @internal
 * @brief 首项系数预计算(Wang's leading coefficient trick)
 *
 * 在 Hensel 提升前,将原多项式关于主变量的首项系数 lc(f, x_main) 分解为
 * 辅助变量的因子,并按求值后的首项系数值分配给各一元因子,
 * 使提升过程保持正确的首项系数与本原内容.
 *
 * 算法步骤:
 * 1. 计算 lc = lc(f, x_main),为辅助变量的多项式
 * 2. 若 lc 为常数,无需预计算
 * 3. 将 lc 在求值点处求值得到 lc_eval(有理数)
 * 4. 对每个一元因子 fᵢ,其首项系数 lc_i 应整除 lc_eval
 * 5. 按 lc_i 在 lc_eval 中的贡献比例,将 lc 的因子分配给各 fᵢ
 * 6. 将分配的 lc 部分乘入各 fᵢ,使提升后首项系数正确
 *
 * @param[in]     poly              原始多元多项式
 * @param[in]     main_var          主变量名
 * @param[in]     eval_points       辅助变量到求值点的映射
 * @param[in,out] univariate_factors 一元因子列表,函数会修改各因子的首项系数
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 * @see D. Y. Y. Yun, "On Square-Free Decomposition Algorithms,"
 *      Proceedings of SYMSAC 1976.
 */

} // namespace LMCAS
