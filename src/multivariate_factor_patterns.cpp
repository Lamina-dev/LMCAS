/**
 * @file multivariate_factor.cpp
 * @brief 多元多项式因式分解器实现.
 */
#include "multivariate_factor.hpp"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include "internal/multivariate_factor_support.hpp"
namespace LMCAS {
namespace detail {
static std::optional<Rational> rational_square_root(const Rational& coefficient) {
    const BigInt& numerator = coefficient.get_numerator();
    const BigInt& denominator = coefficient.get_denominator();
    if (numerator.is_negative()) { return std::nullopt; }
    BigInt numerator_root = numerator.sqrt();
    if (numerator_root * numerator_root != numerator) { return std::nullopt; }
    BigInt denominator_root = denominator.sqrt();
    if (denominator_root * denominator_root != denominator) { return std::nullopt; }
    return Rational(numerator_root, denominator_root);
}

static MultiPoly homogenize_factor(const Polynomial<Rational>& uf,
                                   const std::vector<std::string>& vars) {
    int d = uf.degree();
    // c*x^k becomes c*x^k*y^(d-k).
    std::vector<MultiPoly::Term> homo_terms;
    for (int k = 0; k <= d; ++k) {
        Rational coeff = (k < static_cast<int>(uf.coeffs.size())) ? uf.coeffs[k] : Rational(0);
        if (coeff.is_zero()) { continue; }
        Monomial mono(2, 0);
        mono[0] = k;
        mono[1] = d - k;
        homo_terms.emplace_back(std::move(mono), coeff);
    }
    MultiPoly homo_factor(std::move(homo_terms), vars);
    homo_factor = homo_factor.make_primitive();
    return homo_factor;
}

static void adjust_homogeneous_constant(const MultiPoly& poly,
                                        MultiFactorResult& result) {
    const auto& vars = poly.variables();
    MultiPoly product(Rational(1), vars);
    for (const auto& f : result.factors) {
        product = product * f;
    }
    if (!product.is_zero() && !poly.is_zero()) {
        try {
            MultiPoly ratio = poly.exact_div(product);
            if (ratio.is_constant()) {
                result.constant = ratio.numeric_content();
                if (!ratio.terms().empty() && ratio.terms()[0].second < Rational(0)) {
                    result.constant = -result.constant;
                }
            }
        } catch (const std::runtime_error&) {
            if (!product.terms().empty() && !poly.terms().empty()) {
                result.constant = poly.terms()[0].second / product.terms()[0].second;
            }
        }
    }
}

/**
 * @brief 检测是否为线性多项式(主变量次数 <= 1)
 */
bool is_linear(const MultiPoly& poly, const std::string& main_var)
{
    return poly.degree(main_var) <= 1;
}

/**
 * @brief 检测差平方模式 a^2 - b^2
 *
 * 判断多项式是否恰好有两项,一正一负,且两项的单项式各分量指数均为偶数.
 * 若匹配,返回 (a, b),其中 a^2 为正项,b^2 为负项.
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (a, b),否则返回 nullopt
 */
std::optional<std::pair<MultiPoly, MultiPoly>>
detect_difference_of_squares(const MultiPoly& poly)
{
    if (poly.num_terms() != 2) { return std::nullopt; }

    const auto& terms = poly.terms();
    const auto& vars = poly.variables();

    /// 确定哪项为正,哪项为负
    int pos_idx = -1, neg_idx = -1;
    if (terms[0].second > Rational(0) && terms[1].second < Rational(0)) {
        pos_idx = 0; neg_idx = 1;
    } else if (terms[1].second > Rational(0) && terms[0].second < Rational(0)) {
        pos_idx = 1; neg_idx = 0;
    } else {
        return std::nullopt;
    }

    /// 检查系数绝对值相等(都应为完全平方的系数)
    /// 对于差平方 a^2 - b^2,正项系数和负项系数的绝对值不必相等
    /// 但两项的系数本身必须是完全平方数
    Rational pos_coeff = terms[pos_idx].second;
    Rational neg_coeff = -terms[neg_idx].second; // 取绝对值

    /// 检查系数是否为完全平方有理数
    /// 分子和分母都必须是完全平方数
    auto positive_root = rational_square_root(pos_coeff);
    if (!positive_root) { return std::nullopt; }
    auto negative_root = rational_square_root(neg_coeff);
    if (!negative_root) { return std::nullopt; }

    /// 检查两项的单项式各分量指数均为偶数
    const Monomial& pos_mono = terms[pos_idx].first;
    const Monomial& neg_mono = terms[neg_idx].first;

    for (size_t i = 0; i < pos_mono.size(); ++i) {
        if (pos_mono[i] % 2 != 0) { return std::nullopt; }
    }
    for (size_t i = 0; i < neg_mono.size(); ++i) {
        if (neg_mono[i] % 2 != 0) { return std::nullopt; }
    }

    /// 构造 a:正项单项式指数减半,系数为 sqrt(pos_coeff)
    Monomial a_mono(vars.size(), 0);
    for (size_t i = 0; i < pos_mono.size(); ++i) {
        a_mono[i] = pos_mono[i] / 2;
    }
    const Rational& a_coeff = *positive_root;
    std::vector<MultiPoly::Term> a_terms = {{a_mono, a_coeff}};
    MultiPoly a(std::move(a_terms), vars);

    /// 构造 b:负项单项式指数减半,系数为 sqrt(neg_coeff)
    Monomial b_mono(vars.size(), 0);
    for (size_t i = 0; i < neg_mono.size(); ++i) {
        b_mono[i] = neg_mono[i] / 2;
    }
    const Rational& b_coeff = *negative_root;
    std::vector<MultiPoly::Term> b_terms = {{b_mono, b_coeff}};
    MultiPoly b(std::move(b_terms), vars);

    return std::make_pair(a, b);
}

/**
 * @brief 检测二项式幂模式 xⁿ - yⁿ
 *
 * 判断多项式是否为两个不同变量的纯幂之差,系数分别为 +1 和 -1,
 * 且幂次相同.若匹配,返回 (变量1, 变量2, 幂次).
 *
 * @param[in] poly 输入多项式
 * @return 若匹配则返回 (var1, var2, n),否则返回 nullopt
 */
std::optional<std::tuple<std::string, std::string, int>>
detect_binomial_power(const MultiPoly& poly)
{
    if (poly.num_terms() != 2) { return std::nullopt; }

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

    /// 检查每项的单项式是否为单变量的纯幂(恰好一个非零指数)
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
        if (nonzero_count != 1 || exponent <= 1) { return std::nullopt; }
        return std::make_pair(var_index, exponent);
    };

    auto pos_info = get_single_var_power(terms[pos_idx].first);
    if (!pos_info) { return std::nullopt; }
    auto neg_info = get_single_var_power(terms[neg_idx].first);
    if (!neg_info) { return std::nullopt; }

    auto [pos_var_idx, pos_exp] = *pos_info;
    auto [neg_var_idx, neg_exp] = *neg_info;

    /// 幂次必须相同,变量必须不同
    if (pos_exp != neg_exp) { return std::nullopt; }
    if (pos_var_idx == neg_var_idx) { return std::nullopt; }

    return std::make_tuple(vars[pos_var_idx], vars[neg_var_idx], pos_exp);
}

/**
 * @brief 提取公因子单项式
 *
 * 计算所有项单项式的逐分量最小值(GCD 单项式),将其从每项中减去.
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
            if (val < gcd_mono[i]) { gcd_mono[i] = val; }
        }
    }

    /// 检查是否为平凡(全零)
    bool is_trivial = true;
    for (size_t i = 0; i < gcd_mono.size(); ++i) {
        if (gcd_mono[i] != 0) { is_trivial = false; break; }
    }
    if (is_trivial) { return {identity, poly}; }

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
 * 对齐次二元多项式 f(x, y),令 y=1 得到一元多项式 f(x, 1),
 * 对其进行一元因式分解,然后将每个因子重新齐次化.
 *
 * @param[in] poly 输入多项式（须为齐次二元）
 * @param[in,out] context 调用方的计算预算和取消状态
 * @return 不适用时成功返回 nullopt；分解错误原样传播
 */
Result<std::optional<MultiFactorResult>>
factor_homogeneous_bivariate(
    const MultiPoly& poly,
    ComputationContext& context)
{
    if (!poly.is_homogeneous() || poly.num_vars() != 2) {
        return std::optional<MultiFactorResult>{};
    }

    const auto& vars = poly.variables();
    std::string var_x = vars[0];
    std::string var_y = vars[1];

    /// 去齐次化:令 y = 1,得到关于 x 的一元多项式
    MultiPoly dehomogenized = poly.eval(var_y, Rational(1));

    /// 转换为一元多项式
    Polynomial<Rational> uni_poly;
    try {
        uni_poly = dehomogenized.to_univariate();
    } catch (const std::invalid_argument&) {
        return std::optional<MultiFactorResult>{};
    }

    /// 一元因式分解
    auto checked_factors =
        factor_univariate_bridge_checked(uni_poly, context);
    if (!checked_factors) {
        return Result<std::optional<MultiFactorResult>>::failure(
            checked_factors.error());
    }
    if (checked_factors.value().completeness != Completeness::Complete) {
        return std::optional<MultiFactorResult>{};
    }
    std::vector<Polynomial<Rational>> uni_factors =
        std::move(checked_factors.value().value);
    if (uni_factors.size() <= 1) {
        return std::optional<MultiFactorResult>{};
    }

    /// 重新齐次化每个因子
    MultiFactorResult result;
    result.constant = poly.numeric_content();

    for (const auto& uf : uni_factors) {
        auto homo_factor = homogenize_factor(uf, vars);
        result.factors.push_back(std::move(homo_factor));
        result.multiplicities.push_back(1);
    }

    adjust_homogeneous_constant(poly, result);

    return std::optional<MultiFactorResult>(std::move(result));
}

} // namespace detail
} // namespace LMCAS
