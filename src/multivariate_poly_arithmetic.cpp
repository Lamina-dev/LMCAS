#include "multivariate_poly.hpp"
#include <algorithm>
#include <optional>
#include <stdexcept>

namespace LMCAS {

void MultiPoly::require_compatible_variables(const MultiPoly& other) const
{
    if (vars_ != other.vars_ && !is_zero() && !other.is_zero()) {
        throw std::invalid_argument("MultiPoly: variable list mismatch");
    }
}

/**
 * @brief 多项式相加并合并同类项、去零、排序。
 * 非零操作数须共享变量表；零元素采用非零侧的变量表和单项式序。
 * @param[in] other 加数。
 * @return 和多项式。
 */
MultiPoly MultiPoly::operator+(const MultiPoly& other) const
{
    if (other.is_zero()) return *this;
    if (is_zero()) return other;
    require_compatible_variables(other);

    std::vector<Term> result_terms;
    result_terms.reserve(terms_.size() + other.terms_.size());
    result_terms.insert(result_terms.end(), terms_.begin(), terms_.end());
    result_terms.insert(result_terms.end(), other.terms_.begin(), other.terms_.end());

    MultiPoly result;
    result.terms_ = std::move(result_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    result.normalize();
    return result;
}

/**
 * @brief 将 other 的各项系数取负后与自身相加。
 * @param[in] other 减数。
 * @return 差多项式。
 */
MultiPoly MultiPoly::operator-(const MultiPoly& other) const
{
    if (other.is_zero()) return *this;
    if (is_zero()) return -other;
    require_compatible_variables(other);

    std::vector<Term> result_terms;
    result_terms.reserve(terms_.size() + other.terms_.size());
    result_terms.insert(result_terms.end(), terms_.begin(), terms_.end());

    for (const auto& term : other.terms_) {
        result_terms.emplace_back(term.first, -term.second);
    }

    MultiPoly result;
    result.terms_ = std::move(result_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    result.normalize();
    return result;
}

/**
 * @brief 多项式相乘，指数逐分量相加、系数相乘，再合并同类项。
 * @param[in] other 乘数。
 * @return 积多项式。
 */
MultiPoly MultiPoly::operator*(const MultiPoly& other) const
{
    if (terms_.empty() || other.terms_.empty()) {
        MultiPoly zero;
        const MultiPoly& ring = is_zero() && !other.is_zero() ? other : *this;
        zero.vars_ = ring.vars_;
        zero.order_ = ring.order_;
        return zero;
    }
    require_compatible_variables(other);

    std::vector<Term> result_terms;
    result_terms.reserve(terms_.size() * other.terms_.size());

    size_t n_vars = vars_.size();

    for (const auto& term_a : terms_) {
        for (const auto& term_b : other.terms_) {
            Monomial product_mono(n_vars, 0);
            for (size_t i = 0; i < n_vars; ++i) {
                int exp_a = (i < term_a.first.size()) ? term_a.first[i] : 0;
                int exp_b = (i < term_b.first.size()) ? term_b.first[i] : 0;
                product_mono[i] = exp_a + exp_b;
            }
            Rational product_coeff = term_a.second * term_b.second;
            result_terms.emplace_back(std::move(product_mono), std::move(product_coeff));
        }
    }

    MultiPoly result;
    result.terms_ = std::move(result_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    result.normalize();
    return result;
}

/**
 * @brief 将多项式各项系数取负。
 * @return 取负后的多项式。
 */
MultiPoly MultiPoly::operator-() const
{
    std::vector<Term> negated_terms;
    negated_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        negated_terms.emplace_back(term.first, -term.second);
    }

    MultiPoly result;
    result.terms_ = std::move(negated_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    return result;
}

/**
 * @brief 比较多项式是否相等。
 * 零多项式彼此相等；非零多项式须共享变量表且各单项式系数相同，排序方式不限。
 * @param[in] other 比较对象。
 * @return 相等时返回 true。
 */
bool MultiPoly::operator==(const MultiPoly& other) const
{
    if (terms_.size() != other.terms_.size()) return false;
    if (is_zero()) return true;
    if (vars_ != other.vars_) return false;
    if (order_ != other.order_) {
        const MonomialOrder order(other.order_);
        for (const auto& term : terms_) {
            auto found = std::lower_bound(
                other.terms_.begin(), other.terms_.end(), term.first,
                [&order](const Term& candidate, const Monomial& monomial) {
                    return order(candidate.first, monomial);
                });
            if (found == other.terms_.end() ||
                found->first != term.first || found->second != term.second) {
                return false;
            }
        }
        return true;
    }
    for (size_t i = 0; i < terms_.size(); ++i) {
        if (terms_[i].first != other.terms_[i].first) return false;
        if (terms_[i].second != other.terms_[i].second) return false;
    }
    return true;
}

/**
 * @brief 比较多项式是否不等。
 * @param[in] other 比较对象。
 * @return 不等时返回 true。
 */
bool MultiPoly::operator!=(const MultiPoly& other) const
{
    return !(*this == other);
}

/**
 * @brief 各项系数乘以标量后规范化；零标量生成零多项式。
 * @param[in] scalar 有理数标量。
 * @return 标量乘积。
 */
MultiPoly MultiPoly::operator*(const Rational& scalar) const
{
    if (scalar.is_zero()) {
        MultiPoly zero;
        zero.vars_ = vars_;
        zero.order_ = order_;
        return zero;
    }

    std::vector<Term> scaled_terms;
    scaled_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        scaled_terms.emplace_back(term.first, term.second * scalar);
    }

    MultiPoly result;
    result.terms_ = std::move(scaled_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    result.normalize();
    return result;
}


/**
 * @brief 通过首项消去计算多元多项式的精确商。
 * 每步以余式首项除以除数首项得到商项，再减去商项与除数的乘积；
 * 除数首项单项式无法整除余式首项单项式时，判定为非精确除法。
 * @param[in] divisor 除数多项式。
 * @return 商多项式。
 * @throw std::runtime_error 除数为零或除法不精确。
 * @throw std::invalid_argument 两个非零操作数的变量表或排列不同。
 */
MultiPoly MultiPoly::exact_div(const MultiPoly& divisor) const
{
    if (divisor.is_zero()) {
        throw std::runtime_error("exact division by zero polynomial");
    }

    if (is_zero()) {
        MultiPoly zero;
        zero.vars_ = divisor.vars_;
        zero.order_ = divisor.order_;
        return zero;
    }
    require_compatible_variables(divisor);

    std::optional<MultiPoly> reordered_divisor;
    if (order_ != divisor.order_) {
        reordered_divisor.emplace(divisor);
        reordered_divisor->order_ = order_;
        reordered_divisor->normalize();
    }
    const MultiPoly& ordered_divisor = reordered_divisor ? *reordered_divisor : divisor;
    const Monomial& divisor_lt_mono = ordered_divisor.terms_[0].first;
    const Rational& divisor_lt_coeff = ordered_divisor.terms_[0].second;

    size_t n_vars = vars_.size();

    std::vector<Term> quotient_terms;
    MultiPoly remainder = *this;

    while (!remainder.is_zero()) {
        const Monomial& rem_lt_mono = remainder.terms_[0].first;
        const Rational& rem_lt_coeff = remainder.terms_[0].second;

        if (!divides_monomial(divisor_lt_mono, rem_lt_mono)) {
            throw std::runtime_error("exact division failed");
        }

        Monomial q_mono(n_vars, 0);
        for (size_t i = 0; i < n_vars; ++i) {
            int rem_exp = (i < rem_lt_mono.size()) ? rem_lt_mono[i] : 0;
            int div_exp = (i < divisor_lt_mono.size()) ? divisor_lt_mono[i] : 0;
            q_mono[i] = rem_exp - div_exp;
        }

        Rational q_coeff = rem_lt_coeff / divisor_lt_coeff;

        quotient_terms.emplace_back(q_mono, q_coeff);

        MultiPoly q_term_poly(std::vector<Term>{{q_mono, q_coeff}}, vars_, order_);
        remainder = remainder - (q_term_poly * ordered_divisor);
    }

    return MultiPoly(std::move(quotient_terms), vars_, order_);
}

/**
 * @brief 提取所有系数的有理数 GCD 作为数值内容。
 * 有理数 GCD 满足 gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)。
 * @return 非零多项式返回正有理数；零多项式返回 Rational(0)。
 */
Rational MultiPoly::numeric_content() const
{
    if (terms_.empty()) {
        return Rational(0);
    }

    BigInt num_gcd = terms_[0].second.get_numerator().abs();
    BigInt den_lcm = terms_[0].second.get_denominator();

    for (size_t i = 1; i < terms_.size(); ++i) {
        BigInt num_i = terms_[i].second.get_numerator().abs();
        BigInt den_i = terms_[i].second.get_denominator();

        num_gcd = BigInt::gcd(num_gcd, num_i);
        BigInt den_gcd = BigInt::gcd(den_lcm, den_i);
        den_lcm = den_lcm * den_i / den_gcd;
    }

    if (!num_gcd) {
        return Rational(0);
    }

    return Rational(num_gcd, den_lcm);
}

/**
 * @brief 除以数值内容并调整符号，使首项系数为正。
 * 零多项式原样返回；内容为 1 时仅调整首项符号。
 * @return 零多项式或数值内容为 1、首项系数为正的本原多项式。
 */
MultiPoly MultiPoly::make_primitive() const
{
    if (terms_.empty()) {
        return *this;
    }

    Rational content = numeric_content();

    if (content.is_zero() || content == Rational(1)) {
        if (!terms_.empty() && terms_[0].second < Rational(0)) {
            return (*this) * Rational(-1);
        }
        return *this;
    }

    std::vector<Term> prim_terms;
    prim_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        prim_terms.emplace_back(term.first, term.second / content);
    }

    if (!prim_terms.empty() && prim_terms[0].second < Rational(0)) {
        for (auto& term : prim_terms) {
            term.second = -term.second;
        }
    }

    MultiPoly result;
    result.terms_ = std::move(prim_terms);
    result.vars_ = vars_;
    result.order_ = order_;
    return result;
}

}
