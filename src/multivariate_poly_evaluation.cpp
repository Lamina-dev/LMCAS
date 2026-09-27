#include "multivariate_poly.hpp"
#include <stdexcept>

namespace LMCAS {

/**
 * @brief 将指定变量代入有理数值，返回降维后的多项式
 *
 * 对每一项，将 val 的对应指数次幂乘入系数，然后从单项式中移除该变量维度。
 * 结果多项式的变量数减一（若变量存在于列表中）。
 *
 * @param[in] var 待代入的变量名
 * @param[in] val 代入值
 * @return 代入后的多项式（变量数减一或不变）
 */
MultiPoly MultiPoly::eval(const std::string& var, const Rational& val) const
{
    int var_idx = -1;
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (vars_[i] == var) {
            var_idx = static_cast<int>(i);
            break;
        }
    }
    if (var_idx < 0) {
        return *this;
    }

    std::vector<std::string> new_vars;
    new_vars.reserve(vars_.size() - 1);
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (static_cast<int>(i) != var_idx) {
            new_vars.push_back(vars_[i]);
        }
    }

    std::vector<Term> new_terms;
    new_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        int exp = (static_cast<size_t>(var_idx) < term.first.size())
                      ? term.first[var_idx] : 0;

        Rational new_coeff = term.second;
        if (exp > 0) {
            if (val.is_zero()) {
                continue;
            }
            new_coeff = new_coeff * val.power(BigInt(exp));
        }
        Monomial reduced_mono;
        reduced_mono.reserve(vars_.size() - 1);
        for (size_t i = 0; i < term.first.size(); ++i) {
            if (static_cast<int>(i) != var_idx) {
                reduced_mono.push_back(term.first[i]);
            }
        }

        new_terms.emplace_back(std::move(reduced_mono), std::move(new_coeff));
    }

    return MultiPoly(std::move(new_terms), std::move(new_vars), order_);
}

/**
 * @brief 将多个变量同时代入有理数值
 *
 * 依次对每个 (var, val) 对调用单变量 eval，逐步降维。
 *
 * @param[in] substitution 变量名到值的映射
 * @return 代入后的多项式
 */
MultiPoly MultiPoly::eval(const std::map<std::string, Rational>& substitution) const
{
    MultiPoly result = *this;
    for (const auto& [var, val] : substitution) {
        result = result.eval(var, val);
    }
    return result;
}


/**
 * @brief 将多元多项式转换为一元 Polynomial<Rational>
 *
 * 仅当多项式实际只含一个有效变量（或为常数/零多项式）时有效。
 * 将单项式中该变量的指数映射为系数向量的下标。
 *
 * @return 等价的一元多项式
 * @throw std::invalid_argument 含多个有效变量时抛出
 */
Polynomial<Rational> MultiPoly::to_univariate() const
{
    if (!is_univariate()) {
        throw std::invalid_argument(
            "MultiPoly::to_univariate: polynomial is not univariate");
    }

    if (terms_.empty()) {
        return Polynomial<Rational>(vars_.empty() ? "x" : vars_[0]);
    }
    int active_idx = -1;
    for (size_t vi = 0; vi < vars_.size(); ++vi) {
        for (const auto& term : terms_) {
            if (vi < term.first.size() && term.first[vi] != 0) {
                active_idx = static_cast<int>(vi);
                break;
            }
        }
        if (active_idx >= 0) break;
    }
    std::string var_name = vars_.empty() ? "x" : vars_[0];
    if (active_idx < 0) {
        Rational c = terms_[0].second;
        return Polynomial<Rational>(c, var_name);
    }

    var_name = vars_[static_cast<size_t>(active_idx)];

    int max_deg = 0;
    for (const auto& term : terms_) {
        int exp = term.first[static_cast<size_t>(active_idx)];
        if (exp > max_deg) max_deg = exp;
    }
    std::vector<Rational> coeffs(static_cast<size_t>(max_deg + 1), Rational(0));
    for (const auto& term : terms_) {
        int exp = term.first[static_cast<size_t>(active_idx)];
        coeffs[static_cast<size_t>(exp)] = term.second;
    }

    return Polynomial<Rational>(coeffs, var_name);
}

/**
 * @brief 从一元多项式构造多元多项式
 *
 * 将 Polynomial<Rational> 的每个非零系数映射为含单变量单项式的项。
 *
 * @param[in] poly 一元多项式
 * @param[in] var 对应的变量名
 * @return 等价的 MultiPoly 表示
 */
MultiPoly MultiPoly::from_univariate(const Polynomial<Rational>& poly,
                                     const std::string& var)
{
    std::vector<std::string> vars = {var};

    if (poly.is_zero()) {
        return MultiPoly(std::vector<Term>{}, vars);
    }

    std::vector<Term> terms;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        if (!poly.coeffs[i].is_zero()) {
            Monomial mono = {static_cast<int>(i)};
            terms.emplace_back(std::move(mono), poly.coeffs[i]);
        }
    }

    return MultiPoly(std::move(terms), std::move(vars));
}

}
