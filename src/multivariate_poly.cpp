/**
 * @file multivariate_poly.cpp
 * @brief 多元多项式类 MultiPoly 的实现。
 *
 * 本文件实现 MultiPoly 的构造函数和规范化逻辑。
 * 规范化保证内部表示满足以下不变量：
 * - 项按 MonomialOrder 严格降序排列
 * - 无零系数项
 * - 无重复单项式（同类项已合并）
 * - 所有单项式长度等于变量数
 */

#include "multivariate_poly.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace LMCAS {


/**
 * @brief 构造零多项式
 *
 * 零多项式的项列表为空，变量列表为空，默认使用 GrevLex 序。
 */
MultiPoly::MultiPoly()
    : terms_(), vars_(), order_(MonomialOrderType::GrevLex)
{
}

/**
 * @brief 从项列表和变量名列表构造多元多项式
 *
 * 构造后自动执行规范化：合并同类项、去除零系数项、按单项式序排序。
 * 若项的单项式长度不足变量数，自动补零至正确长度。
 *
 * @param[in] terms 项列表，每项为 (单项式, 系数) 对
 * @param[in] vars 变量名列表，确定各分量对应的变量
 * @param[in] order 单项式序类型，默认分次逆字典序
 */
MultiPoly::MultiPoly(std::vector<Term> terms, std::vector<std::string> vars,
                     MonomialOrderType order)
    : terms_(std::move(terms)), vars_(std::move(vars)), order_(order)
{
    /// 确保所有单项式长度与变量数一致
    for (auto& term : terms_) {
        if (term.first.size() < vars_.size()) {
            term.first.resize(vars_.size(), 0);
        }
    }
    normalize();
}

/**
 * @brief 从常数构造多元多项式
 *
 * 若常数为零，构造零多项式（空项列表）。
 * 否则构造含单个全零单项式的常数多项式。
 *
 * @param[in] constant 常数值
 * @param[in] vars 变量名列表
 */
MultiPoly::MultiPoly(const Rational& constant, const std::vector<std::string>& vars)
    : terms_(), vars_(vars), order_(MonomialOrderType::GrevLex)
{
    if (!constant.is_zero()) {
        Monomial zero_mono(vars_.size(), 0);
        terms_.emplace_back(std::move(zero_mono), constant);
    }
}


/**
 * @internal
 * @brief 规范化多项式内部表示
 *
 * 执行三步操作：
 * 1. 按 MonomialOrder 降序排列所有项
 * 2. 合并相邻的同类项（单项式相同的项，系数相加）
 * 3. 移除零系数项
 *
 * 排序后同类项必然相邻，因此只需单次线性扫描即可完成合并。
 */
void MultiPoly::normalize()
{
    if (terms_.empty()) return;

    /// 步骤 1：按单项式序降序排列
    MonomialOrder cmp(order_);
    std::sort(terms_.begin(), terms_.end(),
              [&cmp](const Term& a, const Term& b) {
                  return cmp(a.first, b.first);
              });

    /// 步骤 2 & 3：合并同类项并移除零系数项
    std::vector<Term> merged;
    merged.reserve(terms_.size());

    for (size_t i = 0; i < terms_.size(); ) {
        Monomial current_mono = terms_[i].first;
        Rational coeff_sum = terms_[i].second;
        size_t j = i + 1;

        /// 合并所有具有相同单项式的连续项
        while (j < terms_.size() && terms_[j].first == current_mono) {
            coeff_sum = coeff_sum + terms_[j].second;
            ++j;
        }

        /// 仅保留非零系数项
        if (!coeff_sum.is_zero()) {
            merged.emplace_back(std::move(current_mono), std::move(coeff_sum));
        }

        i = j;
    }

    terms_ = std::move(merged);
}


bool MultiPoly::is_zero() const
{
    return terms_.empty();
}

bool MultiPoly::is_constant() const
{
    if (terms_.empty()) return true;
    if (terms_.size() != 1) return false;
    /// 检查唯一项的单项式是否全零
    for (int exp : terms_[0].first) {
        if (exp != 0) return false;
    }
    return true;
}

bool MultiPoly::is_univariate() const
{
    if (terms_.empty()) return true;
    /// 统计哪些变量出现了非零指数
    int active_count = 0;
    for (size_t vi = 0; vi < vars_.size(); ++vi) {
        bool active = false;
        for (const auto& term : terms_) {
            if (vi < term.first.size() && term.first[vi] != 0) {
                active = true;
                break;
            }
        }
        if (active) ++active_count;
        if (active_count > 1) return false;
    }
    return true;
}

bool MultiPoly::is_homogeneous() const
{
    if (terms_.empty()) return true;
    int deg = LMCAS::total_degree(terms_[0].first);
    for (size_t i = 1; i < terms_.size(); ++i) {
        if (LMCAS::total_degree(terms_[i].first) != deg) return false;
    }
    return true;
}

int MultiPoly::num_terms() const
{
    return static_cast<int>(terms_.size());
}

int MultiPoly::num_vars() const
{
    return static_cast<int>(vars_.size());
}

const std::vector<std::string>& MultiPoly::variables() const
{
    return vars_;
}

const std::vector<MultiPoly::Term>& MultiPoly::terms() const
{
    return terms_;
}

int MultiPoly::total_degree() const
{
    if (terms_.empty()) return -1;
    int max_deg = 0;
    for (const auto& term : terms_) {
        int deg = LMCAS::total_degree(term.first);
        if (deg > max_deg) max_deg = deg;
    }
    return max_deg;
}

int MultiPoly::degree(const std::string& var) const
{
    if (terms_.empty()) return -1;
    /// 找到变量在 vars_ 中的索引
    int var_idx = -1;
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (vars_[i] == var) {
            var_idx = static_cast<int>(i);
            break;
        }
    }
    if (var_idx < 0) return 0; // 变量不在列表中，次数为 0

    int max_exp = 0;
    for (const auto& term : terms_) {
        if (static_cast<size_t>(var_idx) < term.first.size()) {
            if (term.first[var_idx] > max_exp) {
                max_exp = term.first[var_idx];
            }
        }
    }
    return max_exp;
}

/**
 * @brief 获取关于指定变量的首项系数
 *
 * 将多项式视为 var 的一元多项式，找到 var 的最高次数，
 * 收集所有具有该最高次数的项，移除 var 维度后返回辅助变量的多项式。
 *
 * @param[in] var 主变量名
 * @return 首项系数多项式（在剩余变量上）
 */
MultiPoly MultiPoly::leading_coeff(const std::string& var) const
{
    if (terms_.empty()) {
        return MultiPoly();
    }

    /// 找到变量在 vars_ 中的索引
    int var_idx = -1;
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (vars_[i] == var) {
            var_idx = static_cast<int>(i);
            break;
        }
    }

    /// 变量不在列表中：整个多项式就是"系数"
    if (var_idx < 0) {
        return *this;
    }

    /// 找到 var 的最高指数
    int max_exp = 0;
    for (const auto& term : terms_) {
        if (static_cast<size_t>(var_idx) < term.first.size()) {
            if (term.first[var_idx] > max_exp) {
                max_exp = term.first[var_idx];
            }
        }
    }

    /// 构建剩余变量列表（移除 var）
    std::vector<std::string> remaining_vars;
    remaining_vars.reserve(vars_.size() - 1);
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (static_cast<int>(i) != var_idx) {
            remaining_vars.push_back(vars_[i]);
        }
    }

    /// 收集具有最高指数的项，移除 var 维度
    std::vector<Term> lc_terms;
    for (const auto& term : terms_) {
        int exp_val = (static_cast<size_t>(var_idx) < term.first.size())
                          ? term.first[var_idx] : 0;
        if (exp_val == max_exp) {
            Monomial reduced_mono;
            reduced_mono.reserve(vars_.size() - 1);
            for (size_t i = 0; i < term.first.size(); ++i) {
                if (static_cast<int>(i) != var_idx) {
                    reduced_mono.push_back(term.first[i]);
                }
            }
            lc_terms.emplace_back(std::move(reduced_mono), term.second);
        }
    }

    return MultiPoly(std::move(lc_terms), std::move(remaining_vars), order_);
}


/**
 * @brief 多项式加法
 *
 * 将两个多项式的项列表合并后规范化（合并同类项、去零、排序）。
 * 使用 this 的变量列表和单项式序。
 *
 * @param[in] other 加数
 * @return 和多项式
 */
MultiPoly MultiPoly::operator+(const MultiPoly& other) const
{
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
 * @brief 多项式减法
 *
 * 将 other 的各项系数取负后与 this 相加。
 *
 * @param[in] other 减数
 * @return 差多项式
 */
MultiPoly MultiPoly::operator-(const MultiPoly& other) const
{
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
 * @brief 多项式乘法
 *
 * 对两个多项式的每对项执行：
 * - 单项式：各分量指数逐元素相加
 * - 系数：有理数相乘
 * 结果经规范化合并同类项。
 *
 * @param[in] other 乘数
 * @return 积多项式
 */
MultiPoly MultiPoly::operator*(const MultiPoly& other) const
{
    if (terms_.empty() || other.terms_.empty()) {
        MultiPoly zero;
        zero.vars_ = vars_;
        zero.order_ = order_;
        return zero;
    }

    std::vector<Term> result_terms;
    result_terms.reserve(terms_.size() * other.terms_.size());

    size_t n_vars = vars_.size();

    for (const auto& term_a : terms_) {
        for (const auto& term_b : other.terms_) {
            /// 单项式逐元素相加
            Monomial product_mono(n_vars, 0);
            for (size_t i = 0; i < n_vars; ++i) {
                int exp_a = (i < term_a.first.size()) ? term_a.first[i] : 0;
                int exp_b = (i < term_b.first.size()) ? term_b.first[i] : 0;
                product_mono[i] = exp_a + exp_b;
            }
            /// 系数相乘
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
 * @brief 一元取负
 *
 * 将所有项的系数取负。
 *
 * @return 取负后的多项式
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
    /// 无需 normalize：仅系数取负不改变排序和唯一性
    return result;
}

/**
 * @brief 判等
 *
 * 两个多项式相等当且仅当规范化后的项列表完全相同
 * （相同单项式和相同系数，相同顺序）。
 *
 * @param[in] other 比较对象
 * @return 相等返回 true
 */
bool MultiPoly::operator==(const MultiPoly& other) const
{
    if (terms_.size() != other.terms_.size()) return false;
    for (size_t i = 0; i < terms_.size(); ++i) {
        if (terms_[i].first != other.terms_[i].first) return false;
        if (terms_[i].second != other.terms_[i].second) return false;
    }
    return true;
}

/**
 * @brief 判不等
 *
 * @param[in] other 比较对象
 * @return 不等返回 true
 */
bool MultiPoly::operator!=(const MultiPoly& other) const
{
    return !(*this == other);
}

/**
 * @brief 标量乘法
 *
 * 将每项系数乘以标量，然后规范化（若标量为零则结果为零多项式）。
 *
 * @param[in] scalar 有理数标量
 * @return 各项系数乘以 scalar 后的多项式
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
 * @brief 多元多项式精确除法
 *
 * 假设 divisor 整除 *this，通过迭代首项消去法计算商：
 * 每步取余式首项除以除数首项得到商项，再从余式中减去该商项与除数的乘积。
 * 若过程中发现除数首项单项式不整除余式首项单项式，则除法不精确。
 *
 * @param[in] divisor 除数多项式
 * @return 商多项式
 * @throw std::runtime_error 除数为零或除法不精确时抛出
 */
MultiPoly MultiPoly::exact_div(const MultiPoly& divisor) const
{
    if (divisor.is_zero()) {
        throw std::runtime_error("exact division by zero polynomial");
    }

    /// 被除数为零，商为零
    if (is_zero()) {
        MultiPoly zero;
        zero.vars_ = vars_;
        zero.order_ = order_;
        return zero;
    }

    /// 除数首项（排序后第一项即为首项）
    const Monomial& divisor_lt_mono = divisor.terms_[0].first;
    const Rational& divisor_lt_coeff = divisor.terms_[0].second;

    size_t n_vars = vars_.size();

    std::vector<Term> quotient_terms;
    MultiPoly remainder = *this;

    while (!remainder.is_zero()) {
        /// 余式首项
        const Monomial& rem_lt_mono = remainder.terms_[0].first;
        const Rational& rem_lt_coeff = remainder.terms_[0].second;

        /// 检查除数首项单项式是否整除余式首项单项式
        if (!divides_monomial(divisor_lt_mono, rem_lt_mono)) {
            throw std::runtime_error("exact division failed");
        }

        /// 计算商项单项式：逐分量相减
        Monomial q_mono(n_vars, 0);
        for (size_t i = 0; i < n_vars; ++i) {
            int rem_exp = (i < rem_lt_mono.size()) ? rem_lt_mono[i] : 0;
            int div_exp = (i < divisor_lt_mono.size()) ? divisor_lt_mono[i] : 0;
            q_mono[i] = rem_exp - div_exp;
        }

        /// 计算商项系数
        Rational q_coeff = rem_lt_coeff / divisor_lt_coeff;

        /// 记录商项
        quotient_terms.emplace_back(q_mono, q_coeff);

        /// 构造商项多项式并从余式中减去 (商项 * 除数)
        MultiPoly q_term_poly(std::vector<Term>{{q_mono, q_coeff}}, vars_, order_);
        remainder = remainder - (q_term_poly * divisor);
    }

    return MultiPoly(std::move(quotient_terms), vars_, order_);
}


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
    /// 找到变量在 vars_ 中的索引
    int var_idx = -1;
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (vars_[i] == var) {
            var_idx = static_cast<int>(i);
            break;
        }
    }

    /// 变量不在列表中，返回自身不变
    if (var_idx < 0) {
        return *this;
    }

    /// 构建新变量列表（移除 var）
    std::vector<std::string> new_vars;
    new_vars.reserve(vars_.size() - 1);
    for (size_t i = 0; i < vars_.size(); ++i) {
        if (static_cast<int>(i) != var_idx) {
            new_vars.push_back(vars_[i]);
        }
    }

    /// 对每一项：系数乘以 val^(该变量的指数)，移除该维度
    std::vector<Term> new_terms;
    new_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        int exp = (static_cast<size_t>(var_idx) < term.first.size())
                      ? term.first[var_idx] : 0;

        /// 计算 val^exp * coefficient
        Rational new_coeff = term.second;
        if (exp > 0) {
            if (val.is_zero()) {
                /// val=0 且 exp>0 → 该项系数变为 0，跳过
                continue;
            }
            new_coeff = new_coeff * val.power(BigInt(exp));
        }
        /// exp == 0 时 val^0 = 1，系数不变

        /// 构建移除 var 维度后的单项式
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

    /// 零多项式
    if (terms_.empty()) {
        return Polynomial<Rational>(vars_.empty() ? "x" : vars_[0]);
    }

    /// 找到有效变量的索引（指数非零的变量）
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

    /// 常数多项式（所有指数为零）
    std::string var_name = vars_.empty() ? "x" : vars_[0];
    if (active_idx < 0) {
        Rational c = terms_[0].second;
        return Polynomial<Rational>(c, var_name);
    }

    var_name = vars_[static_cast<size_t>(active_idx)];

    /// 确定最高次数
    int max_deg = 0;
    for (const auto& term : terms_) {
        int exp = term.first[static_cast<size_t>(active_idx)];
        if (exp > max_deg) max_deg = exp;
    }

    /// 构建系数向量（coeffs[i] = x^i 的系数）
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

    /// 零多项式
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


/**
 * @brief 提取数值内容（所有系数的有理数 GCD）
 *
 * 对所有项系数计算最大公约数。有理数 GCD 定义为：
 * gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)
 * 结果始终为正值。零多项式返回 Rational(0)。
 *
 * @return 所有系数的 GCD（正有理数）
 */
Rational MultiPoly::numeric_content() const
{
    if (terms_.empty()) {
        return Rational(0);
    }

    /// 从第一个系数的绝对值开始
    BigInt num_gcd = terms_[0].second.get_numerator().Abs();
    BigInt den_lcm = terms_[0].second.get_denominator();

    for (size_t i = 1; i < terms_.size(); ++i) {
        BigInt num_i = terms_[i].second.get_numerator().Abs();
        BigInt den_i = terms_[i].second.get_denominator();

        /// gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)
        num_gcd = BigInt::gcd(num_gcd, num_i);

        /// lcm(b, d) = b * d / gcd(b, d)
        BigInt den_gcd = BigInt::gcd(den_lcm, den_i);
        den_lcm = den_lcm * den_i / den_gcd;
    }

    if (!num_gcd) {
        return Rational(0);
    }

    return Rational(num_gcd, den_lcm);
}

/**
 * @brief 使多项式本原化（除以数值内容）
 *
 * 计算数值内容后，将所有系数除以该内容。
 * 若首项系数为负，则对所有系数取负以保证首项系数为正。
 * 零多项式或内容为 1 的多项式直接返回自身。
 *
 * @return 本原多项式（数值内容为 1，首项系数为正）
 */
MultiPoly MultiPoly::make_primitive() const
{
    if (terms_.empty()) {
        return *this;
    }

    Rational content = numeric_content();

    if (content.is_zero() || content == Rational(1)) {
        /// 内容为 1 时仍需检查首项系数符号
        if (!terms_.empty() && terms_[0].second < Rational(0)) {
            return (*this) * Rational(-1);
        }
        return *this;
    }

    /// 除以内容
    std::vector<Term> prim_terms;
    prim_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        prim_terms.emplace_back(term.first, term.second / content);
    }

    /// 若首项系数为负，取负使其为正
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


/**
 * @brief 转换为可读字符串表示
 *
 * 格式规则：
 * - 零多项式输出 "0"
 * - 常数多项式输出数值
 * - 各项以 " + " 或 " - " 连接
 * - 系数为 1 时省略（除非为常数项）
 * - 指数为 0 的变量省略，指数为 1 时省略指数
 *
 * @return 多项式的字符串形式
 */
std::string MultiPoly::to_string() const
{
    if (terms_.empty()) {
        return "0";
    }

    std::ostringstream oss;

    for (size_t i = 0; i < terms_.size(); ++i) {
        const Rational& coeff = terms_[i].second;
        const Monomial& mono = terms_[i].first;

        /// 判断是否为常数项（所有指数为零）
        bool is_const_term = true;
        for (size_t vi = 0; vi < mono.size(); ++vi) {
            if (mono[vi] != 0) {
                is_const_term = false;
                break;
            }
        }

        /// 确定系数的符号和绝对值
        bool negative = coeff < Rational(0);
        Rational abs_coeff = coeff.abs();

        /// 输出符号
        if (i == 0) {
            if (negative) oss << "-";
        } else {
            if (negative) {
                oss << " - ";
            } else {
                oss << " + ";
            }
        }

        /// 常数项：直接输出绝对值
        if (is_const_term) {
            oss << abs_coeff.to_string();
            continue;
        }

        /// 非常数项：系数不为 1 时输出系数
        bool coeff_is_one = (abs_coeff == Rational(1));
        if (!coeff_is_one) {
            oss << abs_coeff.to_string() << "*";
        }

        /// 输出变量部分
        bool first_var = true;
        for (size_t vi = 0; vi < mono.size() && vi < vars_.size(); ++vi) {
            if (mono[vi] == 0) continue;

            if (!first_var) {
                oss << "*";
            }
            first_var = false;

            oss << vars_[vi];
            if (mono[vi] != 1) {
                oss << "^" << mono[vi];
            }
        }
    }

    return oss.str();
}

} // namespace LMCAS
