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

}
