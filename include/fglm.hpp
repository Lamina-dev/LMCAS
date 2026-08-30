/**
 * @file fglm.hpp
 * @brief FGLM 算法：零维理想 Gröbner 基在不同单项式序之间的转换。
 *
 * @see J.-C. Faugère, P. Gianni, D. Lazard, and T. Mora,
 *      “Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering,”
 *      Journal of Symbolic Computation 16(4), 1993.
 */
#pragma once

#include "monomial_order.hpp"
#include "rational.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <set>
#include <stdexcept>
#include <cassert>
#include <functional>

namespace lamina {

/** @brief FGLM 算法中使用的多项式表示，以 (单项式, 系数) 对列表存储 */
struct FGLMPoly {
    std::vector<std::pair<Monomial, Rational>> terms;  ///< 项列表
    size_t num_vars;  ///< 变量个数

    /** @brief 默认构造 */
    FGLMPoly() : num_vars(0) {}

    /**
     * @brief 构造指定变量数的零多项式
     * @param n 变量个数
     */
    explicit FGLMPoly(size_t n) : num_vars(n) {}

    /** @brief 判断是否为零多项式 */
    bool is_zero() const { return terms.empty(); }

    /** @brief 获取首项单项式 */
    Monomial LM() const {
        if (terms.empty()) return Monomial();
        return terms.front().first;
    }

    /** @brief 获取首项系数 */
    Rational LC() const {
        if (terms.empty()) return Rational(0);
        return terms.front().second;
    }

    /**
     * @brief 添加一项
     * @param m 单项式
     * @param c 系数（零系数不添加）
     */
    void add_term(const Monomial& m, const Rational& c) {
        if (!c.is_zero()) {
            terms.emplace_back(m, c);
        }
    }

    /**
     * @brief 按指定单项式序排序各项
     * @param order 单项式序
     */
    void sort_terms(const MonomialOrder& order) {
        std::sort(terms.begin(), terms.end(),
            [&order](const std::pair<Monomial, Rational>& a,
                     const std::pair<Monomial, Rational>& b) {
                return order(a.first, b.first);
            });
    }

    /** @brief 合并同类项并去除零项 */
    void normalize() {

        std::vector<std::pair<Monomial, Rational>> cleaned;
        for (auto& [m, c] : terms) {
            if (c.is_zero()) continue;
            if (!cleaned.empty() && cleaned.back().first == m) {
                cleaned.back().second = cleaned.back().second + c;
                if (cleaned.back().second.is_zero()) {
                    cleaned.pop_back();
                }
            } else {
                cleaned.emplace_back(m, c);
            }
        }
        terms = std::move(cleaned);
    }

    /**
     * @brief 对给定基进行多项式约化
     * @param basis Gröbner 基
     * @param order 单项式序
     * @return 约化后的余式
     */
    FGLMPoly reduce(const std::vector<FGLMPoly>& basis,
                    const MonomialOrder& order) const {
        FGLMPoly r(num_vars);
        FGLMPoly f = *this;
        f.sort_terms(order);
        f.normalize();

        int max_steps = 10000;
        while (!f.is_zero() && max_steps-- > 0) {
            bool reduced = false;
            for (const auto& g : basis) {
                if (g.is_zero()) continue;
                Monomial lm_f = f.LM();
                Monomial lm_g = g.LM();
                if (divides_monomial(lm_g, lm_f)) {

                    Monomial quot_mon(num_vars, 0);
                    for (size_t i = 0; i < num_vars; ++i) {
                        int fi = (i < lm_f.size()) ? lm_f[i] : 0;
                        int gi = (i < lm_g.size()) ? lm_g[i] : 0;
                        quot_mon[i] = fi - gi;
                    }
                    Rational quot_coeff = f.LC() / g.LC();

                    FGLMPoly subtracted(num_vars);
                    for (const auto& [gm, gc] : g.terms) {
                        Monomial product_mon(num_vars, 0);
                        for (size_t i = 0; i < num_vars; ++i) {
                            int qi = (i < quot_mon.size()) ? quot_mon[i] : 0;
                            int gmi = (i < gm.size()) ? gm[i] : 0;
                            product_mon[i] = qi + gmi;
                        }
                        subtracted.add_term(product_mon, quot_coeff * gc);
                    }

                    for (const auto& [sm, sc] : subtracted.terms) {
                        f.add_term(sm, -sc);
                    }
                    f.sort_terms(order);
                    f.normalize();
                    reduced = true;
                    break;
                }
            }
            if (!reduced) {

                r.add_term(f.LM(), f.LC());
                f.terms.erase(f.terms.begin());
                f.normalize();
            }
        }
        r.sort_terms(order);
        r.normalize();
        return r;
    }

    /**
     * @brief 从单项式构造多项式（系数为 1）
     * @param m 单项式
     * @param num_vars 变量个数
     * @return 仅含一项的多项式
     */
    static FGLMPoly from_monomial(const Monomial& m, size_t num_vars) {
        FGLMPoly p(num_vars);
        p.add_term(m, Rational(1));
        return p;
    }

    /**
     * @brief 构造零多项式
     * @param num_vars 变量个数
     * @return 零多项式
     */
    static FGLMPoly zero(size_t num_vars) {
        return FGLMPoly(num_vars);
    }
};

/**
 * @brief 计算多项式关于基的标准形（余式）
 * @param f 待约化多项式
 * @param basis Gröbner 基
 * @param order 单项式序
 * @return 标准形
 */
inline FGLMPoly normal_form(const FGLMPoly& f,
                            const std::vector<FGLMPoly>& basis,
                            const MonomialOrder& order) {
    return f.reduce(basis, order);
}

/**
 * @brief 判断理想是否为零维（每个变量都有纯幂次首项）
 * @param basis Gröbner 基
 * @param num_vars 变量个数
 * @return 零维返回 true
 */
inline bool is_zero_dimensional(const std::vector<FGLMPoly>& basis, size_t num_vars) {
    if (basis.empty()) return false;

    for (size_t var = 0; var < num_vars; ++var) {
        bool found = false;
        for (const auto& g : basis) {
            if (g.is_zero()) continue;
            Monomial lm = g.LM();

            bool is_pure = true;
            for (size_t i = 0; i < num_vars; ++i) {
                int exp_i = (i < lm.size()) ? lm[i] : 0;
                if (i == var) {
                    if (exp_i == 0) { is_pure = false; break; }
                } else {
                    if (exp_i != 0) { is_pure = false; break; }
                }
            }
            if (is_pure) { found = true; break; }
        }
        if (!found) return false;
    }
    return true;
}

/**
 * @brief 计算商空间维数（标准单项式个数）
 * @param basis Gröbner 基
 * @param num_vars 变量个数
 * @param max_degree 搜索的最大次数上限
 * @return 商空间维数；-1 表示当前搜索界内维数保持未知
 */
inline int quotient_dimension(const std::vector<FGLMPoly>& basis,
                              size_t num_vars,
                              int max_degree = 50) {
    if (basis.empty()) return -1;

    std::vector<Monomial> lead_mons;
    for (const auto& g : basis) {
        if (!g.is_zero()) {
            lead_mons.push_back(g.LM());
        }
    }

    int count = 0;
    const int limit = 1000;

    for (int d = 0; d <= max_degree && count <= limit; ++d) {

        std::vector<Monomial> degree_mons;
        Monomial current(num_vars, 0);

        std::function<void(size_t, int)> generate =
            [&](size_t var_idx, int remaining_deg) {
                if (var_idx == num_vars - 1) {
                    current[var_idx] = remaining_deg;
                    degree_mons.push_back(current);
                    current[var_idx] = 0;
                    return;
                }
                for (int e = 0; e <= remaining_deg; ++e) {
                    current[var_idx] = e;
                    generate(var_idx + 1, remaining_deg - e);
                }
                current[var_idx] = 0;
            };

        generate(0, d);

        for (const auto& m : degree_mons) {

            bool divisible = false;
            for (const auto& lm : lead_mons) {
                if (divides_monomial(lm, m)) {
                    divisible = true;
                    break;
                }
            }
            if (!divisible) {
                ++count;
                if (count > limit) return count;
            }
        }

        if (degree_mons.size() > 0) {
            bool any_standard = false;
            for (const auto& m : degree_mons) {
                bool divisible = false;
                for (const auto& lm : lead_mons) {
                    if (divides_monomial(lm, m)) { divisible = true; break; }
                }
                if (!divisible) { any_standard = true; break; }
            }
            if (!any_standard && d > 0) break;
        }
    }
    return count;
}

namespace detail {

/** @brief 单项式比较器（按给定序的逆序排列，用于最小堆） */
struct MonomialLessUnder {
    const MonomialOrder* order;
    MonomialLessUnder(const MonomialOrder* o) : order(o) {}
    bool operator()(const Monomial& a, const Monomial& b) const {
        return order->operator()(b, a);
    }
};

/** @brief 按单项式序递增枚举单项式的迭代器 */
class MonomialEnumerator {
public:
    /**
     * @brief 构造单项式枚举器
     * @param num_vars 变量个数
     * @param order 单项式序
     * @param max_deg 最大次数上限
     */
    MonomialEnumerator(size_t num_vars, const MonomialOrder& order, int max_deg = 100)
        : num_vars_(num_vars), order_(order), max_degree_(max_deg) {

        Monomial one(num_vars, 0);
        push(one);
    }

    /**
     * @brief 获取下一个单项式
     * @param out 输出单项式
     * @return 成功返回 true，枚举结束返回 false
     */
    bool next(Monomial& out) {
        if (heap_.empty()) return false;
        out = heap_.front();
        std::pop_heap(heap_.begin(), heap_.end(), cmp_);
        heap_.pop_back();

        int out_deg = total_degree(out);
        if (out_deg < max_degree_) {
            for (size_t i = 0; i < num_vars_; ++i) {
                Monomial succ = out;
                succ[i] += 1;
                push(succ);
            }
        }
        return true;
    }

private:
    size_t num_vars_;
    MonomialOrder order_;
    int max_degree_;
    std::vector<Monomial> heap_;
    std::set<Monomial> visited_;

    struct HeapCmp {
        const MonomialOrder* order;
        HeapCmp() : order(nullptr) {}
        HeapCmp(const MonomialOrder* o) : order(o) {}
        bool operator()(const Monomial& a, const Monomial& b) const {

            return order->operator()(a, b);
        }
    };
    HeapCmp cmp_{&order_};

    void push(const Monomial& m) {
        if (visited_.count(m)) return;
        visited_.insert(m);
        heap_.push_back(m);
        std::push_heap(heap_.begin(), heap_.end(), cmp_);
    }
};

/** @brief 增量式高斯消元器，用于检测向量线性相关性 */
class GaussianEliminator {
public:
    GaussianEliminator() {}

    /**
     * @brief 尝试将向量加入行空间
     * @param v 待加入的向量
     * @param combination 输出线性组合系数（当线性相关时有效）
     * @return true 表示线性无关（已加入），false 表示线性相关
     */
    bool add_vector(const std::vector<Rational>& v,
                    std::vector<Rational>& combination) {
        size_t n = v.size();

        for (auto& row : rows_) {
            row.resize(n, Rational(0));
        }

        std::vector<Rational> working = v;

        combination.assign(basis_count_, Rational(0));

        for (size_t i = 0; i < pivots_.size(); ++i) {
            size_t col = pivots_[i];
            if (col >= working.size()) continue;
            if (working[col].is_zero()) continue;

            Rational factor = working[col] / rows_[i][col];
            for (size_t j = 0; j < n; ++j) {
                if (j < working.size() && j < rows_[i].size()) {
                    working[j] = working[j] - factor * rows_[i][j];
                }
            }

            for (size_t j = 0; j < combination.size(); ++j) {
                combination[j] = combination[j] - factor * combinations_[i][j];
            }
        }

        size_t pivot_col = n;
        for (size_t j = 0; j < n; ++j) {
            if (!working[j].is_zero()) {
                pivot_col = j;
                break;
            }
        }

        if (pivot_col == n) {

            for (auto& c : combination) {
                c = -c;
            }
            return false;
        }

        Rational pivot_val = working[pivot_col];
        for (size_t j = 0; j < n; ++j) {
            working[j] = working[j] / pivot_val;
        }

        std::vector<Rational> comb_row(basis_count_ + 1, Rational(0));
        for (size_t j = 0; j < combination.size(); ++j) {
            comb_row[j] = combination[j] / pivot_val;
        }

        comb_row[basis_count_] = Rational(1) / pivot_val;

        rows_.push_back(working);
        pivots_.push_back(pivot_col);
        combinations_.push_back(comb_row);
        ++basis_count_;

        return true;
    }

    /** @brief 获取当前秩 */
    size_t rank() const { return basis_count_; }

private:
    std::vector<std::vector<Rational>> rows_;
    std::vector<size_t> pivots_;
    std::vector<std::vector<Rational>> combinations_;
    size_t basis_count_ = 0;
};

/**
 * @brief 将多项式转换为坐标向量（相对于单项式基）
 * @param p 多项式
 * @param all_monomials 单项式基（可能被扩展）
 * @param monomial_index 单项式到索引的映射
 * @return 坐标向量
 */
inline std::vector<Rational> poly_to_vector(
    const FGLMPoly& p,
    std::vector<Monomial>& all_monomials,
    std::map<Monomial, size_t>& monomial_index) {

    for (const auto& [m, c] : p.terms) {
        if (monomial_index.find(m) == monomial_index.end()) {
            size_t idx = all_monomials.size();
            all_monomials.push_back(m);
            monomial_index[m] = idx;
        }
    }

    std::vector<Rational> v(all_monomials.size(), Rational(0));
    for (const auto& [m, c] : p.terms) {
        v[monomial_index[m]] = c;
    }
    return v;
}

}

/**
 * @brief FGLM 算法：将零维理想的 Gröbner 基从源序转换到目标序
 * @param source_basis 源序下的 Gröbner 基
 * @param source_order 源单项式序
 * @param target_order 目标单项式序
 * @param num_vars 变量个数
 * @return 目标序下的 Gröbner 基
 * @throw std::runtime_error 理想非零维或维数过大时抛出
 * @see J.-C. Faugère, P. Gianni, D. Lazard, and T. Mora,
 *      “Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering,”
 *      Journal of Symbolic Computation 16(4), 1993.
 */
inline std::vector<FGLMPoly> fglm_convert(
    const std::vector<FGLMPoly>& source_basis,
    const MonomialOrder& source_order,
    const MonomialOrder& target_order,
    size_t num_vars) {

    if (source_basis.empty()) {
        return {};
    }

    if (!is_zero_dimensional(source_basis, num_vars)) {
        throw std::runtime_error(
            "FGLM: ideal is not zero-dimensional. "
            "FGLM only works for zero-dimensional ideals.");
    }

    int dim = quotient_dimension(source_basis, num_vars);
    if (dim <= 0 || dim > 1000) {
        throw std::runtime_error(
            "FGLM: quotient dimension is too large or could not be determined. "
            "dim = " + std::to_string(dim));
    }

    std::vector<FGLMPoly> target_basis;

    std::vector<Monomial> basis_monomials;

    std::vector<Monomial> standard_monomials;
    std::map<Monomial, size_t> monomial_index;
    {
        std::vector<Monomial> source_lms;
        for (const auto& g : source_basis) {
            if (!g.is_zero()) source_lms.push_back(g.LM());
        }

        for (int d = 0; d <= dim + 2 && (int)standard_monomials.size() < dim; ++d) {
            Monomial current(num_vars, 0);
            std::function<void(size_t, int)> gen = [&](size_t var_idx, int remaining) {
                if (var_idx == num_vars - 1) {
                    current[var_idx] = remaining;

                    bool is_std = true;
                    for (const auto& lm : source_lms) {
                        if (divides_monomial(lm, current)) { is_std = false; break; }
                    }
                    if (is_std) standard_monomials.push_back(current);
                    current[var_idx] = 0;
                    return;
                }
                for (int e = 0; e <= remaining; ++e) {
                    current[var_idx] = e;
                    gen(var_idx + 1, remaining - e);
                }
                current[var_idx] = 0;
            };
            gen(0, d);
        }

        for (size_t i = 0; i < standard_monomials.size(); ++i) {
            monomial_index[standard_monomials[i]] = i;
        }
    }

    auto nf_to_vector = [&](const FGLMPoly& nf) -> std::vector<Rational> {
        std::vector<Rational> v(standard_monomials.size(), Rational(0));
        for (const auto& [m, c] : nf.terms) {
            auto it = monomial_index.find(m);
            if (it != monomial_index.end()) {
                v[it->second] = c;
            }
        }
        return v;
    };

    detail::GaussianEliminator gauss;

    int max_mono_degree = dim + 5;
    detail::MonomialEnumerator enumerator(num_vars, target_order, max_mono_degree);

    std::vector<Monomial> target_lms;

    const int max_iterations = dim * 10 + 50;
    int iterations = 0;

    Monomial t;
    while (enumerator.next(t) && iterations < max_iterations) {
        ++iterations;

        bool skip = false;
        for (const auto& lm : target_lms) {
            if (divides_monomial(lm, t)) {
                skip = true;
                break;
            }
        }
        if (skip) continue;

        FGLMPoly t_poly = FGLMPoly::from_monomial(t, num_vars);
        FGLMPoly nf = normal_form(t_poly, source_basis, source_order);

        std::vector<Rational> v = nf_to_vector(nf);

        std::vector<Rational> combination;
        bool independent = gauss.add_vector(v, combination);

        if (independent) {

            basis_monomials.push_back(t);

            if (static_cast<int>(basis_monomials.size()) >= dim) {

            }
        } else {

            FGLMPoly new_elem(num_vars);
            new_elem.add_term(t, Rational(1));

            for (size_t i = 0; i < combination.size() && i < basis_monomials.size(); ++i) {
                if (!combination[i].is_zero()) {

                    new_elem.add_term(basis_monomials[i], -combination[i]);
                }
            }

            new_elem.sort_terms(target_order);
            new_elem.normalize();

            if (!new_elem.is_zero()) {
                Rational lc = new_elem.LC();
                if (!(lc == Rational(1))) {
                    for (auto& [m, c] : new_elem.terms) {
                        c = c / lc;
                    }
                }
            }

            if (!new_elem.is_zero()) {
                target_lms.push_back(new_elem.LM());
                target_basis.push_back(std::move(new_elem));
            }
        }

        if (static_cast<int>(basis_monomials.size()) >= dim) {

            bool complete = true;
            for (size_t var = 0; var < num_vars; ++var) {
                bool has_gen = false;
                for (const auto& lm : target_lms) {

                    bool is_pure = true;
                    for (size_t i = 0; i < num_vars; ++i) {
                        int exp_i = (i < lm.size()) ? lm[i] : 0;
                        if (i == var) {
                            if (exp_i == 0) { is_pure = false; break; }
                        } else {
                            if (exp_i != 0) { is_pure = false; break; }
                        }
                    }
                    if (is_pure) { has_gen = true; break; }
                }
                if (!has_gen) { complete = false; break; }
            }
            if (complete) break;
        }
    }

    return target_basis;
}

/**
 * @brief 将 GrevLex 序下的 Gröbner 基转换为 Lex 序
 * @param grevlex_basis GrevLex 序下的 Gröbner 基
 * @param num_vars 变量个数
 * @return Lex 序下的 Gröbner 基
 */
inline std::vector<FGLMPoly> grevlex_to_lex(
    const std::vector<FGLMPoly>& grevlex_basis,
    size_t num_vars) {
    return fglm_convert(
        grevlex_basis,
        MonomialOrder::grevlex(),
        MonomialOrder::lex(),
        num_vars);
}

}
