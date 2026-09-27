/**
 * @file fglm.hpp
 * @brief FGLM 算法:零维理想 Gröbner 基在不同单项式序之间的转换.
 *
 * @see J.-C. Faugère, P. Gianni, D. Lazard, and T. Mora,
 *      "Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering,"
 *      Journal of Symbolic Computation 16(4), 1993.
 */
#pragma once

#include "lmcas_export.hpp"
#include "monomial_order.hpp"
#include "rational.hpp"
#include <vector>

namespace LMCAS {

/** @brief FGLM 算法中使用的多项式表示,以 (单项式, 系数) 对列表存储 */
struct LMCAS_API FGLMPoly {
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
    Monomial lead_monomial() const {
        if (terms.empty()) return Monomial();
        return terms.front().first;
    }

    /** @brief 获取首项系数 */
    Rational lead_coeff() const {
        if (terms.empty()) return Rational(0);
        return terms.front().second;
    }

    /**
     * @brief 添加一项
     * @param m 单项式
     * @param c 系数(零系数不添加)
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
    void sort_terms(const MonomialOrder& order);

    /** @brief 合并同类项并去除零项 */
    void normalize();

    /**
     * @brief 对给定基进行多项式约化
     * @param basis Gröbner 基
     * @param order 单项式序
     * @return 约化后的余式
     */
    FGLMPoly reduce(const std::vector<FGLMPoly>& basis,
                    const MonomialOrder& order) const;

    /**
     * @brief 从单项式构造多项式(系数为 1)
     * @param m 单项式
     * @param num_vars 变量个数
     * @return 仅含一项的多项式
     */
    static FGLMPoly from_monomial(const Monomial& m, size_t num_vars);

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
 * @brief 计算多项式关于基的标准形(余式)
 * @param f 待约化多项式
 * @param basis Gröbner 基
 * @param order 单项式序
 * @return 标准形
 */
LMCAS_API FGLMPoly normal_form(const FGLMPoly& f,
                              const std::vector<FGLMPoly>& basis,
                              const MonomialOrder& order);

/**
 * @brief 判断理想是否为零维(每个变量都有纯幂次首项)
 * @param basis Gröbner 基
 * @param num_vars 变量个数
 * @return 零维返回 true
 */
LMCAS_API bool is_zero_dimensional(const std::vector<FGLMPoly>& basis, size_t num_vars);

/**
 * @brief 计算商空间维数(标准单项式个数)
 * @param basis Gröbner 基
 * @param num_vars 变量个数
 * @param max_degree 搜索的最大次数上限
 * @return 商空间维数;-1 表示当前搜索界内维数保持未知
 */
LMCAS_API int quotient_dimension(const std::vector<FGLMPoly>& basis,
                                  size_t num_vars, int max_degree = 50);


/**
 * @brief FGLM 算法:将零维理想的 Gröbner 基从源序转换到目标序
 * @param source_basis 源序下的 Gröbner 基
 * @param source_order 源单项式序
 * @param target_order 目标单项式序
 * @param num_vars 变量个数
 * @return 目标序下的 Gröbner 基
 * @throw std::runtime_error 理想非零维或维数过大时抛出
 * @see J.-C. Faugère, P. Gianni, D. Lazard, and T. Mora,
 *      "Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering,"
 *      Journal of Symbolic Computation 16(4), 1993.
 */
LMCAS_API std::vector<FGLMPoly> fglm_convert(
    const std::vector<FGLMPoly>& source_basis,
    const MonomialOrder& source_order,
    const MonomialOrder& target_order, size_t num_vars);

/**
 * @brief 将 GrevLex 序下的 Gröbner 基转换为 Lex 序
 * @param grevlex_basis GrevLex 序下的 Gröbner 基
 * @param num_vars 变量个数
 * @return Lex 序下的 Gröbner 基
 */
LMCAS_API std::vector<FGLMPoly> grevlex_to_lex(
    const std::vector<FGLMPoly>& grevlex_basis, size_t num_vars);

}
