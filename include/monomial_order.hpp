/**
 * @file monomial_order.hpp
 * @brief 单项式序：Lex, GrevLex, DegLex, DegRevLex 及相关工具函数。
 */
#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

namespace lamina {

/** @brief 单项式类型，以整数向量表示各变量的指数 */
using Monomial = std::vector<int>;

/** @brief 单项式序的类型枚举 */
enum class MonomialOrderType {
    Lex,        ///< 字典序
    GrevLex,    ///< 分次逆字典序
    DegLex,     ///< 分次字典序
    DegRevLex   ///< 分次逆字典序（同 GrevLex）
};

/**
 * @brief 计算单项式的全次数
 * @param m 单项式
 * @return 各分量指数之和
 */
inline int total_degree(const Monomial& m) {
    return std::accumulate(m.begin(), m.end(), 0);
}

/**
 * @brief 计算两个单项式的最小公倍单项式
 * @param a 第一个单项式
 * @param b 第二个单项式
 * @return 各分量取最大值构成的单项式
 */
inline Monomial lcm_monomial(const Monomial& a, const Monomial& b) {
    Monomial result(std::max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = std::max(result[i], a[i]);
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = std::max(result[i], b[i]);
    return result;
}

/**
 * @brief 判断 divisor 是否整除 target
 * @param divisor 除数单项式
 * @param target 被除数单项式
 * @return divisor 的每个分量都不超过 target 对应分量时返回 true
 */
inline bool divides_monomial(const Monomial& divisor, const Monomial& target) {
    if (divisor.size() > target.size()) {
        for (size_t i = target.size(); i < divisor.size(); ++i)
            if (divisor[i] != 0) return false;
    }
    for (size_t i = 0; i < std::min(divisor.size(), target.size()); ++i) {
        if (divisor[i] > target[i]) return false;
    }
    return true;
}

/** @brief 单项式序比较器，支持 Lex、GrevLex、DegLex、DegRevLex */
class MonomialOrder {
public:
    /**
     * @brief 构造指定类型的单项式序
     * @param type 序的类型
     */
    explicit MonomialOrder(MonomialOrderType type) : type_(type) {}

    /**
     * @brief 比较两个单项式的大小
     * @param a 左操作数
     * @param b 右操作数
     * @return a > b 时返回 true（按当前序）
     */
    inline bool operator()(const Monomial& a, const Monomial& b) const {
        switch (type_) {
            case MonomialOrderType::Lex:
                return compare_lex(a, b);
            case MonomialOrderType::GrevLex:
            case MonomialOrderType::DegRevLex:
                return compare_grevlex(a, b);
            case MonomialOrderType::DegLex:
                return compare_deglex(a, b);
        }
        return false;
    }

    /** @brief 获取序的类型 */
    MonomialOrderType type() const { return type_; }

    /** @brief 创建字典序 */
    static MonomialOrder lex() { return MonomialOrder(MonomialOrderType::Lex); }

    /** @brief 创建分次逆字典序 */
    static MonomialOrder grevlex() { return MonomialOrder(MonomialOrderType::GrevLex); }

    /** @brief 创建分次字典序 */
    static MonomialOrder deglex() { return MonomialOrder(MonomialOrderType::DegLex); }

    /** @brief 创建分次逆字典序（DegRevLex） */
    static MonomialOrder degrevlex() { return MonomialOrder(MonomialOrderType::DegRevLex); }

private:
    MonomialOrderType type_;

    static inline bool compare_lex(const Monomial& a, const Monomial& b) {
        size_t n = std::max(a.size(), b.size());
        for (size_t i = 0; i < n; ++i) {
            int ai = (i < a.size()) ? a[i] : 0;
            int bi = (i < b.size()) ? b[i] : 0;
            if (ai != bi) return ai > bi;
        }
        return false;
    }

    static inline bool compare_grevlex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;

        size_t n = std::max(a.size(), b.size());
        for (size_t i = n; i > 0; --i) {
            int ai = (i - 1 < a.size()) ? a[i - 1] : 0;
            int bi = (i - 1 < b.size()) ? b[i - 1] : 0;
            if (ai != bi) return ai < bi;
        }
        return false;
    }

    static inline bool compare_deglex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;
        return compare_lex(a, b);
    }
};

}
