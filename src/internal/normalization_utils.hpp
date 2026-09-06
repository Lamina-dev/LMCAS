/** @file internal/normalization_utils.hpp */
#pragma once
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "assumption_context.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace LMCAS {

/**
 * @brief 计算 AST 节点的多项式次数
 * @param node 输入节点
 * @return 节点对应的多项式次数
 */
inline int get_node_degree_helper(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return 0;
    if (std::dynamic_pointer_cast<const VariableNode>(node)) return 1;
    if (auto p = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (auto e = std::dynamic_pointer_cast<const NumberNode>(p->exponent())) {
             if (std::holds_alternative<BigInt>(e->value())) return (int)std::get<BigInt>(e->value()).to_int();
             if (std::holds_alternative<lmmc_real_t>(e->value())) return (int)std::get<lmmc_real_t>(e->value());
        }
        return 1;
    }
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int d = 0;
        for (auto& op : m->operands()) d += get_node_degree_helper(op);
        return d;
    }
    return 0;
}

inline std::shared_ptr<const SymbolicNode> make_normalized_multiply_node(
    const std::vector<std::shared_ptr<const SymbolicNode>>& ops) {
    if (ops.empty()) return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    if (ops.size() == 1) return ops[0];
    return LMCAS::detail::make_node<MultiplyNode>(ops);
}

/**
 * @brief 将节点树中名为 index_var 的变量替换为指定数值（用于展开有限求和/连乘）。
 *
 * 仅处理求和展开所需的节点类型；嵌套的 Summation/Product 若使用同名指标则不深入其绑定体。
 */
inline std::shared_ptr<const SymbolicNode> norm_subst_index(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& index_var,
    const std::shared_ptr<const SymbolicNode>& value) {
    return LMCAS::substitute_free(node, index_var, value);
}

/** @brief 节点比较器，按多项式次数降序排列，用于同类项合并 */
struct NodeCompare {    bool operator()(const std::shared_ptr<const SymbolicNode>& lhs, const std::shared_ptr<const SymbolicNode>& rhs) const {
        if (!lhs && !rhs) return false;
        if (!lhs) return true;
        if (!rhs) return false;

        int d1 = get_node_degree_helper(lhs);
        int d2 = get_node_degree_helper(rhs);
        if (d1 != d2) return d1 > d2;

        bool isNum1 = std::dynamic_pointer_cast<const NumberNode>(lhs) != nullptr;
        bool isNum2 = std::dynamic_pointer_cast<const NumberNode>(rhs) != nullptr;
        if (isNum1 != isNum2) return isNum1;

        if (lhs->type_priority() != rhs->type_priority()) {
            return lhs->type_priority() < rhs->type_priority();
        }

        return lhs->compare(*rhs) < 0;
    }
};

/**
 * @brief 两个数值节点相加
 * @param a 加数节点
 * @param b 加数节点
 * @return 和的数值节点
 */
inline std::shared_ptr<const NumberNode> add_numbers(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
     if (std::holds_alternative<lmmc_real_t>(a->value()) || std::holds_alternative<lmmc_real_t>(b->value())) {
         lmmc_real_t v1 = std::holds_alternative<lmmc_real_t>(a->value()) ? std::get<lmmc_real_t>(a->value()) :
                     (std::holds_alternative<Rational>(a->value()) ? (lmmc_real_t)std::get<Rational>(a->value()).to_double() : (lmmc_real_t)std::get<BigInt>(a->value()).to_double());
         lmmc_real_t v2 = std::holds_alternative<lmmc_real_t>(b->value()) ? std::get<lmmc_real_t>(b->value()) :
                     (std::holds_alternative<Rational>(b->value()) ? (lmmc_real_t)std::get<Rational>(b->value()).to_double() : (lmmc_real_t)std::get<BigInt>(b->value()).to_double());
         lmmc_real_t sum = v1 + v2;
         return LMCAS::detail::make_node<NumberNode>(sum);
     }

     if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
         Rational r1 = std::holds_alternative<Rational>(a->value()) ? std::get<Rational>(a->value()) :
                       (std::holds_alternative<BigInt>(a->value()) ? Rational(std::get<BigInt>(a->value())) : Rational(0));
         Rational r2 = std::holds_alternative<Rational>(b->value()) ? std::get<Rational>(b->value()) :
                       (std::holds_alternative<BigInt>(b->value()) ? Rational(std::get<BigInt>(b->value())) : Rational(0));
         return LMCAS::detail::make_node<NumberNode>(r1 + r2);
     }

     BigInt i1 = std::get<BigInt>(a->value());
     BigInt i2 = std::get<BigInt>(b->value());
     return LMCAS::detail::make_node<NumberNode>(i1 + i2);
}

/**
 * @brief 两个数值节点相乘
 * @param a 乘数节点
 * @param b 乘数节点
 * @return 积的数值节点
 */
inline std::shared_ptr<const NumberNode> multiply_numbers(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
     if (std::holds_alternative<lmmc_real_t>(a->value()) || std::holds_alternative<lmmc_real_t>(b->value())) {
         lmmc_real_t v1 = std::holds_alternative<lmmc_real_t>(a->value()) ? std::get<lmmc_real_t>(a->value()) :
                     (std::holds_alternative<Rational>(a->value()) ? (lmmc_real_t)std::get<Rational>(a->value()).to_double() : (lmmc_real_t)std::get<BigInt>(a->value()).to_double());
         lmmc_real_t v2 = std::holds_alternative<lmmc_real_t>(b->value()) ? std::get<lmmc_real_t>(b->value()) :
                     (std::holds_alternative<Rational>(b->value()) ? (lmmc_real_t)std::get<Rational>(b->value()).to_double() : (lmmc_real_t)std::get<BigInt>(b->value()).to_double());
         lmmc_real_t prod = v1 * v2;
         return LMCAS::detail::make_node<NumberNode>(prod);
     }

     if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
         Rational r1 = std::holds_alternative<Rational>(a->value()) ? std::get<Rational>(a->value()) :
                       (std::holds_alternative<BigInt>(a->value()) ? Rational(std::get<BigInt>(a->value())) : Rational(1));
         Rational r2 = std::holds_alternative<Rational>(b->value()) ? std::get<Rational>(b->value()) :
                       (std::holds_alternative<BigInt>(b->value()) ? Rational(std::get<BigInt>(b->value())) : Rational(1));
         return LMCAS::detail::make_node<NumberNode>(r1 * r2);
     }

     BigInt i1 = std::get<BigInt>(a->value());
     BigInt i2 = std::get<BigInt>(b->value());
     return LMCAS::detail::make_node<NumberNode>(i1 * i2);
}

/**
 * @brief 检查节点是否为负数或带负系数，若是则输出其正值部分
 * @param arg 待检查的节点
 * @param out_positive 输出参数，存放取正后的节点
 * @return 若节点为负则返回 true
 */
inline bool check_negative_arg(const std::shared_ptr<const SymbolicNode>& arg, std::shared_ptr<const SymbolicNode>& out_positive) {
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(arg)) {
        if (std::holds_alternative<lmmc_real_t>(num->value()) && std::get<lmmc_real_t>(num->value()) < 0) {
             out_positive = LMCAS::detail::make_node<NumberNode>(std::abs(std::get<lmmc_real_t>(num->value())));
             return true;
        }
        if (std::holds_alternative<BigInt>(num->value()) &&
            std::get<BigInt>(num->value()).IsNegative()) {
             out_positive = LMCAS::detail::make_node<NumberNode>(std::get<BigInt>(num->value()) * BigInt(-1));
             return true;
        }
        if (std::holds_alternative<Rational>(num->value()) &&
            std::get<Rational>(num->value()).get_numerator().IsNegative()) {
             out_positive = LMCAS::detail::make_node<NumberNode>(std::get<Rational>(num->value()) * Rational(-1));
             return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(arg)) {
        if (!mul->operands().empty()) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0])) {
                 bool is_neg = false;
                 std::shared_ptr<const NumberNode> pos_num = nullptr;
                 if (std::holds_alternative<lmmc_real_t>(num->value()) && std::get<lmmc_real_t>(num->value()) < 0) {
                     is_neg = true;
                     pos_num = LMCAS::detail::make_node<NumberNode>(std::abs(std::get<lmmc_real_t>(num->value())));
                 } else if (std::holds_alternative<BigInt>(num->value()) &&
                            std::get<BigInt>(num->value()).IsNegative()) {
                     is_neg = true;
                     pos_num = LMCAS::detail::make_node<NumberNode>(std::get<BigInt>(num->value()) * BigInt(-1));
                 } else if (std::holds_alternative<Rational>(num->value()) &&
                            std::get<Rational>(num->value()).get_numerator().IsNegative()) {
                     is_neg = true;
                     pos_num = LMCAS::detail::make_node<NumberNode>(std::get<Rational>(num->value()) * Rational(-1));
                 }

                 if (is_neg) {
                     std::vector<std::shared_ptr<const SymbolicNode>> new_ops = mul->operands();

                     bool is_minus_one = false;
                     if (std::holds_alternative<lmmc_real_t>(num->value())) {
                         int eq_1;
                         lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(num->value()), -1.0, 1e-9, 1e-9, &eq_1);
                         is_minus_one = (eq_1 != 0);
                     }
                     else if (std::holds_alternative<BigInt>(num->value())) is_minus_one = (std::get<BigInt>(num->value()) == BigInt(-1));
                     else if (std::holds_alternative<Rational>(num->value())) is_minus_one = (std::get<Rational>(num->value()) == Rational(-1));

                     if (is_minus_one) {
                         new_ops.erase(new_ops.begin());
                         if (new_ops.empty()) out_positive = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                         else if (new_ops.size() == 1) out_positive = new_ops[0];
                         else out_positive = make_normalized_multiply_node(new_ops);
                     } else {

                         new_ops[0] = pos_num;
                         out_positive = make_normalized_multiply_node(new_ops);
                     }
                     return true;
                 }
            }
        }
    }
    return false;
}

/**
 * @brief 提取节点中 π 的有理系数（如 2π/3 中的 2/3）
 * @param node 待检查的节点
 * @param k 输出参数，存放 π 的有理系数
 * @return 若节点为 k*π 形式则返回 true
 */
inline bool get_pi_coeff(const std::shared_ptr<const SymbolicNode>& node, Rational& k) {
    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == "pi") {
            k = Rational(1);
            return true;
        }
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        bool has_pi = false;
        k = Rational(1);

        for (const auto& op : mul->operands()) {
            if (auto v = std::dynamic_pointer_cast<const VariableNode>(op)) {
                if (v->name() == "pi") {
                    if (has_pi) return false;
                    has_pi = true;
                } else return false;
            } else if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (std::holds_alternative<Rational>(n->value())) k = k * std::get<Rational>(n->value());
                else if (std::holds_alternative<BigInt>(n->value())) k = k * Rational(std::get<BigInt>(n->value()));
                else return false;
            } else return false;
        }
        return has_pi;
    }
    return false;
}

} // namespace LMCAS
