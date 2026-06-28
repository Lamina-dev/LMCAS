/**
 * @file normalization_visitor.hpp
 * @brief 规范化访问器，合并同类项、化简数值运算、处理矩阵算术。
 */
#pragma once

#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "../symbolic_ast.hpp"
#include "../assumption_context.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

/**
 * @brief 计算 AST 节点的多项式次数
 * @param node 输入节点
 * @return 节点对应的多项式次数
 */
inline int get_node_degree_helper(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return 0;
    if (std::dynamic_pointer_cast<VariableNode>(node)) return 1;
    if (auto p = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (auto e = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
             if (std::holds_alternative<BigInt>(e->value)) return (int)std::get<BigInt>(e->value).to_int();
             if (std::holds_alternative<lmmc_real_t>(e->value)) return (int)std::get<lmmc_real_t>(e->value);
        }
        return 1;
    }
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        int d = 0;
        for (auto& op : m->operands) d += get_node_degree_helper(op);
        return d;
    }
    return 0;
}

/**
 * @brief 将节点树中名为 index_var 的变量替换为指定数值（用于展开有限求和/连乘）。
 *
 * 仅处理求和展开所需的节点类型；嵌套的 Summation/Product 若使用同名指标则不深入其绑定体。
 */
inline std::shared_ptr<SymbolicNode> norm_subst_index(
    const std::shared_ptr<SymbolicNode>& node,
    const std::string& index_var,
    const std::shared_ptr<SymbolicNode>& value) {
    if (!node) return node;
    if (auto v = std::dynamic_pointer_cast<VariableNode>(node)) {
        if (v->name == index_var) return value->clone();
        return node->clone();
    }
    if (std::dynamic_pointer_cast<NumberNode>(node)) return node->clone();
    if (auto a = std::dynamic_pointer_cast<AddNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        for (auto& op : a->operands) ops.push_back(norm_subst_index(op, index_var, value));
        return std::make_shared<AddNode>(ops);
    }
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        for (auto& op : m->operands) ops.push_back(norm_subst_index(op, index_var, value));
        return std::make_shared<MultiplyNode>(ops);
    }
    if (auto p = std::dynamic_pointer_cast<PowerNode>(node)) {
        return std::make_shared<PowerNode>(
            norm_subst_index(p->base, index_var, value),
            norm_subst_index(p->exponent, index_var, value));
    }
    if (auto f = std::dynamic_pointer_cast<FunctionNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> args;
        for (auto& arg : f->arguments) args.push_back(norm_subst_index(arg, index_var, value));
        return std::make_shared<FunctionNode>(f->type, args);
    }
    // 其它节点类型：保守地原样克隆
    return node->clone();
}

/** @brief 节点比较器，按多项式次数降序排列，用于同类项合并 */
struct NodeCompare {    bool operator()(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) const {
        if (!lhs && !rhs) return false;
        if (!lhs) return true;
        if (!rhs) return false;

        int d1 = get_node_degree_helper(lhs);
        int d2 = get_node_degree_helper(rhs);
        if (d1 != d2) return d1 > d2;

        bool isNum1 = std::dynamic_pointer_cast<NumberNode>(lhs) != nullptr;
        bool isNum2 = std::dynamic_pointer_cast<NumberNode>(rhs) != nullptr;
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
inline std::shared_ptr<NumberNode> add_numbers(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
     if (std::holds_alternative<lmmc_real_t>(a->value) || std::holds_alternative<lmmc_real_t>(b->value)) {
         lmmc_real_t v1 = std::holds_alternative<lmmc_real_t>(a->value) ? std::get<lmmc_real_t>(a->value) :
                     (std::holds_alternative<Rational>(a->value) ? (lmmc_real_t)std::get<Rational>(a->value).to_double() : (lmmc_real_t)std::get<BigInt>(a->value).to_double());
         lmmc_real_t v2 = std::holds_alternative<lmmc_real_t>(b->value) ? std::get<lmmc_real_t>(b->value) :
                     (std::holds_alternative<Rational>(b->value) ? (lmmc_real_t)std::get<Rational>(b->value).to_double() : (lmmc_real_t)std::get<BigInt>(b->value).to_double());
         lmmc_real_t sum = v1 + v2;
         return std::make_shared<NumberNode>(sum);
     }

     if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
         Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) :
                       (std::holds_alternative<BigInt>(a->value) ? Rational(std::get<BigInt>(a->value)) : Rational(0));
         Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) :
                       (std::holds_alternative<BigInt>(b->value) ? Rational(std::get<BigInt>(b->value)) : Rational(0));
         return std::make_shared<NumberNode>(r1 + r2);
     }

     BigInt i1 = std::get<BigInt>(a->value);
     BigInt i2 = std::get<BigInt>(b->value);
     return std::make_shared<NumberNode>(i1 + i2);
}

/**
 * @brief 两个数值节点相乘
 * @param a 乘数节点
 * @param b 乘数节点
 * @return 积的数值节点
 */
inline std::shared_ptr<NumberNode> multiply_numbers(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
     if (std::holds_alternative<lmmc_real_t>(a->value) || std::holds_alternative<lmmc_real_t>(b->value)) {
         lmmc_real_t v1 = std::holds_alternative<lmmc_real_t>(a->value) ? std::get<lmmc_real_t>(a->value) :
                     (std::holds_alternative<Rational>(a->value) ? (lmmc_real_t)std::get<Rational>(a->value).to_double() : (lmmc_real_t)std::get<BigInt>(a->value).to_double());
         lmmc_real_t v2 = std::holds_alternative<lmmc_real_t>(b->value) ? std::get<lmmc_real_t>(b->value) :
                     (std::holds_alternative<Rational>(b->value) ? (lmmc_real_t)std::get<Rational>(b->value).to_double() : (lmmc_real_t)std::get<BigInt>(b->value).to_double());
         lmmc_real_t prod = v1 * v2;
         return std::make_shared<NumberNode>(prod);
     }

     if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
         Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) :
                       (std::holds_alternative<BigInt>(a->value) ? Rational(std::get<BigInt>(a->value)) : Rational(1));
         Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) :
                       (std::holds_alternative<BigInt>(b->value) ? Rational(std::get<BigInt>(b->value)) : Rational(1));
         return std::make_shared<NumberNode>(r1 * r2);
     }

     BigInt i1 = std::get<BigInt>(a->value);
     BigInt i2 = std::get<BigInt>(b->value);
     return std::make_shared<NumberNode>(i1 * i2);
}

/**
 * @brief 检查节点是否为负数或带负系数，若是则输出其正值部分
 * @param arg 待检查的节点
 * @param out_positive 输出参数，存放取正后的节点
 * @return 若节点为负则返回 true
 */
inline bool check_negative_arg(const std::shared_ptr<SymbolicNode>& arg, std::shared_ptr<SymbolicNode>& out_positive) {
    if (auto num = std::dynamic_pointer_cast<NumberNode>(arg)) {
        if (std::holds_alternative<lmmc_real_t>(num->value) && std::get<lmmc_real_t>(num->value) < 0) {
             out_positive = std::make_shared<NumberNode>(std::abs(std::get<lmmc_real_t>(num->value)));
             return true;
        }
        if (std::holds_alternative<BigInt>(num->value) && std::get<BigInt>(num->value).to_double() < 0) {
             out_positive = std::make_shared<NumberNode>(std::get<BigInt>(num->value) * BigInt(-1));
             return true;
        }
        if (std::holds_alternative<Rational>(num->value) && std::get<Rational>(num->value).to_double() < 0) {
             out_positive = std::make_shared<NumberNode>(std::get<Rational>(num->value) * Rational(-1));
             return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(arg)) {
        if (!mul->operands.empty()) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(mul->operands[0])) {
                 bool is_neg = false;
                 std::shared_ptr<NumberNode> pos_num = nullptr;
                 if (std::holds_alternative<lmmc_real_t>(num->value) && std::get<lmmc_real_t>(num->value) < 0) {
                     is_neg = true;
                     pos_num = std::make_shared<NumberNode>(std::abs(std::get<lmmc_real_t>(num->value)));
                 } else if (std::holds_alternative<BigInt>(num->value) && std::get<BigInt>(num->value).to_double() < 0) {
                     is_neg = true;
                     pos_num = std::make_shared<NumberNode>(std::get<BigInt>(num->value) * BigInt(-1));
                 } else if (std::holds_alternative<Rational>(num->value) && std::get<Rational>(num->value).to_double() < 0) {
                     is_neg = true;
                     pos_num = std::make_shared<NumberNode>(std::get<Rational>(num->value) * Rational(-1));
                 }

                 if (is_neg) {
                     std::vector<std::shared_ptr<SymbolicNode>> new_ops = mul->operands;

                     bool is_minus_one = false;
                     if (std::holds_alternative<lmmc_real_t>(num->value)) {
                         int eq_1;
                         lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(num->value), -1.0, 1e-9, 1e-9, &eq_1);
                         is_minus_one = (eq_1 != 0);
                     }
                     else if (std::holds_alternative<BigInt>(num->value)) is_minus_one = (std::get<BigInt>(num->value) == BigInt(-1));
                     else if (std::holds_alternative<Rational>(num->value)) is_minus_one = (std::get<Rational>(num->value) == Rational(-1));

                     if (is_minus_one) {
                         new_ops.erase(new_ops.begin());
                         if (new_ops.size() == 1) out_positive = new_ops[0];
                         else out_positive = std::make_shared<MultiplyNode>(new_ops);
                     } else {

                         new_ops[0] = pos_num;
                         out_positive = std::make_shared<MultiplyNode>(new_ops);
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
inline bool get_pi_coeff(const std::shared_ptr<SymbolicNode>& node, Rational& k) {
    if (auto v = std::dynamic_pointer_cast<VariableNode>(node)) {
        if (v->name == "pi") {
            k = Rational(1);
            return true;
        }
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        bool has_pi = false;
        k = Rational(1);

        for (const auto& op : mul->operands) {
            if (auto v = std::dynamic_pointer_cast<VariableNode>(op)) {
                if (v->name == "pi") {
                    if (has_pi) return false;
                    has_pi = true;
                } else return false;
            } else if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) {
                if (std::holds_alternative<Rational>(n->value)) k = k * std::get<Rational>(n->value);
                else if (std::holds_alternative<BigInt>(n->value)) k = k * Rational(std::get<BigInt>(n->value));
                else return false;
            } else return false;
        }
        return has_pi;
    }
    return false;
}

/** @brief 规范化访问器，合并同类项、化简数值运算、处理矩阵算术和三角函数特殊值 */
class NormalizationVisitor : public SymbolicVisitor {
public:
    std::shared_ptr<SymbolicNode> result;  ///< 规范化结果节点

    /// Default constructor (backward compatible, no assumptions).
    NormalizationVisitor() : assumptions_(nullptr) {}

    /// Construct with an optional AssumptionContext for assumption-aware simplification.
    explicit NormalizationVisitor(const lamina::AssumptionContext* ctx)
        : assumptions_(ctx) {}

    /**
     * @brief 获取规范化结果
     * @return 化简后的 AST 节点
     */
    std::shared_ptr<SymbolicNode> get_result() const {
        return result;
    }

    /**
     * @brief 展开两个节点的乘积（分配律）
     * @param lhs 左操作数
     * @param rhs 右操作数
     * @return 展开后的 AST 节点
     */
    std::shared_ptr<SymbolicNode> expand_product(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) {

        auto add_lhs = std::dynamic_pointer_cast<AddNode>(lhs);
        auto add_rhs = std::dynamic_pointer_cast<AddNode>(rhs);
        
        auto c_lhs = std::dynamic_pointer_cast<ComplexNode>(lhs);
        auto c_rhs = std::dynamic_pointer_cast<ComplexNode>(rhs);
        
        if (c_lhs && c_rhs) {
            auto ac = expand_product(c_lhs->real, c_rhs->real);
            auto bd = expand_product(c_lhs->imag, c_rhs->imag);
            auto ad = expand_product(c_lhs->real, c_rhs->imag);
            auto bc = expand_product(c_lhs->imag, c_rhs->real);
            auto neg_one = std::make_shared<NumberNode>(BigInt(-1));
            auto neg_bd = expand_product(neg_one, bd);
            auto real_part = SymbolicFactory::create_add({ac, neg_bd});
            auto imag_part = SymbolicFactory::create_add({ad, bc});
            NormalizationVisitor norm(assumptions_);
            real_part->accept(norm); auto norm_r = norm.get_result();
            imag_part->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        } else if (c_lhs) {
            auto nr = expand_product(c_lhs->real, rhs);
            auto ni = expand_product(c_lhs->imag, rhs);
            NormalizationVisitor norm(assumptions_);
            nr->accept(norm); auto norm_r = norm.get_result();
            ni->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        } else if (c_rhs) {
            auto nr = expand_product(lhs, c_rhs->real);
            auto ni = expand_product(lhs, c_rhs->imag);
            NormalizationVisitor norm(assumptions_);
            nr->accept(norm); auto norm_r = norm.get_result();
            ni->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        }

        if (add_lhs && add_rhs) {

            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op1 : add_lhs->operands) {
                for (const auto& op2 : add_rhs->operands) {
                    auto prod = expand_product(op1, op2);
                    new_terms.push_back(prod);
                }
            }

            return std::make_shared<AddNode>(new_terms);
        } else if (add_lhs) {

            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op : add_lhs->operands) {
                 auto prod = expand_product(op, rhs);
                 new_terms.push_back(prod);
            }
            return std::make_shared<AddNode>(new_terms);
        } else if (add_rhs) {

            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op : add_rhs->operands) {
                 auto prod = expand_product(lhs, op);
                 new_terms.push_back(prod);
            }
            return std::make_shared<AddNode>(new_terms);
        }

        std::shared_ptr<NumberNode> const_acc = std::make_shared<NumberNode>(BigInt(1));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> bases;

        auto process_factor = [&](const std::shared_ptr<SymbolicNode>& factor) {
             if (auto num = std::dynamic_pointer_cast<NumberNode>(factor)) {
                 if (num->is_zero()) return false;
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<SymbolicNode> base = factor;
                 std::shared_ptr<NumberNode> exp = std::make_shared<NumberNode>(BigInt(1));

                 if (auto pow = std::dynamic_pointer_cast<PowerNode>(factor)) {
                     base = pow->base;
                     if (auto e_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                         exp = e_num;
                     }
                 }

                 auto it = bases.find(base);
                 if (it == bases.end()) {
                     bases[base] = exp;
                 } else {
                     bases[base] = add_numbers(it->second, exp);
                 }
             }
             return true;
        };

        auto flatten_and_process = [&](const std::shared_ptr<SymbolicNode>& node) {
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
                for(const auto& op : mul->operands) {
                    if (!process_factor(op)) return false;
                }
            } else {
                if (!process_factor(node)) return false;
            }
            return true;
        };

        if (!flatten_and_process(lhs)) return std::make_shared<NumberNode>(BigInt(0));
        if (!flatten_and_process(rhs)) return std::make_shared<NumberNode>(BigInt(0));

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }

        std::vector<std::shared_ptr<SymbolicNode>> var_ops;
        for (auto const& [base, exp] : bases) {
            if (exp->is_zero()) {
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(std::make_shared<PowerNode>(base, exp));
            }
        }

        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());

        if (final_ops.empty()) return std::make_shared<NumberNode>(BigInt(1));
        if (final_ops.size() == 1) return final_ops[0];

        if (final_ops.size() > 1 && std::dynamic_pointer_cast<NumberNode>(final_ops.back())) {
             std::rotate(final_ops.begin(), final_ops.end() - 1, final_ops.end());
        }
        return std::make_shared<MultiplyNode>(final_ops);
    }

    void visit(NumberNode& node) override {
        result = node.clone();
    }

    void visit(VariableNode& node) override {
        result = node.clone();
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> simplified_ops;

        for (const auto& op : node.operands) {
            op->accept(*this);
            if (auto add = std::dynamic_pointer_cast<AddNode>(result)) {
                simplified_ops.insert(simplified_ops.end(), add->operands.begin(), add->operands.end());
            } else {
                simplified_ops.push_back(result);
            }
        }

        // Merge ComplexNodes and NumberNodes
        std::vector<std::shared_ptr<SymbolicNode>> real_parts;
        std::vector<std::shared_ptr<SymbolicNode>> imag_parts;
        std::vector<std::shared_ptr<SymbolicNode>> non_complex_ops;
        bool has_complex = false;
        
        for (const auto& op : simplified_ops) {
            if (auto c = std::dynamic_pointer_cast<ComplexNode>(op)) {
                has_complex = true;
                real_parts.push_back(c->real);
                imag_parts.push_back(c->imag);
            } else if (std::dynamic_pointer_cast<NumberNode>(op)) {
                real_parts.push_back(op);
                imag_parts.push_back(SymbolicFactory::create_number(BigInt(0)));
            } else {
                non_complex_ops.push_back(op);
            }
        }
        
        if (has_complex) {
            auto real_sum = SymbolicFactory::create_add(real_parts);
            auto imag_sum = SymbolicFactory::create_add(imag_parts);
            NormalizationVisitor sub_norm(assumptions_);
            real_sum->accept(sub_norm); auto norm_real = sub_norm.get_result();
            imag_sum->accept(sub_norm); auto norm_imag = sub_norm.get_result();
            auto merged_complex = SymbolicFactory::create_complex(norm_real, norm_imag);
            if (!merged_complex->is_zero()) {
                non_complex_ops.push_back(merged_complex);
            }
            simplified_ops = non_complex_ops;
        }

        if (!simplified_ops.empty() && std::dynamic_pointer_cast<MatrixNode>(simplified_ops[0])) {
            auto first_mat = std::dynamic_pointer_cast<MatrixNode>(simplified_ops[0]);
            size_t rows = first_mat->rows;
            size_t cols = first_mat->cols;
            bool all_matrices = true;
            for (const auto& op : simplified_ops) {
                 auto m = std::dynamic_pointer_cast<MatrixNode>(op);
                 if (!m || m->rows != rows || m->cols != cols) {
                     all_matrices = false; break;
                 }
            }

            if (all_matrices) {
                 std::vector<std::shared_ptr<SymbolicNode>> new_elements;
                 new_elements.reserve(rows * cols);

                 for (size_t i = 0; i < rows * cols; ++i) {
                     std::vector<std::shared_ptr<SymbolicNode>> elem_ops;
                     for (const auto& op : simplified_ops) {
                         auto m = std::dynamic_pointer_cast<MatrixNode>(op);
                         std::shared_ptr<SymbolicNode> val;
                         if (std::holds_alternative<MatrixNode::DenseStorage>(m->storage)) {
                             const auto& dense = std::get<MatrixNode::DenseStorage>(m->storage);
                             if (i < dense.size()) val = dense[i];
                             else val = std::make_shared<NumberNode>(BigInt(0));
                         } else {
                             const auto& sparse = std::get<MatrixNode::SparseStorage>(m->storage);
                             auto it = sparse.find(i);
                             if (it != sparse.end()) val = it->second;
                             else val = std::make_shared<NumberNode>(BigInt(0));
                         }
                         if (!val) val = std::make_shared<NumberNode>(BigInt(0));
                         elem_ops.push_back(val);
                     }

                     auto elem_add = std::make_shared<AddNode>(elem_ops);
                     elem_add->accept(*this);
                     new_elements.push_back(result);
                 }

                 result = std::make_shared<MatrixNode>(rows, cols, new_elements);
                 return;
            }
        }

        std::shared_ptr<NumberNode> constant_acc = std::make_shared<NumberNode>(BigInt(0));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> terms;

        for (const auto& op : simplified_ops) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                constant_acc = add_numbers(constant_acc, num);
            } else {
                std::shared_ptr<SymbolicNode> term_part = op;
                std::shared_ptr<NumberNode> coeff_part = std::make_shared<NumberNode>(BigInt(1));

                if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {

                    std::shared_ptr<NumberNode> coeff = nullptr;
                    int coeff_idx = -1;

                    if (!mul->operands.empty()) {
                        if (auto n_first = std::dynamic_pointer_cast<NumberNode>(mul->operands.front())) {
                            coeff = n_first;
                            coeff_idx = 0;
                        } else if (auto n_last = std::dynamic_pointer_cast<NumberNode>(mul->operands.back())) {
                            coeff = n_last;
                            coeff_idx = (int)mul->operands.size() - 1;
                        }
                    }

                    if (coeff) {
                         coeff_part = coeff;
                         if (mul->operands.size() == 2) {
                             term_part = mul->operands[coeff_idx == 0 ? 1 : 0];
                         } else {
                             std::vector<std::shared_ptr<SymbolicNode>> rest;
                             rest.reserve(mul->operands.size() - 1);
                             for(int k=0; k<(int)mul->operands.size(); ++k) {
                                 if (k != coeff_idx) rest.push_back(mul->operands[k]);
                             }

                             term_part = std::make_shared<MultiplyNode>(rest);
                         }
                    }
                }

                auto it = terms.find(term_part);
                if (it == terms.end()) {
                    terms[term_part] = coeff_part;
                } else {
                    terms[term_part] = add_numbers(it->second, coeff_part);
                }
            }
        }

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;

        if (!constant_acc->is_zero()) {
            final_ops.push_back(constant_acc);
        }

        for (auto const& [term, coeff] : terms) {
            if (coeff->is_zero()) continue;

            if (coeff->is_one()) {
                final_ops.push_back(term);
            } else {
                std::vector<std::shared_ptr<SymbolicNode>> m_ops;
                m_ops.push_back(coeff);
                if (auto m = std::dynamic_pointer_cast<MultiplyNode>(term)) {
                     m_ops.insert(m_ops.end(), m->operands.begin(), m->operands.end());
                } else {
                     m_ops.push_back(term);
                }
                final_ops.push_back(std::make_shared<MultiplyNode>(m_ops));
            }
        }

        if (final_ops.empty()) {
            result = std::make_shared<NumberNode>(BigInt(0));
        } else if (final_ops.size() == 1) {
            result = final_ops[0];
        } else {

            std::sort(final_ops.begin(), final_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
                int d1 = get_node_degree_helper(l);
                int d2 = get_node_degree_helper(r);
                if (d1 != d2) return d1 > d2;

                bool isNum1 = std::dynamic_pointer_cast<NumberNode>(l) != nullptr;
                bool isNum2 = std::dynamic_pointer_cast<NumberNode>(r) != nullptr;
                if (isNum1 != isNum2) return isNum2;

                return l->compare(*r) < 0;
            });
            result = std::make_shared<AddNode>(final_ops);
        }
    }

    void visit(ComplexNode& node) override {
        node.real->accept(*this);
        auto norm_r = result;
        node.imag->accept(*this);
        auto norm_i = result;
        result = SymbolicFactory::create_complex(norm_r, norm_i);
    }

    void visit(MultiplyNode& node) override {

        std::vector<std::shared_ptr<SymbolicNode>> sc;

        for (const auto& op : node.operands) {
            op->accept(*this);
            auto res = result;
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(res)) {

                sc.insert(sc.end(), mul->operands.begin(), mul->operands.end());
            } else {
                sc.push_back(res);
            }
        }

        bool has_add = false;

        for(const auto& op : sc) {
            if(std::dynamic_pointer_cast<AddNode>(op)) {
                has_add = true;
                break;
            }
        }

        if (has_add) {
            if (sc.empty()) {
                result = std::make_shared<NumberNode>(BigInt(1));
                return;
            }
            std::shared_ptr<SymbolicNode> current = sc[0];
            for(size_t i=1; i<sc.size(); ++i) {

                current = expand_product(current, sc[i]);
            }

            result = current;
            return;
        }

        bool has_matrix = false;
        for(const auto& op : sc) {
            if (std::dynamic_pointer_cast<MatrixNode>(op)) { has_matrix = true; break; }
            if (auto p = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (std::dynamic_pointer_cast<MatrixNode>(p->base)) { has_matrix = true; break; }
            }
        }

        if (has_matrix) {
            std::vector<std::shared_ptr<SymbolicNode>> new_ops;
            std::shared_ptr<NumberNode> scalar_part = std::make_shared<NumberNode>(BigInt(1));

            for(const auto& op : sc) {
                if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                    scalar_part = multiply_numbers(scalar_part, num);
                } else {
                    new_ops.push_back(op);
                }
            }

            std::vector<std::shared_ptr<SymbolicNode>> fused_ops;
            if (!new_ops.empty()) fused_ops.push_back(new_ops[0]);

            for(size_t i=1; i<new_ops.size(); ++i) {
                auto left = fused_ops.back();
                auto right = new_ops[i];

                auto m_left = std::dynamic_pointer_cast<MatrixNode>(left);
                auto m_right = std::dynamic_pointer_cast<MatrixNode>(right);

                if (m_left && m_right) {

                    if (m_left->cols == m_right->rows) {

                         if (std::holds_alternative<MatrixNode::DenseStorage>(m_left->storage) &&
                             std::holds_alternative<MatrixNode::DenseStorage>(m_right->storage)) {

                             const auto& d_l = std::get<MatrixNode::DenseStorage>(m_left->storage);
                             const auto& d_r = std::get<MatrixNode::DenseStorage>(m_right->storage);

                             size_t R = m_left->rows;
                             size_t C = m_right->cols;
                             size_t K = m_left->cols;

                             MatrixNode::DenseStorage res_data;
                             res_data.reserve(R*C);

                             for(size_t r=0; r<R; ++r) {
                                 for(size_t c=0; c<C; ++c) {

                                     std::vector<std::shared_ptr<SymbolicNode>> sum_ops;
                                     for(size_t k=0; k<K; ++k) {
                                         std::vector<std::shared_ptr<SymbolicNode>> prod_ops = {
                                             d_l[r*K + k], d_r[k*C + c]
                                         };
                                         sum_ops.push_back(std::make_shared<MultiplyNode>(prod_ops));
                                     }

                                     NormalizationVisitor elem_vis(assumptions_);
                                     auto elem_node = std::make_shared<AddNode>(sum_ops);
                                     elem_node->accept(elem_vis);

                                     auto res_val = elem_vis.get_result();
                                     if (auto num = std::dynamic_pointer_cast<NumberNode>(res_val)) {
                                         if (std::holds_alternative<lmmc_real_t>(num->value)) {
                                              lmmc_real_t v = std::get<lmmc_real_t>(num->value);
                                             int eq_0, eq_1;
                                             lmmc_double_nearly_equal_tol(v, 0.0, 1e-10, 1e-10, &eq_0);
                                             lmmc_double_nearly_equal_tol(v, 1.0, 1e-10, 1e-10, &eq_1);
                                             if (eq_0) {
                                                 res_val = std::make_shared<NumberNode>(BigInt(0));
                                             } else if (eq_1) {
                                                 res_val = std::make_shared<NumberNode>(BigInt(1));
                                             }
                                         }
                                     }
                                     res_data.push_back(res_val);
                                 }
                             }

                             fused_ops.pop_back();
                             fused_ops.push_back(std::make_shared<MatrixNode>(R, C, res_data));
                             continue;
                         }
                    }
                }
                fused_ops.push_back(right);
            }

            if (!scalar_part->is_one()) {
                fused_ops.insert(fused_ops.begin(), scalar_part);
            }

            if (fused_ops.empty()) result = std::make_shared<NumberNode>(BigInt(1));
            else if (fused_ops.size() == 1) result = fused_ops[0];
            else result = std::make_shared<MultiplyNode>(fused_ops);

            return;
        }

        std::shared_ptr<NumberNode> const_acc = std::make_shared<NumberNode>(BigInt(1));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> bases;

        for (const auto& op : sc) {
             if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                 if (num->is_zero()) {
                     result = std::make_shared<NumberNode>(BigInt(0));
                     return;
                 }
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<SymbolicNode> base = op;
                 std::shared_ptr<NumberNode> exp = std::make_shared<NumberNode>(BigInt(1));
                 bool is_number_power = false;

                 if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                     // Only split base/exponent of a PowerNode if the exponent is a NumberNode.
                     // For symbolic exponents (e.g. 2^x), keep the PowerNode atomic so we
                     // don't silently drop the exponent when accumulating into `bases`.
                     auto e_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent);
                     if (e_num) {
                     base = pow->base;
                     exp = e_num;

                         if (auto b_num = std::dynamic_pointer_cast<NumberNode>(base)) {
                             long long exp_val = 0;
                             bool exp_ok = false;
                             bool exp_is_half = false;

                             if (std::holds_alternative<BigInt>(e_num->value)) {
                                 exp_val = (long long)std::get<BigInt>(e_num->value).to_double();
                                 exp_ok = true;
                                 } else if (std::holds_alternative<lmmc_real_t>(e_num->value)) {
                                      lmmc_real_t d = std::get<lmmc_real_t>(e_num->value);
                                     int eq_half;
                                     lmmc_double_nearly_equal_tol(d, 0.5, 1e-9, 1e-9, &eq_half);
                                     if (d == std::floor(d)) {
                                         exp_val = (long long)d;
                                         exp_ok = true;
                                     } else if (eq_half) {
                                         exp_is_half = true;
                                     }
                                 } else if (std::holds_alternative<Rational>(e_num->value)) {
                                 Rational r = std::get<Rational>(e_num->value);
                                 if (r.get_denominator() == BigInt(1)) {
                                     exp_val = (long long)r.get_numerator().to_double();
                                     exp_ok = true;
                                 } else if (r.get_numerator() == BigInt(1) && r.get_denominator() == BigInt(2)) {
                                     exp_is_half = true;
                                 }
                             }

                             if (exp_ok) {
                                 std::shared_ptr<NumberNode> pow_val = nullptr;

                                 if (exp_val == -1) {
                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          const auto& bi = std::get<BigInt>(b_num->value);
                                          if (!bi.is_zero()) pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), bi));
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational r = std::get<Rational>(b_num->value);
                                          if (!r.get_numerator().is_zero()) pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                      } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                                          lmmc_real_t bv = std::get<lmmc_real_t>(b_num->value);
                                          if (bv != 0.0) pow_val = std::make_shared<NumberNode>(1.0 / bv);
                                      }
                                 } else if (exp_val == 0) {
                                      pow_val = std::make_shared<NumberNode>(BigInt(1));
                                 } else if (std::abs(exp_val) > 0 && std::abs(exp_val) < 64) {

                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          BigInt b = std::get<BigInt>(b_num->value);
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational b = std::get<Rational>(b_num->value);
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                                          lmmc_real_t b = std::get<lmmc_real_t>(b_num->value);
                                          lmmc_real_t res = 1.0;
                                          for(int k=0;k<abs_exp;++k) res *= b;
                                           base_pow_val = std::make_shared<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {

                                              if (std::holds_alternative<BigInt>(base_pow_val->value)) {
                                                  const auto& bi = std::get<BigInt>(base_pow_val->value);
                                                  if (!bi.is_zero()) pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), bi));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value)) {
                                                  Rational r = std::get<Rational>(base_pow_val->value);
                                                  if (!r.get_numerator().is_zero()) pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              } else if (std::holds_alternative<lmmc_real_t>(base_pow_val->value)) {
                                                  lmmc_real_t bv = std::get<lmmc_real_t>(base_pow_val->value);
                                                  if (bv != 0.0) pow_val = std::make_shared<NumberNode>(1.0 / bv);
                                              }
                                          }
                                      }
                                 }

                                 if (pow_val) {
                                     const_acc = multiply_numbers(const_acc, pow_val);
                                     is_number_power = true;
                                 }
                             } else if (exp_is_half) {

                                 std::shared_ptr<NumberNode> root_val = nullptr;
                                 if (std::holds_alternative<BigInt>(b_num->value)) {

                                     double d = std::get<BigInt>(b_num->value).to_double();
                                     if (d >= 0) {
                                         double r = std::sqrt(d);
                                         int eq_r;
                                         lmmc_double_nearly_equal_tol(r, std::round(r), 1e-9, 1e-9, &eq_r);
                                         if (eq_r) {
                                             BigInt bi((long long)std::round(r));
                                             if (bi * bi == std::get<BigInt>(b_num->value)) {
                                                 root_val = std::make_shared<NumberNode>(bi);
                                             }
                                         }
                                     }
                                 } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                                     double d = std::get<lmmc_real_t>(b_num->value);
                                     if (d >= 0) root_val = std::make_shared<NumberNode>(std::sqrt(d));
                                 }

                                 if (root_val) {
                                     const_acc = multiply_numbers(const_acc, root_val);
                                     is_number_power = true;
                                 }
                             }
                         }
                     }
                 }

                 if (!is_number_power) {
                     auto it = bases.find(base);
                     if (it == bases.end()) {
                         bases[base] = exp;
                     } else {
                         bases[base] = add_numbers(it->second, exp);
                     }
                 }
             }
        }

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }

        std::vector<std::shared_ptr<SymbolicNode>> var_ops;
        for (auto const& [base, exp] : bases) {
            if (exp->is_zero()) {
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(std::make_shared<PowerNode>(base, exp));
            }
        }

        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());

        if (final_ops.empty()) result = std::make_shared<NumberNode>(BigInt(1));
        else if (final_ops.size() == 1) result = final_ops[0];
        else result = std::make_shared<MultiplyNode>(final_ops);
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto s_base = result;
        node.exponent->accept(*this);
        auto s_exp = result;

        if (s_exp->is_zero()) {
            result = std::make_shared<NumberNode>(BigInt(1));
            return;
        }
        if (s_exp->is_one()) {
            result = s_base;
            return;
        }

        if (s_base->is_zero()) {
             // 0^x = 0 only when x is provably positive; otherwise keep the
             // PowerNode so domain issues (e.g. 0^(-1/2)) are not silently
             // turned into 0. This mirrors SymbolicFactory::create_power.
             if (s_exp->is_positive()) {
                 result = std::make_shared<NumberNode>(BigInt(0));
                 return;
             }
             result = std::make_shared<PowerNode>(s_base, s_exp);
             return;
        }
        if (s_base->is_one()) {
            result = std::make_shared<NumberNode>(BigInt(1));
            return;
        }

        if (auto b_num = std::dynamic_pointer_cast<NumberNode>(s_base)) {
            if (auto e_num = std::dynamic_pointer_cast<NumberNode>(s_exp)) {
                 long long exp_val = 0;
                 bool exp_ok = false;
                 bool exp_is_half = false;

                 if (std::holds_alternative<BigInt>(e_num->value)) {
                     exp_val = (long long)std::get<BigInt>(e_num->value).to_double();
                     exp_ok = true;
                 } else if (std::holds_alternative<lmmc_real_t>(e_num->value)) {
                     lmmc_real_t d = std::get<lmmc_real_t>(e_num->value);
                     if (d == std::floor(d)) {
                         exp_val = (long long)d;
                         exp_ok = true;
                     } else if (std::abs(d - 0.5) < 1e-9) {
                         exp_is_half = true;
                     }
                 } else if (std::holds_alternative<Rational>(e_num->value)) {
                     Rational r = std::get<Rational>(e_num->value);
                     if (r.get_denominator() == BigInt(1)) {
                         exp_val = (long long)r.get_numerator().to_double();
                         exp_ok = true;
                     } else if (r.get_numerator() == BigInt(1) && r.get_denominator() == BigInt(2)) {
                         exp_is_half = true;
                     }
                 }

                 if (exp_ok) {
                     std::shared_ptr<NumberNode> pow_val = nullptr;
                     if (exp_val == -1) {
                          if (std::holds_alternative<BigInt>(b_num->value)) {
                              pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(b_num->value)));
                          } else if (std::holds_alternative<Rational>(b_num->value)) {
                              Rational r = std::get<Rational>(b_num->value);
                              if (!r.get_numerator().is_zero()) pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                          } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                              pow_val = std::make_shared<NumberNode>(1.0 / std::get<lmmc_real_t>(b_num->value));
                          }
                     } else if (exp_val == 0) {
                          pow_val = std::make_shared<NumberNode>(BigInt(1));
                     } else if (std::abs(exp_val) > 0 && std::abs(exp_val) < 64) {

                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          BigInt b = std::get<BigInt>(b_num->value);
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational b = std::get<Rational>(b_num->value);
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                                          lmmc_real_t b = std::get<lmmc_real_t>(b_num->value);
                                          lmmc_real_t res = 1.0;
                                          for(int k=0;k<abs_exp;++k) res *= b;
                                           base_pow_val = std::make_shared<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {

                                              if (std::holds_alternative<BigInt>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(base_pow_val->value)));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value)) {
                                                  Rational r = std::get<Rational>(base_pow_val->value);
                                                  pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              } else if (std::holds_alternative<lmmc_real_t>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(1.0 / std::get<lmmc_real_t>(base_pow_val->value));
                                              }
                                          }
                                      }
                     }

                     if (pow_val) {
                         result = pow_val;
                         return;
                     }
                 } else if (exp_is_half) {
                     if (std::holds_alternative<BigInt>(b_num->value)) {
                         double d = std::get<BigInt>(b_num->value).to_double();
                         if (d >= 0) {
                             double r = std::sqrt(d);
                             if (std::abs(r - std::round(r)) < 1e-9) {
                                 BigInt bi((long long)std::round(r));
                                 if (bi * bi == std::get<BigInt>(b_num->value)) {
                                     result = std::make_shared<NumberNode>(bi);
                                     return;
                                 }
                             }

                             BigInt val = std::get<BigInt>(b_num->value);
                             if (val > BigInt(0) && val < BigInt(1000000)) {
                                 long long v_ll = (long long)val.to_double();
                                 long long root = (long long)std::sqrt(v_ll);
                                 for (long long i = root; i >= 2; --i) {
                                     if (v_ll % (i*i) == 0) {
                                         long long s_ll = i;
                                         long long k_ll = v_ll / (i*i);

                                         auto s_node = std::make_shared<NumberNode>(BigInt(s_ll));
                                         auto k_node = std::make_shared<NumberNode>(BigInt(k_ll));
                                         auto half_node = std::make_shared<NumberNode>(Rational(1, 2));
                                         auto pow_node = std::make_shared<PowerNode>(k_node, half_node);

                                         result = SymbolicFactory::create_multiply({s_node, pow_node});
                                         return;
                                     }
                                 }
                             }
                         }
                     } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                         double d = std::get<lmmc_real_t>(b_num->value);
                         if (d >= 0) {
                             result = std::make_shared<NumberNode>(std::sqrt(d));
                             return;
                         }
                     }
                 }

                 // Simplify powers of imaginary unit: (-1)^(n/2) for odd n
                 // i² = -1, i³ = -i, i⁴ = 1 (handled via (-1)^(n/2) reduction)
                 if (std::holds_alternative<Rational>(e_num->value)) {
                     Rational exp_r = std::get<Rational>(e_num->value);
                     bool base_is_neg_one = false;
                     if (std::holds_alternative<BigInt>(b_num->value)) {
                         base_is_neg_one = (std::get<BigInt>(b_num->value) == BigInt(-1));
                     } else if (std::holds_alternative<Rational>(b_num->value)) {
                         base_is_neg_one = (std::get<Rational>(b_num->value) == Rational(-1));
                     } else if (std::holds_alternative<lmmc_real_t>(b_num->value)) {
                         int eq;
                         lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(b_num->value), -1.0, 1e-12, 1e-12, &eq);
                         base_is_neg_one = (eq != 0);
                     }

                     if (base_is_neg_one && exp_r.get_denominator() == BigInt(2)) {
                         // (-1)^(n/2) where n is the numerator
                         long long n_val = (long long)exp_r.get_numerator().to_double();
                         // Reduce n mod 4
                         long long r_mod = ((n_val % 4) + 4) % 4;
                         // i = (-1)^(1/2)
                         auto i_node = std::make_shared<PowerNode>(
                             std::make_shared<NumberNode>(BigInt(-1)),
                             std::make_shared<NumberNode>(Rational(1, 2)));
                         if (r_mod == 0) {
                             // i⁴ = 1
                             result = std::make_shared<NumberNode>(BigInt(1));
                             return;
                         } else if (r_mod == 1) {
                             // i¹ = i (keep as (-1)^(1/2))
                             result = i_node;
                             return;
                         } else if (r_mod == 2) {
                             // i² = -1
                             result = std::make_shared<NumberNode>(BigInt(-1));
                             return;
                         } else { // r_mod == 3
                             // i³ = -i
                             std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {
                                 std::make_shared<NumberNode>(BigInt(-1)), i_node
                             };
                             result = std::make_shared<MultiplyNode>(mul_ops);
                             return;
                         }
                     }
                 }
            }
        }

        if (auto m_base = std::dynamic_pointer_cast<MultiplyNode>(s_base)) {
            std::vector<std::shared_ptr<SymbolicNode>> new_ops;
            for(auto& op : m_base->operands) {

                auto term_pow = std::make_shared<PowerNode>(op, s_exp);
                term_pow->accept(*this);

                if (auto mul_res = std::dynamic_pointer_cast<MultiplyNode>(result)) {
                     new_ops.insert(new_ops.end(), mul_res->operands.begin(), mul_res->operands.end());
                } else {
                     new_ops.push_back(result);
                }
            }

            auto final_mul = std::make_shared<MultiplyNode>(new_ops);
            final_mul->accept(*this);
            return;
        }

        if (auto p_base = std::dynamic_pointer_cast<PowerNode>(s_base)) {

             std::vector<std::shared_ptr<SymbolicNode>> exp_ops;
             exp_ops.push_back(p_base->exponent);
             exp_ops.push_back(s_exp);

             auto mul_exp = std::make_shared<MultiplyNode>(exp_ops);
             mul_exp->accept(*this);

             auto new_pow = std::make_shared<PowerNode>(p_base->base, result);
             new_pow->accept(*this);
             return;
        }

        result = std::make_shared<PowerNode>(s_base, s_exp);
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> s_args;
        for(auto& a : node.arguments) {
            a->accept(*this);
            s_args.push_back(result);
        }

        if (s_args.size() == 1) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(s_args[0])) {
                lmmc_real_t val = 0.0;
                if (std::holds_alternative<lmmc_real_t>(num->value)) val = std::get<lmmc_real_t>(num->value);
                else if (std::holds_alternative<BigInt>(num->value)) val = (lmmc_real_t)std::get<BigInt>(num->value).to_double();
                else if (std::holds_alternative<Rational>(num->value)) val = (lmmc_real_t)std::get<Rational>(num->value).to_double();

                switch (node.type) {
                    case FunctionNode::FuncType::Sin:
                    {
                        int eq; lmmc_double_nearly_equal_tol(val, 0.0, 1e-12, 1e-12, &eq);
                        if (eq) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                        break;
                    }
                    case FunctionNode::FuncType::Cos:
                    {
                        int eq; lmmc_double_nearly_equal_tol(val, 0.0, 1e-12, 1e-12, &eq);
                        if (eq) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
                        break;
                    }
                    case FunctionNode::FuncType::Tan:
                    {
                        int eq; lmmc_double_nearly_equal_tol(val, 0.0, 1e-12, 1e-12, &eq);
                        if (eq) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                        break;
                    }
                    case FunctionNode::FuncType::Exp:
                    {
                        int eq; lmmc_double_nearly_equal_tol(val, 0.0, 1e-12, 1e-12, &eq);
                        if (eq) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
                        break;
                    }
                    case FunctionNode::FuncType::Ln:
                    {
                        int eq1, eq0;
                        lmmc_double_nearly_equal_tol(val, 1.0, 1e-12, 1e-12, &eq1);
                        if (eq1) { result = std::make_shared<NumberNode>(BigInt(0)); return; }

                        lmmc_double_nearly_equal_tol(val, 0.0, 1e-12, 1e-12, &eq0);
                        if (eq0) {
                             std::vector<std::shared_ptr<SymbolicNode>> inf_args;
                             auto inf = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                             std::vector<std::shared_ptr<SymbolicNode>> m_args = {std::make_shared<NumberNode>(BigInt(-1)), inf};
                             result = std::make_shared<MultiplyNode>(m_args);
                             return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Log:
                        if (s_args.size() == 2) {

                             std::vector<std::shared_ptr<SymbolicNode>> args_x = { s_args[0] };
                             std::vector<std::shared_ptr<SymbolicNode>> args_b = { s_args[1] };
                             auto ln_x = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args_x);
                             auto ln_b = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args_b);
                             auto ln_b_inv = std::make_shared<PowerNode>(ln_b, std::make_shared<NumberNode>(BigInt(-1)));
                             std::vector<std::shared_ptr<SymbolicNode>> m_args = { ln_x, ln_b_inv };
                             auto prod = std::make_shared<MultiplyNode>(m_args);

                             NormalizationVisitor v(assumptions_);
                             prod->accept(v);
                             result = v.get_result();
                             return;
                        }
                        break;
                    case FunctionNode::FuncType::Abs:
                        result = std::make_shared<NumberNode>(std::abs(val));
                        return;
                    case FunctionNode::FuncType::Sqrt:
                    {
                        BigInt n;
                        bool check_int = false;

                        if (std::holds_alternative<BigInt>(num->value)) {
                            n = std::get<BigInt>(num->value);
                            check_int = true;
                        } else if (std::holds_alternative<Rational>(num->value)) {
                            Rational r = std::get<Rational>(num->value);
                            if (r.get_denominator() == BigInt(1)) {
                                n = r.get_numerator();
                                check_int = true;
                            }
                        }

                        if (check_int) {
                            if (n >= BigInt(0)) {
                                BigInt s = n.sqrt();
                                if (s * s == n) {
                                    result = std::make_shared<NumberNode>(s);
                                    return;
                                }

                                BigInt one(1);
                                BigInt coeff = one;
                                BigInt rem = n;

                                long long factors[] = {2, 3, 5, 7};
                                for (long long f : factors) {
                                    BigInt sq(f*f);
                                    while (rem % sq == BigInt(0)) {
                                        rem = rem / sq;
                                        coeff = coeff * BigInt(f);
                                    }
                                }

                                if (coeff > one) {
                                     std::vector<std::shared_ptr<SymbolicNode>> inner_args = { std::make_shared<NumberNode>(rem) };
                                     auto inner_sqrt = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sqrt, inner_args);
                                     std::vector<std::shared_ptr<SymbolicNode>> mul_args = { std::make_shared<NumberNode>(coeff), inner_sqrt };
                                     result = std::make_shared<MultiplyNode>(mul_args);
                                     return;
                                }

                            }
                        }

                        if (val >= 0) {
                             double sq = std::sqrt(val);
                             int eq_sq;
                             lmmc_double_nearly_equal_tol(sq, std::round(sq), 1e-12, 1e-12, &eq_sq);
                             if (eq_sq) {
                                 result = std::make_shared<NumberNode>(BigInt((long long)std::round(sq)));
                             } else {
                                 result = std::make_shared<NumberNode>(sq);
                             }
                             return;
                        }
                    }
                        break;
                     case FunctionNode::FuncType::LambertW:
                     {

                         lmmc_real_t w_res;
                         if (lmmc_lambertw(val, &w_res) == LMMC_STATUS_OK) {
                             result = std::make_shared<NumberNode>(w_res);
                             return;
                         }

                         break;
                     }
                     default: break;
                }
            }
        }

        if (node.type == FunctionNode::FuncType::Ln && s_args.size() == 1) {
             if (auto pow = std::dynamic_pointer_cast<PowerNode>(s_args[0])) {
                  auto y = pow->exponent;
                  auto x = pow->base;

                  std::vector<std::shared_ptr<SymbolicNode>> ln_args = { x };
                  auto ln_x = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, ln_args);

                  std::vector<std::shared_ptr<SymbolicNode>> m_args = { y, ln_x };
                  auto prod = std::make_shared<MultiplyNode>(m_args);

                  NormalizationVisitor v(assumptions_);
                  prod->accept(v);
                  result = v.get_result();
                  return;
             }
             if (auto func = std::dynamic_pointer_cast<FunctionNode>(s_args[0])) {
                 if (func->type == FunctionNode::FuncType::Exp && func->arguments.size() == 1) {
                     result = func->arguments[0];
                     return;
                 }
             }
        }

        if (node.type == FunctionNode::FuncType::Log && s_args.size() == 2) {
             std::vector<std::shared_ptr<SymbolicNode>> args_x = { s_args[0] };
             std::vector<std::shared_ptr<SymbolicNode>> args_b = { s_args[1] };
             auto ln_x = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args_x);
             auto ln_b = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args_b);
             auto ln_b_inv = std::make_shared<PowerNode>(ln_b, std::make_shared<NumberNode>(BigInt(-1)));
             std::vector<std::shared_ptr<SymbolicNode>> m_args = { ln_x, ln_b_inv };
             auto prod = std::make_shared<MultiplyNode>(m_args);
             NormalizationVisitor v(assumptions_);
             prod->accept(v);
             result = v.get_result();
             return;
        }

        if (s_args.size() == 1) {
            Rational k_pi_val;
            if (get_pi_coeff(s_args[0], k_pi_val)) {
                BigInt n = k_pi_val.get_numerator();
                BigInt d = k_pi_val.get_denominator();

                BigInt two_d = d * BigInt(2);
                BigInt reduced_n = n % two_d;
                if (reduced_n.IsNegative()) reduced_n = reduced_n + two_d;

                auto zero = std::make_shared<NumberNode>(BigInt(0));
                auto one = std::make_shared<NumberNode>(BigInt(1));
                auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
                auto half = std::make_shared<NumberNode>(Rational(1, 2));
                auto minus_half = std::make_shared<NumberNode>(Rational(-1, 2));

                auto root2 = std::make_shared<PowerNode>(std::make_shared<NumberNode>(BigInt(2)), std::make_shared<NumberNode>(Rational(1, 2)));
                std::vector<std::shared_ptr<SymbolicNode>> half_root2_args = { half, root2 };
                auto half_root2 = std::make_shared<MultiplyNode>(half_root2_args);
                std::vector<std::shared_ptr<SymbolicNode>> minus_half_root2_args = { minus_half, root2 };
                auto minus_half_root2 = std::make_shared<MultiplyNode>(minus_half_root2_args);

                auto root3 = std::make_shared<PowerNode>(std::make_shared<NumberNode>(BigInt(3)), std::make_shared<NumberNode>(Rational(1, 2)));
                std::vector<std::shared_ptr<SymbolicNode>> half_root3_args = { half, root3 };
                auto half_root3 = std::make_shared<MultiplyNode>(half_root3_args);
                std::vector<std::shared_ptr<SymbolicNode>> minus_half_root3_args = { minus_half, root3 };
                auto minus_half_root3 = std::make_shared<MultiplyNode>(minus_half_root3_args);

                std::vector<std::shared_ptr<SymbolicNode>> minus_root3_args = { minus_one, root3 };
                auto minus_root3 = std::make_shared<MultiplyNode>(minus_root3_args);

                auto third = std::make_shared<NumberNode>(Rational(1, 3));
                std::vector<std::shared_ptr<SymbolicNode>> third_root3_args = { third, root3 };
                auto third_root3 = std::make_shared<MultiplyNode>(third_root3_args);
                auto minus_third = std::make_shared<NumberNode>(Rational(-1, 3));
                std::vector<std::shared_ptr<SymbolicNode>> minus_third_root3_args = { minus_third, root3 };
                auto minus_third_root3 = std::make_shared<MultiplyNode>(minus_third_root3_args);

                if (d == BigInt(1)) {

                    if (node.type == FunctionNode::FuncType::Sin || node.type == FunctionNode::FuncType::Tan) {
                        result = zero; return;
                    } else if (node.type == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(0)) result = one;
                        else result = minus_one;
                        return;
                    }
                } else if (d == BigInt(2)) {

                    if (node.type == FunctionNode::FuncType::Sin) {
                        if (reduced_n == BigInt(1)) result = one;
                        else if (reduced_n == BigInt(3)) result = minus_one;
                        return;
                    } else if (node.type == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(3)) {
                            result = zero; return;
                        }
                    }

                } else if (d == BigInt(3)) {

                    if (node.type == FunctionNode::FuncType::Sin) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(2)) result = half_root3;
                         else if (reduced_n == BigInt(4) || reduced_n == BigInt(5)) result = minus_half_root3;
                         return;
                    } else if (node.type == FunctionNode::FuncType::Cos) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = half;
                         else if (reduced_n == BigInt(2) || reduced_n == BigInt(4)) result = minus_half;
                         return;
                    } else if (node.type == FunctionNode::FuncType::Tan) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(4)) result = root3;
                         else if (reduced_n == BigInt(2) || reduced_n == BigInt(5)) result = minus_root3;
                         return;
                    }
                } else if (d == BigInt(4)) {

                     if (node.type == FunctionNode::FuncType::Sin) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(3)) result = half_root2;
                         else if (reduced_n == BigInt(5) || reduced_n == BigInt(7)) result = minus_half_root2;
                         return;
                    } else if (node.type == FunctionNode::FuncType::Cos) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(7)) result = half_root2;
                         else if (reduced_n == BigInt(3) || reduced_n == BigInt(5)) result = minus_half_root2;
                         return;
                    } else if (node.type == FunctionNode::FuncType::Tan) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = one;
                         else if (reduced_n == BigInt(3) || reduced_n == BigInt(7)) result = minus_one;
                         return;
                    }
                } else if (d == BigInt(6)) {

                    if (node.type == FunctionNode::FuncType::Sin) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = half;
                        else if (reduced_n == BigInt(7) || reduced_n == BigInt(11)) result = minus_half;
                        return;
                    } else if (node.type == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(11)) result = half_root3;
                        else if (reduced_n == BigInt(5) || reduced_n == BigInt(7)) result = minus_half_root3;
                        return;
                    } else if (node.type == FunctionNode::FuncType::Tan) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(7)) result = third_root3;
                        else if (reduced_n == BigInt(5) || reduced_n == BigInt(11)) result = minus_third_root3;
                        return;
                    }
                }
            }
        }

        if (s_args.size() == 1) {
             std::shared_ptr<SymbolicNode> pos_arg = nullptr;

             if (check_negative_arg(s_args[0], pos_arg)) {
                 if (node.type == FunctionNode::FuncType::Sin ||
                     node.type == FunctionNode::FuncType::Tan ||
                     node.type == FunctionNode::FuncType::ArcTan ||
                     node.type == FunctionNode::FuncType::ArcSin) {

                     std::vector<std::shared_ptr<SymbolicNode>> new_args = { pos_arg };
                     auto new_func = std::make_shared<FunctionNode>(node.type, new_args);
                     std::vector<std::shared_ptr<SymbolicNode>> mul_ops = { std::make_shared<NumberNode>(BigInt(-1)), new_func };
                     result = std::make_shared<MultiplyNode>(mul_ops);
                     return;
                 } else if (node.type == FunctionNode::FuncType::Cos) {

                     std::vector<std::shared_ptr<SymbolicNode>> new_args = { pos_arg };
                     result = std::make_shared<FunctionNode>(node.type, new_args);
                     return;
                 }
             }
        }

        if (node.type == FunctionNode::FuncType::Log && s_args.size() == 2) {
             auto val = s_args[0];
             auto base = s_args[1];

             if (val->equals(*base)) {
                 result = std::make_shared<NumberNode>(BigInt(1));
                 return;
             }

             if (val->is_one()) {
                 result = std::make_shared<NumberNode>(BigInt(0));
                 return;
             }

             if (auto pow = std::dynamic_pointer_cast<PowerNode>(val)) {
                 if (pow->base->equals(*base)) {
                     result = pow->exponent;
                     return;
                 }
             }
        }

        // Attempt assumption-based simplification before falling through
        auto func_node = std::make_shared<FunctionNode>(node.type, s_args);
        if (auto simplified = try_assumption_simplify(func_node)) {
            // Recursively normalize the simplified result
            simplified->accept(*this);
            return;
        }
        result = func_node;
    }

    void visit(MatrixNode& node) override {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage)) {
             auto& dense = std::get<MatrixNode::DenseStorage>(node.storage);
             MatrixNode::DenseStorage new_dense;
             for(auto& item : dense) {
                 if(item) {
                    item->accept(*this);
                    new_dense.push_back(result);
                 } else {
                    new_dense.push_back(nullptr);
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_dense);
        } else {
             auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage);
             MatrixNode::SparseStorage new_sparse;
             for(auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     new_sparse[idx] = result;
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_sparse);
        }
    }

    void visit(RelationalNode& node) override {

        std::shared_ptr<SymbolicNode> new_left = nullptr;
        std::shared_ptr<SymbolicNode> new_right = nullptr;

        if (node.left) {
            node.left->accept(*this);
            new_left = result;
        }
        if (node.right) {

            node.right->accept(*this);
            new_right = result;
        }

        if (!new_left) new_left = node.left;
        if (!new_right) new_right = node.right;

        result = std::make_shared<RelationalNode>(new_left, new_right, node.op);
    }

    void visit(LogicalNode& node) override {
        // Implication: A ⇒ B = ¬A ∨ B
        if (node.op == LogicalNode::Op::Implies) {
            auto not_left = std::make_shared<LogicalNode>(node.left, nullptr, LogicalNode::Op::Not);
            auto or_node = std::make_shared<LogicalNode>(not_left, node.right, LogicalNode::Op::Or);
            or_node->accept(*this);
            return;
        }

        // NOT handling: De Morgan's laws and double negation
        if (node.op == LogicalNode::Op::Not) {
            // Normalize the operand first
            std::shared_ptr<SymbolicNode> new_left = nullptr;
            if (node.left) {
                node.left->accept(*this);
                new_left = result;
            }
            if (!new_left) new_left = node.left;

            // Double negation: ¬(¬A) = A
            if (auto inner_logical = std::dynamic_pointer_cast<LogicalNode>(new_left)) {
                if (inner_logical->op == LogicalNode::Op::Not) {
                    result = inner_logical->left;
                    return;
                }
                // De Morgan's law: ¬(A∧B) = ¬A∨¬B
                if (inner_logical->op == LogicalNode::Op::And) {
                    auto not_a = std::make_shared<LogicalNode>(inner_logical->left, nullptr, LogicalNode::Op::Not);
                    auto not_b = std::make_shared<LogicalNode>(inner_logical->right, nullptr, LogicalNode::Op::Not);
                    auto or_node = std::make_shared<LogicalNode>(not_a, not_b, LogicalNode::Op::Or);
                    or_node->accept(*this);
                    return;
                }
                // De Morgan's law: ¬(A∨B) = ¬A∧¬B
                if (inner_logical->op == LogicalNode::Op::Or) {
                    auto not_a = std::make_shared<LogicalNode>(inner_logical->left, nullptr, LogicalNode::Op::Not);
                    auto not_b = std::make_shared<LogicalNode>(inner_logical->right, nullptr, LogicalNode::Op::Not);
                    auto and_node = std::make_shared<LogicalNode>(not_a, not_b, LogicalNode::Op::And);
                    and_node->accept(*this);
                    return;
                }
            }

            result = std::make_shared<LogicalNode>(new_left, nullptr, LogicalNode::Op::Not);
            return;
        }

        // And / Or: normalize operands
        std::shared_ptr<SymbolicNode> new_left = nullptr;
        std::shared_ptr<SymbolicNode> new_right = nullptr;

        if (node.left) {
            node.left->accept(*this);
            new_left = result;
        }
        if (node.right) {
            node.right->accept(*this);
            new_right = result;
        }

        if (!new_left) new_left = node.left;
        if (!new_right) new_right = node.right;

        result = std::make_shared<LogicalNode>(new_left, new_right, node.op);
    }

    void visit(PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> new_branches;
        new_branches.reserve(node.branches.size());

        for (const auto& b : node.branches) {
            // Normalize expression
            b.expression->accept(*this);
            auto new_expr = result;

            // Normalize condition
            b.condition->accept(*this);
            auto new_cond = result;

            // Validate condition is RelationalNode or LogicalNode
            // (keep it regardless, but this ensures normalization is applied)
            new_branches.push_back({new_expr, new_cond});
        }

        std::shared_ptr<SymbolicNode> new_default = nullptr;
        if (node.default_expr) {
            node.default_expr->accept(*this);
            new_default = result;
        }

        result = std::make_shared<PiecewiseNode>(std::move(new_branches), new_default);
    }

    void visit(SummationNode& node) override {
        node.body->accept(*this);
        auto new_body = result;

        node.lower_bound->accept(*this);
        auto new_lower = result;

        node.upper_bound->accept(*this);
        auto new_upper = result;

        // 当上下界均为具体整数且范围较小时，展开求和为显式和。
        auto lo_n = std::dynamic_pointer_cast<NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<NumberNode>(new_upper);
        if (lo_n && hi_n && std::holds_alternative<BigInt>(lo_n->value)
            && std::holds_alternative<BigInt>(hi_n->value)) {
            long long lo = (long long)std::get<BigInt>(lo_n->value).to_int();
            long long hi = (long long)std::get<BigInt>(hi_n->value).to_int();
            if (hi < lo) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
            if (hi - lo < 1000) {
                std::vector<std::shared_ptr<SymbolicNode>> terms;
                for (long long kk = lo; kk <= hi; ++kk) {
                    auto kval = std::make_shared<NumberNode>(BigInt((long long)kk));
                    auto term = norm_subst_index(new_body, node.index_var, kval);
                    NormalizationVisitor inner;
                    term->accept(inner);
                    terms.push_back(inner.get_result());
                }
                if (terms.empty()) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                auto sum_node = std::make_shared<AddNode>(terms);
                sum_node->accept(*this);
                return;
            }
        }

        result = std::make_shared<SummationNode>(new_body, node.index_var, new_lower, new_upper);
    }

    void visit(ProductNode_Op& node) override {
        node.body->accept(*this);
        auto new_body = result;

        node.lower_bound->accept(*this);
        auto new_lower = result;

        node.upper_bound->accept(*this);
        auto new_upper = result;

        auto lo_n = std::dynamic_pointer_cast<NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<NumberNode>(new_upper);
        if (lo_n && hi_n && std::holds_alternative<BigInt>(lo_n->value)
            && std::holds_alternative<BigInt>(hi_n->value)) {
            long long lo = (long long)std::get<BigInt>(lo_n->value).to_int();
            long long hi = (long long)std::get<BigInt>(hi_n->value).to_int();
            if (hi < lo) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
            if (hi - lo < 1000) {
                std::vector<std::shared_ptr<SymbolicNode>> factors;
                for (long long kk = lo; kk <= hi; ++kk) {
                    auto kval = std::make_shared<NumberNode>(BigInt((long long)kk));
                    auto term = norm_subst_index(new_body, node.index_var, kval);
                    NormalizationVisitor inner;
                    term->accept(inner);
                    factors.push_back(inner.get_result());
                }
                if (factors.empty()) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
                auto prod_node = std::make_shared<MultiplyNode>(factors);
                prod_node->accept(*this);
                return;
            }
        }

        result = std::make_shared<ProductNode_Op>(new_body, node.index_var, new_lower, new_upper);
    }

    void visit(TransformNode& node) override {
        node.body->accept(*this);
        auto new_body = result;

        result = std::make_shared<TransformNode>(node.transform_type, new_body, node.source_var, node.target_var);
    }

    void visit(QuantifierNode& node) override {
        node.domain->accept(*this);
        auto new_domain = result;

        node.predicate->accept(*this);
        auto new_predicate = result;

        // Simplify ∀x∈S: true → true
        if (node.quantifier_type == QuantifierNode::Type::ForAll) {
            if (new_predicate->is_one()) {
                result = std::make_shared<NumberNode>(BigInt(1));
                return;
            }
        }

        // Simplify ∃x∈S: false → false
        if (node.quantifier_type == QuantifierNode::Type::Exists) {
            if (new_predicate->is_zero()) {
                result = std::make_shared<NumberNode>(BigInt(0));
                return;
            }
        }

        result = std::make_shared<QuantifierNode>(node.quantifier_type, node.bound_var, new_domain, new_predicate);
    }

    void visit(SetBuilderNode& node) override {
        node.domain->accept(*this);
        auto new_domain = result;

        node.predicate->accept(*this);
        auto new_predicate = result;

        result = std::make_shared<SetBuilderNode>(node.element_var, new_domain, new_predicate);
    }

private:
    const lamina::AssumptionContext* assumptions_ = nullptr;

    /**
     * @brief Attempt assumption-based simplification on a node.
     *
     * Applies the following rules when an AssumptionContext is available:
     * - sqrt(x²) → x when x is NonNegative
     * - sqrt(x²) → abs(x) when x is Real (but not NonNegative)
     * - abs(x) → x when x is Positive
     * - abs(x) → -x when x is Negative
     *
     * @param node The node to attempt simplification on
     * @return Simplified node if a rule applied, or nullptr if no rule matched
     */
    std::shared_ptr<SymbolicNode> try_assumption_simplify(
        const std::shared_ptr<SymbolicNode>& node) {
        if (!assumptions_) return nullptr;

        auto func = std::dynamic_pointer_cast<FunctionNode>(node);
        if (!func || func->arguments.size() != 1) return nullptr;

        const auto& arg = func->arguments[0];

        // Rule: sqrt(x²) → x when x is NonNegative
        // Rule: sqrt(x²) → abs(x) when x is Real (but not NonNegative)
        if (func->type == FunctionNode::FuncType::Sqrt) {
            // Check if argument is a PowerNode with exponent 2
            auto pow = std::dynamic_pointer_cast<PowerNode>(arg);
            if (pow) {
                auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent);
                if (exp_num) {
                    bool is_exp_two = false;
                    if (std::holds_alternative<BigInt>(exp_num->value)) {
                        is_exp_two = (std::get<BigInt>(exp_num->value) == BigInt(2));
                    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                        int eq;
                        lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(exp_num->value), 2.0, 1e-12, 1e-12, &eq);
                        is_exp_two = (eq != 0);
                    } else if (std::holds_alternative<Rational>(exp_num->value)) {
                        is_exp_two = (std::get<Rational>(exp_num->value) == Rational(2));
                    }

                    if (is_exp_two) {
                        // We have sqrt(base²) — query the base's properties
                        SymbolicExpr base_expr(pow->base);

                        lamina::Tribool nonneg = assumptions_->is_nonnegative(base_expr);
                        if (nonneg == lamina::Tribool::True) {
                            // sqrt(x²) → x when x is NonNegative
                            return pow->base;
                        }

                        lamina::Tribool real = assumptions_->is_real(base_expr);
                        if (real == lamina::Tribool::True) {
                            // sqrt(x²) → abs(x) when x is Real (but not NonNegative)
                            std::vector<std::shared_ptr<SymbolicNode>> abs_args = { pow->base };
                            return std::make_shared<FunctionNode>(FunctionNode::FuncType::Abs, abs_args);
                        }
                    }
                }
            }
        }

        // Rule: abs(x) → x when x is Positive
        // Rule: abs(x) → -x when x is Negative
        if (func->type == FunctionNode::FuncType::Abs) {
            SymbolicExpr arg_expr(arg);

            lamina::Tribool pos = assumptions_->is_positive(arg_expr);
            if (pos == lamina::Tribool::True) {
                // abs(x) → x when x is Positive
                return arg;
            }

            lamina::Tribool neg = assumptions_->is_negative(arg_expr);
            if (neg == lamina::Tribool::True) {
                // abs(x) → -x when x is Negative
                std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {
                    std::make_shared<NumberNode>(BigInt(-1)), arg
                };
                return std::make_shared<MultiplyNode>(mul_ops);
            }
        }

        return nullptr;
    }
};
