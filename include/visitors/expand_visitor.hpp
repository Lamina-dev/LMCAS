/**
 * @file expand_visitor.hpp
 * @brief 展开访问器，将乘积和幂展开为多项式形式。
 */
#pragma once

#include "../symbolic_ast.hpp"
#include "normalization_visitor.hpp"

/** @brief 展开访问器，将乘积分配律和整数次幂展开为多项式形式 */
class ExpandVisitor : public SymbolicVisitor {
public:
    std::shared_ptr<SymbolicNode> result;  ///< 展开结果节点

    ExpandVisitor() : result(nullptr) {}

    /**
     * @brief 获取展开结果
     * @return 展开后的 AST 节点
     */
    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override { result = node.clone(); }
    void visit(VariableNode& node) override { result = node.clone(); }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = std::make_shared<AddNode>(new_ops);
    }

    void visit(MultiplyNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> expanded_ops;
        for (auto& op : node.operands) {
            op->accept(*this);
            expanded_ops.push_back(result);
        }

        if (expanded_ops.empty()) {
            result = std::make_shared<NumberNode>(BigInt(1));
            return;
        }

        NormalizationVisitor norm;
        std::shared_ptr<SymbolicNode> current = expanded_ops[0];
        for (size_t i = 1; i < expanded_ops.size(); ++i) {
            current = norm.expand_product(current, expanded_ops[i]);
        }
        result = current;
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto expanded_base = result;
        node.exponent->accept(*this);
        auto expanded_exp = result;

        if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(expanded_exp)) {
            long long e_val = -1;
            if (std::holds_alternative<BigInt>(exp_num->value)) e_val = std::get<BigInt>(exp_num->value).to_int();
            else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) e_val = (long long)std::get<lmmc_real_t>(exp_num->value);
            else if (std::holds_alternative<Rational>(exp_num->value)) e_val = (long long)std::get<Rational>(exp_num->value).to_double();

            if (e_val == 0) {
                result = std::make_shared<NumberNode>(BigInt(1));
                return;
            }
            if (e_val > 0 && e_val < 20) {
                NormalizationVisitor norm;
                std::shared_ptr<SymbolicNode> current = expanded_base;
                for (int i = 1; i < (int)e_val; ++i) {
                    current = norm.expand_product(current, expanded_base);
                }
                result = current;
                return;
            }
            // Exponent too large for expansion (>= 20) or negative — keep as PowerNode
        }
        result = std::make_shared<PowerNode>(expanded_base, expanded_exp);
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (auto& arg : node.arguments) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = std::make_shared<FunctionNode>(node.type, new_args);
    }

    void visit(MatrixNode& node) override {
        result = node.clone();
    }

    void visit(RelationalNode& node) override {
        node.left->accept(*this);
        auto head = result;
        node.right->accept(*this);
        auto tail = result;
        result = std::make_shared<RelationalNode>(head, tail, node.op);
    }

    void visit(LogicalNode& node) override {
        node.left->accept(*this);
        auto head = result;
        node.right->accept(*this);
        auto tail = result;
        result = std::make_shared<LogicalNode>(head, tail, node.op);
    }

    void visit(ComplexNode& node) override {
        node.real->accept(*this);
        auto expanded_r = result;
        node.imag->accept(*this);
        auto expanded_i = result;
        result = SymbolicFactory::create_complex(expanded_r, expanded_i);
    }
};
