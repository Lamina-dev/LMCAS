/**
 * @file expand_visitor.hpp
 * @brief 展开访问器，将乘积和幂展开为多项式形式。
 */
#pragma once

#include "../symbolic_ast.hpp"
#include "normalization_visitor.hpp"

namespace LMCAS {

/** @brief 展开访问器，将乘积分配律和整数次幂展开为多项式形式 */
class ExpandVisitor : public LMCAS::detail::SymbolicVisitor {
public:
    std::shared_ptr<const SymbolicNode> result;  ///< 展开结果节点

    ExpandVisitor() : result(nullptr) {}

    /**
     * @brief 获取展开结果
     * @return 展开后的 AST 节点
     */
    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    void visit(const NumberNode& node) override { result = node.clone(); }
    void visit(const VariableNode& node) override { result = node.clone(); }

    void visit(const AddNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = LMCAS::detail::make_node<AddNode>(new_ops);
    }

    void visit(const MultiplyNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> expanded_ops;
        for (auto& op : node.operands()) {
            op->accept(*this);
            expanded_ops.push_back(result);
        }

        if (expanded_ops.empty()) {
            result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
            return;
        }

        NormalizationVisitor norm;
        std::shared_ptr<const SymbolicNode> current = expanded_ops[0];
        for (size_t i = 1; i < expanded_ops.size(); ++i) {
            current = norm.expand_product(current, expanded_ops[i]);
        }
        result = current;
    }

    void visit(const PowerNode& node) override {
        node.base()->accept(*this);
        auto expanded_base = result;
        node.exponent()->accept(*this);
        auto expanded_exp = result;

        if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(expanded_exp)) {
            long long e_val = -1;
            if (std::holds_alternative<BigInt>(exp_num->value())) e_val = std::get<BigInt>(exp_num->value()).to_int();
            else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) e_val = (long long)std::get<lmmc_real_t>(exp_num->value());
            else if (std::holds_alternative<Rational>(exp_num->value())) e_val = (long long)std::get<Rational>(exp_num->value()).to_double();

            if (e_val == 0) {
                result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                return;
            }
            if (e_val > 0 && e_val < 20) {
                NormalizationVisitor norm;
                std::shared_ptr<const SymbolicNode> current = expanded_base;
                for (int i = 1; i < (int)e_val; ++i) {
                    current = norm.expand_product(current, expanded_base);
                }
                result = current;
                return;
            }
            // Exponent too large for expansion (>= 20) or negative — keep as PowerNode
        }
        result = LMCAS::detail::make_node<PowerNode>(expanded_base, expanded_exp);
    }

    void visit(const FunctionNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& arg : node.arguments()) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = LMCAS::detail::make_node<FunctionNode>(node.type(), new_args);
    }
    void visit(const UninterpretedFunctionNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> arguments;
        arguments.reserve(node.arguments().size());
        for (const auto& argument : node.arguments()) {
            argument->accept(*this);
            arguments.push_back(result);
        }
        result = LMCAS::detail::make_node<UninterpretedFunctionNode>(
            node.name(), std::move(arguments));
    }

    void visit(const MatrixNode& node) override {
        result = node.clone();
    }
    void visit(const IntegralNode& node) override {
        node.body()->accept(*this);
        auto body = result;
        std::shared_ptr<const SymbolicNode> lower;
        std::shared_ptr<const SymbolicNode> upper;
        if (node.lower()) {
            node.lower()->accept(*this);
            lower = result;
            node.upper()->accept(*this);
            upper = result;
        }
        result = LMCAS::detail::make_node<IntegralNode>(
            body, node.variable(), lower, upper);
    }

    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        auto head = result;
        node.right()->accept(*this);
        auto tail = result;
        result = LMCAS::detail::make_node<RelationalNode>(head, tail, node.op());
    }

    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        auto head = result;
        std::shared_ptr<const SymbolicNode> tail = nullptr;
        if (node.right()) {
            node.right()->accept(*this);
            tail = result;
        }
        result = LMCAS::detail::make_node<LogicalNode>(head, tail, node.op());
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            branch.expression->accept(*this);
            auto expression = result;
            branch.condition->accept(*this);
            auto condition = result;
            branches.push_back({expression, condition});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            default_expr = result;
        }
        result = LMCAS::detail::make_node<PiecewiseNode>(std::move(branches), default_expr);
    }

    void visit(const SummationNode& node) override {
        node.body()->accept(*this);
        auto body = result;
        node.lower_bound()->accept(*this);
        auto lower = result;
        node.upper_bound()->accept(*this);
        auto upper = result;
        result = LMCAS::detail::make_node<SummationNode>(body, node.index_var(), lower, upper);
    }

    void visit(const ProductNode& node) override {
        node.body()->accept(*this);
        auto body = result;
        node.lower_bound()->accept(*this);
        auto lower = result;
        node.upper_bound()->accept(*this);
        auto upper = result;
        result = LMCAS::detail::make_node<ProductNode>(body, node.index_var(), lower, upper);
    }

    void visit(const TransformNode& node) override {
        node.body()->accept(*this);
        auto body = result;
        node.target()->accept(*this);
        result = LMCAS::detail::make_node<TransformNode>(
            node.transform_type(), body, node.source_var(), result);
    }

    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        auto domain = result;
        node.predicate()->accept(*this);
        auto predicate = result;
        result = LMCAS::detail::make_node<QuantifierNode>(
            node.quantifier_type(), node.bound_var(), domain, predicate);
    }

    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        auto domain = result;
        node.predicate()->accept(*this);
        auto predicate = result;
        result = LMCAS::detail::make_node<SetBuilderNode>(node.element_var(), domain, predicate);
    }

    void visit(const FiniteSetNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        elements.reserve(node.elements().size());
        for (const auto& element : node.elements()) {
            element->accept(*this);
            elements.push_back(result);
        }
        result = LMCAS::detail::make_node<FiniteSetNode>(std::move(elements));
    }

    void visit(const IntervalNode& node) override {
        node.lower()->accept(*this);
        auto lower = result;
        node.upper()->accept(*this);
        result = LMCAS::detail::make_node<IntervalNode>(
            lower, result, node.lower_closed(), node.upper_closed());
    }

    void visit(const MembershipNode& node) override {
        node.element()->accept(*this);
        auto element = result;
        node.set()->accept(*this);
        result = LMCAS::detail::make_node<MembershipNode>(element, result);
    }

    void visit(const QuantityNode& node) override {
        node.value()->accept(*this);
        result = LMCAS::detail::make_node<QuantityNode>(
            result, node.dimension(), node.scale_to_base(), node.display_unit());
    }

    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        auto expanded_r = result;
        node.imag()->accept(*this);
        auto expanded_i = result;
        result = SymbolicFactory::create_complex(expanded_r, expanded_i);
    }
    void visit(const LimitNode& node) override {
        node.body()->accept(*this);
        auto body = result;
        node.point()->accept(*this);
        result = LMCAS::detail::make_node<LimitNode>(
            body, node.variable(), result, node.direction());
    }

    void visit(const RootOfNode& node) override { result = node.clone(); }
};

} // namespace LMCAS
