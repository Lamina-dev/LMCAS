/**
 * @file expand_visitor.hpp
 * @brief 展开访问器，将乘积和幂展开为多项式形式。
 */
#pragma once

#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {
class ExpandVisitor : public LMCAS::detail::SymbolicVisitor {
public:
    std::shared_ptr<const SymbolicNode> result;  /**< 展开结果节点。 */

    explicit ExpandVisitor(detail::RewriteBudget* budget = nullptr)
        : SymbolicVisitor(budget), result(nullptr),
          local_context_(budget ? std::nullopt
              : std::optional<ComputationContext>(std::in_place)),
          context_(budget ? budget->context() : *local_context_) {}
    explicit ExpandVisitor(ComputationContext& context,
                           detail::RewriteBudget* budget = nullptr)
        : SymbolicVisitor(budget), result(nullptr), context_(context) {}
    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    void visit(const NumberNode& node) override {
        check_child_count(0);
        set_result(node.clone(), 1);
    }
    void visit(const VariableNode& node) override {
        check_child_count(0);
        set_result(node.clone(), 1);
    }

    void visit(const AddNode& node) override {
        std::size_t size = 0;
        auto new_ops = expand_children(node.operands(), size);
        size = normalization_check_arithmetic<AddNode>(rewrite_budget(), new_ops, size);
        set_result(LMCAS::detail::make_node<AddNode>(new_ops), size);
    }

    void visit(const MultiplyNode& node) override {
        std::size_t size = 0;
        auto expanded_ops = expand_children(node.operands(), size);
        if (expanded_ops.empty()) {
            set_result(LMCAS::detail::make_node<NumberNode>(BigInt(1)), 1);
            return;
        }

        NormalizationVisitor norm(context_, detail::no_facts(), Domain::Real, rewrite_budget());
        set_result(expanded_ops[0]);
        for (size_t i = 1; i < expanded_ops.size(); ++i) {
            set_result(norm.expand_product(result, expanded_ops[i]));
        }
    }

    LMCAS_API void visit(const PowerNode& node) override;

    void visit(const FunctionNode& node) override {
        std::size_t size = 1;
        auto new_args = expand_children(node.arguments(), size);
        set_result(LMCAS::detail::make_node<FunctionNode>(node.type(), new_args), size);
    }
    void visit(const UninterpretedFunctionNode& node) override {
        std::size_t size = 1;
        auto arguments = expand_children(node.arguments(), size);
        set_result(LMCAS::detail::make_node<UninterpretedFunctionNode>(
            node.name(), std::move(arguments)), size);
    }

    void visit(const MatrixNode& node) override {
        std::size_t size = 1;
        if (rewrite_budget()) {
            if (const auto* dense = std::get_if<MatrixNode::DenseStorage>(&node.storage())) {
                check_child_count(dense->size());
                for (const auto& element : *dense) {
                    append_child(size, element);
                }
            } else {
                const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
                check_child_count(sparse.size());
                for (const auto& entry : sparse) {
                    append_child(size, entry.second);
                }
            }
        }
        set_result(node.clone(), size);
    }
    void visit(const IntegralNode& node) override {
        std::size_t size = 1;
        auto body = expand_child(node.body(), size);
        std::shared_ptr<const SymbolicNode> lower;
        std::shared_ptr<const SymbolicNode> upper;
        if (node.lower()) {
            lower = expand_child(node.lower(), size);
            upper = expand_child(node.upper(), size);
        }
        set_result(LMCAS::detail::make_node<IntegralNode>(
            body, node.variable(), lower, upper), size);
    }

    void visit(const RelationalNode& node) override {
        std::size_t size = 1;
        auto head = expand_child(node.left(), size);
        auto tail = expand_child(node.right(), size);
        set_result(LMCAS::detail::make_node<RelationalNode>(head, tail, node.op()), size);
    }

    void visit(const LogicalNode& node) override {
        std::size_t size = 1;
        auto head = expand_child(node.left(), size);
        std::shared_ptr<const SymbolicNode> tail = nullptr;
        if (node.right()) {
            tail = expand_child(node.right(), size);
        }
        set_result(LMCAS::detail::make_node<LogicalNode>(head, tail, node.op()), size);
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> branches;
        if (rewrite_budget()) {
            auto minimum = rewrite_budget()->append_size(1, node.branches().size());
            rewrite_budget()->append_size(minimum, node.branches().size());
        }
        branches.reserve(node.branches().size());
        std::size_t size = 1;
        for (const auto& branch : node.branches()) {
            auto expression = expand_child(branch.expression, size);
            auto condition = expand_child(branch.condition, size);
            branches.push_back({expression, condition});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) {
            default_expr = expand_child(node.default_expr(), size);
        }
        set_result(LMCAS::detail::make_node<PiecewiseNode>(
            std::move(branches), default_expr), size);
    }

    void visit(const SummationNode& node) override {
        std::size_t size = 1;
        auto body = expand_child(node.body(), size);
        auto lower = expand_child(node.lower_bound(), size);
        auto upper = expand_child(node.upper_bound(), size);
        set_result(LMCAS::detail::make_node<SummationNode>(
            body, node.index_var(), lower, upper), size);
    }

    void visit(const ProductNode& node) override {
        std::size_t size = 1;
        auto body = expand_child(node.body(), size);
        auto lower = expand_child(node.lower_bound(), size);
        auto upper = expand_child(node.upper_bound(), size);
        set_result(LMCAS::detail::make_node<ProductNode>(
            body, node.index_var(), lower, upper), size);
    }

    void visit(const TransformNode& node) override {
        std::size_t size = 1;
        auto body = expand_child(node.body(), size);
        auto target = expand_child(node.target(), size);
        set_result(LMCAS::detail::make_node<TransformNode>(
            node.transform_type(), body, node.source_var(), target), size);
    }

    void visit(const QuantifierNode& node) override {
        std::size_t size = 1;
        auto domain = expand_child(node.domain(), size);
        auto predicate = expand_child(node.predicate(), size);
        set_result(LMCAS::detail::make_node<QuantifierNode>(
            node.quantifier_type(), node.bound_var(), domain, predicate), size);
    }

    void visit(const SetBuilderNode& node) override {
        std::size_t size = 1;
        auto domain = expand_child(node.domain(), size);
        auto predicate = expand_child(node.predicate(), size);
        set_result(LMCAS::detail::make_node<SetBuilderNode>(
            node.element_var(), domain, predicate), size);
    }

    void visit(const FiniteSetNode& node) override {
        std::size_t size = 1;
        auto elements = expand_children(node.elements(), size);
        set_result(LMCAS::detail::make_node<FiniteSetNode>(std::move(elements)), size);
    }

    void visit(const IntervalNode& node) override {
        std::size_t size = 1;
        auto lower = expand_child(node.lower(), size);
        auto upper = expand_child(node.upper(), size);
        set_result(LMCAS::detail::make_node<IntervalNode>(
            lower, upper, node.lower_closed(), node.upper_closed()), size);
    }

    void visit(const MembershipNode& node) override {
        std::size_t size = 1;
        auto element = expand_child(node.element(), size);
        auto set = expand_child(node.set(), size);
        set_result(LMCAS::detail::make_node<MembershipNode>(element, set), size);
    }

    void visit(const QuantityNode& node) override {
        std::size_t size = 1;
        auto value = expand_child(node.value(), size);
        set_result(LMCAS::detail::make_node<QuantityNode>(
            value, node.dimension(), node.scale_to_base(), node.display_unit()), size);
    }

    void visit(const ComplexNode& node) override {
        std::size_t size = 1;
        auto expanded_r = expand_child(node.real(), size);
        auto expanded_i = expand_child(node.imag(), size);
        set_result(SymbolicFactory::create_complex(expanded_r, expanded_i), size);
    }
    void visit(const LimitNode& node) override {
        std::size_t size = 1;
        auto body = expand_child(node.body(), size);
        auto point = expand_child(node.point(), size);
        set_result(LMCAS::detail::make_node<LimitNode>(
            body, node.variable(), point, node.direction()), size);
    }

    void visit(const RootOfNode& node) override {
        check_child_count(0);
        set_result(node.clone(), 1);
    }

private:
    std::optional<ComputationContext> local_context_;
    ComputationContext& context_;

    void check_child_count(std::size_t count, std::size_t roots = 1) const {
        if (auto* budget = rewrite_budget()) {
            budget->append_size(roots, count);
        }
    }

    void append_child(std::size_t& size,
                      const std::shared_ptr<const SymbolicNode>& child) const {
        if (auto* budget = rewrite_budget()) {
            size = budget->append_size(size, budget->measure(child));
        }
    }

    std::shared_ptr<const SymbolicNode> expand_child(
        const std::shared_ptr<const SymbolicNode>& child, std::size_t& size) {
        child->accept(*this);
        append_child(size, result);
        return result;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> expand_children(
        const std::vector<std::shared_ptr<const SymbolicNode>>& children, std::size_t& size) {
        check_child_count(children.size(), size);
        std::vector<std::shared_ptr<const SymbolicNode>> expanded;
        expanded.reserve(children.size());
        for (const auto& child : children) {
            auto expanded_child = expand_child(child, size);
            expanded.push_back(std::move(expanded_child));
        }
        return expanded;
    }

    void set_result(std::shared_ptr<const SymbolicNode> candidate, std::size_t size) {
        if (auto* budget = rewrite_budget()) {
            budget->require_nodes(size);
        }
        result = std::move(candidate);
    }

    void set_result(std::shared_ptr<const SymbolicNode> candidate) {
        std::size_t size = 0;
        append_child(size, candidate);
        set_result(std::move(candidate), size);
    }
};

}
