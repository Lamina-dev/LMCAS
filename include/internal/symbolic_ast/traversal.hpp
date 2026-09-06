/** @file internal/symbolic_ast/traversal.hpp */
#pragma once

#include "special_forms.hpp"

#include <optional>
#include <string_view>
#include <utility>

namespace LMCAS::detail {

/** Scope information for an AST node that introduces one lexical binder. */
struct BinderView {
    std::string_view bound_name;
    SymbolicNodePtr scoped_body;
    std::vector<SymbolicNodePtr> outside_scope;
};

inline std::optional<BinderView> binder_view(const SymbolicNode& node) {
    if (const auto* sum = dynamic_cast<const SummationNode*>(&node)) {
        return BinderView{sum->index_var(), sum->body(),
                          {sum->lower_bound(), sum->upper_bound()}};
    }
    if (const auto* product = dynamic_cast<const ProductNode*>(&node)) {
        return BinderView{product->index_var(), product->body(),
                          {product->lower_bound(), product->upper_bound()}};
    }
    if (const auto* integral = dynamic_cast<const IntegralNode*>(&node)) {
        std::vector<SymbolicNodePtr> bounds;
        if (integral->lower()) bounds.push_back(integral->lower());
        if (integral->upper()) bounds.push_back(integral->upper());
        return BinderView{integral->variable(), integral->body(), std::move(bounds)};
    }
    if (const auto* transform = dynamic_cast<const TransformNode*>(&node)) {
        return BinderView{transform->source_var(), transform->body(),
                          {transform->target()}};
    }
    if (const auto* quantifier = dynamic_cast<const QuantifierNode*>(&node)) {
        return BinderView{quantifier->bound_var(), quantifier->predicate(),
                          {quantifier->domain()}};
    }
    if (const auto* set_builder = dynamic_cast<const SetBuilderNode*>(&node)) {
        return BinderView{set_builder->element_var(), set_builder->predicate(),
                          {set_builder->domain()}};
    }
    if (const auto* limit = dynamic_cast<const LimitNode*>(&node)) {
        return BinderView{limit->variable(), limit->body(), {limit->point()}};
    }
    return std::nullopt;
}

/** Exhaustive depth-first traversal with correct coverage for every AST child. */
class RecursiveSymbolicVisitor : public SymbolicVisitor {
protected:
    void visit_child(const SymbolicNodePtr& child) {
        if (child) child->accept(*this);
    }

    template <class Range>
    void visit_children(const Range& children) {
        for (const auto& child : children) visit_child(child);
    }

public:
    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}
    void visit(const AddNode& node) override { visit_children(node.operands()); }
    void visit(const MultiplyNode& node) override { visit_children(node.operands()); }
    void visit(const PowerNode& node) override {
        visit_child(node.base());
        visit_child(node.exponent());
    }
    void visit(const FunctionNode& node) override { visit_children(node.arguments()); }
    void visit(const UninterpretedFunctionNode& node) override {
        visit_children(node.arguments());
    }
    void visit(const MatrixNode& node) override {
        if (const auto* dense = std::get_if<MatrixNode::DenseStorage>(&node.storage())) {
            visit_children(*dense);
            return;
        }
        for (const auto& [index, value] :
             std::get<MatrixNode::SparseStorage>(node.storage())) {
            (void)index;
            visit_child(value);
        }
    }
    void visit(const RelationalNode& node) override {
        visit_child(node.left());
        visit_child(node.right());
    }
    void visit(const LogicalNode& node) override {
        visit_child(node.left());
        visit_child(node.right());
    }
    void visit(const PiecewiseNode& node) override {
        for (const auto& branch : node.branches()) {
            visit_child(branch.expression);
            visit_child(branch.condition);
        }
        visit_child(node.default_expr());
    }
    void visit(const SummationNode& node) override {
        visit_child(node.lower_bound());
        visit_child(node.upper_bound());
        visit_child(node.body());
    }
    void visit(const ProductNode& node) override {
        visit_child(node.lower_bound());
        visit_child(node.upper_bound());
        visit_child(node.body());
    }
    void visit(const TransformNode& node) override {
        visit_child(node.target());
        visit_child(node.body());
    }
    void visit(const QuantifierNode& node) override {
        visit_child(node.domain());
        visit_child(node.predicate());
    }
    void visit(const SetBuilderNode& node) override {
        visit_child(node.domain());
        visit_child(node.predicate());
    }
    void visit(const FiniteSetNode& node) override { visit_children(node.elements()); }
    void visit(const IntervalNode& node) override {
        visit_child(node.lower());
        visit_child(node.upper());
    }
    void visit(const MembershipNode& node) override {
        visit_child(node.element());
        visit_child(node.set());
    }
    void visit(const QuantityNode& node) override { visit_child(node.value()); }
    void visit(const ComplexNode& node) override {
        visit_child(node.real());
        visit_child(node.imag());
    }
    void visit(const IntegralNode& node) override {
        visit_child(node.lower());
        visit_child(node.upper());
        visit_child(node.body());
    }
    void visit(const LimitNode& node) override {
        visit_child(node.point());
        visit_child(node.body());
    }
    void visit(const RootOfNode&) override {}
};

template <class Node>
class NodeTypeDetector final : public RecursiveSymbolicVisitor {
public:
    bool found() const noexcept { return found_; }
    void visit(const Node&) override { found_ = true; }

private:
    bool found_ = false;
};

template <class Node>
bool contains_node_type(const SymbolicNodePtr& root) {
    if (!root) return false;
    NodeTypeDetector<Node> detector;
    root->accept(detector);
    return detector.found();
}

/** Immutable, path-copying AST rewrite. Unchanged subtrees retain pointer identity. */
class SymbolicRewriter : public SymbolicVisitor {
public:
    SymbolicNodePtr rewrite(const SymbolicNodePtr& node) {
        if (!node) return nullptr;
        const auto previous_current = current_;
        const auto previous_result = result_;
        current_ = node;
        result_.reset();
        node->accept(*this);
        auto rewritten = result_;
        current_ = previous_current;
        result_ = previous_result;
        if (!rewritten) throw std::logic_error("SymbolicRewriter visit did not set a result");
        return rewritten;
    }

protected:
    SymbolicNodePtr current() const noexcept { return current_; }
    void set_result(SymbolicNodePtr result) { result_ = std::move(result); }

    SymbolicNodePtr rewrite_child(const SymbolicNodePtr& child, bool& changed) {
        auto rewritten = rewrite(child);
        changed = changed || rewritten != child;
        return rewritten;
    }

    template <class Range>
    std::vector<SymbolicNodePtr> rewrite_children(const Range& children, bool& changed) {
        std::vector<SymbolicNodePtr> rewritten;
        rewritten.reserve(children.size());
        for (const auto& child : children) {
            rewritten.push_back(rewrite_child(child, changed));
        }
        return rewritten;
    }

public:
    void visit(const NumberNode&) override { set_result(current()); }
    void visit(const VariableNode&) override { set_result(current()); }
    void visit(const AddNode& node) override {
        bool changed = false;
        auto operands = rewrite_children(node.operands(), changed);
        set_result(changed ? SymbolicFactory::create_add(std::move(operands)) : current());
    }
    void visit(const MultiplyNode& node) override {
        bool changed = false;
        auto operands = rewrite_children(node.operands(), changed);
        set_result(changed ? SymbolicFactory::create_multiply(std::move(operands)) : current());
    }
    void visit(const PowerNode& node) override {
        bool changed = false;
        auto base = rewrite_child(node.base(), changed);
        auto exponent = rewrite_child(node.exponent(), changed);
        set_result(changed ? SymbolicFactory::create_power(base, exponent) : current());
    }
    void visit(const FunctionNode& node) override {
        bool changed = false;
        auto arguments = rewrite_children(node.arguments(), changed);
        set_result(changed ? make_node<FunctionNode>(node.type(), std::move(arguments)) : current());
    }
    void visit(const UninterpretedFunctionNode& node) override {
        bool changed = false;
        auto arguments = rewrite_children(node.arguments(), changed);
        set_result(changed ? make_node<UninterpretedFunctionNode>(node.name(), std::move(arguments))
                           : current());
    }
    void visit(const MatrixNode& node) override {
        bool changed = false;
        if (const auto* dense = std::get_if<MatrixNode::DenseStorage>(&node.storage())) {
            auto entries = rewrite_children(*dense, changed);
            set_result(changed ? make_node<MatrixNode>(node.rows(), node.cols(), std::move(entries))
                               : current());
            return;
        }
        auto entries = std::get<MatrixNode::SparseStorage>(node.storage());
        for (auto& [index, value] : entries) {
            (void)index;
            value = rewrite_child(value, changed);
        }
        set_result(changed ? make_node<MatrixNode>(node.rows(), node.cols(), std::move(entries))
                           : current());
    }
    void visit(const RelationalNode& node) override {
        bool changed = false;
        auto left = rewrite_child(node.left(), changed);
        auto right = rewrite_child(node.right(), changed);
        set_result(changed ? make_node<RelationalNode>(left, right, node.op()) : current());
    }
    void visit(const LogicalNode& node) override {
        bool changed = false;
        auto left = rewrite_child(node.left(), changed);
        auto right = node.right() ? rewrite_child(node.right(), changed) : nullptr;
        set_result(changed ? make_node<LogicalNode>(left, right, node.op()) : current());
    }
    void visit(const PiecewiseNode& node) override {
        bool changed = false;
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            branches.push_back({rewrite_child(branch.expression, changed),
                                rewrite_child(branch.condition, changed)});
        }
        auto fallback = node.default_expr()
            ? rewrite_child(node.default_expr(), changed) : nullptr;
        set_result(changed ? make_node<PiecewiseNode>(std::move(branches), fallback) : current());
    }
    void visit(const SummationNode& node) override {
        bool changed = false;
        auto body = rewrite_child(node.body(), changed);
        auto lower = rewrite_child(node.lower_bound(), changed);
        auto upper = rewrite_child(node.upper_bound(), changed);
        set_result(changed ? make_node<SummationNode>(body, node.index_var(), lower, upper)
                           : current());
    }
    void visit(const ProductNode& node) override {
        bool changed = false;
        auto body = rewrite_child(node.body(), changed);
        auto lower = rewrite_child(node.lower_bound(), changed);
        auto upper = rewrite_child(node.upper_bound(), changed);
        set_result(changed ? make_node<ProductNode>(body, node.index_var(), lower, upper)
                           : current());
    }
    void visit(const TransformNode& node) override {
        bool changed = false;
        auto body = rewrite_child(node.body(), changed);
        auto target = rewrite_child(node.target(), changed);
        set_result(changed ? make_node<TransformNode>(node.transform_type(), body,
                                                      node.source_var(), target)
                           : current());
    }
    void visit(const QuantifierNode& node) override {
        bool changed = false;
        auto domain = rewrite_child(node.domain(), changed);
        auto predicate = rewrite_child(node.predicate(), changed);
        set_result(changed ? make_node<QuantifierNode>(node.quantifier_type(), node.bound_var(),
                                                       domain, predicate)
                           : current());
    }
    void visit(const SetBuilderNode& node) override {
        bool changed = false;
        auto domain = rewrite_child(node.domain(), changed);
        auto predicate = rewrite_child(node.predicate(), changed);
        set_result(changed ? make_node<SetBuilderNode>(node.element_var(), domain, predicate)
                           : current());
    }
    void visit(const FiniteSetNode& node) override {
        bool changed = false;
        auto elements = rewrite_children(node.elements(), changed);
        set_result(changed ? make_node<FiniteSetNode>(std::move(elements)) : current());
    }
    void visit(const IntervalNode& node) override {
        bool changed = false;
        auto lower = rewrite_child(node.lower(), changed);
        auto upper = rewrite_child(node.upper(), changed);
        set_result(changed ? make_node<IntervalNode>(lower, upper,
                                                     node.lower_closed(), node.upper_closed())
                           : current());
    }
    void visit(const MembershipNode& node) override {
        bool changed = false;
        auto element = rewrite_child(node.element(), changed);
        auto set = rewrite_child(node.set(), changed);
        set_result(changed ? make_node<MembershipNode>(element, set) : current());
    }
    void visit(const QuantityNode& node) override {
        bool changed = false;
        auto value = rewrite_child(node.value(), changed);
        set_result(changed ? make_node<QuantityNode>(value, node.dimension(),
                                                     node.scale_to_base(), node.display_unit())
                           : current());
    }
    void visit(const ComplexNode& node) override {
        bool changed = false;
        auto real = rewrite_child(node.real(), changed);
        auto imag = rewrite_child(node.imag(), changed);
        set_result(changed ? SymbolicFactory::create_complex(real, imag) : current());
    }
    void visit(const IntegralNode& node) override {
        bool changed = false;
        auto body = rewrite_child(node.body(), changed);
        auto lower = node.lower() ? rewrite_child(node.lower(), changed) : nullptr;
        auto upper = node.upper() ? rewrite_child(node.upper(), changed) : nullptr;
        set_result(changed ? make_node<IntegralNode>(body, node.variable(), lower, upper)
                           : current());
    }
    void visit(const LimitNode& node) override {
        bool changed = false;
        auto body = rewrite_child(node.body(), changed);
        auto point = rewrite_child(node.point(), changed);
        set_result(changed ? make_node<LimitNode>(
                                 body, node.variable(), point, node.direction())
                           : current());
    }
    void visit(const RootOfNode&) override { set_result(current()); }

private:
    SymbolicNodePtr current_;
    SymbolicNodePtr result_;
};

} // namespace LMCAS::detail
