#include "internal/visitors/limit_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/visitors/differentiation_visitor.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <utility>

namespace LMCAS {

int& LimitVisitor::active_limit_depth() {
    static thread_local int depth = 0;
    return depth;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::make_product_or_one(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors) {
    if (factors.empty()) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    }
    if (factors.size() == 1) {
        return factors[0];
    }
    return LMCAS::detail::make_node<MultiplyNode>(factors);
}

LimitVisitor::EvaluationScope::EvaluationScope() {
    if (LimitVisitor::active_limit_depth() < LimitVisitor::max_active_limit_depth) {
        ++LimitVisitor::active_limit_depth();
        entered = true;
    }
}

LimitVisitor::EvaluationScope::~EvaluationScope() {
    if (entered) --LimitVisitor::active_limit_depth();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::get_result() const { return result; }

void LimitVisitor::visit(const NumberNode& node) { result = node.clone(); }

void LimitVisitor::visit(const VariableNode& node) {
    result = node.name() == var ? point->clone() : node.clone();
}

void LimitVisitor::visit(const MatrixNode& node) {
    if (const auto* dense =
            std::get_if<MatrixNode::DenseStorage>(&node.storage())) {
        MatrixNode::DenseStorage entries;
        entries.reserve(dense->size());
        for (const auto& entry : *dense) {
            entry->accept(*this);
            entries.push_back(result);
        }
        result = LMCAS::detail::make_node<MatrixNode>(
            node.rows(), node.cols(), std::move(entries));
        return;
    }
    MatrixNode::SparseStorage entries;
    for (const auto& [index, entry] :
         std::get<MatrixNode::SparseStorage>(node.storage())) {
        entry->accept(*this);
        entries.emplace(index, result);
    }
    result = LMCAS::detail::make_node<MatrixNode>(
        node.rows(), node.cols(), std::move(entries));
}

void LimitVisitor::visit(const RelationalNode& node) {
    node.left()->accept(*this);
    auto left = result;
    node.right()->accept(*this);
    result = LMCAS::detail::make_node<RelationalNode>(left, result, node.op());
}

void LimitVisitor::visit(const LogicalNode& node) {
    node.left()->accept(*this);
    auto left = result;
    std::shared_ptr<const SymbolicNode> right;
    if (node.right()) {
        node.right()->accept(*this);
        right = result;
    }
    result = LMCAS::detail::make_node<LogicalNode>(left, right, node.op());
}

void LimitVisitor::visit(const SummationNode& node) {
    node.lower_bound()->accept(*this);
    auto lower = result;
    node.upper_bound()->accept(*this);
    auto upper = result;
    auto body = node.body();
    if (node.index_var() != var) {
        body->accept(*this);
        body = result;
    }
    result = LMCAS::detail::make_node<SummationNode>(
        body, node.index_var(), lower, upper);
}
void LimitVisitor::visit(const ProductNode& node) {
    node.lower_bound()->accept(*this);
    auto lower = result;
    node.upper_bound()->accept(*this);
    auto upper = result;
    auto body = node.body();
    if (node.index_var() != var) {
        body->accept(*this);
        body = result;
    }
    result = LMCAS::detail::make_node<ProductNode>(
        body, node.index_var(), lower, upper);
}
void LimitVisitor::visit(const TransformNode& node) {
    node.target()->accept(*this);
    auto target = result;
    auto body = node.body();
    if (node.source_var() != var) {
        body->accept(*this);
        body = result;
    }
    result = LMCAS::detail::make_node<TransformNode>(
        node.transform_type(), body, node.source_var(), target);
}
void LimitVisitor::visit(const QuantifierNode& node) {
    node.domain()->accept(*this);
    auto domain = result;
    auto predicate = node.predicate();
    if (node.bound_var() != var) {
        predicate->accept(*this);
        predicate = result;
    }
    result = LMCAS::detail::make_node<QuantifierNode>(
        node.quantifier_type(), node.bound_var(), domain, predicate);
}
void LimitVisitor::visit(const SetBuilderNode& node) {
    node.domain()->accept(*this);
    auto domain = result;
    auto predicate = node.predicate();
    if (node.element_var() != var) {
        predicate->accept(*this);
        predicate = result;
    }
    result = LMCAS::detail::make_node<SetBuilderNode>(
        node.element_var(), domain, predicate);
}
void LimitVisitor::visit(const FiniteSetNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> elements;
    elements.reserve(node.elements().size());
    for (const auto& element : node.elements()) {
        element->accept(*this);
        elements.push_back(result);
    }
    result = LMCAS::detail::make_node<FiniteSetNode>(std::move(elements));
}
void LimitVisitor::visit(const IntervalNode& node) {
    node.lower()->accept(*this);
    auto lower = result;
    node.upper()->accept(*this);
    result = LMCAS::detail::make_node<IntervalNode>(
        lower, result, node.lower_closed(), node.upper_closed());
}
void LimitVisitor::visit(const MembershipNode& node) {
    node.element()->accept(*this);
    auto element = result;
    node.set()->accept(*this);
    result = LMCAS::detail::make_node<MembershipNode>(element, result);
}
void LimitVisitor::visit(const UninterpretedFunctionNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    arguments.reserve(node.arguments().size());
    for (const auto& argument : node.arguments()) {
        argument->accept(*this);
        arguments.push_back(result);
    }
    result = LMCAS::detail::make_node<UninterpretedFunctionNode>(
        node.name(), std::move(arguments));
}
void LimitVisitor::visit(const IntegralNode& node) {
    std::shared_ptr<const SymbolicNode> lower;
    std::shared_ptr<const SymbolicNode> upper;
    if (node.lower()) {
        node.lower()->accept(*this);
        lower = result;
        node.upper()->accept(*this);
        upper = result;
    }
    auto body = node.body();
    if (node.variable() != var) {
        body->accept(*this);
        body = result;
    }
    result = LMCAS::detail::make_node<IntegralNode>(
        body, node.variable(), lower, upper);
}
void LimitVisitor::visit(const LimitNode& node) {
    node.point()->accept(*this);
    auto nested_point = result;
    auto body = node.body();
    if (node.variable() != var) {
        body->accept(*this);
        body = result;
    }
    result = LMCAS::detail::make_node<LimitNode>(
        body, node.variable(), nested_point, node.direction());
}
void LimitVisitor::visit(const RootOfNode& node) {
    result = node.clone();
}
void LimitVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    result = LMCAS::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit());
}

void LimitVisitor::visit(const ComplexNode& node) {
    node.real()->accept(*this);
    auto real = result;
    node.imag()->accept(*this);
    result = SymbolicFactory::create_complex(real, result);
}

} // namespace LMCAS
