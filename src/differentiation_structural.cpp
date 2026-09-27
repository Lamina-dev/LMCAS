#include "internal/visitors/differentiation_visitor.hpp"

namespace LMCAS {

void DifferentiationVisitor::visit(const MatrixNode& node) {

    if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
         const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
         MatrixNode::DenseStorage new_dense;
         for(const auto& item : dense) {
             if (item) {
                item->accept(*this);
                new_dense.push_back(result);
             } else {
                new_dense.push_back(nullptr);
             }
         }
         result = LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense);
    } else {
         const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
         MatrixNode::SparseStorage new_sparse;
         for(const auto& [idx, item] : sparse) {
             item->accept(*this);
             if (!result->is_zero()) {
                 new_sparse[idx] = result;
             }
         }
         result = LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse);
    }
}

void DifferentiationVisitor::visit(const RelationalNode&) {
    unsupported("RelationalNode");
}

void DifferentiationVisitor::visit(const LogicalNode&) {
    unsupported("LogicalNode");
}

void DifferentiationVisitor::visit(const ComplexNode& node) {
    node.real()->accept(*this);
    auto dr = result;
    node.imag()->accept(*this);
    auto di = result;
    result = SymbolicFactory::create_complex(dr, di);
}

void DifferentiationVisitor::visit(const PiecewiseNode& node) {
    std::vector<PiecewiseNode::Branch> new_branches;
    for (const auto& br : node.branches()) {
        PiecewiseNode::Branch new_br;
        br.expression->accept(*this);
        new_br.expression = result;
        new_br.condition = br.condition;
        new_branches.push_back(new_br);
    }
    std::shared_ptr<const SymbolicNode> new_def = nullptr;
    if (node.default_expr()) {
        node.default_expr()->accept(*this);
        new_def = result;
    }
    result = LMCAS::detail::make_node<PiecewiseNode>(std::move(new_branches), new_def);
}

void DifferentiationVisitor::visit(const UninterpretedFunctionNode&) {
    unsupported("UninterpretedFunctionNode");
}

void DifferentiationVisitor::visit(const SummationNode&) {
    unsupported("SummationNode");
}

void DifferentiationVisitor::visit(const ProductNode&) {
    unsupported("ProductNode");
}

void DifferentiationVisitor::visit(const TransformNode&) {
    unsupported("TransformNode");
}

void DifferentiationVisitor::visit(const QuantifierNode&) {
    unsupported("QuantifierNode");
}

void DifferentiationVisitor::visit(const SetBuilderNode&) {
    unsupported("SetBuilderNode");
}
void DifferentiationVisitor::visit(const FiniteSetNode&) { unsupported("FiniteSetNode"); }
void DifferentiationVisitor::visit(const IntervalNode&) { unsupported("IntervalNode"); }
void DifferentiationVisitor::visit(const MembershipNode&) { unsupported("MembershipNode"); }
void DifferentiationVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    result = LMCAS::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit());
}

void DifferentiationVisitor::visit(const IntegralNode& node) {
    if (!node.is_definite()) {
        if (node.variable() == var) {
            result = node.body()->clone();
            return;
        }
        node.body()->accept(*this);
        result = LMCAS::detail::make_node<IntegralNode>(
            result, node.variable());
        return;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    if (node.variable() != var) {
        node.body()->accept(*this);
        if (!result->is_zero()) {
            terms.push_back(LMCAS::detail::make_node<IntegralNode>(
                result, node.variable(), node.lower(), node.upper()));
        }
    }

    auto body_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::expression_from_node(node.body()));
    auto upper_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::expression_from_node(node.upper()));
    auto lower_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::expression_from_node(node.lower()));

    node.upper()->accept(*this);
    auto upper_derivative = result;
    if (!upper_derivative->is_zero()) {
        auto at_upper = body_expression->substitute(
            node.variable(), upper_expression);
        terms.push_back(SymbolicFactory::create_multiply(
            {LMCAS::detail::node(at_upper), upper_derivative}));
    }

    node.lower()->accept(*this);
    auto lower_derivative = result;
    if (!lower_derivative->is_zero()) {
        auto at_lower = body_expression->substitute(
            node.variable(), lower_expression);
        terms.push_back(SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)),
            LMCAS::detail::node(at_lower),
            lower_derivative}));
    }
    result = SymbolicFactory::create_add(std::move(terms));
}
void DifferentiationVisitor::visit(const LimitNode&) { unsupported("LimitNode"); }
void DifferentiationVisitor::visit(const RootOfNode&) { unsupported("RootOfNode"); }

}
