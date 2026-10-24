#include "internal/visitors/differentiation_visitor.hpp"
#include "internal/expression_transform.hpp"

namespace LMCAS {

void DifferentiationVisitor::visit(const MatrixNode& node) {
    std::size_t size = 1;

    if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
         const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
         MatrixNode::DenseStorage new_dense;
         rewrite_budget()->append_size(size, dense.size());
         new_dense.reserve(dense.size());
         for(const auto& item : dense) {
             if (item) {
                item->accept(*this);
                append_result(new_dense, size, result);
             } else {
                new_dense.push_back(nullptr);
             }
         }
         set_result(LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense));
    } else {
         const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
         MatrixNode::SparseStorage new_sparse;
         for(const auto& [idx, item] : sparse) {
             item->accept(*this);
             if (!result->is_zero()) {
                 size = rewrite_budget()->append_size(size, rewrite_budget()->measure(result));
                 new_sparse[idx] = result;
             }
         }
         set_result(LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse));
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
    set_result(SymbolicFactory::create_complex(dr, di));
}

void DifferentiationVisitor::visit(const PiecewiseNode& node) {
    std::vector<PiecewiseNode::Branch> new_branches;
    std::size_t size = 1;
    auto minimum = rewrite_budget()->append_size(size, node.branches().size());
    rewrite_budget()->append_size(minimum, node.branches().size());
    new_branches.reserve(node.branches().size());
    for (const auto& br : node.branches()) {
        PiecewiseNode::Branch new_br;
        br.expression->accept(*this);
        new_br.expression = result;
        new_br.condition = br.condition;
        size = rewrite_budget()->append_size(size, rewrite_budget()->measure(new_br.expression));
        size = rewrite_budget()->append_size(size, rewrite_budget()->measure(new_br.condition));
        new_branches.push_back(new_br);
    }
    std::shared_ptr<const SymbolicNode> new_def = nullptr;
    if (node.default_expr()) {
        node.default_expr()->accept(*this);
        new_def = result;
        rewrite_budget()->append_size(size, rewrite_budget()->measure(new_def));
    }
    set_result(LMCAS::detail::make_node<PiecewiseNode>(std::move(new_branches), new_def));
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
    set_result(LMCAS::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit()));
}

void DifferentiationVisitor::visit(const IntegralNode& node) {
    if (!node.is_definite()) {
        if (node.variable() == var) {
            set_result(node.body());
            return;
        }
        node.body()->accept(*this);
        set_result(LMCAS::detail::make_node<IntegralNode>(
            result, node.variable()));
        return;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    std::size_t size = 1;
    if (node.variable() != var) {
        node.body()->accept(*this);
        if (!result->is_zero()) {
            append_result(terms, size, LMCAS::detail::make_node<IntegralNode>(
                result, node.variable(), node.lower(), node.upper()));
        }
    }

    node.upper()->accept(*this);
    auto upper_derivative = result;
    if (!upper_derivative->is_zero()) {
        auto at_upper = detail::substitute_node(
            node.body(), node.variable(), node.upper(), *rewrite_budget());
        append_result(terms, size, SymbolicFactory::create_multiply(
            {at_upper, upper_derivative}, rewrite_budget()));
    }

    node.lower()->accept(*this);
    auto lower_derivative = result;
    if (!lower_derivative->is_zero()) {
        auto at_lower = detail::substitute_node(
            node.body(), node.variable(), node.lower(), *rewrite_budget());
        append_result(terms, size, SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)), at_lower,
            lower_derivative}, rewrite_budget()));
    }
    set_result(SymbolicFactory::create_add(std::move(terms), rewrite_budget()));
}
void DifferentiationVisitor::visit(const LimitNode&) { unsupported("LimitNode"); }
void DifferentiationVisitor::visit(const RootOfNode&) { unsupported("RootOfNode"); }

}
