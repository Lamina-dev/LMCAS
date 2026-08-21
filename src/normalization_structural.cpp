#include "visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

void NormalizationVisitor::visit(const MatrixNode& node) {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
             auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
             MatrixNode::DenseStorage new_dense;
             for(auto& item : dense) {
                 if(item) {
                    item->accept(*this);
                    new_dense.push_back(result);
                 } else {
                    new_dense.push_back(nullptr);
                 }
             }
             result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense);
        } else {
             auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
             MatrixNode::SparseStorage new_sparse;
             for(auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     new_sparse[idx] = result;
                 }
             }
             result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse);
        }
    }
void NormalizationVisitor::visit(const RelationalNode& node) {

        std::shared_ptr<const SymbolicNode> new_left = nullptr;
        std::shared_ptr<const SymbolicNode> new_right = nullptr;

        if (node.left()) {
            node.left()->accept(*this);
            new_left = result;
        }
        if (node.right()) {

            node.right()->accept(*this);
            new_right = result;
        }

        if (!new_left) new_left = node.left();
        if (!new_right) new_right = node.right();

        result = lamina::detail::make_node<RelationalNode>(new_left, new_right, node.op());
    }
void NormalizationVisitor::visit(const LogicalNode& node) {
        /// Implication: A ⇒ B = ¬A ∨ B
        if (node.op() == LogicalNode::Op::Implies) {
            auto not_left = lamina::detail::make_node<LogicalNode>(node.left(), nullptr, LogicalNode::Op::Not);
            auto or_node = lamina::detail::make_node<LogicalNode>(not_left, node.right(), LogicalNode::Op::Or);
            or_node->accept(*this);
            return;
        }

        /// NOT handling: De Morgan's laws and double negation
        if (node.op() == LogicalNode::Op::Not) {
            /// Normalize the operand first
            std::shared_ptr<const SymbolicNode> new_left = nullptr;
            if (node.left()) {
                node.left()->accept(*this);
                new_left = result;
            }
            if (!new_left) new_left = node.left();

            /// Double negation: ¬(¬A) = A
            if (auto inner_logical = std::dynamic_pointer_cast<const LogicalNode>(new_left)) {
                if (inner_logical->op() == LogicalNode::Op::Not) {
                    result = inner_logical->left();
                    return;
                }
                /// De Morgan's law: ¬(A∧B) = ¬A∨¬B
                if (inner_logical->op() == LogicalNode::Op::And) {
                    auto not_a = lamina::detail::make_node<LogicalNode>(inner_logical->left(), nullptr, LogicalNode::Op::Not);
                    auto not_b = lamina::detail::make_node<LogicalNode>(inner_logical->right(), nullptr, LogicalNode::Op::Not);
                    auto or_node = lamina::detail::make_node<LogicalNode>(not_a, not_b, LogicalNode::Op::Or);
                    or_node->accept(*this);
                    return;
                }
                /// De Morgan's law: ¬(A∨B) = ¬A∧¬B
                if (inner_logical->op() == LogicalNode::Op::Or) {
                    auto not_a = lamina::detail::make_node<LogicalNode>(inner_logical->left(), nullptr, LogicalNode::Op::Not);
                    auto not_b = lamina::detail::make_node<LogicalNode>(inner_logical->right(), nullptr, LogicalNode::Op::Not);
                    auto and_node = lamina::detail::make_node<LogicalNode>(not_a, not_b, LogicalNode::Op::And);
                    and_node->accept(*this);
                    return;
                }
            }

            result = lamina::detail::make_node<LogicalNode>(new_left, nullptr, LogicalNode::Op::Not);
            return;
        }

        /// And / Or: normalize operands
        std::shared_ptr<const SymbolicNode> new_left = nullptr;
        std::shared_ptr<const SymbolicNode> new_right = nullptr;

        if (node.left()) {
            node.left()->accept(*this);
            new_left = result;
        }
        if (node.right()) {
            node.right()->accept(*this);
            new_right = result;
        }

        if (!new_left) new_left = node.left();
        if (!new_right) new_right = node.right();

        result = lamina::detail::make_node<LogicalNode>(new_left, new_right, node.op());
    }
void NormalizationVisitor::visit(const PiecewiseNode& node) {
        std::vector<PiecewiseNode::Branch> new_branches;
        new_branches.reserve(node.branches().size());

        for (const auto& b : node.branches()) {
            /// Normalize expression
            b.expression->accept(*this);
            auto new_expr = result;

            /// Normalize condition
            b.condition->accept(*this);
            auto new_cond = result;

            /// Validate condition is RelationalNode or LogicalNode
            /// (keep it regardless, but this ensures normalization is applied)
            new_branches.push_back({new_expr, new_cond});
        }

        std::shared_ptr<const SymbolicNode> new_default = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            new_default = result;
        }

        result = lamina::detail::make_node<PiecewiseNode>(std::move(new_branches), new_default);
    }
void NormalizationVisitor::visit(const SummationNode& node) {
        node.body()->accept(*this);
        auto new_body = result;

        node.lower_bound()->accept(*this);
        auto new_lower = result;

        node.upper_bound()->accept(*this);
        auto new_upper = result;

        /// 当上下界均为具体整数且范围较小时，展开求和为显式和。
        auto lo_n = std::dynamic_pointer_cast<const NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<const NumberNode>(new_upper);
        if (lo_n && hi_n && std::holds_alternative<BigInt>(lo_n->value())
            && std::holds_alternative<BigInt>(hi_n->value())) {
            long long lo = (long long)std::get<BigInt>(lo_n->value()).to_int();
            long long hi = (long long)std::get<BigInt>(hi_n->value()).to_int();
            if (hi < lo) { result = lamina::detail::make_node<NumberNode>(BigInt(0)); return; }
            if (hi - lo < 1000) {
                std::vector<std::shared_ptr<const SymbolicNode>> terms;
                for (long long kk = lo; kk <= hi; ++kk) {
                    auto kval = lamina::detail::make_node<NumberNode>(BigInt((long long)kk));
                    auto term = norm_subst_index(new_body, node.index_var(), kval);
                    NormalizationVisitor inner;
                    term->accept(inner);
                    terms.push_back(inner.get_result());
                }
                if (terms.empty()) { result = lamina::detail::make_node<NumberNode>(BigInt(0)); return; }
                auto sum_node = lamina::detail::make_node<AddNode>(terms);
                sum_node->accept(*this);
                return;
            }
        }

        result = lamina::detail::make_node<SummationNode>(new_body, node.index_var(), new_lower, new_upper);
    }
void NormalizationVisitor::visit(const ProductNode& node) {
        node.body()->accept(*this);
        auto new_body = result;

        node.lower_bound()->accept(*this);
        auto new_lower = result;

        node.upper_bound()->accept(*this);
        auto new_upper = result;

        auto lo_n = std::dynamic_pointer_cast<const NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<const NumberNode>(new_upper);
        if (lo_n && hi_n && std::holds_alternative<BigInt>(lo_n->value())
            && std::holds_alternative<BigInt>(hi_n->value())) {
            long long lo = (long long)std::get<BigInt>(lo_n->value()).to_int();
            long long hi = (long long)std::get<BigInt>(hi_n->value()).to_int();
            if (hi < lo) { result = lamina::detail::make_node<NumberNode>(BigInt(1)); return; }
            if (hi - lo < 1000) {
                std::vector<std::shared_ptr<const SymbolicNode>> factors;
                for (long long kk = lo; kk <= hi; ++kk) {
                    auto kval = lamina::detail::make_node<NumberNode>(BigInt((long long)kk));
                    auto term = norm_subst_index(new_body, node.index_var(), kval);
                    NormalizationVisitor inner;
                    term->accept(inner);
                    factors.push_back(inner.get_result());
                }
                if (factors.empty()) { result = lamina::detail::make_node<NumberNode>(BigInt(1)); return; }
                auto prod_node = make_normalized_multiply_node(factors);
                prod_node->accept(*this);
                return;
            }
        }

        result = lamina::detail::make_node<ProductNode>(new_body, node.index_var(), new_lower, new_upper);
    }
void NormalizationVisitor::visit(const TransformNode& node) {
        node.body()->accept(*this);
        auto new_body = result;

        result = lamina::detail::make_node<TransformNode>(node.transform_type(), new_body, node.source_var(), node.target_var());
    }
void NormalizationVisitor::visit(const QuantifierNode& node) {
        node.domain()->accept(*this);
        auto new_domain = result;

        node.predicate()->accept(*this);
        auto new_predicate = result;

        /// Simplify ∀x∈S: true → true
        if (node.quantifier_type() == QuantifierNode::Type::ForAll) {
            if (new_predicate->is_one()) {
                result = lamina::detail::make_node<NumberNode>(BigInt(1));
                return;
            }
        }

        /// Simplify ∃x∈S: false → false
        if (node.quantifier_type() == QuantifierNode::Type::Exists) {
            if (new_predicate->is_zero()) {
                result = lamina::detail::make_node<NumberNode>(BigInt(0));
                return;
            }
        }

        result = lamina::detail::make_node<QuantifierNode>(node.quantifier_type(), node.bound_var(), new_domain, new_predicate);
    }
void NormalizationVisitor::visit(const SetBuilderNode& node) {
        node.domain()->accept(*this);
        auto new_domain = result;

        node.predicate()->accept(*this);
        auto new_predicate = result;

        result = lamina::detail::make_node<SetBuilderNode>(node.element_var(), new_domain, new_predicate);
    }

void NormalizationVisitor::visit(const FiniteSetNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> elements;
    elements.reserve(node.elements().size());
    for (const auto& element : node.elements()) {
        element->accept(*this);
        elements.push_back(result);
    }
    result = lamina::detail::make_node<FiniteSetNode>(std::move(elements));
}

void NormalizationVisitor::visit(const IntervalNode& node) {
    node.lower()->accept(*this);
    auto lower = result;
    node.upper()->accept(*this);
    result = lamina::detail::make_node<IntervalNode>(
        lower, result, node.lower_closed(), node.upper_closed());
}

void NormalizationVisitor::visit(const MembershipNode& node) {
    node.element()->accept(*this);
    auto element = result;
    node.set()->accept(*this);
    auto set = result;
    if (auto finite = std::dynamic_pointer_cast<const FiniteSetNode>(set)) {
        result = lamina::detail::make_node<NumberNode>(
            BigInt(finite->contains(*element) ? 1 : 0));
        return;
    }
    result = lamina::detail::make_node<MembershipNode>(element, set);
}

void NormalizationVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    result = lamina::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit());
}
