#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast/traversal.hpp"

#include <cstdint>

namespace LMCAS {

std::shared_ptr<const SymbolicNode> NormalizationVisitor::normalize_bound_body(
    const SymbolicNode& node) {
    const auto binder = detail::binder_view(node);
    if (&facts_ == &detail::no_facts()) {
        binder->scoped_body->accept(*this);
        return result;
    }
    const detail::ScopedFacts scoped(facts_, binder->bound_name);
    NormalizationVisitor inner(context_, scoped, domain_, rewrite_budget());
    binder->scoped_body->accept(inner);
    return inner.get_result();
}

void NormalizationVisitor::visit(const MatrixNode& node) {
        std::size_t nodes = 1;
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
             auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
             MatrixNode::DenseStorage new_dense;
             for(auto& item : dense) {
                 if(item) {
                    item->accept(*this);
                    normalization_append(rewrite_budget(), nodes, new_dense, result);
                 } else {
                    normalization_append(rewrite_budget(), nodes, new_dense, nullptr);
                 }
             }
             set_result(LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense));
        } else {
             auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
             MatrixNode::SparseStorage new_sparse;
             for(auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     if (rewrite_budget()) {
                         nodes = rewrite_budget()->append_size(nodes, rewrite_budget()->measure(result));
                     }
                     new_sparse[idx] = result;
                 }
             }
             set_result(LMCAS::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse));
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

        normalization_check_children(rewrite_budget(), 1, new_left, new_right);
        set_result(LMCAS::detail::make_node<RelationalNode>(new_left, new_right, node.op()));
    }
void NormalizationVisitor::visit(const LogicalNode& node) {
        /// Implication: A ⇒ B = ¬A ∨ B
        if (node.op() == LogicalNode::Op::Implies) {
            normalization_check_children(rewrite_budget(), 1, node.left());
            auto not_left = LMCAS::detail::make_node<LogicalNode>(node.left(), nullptr, LogicalNode::Op::Not);
            normalization_check_children(rewrite_budget(), 1, not_left, node.right());
            auto or_node = LMCAS::detail::make_node<LogicalNode>(not_left, node.right(), LogicalNode::Op::Or);
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
                    set_result(inner_logical->left());
                    return;
                }
                /// De Morgan's law: ¬(A∧B) = ¬A∨¬B
                if (inner_logical->op() == LogicalNode::Op::And) {
                    normalization_check_children(rewrite_budget(), 1, inner_logical->left());
                    auto not_a = LMCAS::detail::make_node<LogicalNode>(inner_logical->left(), nullptr, LogicalNode::Op::Not);
                    normalization_check_children(rewrite_budget(), 1, inner_logical->right());
                    auto not_b = LMCAS::detail::make_node<LogicalNode>(inner_logical->right(), nullptr, LogicalNode::Op::Not);
                    normalization_check_children(rewrite_budget(), 1, not_a, not_b);
                    auto or_node = LMCAS::detail::make_node<LogicalNode>(not_a, not_b, LogicalNode::Op::Or);
                    or_node->accept(*this);
                    return;
                }
                /// De Morgan's law: ¬(A∨B) = ¬A∧¬B
                if (inner_logical->op() == LogicalNode::Op::Or) {
                    normalization_check_children(rewrite_budget(), 1, inner_logical->left());
                    auto not_a = LMCAS::detail::make_node<LogicalNode>(inner_logical->left(), nullptr, LogicalNode::Op::Not);
                    normalization_check_children(rewrite_budget(), 1, inner_logical->right());
                    auto not_b = LMCAS::detail::make_node<LogicalNode>(inner_logical->right(), nullptr, LogicalNode::Op::Not);
                    normalization_check_children(rewrite_budget(), 1, not_a, not_b);
                    auto and_node = LMCAS::detail::make_node<LogicalNode>(not_a, not_b, LogicalNode::Op::And);
                    and_node->accept(*this);
                    return;
                }
            }

            normalization_check_children(rewrite_budget(), 1, new_left);
            set_result(LMCAS::detail::make_node<LogicalNode>(new_left, nullptr, LogicalNode::Op::Not));
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

        normalization_check_children(rewrite_budget(), 1, new_left, new_right);
        set_result(LMCAS::detail::make_node<LogicalNode>(new_left, new_right, node.op()));
    }
void NormalizationVisitor::visit(const PiecewiseNode& node) {
        std::vector<PiecewiseNode::Branch> new_branches;
        std::size_t nodes = 1;
        normalization_check_count(rewrite_budget(), node.branches().size(), 2);
        new_branches.reserve(node.branches().size());

        for (const auto& b : node.branches()) {
            /// Normalize expression
            b.expression->accept(*this);
            auto new_expr = result;

            /// Normalize condition
            b.condition->accept(*this);
            auto new_cond = result;

            if (rewrite_budget()) {
                nodes = rewrite_budget()->append_size(nodes, rewrite_budget()->measure(new_expr));
                nodes = rewrite_budget()->append_size(nodes, rewrite_budget()->measure(new_cond));
            }
            new_branches.push_back({new_expr, new_cond});
        }

        std::shared_ptr<const SymbolicNode> new_default = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            new_default = result;
            if (rewrite_budget()) {
                nodes = rewrite_budget()->append_size(nodes, rewrite_budget()->measure(new_default));
            }
        }

        set_result(LMCAS::detail::make_node<PiecewiseNode>(std::move(new_branches), new_default));
    }
void NormalizationVisitor::visit(const SummationNode& node) {
        auto new_body = normalize_bound_body(node);

        node.lower_bound()->accept(*this);
        auto new_lower = result;

        node.upper_bound()->accept(*this);
        auto new_upper = result;

        /// 当上下界均为具体整数且范围较小时，展开求和为显式和。
        auto lo_n = std::dynamic_pointer_cast<const NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<const NumberNode>(new_upper);
        BigInt lo;
        BigInt hi;
        if (try_get_integer_value(lo_n, lo) && try_get_integer_value(hi_n, hi)) {
            if (hi < lo) {
                set_result(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
                return;
            }
            const BigInt span = hi - lo;
            if (span < BigInt(1000)) {
                const auto count = span.try_to_uint64();
                std::vector<std::shared_ptr<const SymbolicNode>> terms;
                std::size_t nodes = 0;
                normalization_check_count(rewrite_budget(), static_cast<std::size_t>(*count) + 1, 1, false);
                terms.reserve(static_cast<std::size_t>(*count) + 1);
                BigInt index = lo;
                for (std::uint64_t offset = 0; offset <= *count; ++offset, index += BigInt(1)) {
                    auto kval = LMCAS::detail::make_node<NumberNode>(index);
                    auto term = norm_subst_index(new_body, node.index_var(), kval, rewrite_budget());
                    normalization_check_children(rewrite_budget(), 0, term);
                    NormalizationVisitor inner(context_, facts_, domain_, rewrite_budget());
                    term->accept(inner);
                    normalization_append(rewrite_budget(), nodes, terms, inner.get_result());
                }
                if (terms.empty()) { set_result(LMCAS::detail::make_node<NumberNode>(BigInt(0))); return; }
                normalization_check_arithmetic<AddNode>(rewrite_budget(), terms, nodes);
                auto sum_node = LMCAS::detail::make_node<AddNode>(terms);
                sum_node->accept(*this);
                return;
            }
        }

        normalization_check_children(rewrite_budget(), 1, new_body, new_lower, new_upper);
        set_result(LMCAS::detail::make_node<SummationNode>(new_body, node.index_var(), new_lower, new_upper));
    }
void NormalizationVisitor::visit(const ProductNode& node) {
        auto new_body = normalize_bound_body(node);

        node.lower_bound()->accept(*this);
        auto new_lower = result;

        node.upper_bound()->accept(*this);
        auto new_upper = result;

        auto lo_n = std::dynamic_pointer_cast<const NumberNode>(new_lower);
        auto hi_n = std::dynamic_pointer_cast<const NumberNode>(new_upper);
        BigInt lo;
        BigInt hi;
        if (try_get_integer_value(lo_n, lo) && try_get_integer_value(hi_n, hi)) {
            if (hi < lo) {
                set_result(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
                return;
            }
            const BigInt span = hi - lo;
            if (span < BigInt(1000)) {
                const auto count = span.try_to_uint64();
                std::vector<std::shared_ptr<const SymbolicNode>> factors;
                std::size_t nodes = 0;
                normalization_check_count(rewrite_budget(), static_cast<std::size_t>(*count) + 1, 1, false);
                factors.reserve(static_cast<std::size_t>(*count) + 1);
                BigInt index = lo;
                for (std::uint64_t offset = 0; offset <= *count; ++offset, index += BigInt(1)) {
                    auto kval = LMCAS::detail::make_node<NumberNode>(index);
                    auto term = norm_subst_index(new_body, node.index_var(), kval, rewrite_budget());
                    normalization_check_children(rewrite_budget(), 0, term);
                    NormalizationVisitor inner(context_, facts_, domain_, rewrite_budget());
                    term->accept(inner);
                    normalization_append(rewrite_budget(), nodes, factors, inner.get_result());
                }
                if (factors.empty()) { set_result(LMCAS::detail::make_node<NumberNode>(BigInt(1))); return; }
                auto prod_node = make_normalized_multiply_node(factors, rewrite_budget());
                prod_node->accept(*this);
                return;
            }
        }

        normalization_check_children(rewrite_budget(), 1, new_body, new_lower, new_upper);
        set_result(LMCAS::detail::make_node<ProductNode>(new_body, node.index_var(), new_lower, new_upper));
    }
void NormalizationVisitor::visit(const TransformNode& node) {
        auto new_body = normalize_bound_body(node);
        node.target()->accept(*this);
        auto new_target = result;
        normalization_check_children(rewrite_budget(), 1, new_body, new_target);
        set_result(LMCAS::detail::make_node<TransformNode>(
            node.transform_type(), new_body, node.source_var(), new_target));
    }
void NormalizationVisitor::visit(const QuantifierNode& node) {
        node.domain()->accept(*this);
        auto new_domain = result;

        auto new_predicate = normalize_bound_body(node);

        /// Simplify ∀x∈S: true → true
        if (node.quantifier_type() == QuantifierNode::Type::ForAll) {
            if (new_predicate->is_one()) {
                set_result(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
                return;
            }
        }

        /// Simplify ∃x∈S: false → false
        if (node.quantifier_type() == QuantifierNode::Type::Exists) {
            if (new_predicate->is_zero()) {
                set_result(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
                return;
            }
        }

        normalization_check_children(rewrite_budget(), 1, new_domain, new_predicate);
        set_result(LMCAS::detail::make_node<QuantifierNode>(node.quantifier_type(), node.bound_var(), new_domain, new_predicate));
    }
void NormalizationVisitor::visit(const SetBuilderNode& node) {
        node.domain()->accept(*this);
        auto new_domain = result;

        auto new_predicate = normalize_bound_body(node);

        normalization_check_children(rewrite_budget(), 1, new_domain, new_predicate);
        set_result(LMCAS::detail::make_node<SetBuilderNode>(node.element_var(), new_domain, new_predicate));
    }

void NormalizationVisitor::visit(const FiniteSetNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> elements;
    std::size_t nodes = 1;
    normalization_check_count(rewrite_budget(), node.elements().size());
    elements.reserve(node.elements().size());
    for (const auto& element : node.elements()) {
        element->accept(*this);
        normalization_append(rewrite_budget(), nodes, elements, result);
    }
    set_result(LMCAS::detail::make_node<FiniteSetNode>(std::move(elements)));
}

void NormalizationVisitor::visit(const IntervalNode& node) {
    node.lower()->accept(*this);
    auto lower = result;
    node.upper()->accept(*this);
    normalization_check_children(rewrite_budget(), 1, lower, result);
    set_result(LMCAS::detail::make_node<IntervalNode>(
        lower, result, node.lower_closed(), node.upper_closed()));
}

void NormalizationVisitor::visit(const MembershipNode& node) {
    node.element()->accept(*this);
    auto element = result;
    node.set()->accept(*this);
    auto set = result;
    if (auto finite = std::dynamic_pointer_cast<const FiniteSetNode>(set)) {
        set_result(LMCAS::detail::make_node<NumberNode>(
            BigInt(finite->contains(*element) ? 1 : 0)));
        return;
    }
    normalization_check_children(rewrite_budget(), 1, element, set);
    set_result(LMCAS::detail::make_node<MembershipNode>(element, set));
}

void NormalizationVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    normalization_check_children(rewrite_budget(), 1, result);
    set_result(LMCAS::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit()));
}

void NormalizationVisitor::visit(const IntegralNode& node) {
    auto body = normalize_bound_body(node);
    std::shared_ptr<const SymbolicNode> lower;
    std::shared_ptr<const SymbolicNode> upper;
    if (node.lower()) {
        node.lower()->accept(*this);
        lower = result;
        node.upper()->accept(*this);
        upper = result;
    }
    normalization_check_children(rewrite_budget(), 1, body, lower, upper);
    set_result(LMCAS::detail::make_node<IntegralNode>(
        std::move(body), node.variable(), std::move(lower), std::move(upper)));
}

void NormalizationVisitor::visit(const LimitNode& node) {
    auto body = normalize_bound_body(node);
    node.point()->accept(*this);
    normalization_check_children(rewrite_budget(), 1, body, result);
    set_result(LMCAS::detail::make_node<LimitNode>(
        body, node.variable(), result, node.direction()));
}

void NormalizationVisitor::visit(const RootOfNode& node) {
    if (rewrite_budget()) { rewrite_budget()->require_nodes(1); }
    set_result(node.clone());
}

} // namespace LMCAS
