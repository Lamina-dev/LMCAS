#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/squared_norm.hpp"

namespace LMCAS {
namespace {
bool approximate_norm_operand(const NumberNode* number) {
    return number && std::holds_alternative<lmmc_real_t>(number->value());
}

bool numeric_norm_operands(const NumberNode* first, const NumberNode* second,
                           const std::shared_ptr<const SymbolicNode>& second_node) {
    return first && (!second_node || second);
}

std::shared_ptr<const NumberNode> approximate_magnitude(const NumberNode& first,
                                                       const NumberNode* second) {
    const double a = normalization_number_value(first);
    const double value = second ? std::hypot(a, normalization_number_value(*second)) : std::abs(a);
    if (!std::isfinite(value)) {
        return nullptr;
    }
    return detail::make_node<NumberNode>(value);
}

bool retain_first_square(const std::shared_ptr<const SymbolicNode>& first,
                         const std::shared_ptr<const SymbolicNode>& second) {
    if (!first->is_zero()) {
        return true;
    }
    if (!second) {
        return true;
    }
    return second->is_zero();
}

bool complex_norm_operand(const std::shared_ptr<const SymbolicNode>& first,
                          const std::shared_ptr<const SymbolicNode>& second) {
    return dynamic_cast<const ComplexNode*>(first.get()) ||
           dynamic_cast<const ComplexNode*>(second.get());
}

std::shared_ptr<const SymbolicNode> normalized_norm_argument(
    const std::shared_ptr<const SymbolicNode>& first,
    const std::shared_ptr<const SymbolicNode>& second,
    const std::array<const PowerNode*, 2>& squares, detail::RewriteBudget* budget) {
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    std::size_t nodes = 0;
    terms.reserve(second ? 2 : 1);
    if (retain_first_square(first, second)) {
        normalization_check_children(budget, 1, first, squares[0]->exponent());
        normalization_append(budget, nodes, terms,
                             detail::make_node<PowerNode>(first, squares[0]->exponent()));
    }
    if (second && !second->is_zero()) {
        normalization_check_children(budget, 1, second, squares[1]->exponent());
        normalization_append(budget, nodes, terms,
                             detail::make_node<PowerNode>(second, squares[1]->exponent()));
    }
    if (terms.size() == 1) {
        return terms.front();
    }
    normalization_check_arithmetic<AddNode>(budget, terms, nodes);
    return detail::make_node<AddNode>(std::move(terms));
}
}

bool NormalizationVisitor::try_normalize_squared_norm(
    const SymbolicNode& node, std::shared_ptr<const SymbolicNode>& argument) {
    const auto squares = detail::squared_norm_terms(node);
    if (!squares[0]) {
        return false;
    }
    squares[0]->base()->accept(*this);
    auto first = result;
    std::shared_ptr<const SymbolicNode> second;
    if (squares[1]) {
        squares[1]->base()->accept(*this);
        second = result;
    }
    const auto* a = dynamic_cast<const NumberNode*>(first.get());
    const auto* b = dynamic_cast<const NumberNode*>(second.get());
    const bool numeric = numeric_norm_operands(a, b, second);
    const bool approximate = numeric && (approximate_norm_operand(a) || approximate_norm_operand(b));
    if (approximate) {
        if (auto magnitude = approximate_magnitude(*a, b)) {
            set_result(std::move(magnitude));
            return true;
        }
    }
    argument = normalized_norm_argument(first, second, squares, rewrite_budget());
    return finish_squared_norm(node, argument, (numeric && !approximate) || complex_norm_operand(first, second));
}

bool NormalizationVisitor::finish_squared_norm(
    const SymbolicNode& node, std::shared_ptr<const SymbolicNode>& argument, bool normalize_argument) {
    if (normalize_argument) {
        argument->accept(*this);
        argument = result;
        return false;
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(&node)) {
        normalization_check_children(rewrite_budget(), 1, argument, power->exponent());
        set_result(detail::make_node<PowerNode>(argument, power->exponent()));
    } else {
        normalization_check_children(rewrite_budget(), 1, argument);
        auto norm = detail::make_node<FunctionNode>(FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{argument});
        auto simplified = try_assumption_simplify(norm);
        set_result(simplified ? simplified : norm);
    }
    return true;
}
}
