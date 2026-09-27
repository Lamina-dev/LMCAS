#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {
namespace {
bool sum_term_less(const std::shared_ptr<const SymbolicNode>& left,
                   const std::shared_ptr<const SymbolicNode>& right) {
    const BigInt a = get_node_degree_helper(left);
    const BigInt b = get_node_degree_helper(right);
    if (a != b) {
        return a > b;
    }
    const bool left_number = std::dynamic_pointer_cast<const NumberNode>(left) != nullptr;
    const bool right_number = std::dynamic_pointer_cast<const NumberNode>(right) != nullptr;
    if (left_number != right_number) {
        return right_number;
    }
    return left->compare(*right) < 0;
}

std::shared_ptr<const SymbolicNode> weighted_term(
    const std::shared_ptr<const SymbolicNode>& term,
    const std::shared_ptr<const NumberNode>& coefficient, detail::RewriteBudget* budget) {
    if (coefficient->is_one()) {
        return term;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    std::size_t nodes = 0;
    normalization_append(budget, nodes, operands, coefficient);
    if (const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(term)) {
        for (const auto& child : multiply->operands()) {
            normalization_append(budget, nodes, operands, child);
        }
    } else {
        normalization_append(budget, nodes, operands, term);
    }
    return make_normalized_multiply_node(operands, budget);
}
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::get_result() const {
    return result;
}

void NormalizationVisitor::set_result(std::shared_ptr<const SymbolicNode> candidate) {
    if (auto* budget = rewrite_budget()) {
        budget->measure(candidate);
    }
    result.swap(candidate);
}

void NormalizationVisitor::visit(const NumberNode& node) {
    if (rewrite_budget()) { rewrite_budget()->require_nodes(1); }
    set_result(node.clone());
}

void NormalizationVisitor::visit(const VariableNode& node) {
    if (rewrite_budget()) { rewrite_budget()->require_nodes(1); }
    set_result(node.clone());
}

void NormalizationVisitor::visit(const ComplexNode& node) {
    node.real()->accept(*this);
    auto real = result;
    node.imag()->accept(*this);
    normalization_check_children(rewrite_budget(), 1, real, result);
    set_result(SymbolicFactory::create_complex(real, result));
}

void NormalizationVisitor::merge_complex_terms(
    std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    std::vector<std::shared_ptr<const SymbolicNode>> real;
    std::vector<std::shared_ptr<const SymbolicNode>> imag;
    std::vector<std::shared_ptr<const SymbolicNode>> remaining;
    std::size_t real_nodes = 0, imag_nodes = 0, remaining_nodes = 0;
    bool has_complex = false;
    for (const auto& operand : operands) {
        if (const auto complex = std::dynamic_pointer_cast<const ComplexNode>(operand)) {
            has_complex = true;
            normalization_append(rewrite_budget(), real_nodes, real, complex->real());
            normalization_append(rewrite_budget(), imag_nodes, imag, complex->imag());
        } else if (std::dynamic_pointer_cast<const NumberNode>(operand)) {
            normalization_append(rewrite_budget(), real_nodes, real, operand);
            normalization_append(rewrite_budget(), imag_nodes, imag,
                                 SymbolicFactory::create_number(BigInt(0)));
        } else {
            normalization_append(rewrite_budget(), remaining_nodes, remaining, operand);
        }
    }
    if (!has_complex) {
        return;
    }
    auto real_sum = SymbolicFactory::create_add(real, rewrite_budget());
    auto imag_sum = SymbolicFactory::create_add(imag, rewrite_budget());
    NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
    real_sum->accept(visitor);
    auto normalized_real = visitor.get_result();
    imag_sum->accept(visitor);
    normalization_check_children(rewrite_budget(), 1, normalized_real, visitor.get_result());
    auto merged = SymbolicFactory::create_complex(normalized_real, visitor.get_result());
    if (!merged->is_zero()) {
        normalization_append(rewrite_budget(), remaining_nodes, remaining, merged);
    }
    operands = std::move(remaining);
}

std::shared_ptr<const NumberNode> NormalizationVisitor::extract_term_coefficient(
    std::shared_ptr<const SymbolicNode>& term) {
    auto coefficient = detail::make_node<NumberNode>(BigInt(1));
    const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(term);
    if (!multiply || multiply->operands().empty()) {
        return coefficient;
    }
    const auto& operands = multiply->operands();
    auto number = std::dynamic_pointer_cast<const NumberNode>(operands.front());
    std::size_t index = 0;
    if (!number) {
        number = std::dynamic_pointer_cast<const NumberNode>(operands.back());
        index = operands.size() - 1;
    }
    if (!number) {
        return coefficient;
    }
    coefficient = number;
    if (operands.size() == 2) {
        term = operands[index == 0 ? 1 : 0];
        return coefficient;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> rest;
    std::size_t nodes = 0;
    normalization_check_count(rewrite_budget(), operands.size() - 1, 1, false);
    rest.reserve(operands.size() - 1);
    for (std::size_t k = 0; k < operands.size(); ++k) {
        if (k != index) {
            normalization_append(rewrite_budget(), nodes, rest, operands[k]);
        }
    }
    auto rest_product = make_normalized_multiply_node(rest, rewrite_budget());
    NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
    rest_product->accept(visitor);
    term = visitor.get_result();
    return coefficient;
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::collect_sum(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    auto constant = detail::make_node<NumberNode>(BigInt(0));
    struct TermAccum {
        std::shared_ptr<const NumberNode> coefficient;
        std::vector<std::shared_ptr<const SymbolicNode>> originals;
        std::size_t nodes;
    };
    std::map<std::shared_ptr<const SymbolicNode>, TermAccum, NodeCompare> terms;
    for (const auto& operand : operands) {
        if (const auto number = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            constant = add_numbers(constant, number);
            continue;
        }
        auto term = operand;
        auto coefficient = extract_term_coefficient(term);
        auto found = terms.find(term);
        if (found == terms.end()) {
            const auto nodes = rewrite_budget() ? rewrite_budget()->measure(operand) : 0;
            terms.emplace(term, TermAccum{coefficient, {operand}, nodes});
        } else {
            found->second.coefficient = add_numbers(found->second.coefficient, coefficient);
            normalization_append(rewrite_budget(), found->second.nodes, found->second.originals, operand);
        }
    }
    std::vector<std::shared_ptr<const SymbolicNode>> collected;
    std::size_t nodes = 0;
    if (!constant->is_zero()) {
        normalization_append(rewrite_budget(), nodes, collected, constant);
    }
    for (const auto& [term, accumulated] : terms) {
        if (!accumulated.coefficient->is_zero()) {
            normalization_append(rewrite_budget(), nodes, collected,
                                 weighted_term(term, accumulated.coefficient, rewrite_budget()));
        } else if (!normalization_can_discard(term, facts_, domain_, context_)) {
            for (const auto& original : accumulated.originals) {
                normalization_append(rewrite_budget(), nodes, collected, original);
            }
        }
    }
    if (collected.empty()) {
        return detail::make_node<NumberNode>(BigInt(0));
    }
    if (collected.size() == 1) {
        return collected.front();
    }
    normalization_check_arithmetic<AddNode>(rewrite_budget(), collected, nodes);
    std::sort(collected.begin(), collected.end(), sum_term_less);
    return detail::make_node<AddNode>(collected);
}

void NormalizationVisitor::visit(const AddNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    std::size_t nodes = 0;
    for (const auto& operand : node.operands()) {
        operand->accept(*this);
        if (const auto add = std::dynamic_pointer_cast<const AddNode>(result)) {
            for (const auto& child : add->operands()) {
                normalization_append(rewrite_budget(), nodes, operands, child);
            }
        } else {
            normalization_append(rewrite_budget(), nodes, operands, result);
        }
    }
    merge_complex_terms(operands);
    if (normalize_matrix_sum(operands)) {
        return;
    }
    set_result(collect_sum(operands));
}
}
