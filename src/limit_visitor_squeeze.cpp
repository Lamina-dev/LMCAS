#include "internal/visitors/limit_visitor.hpp"

namespace LMCAS {

std::shared_ptr<const SymbolicNode> LimitVisitor::squeeze_product(const MultiplyNode& node) {
    bool has_bounded = false;
    std::vector<std::shared_ptr<const SymbolicNode>> other_factors;
    for (const auto& operand : node.operands()) {
        if (is_bounded(operand) || is_bounded_expression(operand)) {
            has_bounded = true;
        } else {
            other_factors.push_back(operand);
        }
    }
    if (!has_bounded || other_factors.empty()) {
        return nullptr;
    }
    if (tends_to_zero(make_product_or_one(other_factors))) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(0));
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::squeeze_sum(const AddNode& node) {
    bool has_squeezed_term = false;
    std::vector<std::shared_ptr<const SymbolicNode>> other_terms;
    for (const auto& operand : node.operands()) {
        if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(operand)) {
            auto squeezed = squeeze_product(*product);
            if (squeezed && squeezed->is_zero()) {
                has_squeezed_term = true;
                continue;
            }
        }
        other_terms.push_back(operand);
    }
    if (!has_squeezed_term || other_terms.empty()) {
        return nullptr;
    }
    auto remaining = other_terms.size() == 1 ? other_terms[0] :
        std::static_pointer_cast<const SymbolicNode>(LMCAS::detail::make_node<AddNode>(other_terms));
    auto value = eval_limit(remaining);
    if (value && !is_inf(value)) {
        return value;
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::try_squeeze(
    const std::shared_ptr<const SymbolicNode>& expr) {
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        return squeeze_product(*product);
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(expr)) {
        return squeeze_sum(*sum);
    }
    return nullptr;
}

}
