#include "internal/normalization_utils.hpp"

namespace LMCAS {
namespace {
std::shared_ptr<const NumberNode> positive_negative_number(const NumberNode& number) {
    if (const auto* real = std::get_if<lmmc_real_t>(&number.value())) {
        if (*real < 0) {
            return detail::make_node<NumberNode>(std::abs(*real));
        }
        return nullptr;
    }
    if (const auto* integer = std::get_if<BigInt>(&number.value())) {
        if (integer->is_negative()) {
            return detail::make_node<NumberNode>(*integer * BigInt(-1));
        }
        return nullptr;
    }
    const auto& rational = std::get<Rational>(number.value());
    if (rational.get_numerator().is_negative()) {
        return detail::make_node<NumberNode>(rational * Rational(-1));
    }
    return nullptr;
}

}

bool check_negative_arg(const std::shared_ptr<const SymbolicNode>& argument,
                        std::shared_ptr<const SymbolicNode>& positive, detail::RewriteBudget* budget) {
    if (budget) { budget->require_nodes(1); }
    if (const auto number = std::dynamic_pointer_cast<const NumberNode>(argument)) {
        auto absolute = positive_negative_number(*number);
        if (!absolute) {
            return false;
        }
        positive = std::move(absolute);
        return true;
    }
    const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(argument);
    if (!multiply || multiply->operands().empty()) {
        return false;
    }
    const auto number = std::dynamic_pointer_cast<const NumberNode>(multiply->operands().front());
    if (!number) {
        return false;
    }
    auto absolute = positive_negative_number(*number);
    if (!absolute) {
        return false;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    std::size_t nodes = 0;
    for (const auto& operand : multiply->operands()) {
        normalization_append(budget, nodes, operands, operand);
    }
    if (number->is_negative_one()) {
        operands.erase(operands.begin());
        if (budget) { --nodes; }
    } else {
        operands.front() = std::move(absolute);
    }
    positive = make_normalized_multiply_node(operands, budget);
    return true;
}
}
