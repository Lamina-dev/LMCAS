#pragma once

#include "symbolic_ast.hpp"
#include <array>

namespace LMCAS::detail {

/** Recognize real-numeric norm forms without evaluating their squares. */
inline std::array<const PowerNode*, 2> squared_norm_terms(const SymbolicNode& root) {
    const SymbolicNode* argument = nullptr;
    if (const auto* function = dynamic_cast<const FunctionNode*>(&root)) {
        if (function->type() != FunctionNode::FuncType::Sqrt || function->arguments().size() != 1) {
            return {};
        }
        argument = function->arguments().front().get();
    } else if (const auto* power = dynamic_cast<const PowerNode*>(&root)) {
        const auto* exponent = dynamic_cast<const NumberNode*>(power->exponent().get());
        if (!exponent) return {};
        const auto& v = exponent->value();
        const bool half = std::holds_alternative<Rational>(v)
            ? std::get<Rational>(v) == Rational(1, 2)
            : std::holds_alternative<lmmc_real_t>(v) && std::get<lmmc_real_t>(v) == 0.5;
        if (!half) return {};
        argument = power->base().get();
    } else {
        return {};
    }
    const auto square = [](const SymbolicNode* node) -> const PowerNode* {
        const auto* power = dynamic_cast<const PowerNode*>(node);
        if (!power) return nullptr;
        const auto* exponent = dynamic_cast<const NumberNode*>(power->exponent().get());
        if (!exponent) return nullptr;
        const auto& value = exponent->value();
        const bool is_two = std::holds_alternative<BigInt>(value)
            ? std::get<BigInt>(value) == BigInt(2)
            : std::holds_alternative<Rational>(value)
                ? std::get<Rational>(value) == Rational(2)
                : std::get<lmmc_real_t>(value) == 2.0;
        return is_two ? power : nullptr;
    };
    if (const auto* single = square(argument)) return {single, nullptr};
    const auto* sum = dynamic_cast<const AddNode*>(argument);
    if (!sum || sum->operands().size() != 2) return {};
    const auto* first = square(sum->operands()[0].get());
    const auto* second = square(sum->operands()[1].get());
    return first && second ? std::array<const PowerNode*, 2>{first, second}
                           : std::array<const PowerNode*, 2>{};
}

} // namespace LMCAS::detail
