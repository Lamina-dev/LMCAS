#include "internal/squared_norm.hpp"

namespace LMCAS::detail {
namespace {

const SymbolicNode* square_root_argument(const SymbolicNode& root) {
    if (const auto* function = dynamic_cast<const FunctionNode*>(&root)) {
        if (function->type() == FunctionNode::FuncType::Sqrt && function->arguments().size() == 1) {
            return function->arguments().front().get();
        }
        return nullptr;
    }
    const auto* power = dynamic_cast<const PowerNode*>(&root);
    if (!power) {
        return nullptr;
    }
    const auto* exponent = dynamic_cast<const NumberNode*>(power->exponent().get());
    if (!exponent) {
        return nullptr;
    }
    const auto& value = exponent->value();
    const bool half = std::holds_alternative<Rational>(value)
        ? std::get<Rational>(value) == Rational(1, 2)
        : std::holds_alternative<lmmc_real_t>(value) && std::get<lmmc_real_t>(value) == 0.5;
    return half ? power->base().get() : nullptr;
}

const PowerNode* squared_term(const SymbolicNode* node) {
    const auto* power = dynamic_cast<const PowerNode*>(node);
    if (!power) {
        return nullptr;
    }
    const auto* exponent = dynamic_cast<const NumberNode*>(power->exponent().get());
    if (!exponent) {
        return nullptr;
    }
    const auto& value = exponent->value();
    const bool two = std::holds_alternative<BigInt>(value)
        ? std::get<BigInt>(value) == BigInt(2)
        : std::holds_alternative<Rational>(value)
            ? std::get<Rational>(value) == Rational(2)
            : std::get<lmmc_real_t>(value) == 2.0;
    return two ? power : nullptr;
}

}

std::array<const PowerNode*, 2> squared_norm_terms(const SymbolicNode& root) {
    const auto* argument = square_root_argument(root);
    if (const auto* single = squared_term(argument)) {
        return {single, nullptr};
    }
    const auto* sum = dynamic_cast<const AddNode*>(argument);
    if (!sum || sum->operands().size() != 2) {
        return {};
    }
    const auto* first = squared_term(sum->operands()[0].get());
    const auto* second = squared_term(sum->operands()[1].get());
    return first && second ? std::array<const PowerNode*, 2>{first, second}
                           : std::array<const PowerNode*, 2>{};
}

}
