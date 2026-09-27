#include "internal/normalization_utils.hpp"
#include <cmath>

namespace LMCAS {
namespace {
bool number_equals_integer(
    const std::shared_ptr<const NumberNode>& number, const BigInt& value) {
    BigInt integer;
    return try_get_integer_value(number, integer) && integer == value;
}

std::shared_ptr<const SymbolicNode> function_at_zero(FunctionNode::FuncType type) {
    switch (type) {
        case FunctionNode::FuncType::Sin:
        case FunctionNode::FuncType::Tan:
            return detail::make_node<NumberNode>(BigInt(0));
        case FunctionNode::FuncType::Cos:
        case FunctionNode::FuncType::Exp:
            return detail::make_node<NumberNode>(BigInt(1));
        default:
            return nullptr;
    }
}

std::shared_ptr<const SymbolicNode> numeric_absolute(const NumberNode& number) {
    if (const auto* integer = std::get_if<BigInt>(&number.value())) {
        return detail::make_node<NumberNode>(integer->abs());
    }
    if (const auto* rational = std::get_if<Rational>(&number.value())) {
        if (rational->get_numerator().is_negative()) {
            return detail::make_node<NumberNode>(Rational(
                rational->get_numerator().abs(), rational->get_denominator()));
        }
        return detail::make_node<NumberNode>(*rational);
    }
    return detail::make_node<NumberNode>(std::abs(std::get<lmmc_real_t>(number.value())));
}

std::shared_ptr<const SymbolicNode> numeric_square_root(
    const std::shared_ptr<const NumberNode>& number, detail::RewriteBudget* budget) {
    if (const auto* real = std::get_if<lmmc_real_t>(&number->value())) {
        if (*real >= 0) {
            return detail::make_node<NumberNode>(std::sqrt(*real));
        }
        return nullptr;
    }
    if (auto root = normalization_exact_square_root(*number)) {
        return root;
    }
    normalization_check_children(budget, 1, number);
    return detail::make_node<FunctionNode>(FunctionNode::FuncType::Sqrt,
        std::vector<std::shared_ptr<const SymbolicNode>>{number});
}

std::shared_ptr<const SymbolicNode> numeric_lambert_w(
    const std::shared_ptr<const NumberNode>& number) {
    if (const auto* real = std::get_if<lmmc_real_t>(&number->value())) {
        lmmc_real_t value;
        if (lmmc_lambertw(*real, &value) == LMMC_STATUS_OK && std::isfinite(value)) {
            return detail::make_node<NumberNode>(value);
        }
        return nullptr;
    }
    if (number_equals_integer(number, BigInt(0))) {
        return detail::make_node<NumberNode>(BigInt(0));
    }
    return nullptr;
}
}

std::shared_ptr<const SymbolicNode> normalization_numeric_function(
    FunctionNode::FuncType type, const std::shared_ptr<const NumberNode>& number,
    detail::RewriteBudget* budget) {
    if (budget) { budget->require_nodes(1); }
    if (number_equals_integer(number, BigInt(0))) {
        if (auto value = function_at_zero(type)) {
            return value;
        }
    }
    if (type == FunctionNode::FuncType::Ln &&
        number_equals_integer(number, BigInt(1))) {
        return detail::make_node<NumberNode>(BigInt(0));
    }
    switch (type) {
        case FunctionNode::FuncType::Abs: return numeric_absolute(*number);
        case FunctionNode::FuncType::Sqrt: return numeric_square_root(number, budget);
        case FunctionNode::FuncType::LambertW: return numeric_lambert_w(number);
        default: return nullptr;
    }
}
}
