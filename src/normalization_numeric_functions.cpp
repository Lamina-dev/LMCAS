#include "internal/normalization_utils.hpp"

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
    const auto& rational = std::get<Rational>(number.value());
    if (rational.get_numerator().is_negative()) {
        return detail::make_node<NumberNode>(Rational(
            rational.get_numerator().abs(), rational.get_denominator()));
    }
    return detail::make_node<NumberNode>(rational);
}

std::shared_ptr<const SymbolicNode> numeric_square_root(
    const std::shared_ptr<const NumberNode>& number, detail::RewriteBudget* budget) {
    if (auto root = normalization_exact_square_root(*number)) {
        return root;
    }
    normalization_check_children(budget, 1, number);
    return detail::make_node<FunctionNode>(FunctionNode::FuncType::Sqrt,
        std::vector<std::shared_ptr<const SymbolicNode>>{number});
}

std::shared_ptr<const SymbolicNode> numeric_lambert_w(
    const std::shared_ptr<const NumberNode>& number) {
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
