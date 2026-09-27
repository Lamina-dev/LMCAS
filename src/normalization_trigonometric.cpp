#include "internal/normalization_utils.hpp"

namespace LMCAS {
namespace {
std::shared_ptr<const SymbolicNode> scaled_radical(
    const Rational& coefficient, const BigInt& radicand, detail::RewriteBudget* budget) {
    if (budget) { budget->require_nodes(coefficient == Rational(1) ? 3 : 5); }
    auto root = detail::make_node<PowerNode>(detail::make_node<NumberNode>(radicand),
        detail::make_node<NumberNode>(Rational(1, 2)));
    if (coefficient == Rational(1)) {
        return root;
    }
    if (coefficient == Rational(-1)) {
        return make_normalized_multiply_node({detail::make_node<NumberNode>(BigInt(-1)), root}, budget);
    }
    return make_normalized_multiply_node({detail::make_node<NumberNode>(coefficient), root}, budget);
}

std::shared_ptr<const SymbolicNode> thirds_pi(
    FunctionNode::FuncType type, const BigInt& numerator, detail::RewriteBudget* budget) {
    switch (type) {
        case FunctionNode::FuncType::Sin:
            return scaled_radical(numerator < BigInt(3) ? Rational(1, 2) : Rational(-1, 2), BigInt(3), budget);
        case FunctionNode::FuncType::Cos: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(5);
            return detail::make_node<NumberNode>(positive ? Rational(1, 2) : Rational(-1, 2));
        }
        case FunctionNode::FuncType::Tan: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(4);
            return scaled_radical(positive ? Rational(1) : Rational(-1), BigInt(3), budget);
        }
        default:
            return nullptr;
    }
}

std::shared_ptr<const SymbolicNode> quarters_pi(
    FunctionNode::FuncType type, const BigInt& numerator, detail::RewriteBudget* budget) {
    switch (type) {
        case FunctionNode::FuncType::Sin:
            return scaled_radical(numerator < BigInt(4) ? Rational(1, 2) : Rational(-1, 2), BigInt(2), budget);
        case FunctionNode::FuncType::Cos: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(7);
            return scaled_radical(positive ? Rational(1, 2) : Rational(-1, 2), BigInt(2), budget);
        }
        case FunctionNode::FuncType::Tan: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(5);
            return detail::make_node<NumberNode>(BigInt(positive ? 1 : -1));
        }
        default:
            return nullptr;
    }
}

std::shared_ptr<const SymbolicNode> sixths_pi(
    FunctionNode::FuncType type, const BigInt& numerator, detail::RewriteBudget* budget) {
    switch (type) {
        case FunctionNode::FuncType::Sin:
            return detail::make_node<NumberNode>(numerator < BigInt(6) ? Rational(1, 2) : Rational(-1, 2));
        case FunctionNode::FuncType::Cos: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(11);
            return scaled_radical(positive ? Rational(1, 2) : Rational(-1, 2), BigInt(3), budget);
        }
        case FunctionNode::FuncType::Tan: {
            const bool positive = numerator == BigInt(1) || numerator == BigInt(7);
            return scaled_radical(positive ? Rational(1, 3) : Rational(-1, 3), BigInt(3), budget);
        }
        default:
            return nullptr;
    }
}

std::shared_ptr<const SymbolicNode> halves_pi(FunctionNode::FuncType type, const BigInt& numerator) {
    if (type == FunctionNode::FuncType::Sin) {
        return detail::make_node<NumberNode>(BigInt(numerator == BigInt(1) ? 1 : -1));
    }
    if (type == FunctionNode::FuncType::Cos) {
        return detail::make_node<NumberNode>(BigInt(0));
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> integer_pi(FunctionNode::FuncType type, const BigInt& numerator) {
    if (type == FunctionNode::FuncType::Sin || type == FunctionNode::FuncType::Tan) {
        return detail::make_node<NumberNode>(BigInt(0));
    }
    if (type == FunctionNode::FuncType::Cos) {
        return detail::make_node<NumberNode>(BigInt(numerator.is_zero() ? 1 : -1));
    }
    return nullptr;
}
}

std::shared_ptr<const SymbolicNode> normalization_pi_function(
    FunctionNode::FuncType type, const Rational& coefficient, detail::RewriteBudget* budget) {
    if (budget) { budget->require_nodes(1); }
    const BigInt& denominator = coefficient.get_denominator();
    const BigInt period = denominator * BigInt(2);
    BigInt numerator = coefficient.get_numerator() % period;
    if (numerator.is_negative()) numerator += period;
    if (denominator == BigInt(1)) {
        return integer_pi(type, numerator);
    }
    if (denominator == BigInt(2)) {
        return halves_pi(type, numerator);
    }
    if (denominator == BigInt(3)) {
        return thirds_pi(type, numerator, budget);
    }
    if (denominator == BigInt(4)) {
        return quarters_pi(type, numerator, budget);
    }
    if (denominator == BigInt(6)) {
        return sixths_pi(type, numerator, budget);
    }
    return nullptr;
}
}
