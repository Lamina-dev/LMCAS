#include "internal/numeric_evaluation_support.hpp"
#include "root_of_identity.hpp"
#include <algorithm>
#include <cmath>

namespace LMCAS::detail {

Result<ApproxReal> NumericEvaluator::number(const NumberNode& node) {
    double value = 0.0;
    try {
        if (std::holds_alternative<BigInt>(node.value())) {
            value = std::get<BigInt>(node.value()).to_double();
        } else {
            value = std::get<Rational>(node.value()).to_double();
        }
    } catch (const std::exception& error) {
        return numeric_failure(CasErrc::NumericFailure, error.what());
    }
    if (!std::isfinite(value)) {
        return numeric_failure(CasErrc::NumericFailure,
                       "exact number cannot be represented as a finite double");
    }
    return numeric_approximation(value);
}

Result<ApproxReal> NumericEvaluator::variable(const VariableNode& node) {
    if (node.is_constant() && (node.name() == "pi" || node.name() == "π")) {
        return numeric_approximation(static_cast<double>(LMMC_CONST_PI));
    }
    if (node.is_constant() && node.name() == "e") {
        return numeric_approximation(std::exp(1.0));
    }
    if (node.is_constant() && node.name() == "phi") {
        return numeric_approximation((1.0 + std::sqrt(5.0)) / 2.0);
    }
    auto it = bindings_.find(node.name());
    if (it == bindings_.end()) {
        return numeric_failure(CasErrc::UnboundSymbol,
                       "no numeric binding for symbol '" + node.name() + "'");
    }
    if (!std::isfinite(it->second)) {
        return numeric_failure(CasErrc::NumericFailure, "numeric binding is not finite");
    }
    return Result<ApproxReal>::success(
        ApproxReal{it->second, 0.0, NumericStatus::Finite});
}

Result<ApproxReal> NumericEvaluator::sum(const AddNode& node) {
    double value = 0.0;
    bool infinite = false;
    for (const auto& operand : node.operands()) {
        auto evaluated = evaluate(operand);
        if (!evaluated) return evaluated;
        value += evaluated.value().value;
        infinite = infinite || !evaluated.value().is_finite();
        if (std::isinf(value) && !infinite) {
            return numeric_failure(CasErrc::NumericFailure, "numeric evaluation produced a nonfinite result");
        }
    }
    if (std::isinf(value) && infinite) {
        return Result<ApproxReal>::success(ApproxReal{
            value, 0.0, value > 0 ? NumericStatus::PositiveInfinity : NumericStatus::NegativeInfinity});
    }
    return numeric_approximation(value);
}

Result<ApproxReal> NumericEvaluator::product(const MultiplyNode& node) {
    const bool exact_constants = std::all_of(node.operands().begin(), node.operands().end(),
        [](const auto& operand) {
            if (dynamic_cast<const NumberNode*>(operand.get())) return true;
            const auto* power = dynamic_cast<const PowerNode*>(operand.get());
            const auto* exponent = power
                ? dynamic_cast<const NumberNode*>(power->exponent().get()) : nullptr;
            return exponent && exponent->is_negative_one() &&
                   dynamic_cast<const NumberNode*>(power->base().get());
        });
    if (exact_constants) {
        Rational value(1);
        for (const auto& operand : node.operands()) {
            const auto* power = dynamic_cast<const PowerNode*>(operand.get());
            const auto* number = dynamic_cast<const NumberNode*>(
                power ? power->base().get() : operand.get());
            auto step = context_.consume_steps(power ? 3 : 1, "evaluate_numeric");
            if (!step) return Result<ApproxReal>::failure(step.error());
            Rational factor = std::holds_alternative<BigInt>(number->value())
                ? Rational(std::get<BigInt>(number->value()))
                : std::get<Rational>(number->value());
            if (power) {
                if (factor.is_zero()) {
                    return numeric_failure(CasErrc::DomainError, "zero cannot be raised to a non-positive power");
                }
                factor = Rational(factor.get_denominator(), factor.get_numerator());
            }
            value = value * factor;
            auto numerator_bits = context_.require_integer_bits(
                value.get_numerator().bit_length(), "evaluate_numeric");
            if (!numerator_bits) return Result<ApproxReal>::failure(numerator_bits.error());
            auto denominator_bits = context_.require_integer_bits(
                value.get_denominator().bit_length(), "evaluate_numeric");
            if (!denominator_bits) return Result<ApproxReal>::failure(denominator_bits.error());
        }
        try {
            return numeric_approximation(value.to_double());
        } catch (const std::exception& error) {
            return numeric_failure(CasErrc::NumericFailure, error.what());
        }
    }
    double value = 1.0;
    bool infinite = false;
    for (const auto& operand : node.operands()) {
        auto evaluated = evaluate(operand);
        if (!evaluated) return evaluated;
        const double factor = evaluated.value().value;
        const double multiplied = value * factor;
        if (multiplied == 0.0 && value != 0.0 && factor != 0.0) {
            return numeric_failure(CasErrc::NumericFailure, "numeric multiplication underflow");
        }
        value = multiplied;
        infinite = infinite || !evaluated.value().is_finite();
        if (std::isinf(value) && !infinite) {
            return numeric_failure(CasErrc::NumericFailure, "numeric evaluation produced a nonfinite result");
        }
    }
    if (std::isinf(value) && infinite) {
        return Result<ApproxReal>::success(ApproxReal{
            value, 0.0, value > 0 ? NumericStatus::PositiveInfinity : NumericStatus::NegativeInfinity});
    }
    return numeric_approximation(value);
}

Result<ApproxReal> NumericEvaluator::atan2(const FunctionNode& node) {
        if (node.arguments().size() != 2) {
            return numeric_failure(CasErrc::InvalidArgument, "atan2 requires two arguments");
        }
        auto y = evaluate(node.arguments()[0]);
        if (!y) return y;
        auto x = evaluate(node.arguments()[1]);
        if (!x) return x;
        return numeric_approximation(std::atan2(y.value().value, x.value().value));
}

Result<ApproxReal> NumericEvaluator::extremum(const FunctionNode& node) {
        if (node.arguments().empty()) {
            return numeric_failure(CasErrc::InvalidArgument, "min/max requires an argument");
        }
        auto first = evaluate(node.arguments().front());
        if (!first) return first;
        double value = first.value().value;
        for (std::size_t i = 1; i < node.arguments().size(); ++i) {
            auto next = evaluate(node.arguments()[i]);
            if (!next) return next;
            value = node.type() == FunctionNode::FuncType::Max
                        ? std::max(value, next.value().value)
                        : std::min(value, next.value().value);
        }
        return numeric_approximation(value);
}

Result<ApproxReal> NumericEvaluator::norm(const std::array<const PowerNode*, 2>& squares) {
    auto first = evaluate(squares[0]->base());
    if (!first) return first;
    if (!squares[1]) return numeric_approximation(std::abs(first.value().value));
    auto second = evaluate(squares[1]->base());
    if (!second) return second;
    return numeric_approximation(std::hypot(first.value().value, second.value().value));
}

Result<ApproxReal> NumericEvaluator::root(const std::shared_ptr<const SymbolicNode>& node) {
    auto root_expression = LMCAS::detail::make_expression_ptr(node);
    auto root = rootof_evaluate_checked(root_expression, context_);
    if (!root) return Result<ApproxReal>::failure(root.error());
    return numeric_approximation(root.value());
}

}
