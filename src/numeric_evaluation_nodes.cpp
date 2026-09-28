#include "internal/numeric_evaluation_support.hpp"
#include "root_of_identity.hpp"
#include <algorithm>
#include <cmath>

namespace LMCAS::detail {

Result<ApproxReal> NumericEvaluator::number(const NumberNode& node) {
    double value = 0.0;
    try {
        if (std::holds_alternative<lmmc_real_t>(node.value())) {
            value = static_cast<double>(std::get<lmmc_real_t>(node.value()));
        } else if (std::holds_alternative<BigInt>(node.value())) {
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
    if (std::holds_alternative<lmmc_real_t>(node.value())) {
        return Result<ApproxReal>::success(ApproxReal{value, 0.0, NumericStatus::Finite});
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
    double value = 1.0;
    bool infinite = false;
    for (const auto& operand : node.operands()) {
        auto evaluated = evaluate(operand);
        if (!evaluated) return evaluated;
        value *= evaluated.value().value;
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
