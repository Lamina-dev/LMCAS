#include "numeric_evaluation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

namespace lamina {
namespace {

constexpr const char* kOperation = "evaluate_numeric";

class RecursionScope {
public:
    explicit RecursionScope(ComputationContext& context) : context_(context) {}
    ~RecursionScope() { context_.leave_recursion(); }

private:
    ComputationContext& context_;
};

Result<ApproxReal> failure(CasErrc code, std::string message) {
    return Result<ApproxReal>::failure(code, std::move(message), kOperation);
}

Result<ApproxReal> make_approx(double value) {
    if (std::isnan(value)) {
        return failure(CasErrc::NumericFailure, "numeric evaluation produced NaN");
    }

    ApproxReal result;
    result.value = value;
    if (std::isinf(value)) {
        result.status = value > 0 ? NumericStatus::PositiveInfinity
                                  : NumericStatus::NegativeInfinity;
        result.absolute_error = 0.0;
    } else {
        result.status = NumericStatus::Finite;
        result.absolute_error = std::numeric_limits<double>::epsilon() *
                                std::max(1.0, std::abs(value)) * 4.0;
    }
    return Result<ApproxReal>::success(result);
}

Result<ApproxReal> evaluate_node(const std::shared_ptr<const SymbolicNode>& node,
                                 const NumericBindings& bindings,
                                 ComputationContext& context) {
    auto entered = context.enter_recursion(kOperation);
    if (!entered) return Result<ApproxReal>::failure(entered.error());
    RecursionScope scope(context);

    if (!node) {
        return failure(CasErrc::InvalidArgument, "expression contains a null node");
    }

    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        double value = 0.0;
        try {
            if (std::holds_alternative<lmmc_real_t>(number->value())) {
                value = static_cast<double>(std::get<lmmc_real_t>(number->value()));
            } else if (std::holds_alternative<BigInt>(number->value())) {
                value = std::get<BigInt>(number->value()).to_double();
            } else {
                value = std::get<Rational>(number->value()).to_double();
            }
        } catch (const std::exception& error) {
            return failure(CasErrc::NumericFailure, error.what());
        }
        if (!std::isfinite(value)) {
            return failure(CasErrc::NumericFailure,
                           "exact number cannot be represented as a finite double");
        }
        return make_approx(value);
    }

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (variable->name() == "pi" || variable->name() == "π") {
            return make_approx(static_cast<double>(LMMC_CONST_PI));
        }
        if (variable->name() == "e") {
            return make_approx(std::exp(1.0));
        }
        if (variable->name() == "phi") {
            return make_approx((1.0 + std::sqrt(5.0)) / 2.0);
        }
        auto it = bindings.find(variable->name());
        if (it == bindings.end()) {
            return failure(CasErrc::UnboundSymbol,
                           "no numeric binding for symbol '" + variable->name() + "'");
        }
        return make_approx(it->second);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        double value = 0.0;
        for (const auto& operand : add->operands()) {
            auto evaluated = evaluate_node(operand, bindings, context);
            if (!evaluated) return evaluated;
            value += evaluated.value().value;
        }
        return make_approx(value);
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        double value = 1.0;
        for (const auto& operand : multiply->operands()) {
            auto evaluated = evaluate_node(operand, bindings, context);
            if (!evaluated) return evaluated;
            value *= evaluated.value().value;
        }
        return make_approx(value);
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = evaluate_node(power->base(), bindings, context);
        if (!base) return base;
        auto exponent = evaluate_node(power->exponent(), bindings, context);
        if (!exponent) return exponent;
        if (base.value().value == 0.0 && exponent.value().value <= 0.0) {
            return failure(CasErrc::DomainError,
                           "zero cannot be raised to a non-positive power");
        }
        return make_approx(std::pow(base.value().value, exponent.value().value));
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Infinity &&
            function->arguments().empty()) {
            return make_approx(std::numeric_limits<double>::infinity());
        }

        if (function->type() == FunctionNode::FuncType::Atan2) {
            if (function->arguments().size() != 2) {
                return failure(CasErrc::InvalidArgument, "atan2 requires two arguments");
            }
            auto y = evaluate_node(function->arguments()[0], bindings, context);
            if (!y) return y;
            auto x = evaluate_node(function->arguments()[1], bindings, context);
            if (!x) return x;
            return make_approx(std::atan2(y.value().value, x.value().value));
        }

        if (function->type() == FunctionNode::FuncType::Max ||
            function->type() == FunctionNode::FuncType::Min) {
            if (function->arguments().empty()) {
                return failure(CasErrc::InvalidArgument, "min/max requires an argument");
            }
            auto first = evaluate_node(function->arguments().front(), bindings, context);
            if (!first) return first;
            double value = first.value().value;
            for (std::size_t i = 1; i < function->arguments().size(); ++i) {
                auto next = evaluate_node(function->arguments()[i], bindings, context);
                if (!next) return next;
                value = function->type() == FunctionNode::FuncType::Max
                            ? std::max(value, next.value().value)
                            : std::min(value, next.value().value);
            }
            return make_approx(value);
        }

        if (function->arguments().size() != 1) {
            return failure(CasErrc::UnsupportedExpression,
                           "function is not supported by real numeric evaluation");
        }
        auto argument = evaluate_node(function->arguments()[0], bindings, context);
        if (!argument) return argument;
        const double x = argument.value().value;
        double result = 0.0;

        switch (function->type()) {
            case FunctionNode::FuncType::Sin: result = std::sin(x); break;
            case FunctionNode::FuncType::Cos: result = std::cos(x); break;
            case FunctionNode::FuncType::Tan: result = std::tan(x); break;
            case FunctionNode::FuncType::Cot:
                if (std::sin(x) == 0.0) return failure(CasErrc::DomainError, "cot is undefined");
                result = 1.0 / std::tan(x);
                break;
            case FunctionNode::FuncType::Sec:
                if (std::cos(x) == 0.0) return failure(CasErrc::DomainError, "sec is undefined");
                result = 1.0 / std::cos(x);
                break;
            case FunctionNode::FuncType::Csc:
                if (std::sin(x) == 0.0) return failure(CasErrc::DomainError, "csc is undefined");
                result = 1.0 / std::sin(x);
                break;
            case FunctionNode::FuncType::ArcSin:
                if (x < -1.0 || x > 1.0) return failure(CasErrc::DomainError, "asin real domain is [-1, 1]");
                result = std::asin(x);
                break;
            case FunctionNode::FuncType::ArcCos:
                if (x < -1.0 || x > 1.0) return failure(CasErrc::DomainError, "acos real domain is [-1, 1]");
                result = std::acos(x);
                break;
            case FunctionNode::FuncType::ArcTan: result = std::atan(x); break;
            case FunctionNode::FuncType::Sinh: result = std::sinh(x); break;
            case FunctionNode::FuncType::Cosh: result = std::cosh(x); break;
            case FunctionNode::FuncType::Tanh: result = std::tanh(x); break;
            case FunctionNode::FuncType::Ln:
            case FunctionNode::FuncType::Log:
                if (x <= 0.0) return failure(CasErrc::DomainError, "logarithm requires a positive real argument");
                result = std::log(x);
                break;
            case FunctionNode::FuncType::Abs: result = std::abs(x); break;
            case FunctionNode::FuncType::Sqrt:
                if (x < 0.0) return failure(CasErrc::DomainError, "real square root requires a non-negative argument");
                result = std::sqrt(x);
                break;
            case FunctionNode::FuncType::Exp: result = std::exp(x); break;
            case FunctionNode::FuncType::Sgn: result = (x > 0.0) - (x < 0.0); break;
            case FunctionNode::FuncType::Floor: result = std::floor(x); break;
            case FunctionNode::FuncType::Ceil: result = std::ceil(x); break;
            case FunctionNode::FuncType::Round: result = std::round(x); break;
            case FunctionNode::FuncType::Erf: result = std::erf(x); break;
            case FunctionNode::FuncType::LambertW: {
                lmmc_real_t value = 0.0;
                if (lmmc_lambertw(static_cast<lmmc_real_t>(x), &value) != LMMC_STATUS_OK) {
                    return failure(CasErrc::DomainError, "LambertW evaluation failed on the real branch");
                }
                result = static_cast<double>(value);
                break;
            }
            default:
                return failure(CasErrc::UnsupportedExpression,
                               "function is not supported by real numeric evaluation");
        }
        return make_approx(result);
    }

    return failure(CasErrc::UnsupportedExpression,
                   "expression node is not real-numerically evaluable");
}

} // namespace

Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                    const NumericBindings& bindings,
                                    ComputationContext& context) {
    if (!lamina::detail::node(expression)) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "cannot evaluate an empty expression",
                                           kOperation);
    }
    return evaluate_node(lamina::detail::node(expression), bindings, context);
}

} // namespace lamina
