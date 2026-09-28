#include "internal/numeric_evaluation_support.hpp"
#include "internal/squared_norm.hpp"
#include <cmath>
#include <limits>

namespace LMCAS {
namespace {

constexpr const char* kOperation = "evaluate_numeric";

class RecursionScope {
public:
    explicit RecursionScope(ComputationContext& context) : context_(context) {}
    ~RecursionScope() { context_.leave_recursion(); }

private:
    ComputationContext& context_;
};

}

namespace detail {

Result<ApproxReal> numeric_failure(CasErrc code, std::string message) {
    return Result<ApproxReal>::failure(code, std::move(message), kOperation);
}

Result<ApproxReal> numeric_approximation(double value) {
    if (!std::isfinite(value)) {
        return numeric_failure(CasErrc::NumericFailure, "numeric evaluation produced a nonfinite result");
    }
    return Result<ApproxReal>::success(
        ApproxReal{value, std::numeric_limits<double>::infinity(), NumericStatus::Finite});
}

Result<ApproxReal> NumericEvaluator::evaluate(const std::shared_ptr<const SymbolicNode>& node) {
    auto entered = context_.enter_recursion(kOperation);
    if (!entered) {
        return Result<ApproxReal>::failure(entered.error());
    }
    RecursionScope scope(context_);
    if (!node) {
        return numeric_failure(CasErrc::InvalidArgument, "expression contains a null node");
    }
    if (const auto* value = dynamic_cast<const NumberNode*>(node.get())) {
        return number(*value);
    }
    if (const auto* value = dynamic_cast<const VariableNode*>(node.get())) {
        return variable(*value);
    }
    if (const auto* value = dynamic_cast<const AddNode*>(node.get())) {
        return sum(*value);
    }
    if (const auto* value = dynamic_cast<const MultiplyNode*>(node.get())) {
        return product(*value);
    }
    const auto squares = squared_norm_terms(*node);
    if (squares[0]) {
        return norm(squares);
    }
    if (const auto* value = dynamic_cast<const PowerNode*>(node.get())) {
        return power(*value);
    }
    if (dynamic_cast<const RootOfNode*>(node.get())) {
        return root(node);
    }
    if (const auto* value = dynamic_cast<const FunctionNode*>(node.get())) {
        return function(*value);
    }
    return numeric_failure(CasErrc::UnsupportedExpression,
                           "expression node is not real-numerically evaluable");
}

Result<ApproxReal> NumericEvaluator::function(const FunctionNode& node) {
    if (node.type() == FunctionNode::FuncType::Infinity && node.arguments().empty()) {
        return Result<ApproxReal>::success(ApproxReal{
            std::numeric_limits<double>::infinity(), 0.0, NumericStatus::PositiveInfinity});
    }
    if (node.type() == FunctionNode::FuncType::Atan2) {
        return atan2(node);
    }
    if (node.type() == FunctionNode::FuncType::Max || node.type() == FunctionNode::FuncType::Min) {
        return extremum(node);
    }
    if (node.type() == FunctionNode::FuncType::Log &&
        node.arguments().size() == 2) {
        auto value = evaluate(node.arguments()[0]);
        if (!value) {
            return value;
        }
        auto base = evaluate(node.arguments()[1]);
        if (!base) {
            return base;
        }
        const double numeric_value = value.value().value;
        const double numeric_base = base.value().value;
        if (!(numeric_value > 0.0)) {
            return numeric_failure(
                CasErrc::DomainError,
                "logarithm requires a positive real argument");
        }
        if (!(numeric_base > 0.0) || numeric_base == 1.0) {
            return numeric_failure(
                CasErrc::DomainError,
                "logarithm base must be positive and not equal to one");
        }
        return numeric_approximation(
            std::log(numeric_value) / std::log(numeric_base));
    }
    if (node.arguments().size() != 1) {
        return numeric_failure(CasErrc::UnsupportedExpression,
                               "function is not supported by real numeric evaluation");
    }
    auto argument = evaluate(node.arguments()[0]);
    if (!argument) {
        return argument;
    }
    return evaluate_unary_numeric(node.type(), argument.value().value);
}

}

Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                    const NumericBindings& bindings,
                                    ComputationContext& context) {
    if (!detail::node(expression)) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "cannot evaluate an empty expression", kOperation);
    }
    detail::NumericEvaluator evaluator(bindings, context);
    return evaluator.evaluate(detail::node(expression));
}

}
