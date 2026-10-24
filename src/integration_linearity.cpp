#include "internal/integration_support.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

namespace {
std::shared_ptr<const SymbolicNode> integration_exponential_argument(
    const std::shared_ptr<const SymbolicNode>& factor) {
    const auto* power = dynamic_cast<const PowerNode*>(factor.get());
    const auto& base = power ? power->base() : factor;
    const auto* function = dynamic_cast<const FunctionNode*>(base.get());
    if (!function || function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().size() != 1) {
        return nullptr;
    }
    if (!power) { return function->arguments()[0]; }
    BigInt integer;
    if (!try_get_integer_value(
            std::dynamic_pointer_cast<const NumberNode>(power->exponent()), integer) ||
        integer.is_zero()) {
        return nullptr;
    }
    return make_normalized_multiply_node(
        {power->exponent(), function->arguments()[0]});
}

std::shared_ptr<const SymbolicNode> combine_integration_exponentials(
    const std::vector<std::shared_ptr<const SymbolicNode>>& dependents) {
    if (dependents.empty() || (dependents.size() == 1 &&
        !dynamic_cast<const PowerNode*>(dependents[0].get()))) {
        return nullptr;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> exponents;
    exponents.reserve(dependents.size());
    for (const auto& dependent : dependents) {
        auto argument = integration_exponential_argument(dependent);
        if (!argument) { return nullptr; }
        exponents.push_back(std::move(argument));
    }
    return detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Exp,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            SymbolicFactory::create_add(std::move(exponents))});
}
}

static detail::SymbolicNodePtr scale_integral_field(
    const detail::SymbolicNodePtr& value, const detail::SymbolicNodePtr& coefficient,
    detail::RewriteBudget& budget) {
    detail::RewriteScope scope(&budget);
    if (const auto* piecewise = dynamic_cast<const PiecewiseNode*>(value.get())) {
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(piecewise->branches().size());
        for (const auto& branch : piecewise->branches()) {
            branches.push_back({scale_integral_field(branch.expression, coefficient, budget),
                branch.condition});
        }
        auto fallback = piecewise->default_expr()
            ? scale_integral_field(piecewise->default_expr(), coefficient, budget) : nullptr;
        auto result = detail::make_node<PiecewiseNode>(std::move(branches), fallback);
        budget.measure(result);
        return result;
    }
    return detail::normalize_node(
        SymbolicFactory::create_multiply({coefficient, value}, &budget), budget);
}

static bool zero_remainder(const detail::SymbolicNodePtr& value, detail::RewriteBudget& budget) {
    detail::RewriteScope scope(&budget);
    if (value->is_zero()) { return true; }
    const auto* piecewise = dynamic_cast<const PiecewiseNode*>(value.get());
    if (!piecewise) { return false; }
    for (const auto& branch : piecewise->branches()) {
        if (!zero_remainder(branch.expression, budget)) { return false; }
    }
    return !piecewise->default_expr() || zero_remainder(piecewise->default_expr(), budget);
}

Result<std::optional<detail::IntegralParts>> Integrator::apply_linearity(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth) {
    using Attempt = Result<std::optional<detail::IntegralParts>>;
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.linearity");
    if (const auto* add = dynamic_cast<const AddNode*>(detail::node(expr).get())) {
        std::vector<detail::SymbolicNodePtr> primitives, remainders;
        for (const auto& operand : add->operands()) {
            auto term = detail::expression_from_node(operand);
            auto parts = integrate_recursive(term, var, context, depth + 1, &term);
            if (!parts) { return Attempt::failure(parts.error()); }
            primitives.push_back(detail::node(parts.value().primitive));
            if (!zero_remainder(detail::node(parts.value().remainder), budget)) {
                remainders.push_back(detail::node(parts.value().remainder));
            }
        }
        const auto remainder = remainders.empty() ? detail::node(zero_)
            : remainders.size() == 1 ? remainders.front()
            : detail::normalize_node(
                SymbolicFactory::create_add(std::move(remainders), &budget), budget);
        return Attempt::success(detail::IntegralParts{
            detail::expression_from_node(detail::normalize_node(
                SymbolicFactory::create_add(std::move(primitives), &budget), budget)),
            detail::expression_from_node(remainder)});
    }
    if (const auto* mul = dynamic_cast<const MultiplyNode*>(detail::node(expr).get())) {
        std::vector<detail::SymbolicNodePtr> constants, dependents;
        for (const auto& operand : mul->operands()) {
            if (detail::integration_depends_on_variable(operand, var, &budget)) {
                dependents.push_back(operand);
            }
            else { constants.push_back(operand); }
        }
        auto exponential = combine_integration_exponentials(dependents);
        if (!constants.empty() || exponential) {
            auto dependent = exponential ? detail::expression_from_node(exponential)
                : dependents.empty() ? *SymbolicExpr::number(1)
                : dependents.size() == 1 ? detail::expression_from_node(dependents.front())
                : detail::expression_from_node(detail::make_node<MultiplyNode>(dependents));
            auto parts = integrate_recursive(dependent, var, context, depth + 1, &dependent);
            if (!parts) { return Attempt::failure(parts.error()); }
            if (constants.empty()) { return Attempt::success(std::move(parts.value())); }
            auto coefficient = constants.size() == 1 ? constants.front()
                : detail::make_node<MultiplyNode>(std::move(constants));
            coefficient = detail::normalize_node(coefficient, budget);
            auto scale = [&](const SymbolicExpr& value) {
                return detail::expression_from_node(scale_integral_field(
                    detail::node(value), coefficient, budget));
            };
            return Attempt::success(detail::IntegralParts{
                scale(parts.value().primitive), scale(parts.value().remainder)});
        }
    }
    return Attempt::success(std::nullopt);
}

}
