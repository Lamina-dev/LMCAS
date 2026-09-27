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

Result<std::shared_ptr<SymbolicExpr>> Integrator::apply_linearity(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth) {

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        std::vector<std::shared_ptr<const SymbolicNode>> constants;
        std::vector<std::shared_ptr<const SymbolicNode>> dependents;

        for (auto& op : mul->operands()) {
            auto term = LMCAS::detail::expression_from_node(op);
            if (!depends_on_integration_variable(term, var)) {
                constants.push_back(op);
            } else {
                dependents.push_back(op);
            }
        }

        auto exponential = combine_integration_exponentials(dependents);
        if (!constants.empty() || exponential) {
            SymbolicExpr dep_part = exponential
                ? detail::expression_from_node(exponential)
                : dependents.empty() ? *SymbolicExpr::number(1)
                : dependents.size() == 1 ? detail::expression_from_node(dependents[0])
                : detail::expression_from_node(detail::make_node<MultiplyNode>(dependents));

            auto int_part = integrate_recursive(dep_part, var, context, depth + 1);
            if (!int_part) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    int_part.error());
            }
            if (constants.empty()) { return int_part; }
            auto const_part = constants.size() == 1
                ? detail::expression_from_node(constants[0])
                : detail::expression_from_node(detail::make_node<MultiplyNode>(constants));
            return Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::multiply(
                    LMCAS::detail::make_expression_ptr(const_part),
                    int_part.value()));
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(
            LMCAS::detail::node(expr))) {
        std::vector<std::shared_ptr<const SymbolicNode>> results;
        for (auto& op : add->operands()) {
            auto term = LMCAS::detail::expression_from_node(op);
            auto int_term = integrate_recursive(
                term, var, context, depth + 1);
            if (!int_term) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    int_term.error());
            }
            results.push_back(LMCAS::detail::node(int_term.value()));
        }
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<AddNode>(results)));
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(nullptr);
}

}
