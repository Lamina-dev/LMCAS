#include "internal/equivalence_support.hpp"

#include "poly_utils.hpp"

#include <set>

namespace LMCAS::equivalence_detail {
namespace {

bool collect_polynomial_variable_names(
    const std::shared_ptr<const SymbolicNode>& node,
    std::set<std::string>& variables, detail::RewriteBudget& budget) {
    if (!node) { return true; }
    detail::RewriteScope scope(&budget);
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        variables.insert(variable->name());
        return true;
    }
    if (std::dynamic_pointer_cast<const NumberNode>(node)) {
        return true;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (!collect_polynomial_variable_names(operand, variables, budget)) {
                return false;
            }
        }
        return true;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (!collect_polynomial_variable_names(operand, variables, budget)) {
                return false;
            }
        }
        return true;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return collect_polynomial_variable_names(power->base(), variables, budget) &&
               collect_polynomial_variable_names(power->exponent(), variables, budget);
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return collect_polynomial_variable_names(complex_node->real(), variables, budget) &&
               collect_polynomial_variable_names(complex_node->imag(), variables, budget);
    }
    return false;
}

}

Result<std::optional<bool>> prove_rational_polynomial_equivalence(
    const ExprPtr& difference,
    ComputationContext& context,
    const EqvOptions& options, detail::RewriteBudget& budget) {
    if (!difference) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence difference is null",
            kEquivalentOperation);
    }

    if (options.budget.max_rewrite_steps < 4) {
        return Result<std::optional<bool>>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before polynomial normalization",
            kEquivalentOperation);
    }

    auto step = context.consume_steps(4, kEquivalentOperation);
    if (!step) { return Result<std::optional<bool>>::failure(step.error()); }

    std::set<std::string> variables;
    if (!collect_polynomial_variable_names(
            LMCAS::detail::node(difference), variables, budget)) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    if (variables.size() > 1) {
        detail::RewriteBudget growth_budget(
            context, options.budget.max_rewrite_depth, budget.max_nodes(),
            kEquivalentOperation);
        auto expanded = expand_equivalence(difference, growth_budget);
        return Result<std::optional<bool>>::success(
            expanded->is_zero() ? std::optional<bool>{true} : std::nullopt);
    }
    const std::string variable = variables.empty() ? "x" : *variables.begin();
    auto recognized = recognize_rational_polynomial(*difference, variable, context);
    if (!recognized) {
        return Result<std::optional<bool>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    return Result<std::optional<bool>>::success(recognized.value()->is_zero());
}

}
