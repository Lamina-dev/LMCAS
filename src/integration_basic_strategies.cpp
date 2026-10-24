#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

Result<std::optional<SymbolicExpr>> detail::integrate_power_primitive(
    const SymbolicExpr& base, const SymbolicExpr& exponent,
    const std::string& variable, ComputationContext& context) {
    if (depends_on_integration_variable(exponent, variable)) {
        return std::optional<SymbolicExpr>{};
    }
    RewriteBudget budget(context, context.limits().max_recursion_depth,
                         context.limits().max_ast_nodes, "integrate.power");
    auto exact = exact_rational_value(node(exponent));
    if (exact && *exact == Rational(-1)) {
        return std::optional<SymbolicExpr>{*SymbolicExpr::ln(make_expression_ptr(base))};
    }
    auto next = exact ? SymbolicExpr::number(*exact + Rational(1)) :
        make_expression_ptr(normalize_node(node(SymbolicExpr::add(
            make_expression_ptr(exponent), SymbolicExpr::number(1))), budget));
    if (!exact) {
        std::optional<AssumptionFacts> facts;
        if (context.assumptions()) { facts.emplace(*context.assumptions()); }
        auto nonzero = query_nonzero_value(node(next),
            facts ? static_cast<const FactsQuery&>(*facts) : no_facts(), Domain::Real, context);
        if (!nonzero) {
            return Result<std::optional<SymbolicExpr>>::failure(nonzero.error());
        }
        if (nonzero.value() != Tribool::True) { return std::optional<SymbolicExpr>{}; }
    }
    return std::optional<SymbolicExpr>{*SymbolicExpr::divide(
        SymbolicExpr::power(make_expression_ptr(base), next), next)};
}

Result<std::optional<detail::IntegralParts>> Integrator::try_table_lookup(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int) {
    detail::RewriteBudget budget(computation, computation.limits().max_recursion_depth,
                                 computation.limits().max_ast_nodes, "transform");
    auto primitive = lookup_integral(entries_, expr, var, budget,
        [](std::shared_ptr<SymbolicExpr> candidate) { return candidate; });
    if (!primitive) { return std::optional<detail::IntegralParts>{}; }
    return std::optional<detail::IntegralParts>{{*primitive, zero_}};
}

Result<std::optional<detail::IntegralParts>> Integrator::try_power_rule(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int) {
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(detail::node(expr))) {
        if (!variable->is_constant() && variable->name() == var) {
            return std::optional<detail::IntegralParts>{{*SymbolicExpr::multiply(
                SymbolicExpr::power(detail::make_expression_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2)), zero_}};
        }
    }
    auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(expr));
    if (!power) { return std::optional<detail::IntegralParts>{}; }
    auto variable = std::dynamic_pointer_cast<const VariableNode>(power->base());
    if (!variable || variable->is_constant() || variable->name() != var) {
        return std::optional<detail::IntegralParts>{};
    }
    auto primitive = detail::integrate_power_primitive(detail::expression_from_node(power->base()),
        detail::expression_from_node(power->exponent()), var, context);
    if (!primitive) {
        return Result<std::optional<detail::IntegralParts>>::failure(primitive.error());
    }
    if (!primitive.value()) { return std::optional<detail::IntegralParts>{}; }
    return std::optional<detail::IntegralParts>{{std::move(*primitive.value()), zero_}};
}

}
