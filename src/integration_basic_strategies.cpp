#include "internal/integration_support.hpp"

namespace LMCAS {

Result<std::shared_ptr<SymbolicExpr>> TableLookupStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int) {
    (void)computation;

    auto all_entries = ctx.table().get_all_sorted();

    for (const auto* entry : all_entries) {
        MatchMap bindings;
        if (Matcher::match(entry->pattern, expr, entry->wildcards, bindings)) {

            // The matcher allows partial matches on commutative operations
            // (Add/Multiply): unmatched operands are bound to __Add_REST__ /
            // __Mul_REST__ for use by rewrite rules. For integration table
            // lookup that behavior is wrong - a pattern like (1 + _u^2)^-1
            // would otherwise match (1 + x^2 + x + x^3)^-1 with _u=x and the
            // rest x + x^3 silently dropped, producing arctan(x) for an
            // integrand that has nothing to do with arctan. Require an exact
            // match (no leftover) so the result genuinely equals the
            // integrand under this binding.
            if (bindings.find("__Add_REST__") != bindings.end() ||
                bindings.find("__Mul_REST__") != bindings.end()) {
                continue;
            }
            const bool complete_binding = std::all_of(
                entry->wildcards.begin(), entry->wildcards.end(),
                [&](const std::string& wildcard_name) {
                    return bindings.find(wildcard_name) != bindings.end();
                });
            if (!complete_binding) { continue; }

            if (entry->condition && !entry->condition(bindings, var)) {
                continue;
            }

            SymbolicExpr result = Matcher::replace(entry->result, bindings, false);
            auto simplified = result.simplify();
            if (simplified && !contains_unevaluated_integral(LMCAS::detail::node(simplified))) {
                return simplified;
            }
        }
    }
    return nullptr;
}

Result<std::shared_ptr<SymbolicExpr>> PowerRuleStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext&, int) {

    if (auto v_node = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        if (!v_node->is_constant() && v_node->name() == var) {
            return SymbolicExpr::multiply(
                SymbolicExpr::power(detail::make_expression_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2));
        }
    }

    if (auto p_node = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        auto base = LMCAS::detail::expression_from_node(p_node->base());
        auto exp_expr = LMCAS::detail::expression_from_node(p_node->exponent());
        if (auto b_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(base))) {
            if (!b_var->is_constant() && b_var->name() == var &&
                !depends_on_integration_variable(exp_expr, var)) {
                auto n_plus_1 = SymbolicExpr::add(detail::make_expression_ptr(exp_expr), SymbolicExpr::number(1))->simplify();
                if (n_plus_1->is_zero()) {

                    return SymbolicExpr::ln(detail::make_expression_ptr(base));
                }
                return SymbolicExpr::divide(
                    SymbolicExpr::power(detail::make_expression_ptr(base), n_plus_1),
                    n_plus_1);
            }
        }
    }
    return nullptr;
}

}
