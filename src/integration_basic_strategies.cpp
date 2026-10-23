#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

Result<std::shared_ptr<SymbolicExpr>> TableLookupStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int) try {
    detail::RewriteBudget transform_budget(computation, (computation).limits().max_recursion_depth,
                                            (computation).limits().max_ast_nodes, "transform");
    return lookup_integral(ctx.table(), expr, var, ctx.original_integrand(expr), var, transform_budget,
        [](std::shared_ptr<SymbolicExpr> candidate) { return candidate; });
} catch (const CasError& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, "allocation failed", "try_integrate_raw");
} catch (const std::length_error& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, error.what(), "try_integrate_raw");
} catch (const std::exception& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::InternalInvariant, error.what(), "try_integrate_raw");
}

Result<std::shared_ptr<SymbolicExpr>> PowerRuleStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& context, int) try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");

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
                if (contains_inexact_number(detail::node(exp_expr))) { return nullptr; }
                auto exponent = detail::exact_rational_value(detail::node(exp_expr));
                std::shared_ptr<SymbolicExpr> n_plus_1;
                if (exponent) {
                    const Rational next = *exponent + Rational(1);
                    if (next == Rational(0)) {
                        return SymbolicExpr::ln(detail::make_expression_ptr(base));
                    }
                    n_plus_1 = SymbolicExpr::number(next);
                } else {
                    n_plus_1 = detail::make_expression_ptr(detail::normalize_node(
                        detail::node(SymbolicExpr::add(detail::make_expression_ptr(exp_expr),
                            SymbolicExpr::number(1))), transform_budget));
                    std::optional<detail::AssumptionFacts> facts;
                    if (context.assumptions()) { facts.emplace(*context.assumptions()); }
                    auto nonzero = detail::query_nonzero_value(detail::node(n_plus_1),
                        facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts(), Domain::Real, context);
                    if (!nonzero) {
                        return Result<std::shared_ptr<SymbolicExpr>>::failure(nonzero.error());
                    }
                    if (nonzero.value() != Tribool::True) { return nullptr; }
                }
                return SymbolicExpr::divide(
                    SymbolicExpr::power(detail::make_expression_ptr(base), n_plus_1),
                    n_plus_1);
            }
        }
    }
    return nullptr;
} catch (const CasError& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, "allocation failed", "try_integrate_raw");
} catch (const std::length_error& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, error.what(), "try_integrate_raw");
} catch (const std::exception& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::InternalInvariant, error.what(), "try_integrate_raw");
}

}
