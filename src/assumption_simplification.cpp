#include <stdexcept>
#include "assumption_context.hpp"

#include "internal/assumption_facts.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/expression_transform.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"

namespace LMCAS {
namespace detail {

namespace {
template <typename Transform>
ExpressionResult assumed_transform(const ExprPtr& expression, ComputationContext& context,
                                   const char* operation, const char* name,
                                   Transform&& transform) {
    return checked_transform_expr(
        expression, context, operation, name,
        [&](const SymbolicExpr& value) {
            RewriteBudget budget(context, context.limits().max_recursion_depth,
                                 context.limits().max_ast_nodes, operation);
            std::optional<AssumptionFacts> facts;
            if (context.assumptions()) { facts.emplace(*context.assumptions()); }
            return make_expression_ptr(transform(
                node(value), budget, facts ? static_cast<const FactsQuery&>(*facts) : no_facts()));
        });
}
}

ExpressionResult simplify_expression(const ExprPtr& expression,
                                     ComputationContext& context) {
    return assumed_transform(expression, context, "LMCAS.simplify", "simplify",
        [](const SymbolicNodePtr& value, RewriteBudget& budget, const FactsQuery& facts) {
            return normalize_node(value, budget, facts);
        });
}

ExpressionResult expand_expression(const ExprPtr& expression, ComputationContext& context) {
    return assumed_transform(expression, context, "LMCAS.expand", "expand",
        [](const SymbolicNodePtr& value, RewriteBudget& budget, const FactsQuery& facts) {
            return expand_node(value, budget, facts);
        });
}

ExpressionResult differentiate_expression(const ExprPtr& expression,
                                          const std::string& variable,
                                          ComputationContext& context) {
    return assumed_transform(expression, context, "LMCAS.differentiate", "differentiate",
        [&](const SymbolicNodePtr& value, RewriteBudget& budget, const FactsQuery& facts) {
            return differentiate_node(value, variable, budget, facts);
        });
}

ExpressionResult substitute_expression(const ExprPtr& expression,
                                       const std::string& variable,
                                       const ExprPtr& replacement,
                                       ComputationContext& context) {
    return assumed_transform(expression, context, "LMCAS.substitute", "substitute",
        [&](const SymbolicNodePtr& value, RewriteBudget& budget, const FactsQuery& facts) {
            return substitute_node(value, variable, node(replacement), budget, facts);
        });
}

ExpressionResult substitute_expression(
    const ExprPtr& expression,
    const std::vector<Binding>& replacements,
    ComputationContext& context) {
    return assumed_transform(expression, context, "LMCAS.substitute", "substitute",
        [&](SymbolicNodePtr value, RewriteBudget& budget, const FactsQuery& facts) {
            budget.measure(value);
            for (const auto& replacement : replacements) {
                const auto symbol = std::dynamic_pointer_cast<const VariableNode>(
                    node(replacement.symbol));
                if (!symbol || symbol->is_constant() || !replacement.value ||
                    !node(replacement.value)) {
                    throw CasError{CasErrc::InvalidArgument, "invalid substitution binding",
                                   "LMCAS.substitute"};
                }
                value = substitute_node(value, symbol->name(), node(replacement.value),
                                        budget, facts);
            }
            return value;
        });
}

void AssumptionFacts::collect_names(const AssumptionContext& assumptions,
    std::set<std::string>& names, RewriteBudget& budget) {
    const auto include = [&](const SymbolicExpr& expression) {
        budget.measure(node(expression));
        auto variables = free_variables(node(expression), &budget);
        names.merge(variables);
    };
    for (const auto& scope : assumptions.scope_stack_) {
        auto access = budget.context().consume_steps(1, "assumption.names");
        if (!access) { throw access.error(); }
        for (auto& symbol : scope.properties.get_all_symbols()) {
            access = budget.context().consume_steps(1, "assumption.names");
            if (!access) { throw access.error(); }
            for (const auto& [variable, period] : scope.properties.get_period_decls(symbol)) {
                names.insert(variable);
                include(period);
            }
            for (const auto& declaration : scope.properties.get_monotonicity_decls(symbol)) {
                access = budget.context().consume_steps(1, "assumption.names");
                if (!access) { throw access.error(); }
                names.insert(declaration.variable);
            }
            names.insert(std::move(symbol));
        }
        for (const auto& relation : scope.relations.get_relations()) {
            include(relation.lhs);
            include(relation.rhs);
        }
        for (const auto& conditional : scope.conditionals) {
            include(conditional.condition);
            include(conditional.conclusion);
        }
    }
}

Result<Tribool> AssumptionFacts::is_positive(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.assumptions");
    if (!access) { return Result<Tribool>::failure(access.error()); }
    return engine_->query_positive_checked(expression_from_node(node), context);
}

Result<Tribool> AssumptionFacts::is_negative(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.assumptions");
    if (!access) { return Result<Tribool>::failure(access.error()); }
    return engine_->query_negative_checked(expression_from_node(node), context);
}

Result<Tribool> AssumptionFacts::is_nonnegative(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.assumptions");
    if (!access) { return Result<Tribool>::failure(access.error()); }
    return engine_->query_nonnegative_checked(expression_from_node(node), context);
}

Result<Tribool> AssumptionFacts::is_nonzero(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.assumptions");
    if (!access) { return Result<Tribool>::failure(access.error()); }
    return engine_->query_nonzero_checked(expression_from_node(node), context);
}

Result<Tribool> AssumptionFacts::is_real(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.assumptions");
    if (!access) { return Result<Tribool>::failure(access.error()); }
    return engine_->query_real_checked(expression_from_node(node), context);
}

}


}
