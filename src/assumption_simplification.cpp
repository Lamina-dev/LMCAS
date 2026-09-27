#include "assumption_context.hpp"

#include "internal/assumption_facts.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/expression_transform.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"

namespace LMCAS {
namespace detail {

ExpressionResult simplify_expression(const ExprPtr& expression,
                                     ComputationContext& context) {
    return checked_transform_expr(
        expression, context, "LMCAS.simplify", "simplify",
        [&context](const SymbolicExpr& value) {
            RewriteBudget budget(context, context.limits().max_recursion_depth,
                                 context.limits().max_ast_nodes, "LMCAS.simplify");
            budget.measure(detail::node(value));
            std::optional<AssumptionFacts> facts;
            if (context.assumptions()) { facts.emplace(*context.assumptions()); }
            NormalizationVisitor visitor(context, facts ? *facts : no_facts(), Domain::Real, &budget);
            detail::node(value)->accept(visitor);
            auto result = visitor.get_result();
            budget.measure(result);
            return detail::make_expression_ptr(std::move(result));
        });
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

std::shared_ptr<SymbolicExpr> AssumptionContext::simplify(
    const SymbolicExpr& expression) const {
    const auto& root = detail::node(expression);
    if (!root) return nullptr;
    const detail::AssumptionFacts facts(*this);
    NormalizationVisitor visitor(facts);
    root->accept(visitor);
    return detail::make_expression_ptr(visitor.get_result());
}

}
