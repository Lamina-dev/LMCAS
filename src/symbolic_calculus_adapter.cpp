#include "internal/expression_transform.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

std::shared_ptr<SymbolicExpr> SymbolicExpr::series(
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point,
    int order,
    const LMCAS::AssumptionContext*) const {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!impl_->root || !point) return nullptr;
    if (order < 0) return nullptr;

    auto variable_expr = SymbolicExpr::variable(variable);
    auto delta = SymbolicExpr::add(
        variable_expr,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), point));
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    auto derivative = LMCAS::detail::make_expression_ptr(impl_->root->clone());

    for (int n = 0;; ++n) {
        if (n > 0) {
            derivative = detail::make_expression_ptr(detail::differentiate_node(detail::node(derivative), variable, transform_budget));
            if (!derivative) return nullptr;
            derivative = detail::make_expression_ptr(detail::normalize_node(detail::node(derivative), transform_budget));
        }

        auto coefficient = detail::make_expression_ptr(detail::substitute_node(detail::node(derivative), variable, detail::node(point), transform_budget));
        if (!coefficient) return nullptr;
        coefficient = detail::make_expression_ptr(detail::normalize_node(detail::node(coefficient), transform_budget));

        auto term = coefficient;
        if (n > 0) {
            term = SymbolicExpr::multiply(
                term,
                SymbolicExpr::power(delta, SymbolicExpr::number(n)));
        }
        if (n > 1) {
            auto factorial = BigInt::factorial(BigInt(n));
            term = SymbolicExpr::multiply(
                term,
                SymbolicExpr::number(Rational(BigInt(1), factorial)));
        }
        terms.push_back(LMCAS::detail::node(term));
        if (n == order) break;
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    return detail::make_expression_ptr(detail::normalize_node(LMCAS::detail::make_node<AddNode>(terms), transform_budget));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point) {
    if (!expression || !point) return nullptr;
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<LimitNode>(
            LMCAS::detail::node(expression), variable,
            LMCAS::detail::node(point), LimitDirection::Both));
}

} // namespace LMCAS
