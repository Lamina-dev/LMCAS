#include "internal/expression_transform.hpp"
#include "symbolic_implicit_diff.hpp"

namespace LMCAS {

std::shared_ptr<SymbolicExpr> implicit_diff(
    std::shared_ptr<SymbolicExpr> expression,
    const std::string& independent_variable,
    const std::string& dependent_variable) {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!expression || independent_variable.empty() || dependent_variable.empty()) {
        return nullptr;
    }

    auto numerator = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        detail::make_expression_ptr(detail::differentiate_node(detail::node(expression), independent_variable, transform_budget)));
    auto denominator = detail::make_expression_ptr(detail::differentiate_node(detail::node(expression), dependent_variable, transform_budget));
    return SymbolicExpr::divide(numerator, denominator);
}

} // namespace LMCAS
