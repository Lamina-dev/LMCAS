#include "internal/expression_transform.hpp"
#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_matrices.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "jacobian");
    const std::string operation = "jacobian";
    auto valid = vector_calculus_validate_functions_vars(functions, vars,
                                                        context, operation);
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(functions.size() * vars.size(), operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(jacobian_impl(functions, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return jacobian_checked(functions, vars, context);
}

std::shared_ptr<SymbolicExpr> vector_calculus_detail::jacobian_impl(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
{
    size_t m = functions.size();
    size_t n = vars.size();
    std::size_t nodes = 1;
    transform_budget.require_nodes(nodes);

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid;
    grid.reserve(m);

    for (size_t i = 0; i < m; ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row;
        row.reserve(n);
        for (size_t j = 0; j < n; ++j) {
            auto partial = detail::make_expression_ptr(detail::differentiate_node(detail::node(functions[i]), vars[j], transform_budget));
            partial = detail::make_expression_ptr(detail::normalize_node(detail::node(partial), transform_budget));
            nodes = transform_budget.append_size(
                nodes, transform_budget.measure(detail::node(partial)));
            row.push_back(partial);
        }
        grid.push_back(std::move(row));
    }

    return SymbolicExpr::matrix(grid);
}


VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "hessian");
    const std::string operation = "hessian";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(vars.size() * vars.size(), operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(hessian_impl(f, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return hessian_checked(f, vars, context);
}

std::shared_ptr<SymbolicExpr> vector_calculus_detail::hessian_impl(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
{
    size_t n = vars.size();
    std::size_t nodes = 1;
    transform_budget.require_nodes(nodes);

    std::vector<std::shared_ptr<SymbolicExpr>> first_partials;
    first_partials.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        first_partials.push_back(detail::make_expression_ptr(detail::differentiate_node(detail::node(f), vars[i], transform_budget)));
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            auto second_partial = detail::make_expression_ptr(detail::differentiate_node(detail::node(first_partials[i]), vars[j], transform_budget));
            second_partial = detail::make_expression_ptr(detail::normalize_node(detail::node(second_partial), transform_budget));
            const auto partial_nodes = transform_budget.measure(detail::node(second_partial));
            nodes = transform_budget.append_size(nodes, partial_nodes);
            if (i != j) nodes = transform_budget.append_size(nodes, partial_nodes);
            grid[i][j] = second_partial;
            if (i != j) {
                grid[j][i] = second_partial;
            }
        }
    }

    return SymbolicExpr::matrix(grid);
}


std::shared_ptr<SymbolicExpr> jacobian(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "jacobian");
    return jacobian_impl(functions, vars, transform_budget);
}

std::shared_ptr<SymbolicExpr> hessian(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "hessian");
    return hessian_impl(f, vars, transform_budget);
}

}
