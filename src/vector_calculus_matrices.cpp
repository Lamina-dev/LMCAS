#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_matrices.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
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
        return vector_calculus_wrap_expr(jacobian(functions, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusExprResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return jacobian_checked(functions, vars, context);
}

std::shared_ptr<SymbolicExpr> jacobian(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    size_t m = functions.size();
    size_t n = vars.size();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid;
    grid.reserve(m);

    for (size_t i = 0; i < m; ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row;
        row.reserve(n);
        for (size_t j = 0; j < n; ++j) {
            auto partial = functions[i]->differentiate(vars[j]);
            partial = partial->simplify();
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
        return vector_calculus_wrap_expr(hessian(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusExprResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return hessian_checked(f, vars, context);
}

std::shared_ptr<SymbolicExpr> hessian(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    size_t n = vars.size();

    std::vector<std::shared_ptr<SymbolicExpr>> first_partials;
    first_partials.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        first_partials.push_back(f->differentiate(vars[i]));
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            auto second_partial = first_partials[i]->differentiate(vars[j]);
            second_partial = second_partial->simplify();
            grid[i][j] = second_partial;
            if (i != j) {
                grid[j][i] = second_partial;
            }
        }
    }

    return SymbolicExpr::matrix(grid);
}

}
