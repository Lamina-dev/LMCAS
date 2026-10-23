#include "internal/expression_transform.hpp"
#include "matrix_quadratic_form.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/differentiation_visitor.hpp"

#include <optional>
#include <utility>

namespace LMCAS {
namespace {

bool valid_variables(const std::vector<std::string>& variables) {
    if (variables.empty()) return false;
    for (std::size_t index = 0; index < variables.size(); ++index) {
        if (variables[index].empty()) return false;
        for (std::size_t previous = 0; previous < index; ++previous) {
            if (variables[index] == variables[previous]) return false;
        }
    }
    return true;
}

bool independent_of_variables(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::vector<std::string>& variables) {
    if (!expression || !detail::node(expression)) return false;
    for (const auto& variable : variables) {
        if (expression_depends_on_variable(
                detail::node(expression), variable)) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<SymbolicExpr> reconstruct_quadratic_form(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& matrix,
    const std::vector<std::shared_ptr<SymbolicExpr>>& variables) {
    auto reconstructed = SymbolicExpr::number(0);
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        for (std::size_t column = 0; column < matrix.size(); ++column) {
            if (detail::node(matrix[row][column])->is_zero()) continue;
            auto term = SymbolicExpr::multiply(
                matrix[row][column],
                SymbolicExpr::multiply(variables[row], variables[column]));
            reconstructed = SymbolicExpr::add(reconstructed, term);
        }
    }
    return reconstructed;
}

using SymbolicMatrix =
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>;

std::vector<std::shared_ptr<SymbolicExpr>> make_variables(
    const std::vector<std::string>& names) {
    std::vector<std::shared_ptr<SymbolicExpr>> variables;
    variables.reserve(names.size());
    for (const auto& name : names) {
        variables.push_back(SymbolicExpr::variable(name));
    }
    return variables;
}

std::optional<SymbolicMatrix> quadratic_coefficients(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& variables, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    const auto half = SymbolicExpr::number(Rational(1, 2));
    SymbolicMatrix matrix(
        variables.size(),
        std::vector<std::shared_ptr<SymbolicExpr>>(variables.size()));
    for (std::size_t row = 0; row < variables.size(); ++row) {
        auto first_derivative = detail::make_expression_ptr(detail::differentiate_node(detail::node(expanded), variables[row], transform_budget));
        if (!first_derivative) {
            return std::nullopt;
        }
        for (std::size_t column = 0; column < variables.size(); ++column) {
            auto second_derivative =
                detail::make_expression_ptr(detail::differentiate_node(detail::node(first_derivative), variables[column], transform_budget));
            if (!second_derivative) {
                return std::nullopt;
            }
            auto coefficient =
                detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(half, second_derivative)), transform_budget));
            if (!independent_of_variables(coefficient, variables)) {
                return std::nullopt;
            }
            matrix[row][column] = std::move(coefficient);
        }
    }
    return matrix;
}

bool reconstructs_expression(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const SymbolicMatrix& matrix,
    const std::vector<std::shared_ptr<SymbolicExpr>>& variables, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto reconstructed = reconstruct_quadratic_form(matrix, variables);
    auto residual = SymbolicExpr::add(
        expanded,
        SymbolicExpr::multiply(
            SymbolicExpr::number(-1), reconstructed));
    residual = detail::make_expression_ptr(detail::expand_node(detail::node(residual), transform_budget));
    if (residual) {
        residual = detail::make_expression_ptr(detail::normalize_node(detail::node(residual), transform_budget));
    }
    return residual && detail::node(residual) &&
        detail::node(residual)->is_zero();
}

}

std::shared_ptr<SymbolicExpr> quadratic_form_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& vars) {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!expr || !valid_variables(vars)) {
        return nullptr;
    }

    try {
        auto expanded = detail::make_expression_ptr(detail::expand_node(detail::node(expr), transform_budget));
        if (!expanded) {
            return nullptr;
        }
        auto variables = make_variables(vars);
        auto matrix = quadratic_coefficients(expanded, vars, context);
        if (!matrix ||
            !reconstructs_expression(expanded, *matrix, variables, context)) {
            return nullptr;
        }
        return SymbolicExpr::matrix(*matrix);
    } catch (const CasError& error) {
        if (error.code != CasErrc::UnsupportedExpression) throw;
        return nullptr;
    }
}

}
