#pragma once

#include "internal/symbolic_ast.hpp"
#include "multivariate_factor.hpp"

namespace LMCAS::detail {

inline std::vector<std::shared_ptr<SymbolicExpr>> multipoly_term_factors(
    const Monomial& monomial, const Rational& coefficient,
    const std::vector<std::string>& variables) {
    std::vector<std::shared_ptr<SymbolicExpr>> factors;
    if (coefficient != Rational(1) || total_degree(monomial) == 0) {
        factors.push_back(SymbolicExpr::number(coefficient));
    }
    for (std::size_t index = 0;
         index < variables.size() && index < monomial.size(); ++index) {
        if (monomial[index] == 0) continue;
        auto variable = SymbolicExpr::variable(variables[index]);
        factors.push_back(monomial[index] == 1
            ? variable
            : SymbolicExpr::power(
                  variable, SymbolicExpr::number(monomial[index])));
    }
    return factors;
}

inline std::shared_ptr<SymbolicExpr> multiply_expression_factors(
    const std::vector<std::shared_ptr<SymbolicExpr>>& factors) {
    auto expression = factors.empty()
        ? SymbolicExpr::number(1)
        : factors.front();
    for (std::size_t index = 1; index < factors.size(); ++index) {
        expression = SymbolicExpr::multiply(expression, factors[index]);
    }
    return expression;
}

inline std::shared_ptr<SymbolicExpr> multipoly_term_to_expression(
    const Monomial& monomial, const Rational& coefficient,
    const std::vector<std::string>& variables) {
    return multiply_expression_factors(
        multipoly_term_factors(monomial, coefficient, variables));
}

}
