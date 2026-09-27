#include "solver.hpp"
#include "internal/solver_groebner_internal.hpp"

namespace LMCAS {
namespace {
using namespace groebner_detail;

std::vector<Poly> convert_polynomials(const std::vector<SymbolicExpr>& expressions,
                                      PolyContext& context) {
    std::vector<Poly> polynomials;
    for (const auto& expression : expressions) {
        Poly polynomial = to_poly(expression, context);
        if (!polynomial.is_zero()) polynomials.push_back(polynomial);
    }
    return polynomials;
}

void extend_basis(std::vector<Poly>& basis, size_t num_vars) {
    for (auto& polynomial : basis) extend_variables(polynomial, num_vars);
}

std::vector<SymbolicExpr> reconstruct_basis(
    const std::vector<Poly>& basis, const PolyContext& context) {
    std::vector<SymbolicExpr> expressions;
    for (const auto& polynomial : basis) {
        expressions.push_back(from_poly_ext(polynomial, context));
    }
    return expressions;
}
}

std::vector<SymbolicExpr> Solver::groebner_basis(
    const std::vector<SymbolicExpr>& polynomials,
    const std::vector<std::string>& variables) {
    PolyContext context(variables);
    auto basis = convert_polynomials(polynomials, context);
    extend_basis(basis, context.ext_vars.size());
    buchberger(basis, context.ext_vars.size());
    return reconstruct_basis(basis, context);
}

std::vector<SymbolicExpr> Solver::reduced_groebner_basis(
    const std::vector<SymbolicExpr>& polynomials,
    const std::vector<std::string>& variables) {
    auto expressions = groebner_basis(polynomials, variables);
    if (expressions.empty()) return {};
    PolyContext context(variables);
    auto basis = convert_polynomials(expressions, context);
    extend_basis(basis, context.ext_vars.size());
    reduce_basis(basis);
    return reconstruct_basis(basis, context);
}

bool Solver::ideal_membership(const SymbolicExpr& polynomial,
                              const std::vector<SymbolicExpr>& basis,
                              const std::vector<std::string>& variables) {
    PolyContext context(variables);
    auto converted_basis = convert_polynomials(basis, context);
    auto converted = to_poly(polynomial, context);
    extend_basis(converted_basis, context.ext_vars.size());
    extend_variables(converted, context.ext_vars.size());
    return reduce(converted, converted_basis).is_zero();
}

std::vector<SymbolicExpr> Solver::elimination_ideal(
    const std::vector<SymbolicExpr>& basis,
    const std::vector<std::string>& variables,
    int elim_count)
{
    if (elim_count <= 0) {
        return basis;
    }
    if (elim_count >= static_cast<int>(variables.size())) {
        return {};
    }

    PolyContext ctx(variables);
    std::vector<SymbolicExpr> result;

    for (const auto& expr : basis) {
        Poly p = to_poly(expr, ctx);
        if (p.is_zero()) continue;

        bool involves_eliminated = false;
        for (const auto& [m, c] : p.terms) {
            for (int i = 0; i < elim_count && i < static_cast<int>(m.size()); ++i) {
                if (m[i] != 0) {
                    involves_eliminated = true;
                    break;
                }
            }
            if (involves_eliminated) break;
        }

        if (!involves_eliminated) {
            result.push_back(expr);
        }
    }

    return result;
}

}
