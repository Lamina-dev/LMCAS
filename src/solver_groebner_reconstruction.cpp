#include "internal/solver_groebner_internal.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::groebner_detail {
namespace {
std::shared_ptr<const SymbolicNode> variable_node(size_t index, const PolyContext& context) {
    if (index >= context.num_original_vars) {
        auto found = context.aux_to_node.find(index);
        if (found != context.aux_to_node.end()) return found->second;
    }
    return SymbolicFactory::create_variable(context.ext_vars[index]);
}

std::shared_ptr<const SymbolicNode> term_node(const Monomial& monomial,
                                             const Rational& coefficient,
                                             const PolyContext& context) {
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (coefficient.get_denominator() == BigInt(1)) {
        factors.push_back(SymbolicFactory::create_number(coefficient.get_numerator()));
    } else {
        factors.push_back(SymbolicFactory::create_number(coefficient));
    }
    for (size_t index = 0; index < monomial.size() && index < context.ext_vars.size(); ++index) {
        if (monomial[index] <= 0) continue;
        auto variable = variable_node(index, context);
        if (monomial[index] == 1) {
            factors.push_back(variable);
        } else {
            factors.push_back(SymbolicFactory::create_power(
                variable, SymbolicFactory::create_number(BigInt(monomial[index]))));
        }
    }
    if (factors.size() == 1) return factors[0];
    return SymbolicFactory::create_multiply(factors);
}
}

SymbolicExpr from_poly_ext(
    const Poly& polynomial, const PolyContext& context) {
    if (polynomial.terms.empty()) {
        return LMCAS::detail::expression_from_node(SymbolicFactory::create_number(BigInt(0)));
    }
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    for (const auto& [monomial, coefficient] : polynomial.terms) {
        terms.push_back(term_node(monomial, coefficient, context));
    }
    if (terms.size() == 1) return LMCAS::detail::expression_from_node(terms[0]);
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(terms));
}
}
