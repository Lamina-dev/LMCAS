#include "internal/fglm_internal.hpp"
#include <map>
#include <stdexcept>
#include <string>

namespace LMCAS {
namespace {
using namespace fglm_detail;

void generate_standard_degree(Monomial& current, size_t variable, int remaining,
                              const std::vector<Monomial>& leading,
                              std::vector<Monomial>& standard) {
    if (variable == current.size() - 1) {
        current[variable] = remaining;
        if (is_standard(current, leading)) standard.push_back(current);
        current[variable] = 0;
        return;
    }
    for (int exponent = 0; exponent <= remaining; ++exponent) {
        current[variable] = exponent;
        generate_standard_degree(current, variable + 1, remaining - exponent,
                                 leading, standard);
    }
    current[variable] = 0;
}

std::vector<Monomial> standard_basis(const std::vector<FGLMPoly>& source_basis,
                                     size_t num_vars, int dimension) {
    auto leading = leading_monomials(source_basis);
    std::vector<Monomial> standard;
    for (int degree = 0;
         degree <= dimension + 2 && static_cast<int>(standard.size()) < dimension;
         ++degree) {
        Monomial current(num_vars, 0);
        generate_standard_degree(current, 0, degree, leading, standard);
    }
    return standard;
}

std::vector<Rational> coordinates(const FGLMPoly& polynomial,
                                  const std::map<Monomial, size_t>& index) {
    std::vector<Rational> vector(index.size(), Rational(0));
    for (const auto& [monomial, coefficient] : polynomial.terms) {
        auto found = index.find(monomial);
        if (found != index.end()) vector[found->second] = coefficient;
    }
    return vector;
}

FGLMPoly dependent_polynomial(const Monomial& monomial,
                              const std::vector<Rational>& combination,
                              const std::vector<Monomial>& basis,
                              const MonomialOrder& order, size_t num_vars) {
    FGLMPoly result(num_vars);
    result.add_term(monomial, Rational(1));
    for (size_t index = 0; index < combination.size() && index < basis.size(); ++index) {
        if (!combination[index].is_zero()) {
            result.add_term(basis[index], -combination[index]);
        }
    }
    result.sort_terms(order);
    result.normalize();
    if (!result.is_zero()) {
        Rational leading = result.lead_coeff();
        if (!(leading == Rational(1))) {
            for (auto& term : result.terms) term.second = term.second / leading;
        }
    }
    return result;
}

int checked_dimension(const std::vector<FGLMPoly>& source_basis, size_t num_vars) {
    if (!is_zero_dimensional(source_basis, num_vars)) {
        throw std::runtime_error(
            "FGLM: ideal is not zero-dimensional. "
            "FGLM only works for zero-dimensional ideals.");
    }
    int dimension = quotient_dimension(source_basis, num_vars);
    if (dimension <= 0 || dimension > 1000) {
        throw std::runtime_error(
            "FGLM: quotient dimension is too large or could not be determined. "
            "dim = " + std::to_string(dimension));
    }
    return dimension;
}
}

std::vector<FGLMPoly> fglm_convert(const std::vector<FGLMPoly>& source_basis,
                                  const MonomialOrder& source_order,
                                  const MonomialOrder& target_order,
                                  size_t num_vars) {
    int dimension = checked_dimension(source_basis, num_vars);
    std::vector<FGLMPoly> target_basis;
    std::vector<Monomial> basis_monomials;
    auto standard = standard_basis(source_basis, num_vars, dimension);
    std::map<Monomial, size_t> monomial_index;
    for (size_t index = 0; index < standard.size(); ++index) {
        monomial_index[standard[index]] = index;
    }
    fglm_detail::GaussianEliminator gauss;
    fglm_detail::MonomialEnumerator enumerator(num_vars, target_order, dimension + 5);
    std::vector<Monomial> target_leading;
    const int max_iterations = dimension * 10 + 50;
    int iterations = 0;
    Monomial monomial;
    while (iterations < max_iterations && enumerator.next(monomial)) {
        ++iterations;
        if (!fglm_detail::is_standard(monomial, target_leading)) continue;
        auto polynomial = FGLMPoly::from_monomial(monomial, num_vars);
        auto normal = normal_form(polynomial, source_basis, source_order);
        auto vector = coordinates(normal, monomial_index);
        std::vector<Rational> combination;
        if (gauss.add_vector(vector, combination)) {
            basis_monomials.push_back(monomial);
        } else {
            auto dependent = dependent_polynomial(monomial, combination,
                                                   basis_monomials, target_order,
                                                   num_vars);
            if (!dependent.is_zero()) {
                target_leading.push_back(dependent.lead_monomial());
                target_basis.push_back(std::move(dependent));
            }
        }
        if (static_cast<int>(basis_monomials.size()) == dimension &&
            fglm_detail::has_pure_generators(target_leading, num_vars)) {
            return target_basis;
        }
    }
    throw std::runtime_error(
        "FGLM: enumeration ended before the target basis was complete");
}

std::vector<FGLMPoly> grevlex_to_lex(const std::vector<FGLMPoly>& grevlex_basis,
                                    size_t num_vars) {
    return fglm_convert(grevlex_basis, MonomialOrder::grevlex(),
                        MonomialOrder::lex(), num_vars);
}
}
