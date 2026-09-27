#include "internal/fglm_internal.hpp"

namespace LMCAS::fglm_detail {
namespace {
bool is_pure_power(const Monomial& monomial, size_t variable, size_t num_vars) {
    for (size_t index = 0; index < num_vars; ++index) {
        int exponent = index < monomial.size() ? monomial[index] : 0;
        if (index == variable) {
            if (exponent == 0) return false;
        } else if (exponent != 0) {
            return false;
        }
    }
    return true;
}
}

bool is_standard(const Monomial& monomial, const std::vector<Monomial>& leading) {
    for (const auto& divisor : leading) {
        if (divides_monomial(divisor, monomial)) return false;
    }
    return true;
}

bool has_pure_generators(const std::vector<Monomial>& leading, size_t num_vars) {
    for (size_t variable = 0; variable < num_vars; ++variable) {
        bool found = false;
        for (const auto& monomial : leading) {
            if (is_pure_power(monomial, variable, num_vars)) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

std::vector<Monomial> leading_monomials(const std::vector<FGLMPoly>& basis) {
    std::vector<Monomial> leading;
    for (const auto& polynomial : basis) {
        if (!polynomial.is_zero()) leading.push_back(polynomial.lead_monomial());
    }
    return leading;
}

void generate_degree(Monomial& current, size_t variable, int remaining,
                     std::vector<Monomial>& monomials) {
    if (variable == current.size() - 1) {
        current[variable] = remaining;
        monomials.push_back(current);
        current[variable] = 0;
        return;
    }
    for (int exponent = 0; exponent <= remaining; ++exponent) {
        current[variable] = exponent;
        generate_degree(current, variable + 1, remaining - exponent, monomials);
    }
    current[variable] = 0;
}
}

namespace LMCAS {
bool is_zero_dimensional(const std::vector<FGLMPoly>& basis, size_t num_vars) {
    if (basis.empty()) return false;
    for (size_t variable = 0; variable < num_vars; ++variable) {
        bool found = false;
        for (const auto& polynomial : basis) {
            if (polynomial.is_zero()) continue;
            if (fglm_detail::is_pure_power(polynomial.terms.front().first,
                                           variable, num_vars)) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

int quotient_dimension(const std::vector<FGLMPoly>& basis, size_t num_vars,
                       int max_degree) {
    if (basis.empty()) return -1;
    auto leading = fglm_detail::leading_monomials(basis);
    int count = 0;
    const int limit = 1000;
    for (int degree = 0; degree <= max_degree && count <= limit; ++degree) {
        std::vector<Monomial> monomials;
        Monomial current(num_vars, 0);
        fglm_detail::generate_degree(current, 0, degree, monomials);
        bool any_standard = false;
        for (const auto& monomial : monomials) {
            if (!fglm_detail::is_standard(monomial, leading)) continue;
            any_standard = true;
            if (++count > limit) return count;
        }
        if (!monomials.empty() && !any_standard && degree > 0) break;
    }
    return count;
}
}
