#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "result.hpp"

#include <string>
#include <utility>
#include <vector>

namespace LMCAS::detail {

using RationalInterval = std::pair<Rational, Rational>;
inline Rational strict_cauchy_root_bound(
    const Polynomial<Rational>& polynomial) {
    Rational maximum(0);
    const Rational leading = polynomial.lead_coeff().abs();
    for (int degree = 0; degree < polynomial.degree(); ++degree) {
        const Rational ratio =
            polynomial.coeffs[static_cast<std::size_t>(degree)].abs() / leading;
        if (maximum < ratio) maximum = ratio;
    }
    return Rational(2) + maximum;
}


Result<std::vector<RationalInterval>> isolate_real_roots_exact(
    const Polynomial<Rational>& polynomial,
    ComputationContext& context,
    const std::string& operation);

Result<std::size_t> count_real_roots_exact(
    const Polynomial<Rational>& polynomial,
    const Rational& lower,
    const Rational& upper,
    ComputationContext& context,
    const std::string& operation);

}
