#pragma once
/**
 * Shared constructors for multivariate polynomial tests.
 */

#include "test_common.hpp"
#include "multivariate_factor.hpp"

using namespace LMCAS;

inline MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff)
{
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}
