#pragma once

#include "polynomial.hpp"

namespace lamina {

std::vector<Polynomial<Rational>> factor_univariate_bridge(
    const Polynomial<Rational>& polynomial);

} // namespace lamina
