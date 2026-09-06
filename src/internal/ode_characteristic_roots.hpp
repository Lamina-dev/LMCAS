#pragma once

#include "ode_polynomial_utils.hpp"
#include "result.hpp"

#include <string>
#include <vector>

namespace LMCAS::ode_root_detail {

Result<std::vector<CharRoot>> find_characteristic_roots(
    const std::vector<double>& coefficients,
    const std::string& operation);

} // namespace LMCAS::ode_root_detail
