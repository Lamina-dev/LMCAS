#pragma once

#include "../include/symbolic_vector_geometry.hpp"

#include <string>

namespace LMCAS::detail {

VectorStringResult classify_quadric_impl(
    const SurfaceSymbolic& surf,
    ComputationContext& context,
    const std::string& operation);

} // namespace LMCAS::detail
