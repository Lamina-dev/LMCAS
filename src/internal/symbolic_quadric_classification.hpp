#pragma once

#include "symbolic_geometry_quadric.hpp"

#include <string>

namespace LMCAS::detail {

VectorStringResult classify_quadric_impl(
    const SurfaceSymbolic& surf,
    ComputationContext& context,
    const std::string& operation);

}
