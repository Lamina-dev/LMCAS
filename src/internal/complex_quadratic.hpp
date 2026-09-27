#pragma once

#include "computation_context.hpp"
#include "symbolic.hpp"

#include <array>
#include <memory>

namespace LMCAS::detail {
struct QuadraticCartesianRoot {
    std::shared_ptr<SymbolicExpr> real;
    std::shared_ptr<SymbolicExpr> imag;
};
using CartesianQuadraticRoots = std::array<QuadraticCartesianRoot, 2>;

Result<CartesianQuadraticRoots> real_quadratic_roots_checked(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    ComputationContext& context);

}
