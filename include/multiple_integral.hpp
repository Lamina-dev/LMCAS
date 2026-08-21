/** @file multiple_integral.hpp */ #pragma once
#include "integrator.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lamina {

struct IntegrationStep {
    std::string variable;
    std::shared_ptr<SymbolicExpr> lower;
    std::shared_ptr<SymbolicExpr> upper;
};

LAMINA_API Result<SymbolicExpr> integrate_multiple_checked(
    const SymbolicExpr& integrand,
    const std::vector<IntegrationStep>& steps,
    Integrator& integrator,
    ComputationContext& context);

} // namespace lamina
