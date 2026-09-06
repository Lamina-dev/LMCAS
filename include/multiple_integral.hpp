/** @file multiple_integral.hpp */ #pragma once
#include "integrator.hpp"
#include <memory>
#include <string>
#include <vector>

namespace LMCAS {

struct IntegrationStep {
    std::string variable;
    std::shared_ptr<SymbolicExpr> lower;
    std::shared_ptr<SymbolicExpr> upper;
};

LMCAS_API Result<SymbolicExpr> integrate_multiple_checked(
    const SymbolicExpr& integrand,
    const std::vector<IntegrationStep>& steps,
    Integrator& integrator,
    ComputationContext& context);

} // namespace LMCAS
