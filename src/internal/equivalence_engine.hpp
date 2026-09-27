#pragma once

#include "computation_context.hpp"
#include "equivalence_options.hpp"
#include "symbolic.hpp"

namespace LMCAS {

inline constexpr const char* kEquivalentOperation = "LMCAS.equivalent_core";
inline constexpr const char* kEquivalentProfileOperation =
    "LMCAS.equivalent_core.profile";

LMCAS_API bool structurally_equal(const SymbolicExpr& lhs,
                                 const SymbolicExpr& rhs);
LMCAS_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                      const SymbolicExpr& rhs,
                                      ComputationContext& context);
LMCAS_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                      const SymbolicExpr& rhs,
                                      ComputationContext& context,
                                      const EqvOptions& options);

}
