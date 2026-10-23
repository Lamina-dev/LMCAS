#pragma once

#include "computation_context.hpp"
#include "equivalence_options.hpp"
#include "symbolic.hpp"

namespace LMCAS {

inline constexpr const char* kEquivalentOperation = "LMCAS.equivalent";
inline constexpr const char* kEquivalentProfileOperation =
    "LMCAS.equivalent.profile";

LMCAS_API bool structurally_equal(const SymbolicExpr& lhs,
                                 const SymbolicExpr& rhs);
namespace detail {
Result<bool> equivalent_impl(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context,
                             const EqvOptions& options);
}

}
