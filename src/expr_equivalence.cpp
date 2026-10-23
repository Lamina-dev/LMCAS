#include "expr.hpp"

#include "internal/equivalence_engine.hpp"

namespace LMCAS {

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context,
                        const EqvOptions& options) {
    return detail::equivalent_impl(lhs, rhs, context, options);
}


}
