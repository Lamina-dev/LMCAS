#include "expr.hpp"

#include "internal/equivalence_engine.hpp"

namespace LMCAS {

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context,
                        const EqvOptions& options) {
    return equivalent_core(lhs, rhs, context, options);
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context) {
    return equivalent(lhs, rhs, context, EqvOptions{});
}

}
