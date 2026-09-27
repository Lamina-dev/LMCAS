#include "expr.hpp"

#include "internal/equivalence_engine.hpp"

namespace LMCAS {

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context,
                        const EqvOptions& options) {
    auto checked = equivalent_core(lhs, rhs, context, options);
    if (checked) return checked;

    const auto code = checked.error().code;
    if (code == CasErrc::ResourceLimit ||
        code == CasErrc::Inconclusive ||
        code == CasErrc::UnsupportedExpression) {
        (void)context.add_diagnostic(
            Diagnostic{DiagnosticSeverity::Warning,
                       kEquivalentOperation,
                       checked.error().message});
        return Result<bool>::success(false);
    }
    return Result<bool>::failure(checked.error());
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context) {
    return equivalent(lhs, rhs, context, EqvOptions{});
}

}
