#include "internal/integration_support.hpp"

namespace LMCAS {

static bool validate_integration_steps(const std::vector<IntegrationStep>& steps) {
    if (steps.empty()) return false;

    for (size_t i = 0; i < steps.size(); ++i) {
        const auto& s = steps[i];
        if (s.variable.empty()) return false;

        // Bounds must be either both null (indefinite) or both set (definite).
        const bool has_lower = static_cast<bool>(s.lower);
        const bool has_upper = static_cast<bool>(s.upper);
        if (has_lower != has_upper) return false;

        // Reject duplicate variables.
        for (size_t j = 0; j < i; ++j) {
            if (steps[j].variable == s.variable) return false;
        }
    }
    return true;
}

Result<SymbolicExpr> integrate_multiple_checked(
    const SymbolicExpr& integrand,
    const std::vector<IntegrationStep>& steps,
    Integrator& integrator,
    ComputationContext& context) {

    if (!validate_integration_steps(steps)) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "invalid multiple integration steps",
            "integrate.multiple");
    }

    auto current = LMCAS::detail::make_expression_ptr(integrand);

    for (const auto& step : steps) {
        if (!current) {
            return Result<SymbolicExpr>::failure(
                CasErrc::InternalInvariant, "multiple integration lost its expression",
                "integrate.multiple");
        }


        const bool definite = static_cast<bool>(step.lower) && static_cast<bool>(step.upper);

        if (!definite) {
            auto result = integrator.integrate_checked(*current, step.variable, context);
            if (!result) return result;
            current = LMCAS::detail::make_expression_ptr(result.value());
        } else {
            if (!Integrator::depends_on(*current, step.variable)) {
                auto diff = sym_sub(*step.upper, *step.lower);
                auto product = SymbolicExpr::multiply(current, diff);
                auto simp = product->simplify();
                current = simp ? simp : product;
            } else {
                auto result = integrator.integrate_def_checked(
                    *current, step.variable, *step.lower, *step.upper, context);
                if (!result) return result;
                auto res_ptr = LMCAS::detail::make_expression_ptr(result.value());
                auto simp = res_ptr->simplify();
                current = simp ? simp : res_ptr;
            }
        }
    }

    return Result<SymbolicExpr>::success(*current);
}

} // namespace LMCAS
