#include "internal/expression_transform.hpp"
#include "internal/mixed_transcendental_support.hpp"
#include "internal/numeric_probe.hpp"

#include <cmath>
#include <limits>

namespace LMCAS::detail {

lmmc_real_t mixed_evaluate_at(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x, RewriteBudget& budget) {
    auto substituted = make_expression_ptr(substitute_node(
        node(expr), var, node(SymbolicExpr::number(static_cast<double>(x))), budget));
    auto value = try_finite_numeric(*substituted, budget.context());
    return value ? static_cast<lmmc_real_t>(*value)
                 : std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

lmmc_real_t mixed_evaluate_with_retry(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x,
    lmmc_real_t half_width, RewriteBudget& budget, int max_retries) {
    auto value = mixed_evaluate_at(expr, var, x, budget);
    if (!std::isnan(value)) { return value; }
    auto offset = half_width;
    for (int i = 0; i < max_retries; ++i) {
        offset *= 0.5;
        value = mixed_evaluate_at(expr, var, x + offset, budget);
        if (!std::isnan(value)) { return value; }
        value = mixed_evaluate_at(expr, var, x - offset, budget);
        if (!std::isnan(value)) { return value; }
    }
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

bool consume_mixed_step(
    ComputationContext& context, std::optional<CasError>* failure,
    bool* complete) {
    auto step = context.consume_steps(1, "solve_mixed_transcendental_checked");
    if (step) { return true; }
    if (failure) { *failure = step.error(); }
    if (complete) { *complete = false; }
    return false;
}

}
