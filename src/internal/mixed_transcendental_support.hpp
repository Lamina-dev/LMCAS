#pragma once

#include "solve_mixed_transcendental.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::detail {

lmmc_real_t mixed_recursive_eval(const std::shared_ptr<const SymbolicNode>& node);
lmmc_real_t mixed_evaluate_at(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x);
lmmc_real_t mixed_evaluate_with_retry(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x,
    lmmc_real_t half_width, int max_retries = 3);
bool consume_mixed_step(
    ComputationContext* context, std::optional<CasError>* failure,
    bool* complete);
std::vector<IsolatedInterval> isolate_roots_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const SearchInterval& interval,
    const SolveOptions& opts, ComputationContext* context,
    std::optional<CasError>* failure, bool* complete);
std::optional<NumericRoot> refine_root_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const IsolatedInterval& interval,
    const SolveOptions& opts, ComputationContext* context,
    std::optional<CasError>* failure, bool* complete);
Result<std::vector<NumericRoot>> mixed_numerical_path(
    const std::shared_ptr<SymbolicExpr>& factor,
    const std::string& var, const SearchInterval& interval,
    const SolveOptions& opts, ComputationContext& context, bool& complete);

}
