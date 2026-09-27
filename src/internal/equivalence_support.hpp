#pragma once

#include "internal/equivalence_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/normalization_utils.hpp"

#include <optional>

namespace LMCAS::equivalence_detail {

Result<void> validate_eqv_options(const EqvOptions& options);
bool exact_integer_node(const detail::SymbolicNodePtr& node, int expected);
ExprPtr normalize_equivalence(const ExprPtr& expression, detail::RewriteBudget& budget);
ExprPtr expand_equivalence(const ExprPtr& expression, detail::RewriteBudget& budget);
ExprPtr rewrite_trig_basic_identity(const detail::SymbolicNodePtr& node,
                                     detail::RewriteBudget& budget);
ExprPtr rewrite_exp_log_basic_identity(const detail::SymbolicNodePtr& node,
                                       const AssumptionContext* assumptions,
                                       detail::RewriteBudget& budget);
ExprPtr canonicalize_complex_product(const SymbolicExpr& expression,
                                     detail::RewriteBudget& budget);
Result<std::optional<bool>> prove_rational_polynomial_equivalence(
    const ExprPtr& difference, ComputationContext& context,
    const EqvOptions& options, detail::RewriteBudget& budget);
ExprPtr rewritten_complex(const ExprPtr& real, const ExprPtr& imag,
                           const detail::SymbolicNodePtr& original,
                           detail::RewriteBudget& budget);

}
