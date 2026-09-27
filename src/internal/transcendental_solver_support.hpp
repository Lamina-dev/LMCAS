#pragma once

#include "solve_transcendental.hpp"
#include "internal/symbolic_ast.hpp"
#include "assumption.hpp"

namespace LMCAS::detail {
ExprPtr equation_predicate(const ExprPtr& expression);
ExprPtr substitute_raw(const ExprPtr& expression, const std::string& variable,
                       const ExprPtr& value);
SolveResult finalize_solution_set(
    const ExprPtr& original, const std::string& variable, SolutionSet solutions,
    ComputationContext& context, Domain domain, const SolveOptions& options);
}

namespace LMCAS::detail::transcendental {

struct InversePattern {
    FunctionNode::FuncType func_type;
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
    std::shared_ptr<SymbolicExpr> coefficient;
};

struct LambertWPattern {
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
};

struct ExpBasePattern {
    std::shared_ptr<SymbolicExpr> base;
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
};

bool try_evaluate_numeric(const std::shared_ptr<SymbolicExpr>& expr, lmmc_real_t& out);
std::optional<InversePattern> decompose_trig_exp_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
std::optional<LambertWPattern> decompose_lambert_w_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
std::optional<ExpBasePattern> decompose_exp_base_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
SolveResult invert_function(
    const InversePattern& pattern, const ExprPtr& original, const std::string& var,
    ComputationContext& context, const SolveOptions& options);
SolutionSet conditional_preimage(const ExprPtr& original, const std::string& var);
std::vector<std::shared_ptr<SymbolicExpr>> invert_lambert_w(
    const std::shared_ptr<SymbolicExpr>& c);
std::vector<std::shared_ptr<SymbolicExpr>> invert_exp_base(
    const std::shared_ptr<SymbolicExpr>& base, const std::shared_ptr<SymbolicExpr>& c);

}
