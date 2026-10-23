#pragma once
#include "computation_context.hpp"

#include "transcendental_factor.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

bool tf_is_transcendental_type(FunctionNode::FuncType function_type);
int tf_degree_in(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);
std::shared_ptr<SymbolicExpr> tf_back_substitute(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::vector<TransSubstitution>& mappings, ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors_impl(
    const std::vector<std::shared_ptr<SymbolicExpr>>& factors, ComputationContext& context);
Result<TfPolyBuildResult> tf_build_polynomial_impl(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::vector<std::string>& indeterminates,
    const std::string& original_variable, ComputationContext& context);
bool tf_is_linear_irreducible(
    const TransSubstitutionResult& substitution,
    const std::string& variable);
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_multiplicative_structure(
    const std::shared_ptr<SymbolicExpr>& expression);
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_exponential_separation(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable, ComputationContext& context);
bool tf_contains_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);
std::shared_ptr<SymbolicExpr> tf_simplify_pythagorean(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);

} // namespace LMCAS
