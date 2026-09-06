#pragma once

#include "transcendental_factor.hpp"
#include "symbolic_ast.hpp"

namespace LMCAS {

bool tf_is_transcendental_type(FunctionNode::FuncType function_type);
int tf_degree_in(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);
std::shared_ptr<SymbolicExpr> tf_back_substitute(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::vector<TransSubstitution>& mappings);
bool tf_is_linear_irreducible(
    const TransSubstitutionResult& substitution,
    const std::string& variable);
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_multiplicative_structure(
    const std::shared_ptr<SymbolicExpr>& expression);
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_exponential_separation(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);
bool tf_contains_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);
std::shared_ptr<SymbolicExpr> tf_simplify_pythagorean(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);

} // namespace LMCAS
