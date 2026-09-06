#pragma once

#include "inequality_solver.hpp"
#include "poly_utils.hpp"
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS::detail::inequality_support {

Result<IntervalUnion> solve_exact_affine_inequality(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context);

Result<IntervalUnion> solve_exact_polynomial_inequality(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context);

std::optional<double> try_checked_numeric_constant(const SymbolicExpr& expression);
int exact_numeric_sign(const std::shared_ptr<SymbolicExpr>& expression);
int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& polynomial);

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>>
find_roots_with_multiplicity(const std::shared_ptr<SymbolicExpr>& expression,
                             const std::string& variable);

bool root_less_than(const std::shared_ptr<SymbolicExpr>& lhs,
                    const std::shared_ptr<SymbolicExpr>& rhs);
bool roots_equal(const std::shared_ptr<SymbolicExpr>& lhs,
                 const std::shared_ptr<SymbolicExpr>& rhs);

bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expression,
                          const std::vector<std::string>& parameters);

std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const std::string& variable);

} // namespace LMCAS::detail::inequality_support
