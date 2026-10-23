#include "computation_context.hpp"
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

std::optional<double> try_checked_numeric_constant(
    const SymbolicExpr& expression, ComputationContext& context);
int exact_numeric_sign(const std::shared_ptr<SymbolicExpr>& expression,
                       ComputationContext& context);
int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& polynomial,
                           ComputationContext& context);

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>>
find_roots_with_multiplicity(const std::shared_ptr<SymbolicExpr>& expression,
                             const std::string& variable, ComputationContext& context);

bool root_less_than(const std::shared_ptr<SymbolicExpr>& lhs,
                    const std::shared_ptr<SymbolicExpr>& rhs, ComputationContext& context);
bool roots_equal(const std::shared_ptr<SymbolicExpr>& lhs,
                 const std::shared_ptr<SymbolicExpr>& rhs, ComputationContext& context);

bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expression,
                          const std::vector<std::string>& parameters);

std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const std::string& variable, ComputationContext& context);

void include_numeric_root(std::vector<Interval>& intervals,
                          const std::shared_ptr<SymbolicExpr>& root,
                          ComputationContext& context);
IntervalUnion build_parametric_intervals(
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities, int leading_sign, InequalityType type, ComputationContext& context);

}
