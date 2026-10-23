#pragma once
#include "solve_polynomial.hpp"
#include "computation_context.hpp"
#include <optional>

namespace LMCAS::polynomial_solver_detail {
bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr);
std::optional<double> finite_numeric_value(const std::shared_ptr<SymbolicExpr>& expr,
                                         ComputationContext& context);
std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& expr);
std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a,
                                 const std::shared_ptr<SymbolicExpr>& b);
std::shared_ptr<SymbolicExpr> num(int value);
std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, ComputationContext& context);
bool convert_to_rational_poly(const Polynomial<SymbolicPolyCoeff>& poly,
                              Polynomial<Rational>& output, ComputationContext& context);
bool bounded_rational_root_search(const Polynomial<Rational>& polynomial);
struct NumericCubicDepression { double p, q, shift; };
std::optional<NumericCubicDepression> numeric_cubic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> numeric_cubic_roots(
    const NumericCubicDepression& depression, ComputationContext& context);
struct NumericQuarticDepression { double p, q, r, shift; };
std::optional<NumericQuarticDepression> numeric_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e, ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> numeric_quartic_roots(
    const NumericQuarticDepression& depression, const std::string& variable,
    ComputationContext& context);
struct SymbolicQuarticDepression {
    std::shared_ptr<SymbolicExpr> p, q, r, shift;
};
SymbolicQuarticDepression symbolic_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e, ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> symbolic_quartic_roots(
    const SymbolicQuarticDepression& depression, const std::string& variable,
    ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic_impl(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::string& variable, ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic_impl(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e, const std::string& variable,
    ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic_impl(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::string& variable,
    ComputationContext& context);
std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring_impl(
    const Polynomial<SymbolicPolyCoeff>& polynomial, const std::string& variable,
    ComputationContext& context);
}
