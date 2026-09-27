#pragma once
#include "solve_polynomial.hpp"
#include <optional>

namespace LMCAS::polynomial_solver_detail {
bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr);
std::optional<double> finite_numeric_value(const std::shared_ptr<SymbolicExpr>& expr);
std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& expr);
std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a,
                                 const std::shared_ptr<SymbolicExpr>& b);
std::shared_ptr<SymbolicExpr> num(int value);
std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c);
bool convert_to_rational_poly(const Polynomial<SymbolicPolyCoeff>& poly,
                              Polynomial<Rational>& output);
bool bounded_rational_root_search(const Polynomial<Rational>& polynomial);
struct NumericCubicDepression { double p, q, shift; };
std::optional<NumericCubicDepression> numeric_cubic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d);
std::vector<std::shared_ptr<SymbolicExpr>> numeric_cubic_roots(
    const NumericCubicDepression& depression);
struct NumericQuarticDepression { double p, q, r, shift; };
std::optional<NumericQuarticDepression> numeric_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e);
std::vector<std::shared_ptr<SymbolicExpr>> numeric_quartic_roots(
    const NumericQuarticDepression& depression, const std::string& variable);
struct SymbolicQuarticDepression {
    std::shared_ptr<SymbolicExpr> p, q, r, shift;
};
SymbolicQuarticDepression symbolic_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e);
std::vector<std::shared_ptr<SymbolicExpr>> symbolic_quartic_roots(
    const SymbolicQuarticDepression& depression, const std::string& variable);
}
