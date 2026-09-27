#include "internal/symbolic_ast.hpp"
#include "internal/numeric_probe.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS::polynomial_solver_detail {

bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    return !expr || free_variables(LMCAS::detail::node(expr)).empty();
}

std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
}

std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    return SymbolicExpr::add(a, negate(b));
}

std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

std::optional<double> finite_numeric_value(
    const std::shared_ptr<SymbolicExpr>& expr) {
    auto evaluated = detail::try_finite_numeric(expr);
    if (!evaluated) {
        return std::nullopt;
    }
    auto simplified = expr->simplify();
    if (simplified && LMCAS::detail::node(simplified)) {
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(
                LMCAS::detail::node(simplified))) {
            if (std::holds_alternative<Rational>(num->value())) {
                const auto& value = std::get<Rational>(num->value());
                if (!value.get_numerator().is_zero() && *evaluated == 0.0) {
                    return std::nullopt;
                }
            } else if (std::holds_alternative<BigInt>(num->value())) {
                const auto& value = std::get<BigInt>(num->value());
                if (!value.is_zero() && *evaluated == 0.0) {
                    return std::nullopt;
                }
            }
        }
    }
    return evaluated;
}

bool convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly)
{
    out_poly = Polynomial<Rational>(sym_poly.variable_name);
    out_poly.coeffs.resize(sym_poly.coeffs.size(), Rational(0));

    for (size_t i = 0; i < sym_poly.coeffs.size(); ++i) {
        auto coeff_expr = sym_poly.coeffs[i].val;
        if (!coeff_expr) {
            out_poly.coeffs[i] = Rational(0);
            continue;
        }

        auto simplified = coeff_expr->simplify();
        auto num_node = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified));
        if (!num_node) {
            return false;
        }
        if (std::holds_alternative<Rational>(num_node->value())) {
            out_poly.coeffs[i] = std::get<Rational>(num_node->value());
        } else if (std::holds_alternative<BigInt>(num_node->value())) {
            out_poly.coeffs[i] =
                Rational(std::get<BigInt>(num_node->value()));
        } else {
            return false;
        }
    }

    out_poly.trim();
    return true;
}

}
