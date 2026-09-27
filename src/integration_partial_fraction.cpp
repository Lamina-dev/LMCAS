#include "internal/integration_support.hpp"

namespace LMCAS {

namespace {
std::shared_ptr<SymbolicExpr> quadratic_reciprocal_denominator(const SymbolicExpr& expr) {
    std::shared_ptr<SymbolicExpr> den = nullptr;

    if (auto p = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        double exp_val = 0;
        bool is_inv = false;
        if (auto num_node = std::dynamic_pointer_cast<const NumberNode>(p->exponent())) {
            if (std::holds_alternative<lmmc_real_t>(num_node->value()))
                exp_val = std::get<lmmc_real_t>(num_node->value());
            else if (std::holds_alternative<Rational>(num_node->value()))
                exp_val = std::get<Rational>(num_node->value()).to_double();
            else if (std::holds_alternative<BigInt>(num_node->value()))
                exp_val = std::get<BigInt>(num_node->value()).to_double();
            int eq;
            lmmc_double_nearly_equal_tol(exp_val, -1.0, 1e-9, 1e-9, &eq);
            if (eq) is_inv = true;
        }
        if (is_inv) {
            den = detail::make_expression_ptr(LMCAS::detail::expression_from_node(p->base()));
        }
    }
    return den;
}

std::shared_ptr<SymbolicExpr> integrate_numeric_quadratic(
    const Polynomial<SymbolicPolyCoeff>& Q, const std::string& var,
    ComputationContext& computation) {
            SymbolicExpr c_expr = *(Q.coeffs[0].val);
            SymbolicExpr b_expr = *(Q.coeffs[1].val);
            SymbolicExpr a_expr = *(Q.coeffs[2].val);

            if (!a_expr.is_number() || !b_expr.is_number() || !c_expr.is_number()) {
                return nullptr;
            }

            auto a_checked = try_checked_numeric_constant(a_expr, computation);
            auto b_checked = try_checked_numeric_constant(b_expr, computation);
            auto c_checked = try_checked_numeric_constant(c_expr, computation);
            if (!a_checked || !b_checked || !c_checked) {
                return nullptr;
            }

            double a = *a_checked;
            double b = *b_checked;
            double c = *c_checked;

            int eq_a;
            lmmc_double_nearly_equal_tol(a, 0.0, 1e-9, 1e-9, &eq_a);
            if (eq_a) { return nullptr; }

            double delta = b * b - 4 * a * c;
            int eq_delta;
            lmmc_double_nearly_equal_tol(delta, 0.0, 1e-9, 1e-9, &eq_delta);

            if (!eq_delta && delta > 0) {
                double sqrt_delta = std::sqrt(delta);
                auto scalar = SymbolicExpr::number(1.0 / sqrt_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto two_a_x = SymbolicExpr::multiply(detail::make_expression_ptr(*two_a), SymbolicExpr::variable(var));
                auto two_a_x_plus_b = SymbolicExpr::add(detail::make_expression_ptr(*two_a_x), detail::make_expression_ptr(*b_num));
                auto term1_arg = sym_sub(*two_a_x_plus_b, *SymbolicExpr::number(sqrt_delta));
                auto term2_arg = SymbolicExpr::add(detail::make_expression_ptr(*two_a_x_plus_b), SymbolicExpr::number(sqrt_delta));
                auto term1 = SymbolicExpr::ln(detail::make_expression_ptr(*term1_arg));
                auto term2 = SymbolicExpr::ln(detail::make_expression_ptr(*term2_arg));
                return SymbolicExpr::multiply(scalar, sym_sub(*term1, *term2));
            } else if (!eq_delta && delta < 0) {
                double sqrt_neg_delta = std::sqrt(-delta);
                auto scalar = SymbolicExpr::number(2.0 / sqrt_neg_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto num = SymbolicExpr::add(
                    SymbolicExpr::multiply(detail::make_expression_ptr(*two_a), SymbolicExpr::variable(var)),
                    detail::make_expression_ptr(*b_num));
                auto inner = SymbolicExpr::divide(detail::make_expression_ptr(*num), SymbolicExpr::number(sqrt_neg_delta));
                return SymbolicExpr::multiply(scalar, make_arctan(inner));
            }
    return nullptr;
}
}

Result<std::shared_ptr<SymbolicExpr>> PartialFractionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& computation, int) {
    auto den = quadratic_reciprocal_denominator(expr);
    if (!den) { return nullptr; }
    auto Q = symbolic_to_poly<SymbolicPolyCoeff>(den, var);
    if (!Q) {
        if (Q.error().code == CasErrc::UnsupportedExpression) return nullptr;
        return Result<std::shared_ptr<SymbolicExpr>>::failure(Q.error());
    }
    if (Q.value().degree() == 2) {
        return integrate_numeric_quadratic(Q.value(), var, computation);
    }

    return nullptr;
}

}
