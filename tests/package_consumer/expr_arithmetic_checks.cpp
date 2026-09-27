#include "assumption_context.hpp"
#include "expr.hpp"
#include "value.hpp"
#include "symbolic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace LMCAS;

static int check_ordinary_symbol() {
    auto ordinary_i = LMCAS::sym("i");
    if (!ordinary_i || ordinary_i.value()->to_string() != "i") {
        std::cerr << "failed to allow lowercase i as an ordinary symbol\n";
        return 10;
    }
    return 0;
}

static int check_reserved_imaginary_symbol() {
    auto reserved_I = LMCAS::sym("I");
    if (reserved_I ||
        std::string(LMCAS::error_name(reserved_I.error())) !=
            "ImaginaryUnitReserved") {
        std::cerr << "failed to expose LMCAS imaginary unit alias diagnostic\n";
        return 10;
    }
    return 0;
}

static int check_reserved_phi() {
    auto reserved_phi = LMCAS::sym("phi");
    if (reserved_phi) {
        std::cerr << "failed to reserve LMCAS phi constant\n";
        return 10;
    }
    return 0;
}

static int check_reserved_constants() {
    auto reserved_e = LMCAS::sym("e");
    auto reserved_unicode_pi = LMCAS::sym("\xCF\x80");
    if (reserved_e || reserved_unicode_pi) {
        std::cerr << "failed to reserve LMCAS std.math constants\n";
        return 10;
    }
    return 0;
}

static int check_approximate_boundary() {
    auto lsr_approx = LMCAS::approx_real(0.5);
    auto lsr_approx_value =
        lsr_approx ? LMCAS::evalf(*lsr_approx.value())
                   : LMCAS::Result<LMCAS::ApproxReal>::failure(
                         LMCAS::CasErrc::InternalInvariant,
                         "approx_real construction failed", "consumer");
    auto lsr_nan_approx = LMCAS::approx_real(NAN);
    auto lsr_inf_approx = LMCAS::approx_real(INFINITY);
    if (!lsr_approx_value || lsr_approx_value.value().value != 0.5 ||
        lsr_nan_approx || lsr_inf_approx ||
        std::string(LMCAS::error_name(lsr_nan_approx.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_inf_approx.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose finite LMCAS approx_real boundary\n";
        return 10;
    }
    return 0;
}

static int check_pi_constants() {
    auto lsr_pi = LMCAS::pi();
    auto lsr_pi_value = lsr_pi ? LMCAS::evalf(*lsr_pi.value())
                               : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                     LMCAS::CasErrc::InternalInvariant,
                                     "pi construction failed", "consumer");
    auto lsr_unicode_pi_value = LMCAS::evalf(
        *SymbolicExpr::variable("\xCF\x80"));
    if (!lsr_pi_value ||
        !lsr_unicode_pi_value ||
        std::abs(lsr_pi_value.value().value - 3.14159265358979323846) > 1e-15 ||
        std::abs(lsr_unicode_pi_value.value().value - 3.14159265358979323846) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math constants\n";
        return 10;
    }
    return 0;
}

static int check_e_phi_constants() {
    auto lsr_e = LMCAS::e();
    auto lsr_phi = LMCAS::phi();
    auto lsr_e_value = lsr_e ? LMCAS::evalf(*lsr_e.value())
                             : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                   LMCAS::CasErrc::InternalInvariant,
                                   "e construction failed", "consumer");
    auto lsr_phi_value = lsr_phi ? LMCAS::evalf(*lsr_phi.value())
                                 : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                       LMCAS::CasErrc::InternalInvariant,
                                       "phi construction failed", "consumer");
    if (!lsr_e_value ||
        !lsr_phi_value ||
        std::abs(lsr_e_value.value().value - std::exp(1.0)) > 1e-15 ||
        std::abs(lsr_phi_value.value().value -
                 ((1.0 + std::sqrt(5.0)) / 2.0)) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math constants\n";
        return 10;
    }
    return 0;
}

static int check_arithmetic_value() {
    auto expr_two = LMCAS::integer(2);
    auto expr_three = LMCAS::integer(3);
    auto expr_sum =
        LMCAS::add(expr_two.value(), expr_three.value());
    auto expr_product =
        LMCAS::mul(expr_sum.value(), SymbolicExpr::number(4));
    auto expr_quotient =
        LMCAS::div(expr_product.value(), SymbolicExpr::number(2));
    auto expr_difference =
        LMCAS::sub(expr_quotient.value(), SymbolicExpr::number(5));
    auto expr_negated = LMCAS::neg(expr_difference.value());
    auto expr_value =
        expr_negated ? LMCAS::evalf(*expr_negated.value())
                         : LMCAS::Result<LMCAS::ApproxReal>::failure(
                               LMCAS::CasErrc::InternalInvariant,
                               "Expr arithmetic construction failed",
                               "consumer");
    if (!expr_value ||
        std::abs(expr_value.value().value + 5.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS Expr arithmetic wrappers\n";
        return 10;
    }
    return 0;
}

static int check_linear_solution() {
    auto x = SymbolicExpr::variable("x");
    auto expr_polynomial =
        LMCAS::add(x, SymbolicExpr::number(1));
    auto expr_solved =
        expr_polynomial
            ? LMCAS::solve_expr_set(expr_polynomial.value(), "x")
            : LMCAS::Result<LMCAS::ExprSet>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "Expr polynomial construction failed", "consumer");
    if (!expr_solved ||
        !expr_solved.value().contains(*SymbolicExpr::number(-1))) {
        std::cerr << "failed to expose LMCAS Expr arithmetic wrappers\n";
        return 10;
    }
    return 0;
}

static int check_simplification() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_transform_zero = LMCAS::integer(0);
    auto lsr_x_plus_zero =
        LMCAS::add(x, lsr_transform_zero.value());
    auto lsr_simplified =
        lsr_x_plus_zero ? LMCAS::simplify(lsr_x_plus_zero.value())
                        : LMCAS::Result<LMCAS::ExprPtr>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "simplify input failed", "consumer");
    auto lsr_simplified_value =
        lsr_simplified ? LMCAS::evalf(*lsr_simplified.value(),
                                            LMCAS::NumericBindings{{"x", 7.0}})
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "simplify failed", "consumer");
    if (!lsr_simplified_value ||
        std::abs(lsr_simplified_value.value().value - 7.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS Expr transform wrappers\n";
        return 10;
    }
    return 0;
}

static int check_expansion() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_transform_one = LMCAS::integer(1);
    auto lsr_transform_two = LMCAS::integer(2);
    auto lsr_left = LMCAS::add(x, lsr_transform_one.value());
    auto lsr_right = LMCAS::add(x, lsr_transform_two.value());
    auto lsr_product =
        lsr_left && lsr_right
            ? LMCAS::mul(lsr_left.value(), lsr_right.value())
            : LMCAS::Result<LMCAS::ExprPtr>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "expand input failed", "consumer");
    auto lsr_expanded =
        lsr_product ? LMCAS::expand(lsr_product.value())
                    : LMCAS::Result<LMCAS::ExprPtr>::failure(
                          LMCAS::CasErrc::InternalInvariant,
                          "expand input failed", "consumer");
    auto lsr_expanded_value =
        lsr_expanded ? LMCAS::evalf(*lsr_expanded.value(),
                                          LMCAS::NumericBindings{{"x", 3.0}})
                     : LMCAS::Result<LMCAS::ApproxReal>::failure(
                           LMCAS::CasErrc::InternalInvariant,
                           "expand failed", "consumer");
    if (!lsr_expanded_value ||
        std::abs(lsr_expanded_value.value().value - 20.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS Expr transform wrappers\n";
        return 10;
    }
    return 0;
}

static int check_derivative() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_transform_three = LMCAS::integer(3);
    auto lsr_x_cubed = LMCAS::pow(x, lsr_transform_three.value());
    auto lsr_derivative =
        lsr_x_cubed ? LMCAS::differentiate(lsr_x_cubed.value(), "x")
                    : LMCAS::Result<LMCAS::ExprPtr>::failure(
                          LMCAS::CasErrc::InternalInvariant,
                          "differentiate input failed", "consumer");
    auto lsr_derivative_value =
        lsr_derivative ? LMCAS::evalf(*lsr_derivative.value(),
                                            LMCAS::NumericBindings{{"x", 2.0}})
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "differentiate failed", "consumer");
    if (!lsr_derivative_value ||
        std::abs(lsr_derivative_value.value().value - 12.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS Expr transform wrappers\n";
        return 10;
    }
    return 0;
}

static LMCAS::Result<LMCAS::ApproxReal> evaluate_math_expression(
    const LMCAS::ExprResult& expression, const char* failure_message) {
    if (!expression) {
        return LMCAS::Result<LMCAS::ApproxReal>::failure(
            LMCAS::CasErrc::InternalInvariant, failure_message, "consumer");
    }
    return LMCAS::evalf(*expression.value());
}

static int check_elementary_math() {
    auto lsr_math_sin = LMCAS::sin(SymbolicExpr::number(0));
    auto lsr_math_sqrt = LMCAS::sqrt(SymbolicExpr::number(9));
    auto lsr_math_pow =
        LMCAS::pow(SymbolicExpr::number(2), SymbolicExpr::number(4));
    auto lsr_math_sin_value =
        evaluate_math_expression(lsr_math_sin, "sin construction failed");
    auto lsr_math_sqrt_value =
        evaluate_math_expression(lsr_math_sqrt, "sqrt construction failed");
    auto lsr_math_pow_value =
        evaluate_math_expression(lsr_math_pow, "pow construction failed");
    if (!lsr_math_sin_value ||
        !lsr_math_sqrt_value ||
        !lsr_math_pow_value) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    if (std::abs(lsr_math_sin_value.value().value) > 1e-15 ||
        std::abs(lsr_math_sqrt_value.value().value - 3.0) > 1e-15 ||
        std::abs(lsr_math_pow_value.value().value - 16.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    return 0;
}

static int check_inverse_log_math() {
    auto lsr_math_asin = LMCAS::asin(SymbolicExpr::number(0.5));
    auto lsr_math_log10 = LMCAS::log10(SymbolicExpr::number(100));
    auto lsr_math_asin_value =
        evaluate_math_expression(lsr_math_asin, "asin construction failed");
    auto lsr_math_log10_value =
        evaluate_math_expression(lsr_math_log10, "log10 construction failed");
    if (!lsr_math_asin_value ||
        !lsr_math_log10_value ||
        std::abs(lsr_math_asin_value.value().value - std::asin(0.5)) > 1e-12 ||
        std::abs(lsr_math_log10_value.value().value - 2.0) > 1e-12) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    return 0;
}

static int check_rounding() {
    auto lsr_math_floor = LMCAS::floor(SymbolicExpr::number(2.75));
    auto lsr_math_ceil = LMCAS::ceil(SymbolicExpr::number(2.25));
    auto lsr_math_round = LMCAS::round(SymbolicExpr::number(-2.5));
    auto lsr_math_floor_value =
        evaluate_math_expression(lsr_math_floor, "floor construction failed");
    auto lsr_math_ceil_value =
        evaluate_math_expression(lsr_math_ceil, "ceil construction failed");
    auto lsr_math_round_value =
        evaluate_math_expression(lsr_math_round, "round construction failed");
    if (!lsr_math_floor_value ||
        !lsr_math_ceil_value ||
        !lsr_math_round_value) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    if (std::abs(lsr_math_floor_value.value().value - 2.0) > 1e-15 ||
        std::abs(lsr_math_ceil_value.value().value - 3.0) > 1e-15 ||
        std::abs(lsr_math_round_value.value().value + 3.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    return 0;
}

static int check_clamp() {
    auto lsr_math_clamp = LMCAS::clamp(SymbolicExpr::number(7),
                                             SymbolicExpr::number(0),
                                             SymbolicExpr::number(5));
    auto lsr_math_clamp_value =
        lsr_math_clamp ? LMCAS::evalf(*lsr_math_clamp.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "clamp construction failed", "consumer");
    if (!lsr_math_clamp_value ||
        std::abs(lsr_math_clamp_value.value().value - 5.0) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    return 0;
}

int run_expr_arithmetic_consumer_checks() {
    const auto checks = {
        check_ordinary_symbol,
        check_reserved_imaginary_symbol,
        check_reserved_phi,
        check_reserved_constants,
        check_approximate_boundary,
        check_pi_constants,
        check_e_phi_constants,
        check_arithmetic_value,
        check_linear_solution,
        check_simplification,
        check_expansion,
        check_derivative,
        check_elementary_math,
        check_inverse_log_math,
        check_rounding,
        check_clamp,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
