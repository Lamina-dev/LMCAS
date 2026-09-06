#include "assumption_context.hpp"
#include "expr.hpp"
#include "value.hpp"
#include "symbolic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace LMCAS;

int run_expr_consumer_checks() {
    auto x = SymbolicExpr::variable("x");
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto upper_i = LMCAS::I();
    if (!upper_i ||
        !LMCAS::structurally_equal(*i.value(), *upper_i.value())) {
        std::cerr << "failed to expose std.math.I\n";
        return 10;
    }
    auto null_real_complex =
        LMCAS::complex(nullptr, SymbolicExpr::number(1));
    auto null_imag_complex =
        LMCAS::complex(SymbolicExpr::number(0), nullptr);
    if (null_real_complex || null_imag_complex ||
        std::string(LMCAS::error_name(null_real_complex.error())) !=
            "ComplexTypeMismatch" ||
        std::string(LMCAS::error_name(null_imag_complex.error())) !=
            "ComplexTypeMismatch") {
        std::cerr << "failed to expose LMCAS complex type diagnostics\n";
        return 10;
    }
    auto ordinary_i = LMCAS::sym("i");
    if (!ordinary_i || ordinary_i.value()->to_string() != "i") {
        std::cerr << "failed to allow lowercase i as an ordinary symbol\n";
        return 10;
    }
    auto reserved_I = LMCAS::sym("I");
    if (reserved_I ||
        std::string(LMCAS::error_name(reserved_I.error())) !=
            "ImaginaryUnitReserved") {
        std::cerr << "failed to expose LMCAS imaginary unit alias diagnostic\n";
        return 10;
    }
    auto reserved_phi = LMCAS::sym("phi");
    if (reserved_phi) {
        std::cerr << "failed to reserve LMCAS phi constant\n";
        return 10;
    }
    auto reserved_e = LMCAS::sym("e");
    auto reserved_unicode_pi = LMCAS::sym("\xCF\x80");
    if (reserved_e || reserved_unicode_pi) {
        std::cerr << "failed to reserve LMCAS std.math constants\n";
        return 10;
    }
    auto lsr_pi = LMCAS::pi();
    auto lsr_e = LMCAS::e();
    auto lsr_phi = LMCAS::phi();
    auto lsr_pi_value = lsr_pi ? LMCAS::evalf(*lsr_pi.value())
                               : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                     LMCAS::CasErrc::InternalInvariant,
                                     "pi construction failed", "consumer");
    auto lsr_e_value = lsr_e ? LMCAS::evalf(*lsr_e.value())
                             : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                   LMCAS::CasErrc::InternalInvariant,
                                   "e construction failed", "consumer");
    auto lsr_phi_value = lsr_phi ? LMCAS::evalf(*lsr_phi.value())
                                 : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                       LMCAS::CasErrc::InternalInvariant,
                                       "phi construction failed", "consumer");
    auto lsr_unicode_pi_value = LMCAS::evalf(
        *SymbolicExpr::variable("\xCF\x80"));
    if (!lsr_pi_value || !lsr_e_value || !lsr_phi_value ||
        !lsr_unicode_pi_value ||
        std::abs(lsr_pi_value.value().value - 3.14159265358979323846) > 1e-15 ||
        std::abs(lsr_unicode_pi_value.value().value - 3.14159265358979323846) > 1e-15 ||
        std::abs(lsr_e_value.value().value - std::exp(1.0)) > 1e-15 ||
        std::abs(lsr_phi_value.value().value -
                 ((1.0 + std::sqrt(5.0)) / 2.0)) > 1e-15) {
        std::cerr << "failed to expose LMCAS std.math constants\n";
        return 10;
    }
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
    auto expr_polynomial =
        LMCAS::add(x, SymbolicExpr::number(1));
    auto expr_equation =
        LMCAS::eq(expr_polynomial.value(), SymbolicExpr::number(0));
    auto expr_solved =
        expr_polynomial
            ? LMCAS::solve_expr_set(expr_polynomial.value(), "x")
            : LMCAS::Result<LMCAS::ExprSet>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "Expr polynomial construction failed", "consumer");
    auto expr_null = LMCAS::add(nullptr, expr_two.value());
    if (!expr_value || std::abs(expr_value.value().value + 5.0) > 1e-15 ||
        !expr_solved ||
        !expr_equation || expr_equation.value()->to_string().empty() ||
        !expr_solved.value().contains(*SymbolicExpr::number(-1)) ||
        expr_null ||
        std::string(LMCAS::error_name(expr_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS Expr arithmetic wrappers\n";
        return 10;
    }
    auto lsr_transform_zero = LMCAS::integer(0);
    auto lsr_transform_one = LMCAS::integer(1);
    auto lsr_transform_two = LMCAS::integer(2);
    auto lsr_transform_three = LMCAS::integer(3);
    auto lsr_x_plus_zero =
        LMCAS::add(x, lsr_transform_zero.value());
    auto lsr_simplified =
        lsr_x_plus_zero ? LMCAS::simplify(lsr_x_plus_zero.value())
                        : LMCAS::Result<LMCAS::ExprPtr>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "simplify input failed", "consumer");
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
    auto lsr_x_cubed = LMCAS::pow(x, lsr_transform_three.value());
    auto lsr_derivative =
        lsr_x_cubed ? LMCAS::differentiate(lsr_x_cubed.value(), "x")
                    : LMCAS::Result<LMCAS::ExprPtr>::failure(
                          LMCAS::CasErrc::InternalInvariant,
                          "differentiate input failed", "consumer");
    auto lsr_simplified_value =
        lsr_simplified ? LMCAS::evalf(*lsr_simplified.value(),
                                            LMCAS::NumericBindings{{"x", 7.0}})
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "simplify failed", "consumer");
    auto lsr_expanded_value =
        lsr_expanded ? LMCAS::evalf(*lsr_expanded.value(),
                                          LMCAS::NumericBindings{{"x", 3.0}})
                     : LMCAS::Result<LMCAS::ApproxReal>::failure(
                           LMCAS::CasErrc::InternalInvariant,
                           "expand failed", "consumer");
    auto lsr_derivative_value =
        lsr_derivative ? LMCAS::evalf(*lsr_derivative.value(),
                                            LMCAS::NumericBindings{{"x", 2.0}})
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "differentiate failed", "consumer");
    auto lsr_transform_null = LMCAS::expand(nullptr);
    if (!lsr_simplified_value || !lsr_expanded_value ||
        !lsr_derivative_value ||
        std::abs(lsr_simplified_value.value().value - 7.0) > 1e-15 ||
        std::abs(lsr_expanded_value.value().value - 20.0) > 1e-15 ||
        std::abs(lsr_derivative_value.value().value - 12.0) > 1e-15 ||
        lsr_transform_null ||
        std::string(LMCAS::error_name(lsr_transform_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS Expr transform wrappers\n";
        return 10;
    }
    auto lsr_math_sin = LMCAS::sin(SymbolicExpr::number(0));
    auto lsr_math_sqrt = LMCAS::sqrt(SymbolicExpr::number(9));
    auto lsr_math_pow =
        LMCAS::pow(SymbolicExpr::number(2), SymbolicExpr::number(4));
    auto lsr_math_asin = LMCAS::asin(SymbolicExpr::number(0.5));
    auto lsr_math_log10 = LMCAS::log10(SymbolicExpr::number(100));
    auto lsr_math_floor = LMCAS::floor(SymbolicExpr::number(2.75));
    auto lsr_math_ceil = LMCAS::ceil(SymbolicExpr::number(2.25));
    auto lsr_math_round = LMCAS::round(SymbolicExpr::number(-2.5));
    auto lsr_math_clamp = LMCAS::clamp(SymbolicExpr::number(7),
                                             SymbolicExpr::number(0),
                                             SymbolicExpr::number(5));
    auto lsr_math_sin_value =
        lsr_math_sin ? LMCAS::evalf(*lsr_math_sin.value())
                     : LMCAS::Result<LMCAS::ApproxReal>::failure(
                           LMCAS::CasErrc::InternalInvariant,
                           "sin construction failed", "consumer");
    auto lsr_math_sqrt_value =
        lsr_math_sqrt ? LMCAS::evalf(*lsr_math_sqrt.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "sqrt construction failed", "consumer");
    auto lsr_math_pow_value =
        lsr_math_pow ? LMCAS::evalf(*lsr_math_pow.value())
                     : LMCAS::Result<LMCAS::ApproxReal>::failure(
                           LMCAS::CasErrc::InternalInvariant,
                           "pow construction failed", "consumer");
    auto lsr_math_asin_value =
        lsr_math_asin ? LMCAS::evalf(*lsr_math_asin.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "asin construction failed", "consumer");
    auto lsr_math_log10_value =
        lsr_math_log10 ? LMCAS::evalf(*lsr_math_log10.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "log10 construction failed", "consumer");
    auto lsr_math_floor_value =
        lsr_math_floor ? LMCAS::evalf(*lsr_math_floor.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "floor construction failed", "consumer");
    auto lsr_math_ceil_value =
        lsr_math_ceil ? LMCAS::evalf(*lsr_math_ceil.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "ceil construction failed", "consumer");
    auto lsr_math_round_value =
        lsr_math_round ? LMCAS::evalf(*lsr_math_round.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "round construction failed", "consumer");
    auto lsr_math_clamp_value =
        lsr_math_clamp ? LMCAS::evalf(*lsr_math_clamp.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "clamp construction failed", "consumer");
    auto lsr_math_null = LMCAS::sin(nullptr);
    auto lsr_math_null_clamp = LMCAS::clamp(
        SymbolicExpr::number(1), nullptr, SymbolicExpr::number(2));
    if (!lsr_math_sin_value || !lsr_math_sqrt_value ||
        !lsr_math_pow_value || !lsr_math_asin_value ||
        !lsr_math_log10_value || !lsr_math_floor_value ||
        !lsr_math_ceil_value || !lsr_math_round_value ||
        !lsr_math_clamp_value ||
        std::abs(lsr_math_sin_value.value().value) > 1e-15 ||
        std::abs(lsr_math_sqrt_value.value().value - 3.0) > 1e-15 ||
        std::abs(lsr_math_pow_value.value().value - 16.0) > 1e-15 ||
        std::abs(lsr_math_asin_value.value().value - std::asin(0.5)) > 1e-12 ||
        std::abs(lsr_math_log10_value.value().value - 2.0) > 1e-12 ||
        std::abs(lsr_math_floor_value.value().value - 2.0) > 1e-15 ||
        std::abs(lsr_math_ceil_value.value().value - 3.0) > 1e-15 ||
        std::abs(lsr_math_round_value.value().value + 3.0) > 1e-15 ||
        std::abs(lsr_math_clamp_value.value().value - 5.0) > 1e-15 ||
        lsr_math_null ||
        lsr_math_null_clamp ||
        std::string(LMCAS::error_name(lsr_math_null.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_math_null_clamp.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
    LMCAS::ComputationContext lsr_context;
    auto i_rule = LMCAS::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), lsr_context);
    if (!i_rule || !i_rule.value()) {
        std::cerr << "failed to prove LMCAS i*i == -1\n";
        return 11;
    }
    auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
    LMCAS::ComputationContext lsr_power_context;
    auto i_power_rule = LMCAS::equivalent_core(
        *i_power_two, *SymbolicExpr::number(-1), lsr_power_context);
    if (!i_power_rule || !i_power_rule.value()) {
        std::cerr << "failed to prove LMCAS i^2 == -1\n";
        return 11;
    }
    auto legacy_i = SymbolicExpr::variable("i");
    auto legacy_i_squared = SymbolicExpr::multiply(legacy_i, legacy_i);
    LMCAS::ComputationContext legacy_i_context;
    auto legacy_i_rule = LMCAS::equivalent_core(
        *legacy_i_squared, *SymbolicExpr::number(-1), legacy_i_context);
    if (!legacy_i_rule || legacy_i_rule.value()) {
        std::cerr << "ordinary Expr i was treated as the imaginary unit\n";
        return 11;
    }
    LMCAS::EqvOptions valid_budget_options;
    if (!LMCAS::set_eqv_budget(valid_budget_options, 256, 64, 4)) {
        std::cerr << "failed to configure LMCAS equivalence budget\n";
        return 11;
    }
    LMCAS::EqvOptions invalid_budget_options;
    auto invalid_budget =
        LMCAS::set_eqv_budget(invalid_budget_options, 0, 64, 4);
    LMCAS::EqvOptions invalid_depth_budget_options;
    auto invalid_depth_budget =
        LMCAS::set_eqv_budget(invalid_depth_budget_options, 256, 0, 4);
    LMCAS::EqvOptions invalid_growth_budget_options;
    auto invalid_growth_budget =
        LMCAS::set_eqv_budget(invalid_growth_budget_options, 256, 64, 0);
    if (invalid_budget || invalid_depth_budget || invalid_growth_budget ||
        std::string(LMCAS::error_name(invalid_budget.error())) !=
            "EqvBudgetExceeded" ||
        std::string(LMCAS::error_name(invalid_depth_budget.error())) !=
            "EqvBudgetExceeded" ||
        std::string(LMCAS::error_name(invalid_growth_budget.error())) !=
            "EqvBudgetExceeded") {
        std::cerr << "failed to expose LMCAS equivalence budget setter diagnostics\n";
        return 11;
    }
    LMCAS::EqvOptions exhausted_eqv_options;
    (void)LMCAS::set_eqv_budget(exhausted_eqv_options, 1, 64, 4);
    exhausted_eqv_options.budget.max_rewrite_steps = 0;
    LMCAS::ComputationContext exhausted_eqv_context;
    auto exhausted_eqv = LMCAS::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), exhausted_eqv_context,
        exhausted_eqv_options);
    if (exhausted_eqv ||
        std::string(LMCAS::error_name(exhausted_eqv.error())) !=
            "EqvBudgetExceeded") {
        std::cerr << "failed to expose LMCAS equivalence budget diagnostics\n";
        return 11;
    }
    LMCAS::ComputationContext lsr_exhausted_eqv_context;
    auto lsr_exhausted_eqv = LMCAS::equivalent(
        *i_squared, *SymbolicExpr::number(-1), lsr_exhausted_eqv_context,
        exhausted_eqv_options);
    if (!lsr_exhausted_eqv || lsr_exhausted_eqv.value()) {
        std::cerr << "failed to return false for exhausted LMCAS equivalence\n";
        return 11;
    }
    auto x_plus_zero = SymbolicExpr::add(x, SymbolicExpr::number(0));
    LMCAS::ComputationContext identity_eqv_context;
    auto identity_eqv =
        LMCAS::equivalent_core(*x_plus_zero, *x, identity_eqv_context);
    if (!identity_eqv || !identity_eqv.value()) {
        std::cerr << "failed to prove LMCAS Core identity example\n";
        return 11;
    }
    auto x_times_one = SymbolicExpr::multiply(x, SymbolicExpr::number(1));
    auto x_times_zero = SymbolicExpr::multiply(x, SymbolicExpr::number(0));
    auto x_minus_x = SymbolicExpr::add(
        x, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    LMCAS::ComputationContext multiply_identity_context;
    auto multiply_identity = LMCAS::equivalent_core(
        *x_times_one, *x, multiply_identity_context);
    LMCAS::ComputationContext multiply_zero_context;
    auto multiply_zero = LMCAS::equivalent_core(
        *x_times_zero, *SymbolicExpr::number(0), multiply_zero_context);
    LMCAS::ComputationContext subtract_self_context;
    auto subtract_self = LMCAS::equivalent_core(
        *x_minus_x, *SymbolicExpr::number(0), subtract_self_context);
    if (!multiply_identity || !multiply_identity.value() ||
        !multiply_zero || !multiply_zero.value() ||
        !subtract_self || !subtract_self.value()) {
        std::cerr << "failed to prove LMCAS Core algebra identities\n";
        return 11;
    }
    auto x_plus_one = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto x_plus_one_squared = SymbolicExpr::power(x_plus_one, SymbolicExpr::number(2));
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
        SymbolicExpr::number(1));
    LMCAS::ComputationContext polynomial_eqv_context;
    auto polynomial_eqv = LMCAS::equivalent_core(
        *x_plus_one_squared, *expanded_square, polynomial_eqv_context);
    if (!polynomial_eqv || !polynomial_eqv.value()) {
        std::cerr << "failed to prove LMCAS polynomial equivalence example\n";
        return 11;
    }
    auto trig_identity = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        SymbolicExpr::power(SymbolicExpr::cos(x), SymbolicExpr::number(2)));
    LMCAS::EqvOptions trig_eqv_options;
    if (!LMCAS::set_eqv_profile(trig_eqv_options, "Trig-Basic")) {
        std::cerr << "failed to configure LMCAS Trig-Basic profile\n";
        return 11;
    }
    LMCAS::ComputationContext trig_eqv_context;
    auto trig_eqv = LMCAS::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), trig_eqv_context,
        trig_eqv_options);
    if (!trig_eqv || !trig_eqv.value()) {
        std::cerr << "failed to prove LMCAS Trig-Basic equivalence example\n";
        return 11;
    }
    LMCAS::ComputationContext lsr_trig_eqv_context;
    auto lsr_trig_eqv = LMCAS::equivalent(
        *trig_identity, *SymbolicExpr::number(1), lsr_trig_eqv_context,
        trig_eqv_options);
    if (!lsr_trig_eqv || !lsr_trig_eqv.value()) {
        std::cerr << "failed to prove LMCAS equivalent Trig-Basic example\n";
        return 11;
    }
    LMCAS::ComputationContext core_trig_eqv_context;
    auto core_trig_eqv = LMCAS::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), core_trig_eqv_context);
    if (!core_trig_eqv || core_trig_eqv.value()) {
        std::cerr << "enabled LMCAS Trig-Basic rules in Core profile\n";
        return 11;
    }
    LMCAS::EqvOptions exp_log_eqv_options;
    if (!LMCAS::set_eqv_profile(exp_log_eqv_options, "ExpLog-Basic")) {
        std::cerr << "failed to configure LMCAS ExpLog-Basic profile\n";
        return 11;
    }
    if (LMCAS::set_eqv_profile(exp_log_eqv_options,
                                     "Richardson-Complete")) {
        std::cerr << "accepted unsupported LMCAS equivalence profile\n";
        return 11;
    }
    LMCAS::ComputationContext exp_log_eqv_context;
    auto exp_log_eqv = LMCAS::equivalent_core(
        *SymbolicExpr::exp(SymbolicExpr::number(0)), *SymbolicExpr::number(1),
        exp_log_eqv_context, exp_log_eqv_options);
    if (!exp_log_eqv || !exp_log_eqv.value()) {
        std::cerr << "failed to prove LMCAS ExpLog-Basic equivalence example\n";
        return 11;
    }
    LMCAS::ComputationContext ln_one_eqv_context;
    auto ln_one_eqv = LMCAS::equivalent_core(
        *SymbolicExpr::ln(SymbolicExpr::number(1)), *SymbolicExpr::number(0),
        ln_one_eqv_context, exp_log_eqv_options);
    if (!ln_one_eqv || !ln_one_eqv.value()) {
        std::cerr << "failed to prove LMCAS ExpLog-Basic ln(1) example\n";
        return 11;
    }
    auto exp_ln_x = SymbolicExpr::exp(SymbolicExpr::ln(x));
    LMCAS::ComputationContext exp_ln_unproven_context;
    auto exp_ln_unproven = LMCAS::equivalent_core(
        *exp_ln_x, *x, exp_ln_unproven_context, exp_log_eqv_options);
    if (!exp_ln_unproven || exp_ln_unproven.value()) {
        std::cerr << "proved LMCAS exp(ln(x)) without domain evidence\n";
        return 11;
    }
    auto positive_assumptions = std::make_shared<LMCAS::AssumptionContext>();
    auto positive_assumption =
        positive_assumptions->assume_sign("x", LMCAS::Sign::Positive);
    if (!positive_assumption) {
        std::cerr << "failed to create positive LMCAS assumption\n";
        return 11;
    }
    LMCAS::ComputationContext exp_ln_positive_context;
    if (!exp_ln_positive_context.set_assumptions(positive_assumptions)) {
        std::cerr << "failed to attach LMCAS equivalence assumptions\n";
        return 11;
    }
    auto exp_ln_positive = LMCAS::equivalent_core(
        *exp_ln_x, *x, exp_ln_positive_context, exp_log_eqv_options);
    if (!exp_ln_positive || !exp_ln_positive.value()) {
        std::cerr << "failed to prove LMCAS exp(ln(x)) under positive assumption\n";
        return 11;
    }

    auto lsr_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_roots || lsr_roots.value().size() != 2) {
        std::cerr << "failed to lower LMCAS solve result to set<Expr>\n";
        return 12;
    }
    auto lsr_named_roots = LMCAS::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    auto lsr_named_solve = LMCAS::solve(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_named_roots || lsr_named_roots.value().size() != 2 ||
        !lsr_named_solve || lsr_named_solve.value().size() != 2) {
        std::cerr << "failed to call LMCAS roots/solve set<Expr> aliases\n";
        return 12;
    }
    auto lsr_repeated_roots = LMCAS::roots(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), "x");
    if (!lsr_repeated_roots || lsr_repeated_roots.value().size() != 1 ||
        !lsr_repeated_roots.value().contains(*SymbolicExpr::number(0))) {
        std::cerr << "failed to lower repeated LMCAS roots to set<Expr>\n";
        return 12;
    }
    auto lsr_cubic_roots = LMCAS::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(3)),
                          SymbolicExpr::number(-2)),
        "x");
    if (!lsr_cubic_roots || lsr_cubic_roots.value().size() != 3) {
        std::cerr << "failed to lower RootOf LMCAS roots to set<Expr>\n";
        return 12;
    }
    for (const auto& root : lsr_cubic_roots.value().elements()) {
        if (!root) {
            std::cerr << "failed to preserve LMCAS finite Expr roots\n";
            return 12;
        }
    }
    auto lsr_complex_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(1)),
        "x");
    auto negative_i = LMCAS::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(-1));
    if (!lsr_complex_roots || lsr_complex_roots.value().size() != 2 ||
        !lsr_complex_roots.value().contains(*i.value()) ||
        !negative_i ||
        !lsr_complex_roots.value().contains(*negative_i.value())) {
        std::cerr << "failed to lower LMCAS complex roots to set<Expr>\n";
        return 12;
    }
    auto lsr_complex_roots_subset_c = LMCAS::expr_set_subset_domain(
        lsr_complex_roots.value(), LMCAS::complexes());
    if (!lsr_complex_roots_subset_c ||
        !lsr_complex_roots_subset_c.value()) {
        std::cerr << "failed to prove installed LMCAS complex roots subset C\n";
        return 12;
    }
    auto lsr_shifted_complex_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                              SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
            SymbolicExpr::number(2)),
        "x");
    if (!lsr_shifted_complex_roots ||
        lsr_shifted_complex_roots.value().size() != 2) {
        std::cerr << "failed to lower shifted LMCAS complex roots to set<Expr>\n";
        return 12;
    }
    bool saw_negative_one_plus_i = false;
    bool saw_negative_one_minus_i = false;
    for (const auto& root : lsr_shifted_complex_roots.value().elements()) {
        auto lowered_root = LMCAS::eval_complex(*root);
        if (!lowered_root || !lowered_root.value().is_finite()) {
            std::cerr << "failed to evaluate shifted LMCAS complex root\n";
            return 12;
        }
        if (lowered_root.value().real.value == -1.0 &&
            lowered_root.value().imag.value == 1.0) {
            saw_negative_one_plus_i = true;
        }
        if (lowered_root.value().real.value == -1.0 &&
            lowered_root.value().imag.value == -1.0) {
            saw_negative_one_minus_i = true;
        }
    }
    if (!saw_negative_one_plus_i || !saw_negative_one_minus_i) {
        std::cerr << "failed to preserve shifted LMCAS complex root components\n";
        return 12;
    }
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    if (!lsr_set_left || !lsr_set_right ||
        lsr_set_left.value().size() != 2 ||
        !lsr_set_left.value().contains(*SymbolicExpr::number(1))) {
        std::cerr << "failed to construct LMCAS set<Expr>\n";
        return 12;
    }
    auto lsr_empty_set = LMCAS::expr_set({});
    if (!lsr_empty_set || !lsr_empty_set.value().empty()) {
        std::cerr << "failed to construct empty LMCAS set<Expr>\n";
        return 12;
    }
    auto lsr_set_union = lsr_set_left.value().set_union(lsr_set_right.value());
    auto lsr_set_intersection =
        lsr_set_left.value().intersection(lsr_set_right.value());
    auto lsr_set_difference =
        lsr_set_left.value().difference(lsr_set_right.value());
    auto lsr_set_symmetric =
        lsr_set_left.value().symmetric_difference(lsr_set_right.value());
    auto lsr_facade_contains = LMCAS::expr_set_contains(
        lsr_set_left.value(), SymbolicExpr::number(1));
    auto lsr_facade_not_contains = LMCAS::expr_set_not_contains(
        lsr_set_left.value(), SymbolicExpr::number(3));
    auto lsr_facade_union = LMCAS::expr_set_union(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_intersection = LMCAS::expr_set_intersection(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_difference = LMCAS::expr_set_difference(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_symmetric = LMCAS::expr_set_symmetric_difference(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_subset = LMCAS::expr_set_subset(
        lsr_set_intersection, lsr_set_union);
    auto lsr_facade_null_membership = LMCAS::expr_set_contains(
        lsr_set_left.value(), nullptr);
    if (lsr_set_union.size() != 3 ||
        lsr_set_intersection.size() != 1 ||
        !lsr_set_intersection.contains(*SymbolicExpr::number(2)) ||
        lsr_set_difference.size() != 1 ||
        !lsr_set_difference.contains(*SymbolicExpr::number(1)) ||
        lsr_set_symmetric.size() != 2 ||
        !lsr_set_symmetric.contains(*SymbolicExpr::number(1)) ||
        !lsr_set_symmetric.contains(*SymbolicExpr::number(3)) ||
        !lsr_set_intersection.subset_of(lsr_set_union) ||
        !lsr_empty_set.value().subset_of(lsr_set_left.value()) ||
        lsr_set_left.value().set_union(lsr_empty_set.value()).size() !=
            lsr_set_left.value().size() ||
        !lsr_set_left.value().intersection(lsr_empty_set.value()).empty() ||
        lsr_set_left.value().difference(lsr_empty_set.value()).size() !=
            lsr_set_left.value().size() ||
        !lsr_facade_contains || !lsr_facade_contains.value() ||
        !lsr_facade_not_contains || !lsr_facade_not_contains.value() ||
        !lsr_facade_union || lsr_facade_union.value().size() != 3 ||
        !lsr_facade_intersection ||
        lsr_facade_intersection.value().size() != 1 ||
        !lsr_facade_difference ||
        lsr_facade_difference.value().size() != 1 ||
        !lsr_facade_symmetric ||
        lsr_facade_symmetric.value().size() != 2 ||
        !lsr_facade_subset || !lsr_facade_subset.value() ||
        lsr_facade_null_membership ||
        std::string(LMCAS::error_name(
            lsr_facade_null_membership.error())) != "SetElementTypeMismatch") {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    auto lsr_domain_z = LMCAS::integers();
    auto lsr_domain_q = LMCAS::rationals();
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_expr = LMCAS::expressions();
    auto lsr_domain_chain_left =
        LMCAS::domain_subset(lsr_domain_z, lsr_domain_q);
    auto lsr_domain_chain_right =
        LMCAS::domain_subset(lsr_domain_r, lsr_domain_c);
    auto lsr_domain_chain_expr =
        LMCAS::domain_subset(lsr_domain_c, lsr_domain_expr);
    auto lsr_domain_reverse =
        LMCAS::domain_subset(lsr_domain_c, lsr_domain_r);
    auto lsr_domain_exact = LMCAS::domain_contains(
        lsr_domain_z, SymbolicExpr::number(2));
    auto lsr_domain_real = LMCAS::domain_contains(
        lsr_domain_r, SymbolicExpr::number(0.25));
    auto lsr_domain_complex = LMCAS::domain_contains(
        lsr_domain_c, i.value());
    auto lsr_domain_legacy_i = LMCAS::domain_contains(
        lsr_domain_c, SymbolicExpr::variable("i"));
    auto lsr_domain_legacy_i_not_real = LMCAS::domain_contains(
        lsr_domain_r, SymbolicExpr::variable("i"));
    auto lsr_legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto lsr_legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), lsr_legacy_four_i);
    auto lsr_domain_legacy_complex_arithmetic =
        LMCAS::domain_contains(lsr_domain_c,
                                     lsr_legacy_three_plus_four_i);
    auto lsr_domain_unknown = LMCAS::domain_contains(
        lsr_domain_r, x);
    auto lsr_domain_expr_symbol = LMCAS::domain_contains(
        lsr_domain_expr, x);
    auto lsr_set_subset_r = LMCAS::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_r);
    auto lsr_set_subset_c = LMCAS::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_c);
    auto lsr_complex_set = LMCAS::expr_set({i.value()});
    auto lsr_complex_set_subset_r =
        lsr_complex_set ? LMCAS::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_r)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    auto lsr_complex_set_subset_c =
        lsr_complex_set ? LMCAS::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_c)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    auto lsr_unknown_set = LMCAS::expr_set({x});
    auto lsr_legacy_i_set =
        LMCAS::expr_set({SymbolicExpr::variable("i")});
    auto lsr_legacy_i_set_subset_c =
        lsr_legacy_i_set ? LMCAS::expr_set_subset_domain(
                               lsr_legacy_i_set.value(), lsr_domain_c)
                         : LMCAS::Result<bool>::failure(
                               LMCAS::CasErrc::InternalInvariant,
                               "legacy i set construction failed", "consumer");
    auto lsr_legacy_complex_arithmetic_set =
        LMCAS::expr_set({lsr_legacy_three_plus_four_i});
    auto lsr_legacy_complex_arithmetic_set_subset_c =
        lsr_legacy_complex_arithmetic_set
            ? LMCAS::expr_set_subset_domain(
                  lsr_legacy_complex_arithmetic_set.value(), lsr_domain_c)
            : LMCAS::Result<bool>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "legacy complex arithmetic set construction failed",
                  "consumer");
    auto lsr_unknown_set_subset_r =
        lsr_unknown_set ? LMCAS::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_r)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "unknown set construction failed", "consumer");
    auto lsr_unknown_set_subset_expr =
        lsr_unknown_set ? LMCAS::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_expr)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "unknown set construction failed", "consumer");
    if (std::string(lsr_domain_z.name()) != "Z" ||
        std::string(lsr_domain_q.name()) != "Q" ||
        std::string(lsr_domain_r.name()) != "R" ||
        std::string(lsr_domain_c.name()) != "C" ||
        std::string(lsr_domain_expr.name()) != "Expr" ||
        !lsr_domain_chain_left || !lsr_domain_chain_left.value() ||
        !lsr_domain_chain_right || !lsr_domain_chain_right.value() ||
        !lsr_domain_chain_expr || !lsr_domain_chain_expr.value() ||
        !lsr_domain_reverse || lsr_domain_reverse.value() ||
        !lsr_domain_exact || !lsr_domain_exact.value() ||
        !lsr_domain_real || !lsr_domain_real.value() ||
        !lsr_domain_complex || !lsr_domain_complex.value() ||
        lsr_domain_legacy_i ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_i.error())) != "Inconclusive" ||
        lsr_domain_legacy_i_not_real ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_i_not_real.error())) != "Inconclusive" ||
        lsr_domain_legacy_complex_arithmetic ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_complex_arithmetic.error())) != "Inconclusive" ||
        !lsr_domain_expr_symbol || !lsr_domain_expr_symbol.value() ||
        !lsr_set_subset_r || !lsr_set_subset_r.value() ||
        !lsr_set_subset_c || !lsr_set_subset_c.value() ||
        !lsr_complex_set_subset_r || lsr_complex_set_subset_r.value() ||
        !lsr_complex_set_subset_c || !lsr_complex_set_subset_c.value() ||
        lsr_legacy_i_set_subset_c ||
        std::string(LMCAS::error_name(
            lsr_legacy_i_set_subset_c.error())) != "Inconclusive" ||
        lsr_legacy_complex_arithmetic_set_subset_c ||
        std::string(LMCAS::error_name(
            lsr_legacy_complex_arithmetic_set_subset_c.error())) !=
            "Inconclusive" ||
        !lsr_unknown_set_subset_expr ||
        !lsr_unknown_set_subset_expr.value() ||
        lsr_unknown_set_subset_r ||
        std::string(LMCAS::error_name(
            lsr_unknown_set_subset_r.error())) != "Inconclusive" ||
        lsr_domain_unknown ||
        std::string(LMCAS::error_name(lsr_domain_unknown.error())) !=
            "Inconclusive") {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    auto lsr_null_solve_set = LMCAS::solve_set(nullptr, "x");
    auto lsr_null_solve_expr_set = LMCAS::solve_expr_set(nullptr, "x");
    auto lsr_null_roots = LMCAS::roots(nullptr, "x");
    auto lsr_null_solve = LMCAS::solve(nullptr, "x");
    if (lsr_null_solve_set || lsr_null_solve_expr_set || lsr_null_roots ||
        lsr_null_solve ||
        std::string(LMCAS::error_name(lsr_null_solve_set.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(
            lsr_null_solve_expr_set.error())) != "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_roots.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_solve.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to reject null LMCAS set solve inputs\n";
        return 12;
    }
    auto lsr_empty_solve_set_variable = LMCAS::solve_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_variable = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_roots_variable = LMCAS::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_solve_variable = LMCAS::solve(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    if (lsr_empty_solve_set_variable || lsr_empty_variable || lsr_empty_roots_variable ||
        lsr_empty_solve_variable ||
        std::string(LMCAS::error_name(
            lsr_empty_solve_set_variable.error())) != "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_empty_variable.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(
            lsr_empty_roots_variable.error())) != "InvalidArgument" ||
        std::string(LMCAS::error_name(
            lsr_empty_solve_variable.error())) != "InvalidArgument") {
        std::cerr << "failed to reject empty LMCAS set solve variables\n";
        return 12;
    }

    auto lsr_nonfinite_binding = LMCAS::evalf(
        *x, LMCAS::NumericBindings{{"x", INFINITY}});
    auto lsr_nonfinite_expression =
        LMCAS::evalf(*SymbolicExpr::infinity());
    if (lsr_nonfinite_binding ||
        std::string(LMCAS::error_name(
            lsr_nonfinite_binding.error())) != "NumericFailure" ||
        lsr_nonfinite_expression ||
        std::string(LMCAS::error_name(
            lsr_nonfinite_expression.error())) != "NumericFailure") {
        std::cerr << "failed to reject non-finite LMCAS evalf results\n";
        return 12;
    }
    LMCAS::ResourceLimits exhausted_evalf_limits;
    exhausted_evalf_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_evalf_context(exhausted_evalf_limits);
    auto lsr_exhausted_evalf =
        LMCAS::evalf(*x, LMCAS::NumericBindings{{"x", 1.0}},
                           exhausted_evalf_context);
    if (lsr_exhausted_evalf ||
        std::string(LMCAS::error_name(
            lsr_exhausted_evalf.error())) != "ResourceLimit") {
        std::cerr << "failed to expose LMCAS evalf resource limits\n";
        return 12;
    }
    auto lsr_substituted = LMCAS::substitute(
        x_plus_one, "x", SymbolicExpr::number(4));
    auto lsr_substituted_value =
        lsr_substituted ? LMCAS::evalf(*lsr_substituted.value())
                        : LMCAS::Result<LMCAS::ApproxReal>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "substitution failed", "consumer");
    auto lsr_substitute_empty_var =
        LMCAS::substitute(x_plus_one, "", SymbolicExpr::number(4));
    auto lsr_substitute_null_value =
        LMCAS::substitute(x_plus_one, "x", nullptr);
    if (!lsr_substituted_value ||
        lsr_substituted_value.value().value != 5.0 ||
        lsr_substitute_empty_var || lsr_substitute_null_value ||
        std::string(LMCAS::error_name(
            lsr_substitute_empty_var.error())) != "InvalidArgument" ||
        std::string(LMCAS::error_name(
            lsr_substitute_null_value.error())) != "InvalidArgument") {
        std::cerr << "failed to expose LMCAS substitution facade\n";
        return 12;
    }

    auto lsr_match_pattern = SymbolicExpr::add(SymbolicExpr::variable("A"),
                                               SymbolicExpr::number(1));
    auto lsr_match_target = SymbolicExpr::add(SymbolicExpr::number(1), x);
    auto lsr_match =
        LMCAS::expr_match(lsr_match_pattern, lsr_match_target, {"A"});
    auto lsr_nonmatch = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("A")),
        SymbolicExpr::cos(x), {"A"});
    auto lsr_power_match = LMCAS::expr_match(
        SymbolicExpr::power(SymbolicExpr::variable("U"),
                            SymbolicExpr::variable("N")),
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        {"U", "N"});
    auto lsr_function_match = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("U")),
        SymbolicExpr::sin(SymbolicExpr::add(x, SymbolicExpr::number(1))),
        {"U"});
    auto lsr_repeated_match = LMCAS::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(x, x), {"A"});
    auto lsr_invalid_match =
        LMCAS::expr_match(lsr_match_pattern, lsr_match_target, {""});
    if (!lsr_match || !lsr_match.value().matched ||
        lsr_match.value().bindings.size() != 1 ||
        lsr_match.value().bindings[0].name != "A" ||
        !lsr_match.value().bindings[0].value ||
        !LMCAS::structurally_equal(
            *lsr_match.value().bindings[0].value, *x) ||
        !lsr_power_match || !lsr_power_match.value().matched ||
        lsr_power_match.value().bindings.size() != 2 ||
        !lsr_function_match || !lsr_function_match.value().matched ||
        lsr_function_match.value().bindings.size() != 1 ||
        !lsr_repeated_match || !lsr_repeated_match.value().matched ||
        !lsr_nonmatch || lsr_nonmatch.value().matched ||
        lsr_invalid_match ||
        std::string(LMCAS::error_name(
            lsr_invalid_match.error())) != "InvalidArgument") {
        std::cerr << "failed to expose LMCAS expression matching facade\n";
        return 12;
    }

    auto four_i = LMCAS::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = LMCAS::eval_complex(*three_plus_four_i);
    auto lowered_ordinary_i = LMCAS::eval_complex(*legacy_i);
    if (!lowered_complex || !lowered_complex.value().is_finite() ||
        lowered_complex.value().real.value != 3.0 ||
        lowered_complex.value().imag.value != 4.0) {
        std::cerr << "failed to explicitly lower LMCAS Expr to complex\n";
        return 13;
    }
    if (lowered_ordinary_i ||
        lowered_ordinary_i.error().code != LMCAS::CasErrc::UnboundSymbol) {
        std::cerr << "ordinary lowercase i complex lowering returned ";
        if (lowered_ordinary_i) {
            std::cerr << "a value\n";
        } else {
            std::cerr << LMCAS::error_name(lowered_ordinary_i.error())
                      << "\n";
        }
        return 13;
    }
    auto ordinary_multiply_complex = SymbolicExpr::add(
        SymbolicExpr::number(3),
        SymbolicExpr::multiply(SymbolicExpr::number(4), i.value()));
    auto lowered_ordinary_multiply =
        LMCAS::eval_complex(*ordinary_multiply_complex);
    if (!lowered_ordinary_multiply ||
        !lowered_ordinary_multiply.value().is_finite() ||
        lowered_ordinary_multiply.value().real.value != 3.0 ||
        lowered_ordinary_multiply.value().imag.value != 4.0) {
        std::cerr << "failed to lower LMCAS 3 + 4 * i form to complex\n";
        return 13;
    }
    auto zero_inverse_complex = LMCAS::eval_complex(
        *SymbolicExpr::power(SymbolicExpr::number(0), SymbolicExpr::number(-1)));
    if (zero_inverse_complex ||
        std::string(LMCAS::error_name(zero_inverse_complex.error())) !=
            "DomainError") {
        std::cerr << "failed to reject LMCAS complex reciprocal of zero\n";
        return 13;
    }
    auto nonfinite_complex_power = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), SymbolicExpr::infinity()));
    if (nonfinite_complex_power ||
        std::string(LMCAS::error_name(
            nonfinite_complex_power.error())) != "NumericFailure") {
        std::cerr << "failed to reject LMCAS complex non-finite exponent\n";
        return 13;
    }
    LMCAS::ResourceLimits exhausted_complex_limits;
    exhausted_complex_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_complex_context(exhausted_complex_limits);
    auto exhausted_complex = LMCAS::eval_complex(
        *three_plus_four_i, {}, exhausted_complex_context);
    if (exhausted_complex ||
        std::string(LMCAS::error_name(exhausted_complex.error())) !=
            "ResourceLimit") {
        std::cerr << "failed to expose LMCAS complex resource limits\n";
        return 13;
    }
    auto lsr_real = LMCAS::real(three_plus_four_i);
    auto lsr_imag = LMCAS::imag(three_plus_four_i);
    auto lsr_conj = LMCAS::conj(three_plus_four_i);
    auto lsr_abs = LMCAS::abs(three_plus_four_i);
    auto expected_conj = LMCAS::complex(SymbolicExpr::number(3),
                                              SymbolicExpr::number(-4));
    if (!lsr_real || !lsr_imag || !lsr_conj || !lsr_abs || !expected_conj ||
        !LMCAS::structurally_equal(*lsr_real.value(),
                                         *SymbolicExpr::number(3)) ||
        !LMCAS::structurally_equal(*lsr_imag.value(),
                                         *SymbolicExpr::number(4)) ||
        !LMCAS::structurally_equal(*lsr_conj.value(),
                                         *expected_conj.value())) {
        std::cerr << "failed to call LMCAS complex part facade\n";
        return 14;
    }
    auto lsr_null_real = LMCAS::real(nullptr);
    auto lsr_null_imag = LMCAS::imag(nullptr);
    auto lsr_null_conj = LMCAS::conj(nullptr);
    auto lsr_null_abs = LMCAS::abs(nullptr);
    if (lsr_null_real || lsr_null_imag || lsr_null_conj || lsr_null_abs ||
        std::string(LMCAS::error_name(lsr_null_real.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_imag.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_conj.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_abs.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to reject null LMCAS complex facade inputs\n";
        return 14;
    }
    auto lsr_real_value = SymbolicExpr::number(-5);
    auto lsr_real_value_real = LMCAS::real(lsr_real_value);
    auto lsr_real_value_imag = LMCAS::imag(lsr_real_value);
    auto lsr_real_value_conj = LMCAS::conj(lsr_real_value);
    auto lsr_real_value_abs = LMCAS::abs(lsr_real_value);
    auto lsr_real_value_abs_eval =
        lsr_real_value_abs ? LMCAS::evalf(*lsr_real_value_abs.value())
                           : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "abs(-5) construction failed", "consumer");
    if (!lsr_real_value_real || !lsr_real_value_imag ||
        !lsr_real_value_conj || !lsr_real_value_abs_eval ||
        !LMCAS::structurally_equal(*lsr_real_value_real.value(),
                                         *lsr_real_value) ||
        !LMCAS::structurally_equal(*lsr_real_value_imag.value(),
                                         *SymbolicExpr::number(0)) ||
        !LMCAS::structurally_equal(*lsr_real_value_conj.value(),
                                         *lsr_real_value) ||
        lsr_real_value_abs_eval.value().value < 4.999999999999 ||
        lsr_real_value_abs_eval.value().value > 5.000000000001) {
        std::cerr << "failed to promote real values through LMCAS complex facade\n";
        return 14;
    }

    return 0;
}
