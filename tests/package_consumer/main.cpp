#include "assumption_context.hpp"
#include "poly_utils.hpp"
#include "property_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "lsr_expr.hpp"
#include "solve_strategies.hpp"
#include "lmmc/lsr_stdlib.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(x, SymbolicExpr::number(1));
    if (!expr) {
        std::cerr << "failed to construct expression\n";
        return 1;
    }

    lamina::ComputationContext context;
    auto solved = lamina::solve_dispatch_checked(
        expr, "x", context, lamina::SolveOptions{});
    if (!solved || solved.value().kind() != lamina::SolutionSet::Kind::Finite ||
        solved.value().finite_solutions().size() != 1) {
        std::cerr << "failed to solve expression\n";
        return 2;
    }

    auto polynomial = lamina::symbolic_to_poly<Rational>(expr, "x");
    if (polynomial.coeffs.size() != 2 ||
        polynomial.coeffs[0] != Rational(1) ||
        polynomial.coeffs[1] != Rational(1)) {
        std::cerr << "failed to convert expression to polynomial\n";
        return 3;
    }

    auto interval = lamina::Interval::point(SymbolicExpr::number(0));
    lamina::ComputationContext interval_context;
    auto interval_union = lamina::IntervalUnion::from_intervals_checked(
        {interval}, interval_context);
    if (!interval_union || interval_union.value().intervals().size() != 1) {
        std::cerr << "failed to construct interval union\n";
        return 4;
    }

    lamina::PropertyStore properties;
    lamina::ComputationContext property_context;
    auto declared = properties.declare_continuous_checked(
        "f", interval, property_context);
    if (!declared) {
        std::cerr << "failed to declare property\n";
        return 5;
    }
    lamina::ComputationContext query_context;
    auto continuous = properties.is_continuous_checked(
        "f", interval, query_context);
    if (!continuous || !continuous.value()) {
        std::cerr << "failed to query property\n";
        return 6;
    }

    lamina::AssumptionContext assumptions;
    lamina::QueryInterface queries(assumptions);
    auto one = SymbolicExpr::number(1);
    auto positive = queries.query_positive_checked(*one);
    if (!positive || positive.value() != lamina::Tribool::True) {
        std::cerr << "failed to query positivity\n";
        return 7;
    }

    auto assumed_continuous = assumptions.current_properties().declare_continuous_checked(
        "g", interval);
    if (!assumed_continuous) {
        std::cerr << "failed to declare assumption property\n";
        return 8;
    }
    auto context_continuous = assumptions.is_continuous_checked("g", interval);
    if (!context_continuous || context_continuous.value() != lamina::Tribool::True) {
        std::cerr << "failed to query assumption property\n";
        return 9;
    }

    auto i = lamina::lsr::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LSR imaginary unit\n";
        return 10;
    }
    auto lsr_lower_i = lamina::lsr::i();
    auto lsr_upper_i = lamina::lsr::I();
    if (!lsr_lower_i || !lsr_upper_i ||
        !lamina::lsr::structurally_equal(*i.value(), *lsr_lower_i.value()) ||
        !lamina::lsr::structurally_equal(*i.value(), *lsr_upper_i.value())) {
        std::cerr << "failed to expose LSR imaginary unit aliases\n";
        return 10;
    }
    auto null_real_complex =
        lamina::lsr::complex(nullptr, SymbolicExpr::number(1));
    auto null_imag_complex =
        lamina::lsr::complex(SymbolicExpr::number(0), nullptr);
    if (null_real_complex || null_imag_complex ||
        std::string(lamina::lsr::error_name(null_real_complex.error())) !=
            "ComplexTypeMismatch" ||
        std::string(lamina::lsr::error_name(null_imag_complex.error())) !=
            "ComplexTypeMismatch") {
        std::cerr << "failed to expose LSR complex type diagnostics\n";
        return 10;
    }
    auto reserved_i = lamina::lsr::sym("i");
    if (reserved_i ||
        std::string(lamina::lsr::error_name(reserved_i.error())) !=
            "ImaginaryUnitReserved") {
        std::cerr << "failed to expose LSR diagnostic names\n";
        return 10;
    }
    auto reserved_I = lamina::lsr::sym("I");
    if (reserved_I ||
        std::string(lamina::lsr::error_name(reserved_I.error())) !=
            "ImaginaryUnitReserved") {
        std::cerr << "failed to expose LSR imaginary unit alias diagnostic\n";
        return 10;
    }
    auto reserved_phi = lamina::lsr::sym("phi");
    if (reserved_phi) {
        std::cerr << "failed to reserve LSR phi constant\n";
        return 10;
    }
    auto reserved_e = lamina::lsr::sym("e");
    auto reserved_unicode_pi = lamina::lsr::sym("\xCF\x80");
    if (reserved_e || reserved_unicode_pi) {
        std::cerr << "failed to reserve LSR std.math constants\n";
        return 10;
    }
    auto lsr_pi = lamina::lsr::pi();
    auto lsr_e = lamina::lsr::e();
    auto lsr_phi = lamina::lsr::phi();
    auto lsr_pi_value = lsr_pi ? lamina::lsr::evalf(*lsr_pi.value())
                               : lamina::Result<lamina::ApproxReal>::failure(
                                     lamina::CasErrc::InternalInvariant,
                                     "pi construction failed", "consumer");
    auto lsr_e_value = lsr_e ? lamina::lsr::evalf(*lsr_e.value())
                             : lamina::Result<lamina::ApproxReal>::failure(
                                   lamina::CasErrc::InternalInvariant,
                                   "e construction failed", "consumer");
    auto lsr_phi_value = lsr_phi ? lamina::lsr::evalf(*lsr_phi.value())
                                 : lamina::Result<lamina::ApproxReal>::failure(
                                       lamina::CasErrc::InternalInvariant,
                                       "phi construction failed", "consumer");
    auto lsr_unicode_pi_value = lamina::lsr::evalf(
        *SymbolicExpr::variable("\xCF\x80"));
    if (!lsr_pi_value || !lsr_e_value || !lsr_phi_value ||
        !lsr_unicode_pi_value ||
        std::abs(lsr_pi_value.value().value - 3.14159265358979323846) > 1e-15 ||
        std::abs(lsr_unicode_pi_value.value().value - 3.14159265358979323846) > 1e-15 ||
        std::abs(lsr_e_value.value().value - std::exp(1.0)) > 1e-15 ||
        std::abs(lsr_phi_value.value().value -
                 ((1.0 + std::sqrt(5.0)) / 2.0)) > 1e-15) {
        std::cerr << "failed to expose LSR std.math constants\n";
        return 10;
    }
    auto lsr_approx = lamina::lsr::approx_real(0.5);
    auto lsr_approx_value =
        lsr_approx ? lamina::lsr::evalf(*lsr_approx.value())
                   : lamina::Result<lamina::ApproxReal>::failure(
                         lamina::CasErrc::InternalInvariant,
                         "approx_real construction failed", "consumer");
    auto lsr_nan_approx = lamina::lsr::approx_real(NAN);
    auto lsr_inf_approx = lamina::lsr::approx_real(INFINITY);
    if (!lsr_approx_value || lsr_approx_value.value().value != 0.5 ||
        lsr_nan_approx || lsr_inf_approx ||
        std::string(lamina::lsr::error_name(lsr_nan_approx.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_inf_approx.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose finite LSR approx_real boundary\n";
        return 10;
    }
    auto lsr_expr_two = lamina::lsr::integer(2);
    auto lsr_expr_three = lamina::lsr::integer(3);
    auto lsr_expr_sum =
        lamina::lsr::add(lsr_expr_two.value(), lsr_expr_three.value());
    auto lsr_expr_product =
        lamina::lsr::mul(lsr_expr_sum.value(), SymbolicExpr::number(4));
    auto lsr_expr_quotient =
        lamina::lsr::div(lsr_expr_product.value(), SymbolicExpr::number(2));
    auto lsr_expr_difference =
        lamina::lsr::sub(lsr_expr_quotient.value(), SymbolicExpr::number(5));
    auto lsr_expr_negated = lamina::lsr::neg(lsr_expr_difference.value());
    auto lsr_expr_value =
        lsr_expr_negated ? lamina::lsr::evalf(*lsr_expr_negated.value())
                         : lamina::Result<lamina::ApproxReal>::failure(
                               lamina::CasErrc::InternalInvariant,
                               "Expr arithmetic construction failed",
                               "consumer");
    auto lsr_expr_polynomial =
        lamina::lsr::add(x, SymbolicExpr::number(1));
    auto lsr_expr_equation =
        lamina::lsr::eq(lsr_expr_polynomial.value(), SymbolicExpr::number(0));
    auto lsr_expr_solved =
        lsr_expr_polynomial
            ? lamina::lsr::solve_expr_set(lsr_expr_polynomial.value(), "x")
            : lamina::Result<lamina::lsr::ExprSet>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "Expr polynomial construction failed", "consumer");
    auto lsr_expr_null = lamina::lsr::add(nullptr, lsr_expr_two.value());
    if (!lsr_expr_value || std::abs(lsr_expr_value.value().value + 5.0) > 1e-15 ||
        !lsr_expr_solved ||
        !lsr_expr_equation || lsr_expr_equation.value()->to_string().empty() ||
        !lsr_expr_solved.value().contains(*SymbolicExpr::number(-1)) ||
        lsr_expr_null ||
        std::string(lamina::lsr::error_name(lsr_expr_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LSR Expr arithmetic wrappers\n";
        return 10;
    }
    auto lsr_transform_zero = lamina::lsr::integer(0);
    auto lsr_transform_one = lamina::lsr::integer(1);
    auto lsr_transform_two = lamina::lsr::integer(2);
    auto lsr_transform_three = lamina::lsr::integer(3);
    auto lsr_x_plus_zero =
        lamina::lsr::add(x, lsr_transform_zero.value());
    auto lsr_simplified =
        lsr_x_plus_zero ? lamina::lsr::simplify(lsr_x_plus_zero.value())
                        : lamina::Result<lamina::lsr::ExprPtr>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "simplify input failed", "consumer");
    auto lsr_left = lamina::lsr::add(x, lsr_transform_one.value());
    auto lsr_right = lamina::lsr::add(x, lsr_transform_two.value());
    auto lsr_product =
        lsr_left && lsr_right
            ? lamina::lsr::mul(lsr_left.value(), lsr_right.value())
            : lamina::Result<lamina::lsr::ExprPtr>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "expand input failed", "consumer");
    auto lsr_expanded =
        lsr_product ? lamina::lsr::expand(lsr_product.value())
                    : lamina::Result<lamina::lsr::ExprPtr>::failure(
                          lamina::CasErrc::InternalInvariant,
                          "expand input failed", "consumer");
    auto lsr_x_cubed = lamina::lsr::pow(x, lsr_transform_three.value());
    auto lsr_derivative =
        lsr_x_cubed ? lamina::lsr::differentiate(lsr_x_cubed.value(), "x")
                    : lamina::Result<lamina::lsr::ExprPtr>::failure(
                          lamina::CasErrc::InternalInvariant,
                          "differentiate input failed", "consumer");
    auto lsr_simplified_value =
        lsr_simplified ? lamina::lsr::evalf(*lsr_simplified.value(),
                                            lamina::NumericBindings{{"x", 7.0}})
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "simplify failed", "consumer");
    auto lsr_expanded_value =
        lsr_expanded ? lamina::lsr::evalf(*lsr_expanded.value(),
                                          lamina::NumericBindings{{"x", 3.0}})
                     : lamina::Result<lamina::ApproxReal>::failure(
                           lamina::CasErrc::InternalInvariant,
                           "expand failed", "consumer");
    auto lsr_derivative_value =
        lsr_derivative ? lamina::lsr::evalf(*lsr_derivative.value(),
                                            lamina::NumericBindings{{"x", 2.0}})
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "differentiate failed", "consumer");
    auto lsr_transform_null = lamina::lsr::expand(nullptr);
    if (!lsr_simplified_value || !lsr_expanded_value ||
        !lsr_derivative_value ||
        std::abs(lsr_simplified_value.value().value - 7.0) > 1e-15 ||
        std::abs(lsr_expanded_value.value().value - 20.0) > 1e-15 ||
        std::abs(lsr_derivative_value.value().value - 12.0) > 1e-15 ||
        lsr_transform_null ||
        std::string(lamina::lsr::error_name(lsr_transform_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LSR Expr transform wrappers\n";
        return 10;
    }
    auto lsr_math_sin = lamina::lsr::sin(SymbolicExpr::number(0));
    auto lsr_math_sqrt = lamina::lsr::sqrt(SymbolicExpr::number(9));
    auto lsr_math_pow =
        lamina::lsr::pow(SymbolicExpr::number(2), SymbolicExpr::number(4));
    auto lsr_math_asin = lamina::lsr::asin(SymbolicExpr::number(0.5));
    auto lsr_math_log10 = lamina::lsr::log10(SymbolicExpr::number(100));
    auto lsr_math_floor = lamina::lsr::floor(SymbolicExpr::number(2.75));
    auto lsr_math_ceil = lamina::lsr::ceil(SymbolicExpr::number(2.25));
    auto lsr_math_round = lamina::lsr::round(SymbolicExpr::number(-2.5));
    auto lsr_math_clamp = lamina::lsr::clamp(SymbolicExpr::number(7),
                                             SymbolicExpr::number(0),
                                             SymbolicExpr::number(5));
    auto lsr_math_sin_value =
        lsr_math_sin ? lamina::lsr::evalf(*lsr_math_sin.value())
                     : lamina::Result<lamina::ApproxReal>::failure(
                           lamina::CasErrc::InternalInvariant,
                           "sin construction failed", "consumer");
    auto lsr_math_sqrt_value =
        lsr_math_sqrt ? lamina::lsr::evalf(*lsr_math_sqrt.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "sqrt construction failed", "consumer");
    auto lsr_math_pow_value =
        lsr_math_pow ? lamina::lsr::evalf(*lsr_math_pow.value())
                     : lamina::Result<lamina::ApproxReal>::failure(
                           lamina::CasErrc::InternalInvariant,
                           "pow construction failed", "consumer");
    auto lsr_math_asin_value =
        lsr_math_asin ? lamina::lsr::evalf(*lsr_math_asin.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "asin construction failed", "consumer");
    auto lsr_math_log10_value =
        lsr_math_log10 ? lamina::lsr::evalf(*lsr_math_log10.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "log10 construction failed", "consumer");
    auto lsr_math_floor_value =
        lsr_math_floor ? lamina::lsr::evalf(*lsr_math_floor.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "floor construction failed", "consumer");
    auto lsr_math_ceil_value =
        lsr_math_ceil ? lamina::lsr::evalf(*lsr_math_ceil.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "ceil construction failed", "consumer");
    auto lsr_math_round_value =
        lsr_math_round ? lamina::lsr::evalf(*lsr_math_round.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "round construction failed", "consumer");
    auto lsr_math_clamp_value =
        lsr_math_clamp ? lamina::lsr::evalf(*lsr_math_clamp.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "clamp construction failed", "consumer");
    auto lsr_math_null = lamina::lsr::sin(nullptr);
    auto lsr_math_null_clamp = lamina::lsr::clamp(
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
        std::string(lamina::lsr::error_name(lsr_math_null.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_math_null_clamp.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LSR std.math Expr wrappers\n";
        return 10;
    }
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
    lamina::ComputationContext lsr_context;
    auto i_rule = lamina::lsr::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), lsr_context);
    if (!i_rule || !i_rule.value()) {
        std::cerr << "failed to prove LSR i*i == -1\n";
        return 11;
    }
    auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
    lamina::ComputationContext lsr_power_context;
    auto i_power_rule = lamina::lsr::equivalent_core(
        *i_power_two, *SymbolicExpr::number(-1), lsr_power_context);
    if (!i_power_rule || !i_power_rule.value()) {
        std::cerr << "failed to prove LSR i^2 == -1\n";
        return 11;
    }
    auto legacy_i = SymbolicExpr::variable("i");
    auto legacy_i_squared = SymbolicExpr::multiply(legacy_i, legacy_i);
    lamina::ComputationContext legacy_i_context;
    auto legacy_i_rule = lamina::lsr::equivalent_core(
        *legacy_i_squared, *SymbolicExpr::number(-1), legacy_i_context);
    if (!legacy_i_rule || !legacy_i_rule.value()) {
        std::cerr << "failed to normalize legacy Expr i in LSR equivalence\n";
        return 11;
    }
    lamina::lsr::EqvOptions valid_budget_options;
    if (!lamina::lsr::set_eqv_budget(valid_budget_options, 256, 64, 4)) {
        std::cerr << "failed to configure LSR equivalence budget\n";
        return 11;
    }
    lamina::lsr::EqvOptions invalid_budget_options;
    auto invalid_budget =
        lamina::lsr::set_eqv_budget(invalid_budget_options, 0, 64, 4);
    lamina::lsr::EqvOptions invalid_depth_budget_options;
    auto invalid_depth_budget =
        lamina::lsr::set_eqv_budget(invalid_depth_budget_options, 256, 0, 4);
    lamina::lsr::EqvOptions invalid_growth_budget_options;
    auto invalid_growth_budget =
        lamina::lsr::set_eqv_budget(invalid_growth_budget_options, 256, 64, 0);
    if (invalid_budget || invalid_depth_budget || invalid_growth_budget ||
        std::string(lamina::lsr::error_name(invalid_budget.error())) !=
            "EqvBudgetExceeded" ||
        std::string(lamina::lsr::error_name(invalid_depth_budget.error())) !=
            "EqvBudgetExceeded" ||
        std::string(lamina::lsr::error_name(invalid_growth_budget.error())) !=
            "EqvBudgetExceeded") {
        std::cerr << "failed to expose LSR equivalence budget setter diagnostics\n";
        return 11;
    }
    lamina::lsr::EqvOptions exhausted_eqv_options;
    (void)lamina::lsr::set_eqv_budget(exhausted_eqv_options, 1, 64, 4);
    exhausted_eqv_options.budget.max_rewrite_steps = 0;
    lamina::ComputationContext exhausted_eqv_context;
    auto exhausted_eqv = lamina::lsr::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), exhausted_eqv_context,
        exhausted_eqv_options);
    if (exhausted_eqv ||
        std::string(lamina::lsr::error_name(exhausted_eqv.error())) !=
            "EqvBudgetExceeded") {
        std::cerr << "failed to expose LSR equivalence budget diagnostics\n";
        return 11;
    }
    auto x_plus_zero = SymbolicExpr::add(x, SymbolicExpr::number(0));
    lamina::ComputationContext identity_eqv_context;
    auto identity_eqv =
        lamina::lsr::equivalent_core(*x_plus_zero, *x, identity_eqv_context);
    if (!identity_eqv || !identity_eqv.value()) {
        std::cerr << "failed to prove LSR Core identity example\n";
        return 11;
    }
    auto x_times_one = SymbolicExpr::multiply(x, SymbolicExpr::number(1));
    auto x_times_zero = SymbolicExpr::multiply(x, SymbolicExpr::number(0));
    auto x_minus_x = SymbolicExpr::add(
        x, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    lamina::ComputationContext multiply_identity_context;
    auto multiply_identity = lamina::lsr::equivalent_core(
        *x_times_one, *x, multiply_identity_context);
    lamina::ComputationContext multiply_zero_context;
    auto multiply_zero = lamina::lsr::equivalent_core(
        *x_times_zero, *SymbolicExpr::number(0), multiply_zero_context);
    lamina::ComputationContext subtract_self_context;
    auto subtract_self = lamina::lsr::equivalent_core(
        *x_minus_x, *SymbolicExpr::number(0), subtract_self_context);
    if (!multiply_identity || !multiply_identity.value() ||
        !multiply_zero || !multiply_zero.value() ||
        !subtract_self || !subtract_self.value()) {
        std::cerr << "failed to prove LSR Core algebra identities\n";
        return 11;
    }
    auto x_plus_one = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto x_plus_one_squared = SymbolicExpr::power(x_plus_one, SymbolicExpr::number(2));
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
        SymbolicExpr::number(1));
    lamina::ComputationContext polynomial_eqv_context;
    auto polynomial_eqv = lamina::lsr::equivalent_core(
        *x_plus_one_squared, *expanded_square, polynomial_eqv_context);
    if (!polynomial_eqv || !polynomial_eqv.value()) {
        std::cerr << "failed to prove LSR polynomial equivalence example\n";
        return 11;
    }
    auto trig_identity = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        SymbolicExpr::power(SymbolicExpr::cos(x), SymbolicExpr::number(2)));
    lamina::lsr::EqvOptions trig_eqv_options;
    if (!lamina::lsr::set_eqv_profile(trig_eqv_options, "Trig-Basic")) {
        std::cerr << "failed to configure LSR Trig-Basic profile\n";
        return 11;
    }
    lamina::ComputationContext trig_eqv_context;
    auto trig_eqv = lamina::lsr::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), trig_eqv_context,
        trig_eqv_options);
    if (!trig_eqv || !trig_eqv.value()) {
        std::cerr << "failed to prove LSR Trig-Basic equivalence example\n";
        return 11;
    }
    lamina::ComputationContext core_trig_eqv_context;
    auto core_trig_eqv = lamina::lsr::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), core_trig_eqv_context);
    if (!core_trig_eqv || core_trig_eqv.value()) {
        std::cerr << "enabled LSR Trig-Basic rules in Core profile\n";
        return 11;
    }
    lamina::lsr::EqvOptions exp_log_eqv_options;
    if (!lamina::lsr::set_eqv_profile(exp_log_eqv_options, "ExpLog-Basic")) {
        std::cerr << "failed to configure LSR ExpLog-Basic profile\n";
        return 11;
    }
    if (lamina::lsr::set_eqv_profile(exp_log_eqv_options,
                                     "Richardson-Complete")) {
        std::cerr << "accepted unsupported LSR equivalence profile\n";
        return 11;
    }
    lamina::ComputationContext exp_log_eqv_context;
    auto exp_log_eqv = lamina::lsr::equivalent_core(
        *SymbolicExpr::exp(SymbolicExpr::number(0)), *SymbolicExpr::number(1),
        exp_log_eqv_context, exp_log_eqv_options);
    if (!exp_log_eqv || !exp_log_eqv.value()) {
        std::cerr << "failed to prove LSR ExpLog-Basic equivalence example\n";
        return 11;
    }
    lamina::ComputationContext ln_one_eqv_context;
    auto ln_one_eqv = lamina::lsr::equivalent_core(
        *SymbolicExpr::ln(SymbolicExpr::number(1)), *SymbolicExpr::number(0),
        ln_one_eqv_context, exp_log_eqv_options);
    if (!ln_one_eqv || !ln_one_eqv.value()) {
        std::cerr << "failed to prove LSR ExpLog-Basic ln(1) example\n";
        return 11;
    }
    auto exp_ln_x = SymbolicExpr::exp(SymbolicExpr::ln(x));
    lamina::ComputationContext exp_ln_unproven_context;
    auto exp_ln_unproven = lamina::lsr::equivalent_core(
        *exp_ln_x, *x, exp_ln_unproven_context, exp_log_eqv_options);
    if (!exp_ln_unproven || exp_ln_unproven.value()) {
        std::cerr << "proved LSR exp(ln(x)) without domain evidence\n";
        return 11;
    }
    auto positive_assumptions = std::make_shared<lamina::AssumptionContext>();
    positive_assumptions->assume_sign("x", lamina::Sign::Positive);
    lamina::ComputationContext exp_ln_positive_context;
    if (!exp_ln_positive_context.set_assumptions(positive_assumptions)) {
        std::cerr << "failed to attach LSR equivalence assumptions\n";
        return 11;
    }
    auto exp_ln_positive = lamina::lsr::equivalent_core(
        *exp_ln_x, *x, exp_ln_positive_context, exp_log_eqv_options);
    if (!exp_ln_positive || !exp_ln_positive.value()) {
        std::cerr << "failed to prove LSR exp(ln(x)) under positive assumption\n";
        return 11;
    }

    auto lsr_roots = lamina::lsr::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_roots || lsr_roots.value().size() != 2) {
        std::cerr << "failed to lower LSR solve result to set<Expr>\n";
        return 12;
    }
    auto lsr_named_roots = lamina::lsr::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    auto lsr_named_solve = lamina::lsr::solve(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_named_roots || lsr_named_roots.value().size() != 2 ||
        !lsr_named_solve || lsr_named_solve.value().size() != 2) {
        std::cerr << "failed to call LSR roots/solve set<Expr> aliases\n";
        return 12;
    }
    auto lsr_repeated_roots = lamina::lsr::roots(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), "x");
    if (!lsr_repeated_roots || lsr_repeated_roots.value().size() != 1 ||
        !lsr_repeated_roots.value().contains(*SymbolicExpr::number(0))) {
        std::cerr << "failed to lower repeated LSR roots to set<Expr>\n";
        return 12;
    }
    auto lsr_cubic_roots = lamina::lsr::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(3)),
                          SymbolicExpr::number(-2)),
        "x");
    if (!lsr_cubic_roots || lsr_cubic_roots.value().size() != 3) {
        std::cerr << "failed to lower RootOf LSR roots to set<Expr>\n";
        return 12;
    }
    for (const auto& root : lsr_cubic_roots.value().elements()) {
        if (!root) {
            std::cerr << "failed to preserve LSR finite Expr roots\n";
            return 12;
        }
    }
    auto lsr_complex_roots = lamina::lsr::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(1)),
        "x");
    auto negative_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(-1));
    if (!lsr_complex_roots || lsr_complex_roots.value().size() != 2 ||
        !lsr_complex_roots.value().contains(*i.value()) ||
        !negative_i ||
        !lsr_complex_roots.value().contains(*negative_i.value())) {
        std::cerr << "failed to lower LSR complex roots to set<Expr>\n";
        return 12;
    }
    auto lsr_complex_roots_subset_c = lamina::lsr::expr_set_subset_domain(
        lsr_complex_roots.value(), lamina::lsr::complexes());
    if (!lsr_complex_roots_subset_c ||
        !lsr_complex_roots_subset_c.value()) {
        std::cerr << "failed to prove installed LSR complex roots subset C\n";
        return 12;
    }
    auto lsr_shifted_complex_roots = lamina::lsr::solve_expr_set(
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                              SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
            SymbolicExpr::number(2)),
        "x");
    if (!lsr_shifted_complex_roots ||
        lsr_shifted_complex_roots.value().size() != 2) {
        std::cerr << "failed to lower shifted LSR complex roots to set<Expr>\n";
        return 12;
    }
    bool saw_negative_one_plus_i = false;
    bool saw_negative_one_minus_i = false;
    for (const auto& root : lsr_shifted_complex_roots.value().elements()) {
        auto lowered_root = lamina::lsr::eval_complex(*root);
        if (!lowered_root || !lowered_root.value().is_finite()) {
            std::cerr << "failed to evaluate shifted LSR complex root\n";
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
        std::cerr << "failed to preserve shifted LSR complex root components\n";
        return 12;
    }
    auto lsr_set_left = lamina::lsr::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = lamina::lsr::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    if (!lsr_set_left || !lsr_set_right ||
        lsr_set_left.value().size() != 2 ||
        !lsr_set_left.value().contains(*SymbolicExpr::number(1))) {
        std::cerr << "failed to construct LSR set<Expr>\n";
        return 12;
    }
    auto lsr_empty_set = lamina::lsr::expr_set({});
    if (!lsr_empty_set || !lsr_empty_set.value().empty()) {
        std::cerr << "failed to construct empty LSR set<Expr>\n";
        return 12;
    }
    auto lsr_set_union = lsr_set_left.value().set_union(lsr_set_right.value());
    auto lsr_set_intersection =
        lsr_set_left.value().intersection(lsr_set_right.value());
    auto lsr_set_difference =
        lsr_set_left.value().difference(lsr_set_right.value());
    auto lsr_set_symmetric =
        lsr_set_left.value().symmetric_difference(lsr_set_right.value());
    auto lsr_facade_contains = lamina::lsr::expr_set_contains(
        lsr_set_left.value(), SymbolicExpr::number(1));
    auto lsr_facade_not_contains = lamina::lsr::expr_set_not_contains(
        lsr_set_left.value(), SymbolicExpr::number(3));
    auto lsr_facade_union = lamina::lsr::expr_set_union(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_intersection = lamina::lsr::expr_set_intersection(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_difference = lamina::lsr::expr_set_difference(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_symmetric = lamina::lsr::expr_set_symmetric_difference(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_subset = lamina::lsr::expr_set_subset(
        lsr_set_intersection, lsr_set_union);
    auto lsr_facade_null_membership = lamina::lsr::expr_set_contains(
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
        std::string(lamina::lsr::error_name(
            lsr_facade_null_membership.error())) != "SetElementTypeMismatch") {
        std::cerr << "failed to call LSR set<Expr> operations\n";
        return 12;
    }
    auto lsr_domain_z = lamina::lsr::integers();
    auto lsr_domain_q = lamina::lsr::rationals();
    auto lsr_domain_r = lamina::lsr::reals();
    auto lsr_domain_c = lamina::lsr::complexes();
    auto lsr_domain_expr = lamina::lsr::expressions();
    auto lsr_domain_chain_left =
        lamina::lsr::domain_subset(lsr_domain_z, lsr_domain_q);
    auto lsr_domain_chain_right =
        lamina::lsr::domain_subset(lsr_domain_r, lsr_domain_c);
    auto lsr_domain_chain_expr =
        lamina::lsr::domain_subset(lsr_domain_c, lsr_domain_expr);
    auto lsr_domain_reverse =
        lamina::lsr::domain_subset(lsr_domain_c, lsr_domain_r);
    auto lsr_domain_exact = lamina::lsr::domain_contains(
        lsr_domain_z, SymbolicExpr::number(2));
    auto lsr_domain_real = lamina::lsr::domain_contains(
        lsr_domain_r, SymbolicExpr::number(0.25));
    auto lsr_domain_complex = lamina::lsr::domain_contains(
        lsr_domain_c, i.value());
    auto lsr_domain_legacy_i = lamina::lsr::domain_contains(
        lsr_domain_c, SymbolicExpr::variable("i"));
    auto lsr_domain_legacy_i_not_real = lamina::lsr::domain_contains(
        lsr_domain_r, SymbolicExpr::variable("i"));
    auto lsr_legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto lsr_legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), lsr_legacy_four_i);
    auto lsr_domain_legacy_complex_arithmetic =
        lamina::lsr::domain_contains(lsr_domain_c,
                                     lsr_legacy_three_plus_four_i);
    auto lsr_domain_unknown = lamina::lsr::domain_contains(
        lsr_domain_r, x);
    auto lsr_domain_expr_symbol = lamina::lsr::domain_contains(
        lsr_domain_expr, x);
    auto lsr_set_subset_r = lamina::lsr::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_r);
    auto lsr_set_subset_c = lamina::lsr::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_c);
    auto lsr_complex_set = lamina::lsr::expr_set({i.value()});
    auto lsr_complex_set_subset_r =
        lsr_complex_set ? lamina::lsr::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_r)
                        : lamina::Result<bool>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    auto lsr_complex_set_subset_c =
        lsr_complex_set ? lamina::lsr::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_c)
                        : lamina::Result<bool>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    auto lsr_unknown_set = lamina::lsr::expr_set({x});
    auto lsr_legacy_i_set =
        lamina::lsr::expr_set({SymbolicExpr::variable("i")});
    auto lsr_legacy_i_set_subset_c =
        lsr_legacy_i_set ? lamina::lsr::expr_set_subset_domain(
                               lsr_legacy_i_set.value(), lsr_domain_c)
                         : lamina::Result<bool>::failure(
                               lamina::CasErrc::InternalInvariant,
                               "legacy i set construction failed", "consumer");
    auto lsr_legacy_complex_arithmetic_set =
        lamina::lsr::expr_set({lsr_legacy_three_plus_four_i});
    auto lsr_legacy_complex_arithmetic_set_subset_c =
        lsr_legacy_complex_arithmetic_set
            ? lamina::lsr::expr_set_subset_domain(
                  lsr_legacy_complex_arithmetic_set.value(), lsr_domain_c)
            : lamina::Result<bool>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "legacy complex arithmetic set construction failed",
                  "consumer");
    auto lsr_unknown_set_subset_r =
        lsr_unknown_set ? lamina::lsr::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_r)
                        : lamina::Result<bool>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "unknown set construction failed", "consumer");
    auto lsr_unknown_set_subset_expr =
        lsr_unknown_set ? lamina::lsr::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_expr)
                        : lamina::Result<bool>::failure(
                              lamina::CasErrc::InternalInvariant,
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
        !lsr_domain_legacy_i || !lsr_domain_legacy_i.value() ||
        !lsr_domain_legacy_i_not_real ||
        lsr_domain_legacy_i_not_real.value() ||
        !lsr_domain_legacy_complex_arithmetic ||
        !lsr_domain_legacy_complex_arithmetic.value() ||
        !lsr_domain_expr_symbol || !lsr_domain_expr_symbol.value() ||
        !lsr_set_subset_r || !lsr_set_subset_r.value() ||
        !lsr_set_subset_c || !lsr_set_subset_c.value() ||
        !lsr_complex_set_subset_r || lsr_complex_set_subset_r.value() ||
        !lsr_complex_set_subset_c || !lsr_complex_set_subset_c.value() ||
        !lsr_legacy_i_set_subset_c ||
        !lsr_legacy_i_set_subset_c.value() ||
        !lsr_legacy_complex_arithmetic_set_subset_c ||
        !lsr_legacy_complex_arithmetic_set_subset_c.value() ||
        !lsr_unknown_set_subset_expr ||
        !lsr_unknown_set_subset_expr.value() ||
        lsr_unknown_set_subset_r ||
        std::string(lamina::lsr::error_name(
            lsr_unknown_set_subset_r.error())) != "Inconclusive" ||
        lsr_domain_unknown ||
        std::string(lamina::lsr::error_name(lsr_domain_unknown.error())) !=
            "Inconclusive") {
        std::cerr << "failed to call LSR predefined number domain sets\n";
        return 12;
    }
    auto lsr_null_solve_set = lamina::lsr::solve_set(nullptr, "x");
    auto lsr_null_solve_expr_set = lamina::lsr::solve_expr_set(nullptr, "x");
    auto lsr_null_roots = lamina::lsr::roots(nullptr, "x");
    auto lsr_null_solve = lamina::lsr::solve(nullptr, "x");
    if (lsr_null_solve_set || lsr_null_solve_expr_set || lsr_null_roots ||
        lsr_null_solve ||
        std::string(lamina::lsr::error_name(lsr_null_solve_set.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(
            lsr_null_solve_expr_set.error())) != "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_null_roots.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_null_solve.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to reject null LSR set solve inputs\n";
        return 12;
    }
    auto lsr_empty_solve_set_variable = lamina::lsr::solve_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_variable = lamina::lsr::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_roots_variable = lamina::lsr::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    auto lsr_empty_solve_variable = lamina::lsr::solve(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "");
    if (lsr_empty_solve_set_variable || lsr_empty_variable || lsr_empty_roots_variable ||
        lsr_empty_solve_variable ||
        std::string(lamina::lsr::error_name(
            lsr_empty_solve_set_variable.error())) != "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_empty_variable.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(
            lsr_empty_roots_variable.error())) != "InvalidArgument" ||
        std::string(lamina::lsr::error_name(
            lsr_empty_solve_variable.error())) != "InvalidArgument") {
        std::cerr << "failed to reject empty LSR set solve variables\n";
        return 12;
    }

    auto lsr_nonfinite_binding = lamina::lsr::evalf(
        *x, lamina::NumericBindings{{"x", INFINITY}});
    auto lsr_nonfinite_expression =
        lamina::lsr::evalf(*SymbolicExpr::infinity());
    if (lsr_nonfinite_binding ||
        std::string(lamina::lsr::error_name(
            lsr_nonfinite_binding.error())) != "NumericFailure" ||
        lsr_nonfinite_expression ||
        std::string(lamina::lsr::error_name(
            lsr_nonfinite_expression.error())) != "NumericFailure") {
        std::cerr << "failed to reject non-finite LSR evalf results\n";
        return 12;
    }
    lamina::ResourceLimits exhausted_evalf_limits;
    exhausted_evalf_limits.max_steps = 0;
    lamina::ComputationContext exhausted_evalf_context(exhausted_evalf_limits);
    auto lsr_exhausted_evalf =
        lamina::lsr::evalf(*x, lamina::NumericBindings{{"x", 1.0}},
                           exhausted_evalf_context);
    if (lsr_exhausted_evalf ||
        std::string(lamina::lsr::error_name(
            lsr_exhausted_evalf.error())) != "ResourceLimit") {
        std::cerr << "failed to expose LSR evalf resource limits\n";
        return 12;
    }
    auto lsr_substituted = lamina::lsr::substitute(
        x_plus_one, "x", SymbolicExpr::number(4));
    auto lsr_substituted_value =
        lsr_substituted ? lamina::lsr::evalf(*lsr_substituted.value())
                        : lamina::Result<lamina::ApproxReal>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "substitution failed", "consumer");
    auto lsr_substitute_empty_var =
        lamina::lsr::substitute(x_plus_one, "", SymbolicExpr::number(4));
    auto lsr_substitute_null_value =
        lamina::lsr::substitute(x_plus_one, "x", nullptr);
    if (!lsr_substituted_value ||
        lsr_substituted_value.value().value != 5.0 ||
        lsr_substitute_empty_var || lsr_substitute_null_value ||
        std::string(lamina::lsr::error_name(
            lsr_substitute_empty_var.error())) != "InvalidArgument" ||
        std::string(lamina::lsr::error_name(
            lsr_substitute_null_value.error())) != "InvalidArgument") {
        std::cerr << "failed to expose LSR substitution facade\n";
        return 12;
    }

    auto lsr_match_pattern = SymbolicExpr::add(SymbolicExpr::variable("A"),
                                               SymbolicExpr::number(1));
    auto lsr_match_target = SymbolicExpr::add(SymbolicExpr::number(1), x);
    auto lsr_match =
        lamina::lsr::expr_match(lsr_match_pattern, lsr_match_target, {"A"});
    auto lsr_nonmatch = lamina::lsr::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("A")),
        SymbolicExpr::cos(x), {"A"});
    auto lsr_power_match = lamina::lsr::expr_match(
        SymbolicExpr::power(SymbolicExpr::variable("U"),
                            SymbolicExpr::variable("N")),
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        {"U", "N"});
    auto lsr_function_match = lamina::lsr::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("U")),
        SymbolicExpr::sin(SymbolicExpr::add(x, SymbolicExpr::number(1))),
        {"U"});
    auto lsr_repeated_match = lamina::lsr::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(x, x), {"A"});
    auto lsr_invalid_match =
        lamina::lsr::expr_match(lsr_match_pattern, lsr_match_target, {""});
    if (!lsr_match || !lsr_match.value().matched ||
        lsr_match.value().bindings.size() != 1 ||
        lsr_match.value().bindings[0].name != "A" ||
        !lsr_match.value().bindings[0].value ||
        !lamina::lsr::structurally_equal(
            *lsr_match.value().bindings[0].value, *x) ||
        !lsr_power_match || !lsr_power_match.value().matched ||
        lsr_power_match.value().bindings.size() != 2 ||
        !lsr_function_match || !lsr_function_match.value().matched ||
        lsr_function_match.value().bindings.size() != 1 ||
        !lsr_repeated_match || !lsr_repeated_match.value().matched ||
        !lsr_nonmatch || lsr_nonmatch.value().matched ||
        lsr_invalid_match ||
        std::string(lamina::lsr::error_name(
            lsr_invalid_match.error())) != "InvalidArgument") {
        std::cerr << "failed to expose LSR expression matching facade\n";
        return 12;
    }

    auto four_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = lamina::lsr::eval_complex(*three_plus_four_i);
    auto lowered_legacy_i = lamina::lsr::eval_complex(*legacy_i);
    if (!lowered_complex || !lowered_complex.value().is_finite() ||
        lowered_complex.value().real.value != 3.0 ||
        lowered_complex.value().imag.value != 4.0 ||
        !lowered_legacy_i || !lowered_legacy_i.value().is_finite() ||
        lowered_legacy_i.value().real.value != 0.0 ||
        lowered_legacy_i.value().imag.value != 1.0) {
        std::cerr << "failed to explicitly lower LSR Expr to complex\n";
        return 13;
    }
    auto ordinary_multiply_complex = SymbolicExpr::add(
        SymbolicExpr::number(3),
        SymbolicExpr::multiply(SymbolicExpr::number(4), i.value()));
    auto lowered_ordinary_multiply =
        lamina::lsr::eval_complex(*ordinary_multiply_complex);
    if (!lowered_ordinary_multiply ||
        !lowered_ordinary_multiply.value().is_finite() ||
        lowered_ordinary_multiply.value().real.value != 3.0 ||
        lowered_ordinary_multiply.value().imag.value != 4.0) {
        std::cerr << "failed to lower LSR 3 + 4 * i form to complex\n";
        return 13;
    }
    auto zero_inverse_complex = lamina::lsr::eval_complex(
        *SymbolicExpr::power(SymbolicExpr::number(0), SymbolicExpr::number(-1)));
    if (zero_inverse_complex ||
        std::string(lamina::lsr::error_name(zero_inverse_complex.error())) !=
            "DomainError") {
        std::cerr << "failed to reject LSR complex reciprocal of zero\n";
        return 13;
    }
    auto nonfinite_complex_power = lamina::lsr::eval_complex(
        *SymbolicExpr::power(i.value(), SymbolicExpr::infinity()));
    if (nonfinite_complex_power ||
        std::string(lamina::lsr::error_name(
            nonfinite_complex_power.error())) != "NumericFailure") {
        std::cerr << "failed to reject LSR complex non-finite exponent\n";
        return 13;
    }
    lamina::ResourceLimits exhausted_complex_limits;
    exhausted_complex_limits.max_steps = 0;
    lamina::ComputationContext exhausted_complex_context(exhausted_complex_limits);
    auto exhausted_complex = lamina::lsr::eval_complex(
        *three_plus_four_i, {}, exhausted_complex_context);
    if (exhausted_complex ||
        std::string(lamina::lsr::error_name(exhausted_complex.error())) !=
            "ResourceLimit") {
        std::cerr << "failed to expose LSR complex resource limits\n";
        return 13;
    }
    auto lsr_real = lamina::lsr::real(three_plus_four_i);
    auto lsr_imag = lamina::lsr::imag(three_plus_four_i);
    auto lsr_conj = lamina::lsr::conj(three_plus_four_i);
    auto lsr_abs = lamina::lsr::abs(three_plus_four_i);
    auto expected_conj = lamina::lsr::complex(SymbolicExpr::number(3),
                                              SymbolicExpr::number(-4));
    if (!lsr_real || !lsr_imag || !lsr_conj || !lsr_abs || !expected_conj ||
        !lamina::lsr::structurally_equal(*lsr_real.value(),
                                         *SymbolicExpr::number(3)) ||
        !lamina::lsr::structurally_equal(*lsr_imag.value(),
                                         *SymbolicExpr::number(4)) ||
        !lamina::lsr::structurally_equal(*lsr_conj.value(),
                                         *expected_conj.value())) {
        std::cerr << "failed to call LSR complex part facade\n";
        return 14;
    }
    auto lsr_null_real = lamina::lsr::real(nullptr);
    auto lsr_null_imag = lamina::lsr::imag(nullptr);
    auto lsr_null_conj = lamina::lsr::conj(nullptr);
    auto lsr_null_abs = lamina::lsr::abs(nullptr);
    if (lsr_null_real || lsr_null_imag || lsr_null_conj || lsr_null_abs ||
        std::string(lamina::lsr::error_name(lsr_null_real.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_null_imag.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_null_conj.error())) !=
            "InvalidArgument" ||
        std::string(lamina::lsr::error_name(lsr_null_abs.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to reject null LSR complex facade inputs\n";
        return 14;
    }
    auto lsr_real_value = SymbolicExpr::number(-5);
    auto lsr_real_value_real = lamina::lsr::real(lsr_real_value);
    auto lsr_real_value_imag = lamina::lsr::imag(lsr_real_value);
    auto lsr_real_value_conj = lamina::lsr::conj(lsr_real_value);
    auto lsr_real_value_abs = lamina::lsr::abs(lsr_real_value);
    auto lsr_real_value_abs_eval =
        lsr_real_value_abs ? lamina::lsr::evalf(*lsr_real_value_abs.value())
                           : lamina::Result<lamina::ApproxReal>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "abs(-5) construction failed", "consumer");
    if (!lsr_real_value_real || !lsr_real_value_imag ||
        !lsr_real_value_conj || !lsr_real_value_abs_eval ||
        !lamina::lsr::structurally_equal(*lsr_real_value_real.value(),
                                         *lsr_real_value) ||
        !lamina::lsr::structurally_equal(*lsr_real_value_imag.value(),
                                         *SymbolicExpr::number(0)) ||
        !lamina::lsr::structurally_equal(*lsr_real_value_conj.value(),
                                         *lsr_real_value) ||
        lsr_real_value_abs_eval.value().value < 4.999999999999 ||
        lsr_real_value_abs_eval.value().value > 5.000000000001) {
        std::cerr << "failed to promote real values through LSR complex facade\n";
        return 14;
    }

    lmmc_real_t lmmc_out = 0.0;
    lmmc_complex_t lmmc_z = {};
    lmmc_complex_t lmmc_w = {};
    uint64_t lmmc_hash_z = 0;
    uint64_t lmmc_hash_w = 0;
    int lmmc_complex_equal = 0;
    int lmmc_num_equal = 0;
    int lmmc_bool_equal = 0;
    int lmmc_text_equal = 0;
    const char* lmmc_constant_name = nullptr;
    const char* lmmc_constant_unit = nullptr;
    lmmc_real_t lmmc_constant_value = 0.0;
    if (lmmc_lsr_num_equal(2.0, 2.0, &lmmc_num_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_num_equal ||
        lmmc_lsr_num_hash(0.0, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_lsr_num_hash(-0.0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC LSR num key adapters\n";
        return 15;
    }
    {
        lmmc_real_t set_a_values[] = {1.0, 2.0, 2.0, -0.0};
        lmmc_real_t set_b_values[] = {2.0, 3.0, 0.0};
        lmmc_lsr_num_set_t set_a = {};
        lmmc_lsr_num_set_t set_b = {};
        lmmc_lsr_num_set_t set_union = {};
        lmmc_lsr_num_set_t set_difference = {};
        int contains = 0;
        int subset = 0;
        auto cleanup_sets = [&]() {
            lmmc_lsr_num_set_destroy(&set_a);
            lmmc_lsr_num_set_destroy(&set_b);
            lmmc_lsr_num_set_destroy(&set_union);
            lmmc_lsr_num_set_destroy(&set_difference);
        };
        if (lmmc_lsr_num_set_make(set_a_values, 4, &set_a) !=
                LMMC_STATUS_OK ||
            lmmc_lsr_num_set_make(set_b_values, 3, &set_b) !=
                LMMC_STATUS_OK ||
            set_a.size != 3 ||
            lmmc_lsr_num_set_contains(&set_a, 0.0, &contains) !=
                LMMC_STATUS_OK ||
            !contains ||
            lmmc_lsr_num_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 4 ||
            lmmc_lsr_num_set_difference(&set_union,
                                        &set_b,
                                        &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_lsr_num_set_subset(&set_difference, &set_union, &subset) !=
                LMMC_STATUS_OK ||
            !subset ||
            lmmc_lsr_num_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr << "failed to call installed LMMC LSR num set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    if (lmmc_lsr_bool_equal(1, 1, &lmmc_bool_equal) != LMMC_STATUS_OK ||
        !lmmc_bool_equal ||
        lmmc_lsr_bool_hash(1, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_lsr_bool_hash(0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z == lmmc_hash_w ||
        lmmc_lsr_text_equal("alpha", "alpha", &lmmc_text_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_text_equal ||
        lmmc_lsr_text_hash("alpha", &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_lsr_text_hash("alpha", &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC LSR bool/text key adapters\n";
        return 15;
    }
    {
        const char* text_a_values[] = {"alpha", "beta", "alpha"};
        const char* text_b_values[] = {"beta", "gamma"};
        lmmc_lsr_text_set_t set_a = {};
        lmmc_lsr_text_set_t set_b = {};
        lmmc_lsr_text_set_t set_union = {};
        lmmc_lsr_text_set_t set_difference = {};
        int contains = 0;
        int subset = 0;
        auto cleanup_sets = [&]() {
            lmmc_lsr_text_set_destroy(&set_a);
            lmmc_lsr_text_set_destroy(&set_b);
            lmmc_lsr_text_set_destroy(&set_union);
            lmmc_lsr_text_set_destroy(&set_difference);
        };
        if (lmmc_lsr_text_set_make(text_a_values, 3, &set_a) !=
                LMMC_STATUS_OK ||
            lmmc_lsr_text_set_make(text_b_values, 2, &set_b) !=
                LMMC_STATUS_OK ||
            set_a.size != 2 ||
            lmmc_lsr_text_set_contains(&set_a, "alpha", &contains) !=
                LMMC_STATUS_OK ||
            !contains ||
            lmmc_lsr_text_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 3 ||
            lmmc_lsr_text_set_difference(&set_union,
                                         &set_b,
                                         &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_lsr_text_set_subset(&set_difference,
                                     &set_union,
                                     &subset) != LMMC_STATUS_OK ||
            !subset ||
            lmmc_lsr_text_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr << "failed to call installed LMMC LSR text set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    {
        int bool_a_values[] = {1, 1, 0};
        int bool_b_values[] = {0};
        lmmc_lsr_bool_set_t set_a = {};
        lmmc_lsr_bool_set_t set_b = {};
        lmmc_lsr_bool_set_t set_difference = {};
        int contains = 0;
        int subset = 0;
        auto cleanup_sets = [&]() {
            lmmc_lsr_bool_set_destroy(&set_a);
            lmmc_lsr_bool_set_destroy(&set_b);
            lmmc_lsr_bool_set_destroy(&set_difference);
        };
        if (lmmc_lsr_bool_set_make(bool_a_values, 3, &set_a) !=
                LMMC_STATUS_OK ||
            lmmc_lsr_bool_set_make(bool_b_values, 1, &set_b) !=
                LMMC_STATUS_OK ||
            set_a.size != 2 ||
            lmmc_lsr_bool_set_contains(&set_a, 1, &contains) !=
                LMMC_STATUS_OK ||
            !contains ||
            lmmc_lsr_bool_set_difference(&set_a,
                                         &set_b,
                                         &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_lsr_bool_set_subset(&set_difference,
                                     &set_a,
                                     &subset) != LMMC_STATUS_OK ||
            !subset ||
            lmmc_lsr_bool_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr << "failed to call installed LMMC LSR bool set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    if (lmmc_lsr_math_complex(3.0, 4.0, &lmmc_z) != LMMC_STATUS_OK ||
        lmmc_lsr_math_complex(3.0, 4.0, &lmmc_w) != LMMC_STATUS_OK ||
        lmmc_lsr_math_complex_equal(&lmmc_z, &lmmc_w,
                                    &lmmc_complex_equal) != LMMC_STATUS_OK ||
        !lmmc_complex_equal ||
        lmmc_lsr_math_complex_hash(&lmmc_z, &lmmc_hash_z) !=
            LMMC_STATUS_OK ||
        lmmc_lsr_math_complex_hash(&lmmc_w, &lmmc_hash_w) !=
            LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC LSR complex key adapters\n";
        return 15;
    }
    {
        lmmc_complex_t complex_a_values[] = {
            {1.0, 2.0}, {1.0, 2.0}, {-0.0, 0.0}};
        lmmc_complex_t complex_b_values[] = {{0.0, -0.0}, {3.0, 4.0}};
        lmmc_lsr_complex_set_t set_a = {};
        lmmc_lsr_complex_set_t set_b = {};
        lmmc_lsr_complex_set_t set_union = {};
        lmmc_lsr_complex_set_t set_difference = {};
        lmmc_complex_t zero_complex = {0.0, 0.0};
        int contains = 0;
        int subset = 0;
        auto cleanup_sets = [&]() {
            lmmc_lsr_complex_set_destroy(&set_a);
            lmmc_lsr_complex_set_destroy(&set_b);
            lmmc_lsr_complex_set_destroy(&set_union);
            lmmc_lsr_complex_set_destroy(&set_difference);
        };
        if (lmmc_lsr_complex_set_make(complex_a_values, 3, &set_a) !=
                LMMC_STATUS_OK ||
            lmmc_lsr_complex_set_make(complex_b_values, 2, &set_b) !=
                LMMC_STATUS_OK ||
            set_a.size != 2 ||
            lmmc_lsr_complex_set_contains(&set_a,
                                          &zero_complex,
                                          &contains) != LMMC_STATUS_OK ||
            !contains ||
            lmmc_lsr_complex_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 3 ||
            lmmc_lsr_complex_set_difference(&set_union,
                                            &set_b,
                                            &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_lsr_complex_set_subset(&set_difference,
                                        &set_union,
                                        &subset) != LMMC_STATUS_OK ||
            !subset ||
            lmmc_lsr_complex_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr
                << "failed to call installed LMMC LSR complex set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    {
        lmmc_real_t diagonal_values[] = {1.0, 2.0, 3.0};
        lmmc_vec_t diagonal = {3, diagonal_values, 0};
        lmmc_mat_t eye = {};
        lmmc_mat_t diag = {};
        auto cleanup_mats = [&]() {
            lmmc_mat_destroy(&eye);
            lmmc_mat_destroy(&diag);
        };
        if (lmmc_lsr_linalg_eye(3, &eye) != LMMC_STATUS_OK ||
            eye.rows != 3 || eye.cols != 3 ||
            std::abs(eye.data[0 * eye.stride + 0] - 1.0) > 1e-12 ||
            std::abs(eye.data[1 * eye.stride + 1] - 1.0) > 1e-12 ||
            std::abs(eye.data[0 * eye.stride + 2]) > 1e-12 ||
            lmmc_lsr_linalg_diag(&diagonal, &diag) != LMMC_STATUS_OK ||
            diag.rows != 3 || diag.cols != 3 ||
            std::abs(diag.data[0 * diag.stride + 0] - 1.0) > 1e-12 ||
            std::abs(diag.data[1 * diag.stride + 1] - 2.0) > 1e-12 ||
            std::abs(diag.data[2 * diag.stride + 2] - 3.0) > 1e-12 ||
            std::abs(diag.data[0 * diag.stride + 1]) > 1e-12 ||
            lmmc_lsr_linalg_eye(0, &eye) != LMMC_STATUS_INVALID_ARGUMENT ||
            lmmc_lsr_linalg_diag(nullptr, &diag) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_mats();
            std::cerr
                << "failed to call installed LMMC LSR matrix constructors\n";
            return 15;
        }
        cleanup_mats();
    }
    if (lmmc_lsr_constants_count() == 0 ||
        lmmc_lsr_constants_get("C", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.99792458e8) > 1e-3 ||
        std::string(lmmc_lsr_constants_unit("C")) != "m*s^-1" ||
        lmmc_lsr_constants_entry(7,
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) != LMMC_STATUS_OK ||
        std::string(lmmc_constant_name ? lmmc_constant_name : "") != "C" ||
        std::abs(lmmc_constant_value - 2.99792458e8) > 1e-3 ||
        std::string(lmmc_constant_unit ? lmmc_constant_unit : "") !=
            "m*s^-1" ||
        lmmc_lsr_constants_get("NO_SUCH_CONSTANT", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT ||
        lmmc_lsr_constants_unit("NO_SUCH_CONSTANT") != nullptr ||
        lmmc_lsr_constants_entry(lmmc_lsr_constants_count(),
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        std::cerr << "failed to call installed LMMC LSR constants adapters\n";
        return 15;
    }

    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t stats_scaled[] = {2.0, 4.0, 6.0, 8.0};
    if (lmmc_lsr_stats_mean(stats_values, 4, &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_lsr_stats_median(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_lsr_stats_quantile(stats_values, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_lsr_stats_cov(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 3.333333333333333) > 1e-12 ||
        lmmc_lsr_stats_corr(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 1.0) > 1e-12 ||
        lmmc_lsr_stats_normal_pdf(0.0, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.3989422804014327) > 1e-12 ||
        lmmc_lsr_stats_binomial_pmf(2, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.375) > 1e-12 ||
        lmmc_lsr_stats_mean(stats_values, 0, &lmmc_out) !=
            LMMC_STATUS_EMPTY_INPUT) {
        std::cerr << "failed to call installed LMMC LSR stats adapters\n";
        return 15;
    }

    lmmc_rng_t* lmmc_rng = nullptr;
    lmmc_real_t lmmc_rand_a = 0.0;
    lmmc_real_t lmmc_rand_b = 0.0;
    int64_t lmmc_rand_int = 0;
    if (lmmc_rng_create(&lmmc_rng) != LMMC_STATUS_OK ||
        lmmc_lsr_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_lsr_random_rand(lmmc_rng, &lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_rand_a < 0.0 || lmmc_rand_a >= 1.0 ||
        lmmc_lsr_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_lsr_random_rand(lmmc_rng, &lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b ||
        lmmc_lsr_random_randint(lmmc_rng, 1, 3, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 1 || lmmc_rand_int > 3 ||
        lmmc_lsr_random_normal(lmmc_rng, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out) ||
        lmmc_lsr_random_choice(lmmc_rng, stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0)) {
        if (lmmc_rng) lmmc_rng_destroy(lmmc_rng);
        std::cerr << "failed to call installed LMMC LSR random adapters\n";
        return 15;
    }
    lmmc_rng_destroy(lmmc_rng);
    if (lmmc_lsr_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_lsr_random_default_rand(&lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_lsr_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_lsr_random_default_rand(&lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b ||
        lmmc_lsr_random_default_randint(2, 4, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 2 || lmmc_rand_int > 4 ||
        lmmc_lsr_random_default_normal(0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out) ||
        lmmc_lsr_random_default_choice(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0)) {
        lmmc_lsr_random_default_deinit();
        std::cerr << "failed to call installed LMMC LSR default random adapters\n";
        return 15;
    }
    lmmc_lsr_random_default_deinit();

    int lmmc_dimensionless = 0;
    if (lmmc_lsr_units_strip_num(10.0, "km", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10000.0) > 1e-12 ||
        lmmc_lsr_units_convert_from_si(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12 ||
        lmmc_lsr_units_convert_num(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12 ||
        lmmc_lsr_units_strip_scalar(10.0, &lmmc_out) != LMMC_STATUS_OK ||
        lmmc_out != 10.0 ||
        lmmc_lsr_units_is_dimensionless_num(10.0, &lmmc_dimensionless) !=
            LMMC_STATUS_OK ||
        !lmmc_dimensionless) {
        std::cerr << "failed to call installed LMMC LSR unit adapters\n";
        return 15;
    }

    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t cross = {0, nullptr, 0};
    lmmc_vec_t matvec = {0, nullptr, 0};
    lmmc_vec_t vec_add = {0, nullptr, 0};
    lmmc_vec_t vec_add_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_sub = {0, nullptr, 0};
    lmmc_vec_t vec_sub_scalar = {0, nullptr, 0};
    lmmc_vec_t scalar_sub_vec = {0, nullptr, 0};
    lmmc_vec_t vec_mul = {0, nullptr, 0};
    lmmc_vec_t vec_mul_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_div = {0, nullptr, 0};
    lmmc_vec_t vec_div_scalar = {0, nullptr, 0};
    lmmc_vec_t scalar_div_vec = {0, nullptr, 0};
    lmmc_vec_t vec_pow = {0, nullptr, 0};
    lmmc_vec_t vec_pow_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_scale = {0, nullptr, 0};
    lmmc_vec_t shape_vec = {0, nullptr, 0};
    lmmc_lsr_bool_vec_t vec_cmp = {0, nullptr, 0};
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_real_t short_values[] = {1.0, 2.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_vec_t short_vec = {2, short_values, 0};
    lmmc_mat_t matmul = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_add = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_add_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_sub = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_sub_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t scalar_sub_mat = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_mul_elem = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_mul_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_div = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_div_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t scalar_div_mat = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_pow_elem = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_pow_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_scale = {0, 0, 0, nullptr, 0};
    lmmc_lsr_eig_table_t eig_table = {};
    lmmc_lsr_svd_table_t svd_table = {};
    lmmc_lsr_bool_mat_t mat_cmp = {0, 0, 0, nullptr, 0};
    lmmc_real_t dot = 0.0;
    lmmc_real_t norm = 0.0;
    if (lmmc_lsr_linalg_dot(&a_vec, &b_vec, &dot) != LMMC_STATUS_OK ||
        dot != 32.0 ||
        lmmc_lsr_linalg_norm(&a_vec, &norm) != LMMC_STATUS_OK ||
        std::abs(norm - std::sqrt(14.0)) > 1e-12 ||
        lmmc_lsr_linalg_cross(&a_vec, &b_vec, &cross) != LMMC_STATUS_OK ||
        cross.size != 3 || !cross.data ||
        cross.data[0] != -3.0 || cross.data[1] != 6.0 ||
        cross.data[2] != -3.0 ||
        lmmc_lsr_linalg_vec_add(&a_vec, &b_vec, &vec_add) !=
            LMMC_STATUS_OK ||
        vec_add.size != 3 || !vec_add.data ||
        vec_add.data[0] != 5.0 || vec_add.data[1] != 7.0 ||
        vec_add.data[2] != 9.0 ||
        lmmc_lsr_linalg_vec_add_scalar(&a_vec, 10.0, &vec_add_scalar) !=
            LMMC_STATUS_OK ||
        vec_add_scalar.size != 3 || !vec_add_scalar.data ||
        vec_add_scalar.data[0] != 11.0 ||
        vec_add_scalar.data[1] != 12.0 ||
        vec_add_scalar.data[2] != 13.0 ||
        lmmc_lsr_linalg_vec_sub(&b_vec, &a_vec, &vec_sub) !=
            LMMC_STATUS_OK ||
        vec_sub.size != 3 || !vec_sub.data ||
        vec_sub.data[0] != 3.0 || vec_sub.data[1] != 3.0 ||
        vec_sub.data[2] != 3.0 ||
        lmmc_lsr_linalg_vec_sub_scalar(&a_vec, 1.0, &vec_sub_scalar) !=
            LMMC_STATUS_OK ||
        vec_sub_scalar.size != 3 || !vec_sub_scalar.data ||
        vec_sub_scalar.data[0] != 0.0 ||
        vec_sub_scalar.data[1] != 1.0 ||
        vec_sub_scalar.data[2] != 2.0 ||
        lmmc_lsr_linalg_scalar_sub_vec(10.0, &a_vec, &scalar_sub_vec) !=
            LMMC_STATUS_OK ||
        scalar_sub_vec.size != 3 || !scalar_sub_vec.data ||
        scalar_sub_vec.data[0] != 9.0 ||
        scalar_sub_vec.data[1] != 8.0 ||
        scalar_sub_vec.data[2] != 7.0 ||
        lmmc_lsr_linalg_vec_mul(&a_vec, &b_vec, &vec_mul) !=
            LMMC_STATUS_OK ||
        vec_mul.size != 3 || !vec_mul.data ||
        vec_mul.data[0] != 4.0 || vec_mul.data[1] != 10.0 ||
        vec_mul.data[2] != 18.0 ||
        lmmc_lsr_linalg_vec_mul_scalar(&a_vec, 2.0, &vec_mul_scalar) !=
            LMMC_STATUS_OK ||
        vec_mul_scalar.size != 3 || !vec_mul_scalar.data ||
        vec_mul_scalar.data[0] != 2.0 ||
        vec_mul_scalar.data[1] != 4.0 ||
        vec_mul_scalar.data[2] != 6.0 ||
        lmmc_lsr_linalg_vec_div(&b_vec, &a_vec, &vec_div) !=
            LMMC_STATUS_OK ||
        vec_div.size != 3 || !vec_div.data ||
        vec_div.data[0] != 4.0 || vec_div.data[1] != 2.5 ||
        vec_div.data[2] != 2.0 ||
        lmmc_lsr_linalg_vec_div_scalar(&b_vec, 2.0, &vec_div_scalar) !=
            LMMC_STATUS_OK ||
        vec_div_scalar.size != 3 || !vec_div_scalar.data ||
        vec_div_scalar.data[0] != 2.0 ||
        vec_div_scalar.data[1] != 2.5 ||
        vec_div_scalar.data[2] != 3.0 ||
        lmmc_lsr_linalg_scalar_div_vec(12.0, &a_vec, &scalar_div_vec) !=
            LMMC_STATUS_OK ||
        scalar_div_vec.size != 3 || !scalar_div_vec.data ||
        scalar_div_vec.data[0] != 12.0 ||
        scalar_div_vec.data[1] != 6.0 ||
        scalar_div_vec.data[2] != 4.0 ||
        lmmc_lsr_linalg_vec_pow(&a_vec, &short_vec, &vec_pow) !=
            LMMC_STATUS_DIMENSION_MISMATCH ||
        lmmc_lsr_linalg_vec_pow_scalar(&a_vec, 2.0, &vec_pow_scalar) !=
            LMMC_STATUS_OK ||
        vec_pow_scalar.size != 3 || !vec_pow_scalar.data ||
        vec_pow_scalar.data[0] != 1.0 ||
        vec_pow_scalar.data[1] != 4.0 ||
        vec_pow_scalar.data[2] != 9.0 ||
        lmmc_lsr_linalg_vec_compare_scalar(&a_vec,
                                           LMMC_LSR_COMPARE_GT,
                                           1.0,
                                           &vec_cmp) != LMMC_STATUS_OK ||
        vec_cmp.size != 3 || !vec_cmp.data ||
        vec_cmp.data[0] != 0 || vec_cmp.data[1] != 1 ||
        vec_cmp.data[2] != 1 ||
        lmmc_lsr_linalg_vec_scale(&a_vec, 2.0, &vec_scale) !=
            LMMC_STATUS_OK ||
        vec_scale.size != 3 || !vec_scale.data ||
        vec_scale.data[0] != 2.0 || vec_scale.data[1] != 4.0 ||
        vec_scale.data[2] != 6.0 ||
        lmmc_lsr_linalg_matvec(&matrix, &short_vec, &matvec) !=
            LMMC_STATUS_OK ||
        matvec.size != 2 || !matvec.data ||
        matvec.data[0] != 5.0 || matvec.data[1] != 11.0 ||
        lmmc_lsr_linalg_shape_vec(&matrix, &shape_vec) !=
            LMMC_STATUS_OK ||
        shape_vec.size != 2 || !shape_vec.data ||
        shape_vec.data[0] != 2.0 || shape_vec.data[1] != 2.0 ||
        lmmc_lsr_linalg_matmul(&matrix, &rhs_matrix, &matmul) !=
            LMMC_STATUS_OK ||
        matmul.rows != 2 || matmul.cols != 2 || !matmul.data ||
        matmul.data[0] != 19.0 || matmul.data[1] != 22.0 ||
        matmul.data[matmul.stride] != 43.0 ||
        matmul.data[matmul.stride + 1] != 50.0 ||
        lmmc_lsr_linalg_mat_add(&matrix, &rhs_matrix, &mat_add) !=
            LMMC_STATUS_OK ||
        mat_add.rows != 2 || mat_add.cols != 2 || !mat_add.data ||
        mat_add.data[0] != 6.0 || mat_add.data[1] != 8.0 ||
        mat_add.data[mat_add.stride] != 10.0 ||
        mat_add.data[mat_add.stride + 1] != 12.0 ||
        lmmc_lsr_linalg_mat_add_scalar(&matrix, 10.0, &mat_add_scalar) !=
            LMMC_STATUS_OK ||
        mat_add_scalar.rows != 2 || mat_add_scalar.cols != 2 ||
        !mat_add_scalar.data ||
        mat_add_scalar.data[0] != 11.0 ||
        mat_add_scalar.data[1] != 12.0 ||
        mat_add_scalar.data[mat_add_scalar.stride] != 13.0 ||
        mat_add_scalar.data[mat_add_scalar.stride + 1] != 14.0 ||
        lmmc_lsr_linalg_mat_sub(&rhs_matrix, &matrix, &mat_sub) !=
            LMMC_STATUS_OK ||
        mat_sub.rows != 2 || mat_sub.cols != 2 || !mat_sub.data ||
        mat_sub.data[0] != 4.0 || mat_sub.data[1] != 4.0 ||
        mat_sub.data[mat_sub.stride] != 4.0 ||
        mat_sub.data[mat_sub.stride + 1] != 4.0 ||
        lmmc_lsr_linalg_mat_sub_scalar(&matrix, 1.0, &mat_sub_scalar) !=
            LMMC_STATUS_OK ||
        mat_sub_scalar.rows != 2 || mat_sub_scalar.cols != 2 ||
        !mat_sub_scalar.data ||
        mat_sub_scalar.data[0] != 0.0 ||
        mat_sub_scalar.data[1] != 1.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride] != 2.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride + 1] != 3.0 ||
        lmmc_lsr_linalg_scalar_sub_mat(10.0, &matrix, &scalar_sub_mat) !=
            LMMC_STATUS_OK ||
        scalar_sub_mat.rows != 2 || scalar_sub_mat.cols != 2 ||
        !scalar_sub_mat.data ||
        scalar_sub_mat.data[0] != 9.0 ||
        scalar_sub_mat.data[1] != 8.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride] != 7.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride + 1] != 6.0 ||
        lmmc_lsr_linalg_mat_mul_elem(&matrix, &rhs_matrix, &mat_mul_elem) !=
            LMMC_STATUS_OK ||
        mat_mul_elem.rows != 2 || mat_mul_elem.cols != 2 ||
        !mat_mul_elem.data ||
        mat_mul_elem.data[0] != 5.0 || mat_mul_elem.data[1] != 12.0 ||
        mat_mul_elem.data[mat_mul_elem.stride] != 21.0 ||
        mat_mul_elem.data[mat_mul_elem.stride + 1] != 32.0 ||
        lmmc_lsr_linalg_mat_mul_scalar(&matrix, 2.0, &mat_mul_scalar) !=
            LMMC_STATUS_OK ||
        mat_mul_scalar.rows != 2 || mat_mul_scalar.cols != 2 ||
        !mat_mul_scalar.data ||
        mat_mul_scalar.data[0] != 2.0 ||
        mat_mul_scalar.data[1] != 4.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride] != 6.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride + 1] != 8.0 ||
        lmmc_lsr_linalg_mat_div(&rhs_matrix, &matrix, &mat_div) !=
            LMMC_STATUS_OK ||
        mat_div.rows != 2 || mat_div.cols != 2 || !mat_div.data ||
        mat_div.data[0] != 5.0 || mat_div.data[1] != 3.0 ||
        std::abs(mat_div.data[mat_div.stride] - 7.0 / 3.0) > 1e-12 ||
        mat_div.data[mat_div.stride + 1] != 2.0 ||
        lmmc_lsr_linalg_mat_div_scalar(&rhs_matrix, 2.0, &mat_div_scalar) !=
            LMMC_STATUS_OK ||
        mat_div_scalar.rows != 2 || mat_div_scalar.cols != 2 ||
        !mat_div_scalar.data ||
        mat_div_scalar.data[0] != 2.5 ||
        mat_div_scalar.data[1] != 3.0 ||
        mat_div_scalar.data[mat_div_scalar.stride] != 3.5 ||
        mat_div_scalar.data[mat_div_scalar.stride + 1] != 4.0 ||
        lmmc_lsr_linalg_scalar_div_mat(12.0, &matrix, &scalar_div_mat) !=
            LMMC_STATUS_OK ||
        scalar_div_mat.rows != 2 || scalar_div_mat.cols != 2 ||
        !scalar_div_mat.data ||
        scalar_div_mat.data[0] != 12.0 ||
        scalar_div_mat.data[1] != 6.0 ||
        scalar_div_mat.data[scalar_div_mat.stride] != 4.0 ||
        scalar_div_mat.data[scalar_div_mat.stride + 1] != 3.0 ||
        lmmc_lsr_linalg_mat_pow_elem(&matrix, &rhs_matrix, &mat_pow_elem) !=
            LMMC_STATUS_OK ||
        mat_pow_elem.rows != 2 || mat_pow_elem.cols != 2 ||
        !mat_pow_elem.data ||
        mat_pow_elem.data[0] != 1.0 ||
        mat_pow_elem.data[1] != 64.0 ||
        mat_pow_elem.data[mat_pow_elem.stride] != 2187.0 ||
        mat_pow_elem.data[mat_pow_elem.stride + 1] != 65536.0 ||
        lmmc_lsr_linalg_mat_pow_scalar(&matrix, 2.0, &mat_pow_scalar) !=
            LMMC_STATUS_OK ||
        mat_pow_scalar.rows != 2 || mat_pow_scalar.cols != 2 ||
        !mat_pow_scalar.data ||
        mat_pow_scalar.data[0] != 1.0 ||
        mat_pow_scalar.data[1] != 4.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride] != 9.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride + 1] != 16.0 ||
        lmmc_lsr_linalg_mat_compare_scalar(&matrix,
                                           LMMC_LSR_COMPARE_GE,
                                           3.0,
                                           &mat_cmp) != LMMC_STATUS_OK ||
        mat_cmp.rows != 2 || mat_cmp.cols != 2 || !mat_cmp.data ||
        mat_cmp.data[0] != 0 || mat_cmp.data[1] != 0 ||
        mat_cmp.data[mat_cmp.stride] != 1 ||
        mat_cmp.data[mat_cmp.stride + 1] != 1 ||
        lmmc_lsr_linalg_mat_scale(&matrix, 2.0, &mat_scale) !=
            LMMC_STATUS_OK ||
        mat_scale.rows != 2 || mat_scale.cols != 2 || !mat_scale.data ||
        mat_scale.data[0] != 2.0 || mat_scale.data[1] != 4.0 ||
        mat_scale.data[mat_scale.stride] != 6.0 ||
        mat_scale.data[mat_scale.stride + 1] != 8.0 ||
        lmmc_lsr_linalg_eig_table(&matrix, &eig_table) !=
            LMMC_STATUS_OK ||
        lmmc_lsr_eig_table_count(&eig_table) != 4 ||
        std::string(lmmc_lsr_eig_table_key(&eig_table, 0)
                        ? lmmc_lsr_eig_table_key(&eig_table, 0)
                        : "") != "values_real" ||
        std::string(lmmc_lsr_eig_table_key(&eig_table, 3)
                        ? lmmc_lsr_eig_table_key(&eig_table, 3)
                        : "") != "vectors_imag" ||
        lmmc_lsr_eig_table_key(&eig_table, 4) != nullptr ||
        lmmc_lsr_eig_table_get(&eig_table, "values_real") == nullptr ||
        lmmc_lsr_linalg_svd_table(&matrix, &svd_table) !=
            LMMC_STATUS_OK ||
        lmmc_lsr_svd_table_count(&svd_table) != 3 ||
        std::string(lmmc_lsr_svd_table_key(&svd_table, 0)
                        ? lmmc_lsr_svd_table_key(&svd_table, 0)
                        : "") != "U" ||
        std::string(lmmc_lsr_svd_table_key(&svd_table, 2)
                        ? lmmc_lsr_svd_table_key(&svd_table, 2)
                        : "") != "Vt" ||
        lmmc_lsr_svd_table_key(&svd_table, 3) != nullptr ||
        lmmc_lsr_svd_table_get(&svd_table, "S") == nullptr ||
        std::string(lmmc_lsr_error_name(LMMC_STATUS_DIMENSION_MISMATCH)) !=
            "DimensionMismatch" ||
        std::string(lmmc_lsr_error_name(
            LMMC_STATUS_UNIT_STRIP_INVALID)) != "UnitStripInvalid" ||
        lmmc_lsr_units_strip_num(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_INVALID ||
        lmmc_lsr_units_strip_num(1.0, "num<m>", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_LEGACY_SYNTAX ||
        lmmc_lsr_units_convert_from_si(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        lmmc_vec_destroy(&cross);
        lmmc_vec_destroy(&matvec);
        lmmc_vec_destroy(&vec_add);
        lmmc_vec_destroy(&vec_add_scalar);
        lmmc_vec_destroy(&vec_sub);
        lmmc_vec_destroy(&vec_sub_scalar);
        lmmc_vec_destroy(&scalar_sub_vec);
        lmmc_vec_destroy(&vec_mul);
        lmmc_vec_destroy(&vec_mul_scalar);
        lmmc_vec_destroy(&vec_div);
        lmmc_vec_destroy(&vec_div_scalar);
        lmmc_vec_destroy(&scalar_div_vec);
        lmmc_vec_destroy(&vec_pow);
        lmmc_vec_destroy(&vec_pow_scalar);
        lmmc_vec_destroy(&vec_scale);
        lmmc_vec_destroy(&shape_vec);
        lmmc_lsr_bool_vec_destroy(&vec_cmp);
        lmmc_mat_destroy(&matmul);
        lmmc_mat_destroy(&mat_add);
        lmmc_mat_destroy(&mat_add_scalar);
        lmmc_mat_destroy(&mat_sub);
        lmmc_mat_destroy(&mat_sub_scalar);
        lmmc_mat_destroy(&scalar_sub_mat);
        lmmc_mat_destroy(&mat_mul_elem);
        lmmc_mat_destroy(&mat_mul_scalar);
        lmmc_mat_destroy(&mat_div);
        lmmc_mat_destroy(&mat_div_scalar);
        lmmc_mat_destroy(&scalar_div_mat);
        lmmc_mat_destroy(&mat_pow_elem);
        lmmc_mat_destroy(&mat_pow_scalar);
        lmmc_mat_destroy(&mat_scale);
        lmmc_lsr_eig_table_destroy(&eig_table);
        lmmc_lsr_svd_table_destroy(&svd_table);
        lmmc_lsr_bool_mat_destroy(&mat_cmp);
        std::cerr << "failed to call installed LMMC LSR linalg adapters\n";
        return 15;
    }
    lmmc_vec_destroy(&cross);
    lmmc_vec_destroy(&matvec);
    lmmc_vec_destroy(&vec_add);
    lmmc_vec_destroy(&vec_add_scalar);
    lmmc_vec_destroy(&vec_sub);
    lmmc_vec_destroy(&vec_sub_scalar);
    lmmc_vec_destroy(&scalar_sub_vec);
    lmmc_vec_destroy(&vec_mul);
    lmmc_vec_destroy(&vec_mul_scalar);
    lmmc_vec_destroy(&vec_div);
    lmmc_vec_destroy(&vec_div_scalar);
    lmmc_vec_destroy(&scalar_div_vec);
    lmmc_vec_destroy(&vec_pow);
    lmmc_vec_destroy(&vec_pow_scalar);
    lmmc_vec_destroy(&vec_scale);
    lmmc_vec_destroy(&shape_vec);
    lmmc_lsr_bool_vec_destroy(&vec_cmp);
    lmmc_mat_destroy(&matmul);
    lmmc_mat_destroy(&mat_add);
    lmmc_mat_destroy(&mat_add_scalar);
    lmmc_mat_destroy(&mat_sub);
    lmmc_mat_destroy(&mat_sub_scalar);
    lmmc_mat_destroy(&scalar_sub_mat);
    lmmc_mat_destroy(&mat_mul_elem);
    lmmc_mat_destroy(&mat_mul_scalar);
    lmmc_mat_destroy(&mat_div);
    lmmc_mat_destroy(&mat_div_scalar);
    lmmc_mat_destroy(&scalar_div_mat);
    lmmc_mat_destroy(&mat_pow_elem);
    lmmc_mat_destroy(&mat_pow_scalar);
    lmmc_mat_destroy(&mat_scale);
    lmmc_lsr_eig_table_destroy(&eig_table);
    lmmc_lsr_svd_table_destroy(&svd_table);
    lmmc_lsr_bool_mat_destroy(&mat_cmp);

    std::cout << expr->to_string() << '\n';
    if (expr->to_string().empty()) {
        std::cerr << "expression string is empty\n";
        return 21;
    }
    return 0;
}
