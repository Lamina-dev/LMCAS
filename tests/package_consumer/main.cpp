#include "assumption_context.hpp"
#include "poly_utils.hpp"
#include "property_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "lsr_expr.hpp"
#include "solve_strategies.hpp"

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
            lsr_set_left.value().size()) {
        std::cerr << "failed to call LSR set<Expr> operations\n";
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

    auto four_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = lamina::lsr::eval_complex(*three_plus_four_i);
    if (!lowered_complex || !lowered_complex.value().is_finite() ||
        lowered_complex.value().real.value != 3.0 ||
        lowered_complex.value().imag.value != 4.0) {
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

    std::cout << expr->to_string() << '\n';
    if (expr->to_string().empty()) {
        std::cerr << "expression string is empty\n";
        return 21;
    }
    return 0;
}
