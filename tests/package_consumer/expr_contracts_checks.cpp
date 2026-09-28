#include "assumption_context.hpp"
#include "expr.hpp"
#include "value.hpp"
#include "symbolic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace LMCAS;

static int check_expression_inspection() {
    const auto x = SymbolicExpr::variable("inspection_variable");
    const auto one = SymbolicExpr::number(1);
    const auto sum = SymbolicExpr::add(x, one);
    const auto name = LMCAS::symbol_name(x);
    if (!name || *name != "inspection_variable") {
        std::cerr << "public symbol inspection lost the variable name\n";
        return 10;
    }
    if (LMCAS::symbol_name(nullptr) || LMCAS::symbol_name(sum)) {
        std::cerr << "public symbol inspection accepted a non-symbol\n";
        return 10;
    }
    if (LMCAS::relation_op(nullptr) || LMCAS::relation_op(sum)) {
        std::cerr << "public relation inspection accepted a non-relation\n";
        return 10;
    }
    for (const auto op : {RelationOp::EQ, RelationOp::NEQ, RelationOp::LT,
                          RelationOp::GT, RelationOp::LEQ, RelationOp::GEQ}) {
        const auto relation = LMCAS::relation(x, one, op);
        if (!relation || LMCAS::relation_op(relation.value()) != op ||
            LMCAS::symbol_name(relation.value())) {
            std::cerr << "public relation inspection lost the relation operator\n";
            return 10;
        }
    }
    return 0;
}

static int check_null_arithmetic() {
    auto expr_two = LMCAS::integer(2);
    auto expr_null = LMCAS::add(nullptr, expr_two.value());
    if (expr_null ||
        std::string(LMCAS::error_name(expr_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS Expr arithmetic wrappers\n";
        return 10;
    }
    return 0;
}

static int check_null_expansion() {
    auto lsr_transform_null = LMCAS::expand(nullptr);
    if (lsr_transform_null ||
        std::string(LMCAS::error_name(lsr_transform_null.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS Expr transform wrappers\n";
        return 10;
    }
    return 0;
}

static int check_null_math() {
    auto lsr_math_null = LMCAS::sin(nullptr);
    auto lsr_math_null_clamp = LMCAS::clamp(
        SymbolicExpr::number(1), nullptr, SymbolicExpr::number(2));
    if (lsr_math_null ||
        lsr_math_null_clamp ||
        std::string(LMCAS::error_name(lsr_math_null.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_math_null_clamp.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to expose LMCAS std.math Expr wrappers\n";
        return 10;
    }
    return 0;
}

static int check_valid_equivalence_budget() {
    LMCAS::EqvOptions valid_budget_options;
    if (!LMCAS::set_eqv_budget(valid_budget_options, 256, 64, 4)) {
        std::cerr << "failed to configure LMCAS equivalence budget\n";
        return 11;
    }
    auto chain = SymbolicExpr::variable("x");
    for (int depth = 0; depth < 8; ++depth) {
        chain = SymbolicExpr::sin(chain);
    }
    valid_budget_options.budget.max_rewrite_depth = 8;
    ComputationContext context;
    auto exhausted = equivalent_core(*chain, *chain, context, valid_budget_options);
    if (exhausted || exhausted.error().code != CasErrc::ResourceLimit ||
        std::string(error_name(exhausted.error())) != "EqvBudgetExceeded") {
        std::cerr << "failed to exhaust nonzero LMCAS rewrite depth budget\n";
        return 11;
    }
    return 0;
}

static int check_invalid_equivalence_budgets() {
    for (const EqvBudget& budget :
         {EqvBudget{0, 64, 4}, EqvBudget{256, 0, 4}, EqvBudget{256, 64, 0}}) {
        EqvOptions options;
        auto invalid = set_eqv_budget(options, budget.max_rewrite_steps,
                                      budget.max_rewrite_depth, budget.max_node_growth_factor);
        if (invalid || std::string(error_name(invalid.error())) != "EqvBudgetExceeded") {
            std::cerr << "failed to expose LMCAS equivalence budget setter diagnostics\n";
            return 11;
        }
    }
    return 0;
}

static int check_exhausted_core_equivalence() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
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
    return 0;
}


static int check_additive_identity() {
    auto x = SymbolicExpr::variable("x");
    auto x_plus_zero = SymbolicExpr::add(x, SymbolicExpr::number(0));
    LMCAS::ComputationContext identity_eqv_context;
    auto identity_eqv =
        LMCAS::equivalent_core(*x_plus_zero, *x, identity_eqv_context);
    if (!identity_eqv || !identity_eqv.value()) {
        std::cerr << "failed to prove LMCAS Core identity example\n";
        return 11;
    }
    return 0;
}

static int check_multiplicative_identities() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_polynomial_equivalence() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_trigonometric_core_equivalence() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_trigonometric_equivalence() {
    auto x = SymbolicExpr::variable("x");
    auto trig_identity = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        SymbolicExpr::power(SymbolicExpr::cos(x), SymbolicExpr::number(2)));
    LMCAS::EqvOptions trig_eqv_options;
    if (!LMCAS::set_eqv_profile(trig_eqv_options, "Trig-Basic")) {
        std::cerr << "failed to configure LMCAS Trig-Basic profile\n";
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
    return 0;
}


static int check_exponential_profile() {
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
    return 0;
}

static int check_logarithmic_profile() {
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
    LMCAS::ComputationContext ln_one_eqv_context;
    auto ln_one_eqv = LMCAS::equivalent_core(
        *SymbolicExpr::ln(SymbolicExpr::number(1)), *SymbolicExpr::number(0),
        ln_one_eqv_context, exp_log_eqv_options);
    if (!ln_one_eqv || !ln_one_eqv.value()) {
        std::cerr << "failed to prove LMCAS ExpLog-Basic ln(1) example\n";
        return 11;
    }
    return 0;
}

static int check_unproven_log_domain() {
    auto x = SymbolicExpr::variable("x");
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
    auto exp_ln_x = SymbolicExpr::exp(SymbolicExpr::ln(x));
    LMCAS::ComputationContext exp_ln_unproven_context;
    auto exp_ln_unproven = LMCAS::equivalent_core(
        *exp_ln_x, *x, exp_ln_unproven_context, exp_log_eqv_options);
    if (exp_ln_unproven ||
        exp_ln_unproven.error().code != LMCAS::CasErrc::Inconclusive) {
        std::cerr << "failed to preserve unproved logarithm domain\n";
        return 11;
    }
    return 0;
}

static int check_positive_log_domain() {
    auto x = SymbolicExpr::variable("x");
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
    auto exp_ln_x = SymbolicExpr::exp(SymbolicExpr::ln(x));
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
    return 0;
}

static int check_null_solvers() {
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
    return 0;
}

static int check_empty_solve_variables() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_nonfinite_evaluation() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_evaluation_budget() {
    auto x = SymbolicExpr::variable("x");
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
    return 0;
}

static int check_substitution() {
    auto x = SymbolicExpr::variable("x");
    auto x_plus_one = SymbolicExpr::add(x, SymbolicExpr::number(1));
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
    return 0;
}

static int check_commutative_matching() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_match_pattern = SymbolicExpr::add(SymbolicExpr::variable("A"),
                                               SymbolicExpr::number(1));
    auto lsr_match_target = SymbolicExpr::add(SymbolicExpr::number(1), x);
    auto lsr_match =
        LMCAS::expr_match(lsr_match_pattern, lsr_match_target, {"A"});
    if (!lsr_match ||
        !lsr_match.value().matched ||
        lsr_match.value().bindings.size() != 1) {
        std::cerr << "failed to expose LMCAS expression matching facade\n";
        return 12;
    }
    if (lsr_match.value().bindings[0].name != "A" ||
        !lsr_match.value().bindings[0].value ||
        !LMCAS::structurally_equal(
            *lsr_match.value().bindings[0].value, *x)) {
        std::cerr << "failed to expose LMCAS expression matching facade\n";
        return 12;
    }
    return 0;
}

static int check_power_function_matching() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_power_match = LMCAS::expr_match(
        SymbolicExpr::power(SymbolicExpr::variable("U"),
                            SymbolicExpr::variable("N")),
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        {"U", "N"});
    auto lsr_function_match = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("U")),
        SymbolicExpr::sin(SymbolicExpr::add(x, SymbolicExpr::number(1))),
        {"U"});
    if (!lsr_power_match ||
        !lsr_power_match.value().matched ||
        lsr_power_match.value().bindings.size() != 2 ||
        !lsr_function_match ||
        !lsr_function_match.value().matched ||
        lsr_function_match.value().bindings.size() != 1) {
        std::cerr << "failed to expose LMCAS expression matching facade\n";
        return 12;
    }
    return 0;
}

static int check_matching_boundaries() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_match_pattern = SymbolicExpr::add(SymbolicExpr::variable("A"),
                                               SymbolicExpr::number(1));
    auto lsr_match_target = SymbolicExpr::add(SymbolicExpr::number(1), x);
    auto lsr_nonmatch = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("A")),
        SymbolicExpr::cos(x), {"A"});
    auto lsr_repeated_match = LMCAS::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(x, x), {"A"});
    auto lsr_invalid_match =
        LMCAS::expr_match(lsr_match_pattern, lsr_match_target, {""});
    if (!lsr_repeated_match ||
        !lsr_repeated_match.value().matched ||
        !lsr_nonmatch ||
        lsr_nonmatch.value().matched ||
        lsr_invalid_match ||
        std::string(LMCAS::error_name(
            lsr_invalid_match.error())) != "InvalidArgument") {
        std::cerr << "failed to expose LMCAS expression matching facade\n";
        return 12;
    }
    return 0;
}

int run_expr_contracts_consumer_checks() {
    const auto checks = {
        check_expression_inspection,
        check_null_arithmetic,
        check_null_expansion,
        check_null_math,
        check_valid_equivalence_budget,
        check_invalid_equivalence_budgets,
        check_exhausted_core_equivalence,
        check_additive_identity,
        check_multiplicative_identities,
        check_polynomial_equivalence,
        check_trigonometric_core_equivalence,
        check_trigonometric_equivalence,
        check_exponential_profile,
        check_logarithmic_profile,
        check_unproven_log_domain,
        check_positive_log_domain,
        check_null_solvers,
        check_empty_solve_variables,
        check_nonfinite_evaluation,
        check_evaluation_budget,
        check_substitution,
        check_commutative_matching,
        check_power_function_matching,
        check_matching_boundaries,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
