#include "expr.hpp"
#include "assumption_context.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprEquivalence, CommutativeAdditiveIdentity) {
    auto x = LMCAS::sym("x");
    auto one_plus_x = SymbolicExpr::add(SymbolicExpr::number(1), x.value());
    auto x_plus_one = SymbolicExpr::add(x.value(), SymbolicExpr::number(1));
    LMCAS::ComputationContext equivalent_context;
    auto equivalent = LMCAS::equivalent_core(*one_plus_x,
                                             *x_plus_one,
                                             equivalent_context);
    EXPECT_TRUE((equivalent && equivalent.value())) << "equivalent_core proves normalized additive equality";
}

TEST(ExprEquivalence, ArithmeticIdentities) {
    auto x = LMCAS::sym("x");
    auto x_times_one = SymbolicExpr::multiply(x.value(), SymbolicExpr::number(1));
    auto x_times_zero = SymbolicExpr::multiply(x.value(), SymbolicExpr::number(0));
    auto x_minus_x = SymbolicExpr::add(
        x.value(), SymbolicExpr::multiply(SymbolicExpr::number(-1), x.value()));
    LMCAS::ComputationContext multiply_identity_context;
    auto multiply_identity = LMCAS::equivalent_core(
        *x_times_one, *x.value(), multiply_identity_context);
    LMCAS::ComputationContext multiply_zero_context;
    auto multiply_zero = LMCAS::equivalent_core(
        *x_times_zero, *SymbolicExpr::number(0), multiply_zero_context);
    LMCAS::ComputationContext subtract_self_context;
    auto subtract_self = LMCAS::equivalent_core(
        *x_minus_x, *SymbolicExpr::number(0), subtract_self_context);
    EXPECT_TRUE((multiply_identity && multiply_identity.value())) << "equivalent_core proves Core x * 1 identity";
    EXPECT_TRUE((multiply_zero && multiply_zero.value())) << "equivalent_core proves Core x * 0 identity";
    EXPECT_TRUE((subtract_self && subtract_self.value())) << "equivalent_core proves Core x - x identity";
}

TEST(ExprEquivalence, PolynomialIdentity) {
    auto x = LMCAS::sym("x");
    auto x_plus_one = SymbolicExpr::add(x.value(), SymbolicExpr::number(1));
    auto x_plus_one_squared = SymbolicExpr::power(x_plus_one, SymbolicExpr::number(2));
    auto x_squared = SymbolicExpr::power(x.value(), SymbolicExpr::number(2));
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x.value());
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(x_squared, two_x), SymbolicExpr::number(1));
    LMCAS::ComputationContext polynomial_eqv_context;
    auto polynomial_equivalent = LMCAS::equivalent_core(
        *x_plus_one_squared, *expanded_square, polynomial_eqv_context);
    EXPECT_TRUE((polynomial_equivalent && polynomial_equivalent.value())) << "equivalent_core proves the LMCAS polynomial Core example";
}

TEST(ExprEquivalence, RenamedPolynomialIdentity) {
    auto y = LMCAS::sym("y");
    auto y_plus_one = SymbolicExpr::add(y.value(), SymbolicExpr::number(1));
    auto y_plus_one_squared = SymbolicExpr::power(y_plus_one, SymbolicExpr::number(2));
    auto y_squared = SymbolicExpr::power(y.value(), SymbolicExpr::number(2));
    auto two_y = SymbolicExpr::multiply(SymbolicExpr::number(2), y.value());
    auto expanded_y_square = SymbolicExpr::add(
        SymbolicExpr::add(y_squared, two_y), SymbolicExpr::number(1));
    LMCAS::ComputationContext y_polynomial_eqv_context;
    auto y_polynomial_equivalent = LMCAS::equivalent_core(
        *y_plus_one_squared, *expanded_y_square, y_polynomial_eqv_context);
    EXPECT_TRUE((y_polynomial_equivalent && y_polynomial_equivalent.value())) << "equivalent_core polynomial proof is not hard-coded to x";
}

TEST(ExprEquivalence, InvalidRewriteBudgets) {
    LMCAS::EqvOptions no_budget;
    auto invalid_budget = LMCAS::set_eqv_budget(no_budget, 0, 64, 4);
    LMCAS::EqvOptions no_depth_budget;
    auto invalid_depth_budget =
        LMCAS::set_eqv_budget(no_depth_budget, 256, 0, 4);
    LMCAS::EqvOptions no_growth_budget;
    auto invalid_growth_budget =
        LMCAS::set_eqv_budget(no_growth_budget, 256, 64, 0);
    bool invalid_budget_resource = false;
    bool invalid_budget_name = false;
    if (!invalid_budget) {
        invalid_budget_resource =
            invalid_budget.error().code == LMCAS::CasErrc::ResourceLimit;
        invalid_budget_name =
            std::string(LMCAS::error_name(invalid_budget.error())) == "EqvBudgetExceeded";
    }
    bool invalid_depth_resource = false;
    bool invalid_depth_name = false;
    if (!invalid_depth_budget) {
        invalid_depth_resource =
            invalid_depth_budget.error().code == LMCAS::CasErrc::ResourceLimit;
        invalid_depth_name =
            std::string(LMCAS::error_name(invalid_depth_budget.error())) == "EqvBudgetExceeded";
    }
    bool invalid_growth_resource = false;
    bool invalid_growth_name = false;
    if (!invalid_growth_budget) {
        invalid_growth_resource =
            invalid_growth_budget.error().code == LMCAS::CasErrc::ResourceLimit;
        invalid_growth_name =
            std::string(LMCAS::error_name(invalid_growth_budget.error())) == "EqvBudgetExceeded";
    }
    EXPECT_TRUE((invalid_budget_resource)) << "set_eqv_budget rejects zero rewrite steps";
    EXPECT_TRUE((invalid_depth_resource)) << "set_eqv_budget rejects zero rewrite depth";
    EXPECT_TRUE((invalid_growth_resource)) << "set_eqv_budget rejects zero node growth factor";
    EXPECT_TRUE((invalid_budget_name)) << "set_eqv_budget exposes EqvBudgetExceeded for invalid budgets";
    EXPECT_TRUE((invalid_depth_name)) << "zero rewrite depth exposes EqvBudgetExceeded";
    EXPECT_TRUE((invalid_growth_name)) << "zero node growth factor exposes EqvBudgetExceeded";
}

void expect_rewrite_budget_boundary(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    EqvProfile profile, const char* profile_name) {
    SCOPED_TRACE(profile_name);
    EqvOptions limited_options;
    limited_options.profile = profile;
    limited_options.budget.max_rewrite_steps = 8;
    ComputationContext limited_context;
    auto limited = equivalent_core(
        *lhs, *rhs, limited_context, limited_options);
    ASSERT_FALSE(limited);
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(limited.error().operation, "LMCAS.equivalent_core");

    EqvOptions sufficient_options = limited_options;
    sufficient_options.budget.max_rewrite_steps = 256;
    ComputationContext sufficient_context;
    auto sufficient = equivalent_core(
        *lhs, *rhs, sufficient_context, sufficient_options);
    ASSERT_TRUE(sufficient);
    EXPECT_TRUE(sufficient.value());
}

TEST(ExprEquivalence, RewriteStepBudgetCountsBothProofProfiles) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto trig_sum = SymbolicExpr::number(0);
    for (int term = 0; term < 4; ++term) {
        trig_sum = SymbolicExpr::add(
            trig_sum,
            SymbolicExpr::add(
                SymbolicExpr::power(SymbolicExpr::sin(x), two),
                SymbolicExpr::power(SymbolicExpr::cos(x), two)));
    }
    expect_rewrite_budget_boundary(
        trig_sum, SymbolicExpr::number(4),
        EqvProfile::TrigBasic, "Trig-Basic");

    auto exp_log_sum = SymbolicExpr::number(0);
    for (int term = 0; term < 5; ++term) {
        exp_log_sum = SymbolicExpr::add(
            exp_log_sum,
            SymbolicExpr::exp(SymbolicExpr::number(0)));
    }
    expect_rewrite_budget_boundary(
        exp_log_sum, SymbolicExpr::number(5),
        EqvProfile::ExpLogBasic, "ExpLog-Basic");
}

TEST(ExprEquivalence, CancellationPrecedesLocalRewriteExhaustion) {
    auto x = SymbolicExpr::variable("x");
    auto identity = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        SymbolicExpr::power(SymbolicExpr::cos(x), SymbolicExpr::number(2)));
    for (auto profile : {EqvProfile::TrigBasic, EqvProfile::ExpLogBasic}) {
        EqvOptions options;
        options.profile = profile;
        options.budget.max_rewrite_steps = 1;
        CancellationToken cancellation;
        cancellation.cancel();
        ComputationContext context({}, cancellation);
        auto cancelled = equivalent_core(
            *identity, *SymbolicExpr::number(1), context, options);
        ASSERT_FALSE(cancelled);
        EXPECT_EQ(cancelled.error().code, CasErrc::Cancelled);
        EXPECT_EQ(cancelled.error().operation, "LMCAS.equivalent_core");
    }
}

TEST(ExprEquivalence, ComputationResourceErrorRetainsOperation) {
    auto x = SymbolicExpr::variable("x");
    auto identity = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)),
        SymbolicExpr::power(SymbolicExpr::cos(x), SymbolicExpr::number(2)));
    EqvOptions options;
    options.profile = EqvProfile::TrigBasic;
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);
    auto exhausted = equivalent_core(
        *identity, *SymbolicExpr::number(1), context, options);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(exhausted.error().operation, "LMCAS.equivalent_core");
}

std::shared_ptr<SymbolicExpr> binomial_product(std::size_t factors) {
    const char *names[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"};
    auto product = SymbolicExpr::add(
        SymbolicExpr::variable(names[0]), SymbolicExpr::variable(names[1]));
    for (std::size_t index = 1; index < factors; ++index) {
        auto sum = SymbolicExpr::add(SymbolicExpr::variable(names[2 * index]),
                                     SymbolicExpr::variable(names[2 * index + 1]));
        product = SymbolicExpr::multiply(product, sum);
    }
    return product;
}

TEST(ExprEquivalence, RuntimeRewriteDepth) {
    auto chain = SymbolicExpr::variable("x");
    for (int index = 0; index < 8; ++index) {
        chain = SymbolicExpr::sin(chain);
    }
    EqvOptions options;
    options.budget.max_rewrite_depth = 8;
    ComputationContext limited_context;
    auto limited = equivalent_core(*chain, *chain, limited_context, options);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "eight nested sines exceed rewrite depth eight";
    EXPECT_TRUE((!limited && std::string(error_name(limited.error())) == "EqvBudgetExceeded")) << "runtime depth exhaustion retains the equivalence error name";

    options.budget.max_rewrite_depth = 32;
    ComputationContext sufficient_context;
    auto sufficient = equivalent_core(*chain, *chain, sufficient_context, options);
    EXPECT_TRUE((sufficient && sufficient.value())) << "the same nested expression is equivalent with sufficient depth";
}

TEST(ExprEquivalence, RuntimeNodeGrowth) {
    // Product construction flattens the left side to 16 nodes; rhs contributes one.
    auto product = binomial_product(5);
    auto zero = SymbolicExpr::number(0);
    EqvOptions options;
    options.budget.max_node_growth_factor = 4;
    ComputationContext limited_context;
    auto limited = equivalent_core(*product, *zero, limited_context, options);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "five binomial factors exceed fourfold growth while expanding";

    ComputationContext facade_context;
    auto facade = equivalent(*product, *zero, facade_context, options);
    ASSERT_FALSE(facade);
    EXPECT_EQ(facade.error().code, CasErrc::ResourceLimit);

}

TEST(ExprEquivalence, FlattenedProductsFitTheNodeBudget) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto product = SymbolicExpr::multiply(SymbolicExpr::multiply(x, y), z);
    auto zero = SymbolicExpr::number(0);
    EqvOptions options;
    options.budget.max_node_growth_factor = 2;
    ComputationContext context;
    auto result = equivalent_core(*product, *zero, context, options);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive);
}


TEST(ExprEquivalence, SaturatedGrowthAndCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto identity = SymbolicExpr::multiply(SymbolicExpr::number(1), x);
    EqvOptions options;
    options.budget.max_node_growth_factor = std::numeric_limits<std::size_t>::max();
    ComputationContext context;
    auto result = equivalent_core(*identity, *x, context, options);
    EXPECT_TRUE((result && result.value())) << "SIZE_MAX growth saturates instead of wrapping the node allowance";

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = equivalent_core(*identity, *x, cancelled_context, options);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << "equivalence still returns cancellation with a saturated growth allowance";
}

TEST(ExprEquivalence, TrigonometricIdentity) {
    auto y = LMCAS::sym("y");
    auto sin_y_squared = SymbolicExpr::power(
        SymbolicExpr::sin(y.value()), SymbolicExpr::number(2));
    auto cos_y_squared = SymbolicExpr::power(
        SymbolicExpr::cos(y.value()), SymbolicExpr::number(2));
    auto trig_identity = SymbolicExpr::add(sin_y_squared, cos_y_squared);
    LMCAS::EqvOptions trig_profile;
    auto trig_configured = LMCAS::set_eqv_profile(trig_profile, "Trig-Basic");
    EXPECT_TRUE((trig_configured.has_value())) << "set_eqv_profile accepts the LMCAS Trig-Basic profile name";
    LMCAS::ComputationContext trig_profile_context;
    auto trig_equivalent = LMCAS::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), trig_profile_context,
        trig_profile);
    EXPECT_TRUE((trig_equivalent && trig_equivalent.value())) << "equivalent_core proves the LMCAS Trig-Basic sin^2+cos^2 identity";
}

TEST(ExprEquivalence, TrigonometricParity) {
    auto y = LMCAS::sym("y");
    LMCAS::EqvOptions trig_profile;
    auto configured = LMCAS::set_eqv_profile(trig_profile, "Trig-Basic");
    EXPECT_TRUE((configured.has_value())) << "configure trigonometric parity profile";
    auto negative_y = SymbolicExpr::multiply(SymbolicExpr::number(-1), y.value());
    auto sin_negative_y = SymbolicExpr::sin(negative_y);
    auto negative_sin_y = SymbolicExpr::multiply(
        SymbolicExpr::number(-1), SymbolicExpr::sin(y.value()));
    LMCAS::ComputationContext trig_sin_odd_context;
    auto sin_odd_equivalent = LMCAS::equivalent_core(
        *sin_negative_y, *negative_sin_y, trig_sin_odd_context, trig_profile);
    EXPECT_TRUE((sin_odd_equivalent && sin_odd_equivalent.value())) << "equivalent_core proves the LMCAS Trig-Basic sin odd identity";

    auto cos_negative_y = SymbolicExpr::cos(negative_y);
    auto cos_y = SymbolicExpr::cos(y.value());
    LMCAS::ComputationContext trig_cos_even_context;
    auto cos_even_equivalent = LMCAS::equivalent_core(
        *cos_negative_y, *cos_y, trig_cos_even_context, trig_profile);
    EXPECT_TRUE((cos_even_equivalent && cos_even_equivalent.value())) << "equivalent_core proves the LMCAS Trig-Basic cos even identity";
}

TEST(ExprEquivalence, TrigonometricProfileBoundary) {
    auto y = LMCAS::sym("y");
    auto sin_y_squared = SymbolicExpr::power(
        SymbolicExpr::sin(y.value()), SymbolicExpr::number(2));
    auto cos_y_squared = SymbolicExpr::power(
        SymbolicExpr::cos(y.value()), SymbolicExpr::number(2));
    auto trig_identity = SymbolicExpr::add(sin_y_squared, cos_y_squared);
    LMCAS::EqvOptions trig_profile;
    auto trig_configured = LMCAS::set_eqv_profile(trig_profile, "Trig-Basic");
    LMCAS::ComputationContext core_trig_context;
    auto core_trig_equivalent = LMCAS::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), core_trig_context);
    ASSERT_FALSE(core_trig_equivalent);
    EXPECT_EQ(core_trig_equivalent.error().code, CasErrc::Inconclusive);
    LMCAS::ComputationContext lsr_trig_profile_context;
    auto lsr_trig_equivalent = LMCAS::equivalent(
        *trig_identity, *SymbolicExpr::number(1), lsr_trig_profile_context,
        trig_profile);
    EXPECT_TRUE((lsr_trig_equivalent && lsr_trig_equivalent.value())) << "LMCAS equivalent proves enabled Trig-Basic rules";
}

TEST(ExprEquivalence, ExponentialLogarithmProfiles) {
    LMCAS::EqvOptions exp_log_profile;
    auto exp_log_configured =
        LMCAS::set_eqv_profile(exp_log_profile, "ExpLog-Basic");
    EXPECT_TRUE((exp_log_configured.has_value())) << "set_eqv_profile accepts the LMCAS ExpLog-Basic profile name";
    auto unsupported_profile =
        LMCAS::eqv_profile_from_name("Richardson-Complete");
    ASSERT_FALSE(unsupported_profile);
    const auto unsupported_error = unsupported_profile.error();
    EXPECT_EQ(unsupported_error.code, LMCAS::CasErrc::UnsupportedExpression)
        << "eqv_profile_from_name rejects unsupported profile names";
    const std::string unsupported_error_name =
        LMCAS::error_name(unsupported_error);
    EXPECT_EQ(unsupported_error_name, "EqvRuleDisabled")
        << "unsupported equivalence profile exposes EqvRuleDisabled";
    LMCAS::ComputationContext exp_log_profile_context;
    auto exp_zero = SymbolicExpr::exp(SymbolicExpr::number(0));
    auto exp_log_equivalent = LMCAS::equivalent_core(
        *exp_zero, *SymbolicExpr::number(1), exp_log_profile_context,
        exp_log_profile);
    EXPECT_TRUE((exp_log_equivalent && exp_log_equivalent.value())) << "equivalent_core proves the LMCAS ExpLog-Basic exp(0) identity";

    auto ln_one = SymbolicExpr::ln(SymbolicExpr::number(1));
    LMCAS::ComputationContext ln_one_context;
    auto ln_one_equivalent = LMCAS::equivalent_core(
        *ln_one, *SymbolicExpr::number(0), ln_one_context, exp_log_profile);
    EXPECT_TRUE((ln_one_equivalent && ln_one_equivalent.value())) << "equivalent_core proves the LMCAS ExpLog-Basic ln(1) identity";

    auto lambert = SymbolicExpr::lambertw(SymbolicExpr::number(1));
    auto negative_lambert = SymbolicExpr::multiply(
        SymbolicExpr::number(-1), lambert);
    auto exponential = SymbolicExpr::exp(negative_lambert);
    LMCAS::ComputationContext lambert_context;
    auto lambert_equivalent = LMCAS::equivalent_core(
        *exponential, *lambert, lambert_context, exp_log_profile);
    EXPECT_TRUE((lambert_equivalent && lambert_equivalent.value()))
        << "ExpLog-Basic proves the defining LambertW identity at one";
    auto lambert_residual = SymbolicExpr::add(exponential, negative_lambert);
    LMCAS::ComputationContext lambert_residual_context;
    auto lambert_residual_zero = LMCAS::equivalent_core(
        *lambert_residual, *SymbolicExpr::number(0),
        lambert_residual_context, exp_log_profile);
    EXPECT_TRUE((lambert_residual_zero && lambert_residual_zero.value()))
        << "ExpLog-Basic reduces the complete LambertW residual";
}

TEST(ExprEquivalence, ExponentialLogarithmDomainEvidence) {
    auto y = LMCAS::sym("y");
    LMCAS::EqvOptions exp_log_profile;
    auto configured = LMCAS::set_eqv_profile(exp_log_profile, "ExpLog-Basic");
    EXPECT_TRUE((configured.has_value())) << "configure domain-dependent equivalence profile";
    auto exp_ln_y = SymbolicExpr::exp(SymbolicExpr::ln(y.value()));
    LMCAS::ComputationContext exp_ln_unproven_context;
    auto exp_ln_unproven = LMCAS::equivalent_core(
        *exp_ln_y, *y.value(), exp_ln_unproven_context, exp_log_profile);
    ASSERT_FALSE(exp_ln_unproven);
    EXPECT_EQ(exp_ln_unproven.error().code, CasErrc::Inconclusive);

    auto positive_assumptions = std::make_shared<LMCAS::AssumptionContext>();
    auto positive_assumption =
        positive_assumptions->assume_sign("y", LMCAS::Sign::Positive);
    EXPECT_TRUE((positive_assumption.has_value())) << "positive equivalence assumption is accepted";
    LMCAS::ComputationContext exp_ln_positive_context;
    auto set_positive_assumptions =
        exp_ln_positive_context.set_assumptions(positive_assumptions);
    EXPECT_TRUE((set_positive_assumptions.has_value())) << "equivalence context accepts positive assumptions";
    auto exp_ln_positive = LMCAS::equivalent_core(
        *exp_ln_y, *y.value(), exp_ln_positive_context, exp_log_profile);
    EXPECT_TRUE((exp_ln_positive && exp_ln_positive.value())) << "equivalent_core proves ExpLog-Basic exp(ln(y)) for positive y";
}

TEST(ExprEquivalence, ProvesDistinctPolynomialAndPropagatesUnknownProof) {
    auto x = SymbolicExpr::variable("x");
    auto x_plus_one = SymbolicExpr::add(x, SymbolicExpr::number(1));
    ComputationContext context;
    auto distinct = equivalent(*x, *x_plus_one, context);
    ASSERT_TRUE(distinct);
    EXPECT_FALSE(distinct.value());

    auto same = equivalent(*x, *x, context);
    ASSERT_TRUE(same);
    EXPECT_TRUE(same.value());

    auto unknown = equivalent(*SymbolicExpr::sin(x), *x, context);
    ASSERT_FALSE(unknown);
    EXPECT_EQ(unknown.error().code, CasErrc::Inconclusive);
}

} // namespace
