#include "test_mixed_transcendental_support.hpp"
#include "numeric_evaluation.hpp"
#include "internal/mixed_transcendental_support.hpp"

TEST(MixedTranscendentalPipelineContracts, DivisionByZeroRegionReturnsOnlyValidRoots) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(-1)),
        SymbolicExpr::sin(x));
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(result.value().value.empty());
    for (const auto &root : result.value().value) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        EXPECT_TRUE(std::isfinite(*value));
        EXPECT_NE(*value, 0.0);
        EXPECT_LT(std::abs(1.0 / *value + std::sin(*value)), opts.tolerance);
    }
}

TEST(MixedTranscendentalPipelineContracts, LnDomainRootIsPositiveAndSatisfiesResidual) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::ln(x), x);
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().value.size(), 1U);
    auto root = test_numeric_eval(result.value().value.front());
    ASSERT_TRUE(root.has_value());
    EXPECT_GT(*root, 0.0);
    EXPECT_NEAR(*root, 0.5671432904097838, 1e-10);
    EXPECT_LT(std::abs(std::log(*root) + *root), opts.tolerance);
}

TEST(MixedTranscendentalPipelineContracts, TanSingularitiesDoNotBecomeRoots) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::tan(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(result.value().value.empty());
    for (const auto &candidate : result.value().value) {
        auto root = test_numeric_eval(candidate);
        ASSERT_TRUE(root.has_value());
        ASSERT_TRUE(std::isfinite(*root));
        EXPECT_LT(std::abs(std::tan(*root) - *root), opts.tolerance);
    }
}

TEST(MixedTranscendentalPipelineContracts, NullExpressionIsInvalid) {
    SolveOptions opts;
    opts.allow_numeric = true;
    auto result = solve_mixed_transcendental_checked(nullptr, "x", opts);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
}

TEST(MixedTranscendentalPipelineContracts, ExtremeExponentialRootRemainsFinite) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(-1e300));
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.has_search_interval = true;
    opts.search_lo = -100.0;
    opts.search_hi = 1000.0;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(result.value().value.empty());
    EXPECT_EQ(result.value().completeness, Completeness::Complete);
    bool found = false;
    const double expected = std::log(1e300);
    auto payload = SymbolicExpr::number(Rational("1e300"));
    EqvOptions proof_options;
    proof_options.profile = EqvProfile::ExpLogBasic;
    ComputationContext proof_context;
    for (const auto &candidate : result.value().value) {
        auto root = test_numeric_eval(candidate);
        ASSERT_TRUE(root.has_value());
        ASSERT_TRUE(std::isfinite(*root));
        found = found || std::abs(*root - expected) < 1e-8;
        auto preimage = equivalent(
            *SymbolicExpr::exp(candidate), *payload, proof_context, proof_options);
        ASSERT_TRUE(preimage) << preimage.error().message;
        EXPECT_TRUE(preimage.value());
    }
    EXPECT_TRUE(found);
}

TEST(MixedTranscendentalPipelineContracts, NegativeExponentialRetainsSymbolicPreimage) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::exp(x)),
        SymbolicExpr::number(1e300));
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.has_search_interval = true;
    opts.search_lo = -100.0;
    opts.search_hi = 1000.0;

    auto result = solve_mixed_transcendental_checked(expression, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().value.size(), 1U);
    EXPECT_EQ(result.value().completeness, Completeness::Complete);
    auto payload = SymbolicExpr::number(Rational("1e300"));
    EqvOptions proof_options;
    proof_options.profile = EqvProfile::ExpLogBasic;
    ComputationContext proof_context;
    auto preimage = equivalent(*SymbolicExpr::exp(result.value().value.front()),
                              *payload, proof_context, proof_options);
    ASSERT_TRUE(preimage) << preimage.error().message;
    EXPECT_TRUE(preimage.value());
    auto root = test_numeric_eval(result.value().value.front());
    ASSERT_TRUE(root.has_value());
    EXPECT_NEAR(*root, std::log(1e300), 1e-8);
}

TEST(MixedTranscendentalPipelineContracts, NoVariableDependenceReturnsEmpty) {
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(SymbolicExpr::number(3)),
        SymbolicExpr::number(5));
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(result.value().value.empty());
    EXPECT_EQ(
        result.value().completeness, Completeness::Complete);
}

TEST(MixedTranscendentalPipelineContracts, NoVariableDependenceConstantExpression) {
    SolveOptions opts;
    opts.allow_numeric = true;
    auto result = solve_mixed_transcendental_checked(
        SymbolicExpr::number(42), "x", opts);

    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(result.value().value.empty());
    EXPECT_EQ(
        result.value().completeness, Completeness::Complete);
}

TEST(MixedTranscendentalPipelineContracts, IdenticallyZeroEquationIsInconclusive) {
    SolveOptions opts;
    opts.allow_numeric = true;
    auto result = solve_mixed_transcendental_checked(
        SymbolicExpr::number(0), "x", opts);

    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(result.value().value.empty());
    EXPECT_EQ(
        result.value().completeness, Completeness::Inconclusive);
    EXPECT_FALSE(result.value().reason.empty());
}

TEST(MixedTranscendentalPipelineContracts, NoVariableDependenceOtherVariable) {
    auto y = SymbolicExpr::variable("y");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(y), y);
    SolveOptions opts;
    opts.allow_numeric = true;

    auto result = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(result.value().value.empty());
    EXPECT_EQ(result.value().completeness, Completeness::Inconclusive);
    EXPECT_FALSE(result.value().reason.empty());
}

TEST(MixedTranscendentalPipelineContracts, CheckedMixedCandidates) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::sin(x);
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.has_search_interval = true;
    opts.search_lo = -1.0;
    opts.search_hi = 1.0;
    opts.tolerance = 1e-8;
    auto candidates =
        solve_mixed_transcendental_checked(expression, "x", opts);
    ASSERT_TRUE(candidates) << candidates.error().message;
    EXPECT_EQ(candidates.value().completeness, Completeness::Inconclusive);
    for (const auto &root : candidates.value().value) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        ASSERT_TRUE(std::isfinite(*value));
        EXPECT_GE(*value, opts.search_lo);
        EXPECT_LE(*value, opts.search_hi);
        EXPECT_LE(std::fabs(std::sin(*value)), opts.tolerance);
    }
}

static void expect_cubic_real_candidate(const std::shared_ptr<SymbolicExpr> &expression) {
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.has_search_interval = true;
    opts.search_lo = 0.0;
    opts.search_hi = 2.0;
    opts.tolerance = 1e-10;

    auto result = solve_mixed_transcendental_checked(expression, "x", opts);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().value.size(), 1U);

    auto evaluated = evaluate_numeric(*result.value().value.front());
    ASSERT_TRUE(evaluated);
    ASSERT_TRUE(evaluated.value().is_finite());
    const auto value = evaluated.value().value;
    EXPECT_NEAR(value, std::cbrt(2.0), opts.tolerance);
    auto residual = eval_at(expression, "x", value);
    ASSERT_TRUE(residual.has_value());
    EXPECT_LE(std::abs(*residual), opts.tolerance);
}

TEST(MixedTranscendentalPipelineContracts, CheckedMixedExactCubicRoots) {
    auto x = SymbolicExpr::variable("x");
    auto polynomial = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::number(-2));
    expect_cubic_real_candidate(polynomial);
}

TEST(MixedTranscendentalPipelineContracts, CheckedMixedFactoredExactCubicRoots) {
    auto x = SymbolicExpr::variable("x");
    auto polynomial = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::number(-2));
    expect_cubic_real_candidate(
        SymbolicExpr::multiply(polynomial, SymbolicExpr::exp(x)));
}

TEST(MixedTranscendentalPipelineContracts, CheckedMixedInvalidTolerance) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::sin(x);
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.has_search_interval = true;
    opts.search_lo = -1.0;
    opts.search_hi = 1.0;
    opts.tolerance = 1e-8;
    SolveOptions invalid = opts;
    invalid.tolerance = 0.0;
    auto invalid_result =
        solve_mixed_transcendental_checked(expression, "x", invalid);
    ASSERT_FALSE(invalid_result);
    EXPECT_EQ(invalid_result.error().code, CasErrc::InvalidArgument);
}

TEST(MixedTranscendentalPipelineContracts, CheckedMixedStepBudget) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::sin(x);
    SolveOptions opts;
    opts.allow_numeric = true;
    opts.has_search_interval = true;
    opts.search_lo = -1.0;
    opts.search_hi = 1.0;
    opts.tolerance = 1e-8;
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);
    auto limited = solve_mixed_transcendental_checked(
        expression, "x", opts, context);
    ASSERT_FALSE(limited);
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
}

TEST(MixedTranscendentalPipelineContracts, NumericOverflowPreservesOriginalError) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::exp(x), SymbolicExpr::multiply(x, SymbolicExpr::sin(x)));
    auto reference = evaluate_numeric(*SymbolicExpr::exp(SymbolicExpr::number(1000)));
    ASSERT_FALSE(reference);
    ASSERT_EQ(reference.error().code, CasErrc::NumericFailure);

    SolveOptions opts;
    opts.allow_numeric = true;
    opts.has_search_interval = true;
    opts.search_lo = 1000.0;
    opts.search_hi = 1001.0;
    ComputationContext context;
    auto solved = solve_mixed_transcendental_checked(expression, "x", opts, context);
    ASSERT_FALSE(solved);
    EXPECT_EQ(solved.error().code, reference.error().code);
    EXPECT_EQ(solved.error().message, reference.error().message);
    EXPECT_EQ(solved.error().operation, reference.error().operation);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(MixedTranscendentalPipelineContracts, ConstantNumericOverflowPreservesOriginalError) {
    auto expression = SymbolicExpr::exp(SymbolicExpr::number(1000));
    auto reference = evaluate_numeric(*expression);
    ASSERT_FALSE(reference);
    ASSERT_EQ(reference.error().code, CasErrc::NumericFailure);

    SolveOptions opts;
    opts.allow_numeric = true;
    ComputationContext context;
    auto solved = solve_mixed_transcendental_checked(expression, "x", opts, context);
    ASSERT_FALSE(solved);
    EXPECT_EQ(solved.error().code, reference.error().code);
    EXPECT_EQ(solved.error().message, reference.error().message);
    EXPECT_EQ(solved.error().operation, reference.error().operation);
}

TEST(MixedTranscendentalPipelineContracts, SamplesUseRemainingBudget) {
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);
    ASSERT_TRUE(context.consume_steps(1, "previous"));
    detail::RewriteBudget budget(context, limits.max_recursion_depth,
                                 limits.max_ast_nodes, "mixed.sample");
    auto x = SymbolicExpr::variable("x");
    try {
        (void)detail::mixed_evaluate_with_retry(x, "x", 0.0, 0.1, budget);
        FAIL() << "numeric sampling must not create a fresh context";
    } catch (const CasError& error) {
        EXPECT_EQ(error.code, CasErrc::ResourceLimit);
        EXPECT_EQ(error.operation, "mixed.sample");
    }
    EXPECT_EQ(context.recursion_depth(), 0u);
}
