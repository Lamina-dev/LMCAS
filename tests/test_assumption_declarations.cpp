
#include "test_interval_support.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "internal/symbolic_ast.hpp"
#include <string>

using namespace LMCAS;

TEST(AssumptionDeclarations, ContinuousOnlyOnDifferentiableOverlapThrows) {
    PropertyStore store;

    // Declare differentiable on [0, 10]
    Interval diff_interval = make_closed_interval(0.0, 10.0);
    EXPECT_TRUE((store.declare_differentiable("f", diff_interval).has_value())) << "differentiability declaration succeeds";

    // Declaring continuous-only on [5, 15] (overlaps [0,10]) should throw
    Interval cont_interval = make_closed_interval(5.0, 15.0);
    auto failure_46 = store.declare_continuous("f", cont_interval);
    ASSERT_FALSE(failure_46.has_value()) << "Declaring continuous-only on differentiable overlap returns failure";
    EXPECT_TRUE((failure_46.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContinuousOnlyOnExactDifferentiableIntervalThrows) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 5.0);
    EXPECT_TRUE((store.declare_differentiable("g", iv).has_value())) << "differentiability declaration succeeds";

    auto failure_65 = store.declare_continuous("g", iv);
    EXPECT_TRUE((!failure_65.has_value())) << "Declaring continuous-only on exact differentiable interval returns failure";
}

TEST(AssumptionDeclarations, DifferentiableOnContinuousOverlapOk) {
    PropertyStore store;

    // Declare continuous-only on [0, 10]
    Interval cont_interval = make_closed_interval(0.0, 10.0);
    EXPECT_TRUE((store.declare_continuous("h", cont_interval).has_value())) << "continuity declaration succeeds";

    // Declaring differentiable on [5, 15] (overlaps) should NOT throw (it's an upgrade)
    Interval diff_interval = make_closed_interval(5.0, 15.0);
    auto success_76 = store.declare_differentiable("h", diff_interval);
    EXPECT_TRUE((success_76.has_value())) << "Declaring differentiable on continuous-only overlap does not throw";
}

TEST(AssumptionDeclarations, ContinuousOnContinuousOverlapOk) {
    PropertyStore store;

    Interval iv1 = make_closed_interval(0.0, 10.0);
    Interval iv2 = make_closed_interval(5.0, 15.0);
    EXPECT_TRUE((store.declare_continuous("k", iv1).has_value())) << "continuity declaration succeeds";

    auto success_93 = store.declare_continuous("k", iv2);
    EXPECT_TRUE((success_93.has_value())) << "Declaring continuous-only on continuous-only overlap does not throw";
}

TEST(AssumptionDeclarations, DifferentiableImpliesContinuous) {
    PropertyStore store;

    Interval iv = make_closed_interval(1.0, 5.0);
    EXPECT_TRUE((store.declare_differentiable("f", iv).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("f", iv).value())) << "Differentiable symbol is also continuous on same interval";
    EXPECT_TRUE((store.is_differentiable("f", iv).value())) << "Differentiable symbol is differentiable on same interval";
}

TEST(AssumptionDeclarations, ContinuousNotDifferentiable) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 3.0);
    EXPECT_TRUE((store.declare_continuous("f", iv).has_value())) << "continuity declaration succeeds";

    EXPECT_TRUE((store.is_continuous("f", iv).value())) << "Continuous symbol is continuous";
    EXPECT_FALSE((store.is_differentiable("f", iv).value())) << "Continuous-only symbol is NOT differentiable";
}

TEST(AssumptionDeclarations, NonOverlappingIntervalsNoConflict) {
    PropertyStore store;

    Interval diff_iv = make_closed_interval(0.0, 5.0);
    Interval cont_iv = make_closed_interval(6.0, 10.0);

    EXPECT_TRUE((store.declare_differentiable("f", diff_iv).has_value())) << "differentiability declaration succeeds";

    auto success_140 = store.declare_continuous("f", cont_iv);
    EXPECT_TRUE((success_140.has_value())) << "Disjoint intervals do not conflict";
}

TEST(AssumptionDeclarations, MonotonicityDeclareAndRetrieveIncreasing) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    EXPECT_TRUE((store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";

    Monotonicity result = store.get_monotonicity("f", "x", iv).value();
    EXPECT_TRUE((result == Monotonicity::Increasing)) << "get_monotonicity returns Increasing for exact interval";
}

TEST(AssumptionDeclarations, MonotonicityDeclareAndRetrieveDecreasing) {
    PropertyStore store;

    Interval iv = make_closed_interval(-5.0, 5.0);
    EXPECT_TRUE((store.declare_monotonicity("g", "t", iv, Monotonicity::Decreasing).has_value())) << "monotonicity declaration succeeds";

    Monotonicity result = store.get_monotonicity("g", "t", iv).value();
    EXPECT_TRUE((result == Monotonicity::Decreasing)) << "get_monotonicity returns Decreasing for exact interval";
}

TEST(AssumptionDeclarations, MonotonicitySubIntervalCoverage) {
    PropertyStore store;

    // Declare increasing on [0, 10]
    Interval large_iv = make_closed_interval(0.0, 10.0);
    EXPECT_TRUE((store.declare_monotonicity("f", "x", large_iv, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";

    // Query on [2, 8] (sub-interval) should return Increasing
    Interval sub_iv = make_closed_interval(2.0, 8.0);
    Monotonicity result = store.get_monotonicity("f", "x", sub_iv).value();
    EXPECT_TRUE((result == Monotonicity::Increasing)) << "Sub-interval query returns Increasing (covered by larger declaration)";
}

TEST(AssumptionDeclarations, MonotonicityUncoveredIntervalReturnsUnknown) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 5.0);
    EXPECT_TRUE((store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";

    // Query on [6, 10] (not covered) should return Unknown
    Interval uncovered = make_closed_interval(6.0, 10.0);
    Monotonicity result = store.get_monotonicity("f", "x", uncovered).value();
    EXPECT_TRUE((result == Monotonicity::Unknown)) << "Uncovered interval returns Unknown";
}

TEST(AssumptionDeclarations, MonotonicityWrongVariableReturnsUnknown) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    EXPECT_TRUE((store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";

    // Query with different variable
    Monotonicity result = store.get_monotonicity("f", "y", iv).value();
    EXPECT_TRUE((result == Monotonicity::Unknown)) << "Query with wrong variable returns Unknown";
}

TEST(AssumptionDeclarations, MonotonicityUndeclaredSymbolReturnsUnknown) {
    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity result = store.get_monotonicity("undeclared", "x", iv).value();
    EXPECT_TRUE((result == Monotonicity::Unknown)) << "Undeclared symbol returns Unknown";
}

TEST(AssumptionDeclarations, MonotonicityMultipleDeclarations) {
    PropertyStore store;

    Interval iv1 = make_closed_interval(0.0, 5.0);
    Interval iv2 = make_closed_interval(5.0, 10.0);

    EXPECT_TRUE((store.declare_monotonicity("f", "x", iv1, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";
    EXPECT_TRUE((store.declare_monotonicity("f", "x", iv2, Monotonicity::Decreasing).has_value())) << "monotonicity declaration succeeds";

    EXPECT_TRUE((store.get_monotonicity("f", "x", iv1).value() == Monotonicity::Increasing)) << "First interval returns Increasing";
    EXPECT_TRUE((store.get_monotonicity("f", "x", iv2).value() == Monotonicity::Decreasing)) << "Second interval returns Decreasing";
}

TEST(AssumptionDeclarations, MonotonicityEntireLine) {
    PropertyStore store;

    Interval entire = Interval::entire_line();
    EXPECT_TRUE((store.declare_monotonicity("exp", "x", entire, Monotonicity::Increasing).has_value())) << "monotonicity declaration succeeds";

    Interval sub = make_closed_interval(-100.0, 100.0);
    EXPECT_TRUE((store.get_monotonicity("exp", "x", sub).value() == Monotonicity::Increasing)) << "Entire line declaration covers any finite sub-interval";
}

TEST(AssumptionDeclarations, PeriodicityDeclareAndRetrieve) {
    PropertyStore store;

    auto period = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(6.283185307)));

    EXPECT_TRUE((store.declare_periodic("sin_x", "x", period).has_value())) << "period declaration succeeds";

    EXPECT_TRUE((store.is_periodic("sin_x", "x"))) << "Symbol is periodic after declaration";

    auto retrieved = store.get_period("sin_x", "x");
    EXPECT_TRUE((retrieved.has_value())) << "get_period returns a value";
}

TEST(AssumptionDeclarations, PeriodicityNotPeriodicByDefault) {
    PropertyStore store;

    EXPECT_FALSE((store.is_periodic("undeclared", "x"))) << "Undeclared symbol is not periodic";

    auto period = store.get_period("undeclared", "x");
    EXPECT_FALSE((period.has_value())) << "get_period returns nullopt for undeclared symbol";
}

TEST(AssumptionDeclarations, PeriodicityOverwritePeriod) {
    PropertyStore store;

    auto period1 = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(3.14159)));
    auto period2 = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(6.28318)));

    EXPECT_TRUE((store.declare_periodic("f", "x", period1).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((store.is_periodic("f", "x"))) << "f is periodic after first declaration";

    EXPECT_TRUE((store.declare_periodic("f", "x", period2).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((store.is_periodic("f", "x"))) << "f is still periodic after second declaration";

    auto retrieved = store.get_period("f", "x");
    ASSERT_TRUE((retrieved.has_value())) << "get_period returns a value after overwrite";
    if (retrieved.has_value()) {
        double val = retrieved->to_numeric();
        EXPECT_NEAR(val, 6.28318, 1e-4) << "Period is updated to new value";
    }
}

TEST(AssumptionDeclarations, PeriodicityDifferentSymbols) {
    PropertyStore store;

    auto period_sin = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(6.28318)));
    auto period_tan = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(3.14159)));

    EXPECT_TRUE((store.declare_periodic("sin_x", "x", period_sin).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((store.declare_periodic("tan_x", "x", period_tan).has_value())) << "period declaration succeeds";

    EXPECT_TRUE((store.is_periodic("sin_x", "x"))) << "sin_x is periodic";
    EXPECT_TRUE((store.is_periodic("tan_x", "x"))) << "tan_x is periodic";

    auto p_sin = store.get_period("sin_x", "x");
    auto p_tan = store.get_period("tan_x", "x");

    ASSERT_TRUE((p_sin.has_value())) << "sin_x has a period";
    EXPECT_TRUE((p_tan.has_value())) << "tan_x has a period";

    if (p_sin.has_value() && p_tan.has_value()) {
        EXPECT_NEAR(p_sin->to_numeric(), 6.28318, 1e-4) << "sin period is ~2pi";
        EXPECT_NEAR(p_tan->to_numeric(), 3.14159, 1e-4) << "tan period is ~pi";
    }
}

TEST(AssumptionDeclarations, ContradictionFiniteDivergent) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    auto failure_355 = store.declare_finiteness("x", Finiteness::Divergent);
    ASSERT_FALSE(failure_355.has_value()) << "Finite + Divergent returns failure";
    EXPECT_TRUE((failure_355.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionDivergentFinite) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("y", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";

    auto failure_373 = store.declare_finiteness("y", Finiteness::Finite);
    ASSERT_FALSE(failure_373.has_value()) << "Divergent + Finite returns failure";
    EXPECT_TRUE((failure_373.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionPositiveDefiniteNegativeDefinite) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("M", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    auto failure_391 = store.declare_definiteness("M", Definiteness::NegativeDefinite);
    ASSERT_FALSE(failure_391.has_value()) << "PositiveDefinite + NegativeDefinite returns failure";
    EXPECT_TRUE((failure_391.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionPositiveDefiniteIndefinite) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("A", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    auto failure_414 = store.declare_definiteness("A", Definiteness::Indefinite);
    EXPECT_TRUE((!failure_414.has_value())) << "PositiveDefinite + Indefinite returns failure";
}

TEST(AssumptionDeclarations, ContradictionNegativeDefiniteIndefinite) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("B", Definiteness::NegativeDefinite).has_value())) << "definiteness declaration succeeds";

    auto failure_429 = store.declare_definiteness("B", Definiteness::Indefinite);
    EXPECT_TRUE((!failure_429.has_value())) << "NegativeDefinite + Indefinite returns failure";
}

TEST(AssumptionDeclarations, ContradictionTranscendentalAlgebraic) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("pi").has_value())) << "transcendental declaration succeeds";

    auto failure_429 = store.declare_domain("pi", Domain::Algebraic);
    ASSERT_FALSE(failure_429.has_value()) << "Transcendental + Algebraic returns failure";
    EXPECT_TRUE((failure_429.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionTranscendentalRational) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("e").has_value())) << "transcendental declaration succeeds";

    auto failure_462 = store.declare_domain("e", Domain::Rational);
    EXPECT_TRUE((!failure_462.has_value())) << "Transcendental + Rational returns failure";
}

TEST(AssumptionDeclarations, ContradictionTranscendentalInteger) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("pi").has_value())) << "transcendental declaration succeeds";

    auto failure_477 = store.declare_domain("pi", Domain::Integer);
    EXPECT_TRUE((!failure_477.has_value())) << "Transcendental + Integer returns failure";
}

TEST(AssumptionDeclarations, ContradictionAlgebraicThenTranscendental) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("sqrt2", Domain::Algebraic).has_value())) << "domain declaration succeeds";

    auto failure_467 = store.declare_transcendental("sqrt2");
    ASSERT_FALSE(failure_467.has_value()) << "Algebraic + Transcendental returns failure";
    EXPECT_TRUE((failure_467.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionBoundedUnbounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    auto failure_485 = store.declare_bounded("x", Boundedness::Unbounded);
    ASSERT_FALSE(failure_485.has_value()) << "Bounded + Unbounded returns failure";
    EXPECT_TRUE((failure_485.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, ContradictionParityEvenOdd) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("n", Parity::Even).has_value())) << "parity declaration succeeds";

    auto failure_503 = store.declare_parity("n", Parity::Odd);
    ASSERT_FALSE(failure_503.has_value()) << "Even + Odd returns failure";
    EXPECT_TRUE((failure_503.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
}

TEST(AssumptionDeclarations, FiniteImpliesBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "Declaring Finite auto-sets Bounded";
}
