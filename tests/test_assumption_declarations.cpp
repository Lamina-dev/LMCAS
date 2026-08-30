
#include "test_common.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <stdexcept>

using namespace lamina;


static Interval make_closed_interval(double lo, double hi) {
    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lower_val);
    iv.upper = Endpoint::closed(upper_val);
    return iv;
}

static Interval make_open_interval(double lo, double hi) {
    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::open(lower_val);
    iv.upper = Endpoint::open(upper_val);
    return iv;
}


static void test_continuous_only_on_differentiable_overlap_throws() {
    TEST_CASE("Continuity overlap: declaring continuous-only on interval overlapping differentiable returns failure");

    PropertyStore store;

    // Declare differentiable on [0, 10]
    Interval diff_interval = make_closed_interval(0.0, 10.0);
    store.declare_differentiable("f", diff_interval);

    // Declaring continuous-only on [5, 15] (overlaps [0,10]) should throw
    Interval cont_interval = make_closed_interval(5.0, 15.0);
    auto failure_46 = store.declare_continuous("f", cont_interval);
    EXPECT_TRUE(!failure_46.has_value(), "Declaring continuous-only on differentiable overlap returns failure");
    EXPECT_TRUE(failure_46.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_46.error().message, {"Contradiction", "f", "differentiable"},
            "Result error message mentions contradiction, symbol, and differentiable");
}

static void test_continuous_only_on_exact_differentiable_interval_throws() {
    TEST_CASE("Continuity overlap: declaring continuous-only on exact differentiable interval returns failure");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 5.0);
    store.declare_differentiable("g", iv);

    auto failure_65 = store.declare_continuous("g", iv);
    EXPECT_TRUE(!failure_65.has_value(), "Declaring continuous-only on exact differentiable interval returns failure");
}

static void test_differentiable_on_continuous_overlap_ok() {
    TEST_CASE("Continuity overlap: declaring differentiable on continuous-only overlap is OK (upgrade)");

    PropertyStore store;

    // Declare continuous-only on [0, 10]
    Interval cont_interval = make_closed_interval(0.0, 10.0);
    store.declare_continuous("h", cont_interval);

    // Declaring differentiable on [5, 15] (overlaps) should NOT throw (it's an upgrade)
    Interval diff_interval = make_closed_interval(5.0, 15.0);
    auto success_76 = store.declare_differentiable("h", diff_interval);
    EXPECT_TRUE(success_76.has_value(), "Declaring differentiable on continuous-only overlap does not throw");
}

static void test_continuous_on_continuous_overlap_ok() {
    TEST_CASE("Continuity overlap: declaring continuous-only on continuous-only overlap is OK (idempotent)");

    PropertyStore store;

    Interval iv1 = make_closed_interval(0.0, 10.0);
    Interval iv2 = make_closed_interval(5.0, 15.0);
    store.declare_continuous("k", iv1);

    auto success_93 = store.declare_continuous("k", iv2);
    EXPECT_TRUE(success_93.has_value(), "Declaring continuous-only on continuous-only overlap does not throw");
}

static void test_differentiable_implies_continuous() {
    TEST_CASE("Differentiability implies continuity");

    PropertyStore store;

    Interval iv = make_closed_interval(1.0, 5.0);
    store.declare_differentiable("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv).value(),
        "Differentiable symbol is also continuous on same interval");
    EXPECT_TRUE(store.is_differentiable("f", iv).value(),
        "Differentiable symbol is differentiable on same interval");
}

static void test_continuous_not_differentiable() {
    TEST_CASE("Continuous-only is not differentiable");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 3.0);
    store.declare_continuous("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv).value(),
        "Continuous symbol is continuous");
    EXPECT_FALSE(store.is_differentiable("f", iv).value(),
        "Continuous-only symbol is NOT differentiable");
}

static void test_non_overlapping_intervals_no_conflict() {
    TEST_CASE("Non-overlapping intervals: differentiable and continuous-only on disjoint intervals OK");

    PropertyStore store;

    Interval diff_iv = make_closed_interval(0.0, 5.0);
    Interval cont_iv = make_closed_interval(6.0, 10.0);

    store.declare_differentiable("f", diff_iv);

    auto success_140 = store.declare_continuous("f", cont_iv);
    EXPECT_TRUE(success_140.has_value(), "Disjoint intervals do not conflict");
}


static void test_monotonicity_declare_and_retrieve_increasing() {
    TEST_CASE("Monotonicity: declare increasing and retrieve");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing);

    Monotonicity result = store.get_monotonicity("f", "x", iv).value();
    EXPECT_TRUE(result == Monotonicity::Increasing,
        "get_monotonicity returns Increasing for exact interval");
}

static void test_monotonicity_declare_and_retrieve_decreasing() {
    TEST_CASE("Monotonicity: declare decreasing and retrieve");

    PropertyStore store;

    Interval iv = make_closed_interval(-5.0, 5.0);
    store.declare_monotonicity("g", "t", iv, Monotonicity::Decreasing);

    Monotonicity result = store.get_monotonicity("g", "t", iv).value();
    EXPECT_TRUE(result == Monotonicity::Decreasing,
        "get_monotonicity returns Decreasing for exact interval");
}

static void test_monotonicity_sub_interval_coverage() {
    TEST_CASE("Monotonicity: sub-interval is covered by larger declaration");

    PropertyStore store;

    // Declare increasing on [0, 10]
    Interval large_iv = make_closed_interval(0.0, 10.0);
    store.declare_monotonicity("f", "x", large_iv, Monotonicity::Increasing);

    // Query on [2, 8] (sub-interval) should return Increasing
    Interval sub_iv = make_closed_interval(2.0, 8.0);
    Monotonicity result = store.get_monotonicity("f", "x", sub_iv).value();
    EXPECT_TRUE(result == Monotonicity::Increasing,
        "Sub-interval query returns Increasing (covered by larger declaration)");
}

static void test_monotonicity_uncovered_interval_returns_unknown() {
    TEST_CASE("Monotonicity: uncovered interval returns Unknown");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 5.0);
    store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing);

    // Query on [6, 10] (not covered) should return Unknown
    Interval uncovered = make_closed_interval(6.0, 10.0);
    Monotonicity result = store.get_monotonicity("f", "x", uncovered).value();
    EXPECT_TRUE(result == Monotonicity::Unknown,
        "Uncovered interval returns Unknown");
}

static void test_monotonicity_wrong_variable_returns_unknown() {
    TEST_CASE("Monotonicity: wrong variable returns Unknown");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing);

    // Query with different variable
    Monotonicity result = store.get_monotonicity("f", "y", iv).value();
    EXPECT_TRUE(result == Monotonicity::Unknown,
        "Query with wrong variable returns Unknown");
}

static void test_monotonicity_undeclared_symbol_returns_unknown() {
    TEST_CASE("Monotonicity: undeclared symbol returns Unknown");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity result = store.get_monotonicity("undeclared", "x", iv).value();
    EXPECT_TRUE(result == Monotonicity::Unknown,
        "Undeclared symbol returns Unknown");
}

static void test_monotonicity_multiple_declarations() {
    TEST_CASE("Monotonicity: multiple declarations for different intervals");

    PropertyStore store;

    Interval iv1 = make_closed_interval(0.0, 5.0);
    Interval iv2 = make_closed_interval(5.0, 10.0);

    store.declare_monotonicity("f", "x", iv1, Monotonicity::Increasing);
    store.declare_monotonicity("f", "x", iv2, Monotonicity::Decreasing);

    EXPECT_TRUE(store.get_monotonicity("f", "x", iv1).value() == Monotonicity::Increasing,
        "First interval returns Increasing");
    EXPECT_TRUE(store.get_monotonicity("f", "x", iv2).value() == Monotonicity::Decreasing,
        "Second interval returns Decreasing");
}

static void test_monotonicity_entire_line() {
    TEST_CASE("Monotonicity: declaration on entire line covers any sub-interval");

    PropertyStore store;

    Interval entire = Interval::entire_line();
    store.declare_monotonicity("exp", "x", entire, Monotonicity::Increasing);

    Interval sub = make_closed_interval(-100.0, 100.0);
    EXPECT_TRUE(store.get_monotonicity("exp", "x", sub).value() == Monotonicity::Increasing,
        "Entire line declaration covers any finite sub-interval");
}


static void test_periodicity_declare_and_retrieve() {
    TEST_CASE("Periodicity: declare periodic and retrieve period");

    PropertyStore store;

    // Create a period expression: 2*pi (represented as a constant for simplicity)
    auto period = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(6.283185307)));

    store.declare_periodic("sin_x", period);

    EXPECT_TRUE(store.is_periodic("sin_x"), "Symbol is periodic after declaration");

    auto retrieved = store.get_period("sin_x");
    EXPECT_TRUE(retrieved.has_value(), "get_period returns a value");
}

static void test_periodicity_not_periodic_by_default() {
    TEST_CASE("Periodicity: undeclared symbol is not periodic");

    PropertyStore store;

    EXPECT_FALSE(store.is_periodic("undeclared"), "Undeclared symbol is not periodic");

    auto period = store.get_period("undeclared");
    EXPECT_FALSE(period.has_value(), "get_period returns nullopt for undeclared symbol");
}

static void test_periodicity_overwrite_period() {
    TEST_CASE("Periodicity: re-declaring period overwrites previous");

    PropertyStore store;

    auto period1 = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(3.14159)));
    auto period2 = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(6.28318)));

    store.declare_periodic("f", period1);
    EXPECT_TRUE(store.is_periodic("f"), "f is periodic after first declaration");

    store.declare_periodic("f", period2);
    EXPECT_TRUE(store.is_periodic("f"), "f is still periodic after second declaration");

    auto retrieved = store.get_period("f");
    EXPECT_TRUE(retrieved.has_value(), "get_period returns a value after overwrite");
    if (retrieved.has_value()) {
        double val = (*retrieved)->to_numeric();
        EXPECT_NEAR(val, 6.28318, 1e-4, "Period is updated to new value");
    }
}

static void test_periodicity_different_symbols() {
    TEST_CASE("Periodicity: different symbols have independent periods");

    PropertyStore store;

    auto period_sin = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(6.28318)));
    auto period_tan = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(3.14159)));

    store.declare_periodic("sin_x", period_sin);
    store.declare_periodic("tan_x", period_tan);

    EXPECT_TRUE(store.is_periodic("sin_x"), "sin_x is periodic");
    EXPECT_TRUE(store.is_periodic("tan_x"), "tan_x is periodic");

    auto p_sin = store.get_period("sin_x");
    auto p_tan = store.get_period("tan_x");

    EXPECT_TRUE(p_sin.has_value(), "sin_x has a period");
    EXPECT_TRUE(p_tan.has_value(), "tan_x has a period");

    if (p_sin.has_value() && p_tan.has_value()) {
        EXPECT_NEAR((*p_sin)->to_numeric(), 6.28318, 1e-4, "sin period is ~2pi");
        EXPECT_NEAR((*p_tan)->to_numeric(), 3.14159, 1e-4, "tan period is ~pi");
    }
}


static void test_contradiction_finite_divergent() {
    TEST_CASE("Contradiction: Finite then Divergent returns failure");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    auto failure_355 = store.declare_finiteness("x", Finiteness::Divergent);
    EXPECT_TRUE(!failure_355.has_value(), "Finite + Divergent returns failure");
    EXPECT_TRUE(failure_355.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_355.error().message, {"x", "Finite", "Divergent"},
            "Result error mentions symbol, Finite, and Divergent");
}

static void test_contradiction_divergent_finite() {
    TEST_CASE("Contradiction: Divergent then Finite returns failure");

    PropertyStore store;
    store.declare_finiteness("y", Finiteness::Divergent);

    auto failure_373 = store.declare_finiteness("y", Finiteness::Finite);
    EXPECT_TRUE(!failure_373.has_value(), "Divergent + Finite returns failure");
    EXPECT_TRUE(failure_373.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_373.error().message, {"y", "Finite", "Divergent"},
            "Result error mentions symbol, Finite, and Divergent");
}

static void test_contradiction_positive_definite_negative_definite() {
    TEST_CASE("Contradiction: PositiveDefinite then NegativeDefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("M", Definiteness::PositiveDefinite);

    auto failure_391 = store.declare_definiteness("M", Definiteness::NegativeDefinite);
    EXPECT_TRUE(!failure_391.has_value(), "PositiveDefinite + NegativeDefinite returns failure");
    EXPECT_TRUE(failure_391.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_391.error().message, {"M", "definiteness"},
            "Result error mentions symbol and definiteness");
}

static void test_contradiction_positive_definite_indefinite() {
    TEST_CASE("Contradiction: PositiveDefinite then Indefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("A", Definiteness::PositiveDefinite);

    auto failure_414 = store.declare_definiteness("A", Definiteness::Indefinite);
    EXPECT_TRUE(!failure_414.has_value(), "PositiveDefinite + Indefinite returns failure");
}

static void test_contradiction_negative_definite_indefinite() {
    TEST_CASE("Contradiction: NegativeDefinite then Indefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("B", Definiteness::NegativeDefinite);

    auto failure_429 = store.declare_definiteness("B", Definiteness::Indefinite);
    EXPECT_TRUE(!failure_429.has_value(), "NegativeDefinite + Indefinite returns failure");
}

static void test_contradiction_transcendental_algebraic() {
    TEST_CASE("Contradiction: Transcendental then Algebraic domain returns failure");

    PropertyStore store;
    store.declare_transcendental("pi");

    auto failure_429 = store.declare_domain("pi", Domain::Algebraic);
    EXPECT_TRUE(!failure_429.has_value(), "Transcendental + Algebraic returns failure");
    EXPECT_TRUE(failure_429.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_429.error().message, {"pi", "Transcendental", "Algebraic"},
            "Result error mentions symbol, Transcendental, and Algebraic");
}

static void test_contradiction_transcendental_rational() {
    TEST_CASE("Contradiction: Transcendental then Rational domain returns failure");

    PropertyStore store;
    store.declare_transcendental("e");

    auto failure_462 = store.declare_domain("e", Domain::Rational);
    EXPECT_TRUE(!failure_462.has_value(), "Transcendental + Rational returns failure");
}

static void test_contradiction_transcendental_integer() {
    TEST_CASE("Contradiction: Transcendental then Integer domain returns failure");

    PropertyStore store;
    store.declare_transcendental("pi");

    auto failure_477 = store.declare_domain("pi", Domain::Integer);
    EXPECT_TRUE(!failure_477.has_value(), "Transcendental + Integer returns failure");
}

static void test_contradiction_algebraic_then_transcendental() {
    TEST_CASE("Contradiction: Algebraic domain then Transcendental returns failure");

    PropertyStore store;
    store.declare_domain("sqrt2", Domain::Algebraic);

    auto failure_467 = store.declare_transcendental("sqrt2");
    EXPECT_TRUE(!failure_467.has_value(), "Algebraic + Transcendental returns failure");
    EXPECT_TRUE(failure_467.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_467.error().message, {"sqrt2", "Transcendental", "Algebraic"},
            "Result error mentions symbol, Transcendental, and Algebraic");
}

static void test_contradiction_bounded_unbounded() {
    TEST_CASE("Contradiction: Bounded then Unbounded returns failure");

    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    auto failure_485 = store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(!failure_485.has_value(), "Bounded + Unbounded returns failure");
    EXPECT_TRUE(failure_485.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_485.error().message, {"x", "boundedness"},
            "Result error mentions symbol and boundedness");
}

static void test_contradiction_parity_even_odd() {
    TEST_CASE("Contradiction: Even then Odd parity returns failure");

    PropertyStore store;
    store.declare_parity("n", Parity::Even);

    auto failure_503 = store.declare_parity("n", Parity::Odd);
    EXPECT_TRUE(!failure_503.has_value(), "Even + Odd returns failure");
    EXPECT_TRUE(failure_503.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    EXPECT_CONTAINS(failure_503.error().message, {"n", "parity"},
            "Result error mentions symbol and parity");
}


static void test_finite_implies_bounded() {
    TEST_CASE("Finiteness: Finite implies Bounded");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "Declaring Finite auto-sets Bounded");
}

static void test_checked_interval_property_contracts() {
    TEST_CASE("Checked interval properties propagate validation, budgets, and conflicts");

    PropertyStore store;
    Interval interval = make_closed_interval(0.0, 10.0);

    auto empty_symbol = store.declare_continuous_checked("", interval);
    EXPECT_TRUE(!empty_symbol && empty_symbol.error().code == CasErrc::InvalidArgument,
                "checked continuity rejects an empty symbol");
    EXPECT_TRUE(store.get_all_symbols().empty(),
                "failed checked continuity creates no property record");

    auto empty_interval = store.declare_continuous_checked("f", Interval::empty());
    EXPECT_TRUE(!empty_interval && empty_interval.error().code == CasErrc::InvalidArgument,
                "checked continuity rejects an empty interval");
    EXPECT_TRUE(store.get_all_symbols().empty(),
                "empty-interval failure is transactional");

    auto differentiable = store.declare_differentiable_checked("f", interval);
    EXPECT_TRUE(differentiable.has_value(),
                "checked differentiability accepts a valid interval");
    const auto declaration_count = store.get_continuity_decls("f").size();
    auto duplicate = store.declare_differentiable_checked("f", interval);
    EXPECT_TRUE(duplicate.has_value(), "exact duplicate differentiability is idempotent");
    EXPECT_TRUE(store.get_continuity_decls("f").size() == declaration_count,
                "idempotent checked declaration does not duplicate state");

    auto downgrade = store.declare_continuous_checked("f", interval);
    EXPECT_TRUE(!downgrade && downgrade.error().code == CasErrc::InvalidArgument,
                "checked continuity rejects a differentiability downgrade");
    EXPECT_TRUE(store.get_continuity_decls("f").size() == declaration_count,
                "failed checked downgrade preserves declaration state");

    PropertyStore boundary_store;
    Interval left = make_closed_interval(0.0, 5.0);
    Interval open_touch{
        Endpoint::open(SymbolicExpr::number(5.0)),
        Endpoint::closed(SymbolicExpr::number(10.0))};
    Interval closed_touch = make_closed_interval(5.0, 10.0);
    EXPECT_TRUE(boundary_store.declare_differentiable_checked("edge", left).has_value(),
                "boundary test stores differentiability");
    auto disjoint_touch = boundary_store.declare_continuous_checked("edge", open_touch);
    EXPECT_TRUE(disjoint_touch.has_value(),
                "open touching boundary is correctly treated as non-overlapping");
    auto overlapping_touch = boundary_store.declare_continuous_checked("edge", closed_touch);
    EXPECT_TRUE(!overlapping_touch &&
                    overlapping_touch.error().code == CasErrc::InvalidArgument,
                "closed touching boundary is correctly treated as overlapping");

    Interval symbolic{
        Endpoint::closed(SymbolicExpr::variable("a")),
        Endpoint::closed(SymbolicExpr::variable("b"))};
    auto symbolic_declaration = store.declare_continuous_checked("g", symbolic);
    EXPECT_TRUE(!symbolic_declaration &&
                    symbolic_declaration.error().code == CasErrc::UnboundSymbol,
                "checked continuity propagates unbound symbolic endpoints");
    auto symbolic_query = store.is_continuous_checked("f", symbolic);
    EXPECT_TRUE(!symbolic_query && symbolic_query.error().code == CasErrc::UnboundSymbol,
                "checked continuity query propagates endpoint errors");

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = store.declare_monotonicity_checked(
        "g", "x", interval, Monotonicity::Increasing, cancelled_context);
    EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                "checked monotonicity observes cancellation");
    EXPECT_TRUE(store.get_monotonicity_decls("g").empty(),
                "cancelled monotonicity declaration stores no state");

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited = store.declare_differentiable_checked(
        "h", interval, limited_context);
    EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                "checked differentiability observes the step budget");

    auto empty_variable = store.declare_monotonicity_checked(
        "f", "", interval, Monotonicity::Increasing);
    EXPECT_TRUE(!empty_variable && empty_variable.error().code == CasErrc::InvalidArgument,
                "checked monotonicity rejects an empty variable");
    auto monotonic = store.declare_monotonicity_checked(
        "f", "x", interval, Monotonicity::Increasing);
    EXPECT_TRUE(monotonic.has_value(), "checked monotonicity stores a valid declaration");
    auto queried = store.get_monotonicity_checked("f", "x", interval);
    EXPECT_TRUE(queried && queried.value() == Monotonicity::Increasing,
                "checked monotonicity query retrieves the proven declaration");
}


int main() {
    // Continuity/differentiability overlap detection
    test_continuous_only_on_differentiable_overlap_throws();
    test_continuous_only_on_exact_differentiable_interval_throws();
    test_differentiable_on_continuous_overlap_ok();
    test_continuous_on_continuous_overlap_ok();
    test_differentiable_implies_continuous();
    test_continuous_not_differentiable();
    test_non_overlapping_intervals_no_conflict();

    // Monotonicity storage and retrieval
    test_monotonicity_declare_and_retrieve_increasing();
    test_monotonicity_declare_and_retrieve_decreasing();
    test_monotonicity_sub_interval_coverage();
    test_monotonicity_uncovered_interval_returns_unknown();
    test_monotonicity_wrong_variable_returns_unknown();
    test_monotonicity_undeclared_symbol_returns_unknown();
    test_monotonicity_multiple_declarations();
    test_monotonicity_entire_line();

    // Periodicity storage and retrieval
    test_periodicity_declare_and_retrieve();
    test_periodicity_not_periodic_by_default();
    test_periodicity_overwrite_period();
    test_periodicity_different_symbols();

    // Contradiction scenarios
    test_contradiction_finite_divergent();
    test_contradiction_divergent_finite();
    test_contradiction_positive_definite_negative_definite();
    test_contradiction_positive_definite_indefinite();
    test_contradiction_negative_definite_indefinite();
    test_contradiction_transcendental_algebraic();
    test_contradiction_transcendental_rational();
    test_contradiction_transcendental_integer();
    test_contradiction_algebraic_then_transcendental();
    test_contradiction_bounded_unbounded();
    test_contradiction_parity_even_odd();

    // Finiteness implication
    test_finite_implies_bounded();
    test_checked_interval_property_contracts();

    return TEST_REPORT();
}
