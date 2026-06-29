/**
 * @file test_assumption_declarations.cpp
 * @brief Unit tests for PropertyStore declaration methods (Task 2.5).
 *
 * Tests:
 * - Continuity/differentiability overlap detection
 * - Monotonicity storage and retrieval
 * - Periodicity storage and retrieval
 * - All contradiction scenarios produce correct exception messages
 *
 * Validates: Requirements 6.4, 7.1, 7.2, 11.2
 */

#include "test_common.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <stdexcept>

using namespace lamina;

// ============================================================
// Helpers: create intervals for testing
// ============================================================

static Interval make_closed_interval(double lo, double hi) {
    auto lower_val = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lower_val);
    iv.upper = Endpoint::closed(upper_val);
    return iv;
}

static Interval make_open_interval(double lo, double hi) {
    auto lower_val = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::open(lower_val);
    iv.upper = Endpoint::open(upper_val);
    return iv;
}

// ============================================================
// Tests: Continuity/Differentiability overlap detection (Req 6.4)
// ============================================================

static void test_continuous_only_on_differentiable_overlap_throws() {
    TEST_CASE("Continuity overlap: declaring continuous-only on interval overlapping differentiable throws");

    PropertyStore store;

    // Declare differentiable on [0, 10]
    Interval diff_interval = make_closed_interval(0.0, 10.0);
    store.declare_differentiable("f", diff_interval);

    // Declaring continuous-only on [5, 15] (overlaps [0,10]) should throw
    Interval cont_interval = make_closed_interval(5.0, 15.0);
    bool threw = false;
    try {
        store.declare_continuous("f", cont_interval);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"Contradiction", "f", "differentiable"},
            "Exception message mentions contradiction, symbol, and differentiable");
    }
    EXPECT_TRUE(threw, "Declaring continuous-only on differentiable overlap throws");
}

static void test_continuous_only_on_exact_differentiable_interval_throws() {
    TEST_CASE("Continuity overlap: declaring continuous-only on exact differentiable interval throws");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 5.0);
    store.declare_differentiable("g", iv);

    bool threw = false;
    try {
        store.declare_continuous("g", iv);
    } catch (const std::invalid_argument& e) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring continuous-only on exact differentiable interval throws");
}

static void test_differentiable_on_continuous_overlap_ok() {
    TEST_CASE("Continuity overlap: declaring differentiable on continuous-only overlap is OK (upgrade)");

    PropertyStore store;

    // Declare continuous-only on [0, 10]
    Interval cont_interval = make_closed_interval(0.0, 10.0);
    store.declare_continuous("h", cont_interval);

    // Declaring differentiable on [5, 15] (overlaps) should NOT throw (it's an upgrade)
    Interval diff_interval = make_closed_interval(5.0, 15.0);
    bool threw = false;
    try {
        store.declare_differentiable("h", diff_interval);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_FALSE(threw, "Declaring differentiable on continuous-only overlap does not throw");
}

static void test_continuous_on_continuous_overlap_ok() {
    TEST_CASE("Continuity overlap: declaring continuous-only on continuous-only overlap is OK (idempotent)");

    PropertyStore store;

    Interval iv1 = make_closed_interval(0.0, 10.0);
    Interval iv2 = make_closed_interval(5.0, 15.0);
    store.declare_continuous("k", iv1);

    bool threw = false;
    try {
        store.declare_continuous("k", iv2);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_FALSE(threw, "Declaring continuous-only on continuous-only overlap does not throw");
}

static void test_differentiable_implies_continuous() {
    TEST_CASE("Differentiability implies continuity");

    PropertyStore store;

    Interval iv = make_closed_interval(1.0, 5.0);
    store.declare_differentiable("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv),
        "Differentiable symbol is also continuous on same interval");
    EXPECT_TRUE(store.is_differentiable("f", iv),
        "Differentiable symbol is differentiable on same interval");
}

static void test_continuous_not_differentiable() {
    TEST_CASE("Continuous-only is not differentiable");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 3.0);
    store.declare_continuous("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv),
        "Continuous symbol is continuous");
    EXPECT_FALSE(store.is_differentiable("f", iv),
        "Continuous-only symbol is NOT differentiable");
}

static void test_non_overlapping_intervals_no_conflict() {
    TEST_CASE("Non-overlapping intervals: differentiable and continuous-only on disjoint intervals OK");

    PropertyStore store;

    Interval diff_iv = make_closed_interval(0.0, 5.0);
    Interval cont_iv = make_closed_interval(6.0, 10.0);

    store.declare_differentiable("f", diff_iv);

    bool threw = false;
    try {
        store.declare_continuous("f", cont_iv);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_FALSE(threw, "Disjoint intervals do not conflict");
}

// ============================================================
// Tests: Monotonicity storage and retrieval (Req 7.1, 7.2)
// ============================================================

static void test_monotonicity_declare_and_retrieve_increasing() {
    TEST_CASE("Monotonicity: declare increasing and retrieve");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing);

    Monotonicity result = store.get_monotonicity("f", "x", iv);
    EXPECT_TRUE(result == Monotonicity::Increasing,
        "get_monotonicity returns Increasing for exact interval");
}

static void test_monotonicity_declare_and_retrieve_decreasing() {
    TEST_CASE("Monotonicity: declare decreasing and retrieve");

    PropertyStore store;

    Interval iv = make_closed_interval(-5.0, 5.0);
    store.declare_monotonicity("g", "t", iv, Monotonicity::Decreasing);

    Monotonicity result = store.get_monotonicity("g", "t", iv);
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
    Monotonicity result = store.get_monotonicity("f", "x", sub_iv);
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
    Monotonicity result = store.get_monotonicity("f", "x", uncovered);
    EXPECT_TRUE(result == Monotonicity::Unknown,
        "Uncovered interval returns Unknown");
}

static void test_monotonicity_wrong_variable_returns_unknown() {
    TEST_CASE("Monotonicity: wrong variable returns Unknown");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    store.declare_monotonicity("f", "x", iv, Monotonicity::Increasing);

    // Query with different variable
    Monotonicity result = store.get_monotonicity("f", "y", iv);
    EXPECT_TRUE(result == Monotonicity::Unknown,
        "Query with wrong variable returns Unknown");
}

static void test_monotonicity_undeclared_symbol_returns_unknown() {
    TEST_CASE("Monotonicity: undeclared symbol returns Unknown");

    PropertyStore store;

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity result = store.get_monotonicity("undeclared", "x", iv);
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

    EXPECT_TRUE(store.get_monotonicity("f", "x", iv1) == Monotonicity::Increasing,
        "First interval returns Increasing");
    EXPECT_TRUE(store.get_monotonicity("f", "x", iv2) == Monotonicity::Decreasing,
        "Second interval returns Decreasing");
}

static void test_monotonicity_entire_line() {
    TEST_CASE("Monotonicity: declaration on entire line covers any sub-interval");

    PropertyStore store;

    Interval entire = Interval::entire_line();
    store.declare_monotonicity("exp", "x", entire, Monotonicity::Increasing);

    Interval sub = make_closed_interval(-100.0, 100.0);
    EXPECT_TRUE(store.get_monotonicity("exp", "x", sub) == Monotonicity::Increasing,
        "Entire line declaration covers any finite sub-interval");
}

// ============================================================
// Tests: Periodicity storage and retrieval (Req 11.2)
// ============================================================

static void test_periodicity_declare_and_retrieve() {
    TEST_CASE("Periodicity: declare periodic and retrieve period");

    PropertyStore store;

    // Create a period expression: 2*pi (represented as a constant for simplicity)
    auto period = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(6.283185307)));

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

    auto period1 = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(3.14159)));
    auto period2 = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(6.28318)));

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

    auto period_sin = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(6.28318)));
    auto period_tan = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(3.14159)));

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

// ============================================================
// Tests: Contradiction scenarios (various requirements)
// ============================================================

static void test_contradiction_finite_divergent() {
    TEST_CASE("Contradiction: Finite then Divergent throws");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    bool threw = false;
    try {
        store.declare_finiteness("x", Finiteness::Divergent);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"x", "Finite", "Divergent"},
            "Exception mentions symbol, Finite, and Divergent");
    }
    EXPECT_TRUE(threw, "Finite + Divergent throws");
}

static void test_contradiction_divergent_finite() {
    TEST_CASE("Contradiction: Divergent then Finite throws");

    PropertyStore store;
    store.declare_finiteness("y", Finiteness::Divergent);

    bool threw = false;
    try {
        store.declare_finiteness("y", Finiteness::Finite);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"y", "Finite", "Divergent"},
            "Exception mentions symbol, Finite, and Divergent");
    }
    EXPECT_TRUE(threw, "Divergent + Finite throws");
}

static void test_contradiction_positive_definite_negative_definite() {
    TEST_CASE("Contradiction: PositiveDefinite then NegativeDefinite throws");

    PropertyStore store;
    store.declare_definiteness("M", Definiteness::PositiveDefinite);

    bool threw = false;
    try {
        store.declare_definiteness("M", Definiteness::NegativeDefinite);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"M", "definiteness"},
            "Exception mentions symbol and definiteness");
    }
    EXPECT_TRUE(threw, "PositiveDefinite + NegativeDefinite throws");
}

static void test_contradiction_positive_definite_indefinite() {
    TEST_CASE("Contradiction: PositiveDefinite then Indefinite throws");

    PropertyStore store;
    store.declare_definiteness("A", Definiteness::PositiveDefinite);

    bool threw = false;
    try {
        store.declare_definiteness("A", Definiteness::Indefinite);
    } catch (const std::invalid_argument& e) {
        threw = true;
    }
    EXPECT_TRUE(threw, "PositiveDefinite + Indefinite throws");
}

static void test_contradiction_negative_definite_indefinite() {
    TEST_CASE("Contradiction: NegativeDefinite then Indefinite throws");

    PropertyStore store;
    store.declare_definiteness("B", Definiteness::NegativeDefinite);

    bool threw = false;
    try {
        store.declare_definiteness("B", Definiteness::Indefinite);
    } catch (const std::invalid_argument& e) {
        threw = true;
    }
    EXPECT_TRUE(threw, "NegativeDefinite + Indefinite throws");
}

static void test_contradiction_transcendental_algebraic() {
    TEST_CASE("Contradiction: Transcendental then Algebraic domain throws");

    PropertyStore store;
    store.declare_transcendental("pi");

    bool threw = false;
    try {
        store.declare_domain("pi", Domain::Algebraic);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"pi", "Transcendental", "Algebraic"},
            "Exception mentions symbol, Transcendental, and Algebraic");
    }
    EXPECT_TRUE(threw, "Transcendental + Algebraic throws");
}

static void test_contradiction_transcendental_rational() {
    TEST_CASE("Contradiction: Transcendental then Rational domain throws");

    PropertyStore store;
    store.declare_transcendental("e");

    bool threw = false;
    try {
        store.declare_domain("e", Domain::Rational);
    } catch (const std::invalid_argument& e) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Transcendental + Rational throws");
}

static void test_contradiction_transcendental_integer() {
    TEST_CASE("Contradiction: Transcendental then Integer domain throws");

    PropertyStore store;
    store.declare_transcendental("pi");

    bool threw = false;
    try {
        store.declare_domain("pi", Domain::Integer);
    } catch (const std::invalid_argument& e) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Transcendental + Integer throws");
}

static void test_contradiction_algebraic_then_transcendental() {
    TEST_CASE("Contradiction: Algebraic domain then Transcendental throws");

    PropertyStore store;
    store.declare_domain("sqrt2", Domain::Algebraic);

    bool threw = false;
    try {
        store.declare_transcendental("sqrt2");
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"sqrt2", "Transcendental", "Algebraic"},
            "Exception mentions symbol, Transcendental, and Algebraic");
    }
    EXPECT_TRUE(threw, "Algebraic + Transcendental throws");
}

static void test_contradiction_bounded_unbounded() {
    TEST_CASE("Contradiction: Bounded then Unbounded throws");

    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    bool threw = false;
    try {
        store.declare_bounded("x", Boundedness::Unbounded);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"x", "boundedness"},
            "Exception mentions symbol and boundedness");
    }
    EXPECT_TRUE(threw, "Bounded + Unbounded throws");
}

static void test_contradiction_parity_even_odd() {
    TEST_CASE("Contradiction: Even then Odd parity throws");

    PropertyStore store;
    store.declare_parity("n", Parity::Even);

    bool threw = false;
    try {
        store.declare_parity("n", Parity::Odd);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::string msg = e.what();
        EXPECT_CONTAINS(msg, {"n", "parity"},
            "Exception mentions symbol and parity");
    }
    EXPECT_TRUE(threw, "Even + Odd throws");
}

// ============================================================
// Tests: Finiteness implies Bounded (Req 9.3)
// ============================================================

static void test_finite_implies_bounded() {
    TEST_CASE("Finiteness: Finite implies Bounded");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "Declaring Finite auto-sets Bounded");
}

// ============================================================
// main
// ============================================================

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

    return TEST_REPORT();
}
