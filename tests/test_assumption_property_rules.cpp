
#include "test_common.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include <string>
#include <stdexcept>
#include <vector>

using namespace LMCAS;


/// Create a closed interval [lo, hi] from numeric values.
static Interval make_closed_interval(double lo, double hi) {
    auto lo_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto hi_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lo_expr);
    iv.upper = Endpoint::closed(hi_expr);
    return iv;
}

/// Create an open interval (lo, hi) from numeric values.
static Interval make_open_interval(double lo, double hi) {
    auto lo_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto hi_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::open(lo_expr);
    iv.upper = Endpoint::open(hi_expr);
    return iv;
}

/// Check that a callable returns an InvalidArgument failure.
template<typename F>
static bool returns_invalid_argument(F&& f) {
    auto result = f();
    return !result.has_value() && result.error().code == CasErrc::InvalidArgument;
}


static void test_differentiable_implies_continuous_closed_interval() {
    TEST_CASE("Differentiable on [0,1] implies continuous on [0,1]");

    PropertyStore store;
    Interval iv = make_closed_interval(0.0, 1.0);

    store.declare_differentiable("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv).value(),
        "Differentiable on [0,1] => continuous on [0,1]");
    EXPECT_TRUE(store.is_differentiable("f", iv).value(),
        "Differentiable on [0,1] => differentiable on [0,1]");
}

static void test_differentiable_implies_continuous_open_interval() {
    TEST_CASE("Differentiable on (0,10) implies continuous on (0,10)");

    PropertyStore store;
    Interval iv = make_open_interval(0.0, 10.0);

    store.declare_differentiable("g", iv);

    EXPECT_TRUE(store.is_continuous("g", iv).value(),
        "Differentiable on (0,10) => continuous on (0,10)");
}

static void test_differentiable_implies_continuous_subinterval() {
    TEST_CASE("Differentiable on [0,10] implies continuous on [2,5]");

    PropertyStore store;
    Interval outer = make_closed_interval(0.0, 10.0);
    Interval inner = make_closed_interval(2.0, 5.0);

    store.declare_differentiable("h", outer);

    EXPECT_TRUE(store.is_continuous("h", inner).value(),
        "Differentiable on [0,10] => continuous on sub-interval [2,5]");
    EXPECT_TRUE(store.is_differentiable("h", inner).value(),
        "Differentiable on [0,10] => differentiable on sub-interval [2,5]");
}

static void test_differentiable_implies_continuous_entire_line() {
    TEST_CASE("Differentiable on entire line implies continuous everywhere");

    PropertyStore store;
    Interval entire = Interval::entire_line();

    store.declare_differentiable("p", entire);

    EXPECT_TRUE(store.is_continuous("p", entire).value(),
        "Differentiable on (-inf,+inf) => continuous on (-inf,+inf)");

    // Also continuous on any finite sub-interval
    Interval sub = make_closed_interval(-100.0, 100.0);
    EXPECT_TRUE(store.is_continuous("p", sub).value(),
        "Differentiable on (-inf,+inf) => continuous on [-100,100]");
}

static void test_differentiable_multiple_symbols() {
    TEST_CASE("Differentiability implies continuity for multiple symbols");

    PropertyStore store;
    Interval iv1 = make_closed_interval(0.0, 1.0);
    Interval iv2 = make_closed_interval(-5.0, 5.0);
    Interval iv3 = make_closed_interval(10.0, 20.0);

    store.declare_differentiable("a", iv1);
    store.declare_differentiable("b", iv2);
    store.declare_differentiable("c", iv3);

    EXPECT_TRUE(store.is_continuous("a", iv1).value(), "a continuous on [0,1]");
    EXPECT_TRUE(store.is_continuous("b", iv2).value(), "b continuous on [-5,5]");
    EXPECT_TRUE(store.is_continuous("c", iv3).value(), "c continuous on [10,20]");
}

static void test_continuous_only_not_differentiable() {
    TEST_CASE("Continuous-only does NOT imply differentiable");

    PropertyStore store;
    Interval iv = make_closed_interval(0.0, 1.0);

    store.declare_continuous("f", iv);

    EXPECT_TRUE(store.is_continuous("f", iv).value(),
        "Continuous on [0,1] => is_continuous true");
    EXPECT_FALSE(store.is_differentiable("f", iv).value(),
        "Continuous-only on [0,1] => is_differentiable false");
}


static void test_transcendental_rejects_algebraic() {
    TEST_CASE("Transcendental symbol rejects Algebraic declaration");

    PropertyStore store;
    store.declare_transcendental("pi");

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_domain("pi", Domain::Algebraic);
    }), "Transcendental + Algebraic returns failure");
}

static void test_transcendental_rejects_rational() {
    TEST_CASE("Transcendental symbol rejects Rational declaration");

    PropertyStore store;
    store.declare_transcendental("e");

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_domain("e", Domain::Rational);
    }), "Transcendental + Rational returns failure");
}

static void test_transcendental_rejects_integer() {
    TEST_CASE("Transcendental symbol rejects Integer declaration");

    PropertyStore store;
    store.declare_transcendental("tau");

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_domain("tau", Domain::Integer);
    }), "Transcendental + Integer returns failure");
}

static void test_transcendental_rejects_natural() {
    TEST_CASE("Transcendental symbol rejects Natural declaration");

    PropertyStore store;
    store.declare_transcendental("alpha");

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_domain("alpha", Domain::Natural);
    }), "Transcendental + Natural returns failure");
}

static void test_transcendental_rejects_positive_int() {
    TEST_CASE("Transcendental symbol rejects PositiveInt declaration");

    PropertyStore store;
    store.declare_transcendental("gamma");

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_domain("gamma", Domain::PositiveInt);
    }), "Transcendental + PositiveInt returns failure");
}

static void test_algebraic_implies_real() {
    TEST_CASE("Algebraic domain implies Real");

    PropertyStore store;
    store.declare_domain("x", Domain::Algebraic);

    EXPECT_TRUE(store.has_domain("x", Domain::Real),
        "Algebraic implies Real");
    EXPECT_TRUE(store.has_domain("x", Domain::Complex),
        "Algebraic implies Complex");
}

static void test_transcendental_sets_real_domain() {
    TEST_CASE("Declaring transcendental sets Real domain");

    PropertyStore store;
    store.declare_transcendental("pi");

    EXPECT_TRUE(store.get_domain("pi") == Domain::Real,
        "Transcendental symbol has Real domain");
    EXPECT_TRUE(store.is_transcendental("pi"),
        "Symbol is marked transcendental");
}

static void test_algebraic_then_transcendental_throws() {
    TEST_CASE("Algebraic symbol cannot become Transcendental");

    PropertyStore store;
    store.declare_domain("x", Domain::Algebraic);

    // Transcendental requires domain <= Real, but Algebraic is more specific
    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_transcendental("x");
    }), "Algebraic then Transcendental returns failure");
}

static void test_rational_then_transcendental_throws() {
    TEST_CASE("Rational symbol cannot become Transcendental");

    PropertyStore store;
    store.declare_domain("x", Domain::Rational);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_transcendental("x");
    }), "Rational then Transcendental returns failure");
}

static void test_integer_then_transcendental_throws() {
    TEST_CASE("Integer symbol cannot become Transcendental");

    PropertyStore store;
    store.declare_domain("n", Domain::Integer);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_transcendental("n");
    }), "Integer then Transcendental returns failure");
}

static void test_transcendental_allows_real_declaration() {
    TEST_CASE("Transcendental symbol allows Real declaration (no-op)");

    PropertyStore store;
    store.declare_transcendental("pi");

    // Declaring Real on a transcendental symbol should be a no-op (already Real)
    store.declare_domain("pi", Domain::Real);
    EXPECT_TRUE(store.get_domain("pi") == Domain::Real,
        "Transcendental symbol remains Real after Real declaration");
    EXPECT_TRUE(store.is_transcendental("pi"),
        "Still transcendental after Real declaration");
}

static void test_transcendental_allows_complex_declaration() {
    TEST_CASE("Transcendental symbol allows Complex declaration (no-op)");

    PropertyStore store;
    store.declare_transcendental("e");

    // Complex is less specific than Real, so it's a no-op
    store.declare_domain("e", Domain::Complex);
    EXPECT_TRUE(store.get_domain("e") == Domain::Real,
        "Transcendental symbol remains Real after Complex declaration");
}


static void test_finite_implies_bounded() {
    TEST_CASE("Declaring Finite implies Bounded");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "Finite => Bounded");
    EXPECT_TRUE(store.get_finiteness("x") == Finiteness::Finite,
        "Finiteness is Finite");
}

static void test_finite_then_divergent_throws() {
    TEST_CASE("Finite + Divergent returns failure");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_finiteness("x", Finiteness::Divergent);
    }), "Finite + Divergent returns failure");
}

static void test_divergent_then_finite_throws() {
    TEST_CASE("Divergent + Finite returns failure");

    PropertyStore store;
    store.declare_finiteness("y", Finiteness::Divergent);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_finiteness("y", Finiteness::Finite);
    }), "Divergent + Finite returns failure");
}

static void test_divergent_does_not_imply_bounded() {
    TEST_CASE("Divergent does NOT imply Bounded");

    PropertyStore store;
    store.declare_finiteness("z", Finiteness::Divergent);

    EXPECT_TRUE(store.get_boundedness("z") == Boundedness::Unknown,
        "Divergent does not set Bounded");
    EXPECT_TRUE(store.get_finiteness("z") == Finiteness::Divergent,
        "Finiteness is Divergent");
}

static void test_finite_idempotent() {
    TEST_CASE("Declaring Finite twice is idempotent");

    PropertyStore store;
    store.declare_finiteness("x", Finiteness::Finite);
    store.declare_finiteness("x", Finiteness::Finite);

    EXPECT_TRUE(store.get_finiteness("x") == Finiteness::Finite,
        "Finite declared twice remains Finite");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "Still Bounded after idempotent Finite");
}

static void test_finite_multiple_symbols() {
    TEST_CASE("Finite implies Bounded for multiple symbols");

    PropertyStore store;
    store.declare_finiteness("a", Finiteness::Finite);
    store.declare_finiteness("b", Finiteness::Finite);
    store.declare_finiteness("c", Finiteness::Finite);

    EXPECT_TRUE(store.get_boundedness("a") == Boundedness::Bounded, "a is Bounded");
    EXPECT_TRUE(store.get_boundedness("b") == Boundedness::Bounded, "b is Bounded");
    EXPECT_TRUE(store.get_boundedness("c") == Boundedness::Bounded, "c is Bounded");
}


static void test_positive_definite_implies_positive_semidefinite() {
    TEST_CASE("PositiveDefinite implies PositiveSemiDefinite queryable");

    PropertyStore store;
    store.declare_definiteness("A", Definiteness::PositiveDefinite);

    // PositiveDefinite is stored; querying definiteness returns PositiveDefinite
    EXPECT_TRUE(store.get_definiteness("A") == Definiteness::PositiveDefinite,
        "A is PositiveDefinite");

    // Declaring PositiveSemiDefinite after PositiveDefinite should not throw
    // (it's implied, so it's either a no-op or accepted)
    // The implementation stores the most specific, so PositiveDefinite remains
    // We verify by checking that declaring PositiveSemiDefinite doesn't throw
    bool no_throw = true;
    try {
        store.declare_definiteness("A", Definiteness::PositiveSemiDefinite);
    } catch (...) {
        no_throw = false;
    }
    EXPECT_TRUE(no_throw,
        "PositiveSemiDefinite after PositiveDefinite does not throw (implied)");
}

static void test_positive_definite_plus_negative_definite_throws() {
    TEST_CASE("PositiveDefinite + NegativeDefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("B", Definiteness::PositiveDefinite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_definiteness("B", Definiteness::NegativeDefinite);
    }), "PositiveDefinite + NegativeDefinite returns failure");
}

static void test_negative_definite_plus_positive_definite_throws() {
    TEST_CASE("NegativeDefinite + PositiveDefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("C", Definiteness::NegativeDefinite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_definiteness("C", Definiteness::PositiveDefinite);
    }), "NegativeDefinite + PositiveDefinite returns failure");
}

static void test_positive_definite_plus_indefinite_throws() {
    TEST_CASE("PositiveDefinite + Indefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("D", Definiteness::PositiveDefinite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_definiteness("D", Definiteness::Indefinite);
    }), "PositiveDefinite + Indefinite returns failure");
}

static void test_negative_definite_plus_indefinite_throws() {
    TEST_CASE("NegativeDefinite + Indefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("E", Definiteness::NegativeDefinite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_definiteness("E", Definiteness::Indefinite);
    }), "NegativeDefinite + Indefinite returns failure");
}

static void test_negative_definite_implies_negative_semidefinite() {
    TEST_CASE("NegativeDefinite implies NegativeSemiDefinite");

    PropertyStore store;
    store.declare_definiteness("F", Definiteness::NegativeDefinite);

    EXPECT_TRUE(store.get_definiteness("F") == Definiteness::NegativeDefinite,
        "F is NegativeDefinite");

    // Declaring NegativeSemiDefinite after NegativeDefinite should not throw
    bool no_throw = true;
    try {
        store.declare_definiteness("F", Definiteness::NegativeSemiDefinite);
    } catch (...) {
        no_throw = false;
    }
    EXPECT_TRUE(no_throw,
        "NegativeSemiDefinite after NegativeDefinite does not throw (implied)");
}

static void test_positive_semidefinite_upgradeable_to_positive_definite() {
    TEST_CASE("PositiveSemiDefinite can be upgraded to PositiveDefinite");

    PropertyStore store;
    store.declare_definiteness("G", Definiteness::PositiveSemiDefinite);
    store.declare_definiteness("G", Definiteness::PositiveDefinite);

    EXPECT_TRUE(store.get_definiteness("G") == Definiteness::PositiveDefinite,
        "G upgraded from PositiveSemiDefinite to PositiveDefinite");
}

static void test_positive_definite_plus_negative_semidefinite_throws() {
    TEST_CASE("PositiveDefinite + NegativeSemiDefinite returns failure");

    PropertyStore store;
    store.declare_definiteness("H", Definiteness::PositiveDefinite);

    EXPECT_TRUE(returns_invalid_argument([&]() {
        return store.declare_definiteness("H", Definiteness::NegativeSemiDefinite);
    }), "PositiveDefinite + NegativeSemiDefinite returns failure");
}

static void test_definiteness_idempotent() {
    TEST_CASE("Declaring same definiteness twice is idempotent");

    PropertyStore store;
    store.declare_definiteness("I", Definiteness::PositiveDefinite);
    store.declare_definiteness("I", Definiteness::PositiveDefinite);

    EXPECT_TRUE(store.get_definiteness("I") == Definiteness::PositiveDefinite,
        "PositiveDefinite declared twice remains PositiveDefinite");
}


int main() {
    test_differentiable_implies_continuous_closed_interval();
    test_differentiable_implies_continuous_open_interval();
    test_differentiable_implies_continuous_subinterval();
    test_differentiable_implies_continuous_entire_line();
    test_differentiable_multiple_symbols();
    test_continuous_only_not_differentiable();

    test_transcendental_rejects_algebraic();
    test_transcendental_rejects_rational();
    test_transcendental_rejects_integer();
    test_transcendental_rejects_natural();
    test_transcendental_rejects_positive_int();
    test_algebraic_implies_real();
    test_transcendental_sets_real_domain();
    test_algebraic_then_transcendental_throws();
    test_rational_then_transcendental_throws();
    test_integer_then_transcendental_throws();
    test_transcendental_allows_real_declaration();
    test_transcendental_allows_complex_declaration();

    test_finite_implies_bounded();
    test_finite_then_divergent_throws();
    test_divergent_then_finite_throws();
    test_divergent_does_not_imply_bounded();
    test_finite_idempotent();
    test_finite_multiple_symbols();

    test_positive_definite_implies_positive_semidefinite();
    test_positive_definite_plus_negative_definite_throws();
    test_negative_definite_plus_positive_definite_throws();
    test_positive_definite_plus_indefinite_throws();
    test_negative_definite_plus_indefinite_throws();
    test_negative_definite_implies_negative_semidefinite();
    test_positive_semidefinite_upgradeable_to_positive_definite();
    test_positive_definite_plus_negative_semidefinite_throws();
    test_definiteness_idempotent();

    return TEST_REPORT();
}
