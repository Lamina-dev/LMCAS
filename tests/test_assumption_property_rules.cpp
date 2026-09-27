
#include "test_interval_support.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include <string>
#include <stdexcept>
#include <vector>

using namespace LMCAS;

/// Check that a callable returns an InvalidArgument failure.
template <typename F>
static bool returns_invalid_argument(F &&f) {
    auto result = f();
    return !result.has_value() && result.error().code == CasErrc::InvalidArgument;
}

TEST(AssumptionPropertyRules, DifferentiableImpliesContinuousClosedInterval) {
    PropertyStore store;
    Interval iv = make_closed_interval(0.0, 1.0);

    EXPECT_TRUE((store.declare_differentiable("f", iv).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("f", iv).value())) << "Differentiable on [0,1] => continuous on [0,1]";
    EXPECT_TRUE((store.is_differentiable("f", iv).value())) << "Differentiable on [0,1] => differentiable on [0,1]";
}

TEST(AssumptionPropertyRules, DifferentiableImpliesContinuousOpenInterval) {
    PropertyStore store;
    Interval iv = make_open_interval(0.0, 10.0);

    EXPECT_TRUE((store.declare_differentiable("g", iv).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("g", iv).value())) << "Differentiable on (0,10) => continuous on (0,10)";
}

TEST(AssumptionPropertyRules, DifferentiableImpliesContinuousSubinterval) {
    PropertyStore store;
    Interval outer = make_closed_interval(0.0, 10.0);
    Interval inner = make_closed_interval(2.0, 5.0);

    EXPECT_TRUE((store.declare_differentiable("h", outer).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("h", inner).value())) << "Differentiable on [0,10] => continuous on sub-interval [2,5]";
    EXPECT_TRUE((store.is_differentiable("h", inner).value())) << "Differentiable on [0,10] => differentiable on sub-interval [2,5]";
}

TEST(AssumptionPropertyRules, DifferentiableImpliesContinuousEntireLine) {
    PropertyStore store;
    Interval entire = Interval::entire_line();

    EXPECT_TRUE((store.declare_differentiable("p", entire).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("p", entire).value())) << "Differentiable on (-inf,+inf) => continuous on (-inf,+inf)";

    // Also continuous on any finite sub-interval
    Interval sub = make_closed_interval(-100.0, 100.0);
    EXPECT_TRUE((store.is_continuous("p", sub).value())) << "Differentiable on (-inf,+inf) => continuous on [-100,100]";
}

TEST(AssumptionPropertyRules, DifferentiableMultipleSymbols) {
    PropertyStore store;
    Interval iv1 = make_closed_interval(0.0, 1.0);
    Interval iv2 = make_closed_interval(-5.0, 5.0);
    Interval iv3 = make_closed_interval(10.0, 20.0);

    EXPECT_TRUE((store.declare_differentiable("a", iv1).has_value())) << "differentiability declaration succeeds";
    EXPECT_TRUE((store.declare_differentiable("b", iv2).has_value())) << "differentiability declaration succeeds";
    EXPECT_TRUE((store.declare_differentiable("c", iv3).has_value())) << "differentiability declaration succeeds";

    EXPECT_TRUE((store.is_continuous("a", iv1).value())) << "a continuous on [0,1]";
    EXPECT_TRUE((store.is_continuous("b", iv2).value())) << "b continuous on [-5,5]";
    EXPECT_TRUE((store.is_continuous("c", iv3).value())) << "c continuous on [10,20]";
}

TEST(AssumptionPropertyRules, ContinuousOnlyNotDifferentiable) {
    PropertyStore store;
    Interval iv = make_closed_interval(0.0, 1.0);

    EXPECT_TRUE((store.declare_continuous("f", iv).has_value())) << "continuity declaration succeeds";

    EXPECT_TRUE((store.is_continuous("f", iv).value())) << "Continuous on [0,1] => is_continuous true";
    EXPECT_FALSE((store.is_differentiable("f", iv).value())) << "Continuous-only on [0,1] => is_differentiable false";
}

TEST(AssumptionPropertyRules, TranscendentalRejectsAlgebraic) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("pi").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_domain("pi", Domain::Algebraic);
    }))) << "Transcendental + Algebraic returns failure";
}

TEST(AssumptionPropertyRules, TranscendentalRejectsRational) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("e").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_domain("e", Domain::Rational);
    }))) << "Transcendental + Rational returns failure";
}

TEST(AssumptionPropertyRules, TranscendentalRejectsInteger) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("tau").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_domain("tau", Domain::Integer);
    }))) << "Transcendental + Integer returns failure";
}

TEST(AssumptionPropertyRules, TranscendentalRejectsNatural) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("alpha").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_domain("alpha", Domain::Natural);
    }))) << "Transcendental + Natural returns failure";
}

TEST(AssumptionPropertyRules, TranscendentalRejectsPositiveInt) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("gamma").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_domain("gamma", Domain::PositiveInt);
    }))) << "Transcendental + PositiveInt returns failure";
}

TEST(AssumptionPropertyRules, AlgebraicImpliesReal) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Algebraic).has_value())) << "domain declaration succeeds";

    EXPECT_TRUE((store.has_domain("x", Domain::Real))) << "Algebraic implies Real";
    EXPECT_TRUE((store.has_domain("x", Domain::Complex))) << "Algebraic implies Complex";
}

TEST(AssumptionPropertyRules, TranscendentalSetsRealDomain) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("pi").has_value())) << "transcendental declaration succeeds";

    EXPECT_TRUE((store.get_domain("pi") == Domain::Real)) << "Transcendental symbol has Real domain";
    EXPECT_TRUE((store.is_transcendental("pi"))) << "Symbol is marked transcendental";
}

TEST(AssumptionPropertyRules, AlgebraicThenTranscendentalThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Algebraic).has_value())) << "domain declaration succeeds";

    // Transcendental requires domain <= Real, but Algebraic is more specific
    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_transcendental("x");
    }))) << "Algebraic then Transcendental returns failure";
}

TEST(AssumptionPropertyRules, RationalThenTranscendentalThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Rational).has_value())) << "domain declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_transcendental("x");
    }))) << "Rational then Transcendental returns failure";
}

TEST(AssumptionPropertyRules, IntegerThenTranscendentalThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("n", Domain::Integer).has_value())) << "domain declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_transcendental("n");
    }))) << "Integer then Transcendental returns failure";
}

TEST(AssumptionPropertyRules, TranscendentalAllowsRealDeclaration) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("pi").has_value())) << "transcendental declaration succeeds";

    // Declaring Real on a transcendental symbol should be a no-op (already Real)
    EXPECT_TRUE((store.declare_domain("pi", Domain::Real).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("pi") == Domain::Real)) << "Transcendental symbol remains Real after Real declaration";
    EXPECT_TRUE((store.is_transcendental("pi"))) << "Still transcendental after Real declaration";
}

TEST(AssumptionPropertyRules, TranscendentalAllowsComplexDeclaration) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_transcendental("e").has_value())) << "transcendental declaration succeeds";

    // Complex is less specific than Real, so it's a no-op
    EXPECT_TRUE((store.declare_domain("e", Domain::Complex).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("e") == Domain::Real)) << "Transcendental symbol remains Real after Complex declaration";
}

TEST(AssumptionPropertyRules, FiniteImpliesBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "Finite => Bounded";
    EXPECT_TRUE((store.get_finiteness("x") == Finiteness::Finite)) << "Finiteness is Finite";
}

TEST(AssumptionPropertyRules, FiniteThenDivergentThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_finiteness("x", Finiteness::Divergent);
    }))) << "Finite + Divergent returns failure";
}

TEST(AssumptionPropertyRules, DivergentThenFiniteThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("y", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_finiteness("y", Finiteness::Finite);
    }))) << "Divergent + Finite returns failure";
}

TEST(AssumptionPropertyRules, DivergentDoesNotImplyBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("z", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("z") == Boundedness::Unknown)) << "Divergent does not set Bounded";
    EXPECT_TRUE((store.get_finiteness("z") == Finiteness::Divergent)) << "Finiteness is Divergent";
}

TEST(AssumptionPropertyRules, FiniteIdempotent) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((store.declare_finiteness("x", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((store.get_finiteness("x") == Finiteness::Finite)) << "Finite declared twice remains Finite";
    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "Still Bounded after idempotent Finite";
}

TEST(AssumptionPropertyRules, FiniteMultipleSymbols) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_finiteness("a", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((store.declare_finiteness("b", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((store.declare_finiteness("c", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("a") == Boundedness::Bounded)) << "a is Bounded";
    EXPECT_TRUE((store.get_boundedness("b") == Boundedness::Bounded)) << "b is Bounded";
    EXPECT_TRUE((store.get_boundedness("c") == Boundedness::Bounded)) << "c is Bounded";
}
