
#include "test_common.hpp"
#include "property_store.hpp"
#include "symbolic.hpp"
#include <stdexcept>

using namespace LMCAS;

void test_positive_implies_nonnegative_and_nonzero() {
    TEST_CASE("Positive implies NonNegative and NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    EXPECT_TRUE(store.has_sign("x", Sign::Positive), "x has Positive");
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative (implied)");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero (implied)");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative), "x does not have Negative");
    EXPECT_FALSE(store.has_sign("x", Sign::NonPositive), "x does not have NonPositive");
    EXPECT_FALSE(store.has_sign("x", Sign::Zero), "x does not have Zero");
}

void test_negative_implies_nonpositive_and_nonzero() {
    TEST_CASE("Negative implies NonPositive and NonZero");
    PropertyStore store;
    store.declare_sign("y", Sign::Negative);

    EXPECT_TRUE(store.has_sign("y", Sign::Negative), "y has Negative");
    EXPECT_TRUE(store.has_sign("y", Sign::NonPositive), "y has NonPositive (implied)");
    EXPECT_TRUE(store.has_sign("y", Sign::NonZero), "y has NonZero (implied)");
    EXPECT_FALSE(store.has_sign("y", Sign::Positive), "y does not have Positive");
    EXPECT_FALSE(store.has_sign("y", Sign::NonNegative), "y does not have NonNegative");
    EXPECT_FALSE(store.has_sign("y", Sign::Zero), "y does not have Zero");
}

void test_zero_implies_nonnegative_nonpositive_and_integer() {
    TEST_CASE("Zero implies NonNegative, NonPositive, and Integer domain");
    PropertyStore store;
    store.declare_sign("z", Sign::Zero);

    EXPECT_TRUE(store.has_sign("z", Sign::Zero), "z has Zero");
    EXPECT_TRUE(store.has_sign("z", Sign::NonNegative), "z has NonNegative (implied)");
    EXPECT_TRUE(store.has_sign("z", Sign::NonPositive), "z has NonPositive (implied)");
    EXPECT_FALSE(store.has_sign("z", Sign::Positive), "z does not have Positive");
    EXPECT_FALSE(store.has_sign("z", Sign::Negative), "z does not have Negative");
    EXPECT_FALSE(store.has_sign("z", Sign::NonZero), "z does not have NonZero");

    // Zero also implies Integer domain
    EXPECT_TRUE(store.has_domain("z", Domain::Integer), "z has Integer domain (implied by Zero)");
}

void test_nonnegative_no_extra_implications() {
    TEST_CASE("NonNegative has no further implications");
    PropertyStore store;
    store.declare_sign("a", Sign::NonNegative);

    EXPECT_TRUE(store.has_sign("a", Sign::NonNegative), "a has NonNegative");
    EXPECT_FALSE(store.has_sign("a", Sign::Positive), "a does not have Positive");
    EXPECT_FALSE(store.has_sign("a", Sign::NonPositive), "a does not have NonPositive");
    EXPECT_FALSE(store.has_sign("a", Sign::NonZero), "a does not have NonZero");
    EXPECT_FALSE(store.has_sign("a", Sign::Zero), "a does not have Zero");
    EXPECT_FALSE(store.has_sign("a", Sign::Negative), "a does not have Negative");
}

void test_nonpositive_no_extra_implications() {
    TEST_CASE("NonPositive has no further implications");
    PropertyStore store;
    store.declare_sign("b", Sign::NonPositive);

    EXPECT_TRUE(store.has_sign("b", Sign::NonPositive), "b has NonPositive");
    EXPECT_FALSE(store.has_sign("b", Sign::Negative), "b does not have Negative");
    EXPECT_FALSE(store.has_sign("b", Sign::NonNegative), "b does not have NonNegative");
    EXPECT_FALSE(store.has_sign("b", Sign::NonZero), "b does not have NonZero");
    EXPECT_FALSE(store.has_sign("b", Sign::Zero), "b does not have Zero");
    EXPECT_FALSE(store.has_sign("b", Sign::Positive), "b does not have Positive");
}

void test_nonzero_no_extra_implications() {
    TEST_CASE("NonZero has no further implications");
    PropertyStore store;
    store.declare_sign("c", Sign::NonZero);

    EXPECT_TRUE(store.has_sign("c", Sign::NonZero), "c has NonZero");
    EXPECT_FALSE(store.has_sign("c", Sign::Positive), "c does not have Positive");
    EXPECT_FALSE(store.has_sign("c", Sign::Negative), "c does not have Negative");
    EXPECT_FALSE(store.has_sign("c", Sign::NonNegative), "c does not have NonNegative");
    EXPECT_FALSE(store.has_sign("c", Sign::NonPositive), "c does not have NonPositive");
    EXPECT_FALSE(store.has_sign("c", Sign::Zero), "c does not have Zero");
}

void test_idempotent_redeclaration() {
    TEST_CASE("Idempotent re-declaration (same sign already present → no-op)");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    // Re-declaring Positive should be a no-op (no exception)
    store.declare_sign("x", Sign::Positive);
    EXPECT_TRUE(store.has_sign("x", Sign::Positive), "x still has Positive after re-declaration");
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x still has NonNegative after re-declaration");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x still has NonZero after re-declaration");

    // Re-declaring an implied sign should also be a no-op
    store.declare_sign("x", Sign::NonNegative);
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x still has NonNegative after implied re-declaration");

    store.declare_sign("x", Sign::NonZero);
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x still has NonZero after implied re-declaration");
}

void test_contradiction_positive_negative() {
    TEST_CASE("Contradiction: Positive + Negative");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    auto failure_112 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_112.has_value(), "Declaring Negative after Positive returns failure");
    // State should be unchanged
    EXPECT_TRUE(store.has_sign("x", Sign::Positive), "x still has Positive after failed declaration");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative), "x does not have Negative after failed declaration");
}

void test_contradiction_positive_zero() {
    TEST_CASE("Contradiction: Positive + Zero");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    auto failure_129 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(!failure_129.has_value(), "Declaring Zero after Positive returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::Zero), "x does not have Zero after failed declaration");
}

void test_contradiction_positive_nonpositive() {
    TEST_CASE("Contradiction: Positive + NonPositive");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    auto failure_144 = store.declare_sign("x", Sign::NonPositive);
    EXPECT_TRUE(!failure_144.has_value(), "Declaring NonPositive after Positive returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::NonPositive), "x does not have NonPositive after failed declaration");
}

void test_contradiction_negative_zero() {
    TEST_CASE("Contradiction: Negative + Zero");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    auto failure_159 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(!failure_159.has_value(), "Declaring Zero after Negative returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::Zero), "x does not have Zero after failed declaration");
}

void test_contradiction_negative_nonnegative() {
    TEST_CASE("Contradiction: Negative + NonNegative");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    auto failure_174 = store.declare_sign("x", Sign::NonNegative);
    EXPECT_TRUE(!failure_174.has_value(), "Declaring NonNegative after Negative returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::NonNegative), "x does not have NonNegative after failed declaration");
}

void test_contradiction_nonnegative_negative() {
    TEST_CASE("Contradiction: NonNegative + Negative");
    PropertyStore store;
    store.declare_sign("x", Sign::NonNegative);

    auto failure_189 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_189.has_value(), "Declaring Negative after NonNegative returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative), "x does not have Negative after failed declaration");
}

void test_contradiction_nonpositive_positive() {
    TEST_CASE("Contradiction: NonPositive + Positive");
    PropertyStore store;
    store.declare_sign("x", Sign::NonPositive);

    auto failure_204 = store.declare_sign("x", Sign::Positive);
    EXPECT_TRUE(!failure_204.has_value(), "Declaring Positive after NonPositive returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::Positive), "x does not have Positive after failed declaration");
}

void test_contradiction_zero_nonzero() {
    TEST_CASE("Contradiction: Zero + NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);

    auto failure_219 = store.declare_sign("x", Sign::NonZero);
    EXPECT_TRUE(!failure_219.has_value(), "Declaring NonZero after Zero returns failure");
    EXPECT_FALSE(store.has_sign("x", Sign::NonZero), "x does not have NonZero after failed declaration");
}

void test_contradiction_via_implied_signs() {
    TEST_CASE("Contradiction via implied signs: Positive (implies NonNegative) then Negative");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    // x now has: Positive, NonNegative, NonZero

    // Declaring Negative should fail because Negative contradicts NonNegative (implied)
    auto failure_237 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_237.has_value(), "Declaring Negative contradicts implied NonNegative from Positive");
}

void test_contradiction_implied_against_existing() {
    TEST_CASE("Contradiction: new sign's implications conflict with existing signs");
    PropertyStore store;
    store.declare_sign("x", Sign::NonZero);
    // x now has: NonZero

    // Declaring Zero should fail because Zero contradicts NonZero
    auto failure_253 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(!failure_253.has_value(), "Declaring Zero contradicts existing NonZero");
}

void test_compatible_signs_can_coexist() {
    TEST_CASE("Compatible signs can coexist: NonNegative then NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::NonNegative);
    store.declare_sign("x", Sign::NonZero);

    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero");
}

void test_positive_after_nonnegative() {
    TEST_CASE("Positive after NonNegative is compatible (refines)");
    PropertyStore store;
    store.declare_sign("x", Sign::NonNegative);
    store.declare_sign("x", Sign::Positive);

    EXPECT_TRUE(store.has_sign("x", Sign::Positive), "x has Positive");
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero (implied by Positive)");
}

void test_get_signs_returns_all_stored() {
    TEST_CASE("get_signs returns all stored signs (explicit + implied)");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    auto signs = store.get_signs("x");
    EXPECT_TRUE(signs.count(Sign::Positive) > 0, "get_signs includes Positive");
    EXPECT_TRUE(signs.count(Sign::NonNegative) > 0, "get_signs includes NonNegative");
    EXPECT_TRUE(signs.count(Sign::NonZero) > 0, "get_signs includes NonZero");
    EXPECT_TRUE(signs.size() == 3, "get_signs has exactly 3 signs for Positive");
}

void test_undeclared_symbol_has_no_signs() {
    TEST_CASE("Undeclared symbol has no signs");
    PropertyStore store;

    EXPECT_FALSE(store.has_sign("unknown", Sign::Positive), "unknown has no Positive");
    EXPECT_FALSE(store.has_sign("unknown", Sign::Negative), "unknown has no Negative");
    EXPECT_FALSE(store.has_sign("unknown", Sign::NonNegative), "unknown has no NonNegative");
    EXPECT_FALSE(store.has_sign("unknown", Sign::NonPositive), "unknown has no NonPositive");
    EXPECT_FALSE(store.has_sign("unknown", Sign::Zero), "unknown has no Zero");
    EXPECT_FALSE(store.has_sign("unknown", Sign::NonZero), "unknown has no NonZero");

    auto signs = store.get_signs("unknown");
    EXPECT_TRUE(signs.empty(), "get_signs returns empty set for undeclared symbol");
}

void test_zero_domain_promotion_does_not_override_more_specific() {
    TEST_CASE("Zero domain promotion to Integer does not override more specific domain");
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);
    store.declare_sign("x", Sign::Zero);

    // Natural is more specific than Integer, so domain should remain Natural
    EXPECT_TRUE(store.get_domain("x") == Domain::Natural, "Domain remains Natural (more specific than Integer)");
    EXPECT_TRUE(store.has_sign("x", Sign::Zero), "x has Zero");
}

void test_checked_property_store_contracts() {
    TEST_CASE("PropertyStore checked declarations: explicit errors and transactional state");

    PropertyStore store;

    auto bad_symbol = store.declare_sign_checked("", Sign::Positive);
    EXPECT_TRUE(!bad_symbol.has_value(), "checked declare_sign rejects empty symbol");
    EXPECT_TRUE(bad_symbol.error().code == CasErrc::InvalidArgument,
                "checked declare_sign reports InvalidArgument for empty symbol");
    EXPECT_FALSE(store.has_sign("", Sign::Positive),
                 "failed checked declare_sign does not create empty-symbol fact");

    auto positive = store.declare_sign_checked("x", Sign::Positive);
    EXPECT_TRUE(positive.has_value(), "checked declare_sign succeeds");
    EXPECT_TRUE(store.has_sign("x", Sign::Positive),
                "checked declare_sign stores declared sign");

    auto contradiction = store.declare_sign_checked("x", Sign::Negative);
    EXPECT_TRUE(!contradiction.has_value(), "checked declare_sign rejects contradiction");
    EXPECT_TRUE(contradiction.error().code == CasErrc::InvalidArgument,
                "checked declare_sign reports InvalidArgument for contradiction");
    EXPECT_TRUE(store.has_sign("x", Sign::Positive),
                "failed checked declare_sign preserves previous sign");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative),
                 "failed checked declare_sign does not apply conflicting sign");

    auto domain = store.declare_domain_checked("n", Domain::Integer);
    EXPECT_TRUE(domain.has_value(), "checked declare_domain succeeds");
    EXPECT_TRUE(store.has_domain("n", Domain::Integer),
                "checked declare_domain stores domain");

    auto parity = store.declare_parity_checked("n", Parity::Even);
    EXPECT_TRUE(parity.has_value(), "checked declare_parity succeeds");
    auto parity_conflict = store.declare_parity_checked("n", Parity::Odd);
    EXPECT_TRUE(!parity_conflict.has_value(), "checked declare_parity rejects contradiction");
    EXPECT_TRUE(parity_conflict.error().code == CasErrc::InvalidArgument,
                "checked declare_parity reports InvalidArgument for contradiction");
    EXPECT_TRUE(store.get_parity("n") == Parity::Even,
                "failed checked declare_parity preserves previous parity");

    auto bounded = store.declare_bounded_checked("f", Boundedness::Bounded);
    EXPECT_TRUE(bounded.has_value(), "checked declare_bounded succeeds");
    auto bounded_conflict = store.declare_bounded_checked("f", Boundedness::Unbounded);
    EXPECT_TRUE(!bounded_conflict.has_value(),
                "checked declare_bounded rejects contradiction");
    EXPECT_TRUE(store.get_boundedness("f") == Boundedness::Bounded,
                "failed checked declare_bounded preserves previous boundedness");

    auto transcendental = store.declare_transcendental_checked("t");
    EXPECT_TRUE(transcendental.has_value(), "checked declare_transcendental succeeds");
    auto algebraic_conflict = store.declare_domain_checked("t", Domain::Algebraic);
    EXPECT_TRUE(!algebraic_conflict.has_value(),
                "checked declare_domain rejects transcendental/algebraic conflict");
    EXPECT_TRUE(store.is_transcendental("t"),
                "failed checked declare_domain preserves transcendental marker");

    auto finite = store.declare_finiteness_checked("seq", Finiteness::Finite);
    EXPECT_TRUE(finite.has_value(), "checked declare_finiteness succeeds");
    auto divergent = store.declare_finiteness_checked("seq", Finiteness::Divergent);
    EXPECT_TRUE(!divergent.has_value(),
                "checked declare_finiteness rejects contradiction");
    EXPECT_TRUE(store.get_finiteness("seq") == Finiteness::Finite,
                "failed checked declare_finiteness preserves previous finiteness");

    auto pd = store.declare_definiteness_checked("A", Definiteness::PositiveDefinite);
    EXPECT_TRUE(pd.has_value(), "checked declare_definiteness succeeds");
    auto indefinite = store.declare_definiteness_checked("A", Definiteness::Indefinite);
    EXPECT_TRUE(!indefinite.has_value(),
                "checked declare_definiteness rejects contradiction");
    EXPECT_TRUE(store.get_definiteness("A") == Definiteness::PositiveDefinite,
                "failed checked declare_definiteness preserves previous definiteness");

    auto null_period = store.declare_periodic_checked("g", nullptr);
    EXPECT_TRUE(!null_period.has_value(), "checked declare_periodic rejects null period");
    EXPECT_TRUE(null_period.error().code == CasErrc::InvalidArgument,
                "checked declare_periodic reports InvalidArgument for null period");
    EXPECT_FALSE(store.is_periodic("g"),
                 "failed checked declare_periodic does not mark symbol periodic");

    auto period = store.declare_periodic_checked("g", SymbolicExpr::number(2));
    EXPECT_TRUE(period.has_value(), "checked declare_periodic succeeds");
    EXPECT_TRUE(store.is_periodic("g"), "checked declare_periodic stores period");
}

void test_legacy_declarations_delegate_transactionally() {
    TEST_CASE("PropertyStore canonical declarations delegate transactionally");

    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    const auto symbols_before = store.get_all_symbols();

    auto contradiction = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!contradiction.has_value(),
                "canonical sign declaration returns contradiction");
    EXPECT_TRUE(contradiction.error().code == CasErrc::InvalidArgument,
                "sign contradiction reports InvalidArgument");
    EXPECT_TRUE(store.has_sign("x", Sign::Positive),
                "canonical sign failure preserves the established sign");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative),
                 "canonical sign failure does not commit a contradictory sign");

    auto empty_symbol = store.declare_domain("", Domain::Real);
    EXPECT_TRUE(!empty_symbol.has_value(),
                "canonical domain declaration rejects an empty symbol");
    EXPECT_TRUE(empty_symbol.error().code == CasErrc::InvalidArgument,
                "empty domain symbol reports InvalidArgument");
    EXPECT_TRUE(store.get_all_symbols() == symbols_before,
                "failed canonical declaration does not create an empty property record");
}

int main() {
    test_positive_implies_nonnegative_and_nonzero();
    test_negative_implies_nonpositive_and_nonzero();
    test_zero_implies_nonnegative_nonpositive_and_integer();
    test_nonnegative_no_extra_implications();
    test_nonpositive_no_extra_implications();
    test_nonzero_no_extra_implications();
    test_idempotent_redeclaration();
    test_contradiction_positive_negative();
    test_contradiction_positive_zero();
    test_contradiction_positive_nonpositive();
    test_contradiction_negative_zero();
    test_contradiction_negative_nonnegative();
    test_contradiction_nonnegative_negative();
    test_contradiction_nonpositive_positive();
    test_contradiction_zero_nonzero();
    test_contradiction_via_implied_signs();
    test_contradiction_implied_against_existing();
    test_compatible_signs_can_coexist();
    test_positive_after_nonnegative();
    test_get_signs_returns_all_stored();
    test_undeclared_symbol_has_no_signs();
    test_zero_domain_promotion_does_not_override_more_specific();
    test_checked_property_store_contracts();
    test_legacy_declarations_delegate_transactionally();

    return TEST_REPORT();
}
