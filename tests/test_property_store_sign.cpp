/**
 * @file test_property_store_sign.cpp
 * @brief Unit tests for PropertyStore sign declaration with implication and contradiction detection.
 *
 * Validates Requirements 2.1, 2.2, 2.3, 2.4, 2.6:
 * - Sign declaration stores the sign and all implied signs
 * - Contradiction pairs are detected and throw std::invalid_argument
 * - Idempotent re-declaration (same sign already present → no-op)
 * - Zero implies Integer domain
 */

#include "test_common.hpp"
#include "property_store.hpp"
#include <stdexcept>

using namespace lamina;

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

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Negative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Negative after Positive throws");
    // State should be unchanged
    EXPECT_TRUE(store.has_sign("x", Sign::Positive), "x still has Positive after failed declaration");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative), "x does not have Negative after failed declaration");
}

void test_contradiction_positive_zero() {
    TEST_CASE("Contradiction: Positive + Zero");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Zero);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Zero after Positive throws");
    EXPECT_FALSE(store.has_sign("x", Sign::Zero), "x does not have Zero after failed declaration");
}

void test_contradiction_positive_nonpositive() {
    TEST_CASE("Contradiction: Positive + NonPositive");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::NonPositive);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring NonPositive after Positive throws");
    EXPECT_FALSE(store.has_sign("x", Sign::NonPositive), "x does not have NonPositive after failed declaration");
}

void test_contradiction_negative_zero() {
    TEST_CASE("Contradiction: Negative + Zero");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Zero);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Zero after Negative throws");
    EXPECT_FALSE(store.has_sign("x", Sign::Zero), "x does not have Zero after failed declaration");
}

void test_contradiction_negative_nonnegative() {
    TEST_CASE("Contradiction: Negative + NonNegative");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::NonNegative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring NonNegative after Negative throws");
    EXPECT_FALSE(store.has_sign("x", Sign::NonNegative), "x does not have NonNegative after failed declaration");
}

void test_contradiction_nonnegative_negative() {
    TEST_CASE("Contradiction: NonNegative + Negative");
    PropertyStore store;
    store.declare_sign("x", Sign::NonNegative);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Negative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Negative after NonNegative throws");
    EXPECT_FALSE(store.has_sign("x", Sign::Negative), "x does not have Negative after failed declaration");
}

void test_contradiction_nonpositive_positive() {
    TEST_CASE("Contradiction: NonPositive + Positive");
    PropertyStore store;
    store.declare_sign("x", Sign::NonPositive);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Positive);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Positive after NonPositive throws");
    EXPECT_FALSE(store.has_sign("x", Sign::Positive), "x does not have Positive after failed declaration");
}

void test_contradiction_zero_nonzero() {
    TEST_CASE("Contradiction: Zero + NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::NonZero);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring NonZero after Zero throws");
    EXPECT_FALSE(store.has_sign("x", Sign::NonZero), "x does not have NonZero after failed declaration");
}

void test_contradiction_via_implied_signs() {
    TEST_CASE("Contradiction via implied signs: Positive (implies NonNegative) then Negative");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    // x now has: Positive, NonNegative, NonZero

    // Declaring Negative should fail because Negative contradicts NonNegative (implied)
    bool threw = false;
    try {
        store.declare_sign("x", Sign::Negative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Negative contradicts implied NonNegative from Positive");
}

void test_contradiction_implied_against_existing() {
    TEST_CASE("Contradiction: new sign's implications conflict with existing signs");
    PropertyStore store;
    store.declare_sign("x", Sign::NonZero);
    // x now has: NonZero

    // Declaring Zero should fail because Zero contradicts NonZero
    bool threw = false;
    try {
        store.declare_sign("x", Sign::Zero);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Declaring Zero contradicts existing NonZero");
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

    return TEST_REPORT();
}
