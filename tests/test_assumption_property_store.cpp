/**
 * @file test_assumption_property_store.cpp
 * @brief Unified property tests for PropertyStore (Properties 1-8).
 *
 * Feature: assumption-system
 * Validates: Requirements 1.1-1.5, 2.1-2.4, 2.6, 3.1-3.7
 *
 * Property 1: Domain declaration round-trip and idempotence
 * Property 2: Domain hierarchy implication
 * Property 3: Domain specificity preservation
 * Property 4: Domain-sign contradiction detection
 * Property 5: Sign declaration with implication
 * Property 6: Sign contradiction detection
 * Property 7: Parity declaration with domain auto-promotion
 * Property 8: Boundedness declaration consistency
 */

#include "test_common.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <utility>

using namespace lamina;

// ============================================================
// Helper: all Domain values for exhaustive iteration
// ============================================================

static const std::vector<Domain> ALL_DOMAINS = {
    Domain::Complex, Domain::Real, Domain::Algebraic, Domain::Rational,
    Domain::Integer, Domain::Natural, Domain::PositiveInt
};

static const std::vector<Sign> ALL_SIGNS = {
    Sign::Positive, Sign::Negative, Sign::NonNegative,
    Sign::NonPositive, Sign::Zero, Sign::NonZero
};

static const std::vector<std::string> TEST_SYMBOLS = {
    "x", "y", "alpha", "longVariableName123", "a_b_c"
};

// Domain specificity helper (mirrors PropertyStore internal logic)
static int domain_specificity(Domain d) {
    switch (d) {
        case Domain::Complex:     return 0;
        case Domain::Real:        return 1;
        case Domain::Algebraic:   return 2;
        case Domain::Rational:    return 3;
        case Domain::Integer:     return 4;
        case Domain::Natural:     return 5;
        case Domain::PositiveInt: return 6;
    }
    return 0;
}

static std::string domain_name(Domain d) {
    switch (d) {
        case Domain::Complex:     return "Complex";
        case Domain::Real:        return "Real";
        case Domain::Algebraic:   return "Algebraic";
        case Domain::Rational:    return "Rational";
        case Domain::Integer:     return "Integer";
        case Domain::Natural:     return "Natural";
        case Domain::PositiveInt: return "PositiveInt";
    }
    return "?";
}

static std::string sign_name(Sign s) {
    switch (s) {
        case Sign::Positive:    return "Positive";
        case Sign::Negative:    return "Negative";
        case Sign::NonNegative: return "NonNegative";
        case Sign::NonPositive: return "NonPositive";
        case Sign::Zero:        return "Zero";
        case Sign::NonZero:     return "NonZero";
    }
    return "?";
}

// ============================================================
// Property 1: Domain declaration round-trip and idempotence
// **Validates: Requirements 1.1, 1.5**
// ============================================================

void test_property1_domain_roundtrip_all() {
    TEST_CASE("Property 1: Domain declaration round-trip for all domains and symbols");
    // For any valid symbol and any Domain, declaring then querying returns
    // at least as specific as declared.
    for (const auto& sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            store.declare_domain(sym, d);
            Domain result = store.get_domain(sym);
            bool at_least_as_specific =
                domain_specificity(result) >= domain_specificity(d);
            EXPECT_TRUE(at_least_as_specific,
                sym + " domain after declaring " + domain_name(d) +
                " is at least as specific");
        }
    }
}

void test_property1_idempotence() {
    TEST_CASE("Property 1: Re-declaring same domain is idempotent");
    for (const auto& sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            store.declare_domain(sym, d);
            Domain first = store.get_domain(sym);
            // Re-declare same domain
            store.declare_domain(sym, d);
            Domain second = store.get_domain(sym);
            EXPECT_TRUE(first == second,
                sym + " domain unchanged after re-declaring " + domain_name(d));
        }
    }
}

// ============================================================
// Property 2: Domain hierarchy implication
// **Validates: Requirements 1.2**
// ============================================================

void test_property2_hierarchy_implication() {
    TEST_CASE("Property 2: Declaring domain D implies all ancestor domains return True");
    // For each domain D, after declaring it, has_domain should return true
    // for D and all less-specific (ancestor) domains.
    for (const auto& sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            store.declare_domain(sym, d);

            for (Domain ancestor : ALL_DOMAINS) {
                if (domain_specificity(ancestor) <= domain_specificity(d)) {
                    // ancestor is same or less specific -> should be true
                    EXPECT_TRUE(store.has_domain(sym, ancestor),
                        sym + " with " + domain_name(d) +
                        " has ancestor " + domain_name(ancestor));
                }
            }
        }
    }
}

void test_property2_non_ancestors_false() {
    TEST_CASE("Property 2: Domains more specific than declared return False");
    for (Domain d : ALL_DOMAINS) {
        PropertyStore store;
        store.declare_domain("x", d);

        for (Domain more_specific : ALL_DOMAINS) {
            if (domain_specificity(more_specific) > domain_specificity(d)) {
                EXPECT_FALSE(store.has_domain("x", more_specific),
                    "x with " + domain_name(d) +
                    " does NOT have " + domain_name(more_specific));
            }
        }
    }
}

// ============================================================
// Property 3: Domain specificity preservation
// **Validates: Requirements 1.3**
// ============================================================

void test_property3_specificity_preservation() {
    TEST_CASE("Property 3: Declaring less-specific domain after more-specific is no-op");
    // For each pair (D_specific, D_ancestor) where D_ancestor < D_specific,
    // declaring D_specific then D_ancestor should leave domain at D_specific.
    for (const auto& sym : TEST_SYMBOLS) {
        for (Domain specific : ALL_DOMAINS) {
            for (Domain ancestor : ALL_DOMAINS) {
                if (domain_specificity(ancestor) < domain_specificity(specific)) {
                    PropertyStore store;
                    store.declare_domain(sym, specific);
                    Domain before = store.get_domain(sym);

                    store.declare_domain(sym, ancestor);
                    Domain after = store.get_domain(sym);

                    EXPECT_TRUE(before == after,
                        sym + ": declaring " + domain_name(ancestor) +
                        " after " + domain_name(specific) + " is no-op");
                }
            }
        }
    }
}

void test_property3_more_specific_upgrades() {
    TEST_CASE("Property 3: Declaring more-specific domain upgrades");
    // Declaring a more-specific domain should upgrade.
    PropertyStore store;
    store.declare_domain("x", Domain::Real);
    EXPECT_TRUE(store.get_domain("x") == Domain::Real, "x starts at Real");

    store.declare_domain("x", Domain::Integer);
    EXPECT_TRUE(store.get_domain("x") == Domain::Integer,
        "x upgraded to Integer");

    store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(store.get_domain("x") == Domain::PositiveInt,
        "x upgraded to PositiveInt");
}

// ============================================================
// Property 4: Domain-sign contradiction detection
// **Validates: Requirements 1.4**
// ============================================================

void test_property4_natural_negative_contradiction() {
    TEST_CASE("Property 4: Natural + Negative sign throws");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    bool threw = false;
    try {
        store.declare_domain("x", Domain::Natural);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Natural domain with Negative sign throws");
    // State unchanged: domain should still be Complex (default)
    EXPECT_TRUE(store.get_domain("x") == Domain::Complex,
        "x domain unchanged after failed Natural declaration");
}

void test_property4_positiveint_negative_contradiction() {
    TEST_CASE("Property 4: PositiveInt + Negative sign throws");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    bool threw = false;
    try {
        store.declare_domain("x", Domain::PositiveInt);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "PositiveInt domain with Negative sign throws");
}

void test_property4_positiveint_zero_contradiction() {
    TEST_CASE("Property 4: PositiveInt + Zero sign throws");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);

    bool threw = false;
    try {
        store.declare_domain("x", Domain::PositiveInt);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "PositiveInt domain with Zero sign throws");
}

void test_property4_positiveint_nonpositive_contradiction() {
    TEST_CASE("Property 4: PositiveInt + NonPositive sign throws");
    PropertyStore store;
    store.declare_sign("x", Sign::NonPositive);

    bool threw = false;
    try {
        store.declare_domain("x", Domain::PositiveInt);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "PositiveInt domain with NonPositive sign throws");
}

void test_property4_natural_implied_negative_contradiction() {
    TEST_CASE("Property 4: Natural + implied Negative (from Negative declaration) throws");
    PropertyStore store;
    // Declaring Negative implies NonPositive and NonZero
    store.declare_sign("y", Sign::Negative);

    bool threw = false;
    try {
        store.declare_domain("y", Domain::Natural);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Natural domain with implied Negative sign throws");
}

void test_property4_compatible_domain_sign_no_throw() {
    TEST_CASE("Property 4: Compatible domain-sign pairs do not throw");
    // Natural + Positive is fine
    {
        PropertyStore store;
        store.declare_sign("a", Sign::Positive);
        store.declare_domain("a", Domain::Natural);
        EXPECT_TRUE(store.get_domain("a") == Domain::Natural,
            "Natural + Positive is compatible");
    }
    // Natural + NonNegative is fine
    {
        PropertyStore store;
        store.declare_sign("b", Sign::NonNegative);
        store.declare_domain("b", Domain::Natural);
        EXPECT_TRUE(store.get_domain("b") == Domain::Natural,
            "Natural + NonNegative is compatible");
    }
    // PositiveInt + Positive is fine
    {
        PropertyStore store;
        store.declare_sign("c", Sign::Positive);
        store.declare_domain("c", Domain::PositiveInt);
        EXPECT_TRUE(store.get_domain("c") == Domain::PositiveInt,
            "PositiveInt + Positive is compatible");
    }
    // Integer + Negative is fine
    {
        PropertyStore store;
        store.declare_sign("d", Sign::Negative);
        store.declare_domain("d", Domain::Integer);
        EXPECT_TRUE(store.get_domain("d") == Domain::Integer,
            "Integer + Negative is compatible");
    }
}

// ============================================================
// Property 5: Sign declaration with implication
// **Validates: Requirements 2.1, 2.2, 2.4**
// ============================================================

void test_property5_positive_implications() {
    TEST_CASE("Property 5: Positive implies NonNegative and NonZero");
    for (const auto& sym : TEST_SYMBOLS) {
        PropertyStore store;
        store.declare_sign(sym, Sign::Positive);
        EXPECT_TRUE(store.has_sign(sym, Sign::Positive),
            sym + " has Positive");
        EXPECT_TRUE(store.has_sign(sym, Sign::NonNegative),
            sym + " has NonNegative (implied)");
        EXPECT_TRUE(store.has_sign(sym, Sign::NonZero),
            sym + " has NonZero (implied)");
    }
}

void test_property5_negative_implications() {
    TEST_CASE("Property 5: Negative implies NonPositive and NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(store.has_sign("x", Sign::Negative), "x has Negative");
    EXPECT_TRUE(store.has_sign("x", Sign::NonPositive), "x has NonPositive (implied)");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero (implied)");
}

void test_property5_zero_implications() {
    TEST_CASE("Property 5: Zero implies NonNegative, NonPositive, and Integer domain");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(store.has_sign("x", Sign::Zero), "x has Zero");
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative (implied)");
    EXPECT_TRUE(store.has_sign("x", Sign::NonPositive), "x has NonPositive (implied)");
    EXPECT_TRUE(store.has_domain("x", Domain::Integer),
        "x has Integer domain (implied by Zero)");
}

void test_property5_nonneg_nonpos_nonzero_no_extra() {
    TEST_CASE("Property 5: NonNegative/NonPositive/NonZero have no further implications");
    {
        PropertyStore store;
        store.declare_sign("a", Sign::NonNegative);
        auto signs = store.get_signs("a");
        EXPECT_TRUE(signs.size() == 1, "NonNegative has no extra implications");
    }
    {
        PropertyStore store;
        store.declare_sign("b", Sign::NonPositive);
        auto signs = store.get_signs("b");
        EXPECT_TRUE(signs.size() == 1, "NonPositive has no extra implications");
    }
    {
        PropertyStore store;
        store.declare_sign("c", Sign::NonZero);
        auto signs = store.get_signs("c");
        EXPECT_TRUE(signs.size() == 1, "NonZero has no extra implications");
    }
}

void test_property5_redeclaration_noop() {
    TEST_CASE("Property 5: Re-declaring same sign is no-op");
    for (Sign s : ALL_SIGNS) {
        PropertyStore store;
        store.declare_sign("x", s);
        auto signs_before = store.get_signs("x");
        // Re-declare
        store.declare_sign("x", s);
        auto signs_after = store.get_signs("x");
        EXPECT_TRUE(signs_before == signs_after,
            "Re-declaring " + sign_name(s) + " is no-op");
    }
}

void test_property5_redeclare_implied_sign_noop() {
    TEST_CASE("Property 5: Re-declaring an implied sign is no-op");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    auto signs_before = store.get_signs("x");
    // NonNegative is implied by Positive; re-declaring should be no-op
    store.declare_sign("x", Sign::NonNegative);
    auto signs_after = store.get_signs("x");
    EXPECT_TRUE(signs_before == signs_after,
        "Re-declaring implied NonNegative after Positive is no-op");
}

// ============================================================
// Property 6: Sign contradiction detection
// **Validates: Requirements 2.3, 2.6**
// ============================================================

// Contradiction pairs to test exhaustively
static const std::vector<std::pair<Sign, Sign>> CONTRADICTION_PAIRS = {
    {Sign::Positive, Sign::Negative},
    {Sign::Positive, Sign::Zero},
    {Sign::Positive, Sign::NonPositive},
    {Sign::Negative, Sign::Zero},
    {Sign::Negative, Sign::NonNegative},
    {Sign::NonNegative, Sign::Negative},
    {Sign::NonPositive, Sign::Positive},
    {Sign::Zero, Sign::NonZero},
    {Sign::NonZero, Sign::Zero},
};

void test_property6_all_contradiction_pairs() {
    TEST_CASE("Property 6: All contradiction pairs throw std::invalid_argument");
    for (const auto& [s1, s2] : CONTRADICTION_PAIRS) {
        PropertyStore store;
        store.declare_sign("x", s1);

        bool threw = false;
        try {
            store.declare_sign("x", s2);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw,
            sign_name(s1) + " + " + sign_name(s2) + " throws");
        // State unchanged: s2 should not be present
        EXPECT_FALSE(store.has_sign("x", s2),
            "x does not have " + sign_name(s2) + " after failed declaration");
    }
}

void test_property6_implied_contradiction() {
    TEST_CASE("Property 6: Contradiction via implied signs");
    // Positive implies NonNegative; then declaring Negative should fail
    // because Negative contradicts NonNegative
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    bool threw = false;
    try {
        store.declare_sign("x", Sign::Negative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Negative contradicts implied NonNegative from Positive");
    EXPECT_TRUE(store.has_sign("x", Sign::Positive),
        "x still has Positive after failed Negative declaration");
}

void test_property6_compatible_signs_no_throw() {
    TEST_CASE("Property 6: Compatible sign pairs do not throw");
    // NonNegative + NonZero are compatible
    {
        PropertyStore store;
        store.declare_sign("x", Sign::NonNegative);
        store.declare_sign("x", Sign::NonZero);
        EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative");
        EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero");
    }
    // Positive after NonNegative (refines)
    {
        PropertyStore store;
        store.declare_sign("y", Sign::NonNegative);
        store.declare_sign("y", Sign::Positive);
        EXPECT_TRUE(store.has_sign("y", Sign::Positive), "y has Positive");
        EXPECT_TRUE(store.has_sign("y", Sign::NonNegative), "y has NonNegative");
    }
}

// ============================================================
// Property 7: Parity declaration with domain auto-promotion
// **Validates: Requirements 3.1, 3.3, 3.4, 3.7**
// ============================================================

void test_property7_even_promotes_to_integer() {
    TEST_CASE("Property 7: Even parity auto-promotes to Integer domain");
    for (const auto& sym : TEST_SYMBOLS) {
        PropertyStore store;
        // Default domain is Complex
        store.declare_parity(sym, Parity::Even);
        EXPECT_TRUE(store.get_parity(sym) == Parity::Even,
            sym + " has Even parity");
        EXPECT_TRUE(store.has_domain(sym, Domain::Integer),
            sym + " promoted to Integer by Even parity");
    }
}

void test_property7_odd_promotes_to_integer() {
    TEST_CASE("Property 7: Odd parity auto-promotes to Integer domain");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);
    EXPECT_TRUE(store.get_parity("x") == Parity::Odd, "x has Odd parity");
    EXPECT_TRUE(store.has_domain("x", Domain::Integer),
        "x promoted to Integer by Odd parity");
}

void test_property7_parity_does_not_demote() {
    TEST_CASE("Property 7: Parity does not demote more-specific domain");
    // Natural is more specific than Integer
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);
    store.declare_parity("x", Parity::Even);
    EXPECT_TRUE(store.get_domain("x") == Domain::Natural,
        "x domain remains Natural (more specific than Integer)");
}

void test_property7_parity_idempotent() {
    TEST_CASE("Property 7: Re-declaring same parity is no-op");
    PropertyStore store;
    store.declare_parity("x", Parity::Even);
    store.declare_parity("x", Parity::Even);
    EXPECT_TRUE(store.get_parity("x") == Parity::Even,
        "x still Even after re-declaration");

    PropertyStore store2;
    store2.declare_parity("y", Parity::Odd);
    store2.declare_parity("y", Parity::Odd);
    EXPECT_TRUE(store2.get_parity("y") == Parity::Odd,
        "y still Odd after re-declaration");
}

void test_property7_parity_contradiction() {
    TEST_CASE("Property 7: Contradictory parity throws");
    // Even then Odd
    {
        PropertyStore store;
        store.declare_parity("x", Parity::Even);
        bool threw = false;
        try {
            store.declare_parity("x", Parity::Odd);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "Even then Odd throws");
        EXPECT_TRUE(store.get_parity("x") == Parity::Even,
            "x parity unchanged after failed Odd declaration");
    }
    // Odd then Even
    {
        PropertyStore store;
        store.declare_parity("y", Parity::Odd);
        bool threw = false;
        try {
            store.declare_parity("y", Parity::Even);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "Odd then Even throws");
        EXPECT_TRUE(store.get_parity("y") == Parity::Odd,
            "y parity unchanged after failed Even declaration");
    }
}

void test_property7_default_parity_unknown() {
    TEST_CASE("Property 7: Default parity is Unknown");
    PropertyStore store;
    EXPECT_TRUE(store.get_parity("undeclared") == Parity::Unknown,
        "Undeclared symbol has Unknown parity");
}

// ============================================================
// Property 8: Boundedness declaration consistency
// **Validates: Requirements 3.2, 3.5, 3.7**
// ============================================================

void test_property8_bounded_stored() {
    TEST_CASE("Property 8: Declaring Bounded stores Bounded");
    for (const auto& sym : TEST_SYMBOLS) {
        PropertyStore store;
        store.declare_bounded(sym, Boundedness::Bounded);
        EXPECT_TRUE(store.get_boundedness(sym) == Boundedness::Bounded,
            sym + " has Bounded");
    }
}

void test_property8_unbounded_stored() {
    TEST_CASE("Property 8: Declaring Unbounded stores Unbounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
        "x has Unbounded");
}

void test_property8_bounded_with_interval() {
    TEST_CASE("Property 8: Bounded with interval stores bounds");
    PropertyStore store;

    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    store.declare_bounded("x", Boundedness::Bounded, bounds);
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "x has Bounded");
    EXPECT_TRUE(store.get_bounds("x").has_value(),
        "x has bounds stored");
}

void test_property8_idempotent() {
    TEST_CASE("Property 8: Re-declaring same boundedness is no-op");
    {
        PropertyStore store;
        store.declare_bounded("x", Boundedness::Bounded);
        store.declare_bounded("x", Boundedness::Bounded);
        EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
            "x still Bounded after re-declaration");
    }
    {
        PropertyStore store;
        store.declare_bounded("y", Boundedness::Unbounded);
        store.declare_bounded("y", Boundedness::Unbounded);
        EXPECT_TRUE(store.get_boundedness("y") == Boundedness::Unbounded,
            "y still Unbounded after re-declaration");
    }
}

void test_property8_contradiction_bounded_unbounded() {
    TEST_CASE("Property 8: Bounded then Unbounded throws");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    bool threw = false;
    try {
        store.declare_bounded("x", Boundedness::Unbounded);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Bounded then Unbounded throws");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "x remains Bounded after failed Unbounded declaration");
}

void test_property8_contradiction_unbounded_bounded() {
    TEST_CASE("Property 8: Unbounded then Bounded throws");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);

    bool threw = false;
    try {
        store.declare_bounded("x", Boundedness::Bounded);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Unbounded then Bounded throws");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
        "x remains Unbounded after failed Bounded declaration");
}

void test_property8_default_unknown() {
    TEST_CASE("Property 8: Default boundedness is Unknown");
    PropertyStore store;
    EXPECT_TRUE(store.get_boundedness("undeclared") == Boundedness::Unknown,
        "Undeclared symbol has Unknown boundedness");
    EXPECT_FALSE(store.get_bounds("undeclared").has_value(),
        "Undeclared symbol has no bounds");
}

// ============================================================
// Cross-property interaction tests
// ============================================================

void test_cross_domain_sign_interaction() {
    TEST_CASE("Cross: Domain Natural then sign Negative throws (reverse direction)");
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Negative);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Negative sign with Natural domain throws");
}

void test_cross_positiveint_then_zero_throws() {
    TEST_CASE("Cross: Domain PositiveInt then sign Zero throws");
    PropertyStore store;
    store.declare_domain("x", Domain::PositiveInt);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::Zero);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "Zero sign with PositiveInt domain throws");
}

void test_cross_positiveint_then_nonpositive_throws() {
    TEST_CASE("Cross: Domain PositiveInt then sign NonPositive throws");
    PropertyStore store;
    store.declare_domain("x", Domain::PositiveInt);

    bool threw = false;
    try {
        store.declare_sign("x", Sign::NonPositive);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "NonPositive sign with PositiveInt domain throws");
}

void test_default_domain_complex() {
    TEST_CASE("Default domain is Complex for undeclared symbol");
    PropertyStore store;
    EXPECT_TRUE(store.get_domain("undeclared") == Domain::Complex,
        "Undeclared symbol has Complex domain");
    EXPECT_TRUE(store.has_domain("undeclared", Domain::Complex),
        "has_domain returns true for Complex on undeclared symbol");
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 1: Domain declaration round-trip and idempotence
    test_property1_domain_roundtrip_all();
    test_property1_idempotence();

    // Property 2: Domain hierarchy implication
    test_property2_hierarchy_implication();
    test_property2_non_ancestors_false();

    // Property 3: Domain specificity preservation
    test_property3_specificity_preservation();
    test_property3_more_specific_upgrades();

    // Property 4: Domain-sign contradiction detection
    test_property4_natural_negative_contradiction();
    test_property4_positiveint_negative_contradiction();
    test_property4_positiveint_zero_contradiction();
    test_property4_positiveint_nonpositive_contradiction();
    test_property4_natural_implied_negative_contradiction();
    test_property4_compatible_domain_sign_no_throw();

    // Property 5: Sign declaration with implication
    test_property5_positive_implications();
    test_property5_negative_implications();
    test_property5_zero_implications();
    test_property5_nonneg_nonpos_nonzero_no_extra();
    test_property5_redeclaration_noop();
    test_property5_redeclare_implied_sign_noop();

    // Property 6: Sign contradiction detection
    test_property6_all_contradiction_pairs();
    test_property6_implied_contradiction();
    test_property6_compatible_signs_no_throw();

    // Property 7: Parity declaration with domain auto-promotion
    test_property7_even_promotes_to_integer();
    test_property7_odd_promotes_to_integer();
    test_property7_parity_does_not_demote();
    test_property7_parity_idempotent();
    test_property7_parity_contradiction();
    test_property7_default_parity_unknown();

    // Property 8: Boundedness declaration consistency
    test_property8_bounded_stored();
    test_property8_unbounded_stored();
    test_property8_bounded_with_interval();
    test_property8_idempotent();
    test_property8_contradiction_bounded_unbounded();
    test_property8_contradiction_unbounded_bounded();
    test_property8_default_unknown();

    // Cross-property interaction tests
    test_cross_domain_sign_interaction();
    test_cross_positiveint_then_zero_throws();
    test_cross_positiveint_then_nonpositive_throws();
    test_default_domain_complex();

    return TEST_REPORT();
}
