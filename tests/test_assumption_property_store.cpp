
#include "test_common.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <utility>

using namespace LMCAS;


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


void test_domain_roundtrip_all() {
    TEST_CASE("Domain declaration round-trip for all domains and symbols");
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

void test_idempotence() {
    TEST_CASE("Re-declaring same domain is idempotent");
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


void test_hierarchy_implication() {
    TEST_CASE("Declaring domain D implies all ancestor domains return True");
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

void test_non_ancestors_false() {
    TEST_CASE("Domains more specific than declared return False");
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


void test_specificity_preservation() {
    TEST_CASE("Declaring less-specific domain after more-specific is no-op");
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

void test_more_specific_upgrades() {
    TEST_CASE("Declaring more-specific domain upgrades");
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


void test_natural_negative_contradiction() {
    TEST_CASE("Natural + Negative sign returns failure");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    auto failure_184 = store.declare_domain("x", Domain::Natural);
    EXPECT_TRUE(!failure_184.has_value(), "Natural domain with Negative sign returns failure");
    // State unchanged: domain should still be Complex (default)
    EXPECT_TRUE(store.get_domain("x") == Domain::Complex,
        "x domain unchanged after failed Natural declaration");
}

void test_positiveint_negative_contradiction() {
    TEST_CASE("PositiveInt + Negative sign returns failure");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);

    auto failure_201 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(!failure_201.has_value(), "PositiveInt domain with Negative sign returns failure");
}

void test_positiveint_zero_contradiction() {
    TEST_CASE("PositiveInt + Zero sign returns failure");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);

    auto failure_215 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(!failure_215.has_value(), "PositiveInt domain with Zero sign returns failure");
}

void test_positiveint_nonpositive_contradiction() {
    TEST_CASE("PositiveInt + NonPositive sign returns failure");
    PropertyStore store;
    store.declare_sign("x", Sign::NonPositive);

    auto failure_229 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(!failure_229.has_value(), "PositiveInt domain with NonPositive sign returns failure");
}

void test_natural_implied_negative_contradiction() {
    TEST_CASE("Natural + implied Negative (from Negative declaration) returns failure");
    PropertyStore store;
    // Declaring Negative implies NonPositive and NonZero
    store.declare_sign("y", Sign::Negative);

    auto failure_244 = store.declare_domain("y", Domain::Natural);
    EXPECT_TRUE(!failure_244.has_value(), "Natural domain with implied Negative sign returns failure");
}

void test_compatible_domain_sign_no_throw() {
    TEST_CASE("Compatible domain-sign pairs do not throw");
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


void test_positive_implications() {
    TEST_CASE("Positive implies NonNegative and NonZero");
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

void test_negative_implications() {
    TEST_CASE("Negative implies NonPositive and NonZero");
    PropertyStore store;
    store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(store.has_sign("x", Sign::Negative), "x has Negative");
    EXPECT_TRUE(store.has_sign("x", Sign::NonPositive), "x has NonPositive (implied)");
    EXPECT_TRUE(store.has_sign("x", Sign::NonZero), "x has NonZero (implied)");
}

void test_zero_implications() {
    TEST_CASE("Zero implies NonNegative, NonPositive, and Integer domain");
    PropertyStore store;
    store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(store.has_sign("x", Sign::Zero), "x has Zero");
    EXPECT_TRUE(store.has_sign("x", Sign::NonNegative), "x has NonNegative (implied)");
    EXPECT_TRUE(store.has_sign("x", Sign::NonPositive), "x has NonPositive (implied)");
    EXPECT_TRUE(store.has_domain("x", Domain::Integer),
        "x has Integer domain (implied by Zero)");
}

void test_nonneg_nonpos_nonzero_no_extra() {
    TEST_CASE("NonNegative/NonPositive/NonZero have no further implications");
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

void test_redeclaration_noop() {
    TEST_CASE("Re-declaring same sign is no-op");
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

void test_redeclare_implied_sign_noop() {
    TEST_CASE("Re-declaring an implied sign is no-op");
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    auto signs_before = store.get_signs("x");
    // NonNegative is implied by Positive; re-declaring should be no-op
    store.declare_sign("x", Sign::NonNegative);
    auto signs_after = store.get_signs("x");
    EXPECT_TRUE(signs_before == signs_after,
        "Re-declaring implied NonNegative after Positive is no-op");
}


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

void test_all_contradiction_pairs() {
    TEST_CASE("All contradiction pairs throw std::invalid_argument");
    for (const auto& [s1, s2] : CONTRADICTION_PAIRS) {
        PropertyStore store;
        store.declare_sign("x", s1);

        auto failure_392 = store.declare_sign("x", s2);
        EXPECT_TRUE(!failure_392.has_value(), sign_name(s1) + " + " + sign_name(s2) + " returns failure");
        // State unchanged: s2 should not be present
        EXPECT_FALSE(store.has_sign("x", s2),
            "x does not have " + sign_name(s2) + " after failed declaration");
    }
}

void test_implied_contradiction() {
    TEST_CASE("Contradiction via implied signs");
    // Positive implies NonNegative; then declaring Negative should fail
    // because Negative contradicts NonNegative
    PropertyStore store;
    store.declare_sign("x", Sign::Positive);
    auto failure_413 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_413.has_value(), "Negative contradicts implied NonNegative from Positive");
    EXPECT_TRUE(store.has_sign("x", Sign::Positive),
        "x still has Positive after failed Negative declaration");
}

void test_compatible_signs_no_throw() {
    TEST_CASE("Compatible sign pairs do not throw");
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


void test_even_promotes_to_integer() {
    TEST_CASE("Even parity auto-promotes to Integer domain");
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

void test_odd_promotes_to_integer() {
    TEST_CASE("Odd parity auto-promotes to Integer domain");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);
    EXPECT_TRUE(store.get_parity("x") == Parity::Odd, "x has Odd parity");
    EXPECT_TRUE(store.has_domain("x", Domain::Integer),
        "x promoted to Integer by Odd parity");
}

void test_parity_does_not_demote() {
    TEST_CASE("Parity does not demote more-specific domain");
    // Natural is more specific than Integer
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);
    store.declare_parity("x", Parity::Even);
    EXPECT_TRUE(store.get_domain("x") == Domain::Natural,
        "x domain remains Natural (more specific than Integer)");
}

void test_parity_idempotent() {
    TEST_CASE("Re-declaring same parity is no-op");
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

void test_parity_contradiction() {
    TEST_CASE("Contradictory parity returns failure");
    // Even then Odd
    {
        PropertyStore store;
        store.declare_parity("x", Parity::Even);
        auto failure_498 = store.declare_parity("x", Parity::Odd);
        EXPECT_TRUE(!failure_498.has_value(), "Even then Odd returns failure");
        EXPECT_TRUE(store.get_parity("x") == Parity::Even,
            "x parity unchanged after failed Odd declaration");
    }
    // Odd then Even
    {
        PropertyStore store;
        store.declare_parity("y", Parity::Odd);
        auto failure_512 = store.declare_parity("y", Parity::Even);
        EXPECT_TRUE(!failure_512.has_value(), "Odd then Even returns failure");
        EXPECT_TRUE(store.get_parity("y") == Parity::Odd,
            "y parity unchanged after failed Even declaration");
    }
}

void test_default_parity_unknown() {
    TEST_CASE("Default parity is Unknown");
    PropertyStore store;
    EXPECT_TRUE(store.get_parity("undeclared") == Parity::Unknown,
        "Undeclared symbol has Unknown parity");
}


void test_bounded_stored() {
    TEST_CASE("Declaring Bounded stores Bounded");
    for (const auto& sym : TEST_SYMBOLS) {
        PropertyStore store;
        store.declare_bounded(sym, Boundedness::Bounded);
        EXPECT_TRUE(store.get_boundedness(sym) == Boundedness::Bounded,
            sym + " has Bounded");
    }
}

void test_unbounded_stored() {
    TEST_CASE("Declaring Unbounded stores Unbounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
        "x has Unbounded");
}

void test_bounded_with_interval() {
    TEST_CASE("Bounded with interval stores bounds");
    PropertyStore store;

    auto lower_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    store.declare_bounded("x", Boundedness::Bounded, bounds);
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "x has Bounded");
    EXPECT_TRUE(store.get_bounds("x").has_value(),
        "x has bounds stored");
}

void test_idempotent() {
    TEST_CASE("Re-declaring same boundedness is no-op");
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

void test_contradiction_bounded_unbounded() {
    TEST_CASE("Bounded then Unbounded returns failure");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    auto failure_592 = store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(!failure_592.has_value(), "Bounded then Unbounded returns failure");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
        "x remains Bounded after failed Unbounded declaration");
}

void test_contradiction_unbounded_bounded() {
    TEST_CASE("Unbounded then Bounded returns failure");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);

    auto failure_608 = store.declare_bounded("x", Boundedness::Bounded);
    EXPECT_TRUE(!failure_608.has_value(), "Unbounded then Bounded returns failure");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
        "x remains Unbounded after failed Bounded declaration");
}

void test_default_unknown() {
    TEST_CASE("Default boundedness is Unknown");
    PropertyStore store;
    EXPECT_TRUE(store.get_boundedness("undeclared") == Boundedness::Unknown,
        "Undeclared symbol has Unknown boundedness");
    EXPECT_FALSE(store.get_bounds("undeclared").has_value(),
        "Undeclared symbol has no bounds");
}


void test_cross_domain_sign_interaction() {
    TEST_CASE("Cross: Domain Natural then sign Negative throws (reverse direction)");
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);

    auto failure_634 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_634.has_value(), "Negative sign with Natural domain returns failure");
}

void test_cross_positiveint_then_zero_throws() {
    TEST_CASE("Cross: Domain PositiveInt then sign Zero returns failure");
    PropertyStore store;
    store.declare_domain("x", Domain::PositiveInt);

    auto failure_648 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE(!failure_648.has_value(), "Zero sign with PositiveInt domain returns failure");
}

void test_cross_positiveint_then_nonpositive_throws() {
    TEST_CASE("Cross: Domain PositiveInt then sign NonPositive returns failure");
    PropertyStore store;
    store.declare_domain("x", Domain::PositiveInt);

    auto failure_662 = store.declare_sign("x", Sign::NonPositive);
    EXPECT_TRUE(!failure_662.has_value(), "NonPositive sign with PositiveInt domain returns failure");
}

void test_default_domain_complex() {
    TEST_CASE("Default domain is Complex for undeclared symbol");
    PropertyStore store;
    EXPECT_TRUE(store.get_domain("undeclared") == Domain::Complex,
        "Undeclared symbol has Complex domain");
    EXPECT_TRUE(store.has_domain("undeclared", Domain::Complex),
        "has_domain returns true for Complex on undeclared symbol");
}


int main() {
    test_domain_roundtrip_all();
    test_idempotence();

    test_hierarchy_implication();
    test_non_ancestors_false();

    test_specificity_preservation();
    test_more_specific_upgrades();

    test_natural_negative_contradiction();
    test_positiveint_negative_contradiction();
    test_positiveint_zero_contradiction();
    test_positiveint_nonpositive_contradiction();
    test_natural_implied_negative_contradiction();
    test_compatible_domain_sign_no_throw();

    test_positive_implications();
    test_negative_implications();
    test_zero_implications();
    test_nonneg_nonpos_nonzero_no_extra();
    test_redeclaration_noop();
    test_redeclare_implied_sign_noop();

    test_all_contradiction_pairs();
    test_implied_contradiction();
    test_compatible_signs_no_throw();

    test_even_promotes_to_integer();
    test_odd_promotes_to_integer();
    test_parity_does_not_demote();
    test_parity_idempotent();
    test_parity_contradiction();
    test_default_parity_unknown();

    test_bounded_stored();
    test_unbounded_stored();
    test_bounded_with_interval();
    test_idempotent();
    test_contradiction_bounded_unbounded();
    test_contradiction_unbounded_bounded();
    test_default_unknown();

    // Cross-property interaction tests
    test_cross_domain_sign_interaction();
    test_cross_positiveint_then_zero_throws();
    test_cross_positiveint_then_nonpositive_throws();
    test_default_domain_complex();

    return TEST_REPORT();
}
