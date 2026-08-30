
#include "test_common.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include <string>
#include <vector>

using namespace lamina;


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
    return -1;
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


static void test_domain_ordering() {
    TEST_CASE("Domain ordering: Complex < Real < Algebraic < Rational < Integer < Natural < PositiveInt");

    // Verify the full specificity chain
    EXPECT_TRUE(domain_specificity(Domain::Complex) < domain_specificity(Domain::Real),
        "Complex < Real");
    EXPECT_TRUE(domain_specificity(Domain::Real) < domain_specificity(Domain::Algebraic),
        "Real < Algebraic");
    EXPECT_TRUE(domain_specificity(Domain::Algebraic) < domain_specificity(Domain::Rational),
        "Algebraic < Rational");
    EXPECT_TRUE(domain_specificity(Domain::Rational) < domain_specificity(Domain::Integer),
        "Rational < Integer");
    EXPECT_TRUE(domain_specificity(Domain::Integer) < domain_specificity(Domain::Natural),
        "Integer < Natural");
    EXPECT_TRUE(domain_specificity(Domain::Natural) < domain_specificity(Domain::PositiveInt),
        "Natural < PositiveInt");
}

static void test_domain_algebraic_position() {
    TEST_CASE("Domain: Algebraic is between Real and Rational");

    // Algebraic is more specific than Real
    EXPECT_TRUE(domain_specificity(Domain::Algebraic) > domain_specificity(Domain::Real),
        "Algebraic more specific than Real");

    // Algebraic is less specific than Rational
    EXPECT_TRUE(domain_specificity(Domain::Algebraic) < domain_specificity(Domain::Rational),
        "Algebraic less specific than Rational");

    // Algebraic is less specific than Integer
    EXPECT_TRUE(domain_specificity(Domain::Algebraic) < domain_specificity(Domain::Integer),
        "Algebraic less specific than Integer");

    // Algebraic is more specific than Complex
    EXPECT_TRUE(domain_specificity(Domain::Algebraic) > domain_specificity(Domain::Complex),
        "Algebraic more specific than Complex");
}

static void test_domain_algebraic_property_store() {
    TEST_CASE("PropertyStore: Algebraic domain declaration and hierarchy");

    PropertyStore store;

    // Declaring Algebraic should work
    store.declare_domain("x", Domain::Algebraic);
    EXPECT_TRUE(store.get_domain("x") == Domain::Algebraic,
        "get_domain returns Algebraic after declaration");

    // Algebraic implies Real (ancestor)
    EXPECT_TRUE(store.has_domain("x", Domain::Real),
        "Algebraic symbol has Real domain");

    // Algebraic implies Complex (ancestor)
    EXPECT_TRUE(store.has_domain("x", Domain::Complex),
        "Algebraic symbol has Complex domain");

    /// Algebraic 与 Rational/Integer 的细化方向相反，因此两项查询均为 false。
    EXPECT_FALSE(store.has_domain("x", Domain::Rational),
        "Algebraic symbol does NOT have Rational domain");

    EXPECT_FALSE(store.has_domain("x", Domain::Integer),
        "Algebraic symbol does NOT have Integer domain");
}

static void test_domain_algebraic_upgrade() {
    TEST_CASE("PropertyStore: Upgrading from Algebraic to more specific domain");

    PropertyStore store;

    // Start with Algebraic
    store.declare_domain("y", Domain::Algebraic);
    EXPECT_TRUE(store.get_domain("y") == Domain::Algebraic,
        "y starts as Algebraic");

    // Upgrade to Rational (more specific)
    store.declare_domain("y", Domain::Rational);
    EXPECT_TRUE(store.get_domain("y") == Domain::Rational,
        "y upgraded to Rational");

    // Rational still implies Algebraic
    EXPECT_TRUE(store.has_domain("y", Domain::Algebraic),
        "Rational symbol still has Algebraic domain");
}

static void test_domain_algebraic_no_downgrade() {
    TEST_CASE("PropertyStore: Declaring less specific domain is a no-op");

    PropertyStore store;

    // Start with Rational
    store.declare_domain("z", Domain::Rational);
    EXPECT_TRUE(store.get_domain("z") == Domain::Rational,
        "z starts as Rational");

    // Declaring Algebraic (less specific) should be a no-op
    store.declare_domain("z", Domain::Algebraic);
    EXPECT_TRUE(store.get_domain("z") == Domain::Rational,
        "z remains Rational after Algebraic declaration (no downgrade)");

    // Declaring Real (even less specific) should also be a no-op
    store.declare_domain("z", Domain::Real);
    EXPECT_TRUE(store.get_domain("z") == Domain::Rational,
        "z remains Rational after Real declaration (no downgrade)");
}


static void test_default_symbol_properties() {
    TEST_CASE("SymbolProperties defaults: new fields have correct initial values");

    PropertyStore store;

    /// 查询首次出现的符号，PropertyStore 为各属性返回默认值。

    // Default domain is Complex
    EXPECT_TRUE(store.get_domain("fresh") == Domain::Complex,
        "Default domain is Complex");

    // Default boundedness is Unknown
    EXPECT_TRUE(store.get_boundedness("fresh") == Boundedness::Unknown,
        "Default boundedness is Unknown");

    // Default parity is Unknown
    EXPECT_TRUE(store.get_parity("fresh") == Parity::Unknown,
        "Default parity is Unknown");

    // Default signs set is empty
    auto signs = store.get_signs("fresh");
    EXPECT_TRUE(signs.empty(),
        "Default signs set is empty");
}

static void test_default_new_fields_via_domain_declaration() {
    TEST_CASE("SymbolProperties defaults: new fields remain default after domain declaration");

    PropertyStore store;

    // Declare a domain to create the symbol entry, then verify new fields
    store.declare_domain("x", Domain::Real);

    // After declaring domain, the new fields should still be at defaults:
    // transcendental = false (symbol is not transcendental by default)
    // We verify this indirectly: a Real symbol is not automatically transcendental
    // The PropertyStore doesn't expose transcendental directly yet (task 2.3),
    // but we can verify the domain is Real and not Algebraic
    EXPECT_TRUE(store.get_domain("x") == Domain::Real,
        "Domain is Real as declared");

    // Boundedness remains Unknown (finiteness=Unknown doesn't change it)
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unknown,
        "Boundedness remains Unknown (finiteness default)");

    // Parity remains Unknown
    EXPECT_TRUE(store.get_parity("x") == Parity::Unknown,
        "Parity remains Unknown");
}


static void test_monotonicity_enum() {
    TEST_CASE("Monotonicity enum: all values accessible");

    // Verify all enum values compile and are distinct
    Monotonicity m1 = Monotonicity::Increasing;
    Monotonicity m2 = Monotonicity::Decreasing;
    Monotonicity m3 = Monotonicity::NonDecreasing;
    Monotonicity m4 = Monotonicity::NonIncreasing;
    Monotonicity m5 = Monotonicity::Unknown;

    EXPECT_TRUE(m1 != m2, "Increasing != Decreasing");
    EXPECT_TRUE(m1 != m3, "Increasing != NonDecreasing");
    EXPECT_TRUE(m1 != m4, "Increasing != NonIncreasing");
    EXPECT_TRUE(m1 != m5, "Increasing != Unknown");
    EXPECT_TRUE(m2 != m3, "Decreasing != NonDecreasing");
    EXPECT_TRUE(m4 != m5, "NonIncreasing != Unknown");
}


static void test_definiteness_enum() {
    TEST_CASE("Definiteness enum: all values accessible");

    Definiteness d1 = Definiteness::PositiveDefinite;
    Definiteness d2 = Definiteness::PositiveSemiDefinite;
    Definiteness d3 = Definiteness::NegativeDefinite;
    Definiteness d4 = Definiteness::NegativeSemiDefinite;
    Definiteness d5 = Definiteness::Indefinite;
    Definiteness d6 = Definiteness::Unknown;

    EXPECT_TRUE(d1 != d2, "PositiveDefinite != PositiveSemiDefinite");
    EXPECT_TRUE(d1 != d3, "PositiveDefinite != NegativeDefinite");
    EXPECT_TRUE(d3 != d4, "NegativeDefinite != NegativeSemiDefinite");
    EXPECT_TRUE(d5 != d6, "Indefinite != Unknown");
    EXPECT_TRUE(d2 != d4, "PositiveSemiDefinite != NegativeSemiDefinite");
}


static void test_finiteness_enum() {
    TEST_CASE("Finiteness enum: all values accessible");

    Finiteness f1 = Finiteness::Finite;
    Finiteness f2 = Finiteness::Divergent;
    Finiteness f3 = Finiteness::Unknown;

    EXPECT_TRUE(f1 != f2, "Finite != Divergent");
    EXPECT_TRUE(f1 != f3, "Finite != Unknown");
    EXPECT_TRUE(f2 != f3, "Divergent != Unknown");
}


static void test_all_domains_distinct() {
    TEST_CASE("Domain enum: all seven values are distinct");

    std::vector<Domain> all = {
        Domain::Complex, Domain::Real, Domain::Algebraic,
        Domain::Rational, Domain::Integer, Domain::Natural, Domain::PositiveInt
    };

    for (size_t i = 0; i < all.size(); ++i) {
        for (size_t j = i + 1; j < all.size(); ++j) {
            EXPECT_TRUE(all[i] != all[j],
                domain_name(all[i]) + " != " + domain_name(all[j]));
        }
    }
}


static void test_specificity_strictly_increasing() {
    TEST_CASE("Domain specificity: values are strictly increasing through hierarchy");

    std::vector<Domain> ordered = {
        Domain::Complex, Domain::Real, Domain::Algebraic,
        Domain::Rational, Domain::Integer, Domain::Natural, Domain::PositiveInt
    };

    for (size_t i = 0; i + 1 < ordered.size(); ++i) {
        int curr = domain_specificity(ordered[i]);
        int next = domain_specificity(ordered[i + 1]);
        EXPECT_TRUE(curr < next,
            domain_name(ordered[i]) + " specificity (" + std::to_string(curr) +
            ") < " + domain_name(ordered[i + 1]) + " specificity (" + std::to_string(next) + ")");
    }
}


int main() {
    test_domain_ordering();
    test_domain_algebraic_position();
    test_domain_algebraic_property_store();
    test_domain_algebraic_upgrade();
    test_domain_algebraic_no_downgrade();
    test_default_symbol_properties();
    test_default_new_fields_via_domain_declaration();
    test_monotonicity_enum();
    test_definiteness_enum();
    test_finiteness_enum();
    test_all_domains_distinct();
    test_specificity_strictly_increasing();

    return TEST_REPORT();
}
