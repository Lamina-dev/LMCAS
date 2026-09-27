
#include "test_common.hpp"
#include "property_store.hpp"
#include <vector>
#include <string>

using namespace LMCAS;

static const std::vector<std::string> TEST_SYMBOLS = {
    "x", "y", "alpha", "longVariableName123", "a_b_c"};

TEST(AssumptionPropertyStore, NaturalNegativeContradiction) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Negative).has_value())) << "sign declaration succeeds";

    auto failure_184 = store.declare_domain("x", Domain::Natural);
    EXPECT_TRUE((!failure_184.has_value())) << "Natural domain with Negative sign returns failure";
    // State unchanged: domain should still be Complex (default)
    EXPECT_TRUE((store.get_domain("x") == Domain::Complex)) << "x domain unchanged after failed Natural declaration";
}

TEST(AssumptionPropertyStore, PositiveintNegativeContradiction) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Negative).has_value())) << "sign declaration succeeds";

    auto failure_201 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE((!failure_201.has_value())) << "PositiveInt domain with Negative sign returns failure";
}

TEST(AssumptionPropertyStore, PositiveintZeroContradiction) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Zero).has_value())) << "sign declaration succeeds";

    auto failure_215 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE((!failure_215.has_value())) << "PositiveInt domain with Zero sign returns failure";
}

TEST(AssumptionPropertyStore, PositiveintNonpositiveContradiction) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonPositive).has_value())) << "sign declaration succeeds";

    auto failure_229 = store.declare_domain("x", Domain::PositiveInt);
    EXPECT_TRUE((!failure_229.has_value())) << "PositiveInt domain with NonPositive sign returns failure";
}

TEST(AssumptionPropertyStore, NaturalImpliedNegativeContradiction) {
    PropertyStore store;
    // Declaring Negative implies NonPositive and NonZero
    EXPECT_TRUE((store.declare_sign("y", Sign::Negative).has_value())) << "sign declaration succeeds";

    auto failure_244 = store.declare_domain("y", Domain::Natural);
    EXPECT_TRUE((!failure_244.has_value())) << "Natural domain with implied Negative sign returns failure";
}

TEST(AssumptionPropertyStore, CompatibleDomainSignNoThrow) {
    // Natural + Positive is fine
    {
        PropertyStore store;
        EXPECT_TRUE((store.declare_sign("a", Sign::Positive).has_value())) << "sign declaration succeeds";
        EXPECT_TRUE((store.declare_domain("a", Domain::Natural).has_value())) << "domain declaration succeeds";
        EXPECT_TRUE((store.get_domain("a") == Domain::Natural)) << "Natural + Positive is compatible";
    }
    // Natural + NonNegative is fine
    {
        PropertyStore store;
        EXPECT_TRUE((store.declare_sign("b", Sign::NonNegative).has_value())) << "sign declaration succeeds";
        EXPECT_TRUE((store.declare_domain("b", Domain::Natural).has_value())) << "domain declaration succeeds";
        EXPECT_TRUE((store.get_domain("b") == Domain::Natural)) << "Natural + NonNegative is compatible";
    }
    // PositiveInt + Positive is fine
    {
        PropertyStore store;
        EXPECT_TRUE((store.declare_sign("c", Sign::Positive).has_value())) << "sign declaration succeeds";
        EXPECT_TRUE((store.declare_domain("c", Domain::PositiveInt).has_value())) << "domain declaration succeeds";
        EXPECT_TRUE((store.get_domain("c") == Domain::PositiveInt)) << "PositiveInt + Positive is compatible";
    }
    // Integer + Negative is fine
    {
        PropertyStore store;
        EXPECT_TRUE((store.declare_sign("d", Sign::Negative).has_value())) << "sign declaration succeeds";
        EXPECT_TRUE((store.declare_domain("d", Domain::Integer).has_value())) << "domain declaration succeeds";
        EXPECT_TRUE((store.get_domain("d") == Domain::Integer)) << "Integer + Negative is compatible";
    }
}

TEST(AssumptionPropertyStore, EvenPromotesToInteger) {
    for (const auto &sym : TEST_SYMBOLS) {
        PropertyStore store;
        // Default domain is Complex
        EXPECT_TRUE((store.declare_parity(sym, Parity::Even).has_value())) << "parity declaration succeeds";
        EXPECT_TRUE((store.get_parity(sym) == Parity::Even)) << sym + " has Even parity";
        EXPECT_TRUE((store.has_domain(sym, Domain::Integer))) << sym + " promoted to Integer by Even parity";
    }
}

TEST(AssumptionPropertyStore, OddPromotesToInteger) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((store.get_parity("x") == Parity::Odd)) << "x has Odd parity";
    EXPECT_TRUE((store.has_domain("x", Domain::Integer))) << "x promoted to Integer by Odd parity";
}

TEST(AssumptionPropertyStore, ParityDoesNotDemote) {
    // Natural is more specific than Integer
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Natural).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((store.get_domain("x") == Domain::Natural)) << "x domain remains Natural (more specific than Integer)";
}

TEST(AssumptionPropertyStore, CrossDomainSignInteraction) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Natural).has_value())) << "domain declaration succeeds";

    auto failure_634 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_634.has_value())) << "Negative sign with Natural domain returns failure";
}

TEST(AssumptionPropertyStore, CrossPositiveintThenZeroThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::PositiveInt).has_value())) << "domain declaration succeeds";

    auto failure_648 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE((!failure_648.has_value())) << "Zero sign with PositiveInt domain returns failure";
}

TEST(AssumptionPropertyStore, CrossPositiveintThenNonpositiveThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::PositiveInt).has_value())) << "domain declaration succeeds";

    auto failure_662 = store.declare_sign("x", Sign::NonPositive);
    EXPECT_TRUE((!failure_662.has_value())) << "NonPositive sign with PositiveInt domain returns failure";
}
