
#include "test_common.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include <string>
#include <vector>

using namespace LMCAS;

TEST(AssumptionEnums, DomainAlgebraicPropertyStore) {
    PropertyStore store;

    // Declaring Algebraic should work
    EXPECT_TRUE((store.declare_domain("x", Domain::Algebraic).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("x") == Domain::Algebraic)) << "get_domain returns Algebraic after declaration";

    // Algebraic implies Real (ancestor)
    EXPECT_TRUE((store.has_domain("x", Domain::Real))) << "Algebraic symbol has Real domain";

    // Algebraic implies Complex (ancestor)
    EXPECT_TRUE((store.has_domain("x", Domain::Complex))) << "Algebraic symbol has Complex domain";

    /// Algebraic 与 Rational/Integer 的细化方向相反,因此两项查询均为 false.
    EXPECT_FALSE((store.has_domain("x", Domain::Rational))) << "Algebraic symbol does NOT have Rational domain";

    EXPECT_FALSE((store.has_domain("x", Domain::Integer))) << "Algebraic symbol does NOT have Integer domain";
}

TEST(AssumptionEnums, DomainAlgebraicUpgrade) {
    PropertyStore store;

    // Start with Algebraic
    EXPECT_TRUE((store.declare_domain("y", Domain::Algebraic).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("y") == Domain::Algebraic)) << "y starts as Algebraic";

    // Upgrade to Rational (more specific)
    EXPECT_TRUE((store.declare_domain("y", Domain::Rational).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("y") == Domain::Rational)) << "y upgraded to Rational";

    // Rational still implies Algebraic
    EXPECT_TRUE((store.has_domain("y", Domain::Algebraic))) << "Rational symbol still has Algebraic domain";
}

TEST(AssumptionEnums, DomainAlgebraicNoDowngrade) {
    PropertyStore store;

    // Start with Rational
    EXPECT_TRUE((store.declare_domain("z", Domain::Rational).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("z") == Domain::Rational)) << "z starts as Rational";

    // Declaring Algebraic (less specific) should be a no-op
    EXPECT_TRUE((store.declare_domain("z", Domain::Algebraic).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("z") == Domain::Rational)) << "z remains Rational after Algebraic declaration (no downgrade)";

    // Declaring Real (even less specific) should also be a no-op
    EXPECT_TRUE((store.declare_domain("z", Domain::Real).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("z") == Domain::Rational)) << "z remains Rational after Real declaration (no downgrade)";
}
