
#include "test_common.hpp"
#include "property_store.hpp"

using namespace LMCAS;

template <typename F>
static bool returns_invalid_argument(F &&f) {
    auto result = f();
    return !result.has_value() && result.error().code == CasErrc::InvalidArgument;
}

TEST(AssumptionDefinitenessRules, NegativeDefinitePlusPositiveDefiniteThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("C", Definiteness::NegativeDefinite).has_value())) << "definiteness declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_definiteness("C", Definiteness::PositiveDefinite);
    }))) << "NegativeDefinite + PositiveDefinite returns failure";
}

TEST(AssumptionDefinitenessRules, PositiveSemidefiniteUpgradeableToPositiveDefinite) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("G", Definiteness::PositiveSemiDefinite).has_value())) << "definiteness declaration succeeds";
    EXPECT_TRUE((store.declare_definiteness("G", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    EXPECT_TRUE((store.get_definiteness("G") == Definiteness::PositiveDefinite)) << "G upgraded from PositiveSemiDefinite to PositiveDefinite";
}

TEST(AssumptionDefinitenessRules, PositiveDefinitePlusNegativeSemidefiniteThrows) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("H", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    EXPECT_TRUE((returns_invalid_argument([&]() {
        return store.declare_definiteness("H", Definiteness::NegativeSemiDefinite);
    }))) << "PositiveDefinite + NegativeSemiDefinite returns failure";
}

TEST(AssumptionDefinitenessRules, DefinitenessIdempotent) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_definiteness("I", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";
    EXPECT_TRUE((store.declare_definiteness("I", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    EXPECT_TRUE((store.get_definiteness("I") == Definiteness::PositiveDefinite)) << "PositiveDefinite declared twice remains PositiveDefinite";
}
