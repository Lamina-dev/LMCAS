
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
    Domain::Integer, Domain::Natural, Domain::PositiveInt};

static const std::vector<std::string> TEST_SYMBOLS = {
    "x", "y", "alpha", "longVariableName123", "a_b_c"};

static int domain_specificity(Domain d) {
    switch (d) {
    case Domain::Complex: {
        return 0;
    }
    case Domain::Real: {
        return 1;
    }
    case Domain::Algebraic: {
        return 2;
    }
    case Domain::Rational: {
        return 3;
    }
    case Domain::Integer: {
        return 4;
    }
    case Domain::Natural: {
        return 5;
    }
    case Domain::PositiveInt: {
        return 6;
    }
    }
    return 0;
}

static std::string domain_name(Domain d) {
    switch (d) {
    case Domain::Complex: {
        return "Complex";
    }
    case Domain::Real: {
        return "Real";
    }
    case Domain::Algebraic: {
        return "Algebraic";
    }
    case Domain::Rational: {
        return "Rational";
    }
    case Domain::Integer: {
        return "Integer";
    }
    case Domain::Natural: {
        return "Natural";
    }
    case Domain::PositiveInt: {
        return "PositiveInt";
    }
    }
    return "?";
}

TEST(AssumptionDomainStore, DomainRoundtripAll) {
    for (const auto &sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            EXPECT_TRUE((store.declare_domain(sym, d).has_value())) << "domain declaration succeeds";
            Domain result = store.get_domain(sym);
            bool at_least_as_specific =
                domain_specificity(result) >= domain_specificity(d);
            EXPECT_TRUE((at_least_as_specific)) << sym + " domain after declaring " + domain_name(d) +
                                                       " is at least as specific";
        }
    }
}

TEST(AssumptionDomainStore, Idempotence) {
    for (const auto &sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            EXPECT_TRUE((store.declare_domain(sym, d).has_value())) << "domain declaration succeeds";
            Domain first = store.get_domain(sym);
            EXPECT_TRUE((store.declare_domain(sym, d).has_value())) << "domain declaration succeeds";
            Domain second = store.get_domain(sym);
            EXPECT_TRUE((first == second)) << sym + " domain unchanged after re-declaring " + domain_name(d);
        }
    }
}

TEST(AssumptionDomainStore, HierarchyImplication) {
    for (const auto &sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            PropertyStore store;
            EXPECT_TRUE((store.declare_domain(sym, d).has_value())) << "domain declaration succeeds";

            for (Domain ancestor : ALL_DOMAINS) {
                if (domain_specificity(ancestor) <= domain_specificity(d)) {
                    EXPECT_TRUE((store.has_domain(sym, ancestor))) << sym + " with " + domain_name(d) +
                                                                          " has ancestor " + domain_name(ancestor);
                }
            }
        }
    }
}

TEST(AssumptionDomainStore, NonAncestorsFalse) {
    for (Domain d : ALL_DOMAINS) {
        PropertyStore store;
        EXPECT_TRUE((store.declare_domain("x", d).has_value())) << "domain declaration succeeds";

        for (Domain more_specific : ALL_DOMAINS) {
            if (domain_specificity(more_specific) > domain_specificity(d)) {
                EXPECT_FALSE((store.has_domain("x", more_specific))) << "x with " + domain_name(d) +
                                                                            " does NOT have " + domain_name(more_specific);
            }
        }
    }
}

TEST(AssumptionDomainStore, SpecificityPreservation) {
    for (const auto &sym : TEST_SYMBOLS) {
        for (Domain specific : ALL_DOMAINS) {
            for (Domain ancestor : ALL_DOMAINS) {
                if (domain_specificity(ancestor) < domain_specificity(specific)) {
                    PropertyStore store;
                    EXPECT_TRUE((store.declare_domain(sym, specific).has_value())) << "domain declaration succeeds";
                    Domain before = store.get_domain(sym);

                    EXPECT_TRUE((store.declare_domain(sym, ancestor).has_value())) << "domain declaration succeeds";
                    Domain after = store.get_domain(sym);

                    EXPECT_TRUE((before == after)) << sym + ": declaring " + domain_name(ancestor) +
                                                          " after " + domain_name(specific) + " is no-op";
                }
            }
        }
    }
}

TEST(AssumptionDomainStore, MoreSpecificUpgrades) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Real).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("x") == Domain::Real)) << "x starts at Real";

    EXPECT_TRUE((store.declare_domain("x", Domain::Integer).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("x") == Domain::Integer)) << "x upgraded to Integer";

    EXPECT_TRUE((store.declare_domain("x", Domain::PositiveInt).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.get_domain("x") == Domain::PositiveInt)) << "x upgraded to PositiveInt";
}
