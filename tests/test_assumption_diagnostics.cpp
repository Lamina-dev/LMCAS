
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

using namespace LMCAS;

/// Generate a random variable name for property tests
static std::string random_var_name() {
    static const std::vector<std::string> prefixes = {"x", "y", "z", "alpha", "beta", "gamma"};
    std::string prefix = *rc::gen::elementOf(prefixes);
    return prefix + "_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
}

/// Convert Domain enum to its expected string representation in error messages

/// Convert Sign enum to its expected string representation in error messages

/// Generate a pair of contradicting domains.
/// Returns (first_domain, second_domain) where declaring both should throw.
/// Transcendental + any sub-Real domain is a contradiction.
struct DomainContradiction {
    Domain first;
    Domain second;
    bool use_transcendental; // If true, first declare transcendental, then second domain
};

static DomainContradiction random_domain_contradiction() {
    // Transcendental contradicts Algebraic, Rational, Integer, Natural, PositiveInt
    static const std::vector<Domain> sub_real_domains = {
        Domain::Algebraic, Domain::Rational, Domain::Integer,
        Domain::Natural, Domain::PositiveInt};

    DomainContradiction result;
    result.use_transcendental = true;
    result.first = Domain::Real; // transcendental implies Real
    result.second = *rc::gen::elementOf(sub_real_domains);
    return result;
}

/// Generate a pair of contradicting signs.
struct SignContradiction {
    Sign first;
    Sign second;
};

static SignContradiction random_sign_contradiction() {
    // Known contradicting sign pairs:
    // Positive + Negative, Positive + Zero, Positive + NonPositive
    // Negative + Zero, Negative + NonNegative
    // Zero + NonZero
    static const std::vector<SignContradiction> contradictions = {
        {Sign::Positive, Sign::Negative},
        {Sign::Positive, Sign::Zero},
        {Sign::Positive, Sign::NonPositive},
        {Sign::Negative, Sign::Zero},
        {Sign::Negative, Sign::NonNegative},
        {Sign::Zero, Sign::NonZero},
    };
    return *rc::gen::elementOf(contradictions);
}

/// Generate a cross-constraint conflict (domain + sign that contradict).
struct CrossConstraint {
    Domain domain;
    Sign sign;
};

static CrossConstraint random_cross_constraint() {
    // Natural domain (non-negative integers) + Negative sign is a contradiction
    // PositiveInt domain (positive integers) + Negative/Zero/NonPositive is a contradiction
    static const std::vector<CrossConstraint> conflicts = {
        {Domain::Natural, Sign::Negative},
        {Domain::PositiveInt, Sign::Negative},
        {Domain::PositiveInt, Sign::Zero},
        {Domain::PositiveInt, Sign::NonPositive},
    };
    return *rc::gen::elementOf(conflicts);
}

TEST(LmcasAssumptionDiagnostics, DomainContradictionError) {
    EXPECT_TRUE(rc::check("For any domain contradiction (Transcendental then sub-Real), "
                          "the Result reports InvalidArgument",
                          []() {
                              std::string var_name = random_var_name();
                              DomainContradiction contradiction = random_domain_contradiction();

                              PropertyStore store;

                              // First declare transcendental (sets Real domain + transcendental flag)
                              EXPECT_TRUE(store.declare_transcendental(var_name).has_value()) << "transcendental declaration succeeds";

                              // Now try to declare a contradicting sub-Real domain
                              auto failure_126 = store.declare_domain(var_name, contradiction.second);
                              RC_ASSERT(!failure_126.has_value());
                              RC_ASSERT(failure_126.error().code == CasErrc::InvalidArgument);
                          }));
}

TEST(LmcasAssumptionDiagnostics, SignContradictionError) {
    EXPECT_TRUE(rc::check("For any sign contradiction (e.g., Positive then Negative), "
                          "the Result reports InvalidArgument",
                          []() {
                              std::string var_name = random_var_name();
                              SignContradiction contradiction = random_sign_contradiction();

                              PropertyStore store;

                              // First declare the initial sign
                              EXPECT_TRUE(store.declare_sign(var_name, contradiction.first).has_value()) << "sign declaration succeeds";

                              // Now try to declare the contradicting sign
                              auto failure_163 = store.declare_sign(var_name, contradiction.second);
                              RC_ASSERT(!failure_163.has_value());
                              RC_ASSERT(failure_163.error().code == CasErrc::InvalidArgument);
                          }));
}

TEST(LmcasAssumptionDiagnostics, CrossConstraintError) {
    EXPECT_TRUE(rc::check("For any cross-constraint conflict (e.g., Natural + Negative), "
                          "the Result reports InvalidArgument",
                          []() {
                              std::string var_name = random_var_name();
                              CrossConstraint conflict = random_cross_constraint();

                              PropertyStore store;

                              // First declare the domain
                              EXPECT_TRUE(store.declare_domain(var_name, conflict.domain).has_value()) << "domain declaration succeeds";

                              // Now try to declare the contradicting sign
                              auto failure_199 = store.declare_sign(var_name, conflict.sign);
                              RC_ASSERT(!failure_199.has_value());
                              RC_ASSERT(failure_199.error().code == CasErrc::InvalidArgument);
                          }));
}

TEST(LmcasAssumptionDiagnostics, FinitenessContradictionError) {
    EXPECT_TRUE(rc::check("For any Finite+Divergent contradiction, "
                          "the Result reports InvalidArgument",
                          []() {
                              std::string var_name = random_var_name();

                              PropertyStore store;

                              // Randomly choose order: Finite then Divergent, or Divergent then Finite
                              bool finite_first = *rc::gen::arbitrary<bool>();

                              EXPECT_TRUE(store.declare_finiteness(var_name, finite_first ? Finiteness::Finite : Finiteness::Divergent).has_value()) << "finiteness declaration succeeds";

                              auto failure_234 = store.declare_finiteness(var_name, finite_first ? Finiteness::Divergent : Finiteness::Finite);
                              RC_ASSERT(!failure_234.has_value());
                              RC_ASSERT(failure_234.error().code == CasErrc::InvalidArgument);
                          }));
}

TEST(LmcasAssumptionDiagnostics, DefinitenessContradictionError) {
    EXPECT_TRUE(rc::check("For any PositiveDefinite+NegativeDefinite contradiction, "
                          "the Result reports InvalidArgument",
                          []() {
                              std::string var_name = random_var_name();

                              PropertyStore store;

                              // Randomly choose order
                              bool pos_first = *rc::gen::arbitrary<bool>();

                              EXPECT_TRUE(store.declare_definiteness(var_name,
                                                                     pos_first ? Definiteness::PositiveDefinite : Definiteness::NegativeDefinite)
                                              .has_value())
                                  << "definiteness declaration succeeds";

                              auto failure_271 = store.declare_definiteness(var_name,
                                                                            pos_first ? Definiteness::NegativeDefinite : Definiteness::PositiveDefinite);
                              RC_ASSERT(!failure_271.has_value());
                              RC_ASSERT(failure_271.error().code == CasErrc::InvalidArgument);
                          }));
}
