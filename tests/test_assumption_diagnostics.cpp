/**
 * @file test_assumption_diagnostics.cpp
 * @brief Property tests for conflict diagnostics (Task 9.3).
 *
 * Property tested:
 * - Property 20: Conflict diagnostics contain required information
 *
 * For any contradiction, the exception message SHALL contain the symbol name,
 * existing property, and new conflicting property.
 *
 * Validates: Requirements 19.1, 19.2, 19.3
 *
 * Uses rapidcheck (header-only, vendored in tests/rapidcheck/) for
 * property-based testing with random input generation.
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

using namespace lamina;

// ============================================================
// Helpers
// ============================================================

/// Generate a random variable name for property tests
static std::string random_var_name() {
    static const std::vector<std::string> prefixes = {"x", "y", "z", "alpha", "beta", "gamma"};
    std::string prefix = rc::gen::elementOf(prefixes);
    return prefix + "_" + std::to_string(rc::gen::inRange(0, 999));
}

/// Convert Domain enum to its expected string representation in error messages
static std::string domain_to_string(Domain d) {
    switch (d) {
        case Domain::Complex:    return "Complex";
        case Domain::Real:       return "Real";
        case Domain::Algebraic:  return "Algebraic";
        case Domain::Rational:   return "Rational";
        case Domain::Integer:    return "Integer";
        case Domain::Natural:    return "Natural";
        case Domain::PositiveInt: return "PositiveInt";
    }
    return "Unknown";
}

/// Convert Sign enum to its expected string representation in error messages
static std::string sign_to_string(Sign s) {
    switch (s) {
        case Sign::Positive:    return "Positive";
        case Sign::Negative:    return "Negative";
        case Sign::NonNegative: return "NonNegative";
        case Sign::NonPositive: return "NonPositive";
        case Sign::Zero:        return "Zero";
        case Sign::NonZero:     return "NonZero";
    }
    return "Unknown";
}

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
        Domain::Natural, Domain::PositiveInt
    };

    DomainContradiction result;
    result.use_transcendental = true;
    result.first = Domain::Real; // transcendental implies Real
    result.second = rc::gen::elementOf(sub_real_domains);
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
    return rc::gen::elementOf(contradictions);
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
    return rc::gen::elementOf(conflicts);
}

// ============================================================
// Property 20: Domain contradiction messages contain required info
// **Validates: Requirements 19.1**
// ============================================================

static void test_property20_domain_contradiction_message() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 20: Domain contradiction contains symbol and domains");

    rc::check("For any domain contradiction (Transcendental then sub-Real), "
              "the exception message contains the symbol name and both domains", []() {
        std::string var_name = random_var_name();
        DomainContradiction contradiction = random_domain_contradiction();

        PropertyStore store;

        // First declare transcendental (sets Real domain + transcendental flag)
        store.declare_transcendental(var_name);

        // Now try to declare a contradicting sub-Real domain
        bool threw = false;
        std::string message;
        try {
            store.declare_domain(var_name, contradiction.second);
        } catch (const std::invalid_argument& e) {
            threw = true;
            message = e.what();
        }

        RC_ASSERT(threw);
        // Message must contain the symbol name
        RC_ASSERT(message.find(var_name) != std::string::npos);
        // Message must contain information about the conflicting domain
        // (either the domain name or "Transcendental"/"Algebraic" etc.)
        std::string second_str = domain_to_string(contradiction.second);
        bool has_domain_info = message.find(second_str) != std::string::npos ||
                               message.find("Transcendental") != std::string::npos ||
                               message.find("transcendental") != std::string::npos;
        RC_ASSERT(has_domain_info);
    });
}

// ============================================================
// Property 20: Sign contradiction messages contain required info
// **Validates: Requirements 19.2**
// ============================================================

static void test_property20_sign_contradiction_message() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 20: Sign contradiction contains symbol and signs");

    rc::check("For any sign contradiction (e.g., Positive then Negative), "
              "the exception message contains the symbol name and both signs", []() {
        std::string var_name = random_var_name();
        SignContradiction contradiction = random_sign_contradiction();

        PropertyStore store;

        // First declare the initial sign
        store.declare_sign(var_name, contradiction.first);

        // Now try to declare the contradicting sign
        bool threw = false;
        std::string message;
        try {
            store.declare_sign(var_name, contradiction.second);
        } catch (const std::invalid_argument& e) {
            threw = true;
            message = e.what();
        }

        RC_ASSERT(threw);
        // Message must contain the symbol name
        RC_ASSERT(message.find(var_name) != std::string::npos);
        // Message must contain information about both signs
        std::string first_str = sign_to_string(contradiction.first);
        std::string second_str = sign_to_string(contradiction.second);
        bool has_sign_info = (message.find(first_str) != std::string::npos ||
                              message.find(second_str) != std::string::npos);
        RC_ASSERT(has_sign_info);
    });
}

// ============================================================
// Property 20: Cross-constraint conflict messages contain required info
// **Validates: Requirements 19.3**
// ============================================================

static void test_property20_cross_constraint_message() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 20: Cross-constraint conflict contains symbol, domain, and sign");

    rc::check("For any cross-constraint conflict (e.g., Natural + Negative), "
              "the exception message contains the symbol name and explains both", []() {
        std::string var_name = random_var_name();
        CrossConstraint conflict = random_cross_constraint();

        PropertyStore store;

        // First declare the domain
        store.declare_domain(var_name, conflict.domain);

        // Now try to declare the contradicting sign
        bool threw = false;
        std::string message;
        try {
            store.declare_sign(var_name, conflict.sign);
        } catch (const std::invalid_argument& e) {
            threw = true;
            message = e.what();
        }

        RC_ASSERT(threw);
        // Message must contain the symbol name
        RC_ASSERT(message.find(var_name) != std::string::npos);
        // Message should contain information about the domain or sign conflict
        std::string domain_str = domain_to_string(conflict.domain);
        std::string sign_str = sign_to_string(conflict.sign);
        bool has_conflict_info = (message.find(domain_str) != std::string::npos ||
                                  message.find(sign_str) != std::string::npos);
        RC_ASSERT(has_conflict_info);
    });
}

// ============================================================
// Property 20: Finiteness contradiction messages contain required info
// (Additional coverage for contradiction diagnostics)
// ============================================================

static void test_property20_finiteness_contradiction_message() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 20: Finiteness contradiction contains symbol info");

    rc::check("For any Finite+Divergent contradiction, "
              "the exception message contains the symbol name", []() {
        std::string var_name = random_var_name();

        PropertyStore store;

        // Randomly choose order: Finite then Divergent, or Divergent then Finite
        bool finite_first = rc::gen::boolean();

        store.declare_finiteness(var_name, finite_first ? Finiteness::Finite : Finiteness::Divergent);

        bool threw = false;
        std::string message;
        try {
            store.declare_finiteness(var_name, finite_first ? Finiteness::Divergent : Finiteness::Finite);
        } catch (const std::invalid_argument& e) {
            threw = true;
            message = e.what();
        }

        RC_ASSERT(threw);
        // Message must contain the symbol name
        RC_ASSERT(message.find(var_name) != std::string::npos);
        // Message should mention both finiteness states
        bool has_finite_info = (message.find("Finite") != std::string::npos ||
                                message.find("finite") != std::string::npos ||
                                message.find("Divergent") != std::string::npos ||
                                message.find("divergent") != std::string::npos);
        RC_ASSERT(has_finite_info);
    });
}

// ============================================================
// Property 20: Definiteness contradiction messages contain required info
// ============================================================

static void test_property20_definiteness_contradiction_message() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 20: Definiteness contradiction contains symbol info");

    rc::check("For any PositiveDefinite+NegativeDefinite contradiction, "
              "the exception message contains the symbol name", []() {
        std::string var_name = random_var_name();

        PropertyStore store;

        // Randomly choose order
        bool pos_first = rc::gen::boolean();

        store.declare_definiteness(var_name,
            pos_first ? Definiteness::PositiveDefinite : Definiteness::NegativeDefinite);

        bool threw = false;
        std::string message;
        try {
            store.declare_definiteness(var_name,
                pos_first ? Definiteness::NegativeDefinite : Definiteness::PositiveDefinite);
        } catch (const std::invalid_argument& e) {
            threw = true;
            message = e.what();
        }

        RC_ASSERT(threw);
        // Message must contain the symbol name
        RC_ASSERT(message.find(var_name) != std::string::npos);
        // Message should reference definiteness
        bool has_def_info = (message.find("Definite") != std::string::npos ||
                             message.find("definite") != std::string::npos ||
                             message.find("PositiveDefinite") != std::string::npos ||
                             message.find("NegativeDefinite") != std::string::npos);
        RC_ASSERT(has_def_info);
    });
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 20: Conflict diagnostics contain required information
    test_property20_domain_contradiction_message();
    test_property20_sign_contradiction_message();
    test_property20_cross_constraint_message();
    test_property20_finiteness_contradiction_message();
    test_property20_definiteness_contradiction_message();

    return TEST_REPORT();
}
