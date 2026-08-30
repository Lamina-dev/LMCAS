
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <vector>
#include <memory>

using namespace lamina;

static AssumptionContext deserialize_success(const std::string& data) {
    auto result = AssumptionContext::deserialize(data);
    RC_ASSERT(result.has_value());
    return result ? std::move(result.value()) : AssumptionContext();
}


/// Generate a random Domain (excluding Complex which is the default).
static Domain random_domain() {
    std::vector<Domain> domains = {
        Domain::Real, Domain::Algebraic, Domain::Rational,
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(domains);
}

/// Generate a random Sign.
static Sign random_sign() {
    std::vector<Sign> signs = {
        Sign::Positive, Sign::Negative, Sign::NonNegative,
        Sign::NonPositive, Sign::NonZero
    };
    return rc::gen::elementOf(signs);
}

/// Generate a random Parity.
static Parity random_parity() {
    std::vector<Parity> parities = {Parity::Even, Parity::Odd};
    return rc::gen::elementOf(parities);
}

/// Generate a random Boundedness.
static Boundedness random_boundedness() {
    std::vector<Boundedness> values = {Boundedness::Bounded, Boundedness::Unbounded};
    return rc::gen::elementOf(values);
}

/// Generate a random Finiteness (non-Unknown).
static Finiteness random_finiteness() {
    std::vector<Finiteness> values = {Finiteness::Finite, Finiteness::Divergent};
    return rc::gen::elementOf(values);
}

/// Generate a random Definiteness (non-Unknown).
static Definiteness random_definiteness() {
    std::vector<Definiteness> values = {
        Definiteness::PositiveDefinite, Definiteness::PositiveSemiDefinite,
        Definiteness::NegativeDefinite, Definiteness::NegativeSemiDefinite,
        Definiteness::Indefinite
    };
    return rc::gen::elementOf(values);
}

/// Generate a unique variable name.
static std::string random_var_name(int prefix_id) {
    return "s" + std::to_string(prefix_id) + "_" + std::to_string(rc::gen::inRange(0, 99));
}


/// Test: Domain round-trip
static void test_serialization_domain_roundtrip() {
    TEST_CASE("Domain serialization round-trip");

    rc::check("Domain declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(0);
        Domain dom = random_domain();

        ctx.assume_domain(var, dom);

        // Serialize and deserialize
        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        // Query results must match
        RC_ASSERT(restored.get_domain(var) == dom);
        RC_ASSERT(restored.has_domain(var, dom));
    });
}

/// Test: Sign round-trip
static void test_serialization_sign_roundtrip() {
    TEST_CASE("Sign serialization round-trip");

    rc::check("Sign declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(1);
        Sign sign = random_sign();

        ctx.assume_sign(var, sign);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.has_sign(var, sign));
    });
}

/// Test: Parity round-trip
static void test_serialization_parity_roundtrip() {
    TEST_CASE("Parity serialization round-trip");

    rc::check("Parity declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(2);
        Parity par = random_parity();

        ctx.current_properties().declare_parity(var, par);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.get_parity(var) == par);
    });
}

/// Test: Boundedness round-trip
static void test_serialization_boundedness_roundtrip() {
    TEST_CASE("Boundedness serialization round-trip");

    rc::check("Boundedness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(3);
        Boundedness bnd = random_boundedness();

        ctx.current_properties().declare_bounded(var, bnd);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.get_boundedness(var) == bnd);
    });
}

/// Test: Transcendental round-trip
static void test_serialization_transcendental_roundtrip() {
    TEST_CASE("Transcendental serialization round-trip");

    rc::check("Transcendental declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(4);

        ctx.current_properties().declare_transcendental(var);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        // Transcendental implies Real domain
        RC_ASSERT(restored.get_domain(var) == Domain::Real);
        RC_ASSERT(restored.current_properties().is_transcendental(var));
    });
}

/// Test: Finiteness round-trip
static void test_serialization_finiteness_roundtrip() {
    TEST_CASE("Finiteness serialization round-trip");

    rc::check("Finiteness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(5);
        Finiteness fin = random_finiteness();

        ctx.current_properties().declare_finiteness(var, fin);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.current_properties().get_finiteness(var) == fin);

        // If Finite, boundedness should also be Bounded
        if (fin == Finiteness::Finite) {
            RC_ASSERT(restored.get_boundedness(var) == Boundedness::Bounded);
        }
    });
}

/// Test: Definiteness round-trip
static void test_serialization_definiteness_roundtrip() {
    TEST_CASE("Definiteness serialization round-trip");

    rc::check("Definiteness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(6);
        Definiteness def = random_definiteness();

        ctx.current_properties().declare_definiteness(var, def);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.current_properties().get_definiteness(var) == def);
    });
}

/// Test: Multi-scope round-trip
static void test_serialization_multi_scope_roundtrip() {
    TEST_CASE("Multi-scope serialization round-trip");

    rc::check("Multi-scope contexts survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;

        // Root scope: declare domain and sign for var_a
        std::string var_a = "a_" + std::to_string(rc::gen::inRange(0, 99));
        Domain dom_a = random_domain();
        ctx.assume_domain(var_a, dom_a);

        // Push scope: declare sign for var_b
        ctx.push();
        std::string var_b = "b_" + std::to_string(rc::gen::inRange(0, 99));
        Sign sign_b = random_sign();
        ctx.assume_sign(var_b, sign_b);

        int depth_before = ctx.depth();

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        // Depth should match
        RC_ASSERT(restored.depth() == depth_before);

        // Root scope properties should be visible (read-through)
        RC_ASSERT(restored.get_domain(var_a) == dom_a);

        // Child scope properties should be visible
        RC_ASSERT(restored.has_sign(var_b, sign_b));
    });
}

/// Test: Combined properties on a single symbol round-trip
static void test_serialization_combined_properties_roundtrip() {
    TEST_CASE("Combined properties round-trip");

    rc::check("Multiple properties on a single symbol survive round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(7);

        // Declare domain (must be compatible with sign)
        // Use Integer domain with Positive sign (compatible)
        ctx.assume_domain(var, Domain::Integer);
        ctx.assume_sign(var, Sign::Positive);
        ctx.current_properties().declare_parity(var, Parity::Odd);
        ctx.current_properties().declare_bounded(var, Boundedness::Bounded);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.get_domain(var) == Domain::Integer);
        RC_ASSERT(restored.has_sign(var, Sign::Positive));
        RC_ASSERT(restored.get_parity(var) == Parity::Odd);
        RC_ASSERT(restored.get_boundedness(var) == Boundedness::Bounded);
    });
}

/// Test: Simple relation round-trip (variable GT/LT 0)
static void test_serialization_simple_relation_roundtrip() {
    TEST_CASE("Simple relation round-trip");

    rc::check("Simple variable-vs-zero relations survive round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(8);

        // Create a simple relation: var > 0
        auto var_node = lamina::detail::make_node<VariableNode>(var);
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = lamina::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT);
        auto rel_expr = lamina::detail::expression_from_node(rel_node);
        ctx.assume(rel_expr);

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        // The sign property derived from the relation should be preserved
        RC_ASSERT(restored.has_sign(var, Sign::Positive));
    });
}

/// Test: Empty context round-trip
static void test_serialization_empty_roundtrip() {
    TEST_CASE("Empty context round-trip");

    rc::check("Empty context survives serialize/deserialize round-trip", []() {
        AssumptionContext ctx;

        std::string serialized = ctx.serialize();
        AssumptionContext restored = deserialize_success(serialized);

        RC_ASSERT(restored.depth() == ctx.depth());
        // No properties should be set
        std::string var = "nonexistent";
        RC_ASSERT(restored.get_domain(var) == Domain::Complex);
    });
}

static void test_checked_deserialization_contracts() {
    TEST_CASE("checked deserialization errors");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    auto rel_node = lamina::detail::make_node<RelationalNode>(
        lamina::detail::make_node<VariableNode>("x"),
        lamina::detail::make_node<NumberNode>(BigInt(0)),
        RelationalNode::Op::GT);
    auto relation = lamina::detail::expression_from_node(rel_node);
    ctx.assume(relation);

    auto ok = AssumptionContext::deserialize_checked(ctx.serialize());
    EXPECT_TRUE(ok.has_value(), "checked deserialize accepts serialized context");
    if (ok) {
        EXPECT_TRUE(ok.value().has_sign("x", Sign::Positive),
                    "checked deserialize preserves positive sign");
    }

    auto malformed = AssumptionContext::deserialize_checked("SCOPE 0\nRELATION x 0\nEND\n");
    EXPECT_TRUE(!malformed.has_value(), "checked deserialize rejects malformed relation");
    EXPECT_TRUE(malformed.error().code == CasErrc::ParseError,
                "checked deserialize reports ParseError for malformed input");

    auto contradictory = AssumptionContext::deserialize_checked(
        "SCOPE 0\nSIGN x Positive\nRELATION x LT 0\nEND\n");
    EXPECT_TRUE(!contradictory.has_value(),
                "checked deserialize rejects contradictory derived relation property");
    EXPECT_TRUE(contradictory.error().code == CasErrc::ParseError,
                "checked deserialize maps relation contradiction to ParseError");

    auto property_contradiction = AssumptionContext::deserialize_checked(
        "SCOPE 0\nDOMAIN n Natural\nSIGN n Negative\nEND\n");
    EXPECT_TRUE(!property_contradiction.has_value(),
                "checked deserialize rejects contradictory property declarations");
    EXPECT_TRUE(property_contradiction.error().code == CasErrc::ParseError,
                "checked deserialize maps property contradiction to ParseError");

    auto canonical = AssumptionContext::deserialize(
        "SCOPE 0\nDOMAIN n Natural\nSIGN n Negative\nEND\n");
    EXPECT_TRUE(!canonical.has_value(),
                "canonical deserialize returns the checked ParseError");
    EXPECT_TRUE(canonical.error().code == CasErrc::ParseError,
                "canonical deserialize preserves ParseError");
}


int main() {
    test_serialization_domain_roundtrip();
    test_serialization_sign_roundtrip();
    test_serialization_parity_roundtrip();
    test_serialization_boundedness_roundtrip();
    test_serialization_transcendental_roundtrip();
    test_serialization_finiteness_roundtrip();
    test_serialization_definiteness_roundtrip();
    test_serialization_multi_scope_roundtrip();
    test_serialization_combined_properties_roundtrip();
    test_serialization_simple_relation_roundtrip();
    test_serialization_empty_roundtrip();
    test_checked_deserialization_contracts();

    return TEST_REPORT();
}
