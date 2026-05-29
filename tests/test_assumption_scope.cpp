/**
 * @file test_assumption_scope.cpp
 * @brief Property tests for AssumptionContext scope management (Properties 10-11).
 *
 * Feature: assumption-system
 * Validates: Requirements 9.1-9.3, 9.6, 9.7
 *
 * Property 10: Scope push/pop round-trip
 * Property 11: Scope shadowing
 */

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "interval.hpp"
#include <stdexcept>
#include <vector>
#include <string>

using namespace lamina;

// ============================================================
// Helpers
// ============================================================

static const std::vector<Domain> ALL_DOMAINS = {
    Domain::Complex, Domain::Real, Domain::Rational,
    Domain::Integer, Domain::Natural, Domain::PositiveInt
};

static const std::vector<Sign> ALL_SIGNS = {
    Sign::Positive, Sign::Negative, Sign::NonNegative,
    Sign::NonPositive, Sign::Zero, Sign::NonZero
};

static const std::vector<std::string> TEST_SYMBOLS = {
    "x", "y", "alpha", "longVar123", "a_b"
};

static std::string domain_name(Domain d) {
    switch (d) {
        case Domain::Complex:     return "Complex";
        case Domain::Real:        return "Real";
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
// Property 10: Scope push/pop round-trip
// **Validates: Requirements 9.1, 9.2, 9.3, 9.7**
//
// For any set of property and relation declarations made within
// a pushed scope, popping that scope should restore all query
// results to exactly the values they returned before the push
// was called; declarations in the child scope should not be
// visible after pop.
// ============================================================

void test_property10_domain_roundtrip() {
    TEST_CASE("Property 10: Domain declarations in child scope not visible after pop");
    // For each domain, declare in child scope, verify visible in child,
    // then pop and verify not visible (reverts to default Complex).
    for (const auto& sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            if (d == Domain::Complex) continue; // Complex is default, skip

            AssumptionContext ctx;
            // Before push: domain should be Complex (default)
            Domain before = ctx.get_domain(sym);
            EXPECT_TRUE(before == Domain::Complex,
                sym + " domain is Complex before push");

            ctx.push();
            ctx.assume_domain(sym, d);

            // In child scope: domain should be at least d
            EXPECT_TRUE(ctx.has_domain(sym, d),
                sym + " has " + domain_name(d) + " in child scope");

            ctx.pop();

            // After pop: domain should be back to Complex
            Domain after = ctx.get_domain(sym);
            EXPECT_TRUE(after == before,
                sym + " domain restored to Complex after pop (was " +
                domain_name(d) + " in child)");
        }
    }
}

void test_property10_sign_roundtrip() {
    TEST_CASE("Property 10: Sign declarations in child scope not visible after pop");
    for (const auto& sym : TEST_SYMBOLS) {
        for (Sign s : ALL_SIGNS) {
            AssumptionContext ctx;
            // Before push: no signs
            auto signs_before = ctx.get_signs(sym);
            EXPECT_TRUE(signs_before.empty(),
                sym + " has no signs before push");

            ctx.push();
            ctx.assume_sign(sym, s);

            // In child scope: sign should be present
            EXPECT_TRUE(ctx.has_sign(sym, s),
                sym + " has " + sign_name(s) + " in child scope");

            ctx.pop();

            // After pop: signs should be empty again
            auto signs_after = ctx.get_signs(sym);
            EXPECT_TRUE(signs_after.empty(),
                sym + " signs empty after pop (was " +
                sign_name(s) + " in child)");
        }
    }
}

void test_property10_parity_roundtrip() {
    TEST_CASE("Property 10: Parity declarations in child scope not visible after pop");
    AssumptionContext ctx;

    // Before push: parity is Unknown
    EXPECT_TRUE(ctx.get_parity("x") == Parity::Unknown,
        "x parity Unknown before push");

    ctx.push();
    ctx.current_properties().declare_parity("x", Parity::Even);
    EXPECT_TRUE(ctx.get_parity("x") == Parity::Even,
        "x parity Even in child scope");

    ctx.pop();
    EXPECT_TRUE(ctx.get_parity("x") == Parity::Unknown,
        "x parity restored to Unknown after pop");
}

void test_property10_boundedness_roundtrip() {
    TEST_CASE("Property 10: Boundedness declarations in child scope not visible after pop");
    AssumptionContext ctx;

    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Unknown,
        "x boundedness Unknown before push");

    ctx.push();
    ctx.current_properties().declare_bounded("x", Boundedness::Bounded);
    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Bounded,
        "x boundedness Bounded in child scope");

    ctx.pop();
    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Unknown,
        "x boundedness restored to Unknown after pop");
}

void test_property10_relation_roundtrip() {
    TEST_CASE("Property 10: Relations in child scope not visible after pop");
    AssumptionContext ctx;

    // Create a simple relation: x > 0
    auto var_x = std::make_shared<SymbolicExpr>(
        std::make_shared<VariableNode>("x"));
    auto zero = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(BigInt(0)));

    // Before push: no relations, no sign for x
    EXPECT_TRUE(ctx.current_relations().get_relations().empty(),
        "No relations before push");
    EXPECT_FALSE(ctx.has_sign("x", Sign::Positive),
        "x not Positive before push");

    ctx.push();
    ctx.current_relations().add_relation(*var_x, *zero,
        RelationalNode::Op::GT, ctx.current_properties());

    // In child scope: relation present, x is Positive
    EXPECT_TRUE(ctx.has_sign("x", Sign::Positive),
        "x is Positive in child scope (from relation x > 0)");

    ctx.pop();

    // After pop: sign should be gone
    EXPECT_FALSE(ctx.has_sign("x", Sign::Positive),
        "x not Positive after pop");
    EXPECT_TRUE(ctx.current_relations().get_relations().empty(),
        "No relations after pop");
}

void test_property10_multiple_declarations_roundtrip() {
    TEST_CASE("Property 10: Multiple declarations in child scope all reverted on pop");
    AssumptionContext ctx;

    // Set up parent state
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("y", Sign::Positive);

    // Record parent state
    Domain x_domain_before = ctx.get_domain("x");
    auto y_signs_before = ctx.get_signs("y");
    auto z_signs_before = ctx.get_signs("z");
    Parity z_parity_before = ctx.get_parity("z");

    ctx.push();

    // Make multiple declarations in child
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_sign("z", Sign::Negative);
    ctx.current_properties().declare_parity("z", Parity::Odd);

    // Verify child state
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer in child");
    EXPECT_TRUE(ctx.has_sign("z", Sign::Negative),
        "z is Negative in child");
    EXPECT_TRUE(ctx.get_parity("z") == Parity::Odd,
        "z is Odd in child");

    ctx.pop();

    // Verify all restored
    EXPECT_TRUE(ctx.get_domain("x") == x_domain_before,
        "x domain restored after pop");
    EXPECT_TRUE(ctx.get_signs("y") == y_signs_before,
        "y signs unchanged after pop");
    EXPECT_TRUE(ctx.get_signs("z") == z_signs_before,
        "z signs restored (empty) after pop");
    EXPECT_TRUE(ctx.get_parity("z") == z_parity_before,
        "z parity restored (Unknown) after pop");
}

void test_property10_parent_declarations_survive_pop() {
    TEST_CASE("Property 10: Parent scope declarations survive child push/pop");
    AssumptionContext ctx;

    // Declare in root scope
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_sign("x", Sign::Positive);

    ctx.push();
    // Child scope doesn't touch x
    ctx.assume_sign("y", Sign::Negative);
    ctx.pop();

    // x should still have its root declarations
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x domain Integer survives child push/pop");
    EXPECT_TRUE(ctx.has_sign("x", Sign::Positive),
        "x sign Positive survives child push/pop");
}

void test_property10_nested_push_pop_roundtrip() {
    TEST_CASE("Property 10: Nested push/pop restores correctly at each level");
    AssumptionContext ctx;

    // Root: x is Real
    ctx.assume_domain("x", Domain::Real);

    ctx.push(); // depth 2
    ctx.assume_domain("x", Domain::Integer);
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer at depth 2");

    ctx.push(); // depth 3
    ctx.assume_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(ctx.get_domain("x") == Domain::PositiveInt,
        "x is PositiveInt at depth 3");

    ctx.pop(); // back to depth 2
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer after popping depth 3");

    ctx.pop(); // back to depth 1 (root)
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Real,
        "x is Real after popping depth 2");
}

void test_property10_depth_changes() {
    TEST_CASE("Property 10: Depth increases on push and decreases on pop");
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.depth() == 1, "Initial depth is 1");

    ctx.push();
    EXPECT_TRUE(ctx.depth() == 2, "Depth is 2 after push");

    ctx.push();
    EXPECT_TRUE(ctx.depth() == 3, "Depth is 3 after second push");

    ctx.pop();
    EXPECT_TRUE(ctx.depth() == 2, "Depth is 2 after pop");

    ctx.pop();
    EXPECT_TRUE(ctx.depth() == 1, "Depth is 1 after second pop");
}

// ============================================================
// Property 11: Scope shadowing
// **Validates: Requirements 9.6**
//
// For any symbol with a property declared in a parent scope,
// declaring a different (non-contradictory) property for the
// same symbol in a child scope should shadow the parent's value
// for queries within the child scope, while the parent's value
// remains accessible after pop.
// ============================================================

void test_property11_domain_shadowing() {
    TEST_CASE("Property 11: Child scope domain shadows parent domain");
    // For each pair of domains where child is more specific than parent,
    // the child's domain should shadow the parent's.
    for (const auto& sym : TEST_SYMBOLS) {
        // Parent declares Real, child declares Integer (more specific)
        AssumptionContext ctx;
        ctx.assume_domain(sym, Domain::Real);

        Domain parent_domain = ctx.get_domain(sym);
        EXPECT_TRUE(parent_domain == Domain::Real,
            sym + " is Real in parent");

        ctx.push();
        ctx.assume_domain(sym, Domain::Integer);

        // Child scope: should see Integer
        EXPECT_TRUE(ctx.get_domain(sym) == Domain::Integer,
            sym + " is Integer in child (shadows Real)");

        ctx.pop();

        // After pop: should see Real again
        EXPECT_TRUE(ctx.get_domain(sym) == Domain::Real,
            sym + " is Real again after pop");
    }
}

void test_property11_sign_shadowing() {
    TEST_CASE("Property 11: Child scope sign shadows parent sign");
    // Parent: x is NonNegative; Child: x is Positive (compatible, more specific)
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    EXPECT_TRUE(ctx.has_sign("x", Sign::NonNegative),
        "x is NonNegative in parent");
    EXPECT_FALSE(ctx.has_sign("x", Sign::Positive),
        "x is NOT Positive in parent (only NonNegative)");

    ctx.push();
    ctx.assume_sign("x", Sign::Positive);

    // Child scope: should see Positive (and its implications)
    EXPECT_TRUE(ctx.has_sign("x", Sign::Positive),
        "x is Positive in child (shadows parent)");
    EXPECT_TRUE(ctx.has_sign("x", Sign::NonNegative),
        "x is NonNegative in child (implied by Positive)");
    EXPECT_TRUE(ctx.has_sign("x", Sign::NonZero),
        "x is NonZero in child (implied by Positive)");

    ctx.pop();

    // After pop: should see only NonNegative again
    EXPECT_TRUE(ctx.has_sign("x", Sign::NonNegative),
        "x is NonNegative after pop (parent value)");
    EXPECT_FALSE(ctx.has_sign("x", Sign::Positive),
        "x is NOT Positive after pop");
    EXPECT_FALSE(ctx.has_sign("x", Sign::NonZero),
        "x is NOT NonZero after pop (was only NonNegative in parent)");
}

void test_property11_parity_shadowing() {
    TEST_CASE("Property 11: Child scope parity shadows parent parity");
    AssumptionContext ctx;

    // Parent: x is Even
    ctx.current_properties().declare_parity("x", Parity::Even);
    EXPECT_TRUE(ctx.get_parity("x") == Parity::Even,
        "x is Even in parent");

    ctx.push();
    // Child: x is Odd (different parity — this is a new scope, no contradiction
    // because child has its own PropertyStore)
    ctx.current_properties().declare_parity("x", Parity::Odd);

    EXPECT_TRUE(ctx.get_parity("x") == Parity::Odd,
        "x is Odd in child (shadows Even)");

    ctx.pop();

    EXPECT_TRUE(ctx.get_parity("x") == Parity::Even,
        "x is Even after pop (parent value restored)");
}

void test_property11_boundedness_shadowing() {
    TEST_CASE("Property 11: Child scope boundedness shadows parent boundedness");
    AssumptionContext ctx;

    // Parent: x is Unbounded
    ctx.current_properties().declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Unbounded,
        "x is Unbounded in parent");

    ctx.push();
    // Child: x is Bounded (different — new scope, no contradiction)
    ctx.current_properties().declare_bounded("x", Boundedness::Bounded);

    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Bounded,
        "x is Bounded in child (shadows Unbounded)");

    ctx.pop();

    EXPECT_TRUE(ctx.get_boundedness("x") == Boundedness::Unbounded,
        "x is Unbounded after pop (parent value restored)");
}

void test_property11_domain_shadowing_all_pairs() {
    TEST_CASE("Property 11: Domain shadowing for various parent/child domain pairs");
    // Test multiple domain pairs where child is different from parent
    struct DomainPair {
        Domain parent;
        Domain child;
    };
    std::vector<DomainPair> pairs = {
        {Domain::Real, Domain::Integer},
        {Domain::Real, Domain::Natural},
        {Domain::Rational, Domain::PositiveInt},
        {Domain::Integer, Domain::PositiveInt},
        {Domain::Real, Domain::PositiveInt},
    };

    for (const auto& p : pairs) {
        AssumptionContext ctx;
        ctx.assume_domain("x", p.parent);

        ctx.push();
        ctx.assume_domain("x", p.child);

        EXPECT_TRUE(ctx.get_domain("x") == p.child,
            "x is " + domain_name(p.child) + " in child (parent was " +
            domain_name(p.parent) + ")");

        ctx.pop();

        EXPECT_TRUE(ctx.get_domain("x") == p.parent,
            "x is " + domain_name(p.parent) + " after pop");
    }
}

void test_property11_sign_shadowing_various() {
    TEST_CASE("Property 11: Sign shadowing for various parent/child sign pairs");
    // Test sign pairs where child is different but non-contradictory within
    // its own scope (each scope has independent PropertyStore)
    struct SignPair {
        Sign parent;
        Sign child;
    };
    std::vector<SignPair> pairs = {
        {Sign::NonNegative, Sign::Positive},
        {Sign::NonPositive, Sign::Negative},
        {Sign::NonZero, Sign::Positive},
        {Sign::NonZero, Sign::Negative},
        {Sign::NonNegative, Sign::Zero},
    };

    for (const auto& p : pairs) {
        AssumptionContext ctx;
        ctx.assume_sign("x", p.parent);

        ctx.push();
        ctx.assume_sign("x", p.child);

        // Child should see child's sign
        EXPECT_TRUE(ctx.has_sign("x", p.child),
            "x has " + sign_name(p.child) + " in child (parent was " +
            sign_name(p.parent) + ")");

        ctx.pop();

        // Parent's sign should be restored
        EXPECT_TRUE(ctx.has_sign("x", p.parent),
            "x has " + sign_name(p.parent) + " after pop");
    }
}

void test_property11_child_does_not_modify_parent() {
    TEST_CASE("Property 11: Child scope declarations do not modify parent scope");
    AssumptionContext ctx;

    // Parent: x is Real, y is Positive
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("y", Sign::Positive);

    ctx.push();

    // Child: override x to Integer, override y to Negative
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_sign("y", Sign::Negative);

    // Verify child sees child values
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer in child");
    EXPECT_TRUE(ctx.has_sign("y", Sign::Negative),
        "y is Negative in child");

    ctx.pop();

    // Parent values should be completely unchanged
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Real,
        "x is still Real in parent (not modified by child)");
    EXPECT_TRUE(ctx.has_sign("y", Sign::Positive),
        "y is still Positive in parent (not modified by child)");
    EXPECT_FALSE(ctx.has_sign("y", Sign::Negative),
        "y is NOT Negative in parent");
}

void test_property11_read_through_undeclared_in_child() {
    TEST_CASE("Property 11: Child scope reads through to parent for undeclared symbols");
    AssumptionContext ctx;

    // Parent: x is Integer, y is Positive
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_sign("y", Sign::Positive);

    ctx.push();

    // Child doesn't declare x or y — should read through to parent
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x reads through to Integer from parent");
    EXPECT_TRUE(ctx.has_sign("y", Sign::Positive),
        "y reads through to Positive from parent");

    // Child declares z — only z is new
    ctx.assume_sign("z", Sign::NonNegative);
    EXPECT_TRUE(ctx.has_sign("z", Sign::NonNegative),
        "z is NonNegative in child");

    ctx.pop();

    // z should not be visible after pop
    EXPECT_FALSE(ctx.has_sign("z", Sign::NonNegative),
        "z not visible after pop");
    // x and y still visible from root
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x still Integer in root");
    EXPECT_TRUE(ctx.has_sign("y", Sign::Positive),
        "y still Positive in root");
}

void test_property11_multi_level_shadowing() {
    TEST_CASE("Property 11: Multi-level shadowing (grandchild shadows child shadows parent)");
    AssumptionContext ctx;

    // Root: x is Real
    ctx.assume_domain("x", Domain::Real);

    ctx.push(); // Level 2
    // Level 2: x is Integer (shadows Real)
    ctx.assume_domain("x", Domain::Integer);
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer at level 2");

    ctx.push(); // Level 3
    // Level 3: x is PositiveInt (shadows Integer)
    ctx.assume_domain("x", Domain::PositiveInt);
    EXPECT_TRUE(ctx.get_domain("x") == Domain::PositiveInt,
        "x is PositiveInt at level 3");

    ctx.pop(); // Back to level 2
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Integer,
        "x is Integer at level 2 after popping level 3");

    ctx.pop(); // Back to root
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Real,
        "x is Real at root after popping level 2");
}

void test_property11_different_symbols_independent() {
    TEST_CASE("Property 11: Shadowing is per-symbol — different symbols are independent");
    AssumptionContext ctx;

    // Root: x is Real, y is Integer
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Integer);

    ctx.push();
    // Child: only shadow x
    ctx.assume_domain("x", Domain::Natural);

    // x is shadowed, y reads through
    EXPECT_TRUE(ctx.get_domain("x") == Domain::Natural,
        "x is Natural in child (shadowed)");
    EXPECT_TRUE(ctx.get_domain("y") == Domain::Integer,
        "y is Integer in child (read-through, not shadowed)");

    ctx.pop();

    EXPECT_TRUE(ctx.get_domain("x") == Domain::Real,
        "x is Real after pop");
    EXPECT_TRUE(ctx.get_domain("y") == Domain::Integer,
        "y is Integer after pop (unchanged)");
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 10: Scope push/pop round-trip
    test_property10_domain_roundtrip();
    test_property10_sign_roundtrip();
    test_property10_parity_roundtrip();
    test_property10_boundedness_roundtrip();
    test_property10_relation_roundtrip();
    test_property10_multiple_declarations_roundtrip();
    test_property10_parent_declarations_survive_pop();
    test_property10_nested_push_pop_roundtrip();
    test_property10_depth_changes();

    // Property 11: Scope shadowing
    test_property11_domain_shadowing();
    test_property11_sign_shadowing();
    test_property11_parity_shadowing();
    test_property11_boundedness_shadowing();
    test_property11_domain_shadowing_all_pairs();
    test_property11_sign_shadowing_various();
    test_property11_child_does_not_modify_parent();
    test_property11_read_through_undeclared_in_child();
    test_property11_multi_level_shadowing();
    test_property11_different_symbols_independent();

    return TEST_REPORT();
}
