/**
 * @file test_assumption_relation_ext.cpp
 * @brief Unit tests for RelationStore extensions: reversed patterns, transitive closure, and cap.
 *
 * Tests:
 *   - Reversed patterns for all five operators (0 LT var → Positive, etc.)
 *   - Transitive chain of 3+ relations (x GT y, y GT z → x GT z stored)
 *   - 64-relation cap behavior (create enough relations to trigger the cap, verify it stops)
 *
 * Validates: Requirements 3.5, 24.1, 24.2, 24.3, 24.4, 24.5
 */

#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

using namespace lamina;

/// Helper: create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string& name) {
    return SymbolicExpr(std::make_shared<VariableNode>(name));
}

/// Helper: create a SymbolicExpr wrapping a NumberNode with value 0.
static SymbolicExpr make_zero() {
    return SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));
}

// ============================================================================
// Reversed pattern tests (Requirement 24)
// ============================================================================

void test_reversed_0_lt_var_positive() {
    TEST_CASE("Reversed: 0 LT var → Positive (Req 24.1)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("x");

    // 0 < x means x > 0 → Positive
    rs.add_relation(zero, var, RelationalNode::Op::LT, ps);

    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
                "0 LT x should derive Positive for x");
    EXPECT_TRUE(ps.has_sign("x", Sign::NonNegative),
                "Positive implies NonNegative");
    EXPECT_TRUE(ps.has_sign("x", Sign::NonZero),
                "Positive implies NonZero");
}

void test_reversed_0_gt_var_negative() {
    TEST_CASE("Reversed: 0 GT var → Negative (Req 24.2)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("y");

    // 0 > y means y < 0 → Negative
    rs.add_relation(zero, var, RelationalNode::Op::GT, ps);

    EXPECT_TRUE(ps.has_sign("y", Sign::Negative),
                "0 GT y should derive Negative for y");
    EXPECT_TRUE(ps.has_sign("y", Sign::NonPositive),
                "Negative implies NonPositive");
    EXPECT_TRUE(ps.has_sign("y", Sign::NonZero),
                "Negative implies NonZero");
}

void test_reversed_0_geq_var_nonpositive() {
    TEST_CASE("Reversed: 0 GEQ var → NonPositive (Req 24.3)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("z");

    // 0 >= z means z <= 0 → NonPositive
    rs.add_relation(zero, var, RelationalNode::Op::GEQ, ps);

    EXPECT_TRUE(ps.has_sign("z", Sign::NonPositive),
                "0 GEQ z should derive NonPositive for z");
}

void test_reversed_0_leq_var_nonnegative() {
    TEST_CASE("Reversed: 0 LEQ var → NonNegative (Req 24.4)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("w");

    // 0 <= w means w >= 0 → NonNegative
    rs.add_relation(zero, var, RelationalNode::Op::LEQ, ps);

    EXPECT_TRUE(ps.has_sign("w", Sign::NonNegative),
                "0 LEQ w should derive NonNegative for w");
}

void test_reversed_0_neq_var_nonzero() {
    TEST_CASE("Reversed: 0 NEQ var → NonZero (Req 24.5)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("a");

    // 0 != a means a != 0 → NonZero
    rs.add_relation(zero, var, RelationalNode::Op::NEQ, ps);

    EXPECT_TRUE(ps.has_sign("a", Sign::NonZero),
                "0 NEQ a should derive NonZero for a");
}

void test_reversed_all_operators_comprehensive() {
    TEST_CASE("Reversed: All five operators mapped correctly");

    struct TestCase {
        RelationalNode::Op op;
        Sign expected_sign;
        std::string desc;
    };

    std::vector<TestCase> cases = {
        {RelationalNode::Op::LT,  Sign::Positive,    "0 LT var → Positive"},
        {RelationalNode::Op::GT,  Sign::Negative,    "0 GT var → Negative"},
        {RelationalNode::Op::GEQ, Sign::NonPositive,  "0 GEQ var → NonPositive"},
        {RelationalNode::Op::LEQ, Sign::NonNegative, "0 LEQ var → NonNegative"},
        {RelationalNode::Op::NEQ, Sign::NonZero,     "0 NEQ var → NonZero"},
    };

    for (const auto& tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr zero = make_zero();
        SymbolicExpr var = make_var("v");

        rs.add_relation(zero, var, tc.op, ps);

        EXPECT_TRUE(ps.has_sign("v", tc.expected_sign), tc.desc);
    }
}

// ============================================================================
// Transitive closure tests (Requirement 3)
// ============================================================================

void test_transitive_chain_3_gt() {
    TEST_CASE("Transitive: x GT y, y GT z → x GT z (Req 3.1)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);

    // Transitive closure should deduce x GT z
    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
                "x GT y, y GT z → x GT z should be deduced");
}

void test_transitive_chain_geq_gt() {
    TEST_CASE("Transitive: x GEQ y, y GT z → x GT z (Req 3.2)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    rs.add_relation(x, y, RelationalNode::Op::GEQ, ps);
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);

    // GEQ + GT → GT
    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
                "x GEQ y, y GT z → x GT z should be deduced");
}

void test_transitive_chain_gt_geq() {
    TEST_CASE("Transitive: x GT y, y GEQ z → x GT z");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    rs.add_relation(y, z, RelationalNode::Op::GEQ, ps);

    // GT + GEQ → GT
    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
                "x GT y, y GEQ z → x GT z should be deduced");
}

void test_transitive_chain_geq_geq() {
    TEST_CASE("Transitive: x GEQ y, y GEQ z → x GEQ z (Req 3.3)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    rs.add_relation(x, y, RelationalNode::Op::GEQ, ps);
    rs.add_relation(y, z, RelationalNode::Op::GEQ, ps);

    // GEQ + GEQ → GEQ
    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GEQ),
                "x GEQ y, y GEQ z → x GEQ z should be deduced");
}

void test_transitive_chain_4_variables() {
    TEST_CASE("Transitive: chain of 4 variables a GT b GT c GT d → a GT d");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    rs.add_relation(a, b, RelationalNode::Op::GT, ps);
    rs.add_relation(b, c, RelationalNode::Op::GT, ps);
    rs.add_relation(c, d, RelationalNode::Op::GT, ps);

    // Should deduce a GT c, a GT d, b GT d
    EXPECT_TRUE(rs.has_relation(a, c, RelationalNode::Op::GT),
                "a GT b, b GT c → a GT c");
    EXPECT_TRUE(rs.has_relation(a, d, RelationalNode::Op::GT),
                "a GT c, c GT d → a GT d");
    EXPECT_TRUE(rs.has_relation(b, d, RelationalNode::Op::GT),
                "b GT c, c GT d → b GT d");
}

void test_transitive_no_closure_for_lt() {
    TEST_CASE("Transitive: LT does not participate in transitive closure directly");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // LT is not a transitive operator in the implementation (only GT/GEQ are)
    rs.add_relation(x, y, RelationalNode::Op::LT, ps);
    rs.add_relation(y, z, RelationalNode::Op::LT, ps);

    // No transitive deduction for LT (the implementation only handles GT/GEQ)
    EXPECT_FALSE(rs.has_relation(x, z, RelationalNode::Op::LT),
                 "LT does not participate in transitive closure");
}

// ============================================================================
// 64-relation cap test (Requirement 3.5)
// ============================================================================

void test_transitive_cap_64() {
    TEST_CASE("Transitive: 64-relation cap stops deduction (Req 3.5)");

    RelationStore rs;
    PropertyStore ps;

    // Create a long chain: v0 GT v1, v1 GT v2, ..., v(N-1) GT vN
    // With N variables in a chain, adding the last relation can trigger many
    // transitive deductions. We need enough variables so that the total
    // deductions from a single add_relation call would exceed 64.
    //
    // Strategy: build a chain of 12 variables first (v0..v11), then add one
    // more relation that connects to the chain. With 12 existing nodes,
    // adding v12 would only produce 12 new relations. Instead, we build
    // a large chain incrementally and check that the total stored relations
    // are bounded.
    //
    // Actually, the cap is per add_relation call. Let's create a scenario
    // where a single add_relation triggers many deductions:
    // First, add many independent chains that share an endpoint.

    // Build a star pattern: many variables all GT than a central variable "center"
    // Then add "center GT bottom" — this should trigger deductions for all
    // star variables GT bottom.
    
    // With 70 star variables + center GT bottom, the single add_relation(center, bottom, GT)
    // would try to deduce 70 new relations (star_i GT bottom), but cap at 64.

    SymbolicExpr center = make_var("center");
    SymbolicExpr bottom = make_var("bottom");

    // Add 70 relations: star_i GT center
    const int num_star = 70;
    for (int i = 0; i < num_star; ++i) {
        SymbolicExpr star_i = make_var("star_" + std::to_string(i));
        rs.add_relation(star_i, center, RelationalNode::Op::GT, ps);
    }

    // Now add: center GT bottom
    // This should trigger backward chaining: for each star_i GT center,
    // deduce star_i GT bottom. But cap at 64.
    size_t relations_before = rs.get_relations().size();
    rs.add_relation(center, bottom, RelationalNode::Op::GT, ps);
    size_t relations_after = rs.get_relations().size();

    // The new relations added should be: 1 (center GT bottom) + at most 64 deduced
    size_t new_relations = relations_after - relations_before;
    EXPECT_TRUE(new_relations <= 65,
                "Should add at most 1 + 64 = 65 new relations (cap at 64 deductions)");

    // Verify that at least some deductions were made
    EXPECT_TRUE(new_relations > 1,
                "Should have deduced at least some transitive relations");

    // Verify the cap was actually hit (70 potential deductions, only 64 allowed)
    // The 1 is for the directly added relation, the rest are deductions
    size_t deductions = new_relations - 1;
    EXPECT_TRUE(deductions <= 64,
                "Deductions should be capped at 64");
    EXPECT_TRUE(deductions == 64,
                "With 70 potential deductions, exactly 64 should be produced");

    // Verify that some star variables DO have the deduced relation
    int found_count = 0;
    for (int i = 0; i < num_star; ++i) {
        SymbolicExpr star_i = make_var("star_" + std::to_string(i));
        if (rs.has_relation(star_i, bottom, RelationalNode::Op::GT)) {
            ++found_count;
        }
    }
    EXPECT_TRUE(found_count == 64,
                "Exactly 64 star variables should have deduced GT bottom relation");

    // Verify that some star variables do NOT have the deduced relation (cap hit)
    EXPECT_TRUE(found_count < num_star,
                "Not all 70 star variables should have the deduced relation (cap hit)");
}

// ============================================================================
// main
// ============================================================================

int main() {
    // Reversed pattern tests (Req 24)
    test_reversed_0_lt_var_positive();
    test_reversed_0_gt_var_negative();
    test_reversed_0_geq_var_nonpositive();
    test_reversed_0_leq_var_nonnegative();
    test_reversed_0_neq_var_nonzero();
    test_reversed_all_operators_comprehensive();

    // Transitive closure tests (Req 3)
    test_transitive_chain_3_gt();
    test_transitive_chain_geq_gt();
    test_transitive_chain_gt_geq();
    test_transitive_chain_geq_geq();
    test_transitive_chain_4_variables();
    test_transitive_no_closure_for_lt();

    // Cap test (Req 3.5)
    test_transitive_cap_64();

    return TEST_REPORT();
}
