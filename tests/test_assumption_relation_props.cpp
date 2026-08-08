
#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <vector>

using namespace lamina;


/// Create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string& name) {
    return lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>(name));
}

/// Create a SymbolicExpr wrapping a NumberNode with value 0.
static SymbolicExpr make_zero() {
    return lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(0)));
}


static void test_transitive_gt_gt_chain() {
    TEST_CASE("Property 3: GT + GT chain deduces GT (x > y, y > z => x > z)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x > y
    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    // Add y > z — should trigger transitive closure: x > z
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);

    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
        "x > y, y > z => x > z (GT+GT => GT)");
}

static void test_transitive_geq_gt_chain() {
    TEST_CASE("Property 3: GEQ + GT chain deduces GT (x >= y, y > z => x > z)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x >= y
    rs.add_relation(x, y, RelationalNode::Op::GEQ, ps);
    // Add y > z — should trigger: x > z (GEQ+GT => GT)
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);

    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
        "x >= y, y > z => x > z (GEQ+GT => GT)");
}

static void test_transitive_gt_geq_chain() {
    TEST_CASE("Property 3: GT + GEQ chain deduces GT (x > y, y >= z => x > z)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x > y
    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    // Add y >= z — should trigger: x > z (GT+GEQ => GT)
    rs.add_relation(y, z, RelationalNode::Op::GEQ, ps);

    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
        "x > y, y >= z => x > z (GT+GEQ => GT)");
}

static void test_transitive_geq_geq_chain() {
    TEST_CASE("Property 3: GEQ + GEQ chain deduces GEQ (x >= y, y >= z => x >= z)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x >= y
    rs.add_relation(x, y, RelationalNode::Op::GEQ, ps);
    // Add y >= z — should trigger: x >= z (GEQ+GEQ => GEQ)
    rs.add_relation(y, z, RelationalNode::Op::GEQ, ps);

    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GEQ),
        "x >= y, y >= z => x >= z (GEQ+GEQ => GEQ)");
}

static void test_transitive_longer_chain() {
    TEST_CASE("Property 3: Longer chain a > b > c > d deduces a > d");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    rs.add_relation(a, b, RelationalNode::Op::GT, ps);
    rs.add_relation(b, c, RelationalNode::Op::GT, ps);
    rs.add_relation(c, d, RelationalNode::Op::GT, ps);

    // After adding a>b, b>c: a>c should be deduced
    EXPECT_TRUE(rs.has_relation(a, c, RelationalNode::Op::GT),
        "a > b, b > c => a > c");

    // After adding c>d: a>d, b>d should be deduced
    EXPECT_TRUE(rs.has_relation(a, d, RelationalNode::Op::GT),
        "a > b > c > d => a > d");
    EXPECT_TRUE(rs.has_relation(b, d, RelationalNode::Op::GT),
        "b > c > d => b > d");
}

static void test_transitive_mixed_chain() {
    TEST_CASE("Property 3: Mixed chain a > b >= c > d deduces a > d");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    rs.add_relation(a, b, RelationalNode::Op::GT, ps);
    rs.add_relation(b, c, RelationalNode::Op::GEQ, ps);
    rs.add_relation(c, d, RelationalNode::Op::GT, ps);

    // a > b >= c > d => a > d (GT+GEQ => GT, then GT+GT => GT)
    EXPECT_TRUE(rs.has_relation(a, d, RelationalNode::Op::GT),
        "a > b >= c > d => a > d");
}

static void test_transitive_backward_chaining() {
    TEST_CASE("Property 3: Backward chaining — adding earlier link deduces relation");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add y > z first
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);
    // Then add x > y — should trigger backward chain: x > z
    rs.add_relation(x, y, RelationalNode::Op::GT, ps);

    EXPECT_TRUE(rs.has_relation(x, z, RelationalNode::Op::GT),
        "y > z then x > y => x > z (backward chaining)");
}

static void test_transitive_no_duplicate_deduction() {
    TEST_CASE("Property 3: Transitive closure does not add duplicate relations");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    rs.add_relation(y, z, RelationalNode::Op::GT, ps);

    // x > z should be deduced once
    size_t count = 0;
    for (const auto& rel : rs.get_relations()) {
        if (lamina::detail::node(rel.lhs) && lamina::detail::node(rel.rhs)) {
            auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(rel.lhs));
            auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(rel.rhs));
            if (lhs_var && rhs_var &&
                lhs_var->name() == "x" && rhs_var->name() == "z" &&
                rel.op == RelationalNode::Op::GT) {
                ++count;
            }
        }
    }
    EXPECT_TRUE(count == 1, "x > z should appear exactly once (no duplicates)");
}

static void test_transitive_lt_leq_not_transitive() {
    TEST_CASE("Property 3: LT and LEQ do NOT participate in transitive closure");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // LT relations should be stored but not trigger transitive closure
    rs.add_relation(x, y, RelationalNode::Op::LT, ps);
    rs.add_relation(y, z, RelationalNode::Op::LT, ps);

    // x < z should NOT be deduced (only GT/GEQ participate)
    EXPECT_FALSE(rs.has_relation(x, z, RelationalNode::Op::LT),
        "LT does not participate in transitive closure");
}

static void test_transitive_cap_at_64() {
    TEST_CASE("Property 3: Transitive closure capped at 64 new relations per add_relation");

    RelationStore rs;
    PropertyStore ps;

    // Create a long chain of variables: v0 > v1 > v2 > ... > v_N
    // Adding each link deduces transitive relations with all previous links.
    // We need enough variables that a single add_relation would try to produce > 64 deductions.
    // With N existing links forming a chain of length N, adding one more link at the end
    // would try to deduce N new relations. So we need N > 64.
    // Build a chain of 70 variables first (v0 > v1 > ... > v69)

    std::vector<SymbolicExpr> vars;
    for (int i = 0; i < 70; ++i) {
        vars.push_back(make_var("v" + std::to_string(i)));
    }

    // Add the first 69 links one by one (this builds up the chain)
    for (int i = 0; i < 69; ++i) {
        rs.add_relation(vars[i], vars[i + 1], RelationalNode::Op::GT, ps);
    }

    // Count relations before adding the 70th variable link
    size_t before_count = rs.get_relations().size();

    // Now add v_all > v0 which would try to chain through all 69 existing links
    SymbolicExpr v_new = make_var("v_new");
    rs.add_relation(v_new, vars[0], RelationalNode::Op::GT, ps);

    size_t after_count = rs.get_relations().size();
    // The new relations added should be at most 64 + 1 (the original relation itself)
    size_t new_relations = after_count - before_count;

    // The directly added relation is 1, plus at most 64 deduced
    EXPECT_TRUE(new_relations <= 65,
        "At most 64 new deduced relations + 1 original per add_relation call (got " +
        std::to_string(new_relations) + ")");
}

static void test_transitive_sign_derivation_from_chain() {
    TEST_CASE("Property 3: Transitive closure derives sign when chain reaches zero");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr zero = make_zero();

    // x > y, y > 0 => x > 0 => x is Positive
    rs.add_relation(x, y, RelationalNode::Op::GT, ps);
    rs.add_relation(y, zero, RelationalNode::Op::GT, ps);

    // y > 0 directly derives Positive for y
    EXPECT_TRUE(ps.has_sign("y", Sign::Positive),
        "y > 0 derives Positive for y");

    // x > 0 should be deduced transitively, deriving Positive for x
    EXPECT_TRUE(rs.has_relation(x, zero, RelationalNode::Op::GT),
        "x > y > 0 => x > 0 deduced");
    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
        "x > y > 0 => x is Positive (sign derived from transitive closure)");
}


static void test_reversed_0_lt_var_positive() {
    TEST_CASE("Property 12: 0 LT var => var is Positive");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    // 0 < x means x > 0, so x is Positive
    rs.add_relation(zero, x, RelationalNode::Op::LT, ps);

    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
        "0 < x => x is Positive");
    EXPECT_TRUE(ps.has_sign("x", Sign::NonNegative),
        "0 < x => x is NonNegative (implied)");
    EXPECT_TRUE(ps.has_sign("x", Sign::NonZero),
        "0 < x => x is NonZero (implied)");
}

static void test_reversed_0_gt_var_negative() {
    TEST_CASE("Property 12: 0 GT var => var is Negative");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr y = make_var("y");

    // 0 > y means y < 0, so y is Negative
    rs.add_relation(zero, y, RelationalNode::Op::GT, ps);

    EXPECT_TRUE(ps.has_sign("y", Sign::Negative),
        "0 > y => y is Negative");
    EXPECT_TRUE(ps.has_sign("y", Sign::NonPositive),
        "0 > y => y is NonPositive (implied)");
    EXPECT_TRUE(ps.has_sign("y", Sign::NonZero),
        "0 > y => y is NonZero (implied)");
}

static void test_reversed_0_geq_var_nonpositive() {
    TEST_CASE("Property 12: 0 GEQ var => var is NonPositive");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr z = make_var("z");

    // 0 >= z means z <= 0, so z is NonPositive
    rs.add_relation(zero, z, RelationalNode::Op::GEQ, ps);

    EXPECT_TRUE(ps.has_sign("z", Sign::NonPositive),
        "0 >= z => z is NonPositive");
}

static void test_reversed_0_leq_var_nonnegative() {
    TEST_CASE("Property 12: 0 LEQ var => var is NonNegative");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr w = make_var("w");

    // 0 <= w means w >= 0, so w is NonNegative
    rs.add_relation(zero, w, RelationalNode::Op::LEQ, ps);

    EXPECT_TRUE(ps.has_sign("w", Sign::NonNegative),
        "0 <= w => w is NonNegative");
}

static void test_reversed_0_neq_var_nonzero() {
    TEST_CASE("Property 12: 0 NEQ var => var is NonZero");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr v = make_var("v");

    // 0 != v means v != 0, so v is NonZero
    rs.add_relation(zero, v, RelationalNode::Op::NEQ, ps);

    EXPECT_TRUE(ps.has_sign("v", Sign::NonZero),
        "0 != v => v is NonZero");
}

static void test_reversed_all_operators_comprehensive() {
    TEST_CASE("Property 12: All reversed operators mapped correctly");

    struct TestCase {
        RelationalNode::Op op;
        Sign expected_sign;
        std::string desc;
    };

    std::vector<TestCase> cases = {
        {RelationalNode::Op::LT,  Sign::Positive,    "0 LT var => Positive"},
        {RelationalNode::Op::GT,  Sign::Negative,    "0 GT var => Negative"},
        {RelationalNode::Op::GEQ, Sign::NonPositive,  "0 GEQ var => NonPositive"},
        {RelationalNode::Op::LEQ, Sign::NonNegative, "0 LEQ var => NonNegative"},
        {RelationalNode::Op::NEQ, Sign::NonZero,     "0 NEQ var => NonZero"},
    };

    for (const auto& tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr zero = make_zero();
        SymbolicExpr var = make_var("t");

        rs.add_relation(zero, var, tc.op, ps);

        EXPECT_TRUE(ps.has_sign("t", tc.expected_sign), tc.desc);
    }
}

static void test_reversed_multiple_variables() {
    TEST_CASE("Property 12: Reversed pattern works for multiple variables");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();

    // 0 < a, 0 > b, 0 >= c, 0 <= d, 0 != e
    rs.add_relation(zero, make_var("a"), RelationalNode::Op::LT, ps);
    rs.add_relation(zero, make_var("b"), RelationalNode::Op::GT, ps);
    rs.add_relation(zero, make_var("c"), RelationalNode::Op::GEQ, ps);
    rs.add_relation(zero, make_var("d"), RelationalNode::Op::LEQ, ps);
    rs.add_relation(zero, make_var("e"), RelationalNode::Op::NEQ, ps);

    EXPECT_TRUE(ps.has_sign("a", Sign::Positive), "a is Positive");
    EXPECT_TRUE(ps.has_sign("b", Sign::Negative), "b is Negative");
    EXPECT_TRUE(ps.has_sign("c", Sign::NonPositive), "c is NonPositive");
    EXPECT_TRUE(ps.has_sign("d", Sign::NonNegative), "d is NonNegative");
    EXPECT_TRUE(ps.has_sign("e", Sign::NonZero), "e is NonZero");
}

static void test_reversed_non_variable_rhs_no_derivation() {
    TEST_CASE("Property 12: 0 op composite_expr does NOT derive sign");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();

    // RHS is a composite expression (x + y), not a single variable
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto y_node = lamina::detail::make_node<VariableNode>("y");
    auto add_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, y_node});
    auto composite = lamina::detail::expression_from_node(add_node);
    rs.add_relation(zero, composite, RelationalNode::Op::LT, ps);

    // Neither x nor y should have sign derived
    EXPECT_FALSE(ps.has_sign("x", Sign::Positive),
        "0 < (x+y) should not derive sign for x");
    EXPECT_FALSE(ps.has_sign("y", Sign::Positive),
        "0 < (x+y) should not derive sign for y");

    // But the relation should still be stored
    EXPECT_TRUE(rs.has_relation(zero, composite, RelationalNode::Op::LT),
        "Relation with composite RHS should still be stored");
}

static void test_reversed_non_zero_lhs_no_derivation() {
    TEST_CASE("Property 12: Non-zero LHS does NOT trigger reversed pattern");

    RelationStore rs;
    PropertyStore ps;

    // LHS is 5, not 0
    auto five = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(5)));
    SymbolicExpr x = make_var("x");

    rs.add_relation(five, x, RelationalNode::Op::LT, ps);

    // x should NOT have Positive derived (only 0 op var triggers)
    EXPECT_FALSE(ps.has_sign("x", Sign::Positive),
        "5 < x should not trigger reversed pattern (non-zero LHS)");
}

static void test_reversed_relation_stored() {
    TEST_CASE("Property 12: Reversed relations are stored in the RelationStore");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    rs.add_relation(zero, x, RelationalNode::Op::LT, ps);

    EXPECT_TRUE(rs.has_relation(zero, x, RelationalNode::Op::LT),
        "0 < x relation should be stored");
}

static void test_reversed_eq_no_sign_derivation() {
    TEST_CASE("Property 12: 0 EQ var does NOT derive sign (EQ not mapped)");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    rs.add_relation(zero, x, RelationalNode::Op::EQ, ps);

    // EQ is not mapped to any sign in the reversed pattern
    EXPECT_FALSE(ps.has_sign("x", Sign::Positive), "0 == x does not derive Positive");
    EXPECT_FALSE(ps.has_sign("x", Sign::Negative), "0 == x does not derive Negative");
    EXPECT_FALSE(ps.has_sign("x", Sign::NonZero), "0 == x does not derive NonZero");
}


int main() {
    test_transitive_gt_gt_chain();
    test_transitive_geq_gt_chain();
    test_transitive_gt_geq_chain();
    test_transitive_geq_geq_chain();
    test_transitive_longer_chain();
    test_transitive_mixed_chain();
    test_transitive_backward_chaining();
    test_transitive_no_duplicate_deduction();
    test_transitive_lt_leq_not_transitive();
    test_transitive_cap_at_64();
    test_transitive_sign_derivation_from_chain();

    test_reversed_0_lt_var_positive();
    test_reversed_0_gt_var_negative();
    test_reversed_0_geq_var_nonpositive();
    test_reversed_0_leq_var_nonnegative();
    test_reversed_0_neq_var_nonzero();
    test_reversed_all_operators_comprehensive();
    test_reversed_multiple_variables();
    test_reversed_non_variable_rhs_no_derivation();
    test_reversed_non_zero_lhs_no_derivation();
    test_reversed_relation_stored();
    test_reversed_eq_no_sign_derivation();

    return TEST_REPORT();
}
