/**
 * @file test_assumption_relation_store.cpp
 * @brief Property tests for RelationStore — Property 9: Relation storage with sign property derivation.
 *
 * Feature: assumption-system, Property 9: Relation storage with sign property derivation
 *
 * For any variable name and comparison operator (GT, GEQ, LT, LEQ, NEQ) against zero,
 * storing that relation should mark the variable with the corresponding sign property
 * (Positive, NonNegative, Negative, NonPositive, NonZero respectively).
 *
 * Validates: Requirements 4.1, 4.2, 4.3, 4.4, 4.5, 4.6
 */

#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

using namespace lamina;

/// Helper: create a SymbolicExpr wrapping a VariableNode with the given name.
static SymbolicExpr make_var_expr(const std::string& name) {
    return SymbolicExpr(std::make_shared<VariableNode>(name));
}

/// Helper: create a SymbolicExpr wrapping a NumberNode with value 0.
static SymbolicExpr make_zero_expr() {
    return SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));
}

/// Helper: check that a PropertyStore has a specific sign for a symbol.
static bool check_has_sign(const PropertyStore& ps, const std::string& symbol, Sign sign) {
    return ps.has_sign(symbol, sign);
}

// ============================================================================
// Property 9: Relation storage with sign property derivation
// ============================================================================

void test_gt_zero_derives_positive() {
    TEST_CASE("Property 9: GT against zero derives Positive sign");

    // Test with multiple variable names
    std::vector<std::string> var_names = {"x", "y", "alpha", "longVariableName", "a1"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationalNode::Op::GT, ps);

        // The variable should now have Positive sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::Positive),
                    name + " > 0 should derive Positive sign");

        // Positive implies NonNegative and NonZero
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonNegative),
                    name + " > 0 should imply NonNegative");
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " > 0 should imply NonZero");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::GT),
                    name + " > 0 relation should be stored");
    }
}

void test_geq_zero_derives_nonnegative() {
    TEST_CASE("Property 9: GEQ against zero derives NonNegative sign");

    std::vector<std::string> var_names = {"x", "beta", "var_2", "Z", "temp"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationalNode::Op::GEQ, ps);

        // The variable should now have NonNegative sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonNegative),
                    name + " >= 0 should derive NonNegative sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::GEQ),
                    name + " >= 0 relation should be stored");
    }
}

void test_lt_zero_derives_negative() {
    TEST_CASE("Property 9: LT against zero derives Negative sign");

    std::vector<std::string> var_names = {"x", "gamma", "n", "val", "q"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationalNode::Op::LT, ps);

        // The variable should now have Negative sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::Negative),
                    name + " < 0 should derive Negative sign");

        // Negative implies NonPositive and NonZero
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonPositive),
                    name + " < 0 should imply NonPositive");
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " < 0 should imply NonZero");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::LT),
                    name + " < 0 relation should be stored");
    }
}

void test_leq_zero_derives_nonpositive() {
    TEST_CASE("Property 9: LEQ against zero derives NonPositive sign");

    std::vector<std::string> var_names = {"x", "delta", "m", "result", "w"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationalNode::Op::LEQ, ps);

        // The variable should now have NonPositive sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonPositive),
                    name + " <= 0 should derive NonPositive sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::LEQ),
                    name + " <= 0 relation should be stored");
    }
}

void test_neq_zero_derives_nonzero() {
    TEST_CASE("Property 9: NEQ against zero derives NonZero sign");

    std::vector<std::string> var_names = {"x", "epsilon", "k", "divisor", "p"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationalNode::Op::NEQ, ps);

        // The variable should now have NonZero sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " != 0 should derive NonZero sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::NEQ),
                    name + " != 0 relation should be stored");
    }
}

void test_all_operators_comprehensive() {
    TEST_CASE("Property 9: All operators mapped correctly for a single variable");

    // Test all 5 operators on the same variable name (each in a fresh store)
    struct TestCase {
        RelationalNode::Op op;
        Sign expected_sign;
        std::string op_str;
    };

    std::vector<TestCase> cases = {
        {RelationalNode::Op::GT,  Sign::Positive,    "GT"},
        {RelationalNode::Op::GEQ, Sign::NonNegative, "GEQ"},
        {RelationalNode::Op::LT,  Sign::Negative,    "LT"},
        {RelationalNode::Op::LEQ, Sign::NonPositive,  "LEQ"},
        {RelationalNode::Op::NEQ, Sign::NonZero,     "NEQ"},
    };

    for (const auto& tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr("x");
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, tc.op, ps);

        EXPECT_TRUE(check_has_sign(ps, "x", tc.expected_sign),
                    "x " + tc.op_str + " 0 should derive expected sign");
    }
}

void test_composite_relation_no_sign_derivation() {
    TEST_CASE("Property 9: Composite LHS (non-variable) does not derive sign");

    RelationStore rs;
    PropertyStore ps;

    // Create a composite expression: x + y (AddNode, not a single VariableNode)
    auto x_node = std::make_shared<VariableNode>("x");
    auto y_node = std::make_shared<VariableNode>("y");
    auto add_node = std::make_shared<AddNode>(
        std::vector<std::shared_ptr<SymbolicNode>>{x_node, y_node});
    SymbolicExpr composite_expr(add_node);
    SymbolicExpr zero_expr = make_zero_expr();

    rs.add_relation(composite_expr, zero_expr, RelationalNode::Op::GT, ps);

    // Neither x nor y should have sign derived (composite LHS)
    EXPECT_FALSE(check_has_sign(ps, "x", Sign::Positive),
                 "Composite LHS should not derive sign for x");
    EXPECT_FALSE(check_has_sign(ps, "y", Sign::Positive),
                 "Composite LHS should not derive sign for y");

    // But the relation should still be stored
    EXPECT_TRUE(rs.has_relation(composite_expr, zero_expr, RelationalNode::Op::GT),
                "Composite relation should still be stored");
}

void test_nonzero_rhs_no_sign_derivation() {
    TEST_CASE("Property 9: Non-zero RHS does not derive sign property");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    // RHS is 5, not 0
    SymbolicExpr five_expr(std::make_shared<NumberNode>(BigInt(5)));

    rs.add_relation(var_expr, five_expr, RelationalNode::Op::GT, ps);

    // x > 5 should NOT derive Positive sign (only x > 0 pattern triggers derivation)
    EXPECT_FALSE(check_has_sign(ps, "x", Sign::Positive),
                 "x > 5 should not derive Positive sign (non-zero RHS)");

    // But the relation should still be stored
    EXPECT_TRUE(rs.has_relation(var_expr, five_expr, RelationalNode::Op::GT),
                "Non-zero RHS relation should still be stored");
}

void test_relation_stored_regardless_of_pattern() {
    TEST_CASE("Property 9: Relations are always stored regardless of pattern");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    // Add multiple relations
    rs.add_relation(var_expr, zero_expr, RelationalNode::Op::GT, ps);

    const auto& relations = rs.get_relations();
    EXPECT_TRUE(relations.size() == 1, "Should have 1 stored relation");

    // Add another relation
    SymbolicExpr y_expr = make_var_expr("y");
    rs.add_relation(y_expr, zero_expr, RelationalNode::Op::LT, ps);

    EXPECT_TRUE(rs.get_relations().size() == 2, "Should have 2 stored relations");
}

void test_clear_removes_all_relations() {
    TEST_CASE("Property 9: clear() removes all stored relations");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    rs.add_relation(var_expr, zero_expr, RelationalNode::Op::GT, ps);
    rs.add_relation(make_var_expr("y"), zero_expr, RelationalNode::Op::LT, ps);

    EXPECT_TRUE(rs.get_relations().size() == 2, "Should have 2 relations before clear");

    rs.clear();

    EXPECT_TRUE(rs.get_relations().empty(), "Should have 0 relations after clear");
    EXPECT_FALSE(rs.has_relation(var_expr, zero_expr, RelationalNode::Op::GT),
                 "has_relation should return false after clear");
}

int main() {
    test_gt_zero_derives_positive();
    test_geq_zero_derives_nonnegative();
    test_lt_zero_derives_negative();
    test_leq_zero_derives_nonpositive();
    test_neq_zero_derives_nonzero();
    test_all_operators_comprehensive();
    test_composite_relation_no_sign_derivation();
    test_nonzero_rhs_no_sign_derivation();
    test_relation_stored_regardless_of_pattern();
    test_clear_removes_all_relations();

    return TEST_REPORT();
}
