
#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <string>
#include <vector>

using namespace LMCAS;

/// Create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string &name) {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
}

/// Create a SymbolicExpr wrapping a NumberNode with value 0.
static SymbolicExpr make_zero() {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
}

TEST(AssumptionRelationProps, TransitiveGtGtChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x > y
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    // Add y > z — should trigger transitive closure: x > z
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x > y, y > z => x > z (GT+GT => GT)";
}

TEST(AssumptionRelationProps, TransitiveGeqGtChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x >= y
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    // Add y > z — should trigger: x > z (GEQ+GT => GT)
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x >= y, y > z => x > z (GEQ+GT => GT)";
}

TEST(AssumptionRelationProps, TransitiveGtGeqChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x > y
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    // Add y >= z — should trigger: x > z (GT+GEQ => GT)
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x > y, y >= z => x > z (GT+GEQ => GT)";
}

TEST(AssumptionRelationProps, TransitiveGeqGeqChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add x >= y
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    // Add y >= z — should trigger: x >= z (GEQ+GEQ => GEQ)
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GEQ))) << "x >= y, y >= z => x >= z (GEQ+GEQ => GEQ)";
}

TEST(AssumptionRelationProps, TransitiveLongerChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    EXPECT_TRUE((rs.add_relation(a, b, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(b, c, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(c, d, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // After adding a>b, b>c: a>c should be deduced
    EXPECT_TRUE((rs.has_relation(a, c, RelationalNode::Op::GT))) << "a > b, b > c => a > c";

    // After adding c>d: a>d, b>d should be deduced
    EXPECT_TRUE((rs.has_relation(a, d, RelationalNode::Op::GT))) << "a > b > c > d => a > d";
    EXPECT_TRUE((rs.has_relation(b, d, RelationalNode::Op::GT))) << "b > c > d => b > d";
}

TEST(AssumptionRelationProps, TransitiveMixedChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    EXPECT_TRUE((rs.add_relation(a, b, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(b, c, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(c, d, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // a > b >= c > d => a > d (GT+GEQ => GT, then GT+GT => GT)
    EXPECT_TRUE((rs.has_relation(a, d, RelationalNode::Op::GT))) << "a > b >= c > d => a > d";
}

TEST(AssumptionRelationProps, TransitiveBackwardChaining) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // Add y > z first
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    // Then add x > y — should trigger backward chain: x > z
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "y > z then x > y => x > z (backward chaining)";
}

TEST(AssumptionRelationProps, TransitiveNoDuplicateDeduction) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // x > z should be deduced once
    size_t count = 0;
    for (const auto &rel : rs.get_relations()) {
        if (!LMCAS::detail::node(rel.lhs) || !LMCAS::detail::node(rel.rhs)) {
            continue;
        }
        auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.lhs));
        auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.rhs));
        if (!lhs_var || !rhs_var) {
            continue;
        }
        if (lhs_var->name() == "x" && rhs_var->name() == "z" &&
            rel.op == RelationalNode::Op::GT) {
            ++count;
        }
    }
    EXPECT_TRUE((count == 1)) << "x > z should appear exactly once (no duplicates)";
}

TEST(AssumptionRelationProps, TransitiveLtLeqNotTransitive) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // LT relations should be stored but not trigger transitive closure
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    // x < z should NOT be deduced (only GT/GEQ participate)
    EXPECT_FALSE((rs.has_relation(x, z, RelationalNode::Op::LT))) << "LT does not participate in transitive closure";
}

TEST(AssumptionRelationProps, TransitiveCapAt64) {
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
        EXPECT_TRUE((rs.add_relation(vars[i], vars[i + 1], RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    }

    // Count relations before adding the 70th variable link
    size_t before_count = rs.get_relations().size();

    // Now add v_all > v0 which would try to chain through all 69 existing links
    SymbolicExpr v_new = make_var("v_new");
    EXPECT_TRUE((rs.add_relation(v_new, vars[0], RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    size_t after_count = rs.get_relations().size();
    // The new relations added should be at most 64 + 1 (the original relation itself)
    size_t new_relations = after_count - before_count;

    // The directly added relation is 1, plus at most 64 deduced
    EXPECT_TRUE((new_relations <= 65)) << "At most 64 new deduced relations + 1 original per add_relation call (got " +
                                              std::to_string(new_relations) + ")";
}

TEST(AssumptionRelationProps, TransitiveSignDerivationFromChain) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr zero = make_zero();

    // x > y, y > 0 => x > 0 => x is Positive
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, zero, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // y > 0 directly derives Positive for y
    EXPECT_TRUE((ps.has_sign("y", Sign::Positive))) << "y > 0 derives Positive for y";

    // x > 0 should be deduced transitively, deriving Positive for x
    EXPECT_TRUE((rs.has_relation(x, zero, RelationalNode::Op::GT))) << "x > y > 0 => x > 0 deduced";
    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "x > y > 0 => x is Positive (sign derived from transitive closure)";
}

TEST(AssumptionRelationProps, Reversed0LtVarPositive) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    // 0 < x means x > 0, so x is Positive
    EXPECT_TRUE((rs.add_relation(zero, x, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "0 < x => x is Positive";
    EXPECT_TRUE((ps.has_sign("x", Sign::NonNegative))) << "0 < x => x is NonNegative (implied)";
    EXPECT_TRUE((ps.has_sign("x", Sign::NonZero))) << "0 < x => x is NonZero (implied)";
}

TEST(AssumptionRelationProps, Reversed0GtVarNegative) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr y = make_var("y");

    // 0 > y means y < 0, so y is Negative
    EXPECT_TRUE((rs.add_relation(zero, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("y", Sign::Negative))) << "0 > y => y is Negative";
    EXPECT_TRUE((ps.has_sign("y", Sign::NonPositive))) << "0 > y => y is NonPositive (implied)";
    EXPECT_TRUE((ps.has_sign("y", Sign::NonZero))) << "0 > y => y is NonZero (implied)";
}

TEST(AssumptionRelationProps, Reversed0GeqVarNonpositive) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr z = make_var("z");

    // 0 >= z means z <= 0, so z is NonPositive
    EXPECT_TRUE((rs.add_relation(zero, z, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("z", Sign::NonPositive))) << "0 >= z => z is NonPositive";
}

TEST(AssumptionRelationProps, Reversed0LeqVarNonnegative) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr w = make_var("w");

    // 0 <= w means w >= 0, so w is NonNegative
    EXPECT_TRUE((rs.add_relation(zero, w, RelationalNode::Op::LEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("w", Sign::NonNegative))) << "0 <= w => w is NonNegative";
}

TEST(AssumptionRelationProps, Reversed0NeqVarNonzero) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr v = make_var("v");

    // 0 != v means v != 0, so v is NonZero
    EXPECT_TRUE((rs.add_relation(zero, v, RelationalNode::Op::NEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("v", Sign::NonZero))) << "0 != v => v is NonZero";
}

TEST(AssumptionRelationProps, ReversedAllOperatorsComprehensive) {
    struct TestCase {
        RelationalNode::Op op;
        Sign expected_sign;
        std::string desc;
    };

    std::vector<TestCase> cases = {
        {RelationalNode::Op::LT, Sign::Positive, "0 LT var => Positive"},
        {RelationalNode::Op::GT, Sign::Negative, "0 GT var => Negative"},
        {RelationalNode::Op::GEQ, Sign::NonPositive, "0 GEQ var => NonPositive"},
        {RelationalNode::Op::LEQ, Sign::NonNegative, "0 LEQ var => NonNegative"},
        {RelationalNode::Op::NEQ, Sign::NonZero, "0 NEQ var => NonZero"},
    };

    for (const auto &tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr zero = make_zero();
        SymbolicExpr var = make_var("t");

        EXPECT_TRUE((rs.add_relation(zero, var, tc.op, ps).has_value())) << "relation insertion succeeds";

        EXPECT_TRUE((ps.has_sign("t", tc.expected_sign))) << tc.desc;
    }
}

TEST(AssumptionRelationProps, ReversedMultipleVariables) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();

    // 0 < a, 0 > b, 0 >= c, 0 <= d, 0 != e
    EXPECT_TRUE((rs.add_relation(zero, make_var("a"), RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(zero, make_var("b"), RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(zero, make_var("c"), RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(zero, make_var("d"), RelationalNode::Op::LEQ, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(zero, make_var("e"), RelationalNode::Op::NEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("a", Sign::Positive))) << "a is Positive";
    EXPECT_TRUE((ps.has_sign("b", Sign::Negative))) << "b is Negative";
    EXPECT_TRUE((ps.has_sign("c", Sign::NonPositive))) << "c is NonPositive";
    EXPECT_TRUE((ps.has_sign("d", Sign::NonNegative))) << "d is NonNegative";
    EXPECT_TRUE((ps.has_sign("e", Sign::NonZero))) << "e is NonZero";
}

TEST(AssumptionRelationProps, ReversedNonVariableRhsNoDerivation) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();

    // RHS is a composite expression (x + y), not a single variable
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto y_node = LMCAS::detail::make_node<VariableNode>("y");
    auto add_node = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, y_node});
    auto composite = LMCAS::detail::expression_from_node(add_node);
    EXPECT_TRUE((rs.add_relation(zero, composite, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    // Neither x nor y should have sign derived
    EXPECT_FALSE((ps.has_sign("x", Sign::Positive))) << "0 < (x+y) should not derive sign for x";
    EXPECT_FALSE((ps.has_sign("y", Sign::Positive))) << "0 < (x+y) should not derive sign for y";

    // But the relation should still be stored
    EXPECT_TRUE((rs.has_relation(zero, composite, RelationalNode::Op::LT))) << "Relation with composite RHS should still be stored";
}

TEST(AssumptionRelationProps, ReversedNonZeroLhsNoDerivation) {
    RelationStore rs;
    PropertyStore ps;

    // LHS is 5, not 0
    auto five = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(5)));
    SymbolicExpr x = make_var("x");

    EXPECT_TRUE((rs.add_relation(five, x, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    // x should NOT have Positive derived (only 0 op var triggers)
    EXPECT_FALSE((ps.has_sign("x", Sign::Positive))) << "5 < x should not trigger reversed pattern (non-zero LHS)";
}

TEST(AssumptionRelationProps, ReversedRelationStored) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    EXPECT_TRUE((rs.add_relation(zero, x, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.has_relation(zero, x, RelationalNode::Op::LT))) << "0 < x relation should be stored";
}

TEST(AssumptionRelationProps, ReversedEqNoSignDerivation) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr x = make_var("x");

    EXPECT_TRUE((rs.add_relation(zero, x, RelationalNode::Op::EQ, ps).has_value())) << "relation insertion succeeds";

    // EQ is not mapped to any sign in the reversed pattern
    EXPECT_FALSE((ps.has_sign("x", Sign::Positive))) << "0 == x does not derive Positive";
    EXPECT_FALSE((ps.has_sign("x", Sign::Negative))) << "0 == x does not derive Negative";
    EXPECT_FALSE((ps.has_sign("x", Sign::NonZero))) << "0 == x does not derive NonZero";
}
