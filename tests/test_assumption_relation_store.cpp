
#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"

using namespace lamina;

/// Helper: create a variable expression through the stable public factory.
static SymbolicExpr make_var_expr(const std::string& name) {
    return *SymbolicExpr::variable(name);
}

/// Helper: create zero through the stable public factory.
static SymbolicExpr make_zero_expr() {
    return *SymbolicExpr::number(0);
}

/// Helper: check that a PropertyStore has a specific sign for a symbol.
static bool check_has_sign(const PropertyStore& ps, const std::string& symbol, Sign sign) {
    return ps.has_sign(symbol, sign);
}


void test_gt_zero_derives_positive() {
    TEST_CASE("GT against zero derives Positive sign");

    // Test with multiple variable names
    std::vector<std::string> var_names = {"x", "y", "alpha", "longVariableName", "a1"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps);

        // The variable should now have Positive sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::Positive),
                    name + " > 0 should derive Positive sign");

        // Positive implies NonNegative and NonZero
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonNegative),
                    name + " > 0 should imply NonNegative");
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " > 0 should imply NonZero");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationOp::GT),
                    name + " > 0 relation should be stored");
    }
}

void test_geq_zero_derives_nonnegative() {
    TEST_CASE("GEQ against zero derives NonNegative sign");

    std::vector<std::string> var_names = {"x", "beta", "var_2", "Z", "temp"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationOp::GEQ, ps);

        // The variable should now have NonNegative sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonNegative),
                    name + " >= 0 should derive NonNegative sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationOp::GEQ),
                    name + " >= 0 relation should be stored");
    }
}

void test_lt_zero_derives_negative() {
    TEST_CASE("LT against zero derives Negative sign");

    std::vector<std::string> var_names = {"x", "gamma", "n", "val", "q"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationOp::LT, ps);

        // The variable should now have Negative sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::Negative),
                    name + " < 0 should derive Negative sign");

        // Negative implies NonPositive and NonZero
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonPositive),
                    name + " < 0 should imply NonPositive");
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " < 0 should imply NonZero");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationOp::LT),
                    name + " < 0 relation should be stored");
    }
}

void test_leq_zero_derives_nonpositive() {
    TEST_CASE("LEQ against zero derives NonPositive sign");

    std::vector<std::string> var_names = {"x", "delta", "m", "result", "w"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationOp::LEQ, ps);

        // The variable should now have NonPositive sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonPositive),
                    name + " <= 0 should derive NonPositive sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationOp::LEQ),
                    name + " <= 0 relation should be stored");
    }
}

void test_neq_zero_derives_nonzero() {
    TEST_CASE("NEQ against zero derives NonZero sign");

    std::vector<std::string> var_names = {"x", "epsilon", "k", "divisor", "p"};

    for (const auto& name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        rs.add_relation(var_expr, zero_expr, RelationOp::NEQ, ps);

        // The variable should now have NonZero sign
        EXPECT_TRUE(check_has_sign(ps, name, Sign::NonZero),
                    name + " != 0 should derive NonZero sign");

        // The relation should be stored
        EXPECT_TRUE(rs.has_relation(var_expr, zero_expr, RelationOp::NEQ),
                    name + " != 0 relation should be stored");
    }
}

void test_all_operators_comprehensive() {
    TEST_CASE("All operators mapped correctly for a single variable");

    // Test all 5 operators on the same variable name (each in a fresh store)
    struct TestCase {
        RelationOp op;
        Sign expected_sign;
        std::string op_str;
    };

    std::vector<TestCase> cases = {
        {RelationOp::GT,  Sign::Positive,    "GT"},
        {RelationOp::GEQ, Sign::NonNegative, "GEQ"},
        {RelationOp::LT,  Sign::Negative,    "LT"},
        {RelationOp::LEQ, Sign::NonPositive,  "LEQ"},
        {RelationOp::NEQ, Sign::NonZero,     "NEQ"},
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
    TEST_CASE("Composite LHS (non-variable) does not derive sign");

    RelationStore rs;
    PropertyStore ps;

    /// 创建复合表达式，覆盖多变量关系存储路径。
    auto composite_expr = *SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::variable("y"));
    SymbolicExpr zero_expr = make_zero_expr();

    rs.add_relation(composite_expr, zero_expr, RelationOp::GT, ps);

    // Neither x nor y should have sign derived (composite LHS)
    EXPECT_FALSE(check_has_sign(ps, "x", Sign::Positive),
                 "Composite LHS should not derive sign for x");
    EXPECT_FALSE(check_has_sign(ps, "y", Sign::Positive),
                 "Composite LHS should not derive sign for y");

    // But the relation should still be stored
    EXPECT_TRUE(rs.has_relation(composite_expr, zero_expr, RelationOp::GT),
                "Composite relation should still be stored");
}

void test_nonzero_rhs_no_sign_derivation() {
    TEST_CASE("Non-zero RHS does not derive sign property");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    // RHS is 5, not 0
    auto five_expr = *SymbolicExpr::number(5);
    rs.add_relation(var_expr, five_expr, RelationOp::GT, ps);

    // x > 5 should NOT derive Positive sign (only x > 0 pattern triggers derivation)
    EXPECT_FALSE(check_has_sign(ps, "x", Sign::Positive),
                 "x > 5 should not derive Positive sign (non-zero RHS)");

    // But the relation should still be stored
    EXPECT_TRUE(rs.has_relation(var_expr, five_expr, RelationOp::GT),
                "Non-zero RHS relation should still be stored");
}

void test_relation_stored_regardless_of_pattern() {
    TEST_CASE("Relations are always stored regardless of pattern");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    // Add multiple relations
    rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps);

    const auto& relations = rs.get_relations();
    EXPECT_TRUE(relations.size() == 1, "Should have 1 stored relation");

    // Add another relation
    SymbolicExpr y_expr = make_var_expr("y");
    rs.add_relation(y_expr, zero_expr, RelationOp::LT, ps);

    EXPECT_TRUE(rs.get_relations().size() == 2, "Should have 2 stored relations");
}

void test_clear_removes_all_relations() {
    TEST_CASE("clear() removes all stored relations");

    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps);
    rs.add_relation(make_var_expr("y"), zero_expr, RelationOp::LT, ps);

    EXPECT_TRUE(rs.get_relations().size() == 2, "Should have 2 relations before clear");

    rs.clear();

    EXPECT_TRUE(rs.get_relations().empty(), "Should have 0 relations after clear");
    EXPECT_FALSE(rs.has_relation(var_expr, zero_expr, RelationOp::GT),
                 "has_relation should return false after clear");
}

void test_checked_add_relation_contracts() {
    TEST_CASE("RelationStore checked add_relation: errors and transactional state");

    RelationStore rs;
    PropertyStore ps;
    SymbolicExpr x = make_var_expr("x");
    SymbolicExpr zero = make_zero_expr();

    auto success = rs.add_relation_checked(x, zero, RelationOp::GT, ps);
    EXPECT_TRUE(success.has_value(), "checked add_relation succeeds for x > 0");
    EXPECT_TRUE(rs.has_relation(x, zero, RelationOp::GT),
                "checked add_relation stores successful relation");
    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
                "checked add_relation commits derived property");

    const auto relation_count = rs.get_relations().size();
    auto conflict = rs.add_relation_checked(x, zero, RelationOp::LT, ps);
    EXPECT_TRUE(!conflict.has_value(),
                "checked add_relation rejects property contradiction");
    EXPECT_TRUE(conflict.error().code == CasErrc::InvalidArgument,
                "checked add_relation reports InvalidArgument for contradiction");
    EXPECT_TRUE(rs.get_relations().size() == relation_count,
                "failed checked add_relation preserves relation count");
    EXPECT_FALSE(rs.has_relation(x, zero, RelationOp::LT),
                 "failed checked add_relation does not store conflicting relation");
    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
                "failed checked add_relation preserves previous property");
    EXPECT_FALSE(ps.has_sign("x", Sign::Negative),
                 "failed checked add_relation does not apply conflicting property");
}

void test_legacy_add_relation_is_transactional() {
    TEST_CASE("RelationStore canonical add_relation delegates transactionally");

    RelationStore rs;
    PropertyStore ps;
    SymbolicExpr x = make_var_expr("x");
    SymbolicExpr zero = make_zero_expr();

    rs.add_relation(x, zero, RelationOp::GT, ps);
    const auto relation_count = rs.get_relations().size();

    auto failure_318 = rs.add_relation(x, zero, RelationOp::LT, ps);
    EXPECT_TRUE(!failure_318.has_value(), "canonical add_relation maps checked contradiction to invalid_argument");
    EXPECT_TRUE(rs.get_relations().size() == relation_count,
                "canonical add_relation does not retain a failed relation");
    EXPECT_FALSE(rs.has_relation(x, zero, RelationOp::LT),
                 "canonical add_relation preserves transactional relation state");
    EXPECT_TRUE(ps.has_sign("x", Sign::Positive),
                "canonical add_relation preserves the previously proven sign");
    EXPECT_FALSE(ps.has_sign("x", Sign::Negative),
                 "canonical add_relation does not commit a contradictory sign");
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
    test_checked_add_relation_contracts();
    test_legacy_add_relation_is_transactional();

    return TEST_REPORT();
}
