
#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "assumption_context.hpp"
#include <utility>
#include "symbolic.hpp"

using namespace LMCAS;

/// Helper: create a variable expression through the stable public factory.
static SymbolicExpr make_var_expr(const std::string &name) {
    return *SymbolicExpr::variable(name);
}

/// Helper: create zero through the stable public factory.
static SymbolicExpr make_zero_expr() {
    return *SymbolicExpr::number(0);
}

TEST(AssumptionRelationStore, GtZeroDerivesPositive) {
    // Test with multiple variable names
    std::vector<std::string> var_names = {"x", "y", "alpha", "longVariableName", "a1"};

    for (const auto &name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";

        // The variable should now have Positive sign
        EXPECT_TRUE((ps.has_sign(name, Sign::Positive))) << name + " > 0 should derive Positive sign";

        // Positive implies NonNegative and NonZero
        EXPECT_TRUE((ps.has_sign(name, Sign::NonNegative))) << name + " > 0 should imply NonNegative";
        EXPECT_TRUE((ps.has_sign(name, Sign::NonZero))) << name + " > 0 should imply NonZero";

        // The relation should be stored
        EXPECT_TRUE((rs.has_relation(var_expr, zero_expr, RelationOp::GT))) << name + " > 0 relation should be stored";
    }
}

TEST(AssumptionRelationStore, GeqZeroDerivesNonnegative) {
    std::vector<std::string> var_names = {"x", "beta", "var_2", "Z", "temp"};

    for (const auto &name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::GEQ, ps).has_value())) << "relation insertion succeeds";

        // The variable should now have NonNegative sign
        EXPECT_TRUE((ps.has_sign(name, Sign::NonNegative))) << name + " >= 0 should derive NonNegative sign";

        // The relation should be stored
        EXPECT_TRUE((rs.has_relation(var_expr, zero_expr, RelationOp::GEQ))) << name + " >= 0 relation should be stored";
    }
}

TEST(AssumptionRelationStore, LtZeroDerivesNegative) {
    std::vector<std::string> var_names = {"x", "gamma", "n", "val", "q"};

    for (const auto &name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::LT, ps).has_value())) << "relation insertion succeeds";

        // The variable should now have Negative sign
        EXPECT_TRUE((ps.has_sign(name, Sign::Negative))) << name + " < 0 should derive Negative sign";

        // Negative implies NonPositive and NonZero
        EXPECT_TRUE((ps.has_sign(name, Sign::NonPositive))) << name + " < 0 should imply NonPositive";
        EXPECT_TRUE((ps.has_sign(name, Sign::NonZero))) << name + " < 0 should imply NonZero";

        // The relation should be stored
        EXPECT_TRUE((rs.has_relation(var_expr, zero_expr, RelationOp::LT))) << name + " < 0 relation should be stored";
    }
}

TEST(AssumptionRelationStore, LeqZeroDerivesNonpositive) {
    std::vector<std::string> var_names = {"x", "delta", "m", "result", "w"};

    for (const auto &name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::LEQ, ps).has_value())) << "relation insertion succeeds";

        // The variable should now have NonPositive sign
        EXPECT_TRUE((ps.has_sign(name, Sign::NonPositive))) << name + " <= 0 should derive NonPositive sign";

        // The relation should be stored
        EXPECT_TRUE((rs.has_relation(var_expr, zero_expr, RelationOp::LEQ))) << name + " <= 0 relation should be stored";
    }
}

TEST(AssumptionRelationStore, NeqZeroDerivesNonzero) {
    std::vector<std::string> var_names = {"x", "epsilon", "k", "divisor", "p"};

    for (const auto &name : var_names) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr(name);
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::NEQ, ps).has_value())) << "relation insertion succeeds";

        // The variable should now have NonZero sign
        EXPECT_TRUE((ps.has_sign(name, Sign::NonZero))) << name + " != 0 should derive NonZero sign";

        // The relation should be stored
        EXPECT_TRUE((rs.has_relation(var_expr, zero_expr, RelationOp::NEQ))) << name + " != 0 relation should be stored";
    }
}

TEST(AssumptionRelationStore, AllOperatorsComprehensive) {
    // Test all 5 operators on the same variable name (each in a fresh store)
    struct TestCase {
        RelationOp op;
        Sign expected_sign;
        std::string op_str;
    };

    std::vector<TestCase> cases = {
        {RelationOp::GT, Sign::Positive, "GT"},
        {RelationOp::GEQ, Sign::NonNegative, "GEQ"},
        {RelationOp::LT, Sign::Negative, "LT"},
        {RelationOp::LEQ, Sign::NonPositive, "LEQ"},
        {RelationOp::NEQ, Sign::NonZero, "NEQ"},
    };

    for (const auto &tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr var_expr = make_var_expr("x");
        SymbolicExpr zero_expr = make_zero_expr();

        EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, tc.op, ps).has_value())) << "relation insertion succeeds";

        EXPECT_TRUE((ps.has_sign("x", tc.expected_sign))) << "x " + tc.op_str + " 0 should derive expected sign";
    }
}

TEST(AssumptionRelationStore, CompositeRelationNoSignDerivation) {
    RelationStore rs;
    PropertyStore ps;

    /// 创建复合表达式,覆盖多变量关系存储路径.
    auto composite_expr = *SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::variable("y"));
    SymbolicExpr zero_expr = make_zero_expr();

    EXPECT_TRUE((rs.add_relation(composite_expr, zero_expr, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";

    // Neither x nor y should have sign derived (composite LHS)
    EXPECT_FALSE((ps.has_sign("x", Sign::Positive))) << "Composite LHS should not derive sign for x";
    EXPECT_FALSE((ps.has_sign("y", Sign::Positive))) << "Composite LHS should not derive sign for y";

    // But the relation should still be stored
    EXPECT_TRUE((rs.has_relation(composite_expr, zero_expr, RelationOp::GT))) << "Composite relation should still be stored";
}

TEST(AssumptionRelationStore, NonzeroRhsNoSignDerivation) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    // RHS is 5, not 0
    auto five_expr = *SymbolicExpr::number(5);
    EXPECT_TRUE((rs.add_relation(var_expr, five_expr, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";

    // x > 5 should NOT derive Positive sign (only x > 0 pattern triggers derivation)
    EXPECT_FALSE((ps.has_sign("x", Sign::Positive))) << "x > 5 should not derive Positive sign (non-zero RHS)";

    // But the relation should still be stored
    EXPECT_TRUE((rs.has_relation(var_expr, five_expr, RelationOp::GT))) << "Non-zero RHS relation should still be stored";
}

TEST(AssumptionRelationStore, RelationStoredRegardlessOfPattern) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    // Add multiple relations
    EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";

    const auto &relations = rs.get_relations();
    EXPECT_TRUE((relations.size() == 1)) << "Should have 1 stored relation";

    // Add another relation
    SymbolicExpr y_expr = make_var_expr("y");
    EXPECT_TRUE((rs.add_relation(y_expr, zero_expr, RelationOp::LT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.get_relations().size() == 2)) << "Should have 2 stored relations";
}

TEST(AssumptionRelationStore, ClearRemovesAllRelations) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr var_expr = make_var_expr("x");
    SymbolicExpr zero_expr = make_zero_expr();

    EXPECT_TRUE((rs.add_relation(var_expr, zero_expr, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(make_var_expr("y"), zero_expr, RelationOp::LT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((rs.get_relations().size() == 2)) << "Should have 2 relations before clear";

    rs.clear();

    EXPECT_TRUE((rs.get_relations().empty())) << "Should have 0 relations after clear";
    EXPECT_FALSE((rs.has_relation(var_expr, zero_expr, RelationOp::GT))) << "has_relation should return false after clear";
}

TEST(AssumptionRelationStore, CheckedAddRelationContracts) {
    RelationStore rs;
    PropertyStore ps;
    SymbolicExpr x = make_var_expr("x");
    SymbolicExpr zero = make_zero_expr();

    auto success = rs.add_relation_checked(x, zero, RelationOp::GT, ps);
    EXPECT_TRUE((success.has_value())) << "checked add_relation succeeds for x > 0";
    EXPECT_TRUE((rs.has_relation(x, zero, RelationOp::GT))) << "checked add_relation stores successful relation";
    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "checked add_relation commits derived property";

    const auto relation_count = rs.get_relations().size();
    auto conflict = rs.add_relation_checked(x, zero, RelationOp::LT, ps);
    EXPECT_TRUE((!conflict.has_value())) << "checked add_relation rejects property contradiction";
    EXPECT_TRUE((conflict.error().code == CasErrc::InvalidArgument)) << "checked add_relation reports InvalidArgument for contradiction";
    EXPECT_TRUE((rs.get_relations().size() == relation_count)) << "failed checked add_relation preserves relation count";
    EXPECT_FALSE((rs.has_relation(x, zero, RelationOp::LT))) << "failed checked add_relation does not store conflicting relation";
    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "failed checked add_relation preserves previous property";
    EXPECT_FALSE((ps.has_sign("x", Sign::Negative))) << "failed checked add_relation does not apply conflicting property";
}

TEST(AssumptionRelationStore, LegacyAddRelationIsTransactional) {
    RelationStore rs;
    PropertyStore ps;
    SymbolicExpr x = make_var_expr("x");
    SymbolicExpr zero = make_zero_expr();

    EXPECT_TRUE((rs.add_relation(x, zero, RelationOp::GT, ps).has_value())) << "relation insertion succeeds";
    const auto relation_count = rs.get_relations().size();

    auto failure_318 = rs.add_relation(x, zero, RelationOp::LT, ps);
    EXPECT_TRUE((!failure_318.has_value())) << "canonical add_relation maps checked contradiction to invalid_argument";
    EXPECT_TRUE((rs.get_relations().size() == relation_count)) << "canonical add_relation does not retain a failed relation";
    EXPECT_FALSE((rs.has_relation(x, zero, RelationOp::LT))) << "canonical add_relation preserves transactional relation state";
    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "canonical add_relation preserves the previously proven sign";
    EXPECT_FALSE((ps.has_sign("x", Sign::Negative))) << "canonical add_relation does not commit a contradictory sign";
}

TEST(AssumptionRelationStore, ExplicitAndDerivedConclusionsShareOneEntry) {
    AssumptionContext ctx;
    const auto a = make_var_expr("a");
    const auto b = make_var_expr("b");
    const auto c = make_var_expr("c");
    auto& relations = ctx.current_relations();
    auto& properties = ctx.current_properties();

    ASSERT_TRUE(relations.add_relation(a, b, RelationOp::GT, properties));
    ASSERT_TRUE(relations.add_relation(b, c, RelationOp::GT, properties));
    ASSERT_TRUE(relations.has_relation(a, c, RelationOp::GT));
    ASSERT_EQ(relations.get_relations().size(), 3u);
    ASSERT_TRUE(relations.add_relation(a, c, RelationOp::GT, properties));
    ASSERT_TRUE(relations.add_relation(make_var_expr("a"), make_var_expr("c"),
                                       RelationOp::GT, properties));
    EXPECT_EQ(relations.get_relations().size(), 3u);

    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    EXPECT_TRUE(ctx.has_relation(a, c, RelationOp::GT));
    ASSERT_TRUE(ctx.pop());
    ctx.current_relations().clear();

    auto& rebuilt = ctx.current_relations();
    auto& rebuilt_properties = ctx.current_properties();
    ASSERT_TRUE(rebuilt.add_relation(a, c, RelationOp::GT, rebuilt_properties));
    ASSERT_TRUE(rebuilt.add_relation(a, b, RelationOp::GT, rebuilt_properties));
    ASSERT_TRUE(rebuilt.add_relation(b, c, RelationOp::GT, rebuilt_properties));
    EXPECT_EQ(rebuilt.get_relations().size(), 3u);
}

TEST(AssumptionRelationStore, RepeatedDeclarationStillDerivesAndChecksProperties) {
    RelationStore relations;
    PropertyStore properties;
    const auto x = make_var_expr("x");
    const auto zero = make_zero_expr();
    ASSERT_TRUE(relations.add_relation(x, zero, RelationOp::GT, properties));

    properties = PropertyStore{};
    ASSERT_TRUE(relations.add_relation(x, zero, RelationOp::GT, properties));
    EXPECT_TRUE(properties.has_sign("x", Sign::Positive));
    EXPECT_EQ(relations.get_relations().size(), 1u);

    properties = PropertyStore{};
    ASSERT_TRUE(properties.declare_sign("x", Sign::Negative));
    const auto revision = relations.revision();
    const auto rejected = relations.add_relation_checked(x, zero, RelationOp::GT, properties);
    ASSERT_FALSE(rejected);
    EXPECT_EQ(rejected.error().code, CasErrc::InvalidArgument);
    EXPECT_EQ(relations.revision(), revision);
    EXPECT_EQ(relations.get_relations().size(), 1u);
    EXPECT_TRUE(properties.has_sign("x", Sign::Negative));
    EXPECT_FALSE(properties.has_sign("x", Sign::Positive));
}

TEST(AssumptionRelationStore, CopyMoveAndClearPreserveAlternativeProofBehavior) {
    RelationStore original;
    PropertyStore properties;
    const auto a = make_var_expr("a");
    const auto b = make_var_expr("b");
    const auto c = make_var_expr("c");
    const auto d = make_var_expr("d");
    ASSERT_TRUE(original.add_relation(a, b, RelationOp::GT, properties));
    ASSERT_TRUE(original.add_relation(b, c, RelationOp::GT, properties));
    ASSERT_TRUE(original.add_relation(a, d, RelationOp::GT, properties));
    ASSERT_TRUE(original.add_relation(d, c, RelationOp::GT, properties));
    ASSERT_EQ(original.get_relations().size(), 5u);

    RelationStore copied(original);
    original.clear();
    AssumptionContext ctx;
    ctx.current_relations() = copied;
    copied.clear();
    const auto check_alternatives = [&]() {
        ctx.push();
        ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
        EXPECT_TRUE(ctx.has_relation(a, c, RelationOp::GT));
        ASSERT_TRUE(ctx.assume_sign("d", Sign::Negative));
        EXPECT_FALSE(ctx.has_relation(a, c, RelationOp::GT));
        ASSERT_TRUE(ctx.pop());
        EXPECT_TRUE(ctx.has_relation(a, c, RelationOp::GT));
    };
    check_alternatives();

    RelationStore moved(std::move(ctx.current_relations()));
    EXPECT_FALSE(ctx.has_relation(a, c, RelationOp::GT));
    ctx.current_relations() = std::move(moved);
    EXPECT_FALSE(moved.has_relation(a, c, RelationOp::GT));
    check_alternatives();

    ctx.current_relations().clear();
    EXPECT_FALSE(ctx.has_relation(a, c, RelationOp::GT));
    ASSERT_TRUE(ctx.current_relations().add_relation(a, b, RelationOp::GT,
                                                    ctx.current_properties()));
    ASSERT_TRUE(ctx.current_relations().add_relation(b, c, RelationOp::GT,
                                                    ctx.current_properties()));
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    EXPECT_FALSE(ctx.has_relation(a, c, RelationOp::GT));
}
