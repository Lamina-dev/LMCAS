
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace LMCAS;

/// Helper: create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var_expr(const std::string &name) {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
}

/// Helper: create a FunctionNode expression (e.g., ln(x), sqrt(x), exp(x)).
static SymbolicExpr make_func_expr(FunctionNode::FuncType type, const std::string &var_name) {
    auto var_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto func_node = LMCAS::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{var_node});
    return LMCAS::detail::expression_from_node(func_node);
}

/// Helper: create a PowerNode expression (var^n).
static SymbolicExpr make_power_expr(const std::string &var_name, int n) {
    auto var_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto exp_node = LMCAS::detail::make_node<NumberNode>(BigInt(n));
    auto pow_node = LMCAS::detail::make_node<PowerNode>(var_node, exp_node);
    return LMCAS::detail::expression_from_node(pow_node);
}

TEST(AssumptionMonotonicity, LnMonotonicityBothPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Add relation x > y
    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    // Apply monotonicity rules
    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that ln(x) > ln(y) was deduced
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "ln(x) > ln(y) should be deduced when both x,y are Positive and x > y";
}

TEST(AssumptionMonotonicity, CheckedMonotonicityRulesSuccess) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    auto seed = ctx.current_relations().add_relation_checked(
        x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());
    EXPECT_TRUE((seed.has_value())) << "checked seed relation succeeds";

    auto result = engine.apply_monotonicity_rules_checked(
        rel, ctx.current_relations(), ctx.current_properties());
    EXPECT_TRUE((result.has_value())) << "checked monotonicity rules report success";

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "checked monotonicity rules deduce ln(x) > ln(y)";
}

TEST(AssumptionMonotonicity, CheckedMonotonicityRulesNoop) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");
    Relation rel{x_expr, y_expr, RelationalNode::Op::LEQ};

    auto result = engine.apply_monotonicity_rules_checked(
        rel, ctx.current_relations(), ctx.current_properties());
    EXPECT_TRUE((result.has_value())) << "checked monotonicity no-op reports success";
    EXPECT_TRUE((ctx.current_relations().get_relations().empty())) << "checked monotonicity no-op leaves relation store unchanged";
}

TEST(AssumptionMonotonicity, SqrtMonotonicityBothPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that sqrt(x) > sqrt(y) was deduced
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT))) << "sqrt(x) > sqrt(y) should be deduced when both x,y are Positive and x > y";
}

TEST(AssumptionMonotonicity, ExpMonotonicityBothReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that exp(x) > exp(y) was deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT))) << "exp(x) > exp(y) should be deduced when both x,y are Real and x > y";
}

TEST(AssumptionMonotonicity, ExpMonotonicityPositiveImpliesReal) {
    AssumptionContext ctx;
    // Positive implies NonNegative and NonZero, but we also need Real domain
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that exp(x) > exp(y) was deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT))) << "exp(x) > exp(y) should be deduced when both are Positive+Real";
}

TEST(AssumptionMonotonicity, LnGuardMissingPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    // y has no sign assumption - guard should prevent ln rule

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // ln(x) > ln(y) should NOT be deduced
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_FALSE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "ln rule should NOT apply when y lacks Positive assumption";
}

TEST(AssumptionMonotonicity, SqrtGuardMissingPositive) {
    AssumptionContext ctx;
    // Neither x nor y has Positive assumption
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // sqrt(x) > sqrt(y) should NOT be deduced
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");

    EXPECT_FALSE((ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT))) << "sqrt rule should NOT apply when variables lack Positive assumption";
}

TEST(AssumptionMonotonicity, ExpGuardMissingReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    // y has no domain assumption (defaults to Complex)

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // exp(x) > exp(y) should NOT be deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_FALSE((ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT))) << "exp rule should NOT apply when y lacks Real domain";
}

TEST(AssumptionMonotonicity, NoRulesForNonGtRelation) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    /// 使用 LT 关系验证反向单调推导.
    Relation rel{x_expr, y_expr, RelationalNode::Op::LT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::LT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // No deduced relations should be added (only the original LT)
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_FALSE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "No monotonicity rules should apply for non-GT relations";
}

TEST(AssumptionMonotonicity, NoRulesForNonVariableOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    // Create a composite LHS: x + 1
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto one_node = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto add_node = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, one_node});
    auto composite_expr = LMCAS::detail::expression_from_node(add_node);
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{composite_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(composite_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // No deduced relations should be added for composite operands
    // (only 1 relation in the store: the original one)
    EXPECT_TRUE((ctx.current_relations().get_relations().size() == 1)) << "No monotonicity rules should apply for non-variable operands";
}

TEST(AssumptionMonotonicity, PowerMonotonicityBothNonnegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // First, add a relation that contains a power expression with x^2
    // so that the exponent 2 is "appearing in expressions"
    SymbolicExpr x_squared = make_power_expr("x", 2);
    auto zero_expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    EXPECT_TRUE((ctx.current_relations().add_relation(x_squared, zero_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    // Now add x > y and apply monotonicity
    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that x^2 > y^2 was deduced
    SymbolicExpr pow_x = make_power_expr("x", 2);
    SymbolicExpr pow_y = make_power_expr("y", 2);

    EXPECT_TRUE((ctx.current_relations().has_relation(pow_x, pow_y, RelationalNode::Op::GT))) << "x^2 > y^2 should be deduced when both are NonNegative and x > y";
}

TEST(AssumptionMonotonicity, PowerGuardMissingNonnegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());
    // No NonNegative assumption

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Add a power expression to provide exponent context
    SymbolicExpr x_squared = make_power_expr("x", 2);
    auto zero_expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    EXPECT_TRUE((ctx.current_relations().add_relation(x_squared, zero_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // x^2 > y^2 should NOT be deduced
    SymbolicExpr pow_x = make_power_expr("x", 2);
    SymbolicExpr pow_y = make_power_expr("y", 2);

    EXPECT_FALSE((ctx.current_relations().has_relation(pow_x, pow_y, RelationalNode::Op::GT))) << "Power rule should NOT apply when variables lack NonNegative";
}

TEST(AssumptionMonotonicity, RecursiveMonotonicityDepthLimit) {
    // With both Positive and Real, applying x > y should produce:
    // Level 0: ln(x) > ln(y), sqrt(x) > sqrt(y), exp(x) > exp(y)
    // The deduced relations have FunctionNode operands (not VariableNodes),
    // so recursion won't produce further deductions (guard: non-variable operands).
    // This test verifies the recursion doesn't crash and the depth limit works.

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    // This should not crash or infinite-loop
    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Verify all three rules were applied
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "ln(x) > ln(y) should be deduced";
    EXPECT_TRUE((ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT))) << "sqrt(x) > sqrt(y) should be deduced";
    EXPECT_TRUE((ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT))) << "exp(x) > exp(y) should be deduced";
}

TEST(AssumptionMonotonicity, AllRulesAppliedTogether) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("b", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr a_expr = make_var_expr("a");
    SymbolicExpr b_expr = make_var_expr("b");

    Relation rel{a_expr, b_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(a_expr, b_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Positive -> ln and sqrt rules apply
    // Real -> exp rule applies
    // Positive implies NonNegative -> power rule applies (if exponents exist)
    SymbolicExpr ln_a = make_func_expr(FunctionNode::FuncType::Ln, "a");
    SymbolicExpr ln_b = make_func_expr(FunctionNode::FuncType::Ln, "b");
    SymbolicExpr sqrt_a = make_func_expr(FunctionNode::FuncType::Sqrt, "a");
    SymbolicExpr sqrt_b = make_func_expr(FunctionNode::FuncType::Sqrt, "b");
    SymbolicExpr exp_a = make_func_expr(FunctionNode::FuncType::Exp, "a");
    SymbolicExpr exp_b = make_func_expr(FunctionNode::FuncType::Exp, "b");

    EXPECT_TRUE((ctx.current_relations().has_relation(ln_a, ln_b, RelationalNode::Op::GT))) << "ln(a) > ln(b) should be deduced";
    EXPECT_TRUE((ctx.current_relations().has_relation(sqrt_a, sqrt_b, RelationalNode::Op::GT))) << "sqrt(a) > sqrt(b) should be deduced";
    EXPECT_TRUE((ctx.current_relations().has_relation(exp_a, exp_b, RelationalNode::Op::GT))) << "exp(a) > exp(b) should be deduced";
}

/// Helper: create a closed interval [lo, hi] from numeric values.

/// Helper: create an open interval (lo, hi) from numeric values.
