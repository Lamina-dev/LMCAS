
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace lamina;

/// Helper: create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var_expr(const std::string& name) {
    return lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>(name));
}

/// Helper: create a FunctionNode expression (e.g., ln(x), sqrt(x), exp(x)).
static SymbolicExpr make_func_expr(FunctionNode::FuncType type, const std::string& var_name) {
    auto var_node = lamina::detail::make_node<VariableNode>(var_name);
    auto func_node = lamina::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{var_node});
    return lamina::detail::expression_from_node(func_node);
}

/// Helper: create a PowerNode expression (var^n).
static SymbolicExpr make_power_expr(const std::string& var_name, int n) {
    auto var_node = lamina::detail::make_node<VariableNode>(var_name);
    auto exp_node = lamina::detail::make_node<NumberNode>(BigInt(n));
    auto pow_node = lamina::detail::make_node<PowerNode>(var_node, exp_node);
    return lamina::detail::expression_from_node(pow_node);
}


void test_ln_monotonicity_both_positive() {
    TEST_CASE("x > y, both Positive → ln(x) > ln(y)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Add relation x > y
    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    // Apply monotonicity rules
    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that ln(x) > ln(y) was deduced
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
                "ln(x) > ln(y) should be deduced when both x,y are Positive and x > y");
}

void test_checked_monotonicity_rules_success() {
    TEST_CASE("Checked monotonicity rules: explicit success result");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    auto seed = ctx.current_relations().add_relation_checked(
        x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());
    EXPECT_TRUE(seed.has_value(), "checked seed relation succeeds");

    auto result = engine.apply_monotonicity_rules_checked(
        rel, ctx.current_relations(), ctx.current_properties());
    EXPECT_TRUE(result.has_value(), "checked monotonicity rules report success");

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
                "checked monotonicity rules deduce ln(x) > ln(y)");
}

void test_checked_monotonicity_rules_noop() {
    TEST_CASE("Checked monotonicity rules: unsupported relation is no-op success");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");
    Relation rel{x_expr, y_expr, RelationalNode::Op::LEQ};

    auto result = engine.apply_monotonicity_rules_checked(
        rel, ctx.current_relations(), ctx.current_properties());
    EXPECT_TRUE(result.has_value(), "checked monotonicity no-op reports success");
    EXPECT_TRUE(ctx.current_relations().get_relations().empty(),
                "checked monotonicity no-op leaves relation store unchanged");
}


void test_sqrt_monotonicity_both_positive() {
    TEST_CASE("x > y, both Positive → sqrt(x) > sqrt(y)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that sqrt(x) > sqrt(y) was deduced
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT),
                "sqrt(x) > sqrt(y) should be deduced when both x,y are Positive and x > y");
}


void test_exp_monotonicity_both_real() {
    TEST_CASE("x > y, both Real → exp(x) > exp(y)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that exp(x) > exp(y) was deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT),
                "exp(x) > exp(y) should be deduced when both x,y are Real and x > y");
}

void test_exp_monotonicity_positive_implies_real() {
    TEST_CASE("x > y, both Positive (implies Real) → exp(x) > exp(y)");

    AssumptionContext ctx;
    // Positive implies NonNegative and NonZero, but we also need Real domain
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that exp(x) > exp(y) was deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT),
                "exp(x) > exp(y) should be deduced when both are Positive+Real");
}


void test_ln_guard_missing_positive() {
    TEST_CASE("ln rule NOT applied when one variable lacks Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    // y has no sign assumption — guard should prevent ln rule

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // ln(x) > ln(y) should NOT be deduced
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_FALSE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
                 "ln rule should NOT apply when y lacks Positive assumption");
}

void test_sqrt_guard_missing_positive() {
    TEST_CASE("sqrt rule NOT applied when one variable lacks Positive");

    AssumptionContext ctx;
    // Neither x nor y has Positive assumption
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // sqrt(x) > sqrt(y) should NOT be deduced
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");

    EXPECT_FALSE(ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT),
                 "sqrt rule should NOT apply when variables lack Positive assumption");
}

void test_exp_guard_missing_real() {
    TEST_CASE("exp rule NOT applied when one variable lacks Real domain");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    // y has no domain assumption (defaults to Complex)

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // exp(x) > exp(y) should NOT be deduced
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_FALSE(ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT),
                 "exp rule should NOT apply when y lacks Real domain");
}

void test_no_rules_for_non_gt_relation() {
    TEST_CASE("No monotonicity rules applied for non-GT relations");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Use LT instead of GT
    Relation rel{x_expr, y_expr, RelationalNode::Op::LT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::LT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // No deduced relations should be added (only the original LT)
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_FALSE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
                 "No monotonicity rules should apply for non-GT relations");
}

void test_no_rules_for_non_variable_operands() {
    TEST_CASE("No monotonicity rules applied for non-variable operands");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);

    InferenceEngine engine(ctx);

    // Create a composite LHS: x + 1
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto one_node = lamina::detail::make_node<NumberNode>(BigInt(1));
    auto add_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, one_node});
    auto composite_expr = lamina::detail::expression_from_node(add_node);
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{composite_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(composite_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // No deduced relations should be added for composite operands
    // (only 1 relation in the store: the original one)
    EXPECT_TRUE(ctx.current_relations().get_relations().size() == 1,
                "No monotonicity rules should apply for non-variable operands");
}


void test_power_monotonicity_both_nonnegative() {
    TEST_CASE("x > y, both NonNegative → x^n > y^n for n in expressions");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);
    ctx.assume_sign("y", Sign::NonNegative);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // First, add a relation that contains a power expression with x^2
    // so that the exponent 2 is "appearing in expressions"
    SymbolicExpr x_squared = make_power_expr("x", 2);
    auto zero_expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(0)));
    ctx.current_relations().add_relation(x_squared, zero_expr, RelationalNode::Op::GT, ctx.current_properties());

    // Now add x > y and apply monotonicity
    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Check that x^2 > y^2 was deduced
    SymbolicExpr pow_x = make_power_expr("x", 2);
    SymbolicExpr pow_y = make_power_expr("y", 2);

    EXPECT_TRUE(ctx.current_relations().has_relation(pow_x, pow_y, RelationalNode::Op::GT),
                "x^2 > y^2 should be deduced when both are NonNegative and x > y");
}

void test_power_guard_missing_nonnegative() {
    TEST_CASE("Power rule NOT applied when variables lack NonNegative");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);
    // No NonNegative assumption

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Add a power expression to provide exponent context
    SymbolicExpr x_squared = make_power_expr("x", 2);
    auto zero_expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(0)));
    ctx.current_relations().add_relation(x_squared, zero_expr, RelationalNode::Op::GT, ctx.current_properties());

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // x^2 > y^2 should NOT be deduced
    SymbolicExpr pow_x = make_power_expr("x", 2);
    SymbolicExpr pow_y = make_power_expr("y", 2);

    EXPECT_FALSE(ctx.current_relations().has_relation(pow_x, pow_y, RelationalNode::Op::GT),
                 "Power rule should NOT apply when variables lack NonNegative");
}


void test_recursive_monotonicity_depth_limit() {
    TEST_CASE("Monotonicity rules applied recursively up to depth 8");

    // With both Positive and Real, applying x > y should produce:
    // Level 0: ln(x) > ln(y), sqrt(x) > sqrt(y), exp(x) > exp(y)
    // The deduced relations have FunctionNode operands (not VariableNodes),
    // so recursion won't produce further deductions (guard: non-variable operands).
    // This test verifies the recursion doesn't crash and the depth limit works.

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    // This should not crash or infinite-loop
    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Verify all three rules were applied
    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");
    SymbolicExpr sqrt_x = make_func_expr(FunctionNode::FuncType::Sqrt, "x");
    SymbolicExpr sqrt_y = make_func_expr(FunctionNode::FuncType::Sqrt, "y");
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
                "ln(x) > ln(y) should be deduced");
    EXPECT_TRUE(ctx.current_relations().has_relation(sqrt_x, sqrt_y, RelationalNode::Op::GT),
                "sqrt(x) > sqrt(y) should be deduced");
    EXPECT_TRUE(ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT),
                "exp(x) > exp(y) should be deduced");
}

void test_all_rules_applied_together() {
    TEST_CASE("All applicable rules applied when both Positive and Real");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);
    ctx.assume_domain("a", Domain::Real);
    ctx.assume_domain("b", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr a_expr = make_var_expr("a");
    SymbolicExpr b_expr = make_var_expr("b");

    Relation rel{a_expr, b_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(a_expr, b_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // Positive → ln and sqrt rules apply
    // Real → exp rule applies
    // Positive implies NonNegative → power rule applies (if exponents exist)
    SymbolicExpr ln_a = make_func_expr(FunctionNode::FuncType::Ln, "a");
    SymbolicExpr ln_b = make_func_expr(FunctionNode::FuncType::Ln, "b");
    SymbolicExpr sqrt_a = make_func_expr(FunctionNode::FuncType::Sqrt, "a");
    SymbolicExpr sqrt_b = make_func_expr(FunctionNode::FuncType::Sqrt, "b");
    SymbolicExpr exp_a = make_func_expr(FunctionNode::FuncType::Exp, "a");
    SymbolicExpr exp_b = make_func_expr(FunctionNode::FuncType::Exp, "b");

    EXPECT_TRUE(ctx.current_relations().has_relation(ln_a, ln_b, RelationalNode::Op::GT),
                "ln(a) > ln(b) should be deduced");
    EXPECT_TRUE(ctx.current_relations().has_relation(sqrt_a, sqrt_b, RelationalNode::Op::GT),
                "sqrt(a) > sqrt(b) should be deduced");
    EXPECT_TRUE(ctx.current_relations().has_relation(exp_a, exp_b, RelationalNode::Op::GT),
                "exp(a) > exp(b) should be deduced");
}


/// Helper: create a closed interval [lo, hi] from numeric values.
static Interval make_closed_interval(double lo, double hi) {
    auto lo_expr = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto hi_expr = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lo_expr);
    iv.upper = Endpoint::closed(hi_expr);
    return iv;
}

/// Helper: create an open interval (lo, hi) from numeric values.
static Interval make_open_interval(double lo, double hi) {
    auto lo_expr = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto hi_expr = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::open(lo_expr);
    iv.upper = Endpoint::open(hi_expr);
    return iv;
}

void test_exp_increasing_on_reals() {
    TEST_CASE("exp is auto-inferred as Increasing on all of R");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // exp(x) should be Increasing on any interval
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval entire = Interval::entire_line();
    Monotonicity mono = engine.infer_monotonicity(exp_x, "x", entire);

    EXPECT_TRUE(mono == Monotonicity::Increasing,
        "exp(x) is Increasing on entire real line");
}

void test_exp_increasing_on_finite_interval() {
    TEST_CASE("exp is Increasing on finite interval [0, 10]");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity mono = engine.infer_monotonicity(exp_x, "x", iv);

    EXPECT_TRUE(mono == Monotonicity::Increasing,
        "exp(x) is Increasing on [0, 10]");
}

void test_ln_increasing_on_positive_reals() {
    TEST_CASE("ln is auto-inferred as Increasing on positive reals");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");

    // ln is increasing on (0, +inf)
    Interval pos_reals = make_open_interval(0.0, 1000.0);
    Monotonicity mono = engine.infer_monotonicity(ln_x, "x", pos_reals);

    EXPECT_TRUE(mono == Monotonicity::Increasing,
        "ln(x) is Increasing on (0, 1000)");
}

void test_ln_increasing_on_closed_positive() {
    TEST_CASE("ln is Increasing on [1, 100]");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");

    Interval iv = make_closed_interval(1.0, 100.0);
    Monotonicity mono = engine.infer_monotonicity(ln_x, "x", iv);

    EXPECT_TRUE(mono == Monotonicity::Increasing,
        "ln(x) is Increasing on [1, 100]");
}

void test_negation_reverses_monotonicity() {
    TEST_CASE("Negation reverses monotonicity (-exp(x) is Decreasing)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Create -exp(x) = MultiplyNode([-1, exp(x)])
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto exp_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Exp,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto neg_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
    auto neg_exp = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, exp_node});
    auto neg_exp_x = lamina::detail::expression_from_node(neg_exp);
    Interval entire = Interval::entire_line();
    Monotonicity mono = engine.infer_monotonicity(neg_exp_x, "x", entire);

    EXPECT_TRUE(mono == Monotonicity::Decreasing,
        "-exp(x) is Decreasing (negation reverses Increasing)");
}

void test_negation_reverses_ln() {
    TEST_CASE("-ln(x) is Decreasing on positive reals");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // Create -ln(x)
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto ln_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto neg_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
    auto neg_ln = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, ln_node});
    auto neg_ln_x = lamina::detail::expression_from_node(neg_ln);
    Interval iv = make_closed_interval(1.0, 100.0);
    Monotonicity mono = engine.infer_monotonicity(neg_ln_x, "x", iv);

    EXPECT_TRUE(mono == Monotonicity::Decreasing,
        "-ln(x) is Decreasing on [1, 100]");
}

void test_declared_monotonicity_deduction() {
    TEST_CASE("Declared monotonically increasing f with x > y deduces f(x) > f(y)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    // Add x > y relation
    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    // Apply monotonicity rules — exp is auto-inferred increasing
    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    // exp(x) > exp(y) should be deduced (exp is increasing on R)
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT),
        "x > y with exp increasing => exp(x) > exp(y)");
}

void test_ln_deduction_from_inequality() {
    TEST_CASE("x > y with both Positive deduces ln(x) > ln(y)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties());

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE(ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT),
        "x > y with both Positive => ln(x) > ln(y)");
}

void test_unknown_for_non_monotone_function() {
    TEST_CASE("Non-monotone function (sin) returns Unknown monotonicity");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr sin_x = make_func_expr(FunctionNode::FuncType::Sin, "x");

    Interval iv = make_closed_interval(0.0, 2.0 * M_PI);
    Monotonicity mono = engine.infer_monotonicity(sin_x, "x", iv);

    EXPECT_TRUE(mono == Monotonicity::Unknown,
        "sin(x) on [0, 2*pi] has Unknown monotonicity (not monotone on full period)");
}

void test_wrong_variable_returns_unknown() {
    TEST_CASE("Querying monotonicity w.r.t. wrong variable returns Unknown");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // exp(x) is increasing w.r.t. x, but Unknown w.r.t. y
    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity mono = engine.infer_monotonicity(exp_x, "y", iv);

    EXPECT_TRUE(mono == Monotonicity::Unknown,
        "exp(x) w.r.t. y returns Unknown (wrong variable)");
}

int main() {
    test_exp_increasing_on_reals();
    test_exp_increasing_on_finite_interval();
    test_ln_increasing_on_positive_reals();
    test_ln_increasing_on_closed_positive();
    test_negation_reverses_monotonicity();
    test_negation_reverses_ln();
    test_declared_monotonicity_deduction();
    test_ln_deduction_from_inequality();
    test_unknown_for_non_monotone_function();
    test_wrong_variable_returns_unknown();

    // Properties 28-31: Monotonicity deduction rules (existing tests)
    test_ln_monotonicity_both_positive();
    test_checked_monotonicity_rules_success();
    test_checked_monotonicity_rules_noop();
    test_sqrt_monotonicity_both_positive();
    test_exp_monotonicity_both_real();
    test_exp_monotonicity_positive_implies_real();
    test_ln_guard_missing_positive();
    test_sqrt_guard_missing_positive();
    test_exp_guard_missing_real();
    test_no_rules_for_non_gt_relation();
    test_no_rules_for_non_variable_operands();
    test_power_monotonicity_both_nonnegative();
    test_power_guard_missing_nonnegative();
    test_recursive_monotonicity_depth_limit();
    test_all_rules_applied_together();

    return TEST_REPORT();
}
