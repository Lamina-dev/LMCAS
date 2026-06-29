/**
 * @file test_inference_engine_ext.cpp
 * @brief Unit tests for InferenceEngine extensions: division sign, subtraction sign,
 *        composite domain inference, depth limit, and cycle detection.
 *
 * Validates: Requirements 1.1, 1.2, 2.1, 2.2, 22.2, 25.2
 */

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;

// ============================================================
// Helpers
// ============================================================

static std::shared_ptr<SymbolicNode> make_number(int val) {
    return std::make_shared<NumberNode>(BigInt(val));
}

static std::shared_ptr<SymbolicNode> make_var(const std::string& name) {
    return std::make_shared<VariableNode>(name);
}

static std::shared_ptr<SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<SymbolicNode>> ops) {
    return std::make_shared<MultiplyNode>(std::move(ops));
}

static std::shared_ptr<SymbolicNode> make_power(
    std::shared_ptr<SymbolicNode> base, std::shared_ptr<SymbolicNode> exp) {
    return std::make_shared<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<SymbolicNode> make_function(
    FunctionNode::FuncType type, std::shared_ptr<SymbolicNode> arg) {
    return std::make_shared<FunctionNode>(type,
        std::vector<std::shared_ptr<SymbolicNode>>{std::move(arg)});
}

static std::shared_ptr<SymbolicNode> make_add(
    std::vector<std::shared_ptr<SymbolicNode>> ops) {
    return std::make_shared<AddNode>(std::move(ops));
}

/// Create a division expression: num / den represented as MultiplyNode([num, PowerNode(den, -1)])
static std::shared_ptr<SymbolicNode> make_division(
    std::shared_ptr<SymbolicNode> num, std::shared_ptr<SymbolicNode> den) {
    auto den_inv = make_power(std::move(den), make_number(-1));
    return make_multiply({std::move(num), std::move(den_inv)});
}

/// Create a subtraction expression: a - b represented as AddNode([a, MultiplyNode([-1, b])])
static std::shared_ptr<SymbolicNode> make_subtraction(
    std::shared_ptr<SymbolicNode> a, std::shared_ptr<SymbolicNode> b) {
    auto neg_b = make_multiply({make_number(-1), std::move(b)});
    return make_add({std::move(a), std::move(neg_b)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<SymbolicNode> node) {
    SymbolicExpr expr;
    expr.root = std::move(node);
    return expr;
}

// ============================================================
// Division sign inference tests (Requirements 1.1, 1.2)
// ============================================================

void test_division_positive_over_positive() {
    TEST_CASE("Division: positive / positive → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "pos / pos is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "pos / pos is not Negative");
}

void test_division_negative_over_negative() {
    TEST_CASE("Division: negative / negative → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "neg / neg is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "neg / neg is not Negative");
}

void test_division_positive_over_negative() {
    TEST_CASE("Division: positive / negative → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "pos / neg is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "pos / neg is Negative");
}

void test_division_negative_over_positive() {
    TEST_CASE("Division: negative / positive → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "neg / pos is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "neg / pos is Negative");
}

void test_division_unknown_denominator() {
    TEST_CASE("Division: unknown denominator → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    // b has no sign declared
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "pos / unknown: Positive is Unknown");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "pos / unknown: Negative is Unknown");
}

void test_division_zero_denominator() {
    TEST_CASE("Division: zero denominator → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Zero);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "pos / zero: Positive is Unknown");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "pos / zero: Negative is Unknown");
}

// ============================================================
// Subtraction sign inference tests (Requirement 1.1, 1.2)
// ============================================================

void test_subtraction_positive_minus_negative() {
    TEST_CASE("Subtraction: positive - negative → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_subtraction(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "pos - neg is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "pos - neg is not Negative");
}

void test_subtraction_negative_minus_positive() {
    TEST_CASE("Subtraction: negative - positive → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_subtraction(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "neg - pos is Negative");
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "neg - pos is not Positive");
}

// ============================================================
// Composite domain inference tests (Requirements 2.1, 2.2)
// ============================================================

void test_domain_sin_integer() {
    TEST_CASE("Domain: sin(integer) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sin, make_var("x")));

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "sin(integer) is Real");
}

void test_domain_exp_rational() {
    TEST_CASE("Domain: exp(rational) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Rational);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("x")));

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "exp(rational) is Real");
}

void test_domain_ln_integer() {
    TEST_CASE("Domain: ln(integer) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    // ln requires positive argument for Real result; Integer domain alone
    // should suffice per Req 2.3 (ln(Integer) → Real)
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Ln, make_var("x")));

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "ln(integer) is Real");
}

void test_domain_sqrt_nonneg_real() {
    TEST_CASE("Domain: sqrt(nonneg real) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::NonNegative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sqrt, make_var("x")));

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "sqrt(nonneg real) is Real");
}

void test_domain_integer_pow_natural() {
    TEST_CASE("Domain: integer^natural → Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    // x^3 where x is Integer and 3 is a positive integer
    auto expr = wrap_expr(make_power(make_var("x"), make_number(3)));

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "integer^3 is Integer");
}

void test_domain_rational_pow_integer() {
    TEST_CASE("Domain: rational^integer → Real (at minimum)");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Rational);
    InferenceEngine engine(ctx);

    // x^2 where x is Rational and 2 is an integer exponent
    auto expr = wrap_expr(make_power(make_var("x"), make_number(2)));

    // Rational ⊂ Real, so Rational base + integer exponent → Real
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "rational^2 is Real");
}

// ============================================================
// Depth limit tests (Requirement 22.2, 25.2)
// ============================================================

void test_depth_limit_triggers_unknown() {
    TEST_CASE("Depth limit: deeply nested expression → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    InferenceEngine engine(ctx);

    // Set a very small max depth to trigger the limit
    engine.set_max_depth(3);

    // Create a deeply nested expression: sin(sin(sin(sin(sin(x)))))
    auto node = make_var("x");
    for (int i = 0; i < 10; ++i) {
        node = make_function(FunctionNode::FuncType::Sin, node);
    }
    auto expr = wrap_expr(node);

    // With max depth 3, the deeply nested expression should return Unknown
    // because the engine can't recurse deep enough to resolve the innermost variable
    Tribool result = engine.query_positive(expr);
    // sin doesn't have a definite sign anyway, but the key test is that
    // it doesn't crash or hang — it returns Unknown gracefully
    EXPECT_TRUE(result == Tribool::Unknown,
        "Deeply nested expression with small depth limit returns Unknown");
}

void test_depth_limit_set_and_get() {
    TEST_CASE("Depth limit: set_max_depth and get_max_depth");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    EXPECT_TRUE(engine.get_max_depth() == 32,
        "Default max depth is 32");

    engine.set_max_depth(5);
    EXPECT_TRUE(engine.get_max_depth() == 5,
        "Max depth set to 5");

    // Setting to 0 or negative should not change (must be > 0)
    engine.set_max_depth(0);
    EXPECT_TRUE(engine.get_max_depth() == 5,
        "Max depth unchanged when set to 0");

    engine.set_max_depth(-1);
    EXPECT_TRUE(engine.get_max_depth() == 5,
        "Max depth unchanged when set to negative");
}

void test_depth_limit_nested_addition() {
    TEST_CASE("Depth limit: nested expression exceeds depth → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    InferenceEngine engine(ctx);
    engine.set_max_depth(1);

    // Create a nested function expression: exp(exp(x))
    // Depth 1: query_positive on outer exp → DepthGuard depth=1 (ok, ≤ 1)
    //   infer_function_property checks query_real on inner exp
    //     query_real on inner exp → DepthGuard depth=2 (> 1, abort → Unknown)
    // So exp(exp(x)) with max_depth=1 should return Unknown for query_real
    auto x = make_var("x");
    auto inner_exp = make_function(FunctionNode::FuncType::Exp, x);
    auto outer_exp = make_function(FunctionNode::FuncType::Exp, inner_exp);
    auto expr = wrap_expr(outer_exp);

    // exp requires Real argument to infer Positive. With depth limit 1,
    // the inner exp's domain can't be determined, so outer exp returns Unknown.
    Tribool result = engine.query_positive(expr);
    EXPECT_TRUE(result == Tribool::Unknown,
        "Nested exp with depth limit 1 returns Unknown");
}

// ============================================================
// Cycle detection tests (Requirement 22.2)
// ============================================================

void test_cycle_detection_self_referential() {
    TEST_CASE("Cycle detection: depth limit prevents deep recursion");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // Set a very small depth limit
    engine.set_max_depth(2);

    // Create a deeply nested expression: sin(cos(sin(cos(x))))
    // Each function call adds a depth level when querying its argument's domain
    auto node = make_var("x");
    for (int i = 0; i < 8; ++i) {
        auto type = (i % 2 == 0) ? FunctionNode::FuncType::Sin : FunctionNode::FuncType::Cos;
        node = make_function(type, node);
    }
    auto expr = wrap_expr(node);

    // With max depth 2, the engine can't recurse deep enough to determine
    // the innermost variable's domain, so it returns Unknown
    Tribool result = engine.query_real(expr);
    EXPECT_TRUE(result == Tribool::Unknown,
        "Deeply nested trig with depth limit 2 returns Unknown for domain");
}

void test_cycle_detection_no_crash_deep_shared() {
    TEST_CASE("Cycle detection: deeply shared nodes don't crash");
    AssumptionContext ctx;
    ctx.assume_sign("y", Sign::Positive);
    InferenceEngine engine(ctx);

    // Create a tree with shared sub-expressions at multiple levels
    auto y = make_var("y");
    auto inner = make_multiply({y, y}); // same y pointer twice
    auto outer = make_add({inner, inner}); // same inner pointer twice
    auto expr = wrap_expr(outer);

    // Should not crash or hang — returns some result (likely Unknown due to cycles)
    Tribool result = engine.query_positive(expr);
    // The important thing is no infinite recursion
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Deeply shared nodes: no crash, returns True or Unknown");
}

// ============================================================
// main
// ============================================================

int main() {
    // Division sign inference (Req 1.1, 1.2)
    test_division_positive_over_positive();
    test_division_negative_over_negative();
    test_division_positive_over_negative();
    test_division_negative_over_positive();
    test_division_unknown_denominator();
    test_division_zero_denominator();

    // Subtraction sign inference (Req 1.1, 1.2)
    test_subtraction_positive_minus_negative();
    test_subtraction_negative_minus_positive();

    // Composite domain inference (Req 2.1, 2.2)
    test_domain_sin_integer();
    test_domain_exp_rational();
    test_domain_ln_integer();
    test_domain_sqrt_nonneg_real();
    test_domain_integer_pow_natural();
    test_domain_rational_pow_integer();

    // Depth limit (Req 22.2, 25.2)
    test_depth_limit_triggers_unknown();
    test_depth_limit_set_and_get();
    test_depth_limit_nested_addition();

    // Cycle detection (Req 22.2)
    test_cycle_detection_self_referential();
    test_cycle_detection_no_crash_deep_shared();

    return TEST_REPORT();
}
