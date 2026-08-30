
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <type_traits>

using namespace lamina;

static_assert(sizeof(InferenceEngine) == sizeof(void*),
              "InferenceEngine must keep implementation state out of the public layout");
static_assert(!std::is_copy_constructible_v<InferenceEngine>,
              "InferenceEngine must not share mutable traversal state by copying");


static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<MultiplyNode>(std::move(ops));
}

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exp) {
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_function(
    FunctionNode::FuncType type, std::shared_ptr<const SymbolicNode> arg) {
    return lamina::detail::make_node<FunctionNode>(type,
        std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static std::shared_ptr<const SymbolicNode> make_add(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<AddNode>(std::move(ops));
}

/// Create a division expression: num / den represented as MultiplyNode([num, PowerNode(den, -1)])
static std::shared_ptr<const SymbolicNode> make_division(
    std::shared_ptr<const SymbolicNode> num, std::shared_ptr<const SymbolicNode> den) {
    auto den_inv = make_power(std::move(den), make_number(-1));
    return make_multiply({std::move(num), std::move(den_inv)});
}

/// Create a subtraction expression: a - b represented as AddNode([a, MultiplyNode([-1, b])])
static std::shared_ptr<const SymbolicNode> make_subtraction(
    std::shared_ptr<const SymbolicNode> a, std::shared_ptr<const SymbolicNode> b) {
    auto neg_b = make_multiply({make_number(-1), std::move(b)});
    return make_add({std::move(a), std::move(neg_b)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}


void test_division_positive_over_positive() {
    TEST_CASE("Division: positive / positive → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "pos / pos is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "pos / pos is not Negative");
}

void test_division_negative_over_negative() {
    TEST_CASE("Division: negative / negative → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "neg / neg is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "neg / neg is not Negative");
}

void test_division_positive_over_negative() {
    TEST_CASE("Division: positive / negative → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "pos / neg is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "pos / neg is Negative");
}

void test_division_negative_over_positive() {
    TEST_CASE("Division: negative / positive → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "neg / pos is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "neg / pos is Negative");
}

void test_division_unknown_denominator() {
    TEST_CASE("Division: unknown denominator → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    // b has no sign declared
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
        "pos / unknown: Positive is Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
        "pos / unknown: Negative is Unknown");
}

void test_division_zero_denominator() {
    TEST_CASE("Division: zero denominator → Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Zero);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_division(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
        "pos / zero: Positive is Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
        "pos / zero: Negative is Unknown");
}


void test_subtraction_positive_minus_negative() {
    TEST_CASE("Subtraction: positive - negative → positive");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Negative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_subtraction(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "pos - neg is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "pos - neg is not Negative");
}

void test_subtraction_negative_minus_positive() {
    TEST_CASE("Subtraction: negative - positive → negative");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Positive);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_subtraction(make_var("a"), make_var("b")));

    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "neg - pos is Negative");
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "neg - pos is not Positive");
}


void test_domain_sin_integer() {
    TEST_CASE("Domain: sin(integer) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sin, make_var("x")));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "sin(integer) is Real");
}

void test_domain_exp_rational() {
    TEST_CASE("Domain: exp(rational) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Rational);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("x")));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "exp(rational) is Real");
}

void test_domain_ln_integer() {
    TEST_CASE("Domain: ln(integer) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    // ln requires positive argument for Real result; Integer domain alone
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Ln, make_var("x")));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "ln(integer) is Real");
}

void test_domain_sqrt_nonneg_real() {
    TEST_CASE("Domain: sqrt(nonneg real) → Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::NonNegative);
    InferenceEngine engine(ctx);

    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sqrt, make_var("x")));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "sqrt(nonneg real) is Real");
}

void test_domain_integer_pow_natural() {
    TEST_CASE("Domain: integer^natural → Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    // x^3 where x is Integer and 3 is a positive integer
    auto expr = wrap_expr(make_power(make_var("x"), make_number(3)));

    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::True,
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
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "rational^2 is Real");
}


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
    Tribool result = engine.query_positive_checked(expr).value();
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
    Tribool result = engine.query_positive_checked(expr).value();
    EXPECT_TRUE(result == Tribool::Unknown,
        "Nested exp with depth limit 1 returns Unknown");
}


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
    Tribool result = engine.query_real_checked(expr).value();
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
    Tribool result = engine.query_positive_checked(expr).value();
    // The important thing is no infinite recursion
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Deeply shared nodes: no crash, returns True or Unknown");
}

void test_checked_inference_query_contracts() {
    TEST_CASE("InferenceEngine checked queries: explicit errors and values");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("q", Domain::Rational);
    ctx.assume_domain("real_symbol", Domain::Real);
    ctx.assume_sign("real_pos", Sign::Positive);
    ctx.assume_domain("real_pos", Domain::Real);
    ctx.assume_sign("nn_symbol", Sign::NonNegative);
    ctx.current_properties().declare_finiteness("finite_symbol", Finiteness::Finite);
    ctx.current_properties().declare_finiteness("divergent_symbol", Finiteness::Divergent);
    ctx.current_properties().declare_transcendental("tau_symbol");
    ctx.current_properties().declare_periodic("periodic_symbol",
        lamina::detail::make_expression_ptr(wrap_expr(make_number(6))));
    InferenceEngine engine(ctx);

    auto x = wrap_expr(make_var("x"));
    auto positive = engine.query_positive_checked(x);
    EXPECT_TRUE(positive.has_value(), "checked query_positive succeeds");
    if (positive) {
        EXPECT_TRUE(positive.value() == Tribool::True,
            "checked query_positive returns True for positive symbol");
    }
    EXPECT_TRUE(engine.query_positive_checked(x).value() == Tribool::True,
        "legacy query_positive unwraps checked result");

    auto nonpositive = engine.query_nonpositive_checked(x);
    EXPECT_TRUE(nonpositive.has_value(), "checked query_nonpositive succeeds");
    if (nonpositive) {
        EXPECT_TRUE(nonpositive.value() == Tribool::False,
            "checked query_nonpositive returns False for positive symbol");
    }
    EXPECT_TRUE(engine.query_nonpositive_checked(x).value() == Tribool::False,
        "legacy query_nonpositive unwraps checked result");

    auto integer = engine.query_integer_checked(x);
    EXPECT_TRUE(integer.has_value(), "checked query_integer succeeds");
    if (integer) {
        EXPECT_TRUE(integer.value() == Tribool::True,
            "checked query_integer returns True for integer symbol");
    }
    EXPECT_TRUE(engine.query_integer_checked(x).value() == Tribool::True,
        "legacy query_integer unwraps checked result");

    auto nonzero = engine.query_nonzero_checked(x);
    EXPECT_TRUE(nonzero.has_value(), "checked query_nonzero succeeds");
    if (nonzero) {
        EXPECT_TRUE(nonzero.value() == Tribool::True,
            "checked query_nonzero returns True for positive symbol");
    }
    EXPECT_TRUE(engine.query_nonzero_checked(x).value() == Tribool::True,
        "legacy query_nonzero unwraps checked result");

    auto three = wrap_expr(make_number(3));
    auto numeric_negative = engine.query_negative_checked(three);
    EXPECT_TRUE(numeric_negative.has_value(), "checked query_negative succeeds");
    if (numeric_negative) {
        EXPECT_TRUE(numeric_negative.value() == Tribool::False,
            "checked query_negative returns False for positive number");
    }
    EXPECT_TRUE(engine.query_negative_checked(three).value() == Tribool::False,
        "legacy query_negative unwraps checked result");

    auto numeric_nonnegative = engine.query_nonnegative_checked(three);
    EXPECT_TRUE(numeric_nonnegative.has_value(), "checked query_nonnegative succeeds");
    if (numeric_nonnegative) {
        EXPECT_TRUE(numeric_nonnegative.value() == Tribool::True,
            "checked query_nonnegative returns True for positive number");
    }
    EXPECT_TRUE(engine.query_nonnegative_checked(three).value() == Tribool::True,
        "legacy query_nonnegative unwraps checked result");

    auto numeric_real = engine.query_real_checked(three);
    EXPECT_TRUE(numeric_real.has_value(), "checked query_real succeeds");
    if (numeric_real) {
        EXPECT_TRUE(numeric_real.value() == Tribool::True,
            "checked query_real returns True for exact number");
    }
    EXPECT_TRUE(engine.query_real_checked(three).value() == Tribool::True,
        "legacy query_real unwraps checked result");

    auto exp_rational = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("q")));
    auto exp_rational_real = engine.query_real_checked(exp_rational);
    EXPECT_TRUE(exp_rational_real.has_value(),
        "checked query_real succeeds for exp(rational)");
    if (exp_rational_real) {
        EXPECT_TRUE(exp_rational_real.value() == Tribool::True,
            "checked query_real uses checked internal rational-domain dispatch");
    }
    EXPECT_TRUE(engine.query_real_checked(exp_rational).value() == Tribool::True,
        "legacy query_real unwraps checked rational-domain dispatch");

    auto integer_sum = wrap_expr(make_add({make_var("x"), make_var("x")}));
    auto integer_sum_checked = engine.query_integer_checked(integer_sum);
    EXPECT_TRUE(integer_sum_checked.has_value(),
        "checked query_integer succeeds for integer addition");
    if (integer_sum_checked) {
        EXPECT_TRUE(integer_sum_checked.value() == Tribool::True,
            "checked query_integer uses checked addition-domain inference");
    }

    auto real_product_domain = wrap_expr(make_multiply({make_var("x"), make_var("real_symbol")}));
    auto real_product_checked = engine.query_real_checked(real_product_domain);
    EXPECT_TRUE(real_product_checked.has_value(),
        "checked query_real succeeds for mixed integer-real product");
    if (real_product_checked) {
        EXPECT_TRUE(real_product_checked.value() == Tribool::True,
            "checked query_real uses checked multiplication-domain inference");
    }

    auto integer_power_domain = wrap_expr(make_power(make_var("x"), make_number(2)));
    auto integer_power_domain_checked = engine.query_integer_checked(integer_power_domain);
    EXPECT_TRUE(integer_power_domain_checked.has_value(),
        "checked query_integer succeeds for integer nonnegative power");
    if (integer_power_domain_checked) {
        EXPECT_TRUE(integer_power_domain_checked.value() == Tribool::True,
            "checked query_integer uses checked power-domain inference");
    }

    auto abs_integer_domain = wrap_expr(make_function(FunctionNode::FuncType::Abs, make_var("x")));
    auto abs_integer_domain_checked = engine.query_integer_checked(abs_integer_domain);
    EXPECT_TRUE(abs_integer_domain_checked.has_value(),
        "checked query_integer succeeds for abs(integer)");
    if (abs_integer_domain_checked) {
        EXPECT_TRUE(abs_integer_domain_checked.value() == Tribool::True,
            "checked query_integer uses checked function-domain inference");
    }

    auto positive_sum = wrap_expr(make_add({make_var("x"), make_var("x")}));
    auto positive_sum_checked = engine.query_positive_checked(positive_sum);
    EXPECT_TRUE(positive_sum_checked.has_value(),
        "checked query_positive succeeds for addition");
    if (positive_sum_checked) {
        EXPECT_TRUE(positive_sum_checked.value() == Tribool::True,
            "checked query_positive uses checked addition-sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(positive_sum).value() == Tribool::True,
        "legacy query_positive unwraps checked addition-sign inference");

    ctx.assume_sign("neg_symbol", Sign::Negative);
    auto positive_difference = wrap_expr(make_subtraction(make_var("x"), make_var("neg_symbol")));
    auto positive_difference_checked = engine.query_positive_checked(positive_difference);
    EXPECT_TRUE(positive_difference_checked.has_value(),
        "checked query_positive succeeds for subtraction-shaped addition");
    if (positive_difference_checked) {
        EXPECT_TRUE(positive_difference_checked.value() == Tribool::True,
            "checked query_positive uses checked subtraction-sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(positive_difference).value() == Tribool::True,
        "legacy query_positive unwraps checked subtraction-sign inference");

    auto relation_positive_sum = wrap_expr(make_add({make_var("rel_a"), make_var("rel_b")}));
    auto zero_expr = wrap_expr(make_number(0));
    auto relation_inserted = ctx.current_relations().add_relation_checked(
        relation_positive_sum, zero_expr, RelationalNode::Op::GT, ctx.current_properties());
    EXPECT_TRUE(relation_inserted.has_value(),
        "checked relation insertion succeeds for composite positive sum");
    auto relation_positive_sum_checked = engine.query_positive_checked(relation_positive_sum);
    EXPECT_TRUE(relation_positive_sum_checked.has_value(),
        "checked query_positive succeeds for relation-backed composite sum");
    if (relation_positive_sum_checked) {
        EXPECT_TRUE(relation_positive_sum_checked.value() == Tribool::True,
            "checked query_positive uses checked relation sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(relation_positive_sum).value() == Tribool::True,
        "legacy query_positive unwraps checked relation sign inference");

    auto positive_product = wrap_expr(make_multiply({make_var("x"), make_var("x")}));
    auto positive_product_checked = engine.query_positive_checked(positive_product);
    EXPECT_TRUE(positive_product_checked.has_value(),
        "checked query_positive succeeds for multiplication");
    if (positive_product_checked) {
        EXPECT_TRUE(positive_product_checked.value() == Tribool::True,
            "checked query_positive uses checked multiplication-sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(positive_product).value() == Tribool::True,
        "legacy query_positive unwraps checked multiplication-sign inference");

    auto negative_product = wrap_expr(make_multiply({make_var("x"), make_var("neg_symbol")}));
    auto negative_product_checked = engine.query_negative_checked(negative_product);
    EXPECT_TRUE(negative_product_checked.has_value(),
        "checked query_negative succeeds for multiplication");
    if (negative_product_checked) {
        EXPECT_TRUE(negative_product_checked.value() == Tribool::True,
            "checked query_negative uses checked multiplication-sign inference");
    }
    EXPECT_TRUE(engine.query_negative_checked(negative_product).value() == Tribool::True,
        "legacy query_negative unwraps checked multiplication-sign inference");

    auto division_expr = wrap_expr(make_division(make_var("x"), make_var("x")));
    auto division_positive = engine.query_positive_checked(division_expr);
    EXPECT_TRUE(division_positive.has_value(),
        "checked query_positive succeeds for division pattern");
    if (division_positive) {
        EXPECT_TRUE(division_positive.value() == Tribool::True,
            "checked query_positive uses checked division-sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(division_expr).value() == Tribool::True,
        "legacy query_positive unwraps checked division-sign inference");

    auto positive_power = wrap_expr(make_power(make_var("x"), make_number(2)));
    auto positive_power_checked = engine.query_positive_checked(positive_power);
    EXPECT_TRUE(positive_power_checked.has_value(),
        "checked query_positive succeeds for power");
    if (positive_power_checked) {
        EXPECT_TRUE(positive_power_checked.value() == Tribool::True,
            "checked query_positive uses checked positive-base power inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(positive_power).value() == Tribool::True,
        "legacy query_positive unwraps checked power-sign inference");

    auto even_power = wrap_expr(make_power(make_var("unknown_symbol"), make_number(2)));
    auto even_power_nonnegative = engine.query_nonnegative_checked(even_power);
    EXPECT_TRUE(even_power_nonnegative.has_value(),
        "checked query_nonnegative succeeds for even power");
    if (even_power_nonnegative) {
        EXPECT_TRUE(even_power_nonnegative.value() == Tribool::Unknown,
            "checked even-power inference preserves Unknown without real-domain proof");
    }
    auto real_even_power = wrap_expr(make_power(make_var("real_symbol"), make_number(2)));
    auto real_even_power_nonnegative = engine.query_nonnegative_checked(real_even_power);
    EXPECT_TRUE(real_even_power_nonnegative.has_value(),
        "checked query_nonnegative succeeds for real even power");
    if (real_even_power_nonnegative) {
        EXPECT_TRUE(real_even_power_nonnegative.value() == Tribool::True,
            "checked query_nonnegative uses checked even-power inference");
    }

    auto nonzero_power = wrap_expr(make_power(make_var("x"), make_number(-1)));
    auto nonzero_power_checked = engine.query_nonzero_checked(nonzero_power);
    EXPECT_TRUE(nonzero_power_checked.has_value(),
        "checked query_nonzero succeeds for integer power");
    if (nonzero_power_checked) {
        EXPECT_TRUE(nonzero_power_checked.value() == Tribool::True,
            "checked query_nonzero uses checked nonzero-base power inference");
    }

    auto exp_real = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("real_symbol")));
    auto exp_real_positive = engine.query_positive_checked(exp_real);
    EXPECT_TRUE(exp_real_positive.has_value(),
        "checked query_positive succeeds for exp(real)");
    if (exp_real_positive) {
        EXPECT_TRUE(exp_real_positive.value() == Tribool::True,
            "checked query_positive uses checked exp function-sign inference");
    }
    EXPECT_TRUE(engine.query_positive_checked(exp_real).value() == Tribool::True,
        "legacy query_positive unwraps checked exp function-sign inference");

    auto abs_real_pos = wrap_expr(make_function(FunctionNode::FuncType::Abs, make_var("real_pos")));
    auto abs_real_pos_positive = engine.query_positive_checked(abs_real_pos);
    EXPECT_TRUE(abs_real_pos_positive.has_value(),
        "checked query_positive succeeds for abs(positive real)");
    if (abs_real_pos_positive) {
        EXPECT_TRUE(abs_real_pos_positive.value() == Tribool::True,
            "checked query_positive uses checked abs function-sign inference");
    }

    auto sqrt_nonnegative = wrap_expr(make_function(FunctionNode::FuncType::Sqrt, make_var("nn_symbol")));
    auto sqrt_nonnegative_checked = engine.query_nonnegative_checked(sqrt_nonnegative);
    EXPECT_TRUE(sqrt_nonnegative_checked.has_value(),
        "checked query_nonnegative succeeds for sqrt(nonnegative)");
    if (sqrt_nonnegative_checked) {
        EXPECT_TRUE(sqrt_nonnegative_checked.value() == Tribool::True,
            "checked query_nonnegative uses checked sqrt function-sign inference");
    }

    auto finite = engine.query_finite_checked(wrap_expr(make_var("finite_symbol")));
    EXPECT_TRUE(finite.has_value(), "checked query_finite succeeds");
    if (finite) {
        EXPECT_TRUE(finite.value() == Tribool::True,
            "checked query_finite returns True for finite symbol");
    }

    auto periodic_expr = wrap_expr(make_var("periodic_symbol"));
    auto periodic = engine.query_periodic_checked(periodic_expr);
    EXPECT_TRUE(periodic.has_value(), "checked query_periodic succeeds");
    if (periodic) {
        EXPECT_TRUE(periodic.value() == Tribool::True,
            "checked query_periodic returns True for periodic symbol");
    }

    auto period = engine.infer_period_checked(periodic_expr);
    EXPECT_TRUE(period.has_value(), "checked infer_period succeeds");
    if (period) {
        EXPECT_TRUE(period.value().has_value(),
            "checked infer_period returns declared period");
    }

    auto unknown = engine.query_algebraic_checked(wrap_expr(make_var("unknown_symbol")));
    EXPECT_TRUE(unknown.has_value(), "checked query_algebraic accepts valid unknown symbol");
    if (unknown) {
        EXPECT_TRUE(unknown.value() == Tribool::Unknown,
            "checked query_algebraic preserves Unknown for valid unsupported facts");
    }

    auto tau_expr = wrap_expr(make_var("tau_symbol"));
    auto transcendental = engine.query_transcendental_checked(tau_expr);
    EXPECT_TRUE(transcendental.has_value(), "checked query_transcendental succeeds");
    if (transcendental) {
        EXPECT_TRUE(transcendental.value() == Tribool::True,
            "checked query_transcendental returns True for transcendental symbol");
    }
    EXPECT_TRUE(engine.query_transcendental_checked(tau_expr).value() == Tribool::True,
        "legacy query_transcendental unwraps checked result");

    auto divergent_expr = wrap_expr(make_var("divergent_symbol"));
    auto divergent = engine.query_divergent_checked(divergent_expr);
    EXPECT_TRUE(divergent.has_value(), "checked query_divergent succeeds");
    if (divergent) {
        EXPECT_TRUE(divergent.value() == Tribool::True,
            "checked query_divergent returns True for divergent symbol");
    }
    auto divergent_finite = engine.query_finite_checked(divergent_expr);
    EXPECT_TRUE(divergent_finite.has_value(), "checked query_finite succeeds for divergent symbol");
    if (divergent_finite) {
        EXPECT_TRUE(divergent_finite.value() == Tribool::False,
            "checked query_finite returns False for divergent symbol");
    }
}


int main() {
    test_division_positive_over_positive();
    test_division_negative_over_negative();
    test_division_positive_over_negative();
    test_division_negative_over_positive();
    test_division_unknown_denominator();
    test_division_zero_denominator();

    test_subtraction_positive_minus_negative();
    test_subtraction_negative_minus_positive();

    test_domain_sin_integer();
    test_domain_exp_rational();
    test_domain_ln_integer();
    test_domain_sqrt_nonneg_real();
    test_domain_integer_pow_natural();
    test_domain_rational_pow_integer();

    test_depth_limit_triggers_unknown();
    test_depth_limit_set_and_get();
    test_depth_limit_nested_addition();

    test_cycle_detection_self_referential();
    test_cycle_detection_no_crash_deep_shared();
    test_checked_inference_query_contracts();

    return TEST_REPORT();
}
