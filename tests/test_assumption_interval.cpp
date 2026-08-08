
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <cmath>

using namespace lamina;


/// Create a VariableNode wrapped in a shared_ptr<SymbolicNode>
static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

/// Create a NumberNode from a double value
static std::shared_ptr<const SymbolicNode> make_num(double val) {
    return lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val));
}

/// Create a NumberNode from an int value (BigInt)
static std::shared_ptr<const SymbolicNode> make_int_num(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

/// Wrap a SymbolicNode into a SymbolicExpr for querying
static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

/// Declare bounded interval [lo, hi] for a variable in the context
static void declare_bounds(AssumptionContext& ctx, const std::string& var, double lo, double hi) {
    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    ctx.current_properties().declare_bounded(var, Boundedness::Bounded, bounds);
}

static std::shared_ptr<SymbolicExpr> expr_from_node(std::shared_ptr<const SymbolicNode> node) {
    return lamina::detail::make_expression_ptr(std::move(node));
}

static void declare_expr_bounds(
    AssumptionContext& ctx,
    const std::string& var,
    std::shared_ptr<SymbolicExpr> lo,
    std::shared_ptr<SymbolicExpr> hi)
{
    Interval bounds;
    bounds.lower = Endpoint::closed(std::move(lo));
    bounds.upper = Endpoint::closed(std::move(hi));
    ctx.current_properties().declare_bounded(var, Boundedness::Bounded, bounds);
}

/// Extract numeric lower bound from an interval
static double get_lower(const Interval& iv) {
    if (iv.lower.is_neg_infinity) return -std::numeric_limits<double>::infinity();
    if (iv.lower.value) return iv.lower.value->to_numeric();
    return 0.0;
}

/// Extract numeric upper bound from an interval
static double get_upper(const Interval& iv) {
    if (iv.upper.is_pos_infinity) return std::numeric_limits<double>::infinity();
    if (iv.upper.value) return iv.upper.value->to_numeric();
    return 0.0;
}


void test_addition_propagation() {
    TEST_CASE("Addition: [1,3] + [2,5] = [3,8]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 3.0);
    declare_bounds(ctx, "y", 2.0, 5.0);

    // x + y
    auto add_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), make_var("y")});
    auto expr = wrap_expr(add_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Addition produces a bounded interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 3.0, 1e-10, "Lower bound is 1+2=3");
        EXPECT_NEAR(get_upper(*result), 8.0, 1e-10, "Upper bound is 3+5=8");
    }
}

void test_addition_negative_bounds() {
    TEST_CASE("Addition: [-2,1] + [-3,4] = [-5,5]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -2.0, 1.0);
    declare_bounds(ctx, "y", -3.0, 4.0);

    auto add_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), make_var("y")});
    auto expr = wrap_expr(add_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Addition with negatives produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), -5.0, 1e-10, "Lower bound is -2+(-3)=-5");
        EXPECT_NEAR(get_upper(*result), 5.0, 1e-10, "Upper bound is 1+4=5");
    }
}

void test_numeric_expression_endpoint_propagation() {
    TEST_CASE("Addition: expression endpoints [1+1, 3+2] + 1 = [3,6]");
    AssumptionContext ctx;
    auto lo = expr_from_node(lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_int_num(1), make_int_num(1)}));
    auto hi = expr_from_node(lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_int_num(3), make_int_num(2)}));
    declare_expr_bounds(ctx, "x", lo, hi);

    auto add_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), make_int_num(1)});
    auto expr = wrap_expr(add_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(),
                "Addition propagates finite numeric expression endpoints");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 3.0, 1e-10, "Lower bound is (1+1)+1=3");
        EXPECT_NEAR(get_upper(*result), 6.0, 1e-10, "Upper bound is (3+2)+1=6");
    }
}


void test_subtraction_propagation() {
    TEST_CASE("Subtraction: [1,5] - [2,3] = [-2,3]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", 2.0, 3.0);

    // x - y is represented as x + (-1)*y
    auto neg_y = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_int_num(-1), make_var("y")});
    auto sub_node = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), neg_y});
    auto expr = wrap_expr(sub_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Subtraction produces a bounded interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), -2.0, 1e-10, "Lower bound is 1-3=-2");
        EXPECT_NEAR(get_upper(*result), 3.0, 1e-10, "Upper bound is 5-2=3");
    }
}


void test_multiplication_propagation() {
    TEST_CASE("Multiplication: [2,3] * [4,5] = [8,15]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 3.0);
    declare_bounds(ctx, "y", 4.0, 5.0);

    auto mul_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Multiplication produces a bounded interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 8.0, 1e-10, "Lower bound is min(8,10,12,15)=8");
        EXPECT_NEAR(get_upper(*result), 15.0, 1e-10, "Upper bound is max(8,10,12,15)=15");
    }
}

void test_multiplication_mixed_signs() {
    TEST_CASE("Multiplication: [-2,3] * [1,4] = [-8,12]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -2.0, 3.0);
    declare_bounds(ctx, "y", 1.0, 4.0);

    auto mul_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Multiplication with mixed signs produces interval");
    if (result.has_value()) {
        // Products: (-2)*1=-2, (-2)*4=-8, 3*1=3, 3*4=12
        EXPECT_NEAR(get_lower(*result), -8.0, 1e-10, "Lower bound is min(-2,-8,3,12)=-8");
        EXPECT_NEAR(get_upper(*result), 12.0, 1e-10, "Upper bound is max(-2,-8,3,12)=12");
    }
}


void test_squaring_nonnegative() {
    TEST_CASE("Squaring: [2,5]^2 = [4,25] (a>=0)");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 5.0);

    auto pow_node = lamina::detail::make_node<PowerNode>(make_var("x"), make_int_num(2));
    auto expr = wrap_expr(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Squaring non-negative produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 4.0, 1e-10, "Lower bound is 2^2=4");
        EXPECT_NEAR(get_upper(*result), 25.0, 1e-10, "Upper bound is 5^2=25");
    }
}

void test_squaring_spanning_zero() {
    TEST_CASE("Squaring: [-3,2]^2 = [0,9] (spans zero)");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -3.0, 2.0);

    auto pow_node = lamina::detail::make_node<PowerNode>(make_var("x"), make_int_num(2));
    auto expr = wrap_expr(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Squaring spanning zero produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 0.0, 1e-10, "Lower bound is 0 (interval spans zero)");
        EXPECT_NEAR(get_upper(*result), 9.0, 1e-10, "Upper bound is max(9,4)=9");
    }
}

void test_squaring_nonpositive() {
    TEST_CASE("Squaring: [-5,-2]^2 = [4,25] (b<=0)");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -5.0, -2.0);

    auto pow_node = lamina::detail::make_node<PowerNode>(make_var("x"), make_int_num(2));
    auto expr = wrap_expr(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Squaring non-positive produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 4.0, 1e-10, "Lower bound is (-2)^2=4");
        EXPECT_NEAR(get_upper(*result), 25.0, 1e-10, "Upper bound is (-5)^2=25");
    }
}


void test_division_positive() {
    TEST_CASE("Division: [2,6] / [1,3] = [2/3, 6]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 6.0);
    declare_bounds(ctx, "y", 1.0, 3.0);

    // x / y is represented as x * y^(-1)
    auto y_inv = lamina::detail::make_node<PowerNode>(make_var("y"), make_int_num(-1));
    auto div_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), y_inv});
    auto expr = wrap_expr(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "Division by positive interval produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 2.0 / 3.0, 1e-10, "Lower bound is 2/3");
        EXPECT_NEAR(get_upper(*result), 6.0, 1e-10, "Upper bound is 6/1=6");
    }
}


void test_division_by_zero_containing() {
    TEST_CASE("Division: [1,5] / [-1,2] = unbounded (divisor contains zero)");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", -1.0, 2.0);

    // x / y = x * y^(-1)
    auto y_inv = lamina::detail::make_node<PowerNode>(make_var("y"), make_int_num(-1));
    auto div_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), y_inv});
    auto expr = wrap_expr(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE(result.has_value(), "Division by zero-containing interval returns nullopt");
}

void test_division_by_zero_at_boundary() {
    TEST_CASE("Division: [1,5] / [0,3] = unbounded (divisor contains zero at boundary)");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", 0.0, 3.0);

    auto y_inv = lamina::detail::make_node<PowerNode>(make_var("y"), make_int_num(-1));
    auto div_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x"), y_inv});
    auto expr = wrap_expr(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE(result.has_value(), "Division by interval containing zero at boundary returns nullopt");
}


void test_sin_propagation() {
    TEST_CASE("sin([0, 6.28]) = [-1, 1]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 0.0, 6.28);

    auto sin_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin,
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x")});
    auto expr = wrap_expr(sin_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "sin of bounded input produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), -1.0, 1e-10, "sin lower bound is -1");
        EXPECT_NEAR(get_upper(*result), 1.0, 1e-10, "sin upper bound is 1");
    }
}

void test_cos_propagation() {
    TEST_CASE("cos([0, 3.14]) = [-1, 1]");
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 0.0, 3.14);

    auto cos_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Cos,
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var("x")});
    auto expr = wrap_expr(cos_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "cos of bounded input produces interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), -1.0, 1e-10, "cos lower bound is -1");
        EXPECT_NEAR(get_upper(*result), 1.0, 1e-10, "cos upper bound is 1");
    }
}


void test_number_node_propagation() {
    TEST_CASE("NumberNode 5.0 produces [5, 5]");
    AssumptionContext ctx;

    auto expr = wrap_expr(make_num(5.0));

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_TRUE(result.has_value(), "NumberNode produces a point interval");
    if (result.has_value()) {
        EXPECT_NEAR(get_lower(*result), 5.0, 1e-10, "Lower bound is 5");
        EXPECT_NEAR(get_upper(*result), 5.0, 1e-10, "Upper bound is 5");
    }
}


void test_unbounded_variable() {
    TEST_CASE("Variable without bounds returns nullopt");
    AssumptionContext ctx;

    auto expr = wrap_expr(make_var("x"));

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE(result.has_value(), "Unbounded variable returns nullopt");
}


int main() {
    // Addition
    test_addition_propagation();
    test_addition_negative_bounds();
    test_numeric_expression_endpoint_propagation();

    // Subtraction
    test_subtraction_propagation();

    // Multiplication
    test_multiplication_propagation();
    test_multiplication_mixed_signs();

    // Squaring
    test_squaring_nonnegative();
    test_squaring_spanning_zero();
    test_squaring_nonpositive();

    // Division
    test_division_positive();
    test_division_by_zero_containing();
    test_division_by_zero_at_boundary();

    // Trig functions
    test_sin_propagation();
    test_cos_propagation();

    // Edge cases
    test_number_node_propagation();
    test_unbounded_variable();

    return TEST_REPORT();
}
