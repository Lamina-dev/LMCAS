
#include "test_assumption_interval_support.hpp"
#include <cmath>
#include <limits>
#include <vector>

/// Declare bounded interval [lo, hi] for a variable in the context
static void declare_bounds(AssumptionContext &ctx, const std::string &var, double lo, double hi) {
    auto lower_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    EXPECT_TRUE((ctx.current_properties().declare_bounded(var, Boundedness::Bounded, bounds).has_value())) << "boundedness declaration succeeds";
}

/// Extract numeric lower bound from an interval
static double get_lower(const Interval &iv) {
    if (iv.lower.is_neg_infinity)
        return -std::numeric_limits<double>::infinity();
    if (iv.lower.value)
        return iv.lower.value->to_numeric();
    return 0.0;
}

/// Extract numeric upper bound from an interval
static double get_upper(const Interval &iv) {
    if (iv.upper.is_pos_infinity)
        return std::numeric_limits<double>::infinity();
    if (iv.upper.value)
        return iv.upper.value->to_numeric();
    return 0.0;
}

TEST(AssumptionInterval, AdditionPropagation) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 3.0);
    declare_bounds(ctx, "y", 2.0, 5.0);

    // x + y
    auto add_node = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(add_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Addition produces a bounded interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {3.0, 8.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, AdditionNegativeBounds) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -2.0, 1.0);
    declare_bounds(ctx, "y", -3.0, 4.0);

    auto add_node = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(add_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Addition with negatives produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {-5.0, 5.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, SubtractionPropagation) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", 2.0, 3.0);

    // x - y is represented as x + (-1)*y
    auto neg_y = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(-1), test_variable_node("y")});
    auto sub_node = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), neg_y});
    auto expr = test_expression_from_node(sub_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Subtraction produces a bounded interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {-2.0, 3.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, MultiplicationPropagation) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 3.0);
    declare_bounds(ctx, "y", 4.0, 5.0);

    auto mul_node = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Multiplication produces a bounded interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {8.0, 15.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, MultiplicationMixedSigns) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -2.0, 3.0);
    declare_bounds(ctx, "y", 1.0, 4.0);

    auto mul_node = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Multiplication with mixed signs produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {-8.0, 12.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, SquaringNonnegative) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 5.0);

    auto pow_node = LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2));
    auto expr = test_expression_from_node(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Squaring non-negative produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {4.0, 25.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, SquaringSpanningZero) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -3.0, 2.0);

    auto pow_node = LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2));
    auto expr = test_expression_from_node(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Squaring spanning zero produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {0.0, 9.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, SquaringNonpositive) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", -5.0, -2.0);

    auto pow_node = LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2));
    auto expr = test_expression_from_node(pow_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Squaring non-positive produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {4.0, 25.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, DivisionPositive) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 2.0, 6.0);
    declare_bounds(ctx, "y", 1.0, 3.0);

    // x / y is represented as x * y^(-1)
    auto y_inv = LMCAS::detail::make_node<PowerNode>(test_variable_node("y"), test_integer_node(-1));
    auto div_node = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), y_inv});
    auto expr = test_expression_from_node(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "Division by positive interval produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {2.0 / 3.0, 6.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, DivisionByZeroContaining) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", -1.0, 2.0);

    // x / y = x * y^(-1)
    auto y_inv = LMCAS::detail::make_node<PowerNode>(test_variable_node("y"), test_integer_node(-1));
    auto div_node = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), y_inv});
    auto expr = test_expression_from_node(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE((result.has_value())) << "Division by zero-containing interval returns nullopt";
}

TEST(AssumptionInterval, DivisionByZeroAtBoundary) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 1.0, 5.0);
    declare_bounds(ctx, "y", 0.0, 3.0);

    auto y_inv = LMCAS::detail::make_node<PowerNode>(test_variable_node("y"), test_integer_node(-1));
    auto div_node = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), y_inv});
    auto expr = test_expression_from_node(div_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE((result.has_value())) << "Division by interval containing zero at boundary returns nullopt";
}

TEST(AssumptionInterval, SinPropagation) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 0.0, 6.28);

    auto sin_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin,
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x")});
    auto expr = test_expression_from_node(sin_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "sin of bounded input produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {-1.0, 1.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, CosPropagation) {
    AssumptionContext ctx;
    declare_bounds(ctx, "x", 0.0, 3.14);

    auto cos_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Cos,
        std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x")});
    auto expr = test_expression_from_node(cos_node);

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "cos of bounded input produces interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {-1.0, 1.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, NumberNodePropagation) {
    AssumptionContext ctx;

    auto expr = test_expression_from_node(make_num(5.0));

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    ASSERT_TRUE(result.has_value()) << "NumberNode produces a point interval";
    const double actual_bounds[] = {get_lower(*result), get_upper(*result)};
    const double expected_bounds[] = {5.0, 5.0};
    for (size_t endpoint = 0; endpoint < 2; ++endpoint) {
        SCOPED_TRACE(endpoint);
        EXPECT_TRUE(std::isfinite(actual_bounds[endpoint]));
        EXPECT_NEAR(actual_bounds[endpoint], expected_bounds[endpoint], 1e-10);
    }
}

TEST(AssumptionInterval, UnboundedVariable) {
    AssumptionContext ctx;

    auto expr = test_expression_from_node(test_variable_node("x"));

    InferenceEngine engine(ctx);
    auto result = engine.propagate_bounds(expr);

    EXPECT_FALSE((result.has_value())) << "Unbounded variable returns nullopt";
}

TEST(AssumptionInterval, ArctanFiniteBounds) {
    AssumptionContext bounded_context;
    EXPECT_TRUE((bounded_context.assume_domain("x", Domain::Real).has_value())) << "有限区间实数域假设应成功";
    declare_bounds(bounded_context, "x", -1.0, 1.0);
    auto bounded_expression = test_expression_from_node(LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcTan,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            test_variable_node("x")}));
    InferenceEngine bounded_engine(bounded_context);
    auto bounded = bounded_engine.propagate_bounds(bounded_expression);
    ASSERT_TRUE((bounded.has_value())) << "有限参数区间应传播到 atan";
    if (bounded) {
        {
            const double actual_value = (get_lower(*bounded));
            const double expected_value = (-std::atan(1.0));
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
        {
            const double actual_value = (get_upper(*bounded));
            const double expected_value = (std::atan(1.0));
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(AssumptionInterval, ArctanGlobalBounds) {
    auto bounded_expression = test_expression_from_node(LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcTan,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            test_variable_node("x")}));

    AssumptionContext real_context;
    EXPECT_TRUE((real_context.assume_domain("x", Domain::Real).has_value())) << "无界实数域假设应成功";
    InferenceEngine real_engine(real_context);
    auto global = real_engine.propagate_bounds(bounded_expression);
    ASSERT_TRUE((global.has_value())) << "实参数 atan 应具有全局值域";
    if (global) {
        EXPECT_TRUE((global->lower.is_open && global->upper.is_open)) << "atan 全局值域两端均为开端点";
        {
            const double actual_value = (get_lower(*global));
            const double expected_value = (-LMMC_CONST_PI / 2.0);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
        {
            const double actual_value = (get_upper(*global));
            const double expected_value = (LMMC_CONST_PI / 2.0);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}
