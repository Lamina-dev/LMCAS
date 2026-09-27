
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>

using namespace LMCAS;

static SymbolicExpr make_func_expr(FunctionNode::FuncType type,
                                   std::shared_ptr<const SymbolicNode> arg) {
    std::vector<std::shared_ptr<const SymbolicNode>> args = {std::move(arg)};
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<FunctionNode>(type, std::move(args)));
}

static std::shared_ptr<const SymbolicNode> var(const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> num(int v) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(v));
}

TEST(AssumptionInferenceTrig, SinRealArgRealDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sin(x) is Real when x is Real";
}

TEST(AssumptionInferenceTrig, SinRealArgSignUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "sin(x) Positive is Unknown (sin can be negative)";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "sin(x) Negative is Unknown (sin can be positive)";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "sin(x) NonNegative is Unknown (sin can be negative)";
}

TEST(AssumptionInferenceTrig, SinIntegerArg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("n"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sin(n) is Real when n is Integer";
}

TEST(AssumptionInferenceTrig, SinNoAssumption) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "sin(x) Real is Unknown when x has no assumptions";
}

TEST(AssumptionInferenceTrig, CosRealArgRealDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "cos(x) is Real when x is Real";
}

TEST(AssumptionInferenceTrig, CosRealArgSignUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "cos(x) Positive is Unknown (cos can be negative)";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "cos(x) Negative is Unknown (cos can be positive)";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "cos(x) NonNegative is Unknown (cos can be negative)";
}

TEST(AssumptionInferenceTrig, CosIntegerArg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("n"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "cos(n) is Real when n is Integer";
}

TEST(AssumptionInferenceTrig, CosNoAssumption) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "cos(x) Real is Unknown when x has no assumptions";
}

TEST(AssumptionInferenceTrig, SinNumericArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, num(1));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sin(1) is Real";
}

TEST(AssumptionInferenceTrig, CosNumericArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, num(0));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "cos(0) is Real";
}

TEST(AssumptionInferenceTrig, TanRealArgRealDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Tan, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "tan(x) is Real when x is Real";
}

TEST(AssumptionInferenceTrig, TanRealArgSignUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Tan, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "tan(x) Positive is Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "tan(x) Negative is Unknown";
}

TEST(AssumptionInferenceTrig, ArctanSignAndRealDomain) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_domain("x", Domain::Real).has_value())) << "实数域假设应成功";
    EXPECT_TRUE((context.assume_sign("x", Sign::Positive).has_value())) << "正号假设应成功";
    InferenceEngine engine(context);
    auto expression =
        make_func_expr(FunctionNode::FuncType::ArcTan, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expression).value() == Tribool::True)) << "正参数的 atan 为正";
    EXPECT_TRUE((engine.query_negative_checked(expression).value() == Tribool::False)) << "正参数的 atan 不为负";
    EXPECT_TRUE((engine.query_nonnegative_checked(expression).value() == Tribool::True)) << "正参数的 atan 非负";
    EXPECT_TRUE((engine.query_nonpositive_checked(expression).value() == Tribool::False)) << "正参数的 atan 不为非正";
    EXPECT_TRUE((engine.query_nonzero_checked(expression).value() == Tribool::True)) << "正参数的 atan 非零";
    EXPECT_TRUE((engine.query_real_checked(expression).value() == Tribool::True)) << "实参数的 atan 为实数";
}
