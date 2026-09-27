#include "expr.hpp"
#include "limit_result.hpp"
#include "numeric_evaluation.hpp"
#include "internal/symbolic_ast.hpp"
#include "test_common.hpp"

using namespace LMCAS;

static void expect_finite_limit(
    const char* source, const ExprPtr& point,
    LimitDirection direction, double expected)
{
    auto expression = parse_expr(source);
    ASSERT_TRUE(expression.has_value()) << source;
    auto result =
        limit_checked(expression.value(), "x", point, direction);
    ASSERT_TRUE(result.has_value()) << source;
    auto finite = std::get_if<FiniteLimit>(&result.value().value);
    ASSERT_TRUE(finite != nullptr) << source;
    auto value = evaluate_numeric(*finite->value);
    ASSERT_TRUE(value.has_value()) << source;
    EXPECT_TRUE(std::isfinite(value.value().value)) << source;
    EXPECT_NEAR(value.value().value, expected, 1e-6) << source;
}

TEST(LmcasIndeterminateLimits, RealDefinednessProtection) {
    using Node = std::shared_ptr<const SymbolicNode>;
    auto zero = detail::node(SymbolicExpr::number(0));
    auto negative = detail::node(SymbolicExpr::number(-1));
    auto half = detail::node(SymbolicExpr::number(Rational(1, 2)));
    const std::vector<Node> invalid_values{
        detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<Node>{negative}),
        detail::make_node<PowerNode>(negative, half),
        detail::make_node<PowerNode>(zero, zero)};
    for (const auto &value : invalid_values) {
        auto result = limit_checked(detail::make_expression_ptr(value), "x", SymbolicExpr::number(1));
        EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << ("default real limits reject undefined constants without complex continuation");
        auto endpoint = limit_checked(SymbolicExpr::variable("x"), "x", detail::make_expression_ptr(value));
        EXPECT_TRUE((!endpoint && endpoint.error().code == CasErrc::Inconclusive)) << ("default real limits reject undefined endpoints");
    }
    auto invalid_neighborhood = parse_expr("sqrt(-x)");
    ASSERT_TRUE(invalid_neighborhood.has_value())
        << "invalid real neighborhood fixture parses";
    auto result = limit_checked(
        invalid_neighborhood.value(), "x", SymbolicExpr::number(1));
    EXPECT_TRUE(
        !result && result.error().code == CasErrc::Inconclusive)
        << "real neighborhood obligations remain required";
}

TEST(LmcasIndeterminateLimits, ZeroTimesLogarithm) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("x*ln(x)", zero, LimitDirection::FromAbove, 0.0);
}

TEST(LmcasIndeterminateLimits, ExponentialDecay) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("x*exp(-x)", infinity, LimitDirection::Both, 0.0);
}

TEST(LmcasIndeterminateLimits, ReciprocalDifference) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("1/x-1/sin(x)", zero, LimitDirection::Both, 0.0);
}

TEST(LmcasIndeterminateLimits, RadicalDifferenceAtInfinity) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("x-sqrt(x^2+x)", infinity, LimitDirection::Both, -0.5);
}

TEST(LmcasIndeterminateLimits, ExponentialPower) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("(1+1/x)^x", infinity, LimitDirection::Both, std::exp(1.0));
}

TEST(LmcasIndeterminateLimits, ZeroToZeroPower) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("x^x", zero, LimitDirection::FromAbove, 1.0);
}

TEST(LmcasIndeterminateLimits, PowerAtOne) {
    expect_finite_limit("x^x", SymbolicExpr::number(1), LimitDirection::Both, 1.0);
}

TEST(LmcasIndeterminateLimits, InfiniteToZeroPower) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("x^(1/x)", infinity, LimitDirection::Both, 1.0);
}

TEST(LmcasIndeterminateLimits, CubicSineResidual) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("(sin(x)-x)/x^3", zero, LimitDirection::Both, -1.0 / 6.0);
}

TEST(LmcasIndeterminateLimits, SineRatio) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("sin(x)/x", zero, LimitDirection::Both, 1.0);
}

TEST(LmcasIndeterminateLimits, ExponentialRatio) {
    auto zero = SymbolicExpr::number(0);
    expect_finite_limit("(exp(x)-1)/x", zero, LimitDirection::Both, 1.0);
}

TEST(LmcasIndeterminateLimits, ScaledRadicalDifference) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("sqrt(4*x^2+3*x)-2*x", infinity, LimitDirection::Both, 0.75);
}

TEST(LmcasIndeterminateLimits, TwoRadicalDifference) {
    auto infinity = SymbolicExpr::infinity();
    expect_finite_limit("sqrt(x^2+3*x)-sqrt(x^2+x)", infinity, LimitDirection::Both, 1.0);
}
