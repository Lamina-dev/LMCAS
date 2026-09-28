/**
 * @file test_calculus_utils_geometry.cpp
 * @brief 曲率、曲面面积与拐点的行为约定。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"
#include "expr.hpp"
#include <string>
#include <variant>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string &name) { return SE::variable(name); }
static std::shared_ptr<SymbolicExpr> bigint_num(const BigInt &n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(
            std::variant<BigInt, Rational, lmmc_real_t>{
                std::in_place_type<BigInt>, n}));
}

TEST(CalculusUtilsGeometry, CurvatureInvalidInputs) {
    auto null_result = LMCAS::curvature_checked(nullptr, "x");
    EXPECT_TRUE((!null_result &&
                 null_result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked curvature rejects null expression";

    auto empty_var = LMCAS::curvature_checked(var("x"), "");
    EXPECT_TRUE((!empty_var &&
                 empty_var.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked curvature rejects empty variable";
}

TEST(CalculusUtilsGeometry, ParametricZeroVelocity) {
    auto checked = LMCAS::curvature_parametric_checked(num(1), num(2), "t");
    EXPECT_TRUE((!checked &&
                 checked.error().code == LMCAS::CasErrc::DomainError))
        << "checked parametric curvature rejects zero velocity";
}

TEST(CalculusUtilsGeometry, CurvatureDerivativeErrors) {
    auto floor_result = LMCAS::floor(var("x"));
    ASSERT_TRUE(floor_result);
    const auto& floor_expression = floor_result.value();
    auto integral = detail::make_expression_ptr(
        detail::make_node<IntegralNode>(detail::node(floor_expression), "x"));
    for (const auto& expression : {floor_expression, integral}) {
        auto explicit_result = curvature_checked(expression, "x");
        ASSERT_FALSE(explicit_result);
        EXPECT_EQ(explicit_result.error().code, CasErrc::UnsupportedExpression);
        EXPECT_EQ(explicit_result.error().operation, "curvature");
        for (bool first_coordinate : {true, false}) {
            auto parametric = curvature_parametric_checked(
                first_coordinate ? expression : var("x"),
                first_coordinate ? var("x") : expression, "x");
            ASSERT_FALSE(parametric);
            EXPECT_EQ(parametric.error().code, CasErrc::UnsupportedExpression);
            EXPECT_EQ(parametric.error().operation, "curvature_parametric");
        }
    }
}

TEST(CalculusUtilsGeometry, InflectionDerivativeErrors) {
    auto floor_result = LMCAS::floor(var("x"));
    ASSERT_TRUE(floor_result);
    const auto& floor_expression = floor_result.value();
    auto integral = detail::make_expression_ptr(
        detail::make_node<IntegralNode>(detail::node(floor_expression), "x"));
    for (const auto& expression : {floor_expression, integral}) {
        auto result = inflection_points_checked(expression, "x");
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().code, CasErrc::UnsupportedExpression);
        EXPECT_EQ(result.error().operation, "inflection_points");
    }
}

TEST(CalculusUtilsGeometry, DerivativeResourceLimit) {
    auto x = var("x");
    auto square = SE::power(x, num(2));
    ResourceLimits limits;
    limits.max_steps = 2;
    ComputationContext explicit_context(limits);
    auto explicit_result = curvature_checked(square, "x", explicit_context);
    ASSERT_FALSE(explicit_result);
    EXPECT_EQ(explicit_result.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(explicit_result.error().operation, "curvature");

    limits.max_steps = 4;
    ComputationContext parametric_context(limits);
    auto parametric = curvature_parametric_checked(
        x, square, "x", parametric_context);
    ASSERT_FALSE(parametric);
    EXPECT_EQ(parametric.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(parametric.error().operation, "curvature_parametric");

    limits.max_steps = 2;
    ComputationContext inflection_context(limits);
    auto inflection = inflection_points_checked(
        SE::power(x, num(3)), "x", inflection_context);
    ASSERT_FALSE(inflection);
    EXPECT_EQ(inflection.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(inflection.error().operation, "inflection_points");
}

TEST(CalculusUtilsGeometry, DerivativeCancellation) {
    CancellationToken cancellation;
    cancellation.cancel();
    auto x = var("x");
    auto square = SE::power(x, num(2));
    ComputationContext explicit_context({}, cancellation);
    auto explicit_result = curvature_checked(square, "x", explicit_context);
    ASSERT_FALSE(explicit_result);
    EXPECT_EQ(explicit_result.error().code, CasErrc::Cancelled);
    EXPECT_EQ(explicit_result.error().operation, "curvature");

    ComputationContext parametric_context({}, cancellation);
    auto parametric = curvature_parametric_checked(
        x, square, "x", parametric_context);
    ASSERT_FALSE(parametric);
    EXPECT_EQ(parametric.error().code, CasErrc::Cancelled);
    EXPECT_EQ(parametric.error().operation, "curvature_parametric");

    ComputationContext inflection_context({}, cancellation);
    auto inflection = inflection_points_checked(
        SE::power(x, num(3)), "x", inflection_context);
    ASSERT_FALSE(inflection);
    EXPECT_EQ(inflection.error().code, CasErrc::Cancelled);
    EXPECT_EQ(inflection.error().operation, "inflection_points");
}

TEST(CalculusUtilsGeometry, SurfaceHugeBounds) {
    auto x = var("x");
    auto unsupported = SE::exp(SE::power(x, num(2)));
    const BigInt huge("1" + std::string(400, '0'));
    auto result = LMCAS::surface_area_revolution_x_checked(
        unsupported, "x", bigint_num(BigInt(0)), bigint_num(huge));
    EXPECT_TRUE((!result)) << "huge exact bounds do not fabricate a surface area";
}

TEST(CalculusUtilsGeometry, SurfaceCheckedContracts) {
    auto x = var("x");
    auto zero = num(0);
    auto one = num(1);

    auto exact = LMCAS::surface_area_revolution_x_checked(num(1), "x", zero, one);
    ASSERT_TRUE((exact.has_value())) << "checked x-axis surface area succeeds for constant radius";
    auto evaluated = LMCAS::evalf(*exact.value());
    ASSERT_TRUE(evaluated);
    EXPECT_TRUE(evaluated.value().is_finite());
    EXPECT_NEAR(evaluated.value().value, 2.0 * std::acos(-1.0), 1e-9)
        << "checked constant surface area is 2*pi";

    auto null_result = LMCAS::surface_area_revolution_x_checked(
        nullptr, "x", zero, one);
    EXPECT_TRUE((!null_result &&
                 null_result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked surface area rejects null expression";

    auto null_bound = LMCAS::surface_area_revolution_y_checked(
        x, "x", zero, nullptr);
    EXPECT_TRUE((!null_bound &&
                 null_bound.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked surface area rejects null bounds";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::surface_area_revolution_x_checked(
        x, "x", zero, one, limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked surface area observes exhausted step budget";
}

TEST(CalculusUtilsGeometry, SurfaceNumericDomain) {
    auto x = var("x");
    auto zero = num(0);
    auto one = num(1);
    auto unsupported = SE::exp(SE::power(x, num(2)));

    auto checked = LMCAS::surface_area_revolution_x_checked(
        unsupported, "x", zero, one);
    EXPECT_TRUE((!checked &&
                 checked.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked surface area rejects implicit numeric fallback";
}

TEST(CalculusUtilsGeometry, InflectionInvalidInputs) {
    auto null_result = LMCAS::inflection_points_checked(nullptr, "x");
    EXPECT_TRUE((!null_result &&
                 null_result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked inflection points rejects null expression";

    auto empty_var = LMCAS::inflection_points_checked(var("x"), "");
    EXPECT_TRUE((!empty_var &&
                 empty_var.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked inflection points rejects empty variable";
}

TEST(CalculusUtilsGeometry, InflectionOutcomes) {
    auto x = var("x");

    auto no_inflection = LMCAS::inflection_points_checked(SE::power(x, num(2)), "x");
    EXPECT_TRUE((no_inflection && no_inflection.value().empty())) << "checked inflection points reports an empty set for f'' = constant nonzero";

    auto cubic = LMCAS::inflection_points_checked(SE::power(x, num(3)), "x");
    EXPECT_TRUE((cubic && cubic.value().size() == 1)) << "checked inflection points reports finite exact candidates for x^3";

    auto exponential = LMCAS::inflection_points_checked(SE::exp(x), "x");
    EXPECT_TRUE((exponential && exponential.value().empty())) << "checked inflection points proves exp(x) has no real inflection point";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::inflection_points_checked(
        SE::power(x, num(3)), "x", limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked inflection points propagates solver budget failure";
}
