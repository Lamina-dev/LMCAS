#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprValues, ResolvedUnitConversion) {
    LMCAS::ComputationContext unit_context;
    LMCAS::UnitDefinition level{
        LMCAS::DimensionSignature::base("user::score"), Rational(100)};
    LMCAS::UnitDefinition score{
        LMCAS::DimensionSignature::base("user::score"), Rational(1)};
    auto three = LMCAS::integer(3);
    auto points = LMCAS::with_unit_definition(
        three.value(), "level", level, unit_context);
    EXPECT_TRUE((points.has_value())) << "attach resolved unit";
    auto converted = LMCAS::convert_to_unit_definition(
        points.value(), "score", score, unit_context);
    EXPECT_TRUE((converted.has_value())) << "convert resolved unit";
    auto display = LMCAS::strip_to_display_value(
        converted.value(), unit_context);
    EXPECT_TRUE((display && display.value()->simplify()->to_string() == "300")) << "preserve target magnitude";
}

TEST(ExprValues, PiConstant) {
    auto constant = LMCAS::pi();
    ASSERT_TRUE(constant.has_value()) << "std.math.pi Expr can be constructed";
    auto value = LMCAS::evalf(*constant.value());
    ASSERT_TRUE(value.has_value());
    EXPECT_TRUE(value.value().is_finite()) << "std.math.pi explicitly evaluates through evalf";
    EXPECT_TRUE(std::isfinite(value.value().value));
    EXPECT_TRUE(std::isfinite(LMMC_CONST_PI));
    EXPECT_NEAR(value.value().value, LMMC_CONST_PI, 1e-15);
}

TEST(ExprValues, EulerConstant) {
    auto constant = LMCAS::e();
    ASSERT_TRUE(constant.has_value()) << "std.math.e Expr can be constructed";
    auto value = LMCAS::evalf(*constant.value());
    ASSERT_TRUE(value.has_value());
    EXPECT_TRUE(value.value().is_finite()) << "std.math.e explicitly evaluates through evalf";
    const double expected = std::exp(1.0);
    EXPECT_TRUE(std::isfinite(value.value().value));
    EXPECT_TRUE(std::isfinite(expected));
    EXPECT_NEAR(value.value().value, expected, 1e-15);
}

TEST(ExprValues, GoldenRatioConstant) {
    auto constant = LMCAS::phi();
    ASSERT_TRUE(constant.has_value()) << "std.math.phi Expr can be constructed";
    auto value = LMCAS::evalf(*constant.value());
    ASSERT_TRUE(value.has_value());
    EXPECT_TRUE(value.value().is_finite()) << "std.math.phi explicitly evaluates through evalf";
    const double expected = (1.0 + std::sqrt(5.0)) / 2.0;
    EXPECT_TRUE(std::isfinite(value.value().value));
    EXPECT_TRUE(std::isfinite(expected));
    EXPECT_NEAR(value.value().value, expected, 1e-15);
}

TEST(ExprValues, UnicodePiVariableIsNotConstant) {
    auto symbol = SymbolicExpr::variable("\xCF\x80");
    auto value = LMCAS::evalf(*symbol);
    ASSERT_FALSE(value);
    EXPECT_EQ(value.error().code, CasErrc::UnboundSymbol);
    auto constant = LMCAS::pi();
    ASSERT_TRUE(constant);
    EXPECT_FALSE(LMCAS::structurally_equal(*symbol, *constant.value()));
}

TEST(ExprValues, ApproximateRealValues) {
    auto approximate_half = LMCAS::approx_real(0.5);
    EXPECT_TRUE((approximate_half.has_value())) << "approx_real constructs an explicit approximate Expr";
    auto approximate_half_value =
        approximate_half ? LMCAS::evalf(*approximate_half.value())
                         : LMCAS::Result<LMCAS::ApproxReal>::failure(
                               LMCAS::CasErrc::InternalInvariant,
                               "approx_real construction failed", "test");
    EXPECT_TRUE((approximate_half_value &&
                 approximate_half_value.value().is_finite()))
        << "explicit approximate Expr can be evaluated with evalf";
    EXPECT_NEAR(approximate_half_value.value().value, 0.5, 0.0) << "approx_real preserves the requested finite value";

    auto nan_approx = LMCAS::approx_real(NAN);
    auto inf_approx = LMCAS::approx_real(INFINITY);
    EXPECT_TRUE((!nan_approx &&
                 nan_approx.error().code == LMCAS::CasErrc::InvalidArgument))
        << "approx_real rejects NaN";
    EXPECT_TRUE((!inf_approx &&
                 inf_approx.error().code == LMCAS::CasErrc::InvalidArgument))
        << "approx_real rejects infinity";
}

} // namespace
