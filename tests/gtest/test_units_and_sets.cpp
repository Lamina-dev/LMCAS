#include <gtest/gtest.h>

#include "expr.hpp"
#include "symbolic_ast.hpp"

using namespace LMCAS;

namespace LMCAS {
namespace {

using namespace LMCAS;

ExprPtr exact(long long value) {
    auto result = integer(value);
    EXPECT_TRUE(result);
    return result.value();
}

TEST(UnitSystemTest, UsesExactDimensionsAndConversionFactors) {
    ComputationContext context;
    auto factor = context.units().conversion_factor("km", "m");
    ASSERT_TRUE(factor);
    EXPECT_EQ(factor.value(), Rational(1000));

    auto mismatch = context.units().conversion_factor("m", "s");
    ASSERT_FALSE(mismatch);
    EXPECT_EQ(mismatch.error().code, CasErrc::DimensionMismatch);

    ASSERT_TRUE(context.units().declare_base_unit("widget"));
    auto duplicate = context.units().declare_base_unit("widget");
    ASSERT_FALSE(duplicate);
}

TEST(QuantityTest, ConvertsAndStripsWithoutApproximation) {
    ComputationContext context;
    auto distance = with_unit(exact(2), "km", context);
    ASSERT_TRUE(distance);
    auto metres = convert_to_unit(distance.value(), "m", context);
    ASSERT_TRUE(metres);

    auto stripped = strip_to_base_value(metres.value(), context);
    ASSERT_TRUE(stripped);
    auto simplified = stripped.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), "2000");

    auto display = strip_to_display_value(distance.value(), context);
    ASSERT_TRUE(display);
    EXPECT_EQ(display.value()->to_string(), "2");

    auto one_metre = with_unit(exact(1), "m", context);
    auto hundred_centimetres = with_unit(exact(100), "cm", context);
    ASSERT_TRUE(one_metre);
    ASSERT_TRUE(hundred_centimetres);
    auto same_length = equivalent_core(*one_metre.value(),
                                       *hundred_centimetres.value(), context);
    ASSERT_TRUE(same_length);
    EXPECT_TRUE(same_length.value());
}

TEST(QuantityTest, AcceptsResolvedCustomUnits) {
    ComputationContext context;
    UnitDefinition level{DimensionSignature::base("user::score"), Rational(100)};
    UnitDefinition score{DimensionSignature::base("user::score"), Rational(1)};

    auto points = with_unit_definition(exact(3), "level", level, context);
    ASSERT_TRUE(points);
    auto converted = convert_to_unit_definition(
        points.value(), "score", score, context);
    ASSERT_TRUE(converted);

    auto display = strip_to_display_value(converted.value(), context);
    ASSERT_TRUE(display);
    auto simplified = display.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), "300");

    auto mismatch = convert_to_unit_definition(
        points.value(), "second",
        UnitDefinition{DimensionSignature::base("s"), Rational(1)}, context);
    ASSERT_FALSE(mismatch);
    EXPECT_EQ(mismatch.error().code, CasErrc::DimensionMismatch);
}

TEST(QuantityTest, EnforcesDimensionRules) {
    ComputationContext context;
    auto metre = with_unit(exact(1), "m", context);
    auto second = with_unit(exact(1), "s", context);
    ASSERT_TRUE(metre);
    ASSERT_TRUE(second);

    auto invalid_sum = add(metre.value(), second.value(), context);
    ASSERT_FALSE(invalid_sum);
    EXPECT_EQ(invalid_sum.error().code, CasErrc::DimensionMismatch);

    auto speed = div(metre.value(), second.value(), context);
    ASSERT_TRUE(speed);
    auto dimension = dimension_of(*speed.value());
    ASSERT_TRUE(dimension);
    DimensionSignature expected = DimensionSignature::base("m").divided_by(
        DimensionSignature::base("s"));
    EXPECT_EQ(dimension.value(), expected);

    auto invalid_sine = sin(metre.value(), context);
    ASSERT_FALSE(invalid_sine);
    EXPECT_EQ(invalid_sine.error().code, CasErrc::DimensionMismatch);

    auto square = mul(metre.value(), metre.value(), context);
    ASSERT_TRUE(square);
    auto root = sqrt(square.value(), context);
    ASSERT_TRUE(root);
    auto root_dimension = dimension_of(*root.value());
    ASSERT_TRUE(root_dimension);
    EXPECT_EQ(root_dimension.value(), DimensionSignature::base("m"));
}

TEST(FiniteSetTest, IsUnorderedDeduplicatedAndSupportsAlgebra) {
    ComputationContext context;
    auto left = finite_set({exact(2), exact(1), exact(1)}, context);
    auto same = finite_set({exact(1), exact(2)}, context);
    ASSERT_TRUE(left);
    ASSERT_TRUE(same);
    EXPECT_TRUE(structurally_equal(*left.value(), *same.value()));
    EXPECT_EQ(detail::node(left.value())->hash(), detail::node(same.value())->hash());

    auto right = finite_set({exact(2), exact(3)}, context);
    ASSERT_TRUE(right);
    auto intersection = finite_set_intersection(*left.value(), *right.value(), context);
    ASSERT_TRUE(intersection);
    EXPECT_EQ(intersection.value()->to_string(), "{2}");

    auto membership = member(exact(1), left.value(), context);
    ASSERT_TRUE(membership);
    EXPECT_EQ(membership.value()->to_string(), "1");
}

TEST(IntervalTest, HonorsOpenAndClosedEndpoints) {
    ComputationContext context;
    auto range = interval(exact(0), exact(1), false, true, context);
    ASSERT_TRUE(range);
    EXPECT_EQ(range.value()->to_string(), "(0, 1]");

    auto lower = member(exact(0), range.value(), context);
    auto upper = member(exact(1), range.value(), context);
    ASSERT_TRUE(lower);
    ASSERT_TRUE(upper);
    EXPECT_EQ(lower.value()->to_string(), "0");
    EXPECT_EQ(upper.value()->to_string(), "1");

    auto reversed = interval(exact(2), exact(1), true, true, context);
    ASSERT_FALSE(reversed);
    EXPECT_EQ(reversed.error().code, CasErrc::DomainError);
}

} // namespace
} // namespace LMCAS
