#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

TEST(SymbolicVectorGeometryVector, DotOrthogonal) {
    auto one = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    std::vector<std::shared_ptr<SymbolicExpr>> a = {one, zero, zero};
    std::vector<std::shared_ptr<SymbolicExpr>> b = {zero, one, zero};
    auto result = LMCAS::vector_dot(a, b);
    auto simplified = result ? result->simplify() : nullptr;
    EXPECT_TRUE(test_expression_text((simplified), ("0"))) << "orthogonal dot product is 0";
}

TEST(SymbolicVectorGeometryVector, DotParallel) {
    auto result = LMCAS::vector_dot(
        {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
        {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)});
    auto simplified = result ? result->simplify() : nullptr;
    EXPECT_TRUE(test_expression_text((simplified), ("28"))) << "parallel dot product is 28";
}

TEST(SymbolicVectorGeometryVector, DotGeneral) {
    auto result = LMCAS::vector_dot(
        {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(2)},
        {SymbolicExpr::number(3), SymbolicExpr::number(4), SymbolicExpr::number(1)});
    auto simplified = result ? result->simplify() : nullptr;
    EXPECT_TRUE(test_expression_text((simplified), ("5"))) << "general dot product is 5";
}

TEST(SymbolicVectorGeometryVector, CrossParallel) {
    auto result = LMCAS::vector_cross(
        {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
        {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)});
    EXPECT_TRUE((result.size() == 3)) << "cross product has 3 components";
    for (int i = 0; i < 3; ++i) {
        auto s = result[i] ? result[i]->simplify() : nullptr;
        EXPECT_TRUE(test_expression_text((s), ("0"))) << "parallel cross component " + std::to_string(i) + " is 0";
    }
}

TEST(SymbolicVectorGeometryVector, CrossBasis) {
    auto one = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    auto result = LMCAS::vector_cross(
        {one, zero, zero},
        {zero, one, zero});
    ASSERT_EQ(result.size(), 3u) << "i x j has 3 components";
    const char *expected[] = {"0", "0", "1"};
    for (size_t component = 0; component < result.size(); ++component) {
        SCOPED_TRACE(component);
        ASSERT_TRUE(result[component]);
        const auto simplified = result[component]->simplify();
        ASSERT_TRUE(simplified);
        EXPECT_EQ(simplified->to_string(), expected[component]) << "i x j = k";
    }
}

TEST(SymbolicVectorGeometryVector, CrossGeneral) {
    auto result = LMCAS::vector_cross(
        {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
        {SymbolicExpr::number(4), SymbolicExpr::number(5), SymbolicExpr::number(6)});
    ASSERT_EQ(result.size(), 3u) << "general cross has 3 components";
    const char *expected[] = {"-3", "6", "-3"};
    for (size_t component = 0; component < result.size(); ++component) {
        SCOPED_TRACE(component);
        ASSERT_TRUE(result[component]);
        const auto simplified = result[component]->simplify();
        ASSERT_TRUE(simplified);
        EXPECT_EQ(simplified->to_string(), expected[component])
            << "(1,2,3) x (4,5,6) = (-3,6,-3)";
    }
}

TEST(SymbolicVectorGeometryVector, AngleOrthogonal) {
    auto one = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    double angle = LMCAS::vector_angle_checked({one, zero, zero}, {zero, one, zero}).value();
    EXPECT_TRUE((std::abs(angle - M_PI / 2.0) < 1e-9)) << "orthogonal angle is pi/2";
}

TEST(SymbolicVectorGeometryVector, AngleParallel) {
    double angle = LMCAS::vector_angle_checked({SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)}, {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}).value();
    EXPECT_TRUE((std::abs(angle - 0.0) < 1e-9)) << "parallel angle is 0";
}

TEST(SymbolicVectorGeometryVector, AngleExtremeScale) {
    const auto large = SymbolicExpr::number(1.0e200);
    const auto negative_large = SymbolicExpr::number(-1.0e200);
    auto angle = LMCAS::vector_angle_checked(
        {large, large}, {large, negative_large});
    ASSERT_TRUE((angle.has_value())) << "extreme-scale finite vectors have a defined angle";
    if (angle) {
        EXPECT_TRUE((std::isfinite(angle.value()))) << "extreme-scale vector angle remains finite";
        EXPECT_TRUE((std::abs(angle.value() - M_PI / 2.0) < 1e-12)) << "extreme-scale orthogonal vectors retain pi/2 angle";
    }
}

TEST(SymbolicVectorGeometryVector, DotCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto dot = LMCAS::vector_dot_checked({one, zero}, {zero, one});
    ASSERT_TRUE((dot.has_value())) << "checked vector_dot succeeds";
    if (dot) {
        EXPECT_TRUE(test_expression_text((dot.value()->simplify()), ("0"))) << "checked dot product of orthogonal vectors is 0";
    }
}

TEST(SymbolicVectorGeometryVector, DotCheckedDimension) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto bad_dot = LMCAS::vector_dot_checked({one}, {one, zero});
    EXPECT_TRUE((!bad_dot &&
                 bad_dot.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked vector_dot rejects dimension mismatch";
    bool legacy_dot_threw = false;
    try {
        (void)LMCAS::vector_dot({one}, {one, zero});
    } catch (const std::invalid_argument &) {
        legacy_dot_threw = true;
    }
    EXPECT_TRUE((legacy_dot_threw)) << "legacy vector_dot preserves invalid_argument for dimension mismatch";
}

TEST(SymbolicVectorGeometryVector, CrossCheckedComponent) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto bad_component = LMCAS::vector_cross_checked({one, nullptr, zero},
                                                     {zero, one, zero});
    EXPECT_TRUE((!bad_component &&
                 bad_component.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked vector_cross rejects null component";
}

TEST(SymbolicVectorGeometryVector, CrossCheckedDimension) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto bad_cross_dim = LMCAS::vector_cross_checked({one, zero},
                                                     {zero, one});
    EXPECT_TRUE((!bad_cross_dim &&
                 bad_cross_dim.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked vector_cross rejects non-3D vectors";
}

TEST(SymbolicVectorGeometryVector, AngleCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto angle = LMCAS::vector_angle_checked({one, zero}, {zero, one});
    ASSERT_TRUE((angle.has_value())) << "checked vector_angle succeeds";
    if (angle) {
        EXPECT_TRUE((std::abs(angle.value() - M_PI / 2.0) < 1e-9)) << "checked vector_angle of orthogonal vectors is pi/2";
    }
}

TEST(SymbolicVectorGeometryVector, AngleCheckedZero) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto zero_angle = LMCAS::vector_angle_checked({zero, zero}, {one, zero});
    EXPECT_TRUE((!zero_angle &&
                 zero_angle.error().code == LMCAS::CasErrc::DomainError))
        << "checked vector_angle rejects zero-length vectors";
}

TEST(SymbolicVectorGeometryVector, AngleCheckedSymbolic) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    auto symbolic_angle = LMCAS::vector_angle_checked({x, zero}, {one, zero});
    EXPECT_TRUE((!symbolic_angle &&
                 symbolic_angle.error().code == LMCAS::CasErrc::NumericFailure))
        << "checked vector_angle rejects symbolic components";
}

TEST(SymbolicVectorGeometryVector, AngleCheckedExpression) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto two_as_expr = SymbolicExpr::add(one, one);
    auto expression_angle = LMCAS::vector_angle_checked(
        {two_as_expr, zero}, {SymbolicExpr::number(2), zero});
    ASSERT_TRUE((expression_angle.has_value())) << "checked vector_angle accepts finite numeric expressions";
    if (expression_angle) {
        EXPECT_TRUE((std::abs(expression_angle.value()) < 1e-9)) << "checked vector_angle evaluates expression components";
    }
}

TEST(SymbolicVectorGeometryVector, DotCheckedCancellation) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::vector_dot_checked({one, zero}, {zero, one}, context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked vector_dot observes cancelled context";
}
