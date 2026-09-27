#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

TEST(SymbolicVectorGeometryLinePlane, IntersectionPoint) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic line;
    line.point = {zero, zero, zero};
    line.direction = {one, one, one};

    LMCAS::PlaneSymbolic plane;
    plane.normal = {one, one, one};
    plane.d = three;

    auto result = LMCAS::line_plane_intersection_checked(line, plane).value();
    EXPECT_TRUE((result.size() == 3)) << "intersection has 3 coordinates";
    for (int i = 0; i < 3; ++i) {
        auto s = result[i] ? result[i]->simplify() : nullptr;
        EXPECT_TRUE(test_expression_text((s), ("1"))) << "intersection coord " + std::to_string(i) + " is 1";
    }
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionTinyScale) {
    auto zero = SymbolicExpr::number(0);
    auto tiny = SymbolicExpr::number(1.0e-200);
    LMCAS::LineSymbolic line{
        {zero, zero, zero}, {tiny, zero, zero}};
    LMCAS::PlaneSymbolic plane{{tiny, zero, zero}, tiny};

    auto result = LMCAS::line_plane_intersection_checked(line, plane);
    ASSERT_TRUE((result.has_value())) << "nonparallel tiny directions retain a unique intersection";
    if (result) {
        auto x = test_numeric_eval(result.value()[0]);
        EXPECT_TRUE((x.has_value() && std::abs(*x - 1.0) < 1e-12)) << "line and plane scales cancel from the intersection";
    }
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionCancellingDot) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto maximum = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto negative_maximum = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    LMCAS::LineSymbolic line{
        {maximum, negative_maximum, zero}, {one, zero, zero}};
    LMCAS::PlaneSymbolic plane{{two, two, zero}, zero};

    auto result = LMCAS::line_plane_intersection_checked(line, plane);
    ASSERT_TRUE((result.has_value())) << "a finite on-plane point survives cancelling dot products";
    if (result) {
        auto x = test_numeric_eval(result.value()[0]);
        auto y = test_numeric_eval(result.value()[1]);
        EXPECT_TRUE((x.has_value() && y.has_value() &&
                     *x == std::numeric_limits<double>::max() &&
                     *y == -std::numeric_limits<double>::max()))
            << "scaled intersection preserves the extreme on-plane point";
    }
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionCancellingUpdate) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto maximum = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto negative_maximum = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    LMCAS::LineSymbolic line{
        {negative_maximum, zero, zero}, {two, zero, zero}};
    LMCAS::PlaneSymbolic plane{{one, zero, zero}, maximum};

    auto result = LMCAS::line_plane_intersection_checked(line, plane);
    ASSERT_TRUE((result.has_value())) << "finite intersection survives an overflowing displacement";
    if (result) {
        auto x = test_numeric_eval(result.value()[0]);
        EXPECT_TRUE((x.has_value() &&
                     *x == std::numeric_limits<double>::max()))
            << "scaled coordinate update preserves the finite endpoint";
    }
}

TEST(SymbolicVectorGeometryLinePlane, PointPlaneKnownDistance) {
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto zero = SymbolicExpr::number(0);

    std::vector<std::shared_ptr<SymbolicExpr>> point = {one, two, three};

    LMCAS::PlaneSymbolic plane;
    plane.normal = {one, one, one};
    plane.d = zero;

    auto result = LMCAS::point_plane_distance_checked(point, plane).value();
    auto val = test_numeric_eval(result ? result->simplify() : nullptr);
    double expected = 6.0 / std::sqrt(3.0);
    ASSERT_TRUE((val.has_value())) << "point-plane distance is numeric";
    if (val) {
        EXPECT_TRUE((std::abs(*val - expected) < 1e-6)) << "point-plane distance is 6/sqrt(3)";
    }
}

TEST(SymbolicVectorGeometryLinePlane, PointPlaneHugeNormal) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto huge = SymbolicExpr::number(1e200);
    LMCAS::PlaneSymbolic plane{{huge, zero, zero}, zero};

    auto result = LMCAS::point_plane_distance_checked(
        {one, zero, zero}, plane);
    ASSERT_TRUE((result.has_value())) << "finite huge plane normal has a representable distance";
    if (result) {
        auto value = test_numeric_eval(result.value());
        EXPECT_TRUE((value.has_value() && std::abs(*value - 1.0) < 1e-12)) << "huge normal scale cancels from point-plane distance";
    }
}

TEST(SymbolicVectorGeometryLinePlane, PointPlaneCancellingDot) {
    auto zero = SymbolicExpr::number(0);
    auto two = SymbolicExpr::number(2);
    auto maximum = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto negative_maximum = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    LMCAS::PlaneSymbolic plane{{two, two, zero}, zero};

    auto result = LMCAS::point_plane_distance_checked(
        {maximum, negative_maximum, zero}, plane);
    ASSERT_TRUE((result.has_value())) << "cancelling extreme dot products define a finite distance";
    if (result) {
        auto value = test_numeric_eval(result.value());
        EXPECT_TRUE((value.has_value() && *value == 0.0)) << "scaled point-normal dot product preserves cancellation";
    }
}

TEST(SymbolicVectorGeometryLinePlane, SkewKnownDistance) {
    /**
     * @brief x 轴与过 (0,1,0)、方向为 (0,0,1) 的异面直线间距为 1。
     * 第一条直线过 (0,0,0)，方向 d1=(1,0,0)；第二条平行 z 轴，方向 d2=(0,0,1)。
     * d1×d2=(0,-1,0)，|d1×d2|=1，a2-a1=(0,1,0)，
     * (a2-a1)·(d1×d2)=-1，距离 |(a2-a1)·(d1×d2)|/|d1×d2|=1。
     */
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    LMCAS::LineSymbolic l1;
    l1.point = {zero, zero, zero};
    l1.direction = {one, zero, zero};

    LMCAS::LineSymbolic l2;
    l2.point = {zero, one, zero};
    l2.direction = {zero, zero, one};

    auto result = LMCAS::skew_lines_distance_checked(l1, l2).value();
    auto simplified = result ? result->simplify() : nullptr;
    auto val = test_numeric_eval(simplified);
    if (val) {
        EXPECT_TRUE((std::abs(*val - 1.0) < 1e-9)) << "skew lines distance is 1";
    } else {
        EXPECT_TRUE(test_expression_text((simplified), ("1"))) << "skew lines distance simplifies to 1";
    }
}

TEST(SymbolicVectorGeometryLinePlane, SkewHugeDirections) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto huge = SymbolicExpr::number(1e200);
    LMCAS::LineSymbolic l1{
        {zero, zero, zero}, {huge, zero, zero}};
    LMCAS::LineSymbolic l2{
        {zero, one, zero}, {zero, zero, huge}};

    auto result = LMCAS::skew_lines_distance_checked(l1, l2);
    ASSERT_TRUE((result.has_value())) << "finite huge directions define a representable skew distance";
    if (result) {
        auto value = test_numeric_eval(result.value());
        EXPECT_TRUE((value.has_value() && std::abs(*value - 1.0) < 1e-12)) << "direction scales cancel from skew-line distance";
    }
}

TEST(SymbolicVectorGeometryLinePlane, SkewExtremeOffsets) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto maximum = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto negative_maximum = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    constexpr double epsilon = 1e-200;
    auto small = SymbolicExpr::number(epsilon);
    LMCAS::LineSymbolic l1{
        {zero, negative_maximum, zero}, {one, zero, zero}};
    LMCAS::LineSymbolic l2{
        {zero, maximum, zero}, {zero, one, small}};

    auto result = LMCAS::skew_lines_distance_checked(l1, l2);
    ASSERT_TRUE((result.has_value())) << "extreme finite point offsets have a representable skew distance";
    if (result) {
        auto value = test_numeric_eval(result.value());
        const double expected =
            (std::numeric_limits<double>::max() * epsilon) * 2.0;
        EXPECT_TRUE((value.has_value() && std::isfinite(*value) &&
                     std::abs(*value / expected - 1.0) < 1e-12))
            << "point-offset scaling preserves skew-line distance";
    }
}

TEST(SymbolicVectorGeometryLinePlane, DihedralPerpendicular) {
    LMCAS::PlaneSymbolic p1{{SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)}, SymbolicExpr::number(0)};
    LMCAS::PlaneSymbolic p2{{SymbolicExpr::number(0), SymbolicExpr::number(1), SymbolicExpr::number(0)}, SymbolicExpr::number(0)};
    auto ang = LMCAS::dihedral_angle_checked(p1, p2).value();
    auto v = test_numeric_eval(ang ? ang->simplify() : nullptr);
    EXPECT_TRUE((v.has_value() && std::abs(*v - LMMC_CONST_PI / 2.0) < 1e-6)) << "dihedral angle of perpendicular planes is pi/2";
}

TEST(SymbolicVectorGeometryLinePlane, DihedralHugeNormals) {
    auto huge = SymbolicExpr::number(1e200);
    LMCAS::PlaneSymbolic p1{{huge, SymbolicExpr::number(0), SymbolicExpr::number(0)}, SymbolicExpr::number(0)};
    LMCAS::PlaneSymbolic p2{{SymbolicExpr::number(0), huge, SymbolicExpr::number(0)}, SymbolicExpr::number(0)};
    auto angle = LMCAS::dihedral_angle_checked(p1, p2);
    ASSERT_TRUE((angle.has_value())) << "finite huge plane normals remain valid";
    if (angle) {
        auto value = test_numeric_eval(angle.value());
        EXPECT_TRUE((value.has_value() &&
                     std::abs(*value - LMMC_CONST_PI / 2.0) < 1e-12))
            << "huge perpendicular normals retain pi/2 angle";
    }
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic line{{zero, zero, zero}, {one, one, one}};
    LMCAS::PlaneSymbolic plane{{one, one, one}, three};

    auto intersection = LMCAS::line_plane_intersection_checked(line, plane);
    ASSERT_TRUE((intersection.has_value())) << "checked line_plane_intersection succeeds";
    if (intersection) {
        EXPECT_TRUE((intersection.value().size() == 3)) << "checked intersection has three coordinates";
        EXPECT_TRUE(test_expression_text((intersection.value()[0]->simplify()), ("1"))) << "checked intersection x = 1";
    }
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionParallel) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic parallel{{zero, zero, zero}, {one, zero, zero}};
    LMCAS::PlaneSymbolic z_plane{{zero, zero, one}, one};
    auto no_unique = LMCAS::line_plane_intersection_checked(parallel, z_plane);
    EXPECT_TRUE((!no_unique &&
                 no_unique.error().code == LMCAS::CasErrc::DomainError))
        << "checked line_plane_intersection rejects parallel line-plane";
}

TEST(SymbolicVectorGeometryLinePlane, IntersectionZeroDirection) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::PlaneSymbolic plane{{one, one, one}, three};

    LMCAS::LineSymbolic zero_direction{{zero, zero, zero}, {zero, zero, zero}};
    auto bad_line = LMCAS::line_plane_intersection_checked(zero_direction, plane);
    EXPECT_TRUE((!bad_line &&
                 bad_line.error().code == LMCAS::CasErrc::DomainError))
        << "checked line_plane_intersection rejects zero direction";
}

TEST(SymbolicVectorGeometryLinePlane, PointPlaneZeroNormal) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::PlaneSymbolic zero_plane{{zero, zero, zero}, one};
    auto bad_distance = LMCAS::point_plane_distance_checked({one, two, three}, zero_plane);
    EXPECT_TRUE((!bad_distance &&
                 bad_distance.error().code == LMCAS::CasErrc::DomainError))
        << "checked point_plane_distance rejects zero plane normal";
}

TEST(SymbolicVectorGeometryLinePlane, SkewCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic l1{{zero, zero, zero}, {one, zero, zero}};
    LMCAS::LineSymbolic l2{{zero, one, zero}, {zero, zero, one}};
    auto skew = LMCAS::skew_lines_distance_checked(l1, l2);
    ASSERT_TRUE((skew.has_value())) << "checked skew_lines_distance succeeds";
    if (skew) {
        auto value = test_numeric_eval(skew.value()->simplify());
        EXPECT_TRUE((value.has_value() && std::abs(*value - 1.0) < 1e-9)) << "checked skew distance is 1";
    }
}

TEST(SymbolicVectorGeometryLinePlane, SkewParallel) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic l1{{zero, zero, zero}, {one, zero, zero}};

    LMCAS::LineSymbolic parallel_l2{{zero, one, zero}, {one, zero, zero}};
    auto parallel_distance = LMCAS::skew_lines_distance_checked(l1, parallel_l2);
    EXPECT_TRUE((!parallel_distance &&
                 parallel_distance.error().code == LMCAS::CasErrc::DomainError))
        << "checked skew_lines_distance rejects parallel directions";
}

TEST(SymbolicVectorGeometryLinePlane, DihedralCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto angle = LMCAS::dihedral_angle_checked(
        LMCAS::PlaneSymbolic{{one, zero, zero}, zero},
        LMCAS::PlaneSymbolic{{zero, one, zero}, zero});
    ASSERT_TRUE((angle.has_value())) << "checked dihedral_angle succeeds";
    if (angle) {
        auto value = test_numeric_eval(angle.value()->simplify());
        EXPECT_TRUE((value.has_value() && std::abs(*value - LMMC_CONST_PI / 2.0) < 1e-6)) << "checked dihedral angle is pi/2";
    }
}

TEST(SymbolicVectorGeometryLinePlane, PointPlaneCancellation) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::PlaneSymbolic plane{{one, one, one}, three};

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::point_plane_distance_checked({one, two, three}, plane, context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked point_plane_distance observes cancellation";
}
