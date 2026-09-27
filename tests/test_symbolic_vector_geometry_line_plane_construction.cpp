#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

TEST(SymbolicVectorGeometryLinePlaneConstruction, LineFromPointsDirection) {
    std::vector<std::shared_ptr<SymbolicExpr>> p1 = {SymbolicExpr::number(0), SymbolicExpr::number(0), SymbolicExpr::number(0)};
    std::vector<std::shared_ptr<SymbolicExpr>> p2 = {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)};
    auto line = LMCAS::line_from_two_points_checked(p1, p2).value();
    EXPECT_TRUE(test_expression_text((line.direction[0]->simplify()), ("1"))) << "dir.x = 1";
    EXPECT_TRUE(test_expression_text((line.direction[1]->simplify()), ("2"))) << "dir.y = 2";
    EXPECT_TRUE(test_expression_text((line.direction[2]->simplify()), ("3"))) << "dir.z = 3";
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, PlaneFromPointsNormal) {
    std::vector<std::shared_ptr<SymbolicExpr>> p1 = {SymbolicExpr::number(0), SymbolicExpr::number(0), SymbolicExpr::number(0)};
    std::vector<std::shared_ptr<SymbolicExpr>> p2 = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};
    std::vector<std::shared_ptr<SymbolicExpr>> p3 = {SymbolicExpr::number(0), SymbolicExpr::number(1), SymbolicExpr::number(0)};
    auto plane = LMCAS::plane_from_three_points_checked(p1, p2, p3).value();
    EXPECT_TRUE(test_expression_text((plane.normal[0]->simplify()), ("0"))) << "n.x = 0";
    EXPECT_TRUE(test_expression_text((plane.normal[1]->simplify()), ("0"))) << "n.y = 0";
    EXPECT_TRUE(test_expression_text((plane.normal[2]->simplify()), ("1"))) << "n.z = 1";
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, LineConstructionCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto constructed_line = LMCAS::line_from_two_points_checked(
        {zero, zero, zero}, {one, two, three});
    ASSERT_TRUE((constructed_line.has_value())) << "checked line_from_two_points succeeds";
    if (constructed_line) {
        EXPECT_TRUE(test_expression_text((constructed_line.value().direction[2]->simplify()), ("3"))) << "checked constructed line z direction = 3";
    }
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, LineIdenticalPoints) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto identical_points = LMCAS::line_from_two_points_checked(
        {one, one, one}, {one, one, one});
    EXPECT_TRUE((!identical_points &&
                 identical_points.error().code == LMCAS::CasErrc::DomainError))
        << "checked line_from_two_points rejects identical points";
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, PlaneConstructionCheckedValue) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto checked_plane = LMCAS::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {zero, one, zero});
    ASSERT_TRUE((checked_plane.has_value())) << "checked plane_from_three_points succeeds";
    if (checked_plane) {
        EXPECT_TRUE(test_expression_text((checked_plane.value().normal[2]->simplify()), ("1"))) << "checked plane normal z = 1";
    }
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, PlaneCollinearPoints) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto collinear_plane = LMCAS::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {two, zero, zero});
    EXPECT_TRUE((!collinear_plane &&
                 collinear_plane.error().code == LMCAS::CasErrc::DomainError))
        << "checked plane_from_three_points rejects collinear points";
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, LineExtremePoints) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto negative_max = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    auto positive_max = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto extreme_line = LMCAS::line_from_two_points_checked(
        {negative_max, zero, zero}, {positive_max, zero, zero});
    ASSERT_TRUE((extreme_line.has_value())) << "opposite extreme finite points define a line";
    if (extreme_line) {
        auto dx = test_numeric_eval(extreme_line.value().direction[0]);
        EXPECT_TRUE((dx.has_value() && std::isfinite(*dx) && *dx > 0.0)) << "line direction is scaled before point subtraction overflows";
    }
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, PlaneHugeEdges) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto huge = SymbolicExpr::number(1.0e200);
    auto huge_plane = LMCAS::plane_from_three_points_checked(
        {zero, zero, zero}, {huge, zero, zero}, {zero, huge, zero});
    ASSERT_TRUE((huge_plane.has_value())) << "finite huge edges define a representable plane direction";
    if (huge_plane) {
        auto nx = test_numeric_eval(huge_plane.value().normal[0]);
        auto ny = test_numeric_eval(huge_plane.value().normal[1]);
        auto nz = test_numeric_eval(huge_plane.value().normal[2]);
        EXPECT_TRUE((nx.has_value() && ny.has_value() && nz.has_value() &&
                     *nx == 0.0 && *ny == 0.0 && *nz > 0.0 &&
                     std::isfinite(*nz)))
            << "plane construction scales edge vectors before crossing";
    }
}

TEST(SymbolicVectorGeometryLinePlaneConstruction, PlaneExtremePoints) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    auto negative_max = SymbolicExpr::number(
        -std::numeric_limits<double>::max());
    auto positive_max = SymbolicExpr::number(
        std::numeric_limits<double>::max());
    auto extreme_plane = LMCAS::plane_from_three_points_checked(
        {negative_max, zero, zero},
        {positive_max, zero, zero},
        {negative_max, positive_max, zero});
    ASSERT_TRUE((extreme_plane.has_value())) << "extreme finite point differences define a plane";
    if (extreme_plane) {
        auto nz = test_numeric_eval(extreme_plane.value().normal[2]);
        EXPECT_TRUE((nz.has_value() && std::isfinite(*nz) && *nz > 0.0)) << "plane edges are formed after common point scaling";
    }
}
