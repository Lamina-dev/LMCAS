#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

TEST(SymbolicVectorGeometryQuadric, QuadricSphere) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::add(SymbolicExpr::multiply(z, z), SymbolicExpr::number(-1)));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};
    std::string c = LMCAS::classify_quadric_checked(surf).value();
    EXPECT_TRUE((c == "sphere")) << "x^2+y^2+z^2=1 classified as sphere";

    auto checked = LMCAS::classify_quadric_checked(surf);
    ASSERT_TRUE((checked.has_value())) << "checked classify_quadric succeeds";
    if (checked) {
        EXPECT_TRUE((checked.value() == "sphere")) << "checked classify_quadric reports sphere";
    }
}

TEST(SymbolicVectorGeometryQuadric, QuadricEquationScaling) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto quadratic = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                          SymbolicExpr::multiply(y, y)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(z, z)),
            SymbolicExpr::number(-1)));
    auto F = SymbolicExpr::multiply(
        SymbolicExpr::number(1e-12), quadratic);
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((classification.has_value() &&
                 classification.value() == "hyperboloid"))
        << "nonzero equation scaling preserves hyperboloid classification";
}

TEST(SymbolicVectorGeometryQuadric, QuadricUnsquaredAxis) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                          SymbolicExpr::multiply(z, z)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((classification.has_value() &&
                 classification.value() == "paraboloid"))
        << "x^2+z^2-y=0 is a paraboloid";
}

TEST(SymbolicVectorGeometryQuadric, QuadricLargeRangeTerm) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                          SymbolicExpr::multiply(y, y)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::number(DBL_MAX), x),
            z));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((classification && classification.value() == "paraboloid")) << "a large range-space term does not hide the null-axis term";
}

TEST(SymbolicVectorGeometryQuadric, QuadricTranslatedCone) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto shifted_x_square = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::multiply(x, x),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(shifted_x_square,
                          SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(z, z)));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((classification.has_value() &&
                 classification.value() == "cone"))
        << "(x-1)^2+y^2-z^2=0 is a translated cone";
}

TEST(SymbolicVectorGeometryQuadric, QuadricRotatedAxes) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto xy = SymbolicExpr::multiply(x, y);
    auto rotated_ellipsoid = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(x, x)),
                SymbolicExpr::multiply(SymbolicExpr::number(2), xy)),
            SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(y, y))),
        SymbolicExpr::add(SymbolicExpr::multiply(z, z), SymbolicExpr::number(-1)));
    SurfaceSymbolic ellipsoid{
        rotated_ellipsoid, {"x", "y", "z"}};

    auto ellipsoid_classification =
        LMCAS::classify_quadric_checked(ellipsoid);
    EXPECT_TRUE((ellipsoid_classification &&
                 ellipsoid_classification.value() == "ellipsoid"))
        << "positive-definite mixed quadric is classified after rotation";

    auto rotated_paraboloid = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::multiply(x, x),
                SymbolicExpr::multiply(SymbolicExpr::number(2), xy)),
            SymbolicExpr::multiply(y, y)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(z, z),
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(-1), x), y)));
    SurfaceSymbolic paraboloid{
        rotated_paraboloid, {"x", "y", "z"}};

    auto paraboloid_classification =
        LMCAS::classify_quadric_checked(paraboloid);
    EXPECT_TRUE((paraboloid_classification &&
                 paraboloid_classification.value() == "paraboloid"))
        << "mixed rank-two quadric follows its rotated null axis";
}

TEST(SymbolicVectorGeometryQuadric, QuadricParabolicCylinder) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto F = SymbolicExpr::add(
        SymbolicExpr::multiply(x, x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((classification && classification.value() == "cylinder")) << "x^2-y=0 is a parabolic cylinder along the z-axis";
}

TEST(SymbolicVectorGeometryQuadric, QuadricUnresolvedEigenvalue) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::multiply(x, x),
            SymbolicExpr::multiply(
                SymbolicExpr::number(1e-16),
                SymbolicExpr::multiply(y, y))),
        SymbolicExpr::add(SymbolicExpr::multiply(z, z), SymbolicExpr::number(-1)));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};

    auto classification = LMCAS::classify_quadric_checked(surf);
    EXPECT_TRUE((!classification &&
                 classification.error().code == LMCAS::CasErrc::Inconclusive))
        << "an eigenvalue inside the backward-error bound is unresolved";
}

TEST(SymbolicVectorGeometryQuadric, QuadricEmptySphere) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto squares = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                          SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(z, z));

    SurfaceSymbolic empty_sphere{
        SymbolicExpr::add(squares, SymbolicExpr::number(1)), {"x", "y", "z"}};
    auto empty = LMCAS::classify_quadric_checked(empty_sphere);
    EXPECT_TRUE((!empty &&
                 empty.error().code == LMCAS::CasErrc::Inconclusive))
        << "x^2+y^2+z^2+1=0 is not reported as a real sphere";
}

TEST(SymbolicVectorGeometryQuadric, QuadricPointLocus) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto point_equation = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
            SymbolicExpr::number(1)),
        SymbolicExpr::add(SymbolicExpr::multiply(y, y),
                          SymbolicExpr::multiply(z, z)));
    SurfaceSymbolic point_surface{
        point_equation, {"x", "y", "z"}};
    auto point = LMCAS::classify_quadric_checked(point_surface);
    EXPECT_TRUE((!point &&
                 point.error().code == LMCAS::CasErrc::Inconclusive))
        << "(x-1)^2+y^2+z^2=0 is a point, not a sphere";
}

TEST(SymbolicVectorGeometryQuadric, QuadricEmptyCylinder) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    SurfaceSymbolic empty_cylinder{
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(y, y)),
            SymbolicExpr::number(1)),
        {"x", "y", "z"}};
    auto empty_rank_two =
        LMCAS::classify_quadric_checked(empty_cylinder);
    EXPECT_TRUE((!empty_rank_two &&
                 empty_rank_two.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "x^2+y^2+1=0 is not reported as a real cylinder";
}

TEST(SymbolicVectorGeometryQuadric, QuadricRealCylinder) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    SurfaceSymbolic real_cylinder{
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(y, y)),
            SymbolicExpr::number(-1)),
        {"x", "y", "z"}};
    auto real_rank_two =
        LMCAS::classify_quadric_checked(real_cylinder);
    EXPECT_TRUE((real_rank_two &&
                 real_rank_two.value() == "cylinder"))
        << "x^2+y^2-1=0 remains a real cylinder";
}

TEST(SymbolicVectorGeometryQuadric, QuadricPlanePair) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    SurfaceSymbolic plane_pair{
        SymbolicExpr::add(
            SymbolicExpr::multiply(x, x),
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1), SymbolicExpr::multiply(y, y))),
        {"x", "y", "z"}};
    auto degenerate_rank_two =
        LMCAS::classify_quadric_checked(plane_pair);
    EXPECT_TRUE((!degenerate_rank_two &&
                 degenerate_rank_two.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "x^2-y^2=0 is a plane pair, not a cylinder";
}

TEST(SymbolicVectorGeometryQuadric, QuadricNullEquation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    SurfaceSymbolic invalid{nullptr, {"x", "y", "z"}};
    auto invalid_result = LMCAS::classify_quadric_checked(invalid);
    EXPECT_TRUE((!invalid_result &&
                 invalid_result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked classify_quadric rejects null surface equation";
}

TEST(SymbolicVectorGeometryQuadric, QuadricLinearEquation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto linear = SymbolicExpr::add(x, y);
    SurfaceSymbolic non_quadric{linear, {"x", "y", "z"}};
    auto unknown = LMCAS::classify_quadric_checked(non_quadric);
    EXPECT_TRUE((!unknown &&
                 unknown.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked classify_quadric reports an unclassified surface as Inconclusive";
}

TEST(SymbolicVectorGeometryQuadric, QuadricSymbolicCoefficient) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto a = SymbolicExpr::variable("a");

    auto symbolic_coeff = SymbolicExpr::add(
        SymbolicExpr::multiply(a, SymbolicExpr::multiply(x, x)),
        SymbolicExpr::add(SymbolicExpr::multiply(y, y),
                          SymbolicExpr::multiply(z, z)));
    SurfaceSymbolic symbolic{symbolic_coeff, {"x", "y", "z"}};
    auto unsupported = LMCAS::classify_quadric_checked(symbolic);
    EXPECT_TRUE((!unsupported &&
                 unsupported.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked classify_quadric rejects unproved symbolic coefficients";
}

TEST(SymbolicVectorGeometryQuadric, QuadricCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto linear = SymbolicExpr::add(x, y);
    SurfaceSymbolic non_quadric{linear, {"x", "y", "z"}};

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::classify_quadric_checked(non_quadric, context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked classify_quadric observes cancellation";
}
