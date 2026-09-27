#include "expr.hpp"
#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

TEST(SymbolicVectorGeometrySurface, SurfaceSphereNormalAndTangent) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::add(SymbolicExpr::multiply(z, z), SymbolicExpr::number(-1)));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> pt = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};

    auto checked_normal = LMCAS::surface_normal_checked(surf, pt);
    ASSERT_TRUE((checked_normal.has_value())) << "checked surface normal succeeds";
    if (checked_normal) {
        EXPECT_TRUE(test_proved_equivalent(checked_normal.value()[0], parse_expr("1").value())) << "unit normal x = 1";
        EXPECT_TRUE(test_proved_equivalent(checked_normal.value()[1], parse_expr("0").value())) << "unit normal y = 0";
        EXPECT_TRUE(test_proved_equivalent(checked_normal.value()[2], parse_expr("0").value())) << "unit normal z = 0";
    }

    auto checked_plane = LMCAS::tangent_plane_checked(surf, pt);
    ASSERT_TRUE((checked_plane.has_value())) << "checked tangent plane succeeds";
    if (checked_plane) {
        EXPECT_TRUE(test_proved_equivalent(checked_plane.value().normal[0], parse_expr("2").value())) << "tangent plane normal x = 2";
        EXPECT_TRUE(test_proved_equivalent(checked_plane.value().d, parse_expr("2").value())) << "tangent plane d = 2";
    }
}

TEST(SymbolicVectorGeometrySurface, SurfaceHugeGradient) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto huge = SymbolicExpr::number(1e200);
    auto F = SymbolicExpr::add(
        SymbolicExpr::multiply(huge, x),
        SymbolicExpr::multiply(huge, y));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> origin = {SymbolicExpr::number(0), SymbolicExpr::number(0), SymbolicExpr::number(0)};

    auto normal = LMCAS::surface_normal_checked(surf, origin);
    ASSERT_TRUE((normal.has_value())) << "finite huge gradient has a representable unit normal";
    if (normal) {
        auto nx = test_numeric_eval(normal.value()[0]);
        auto ny = test_numeric_eval(normal.value()[1]);
        EXPECT_TRUE((nx.has_value() && ny.has_value() &&
                     std::abs(*nx - std::sqrt(0.5)) < 1e-12 &&
                     std::abs(*ny - std::sqrt(0.5)) < 1e-12))
            << "huge gradient normal is scaled before normalization";
    }

    auto plane = LMCAS::tangent_plane_checked(surf, origin);
    ASSERT_TRUE((plane.has_value())) << "finite huge gradient defines a tangent plane";
    if (plane) {
        auto nx = test_numeric_eval(plane.value().normal[0]);
        auto ny = test_numeric_eval(plane.value().normal[1]);
        EXPECT_TRUE((nx.has_value() && ny.has_value() &&
                     *nx == 1e200 && *ny == 1e200))
            << "tangent plane preserves finite huge coefficients";
    }
}

TEST(SymbolicVectorGeometrySurface, TangentPlaneCancellingDot) {
    const double maximum = std::numeric_limits<double>::max();
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto F = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y));
    SurfaceSymbolic surf{F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> point = {
        SymbolicExpr::number(maximum),
        SymbolicExpr::number(-maximum),
        SymbolicExpr::number(0)};

    auto plane = LMCAS::tangent_plane_checked(surf, point);
    ASSERT_TRUE((plane.has_value())) << "finite tangent plane survives a cancelling gradient dot";
    if (plane) {
        auto d = test_numeric_eval(plane.value().d);
        EXPECT_TRUE((d.has_value() && *d == 0.0)) << "cancelling extreme tangent-plane constant is zero";
    }
}

TEST(SymbolicVectorGeometrySurface, SurfaceSingularNormal) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto zero = SymbolicExpr::number(0);

    auto singular_F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(z, z));
    SurfaceSymbolic singular{singular_F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> origin = {zero, zero, zero};

    auto singular_normal = LMCAS::surface_normal_checked(singular, origin);
    EXPECT_TRUE((!singular_normal.has_value())) << "checked surface normal rejects singular point";
    EXPECT_TRUE((singular_normal.error().code == LMCAS::CasErrc::DomainError)) << "checked surface normal reports DomainError at singular point";
}

TEST(SymbolicVectorGeometrySurface, SurfaceSingularTangent) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto zero = SymbolicExpr::number(0);

    auto singular_F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(z, z));
    SurfaceSymbolic singular{singular_F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> origin = {zero, zero, zero};

    auto singular_plane = LMCAS::tangent_plane_checked(singular, origin);
    EXPECT_TRUE((!singular_plane.has_value())) << "checked tangent plane rejects singular point";
    EXPECT_TRUE((singular_plane.error().code == LMCAS::CasErrc::DomainError)) << "checked tangent plane reports DomainError at singular point";
}

TEST(SymbolicVectorGeometrySurface, SurfaceUnsupportedDerivative) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto zero = SymbolicExpr::number(0);

    auto singular_F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(z, z));
    SurfaceSymbolic singular{singular_F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> origin = {zero, zero, zero};

    auto unsupported_F = SymbolicExpr::eq(x, zero);
    SurfaceSymbolic unsupported{unsupported_F, {"x", "y", "z"}};
    auto unsupported_normal = LMCAS::surface_normal_checked(unsupported, origin);
    EXPECT_TRUE((!unsupported_normal.has_value())) << "checked surface normal rejects unsupported derivatives";
    EXPECT_TRUE((unsupported_normal.error().code == LMCAS::CasErrc::Inconclusive)) << "checked surface normal reports Inconclusive for unsupported derivatives";
}

TEST(SymbolicVectorGeometrySurface, SurfaceTangentCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto zero = SymbolicExpr::number(0);

    auto singular_F = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
        SymbolicExpr::multiply(z, z));
    SurfaceSymbolic singular{singular_F, {"x", "y", "z"}};
    std::vector<std::shared_ptr<SymbolicExpr>> origin = {zero, zero, zero};

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::tangent_plane_checked(singular, origin, context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked tangent plane observes cancellation";
    EXPECT_TRUE((cancelled.error().code == LMCAS::CasErrc::Cancelled)) << "checked tangent plane reports Cancelled";
}
