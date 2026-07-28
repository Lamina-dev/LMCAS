#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cmath>

void test_vector_dot() {
    TEST_CASE("Vector Dot: Orthogonal vectors (expect 0)");
    {
        // i . j = 0
        auto one = SymbolicExpr::number(1);
        auto zero = SymbolicExpr::number(0);
        std::vector<std::shared_ptr<SymbolicExpr>> a = {one, zero, zero};
        std::vector<std::shared_ptr<SymbolicExpr>> b = {zero, one, zero};
        auto result = lamina::vector_dot(a, b);
        std::string s = result ? result->to_string() : "null";
        std::cout << "  i . j = " << s << std::endl;
        // Simplify to check for zero
        auto simplified = result ? result->simplify() : nullptr;
        std::string ss = simplified ? simplified->to_string() : "null";
        std::cout << "  simplified: " << ss << std::endl;
        EXPECT_EQ_EXPR_STR(simplified, "0", "orthogonal dot product is 0");
    }

    TEST_CASE("Vector Dot: Parallel vectors");
    {
        // (1,2,3) . (2,4,6) = 2+8+18 = 28
        auto result = lamina::vector_dot(
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}
        );
        auto simplified = result ? result->simplify() : nullptr;
        std::string s = simplified ? simplified->to_string() : "null";
        std::cout << "  (1,2,3).(2,4,6) = " << s << std::endl;
        EXPECT_EQ_EXPR_STR(simplified, "28", "parallel dot product is 28");
    }

    TEST_CASE("Vector Dot: General vectors");
    {
        // (1,0,2) . (3,4,1) = 3+0+2 = 5
        auto result = lamina::vector_dot(
            {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4), SymbolicExpr::number(1)}
        );
        auto simplified = result ? result->simplify() : nullptr;
        std::string s = simplified ? simplified->to_string() : "null";
        std::cout << "  (1,0,2).(3,4,1) = " << s << std::endl;
        EXPECT_EQ_EXPR_STR(simplified, "5", "general dot product is 5");
    }
}

void test_vector_cross() {
    TEST_CASE("Vector Cross: Parallel vectors (zero vector)");
    {
        // (1,2,3) x (2,4,6) = (0,0,0)
        auto result = lamina::vector_cross(
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}
        );
        EXPECT_TRUE(result.size() == 3, "cross product has 3 components");
        for (int i = 0; i < 3; ++i) {
            auto s = result[i] ? result[i]->simplify() : nullptr;
            std::string str = s ? s->to_string() : "null";
            std::cout << "  parallel cross[" << i << "] = " << str << std::endl;
            EXPECT_EQ_EXPR_STR(s, "0", "parallel cross component " + std::to_string(i) + " is 0");
        }
    }

    TEST_CASE("Vector Cross: i x j = k");
    {
        // (1,0,0) x (0,1,0) = (0,0,1)
        auto one = SymbolicExpr::number(1);
        auto zero = SymbolicExpr::number(0);
        auto result = lamina::vector_cross(
            {one, zero, zero},
            {zero, one, zero}
        );
        EXPECT_TRUE(result.size() == 3, "i x j has 3 components");
        auto s0 = result[0] ? result[0]->simplify() : nullptr;
        auto s1 = result[1] ? result[1]->simplify() : nullptr;
        auto s2 = result[2] ? result[2]->simplify() : nullptr;
        std::cout << "  i x j = (" << (s0 ? s0->to_string() : "null") << ", "
                  << (s1 ? s1->to_string() : "null") << ", "
                  << (s2 ? s2->to_string() : "null") << ")" << std::endl;
        EXPECT_EQ_EXPR_STR(s0, "0", "i x j component x is 0");
        EXPECT_EQ_EXPR_STR(s1, "0", "i x j component y is 0");
        EXPECT_EQ_EXPR_STR(s2, "1", "i x j component z is 1");
    }

    TEST_CASE("Vector Cross: General 3D vectors");
    {
        // (1,2,3) x (4,5,6) = (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
        auto result = lamina::vector_cross(
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
            {SymbolicExpr::number(4), SymbolicExpr::number(5), SymbolicExpr::number(6)}
        );
        EXPECT_TRUE(result.size() == 3, "general cross has 3 components");
        auto s0 = result[0] ? result[0]->simplify() : nullptr;
        auto s1 = result[1] ? result[1]->simplify() : nullptr;
        auto s2 = result[2] ? result[2]->simplify() : nullptr;
        std::cout << "  (1,2,3) x (4,5,6) = (" << (s0 ? s0->to_string() : "null") << ", "
                  << (s1 ? s1->to_string() : "null") << ", "
                  << (s2 ? s2->to_string() : "null") << ")" << std::endl;
        EXPECT_EQ_EXPR_STR(s0, "-3", "general cross x = -3");
        EXPECT_EQ_EXPR_STR(s1, "6", "general cross y = 6");
        EXPECT_EQ_EXPR_STR(s2, "-3", "general cross z = -3");
    }
}

void test_vector_angle() {
    TEST_CASE("Vector Angle: Orthogonal vectors (pi/2)");
    {
        // angle between (1,0,0) and (0,1,0) = pi/2
        auto one = SymbolicExpr::number(1);
        auto zero = SymbolicExpr::number(0);
        double angle = lamina::vector_angle(
            {one, zero, zero},
            {zero, one, zero}
        );
        std::cout << "  angle(i, j) = " << angle << " (expected " << M_PI / 2.0 << ")" << std::endl;
        EXPECT_TRUE(std::abs(angle - M_PI / 2.0) < 1e-9, "orthogonal angle is pi/2");
    }

    TEST_CASE("Vector Angle: Parallel vectors (0)");
    {
        // angle between (1,2,3) and (2,4,6) = 0
        double angle = lamina::vector_angle(
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}
        );
        std::cout << "  angle((1,2,3), (2,4,6)) = " << angle << " (expected 0)" << std::endl;
        EXPECT_TRUE(std::abs(angle - 0.0) < 1e-9, "parallel angle is 0");
    }
}

void test_vector_checked_contracts() {
    TEST_CASE("Vector Geometry Checked API Contracts");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto dot = lamina::vector_dot_checked({one, zero}, {zero, one});
    EXPECT_TRUE(dot.has_value(), "checked vector_dot succeeds");
    if (dot) {
        EXPECT_EQ_EXPR_STR(dot.value()->simplify(), "0",
                           "checked dot product of orthogonal vectors is 0");
    }

    auto bad_dot = lamina::vector_dot_checked({one}, {one, zero});
    EXPECT_TRUE(!bad_dot &&
                    bad_dot.error().code == lamina::CasErrc::InvalidArgument,
                "checked vector_dot rejects dimension mismatch");
    bool legacy_dot_threw = false;
    try {
        (void)lamina::vector_dot({one}, {one, zero});
    } catch (const std::invalid_argument&) {
        legacy_dot_threw = true;
    }
    EXPECT_TRUE(legacy_dot_threw,
                "legacy vector_dot preserves invalid_argument for dimension mismatch");

    auto bad_component = lamina::vector_cross_checked({one, nullptr, zero},
                                                     {zero, one, zero});
    EXPECT_TRUE(!bad_component &&
                    bad_component.error().code == lamina::CasErrc::InvalidArgument,
                "checked vector_cross rejects null component");

    auto bad_cross_dim = lamina::vector_cross_checked({one, zero},
                                                     {zero, one});
    EXPECT_TRUE(!bad_cross_dim &&
                    bad_cross_dim.error().code == lamina::CasErrc::InvalidArgument,
                "checked vector_cross rejects non-3D vectors");

    auto angle = lamina::vector_angle_checked({one, zero}, {zero, one});
    EXPECT_TRUE(angle.has_value(), "checked vector_angle succeeds");
    if (angle) {
        EXPECT_TRUE(std::abs(angle.value() - M_PI / 2.0) < 1e-9,
                    "checked vector_angle of orthogonal vectors is pi/2");
    }

    auto zero_angle = lamina::vector_angle_checked({zero, zero}, {one, zero});
    EXPECT_TRUE(!zero_angle &&
                    zero_angle.error().code == lamina::CasErrc::DomainError,
                "checked vector_angle rejects zero-length vectors");
    EXPECT_TRUE(std::isnan(lamina::vector_angle({zero, zero}, {one, zero})),
                "legacy vector_angle unwraps domain failure to NaN");

    auto x = SymbolicExpr::variable("x");
    auto symbolic_angle = lamina::vector_angle_checked({x, zero}, {one, zero});
    EXPECT_TRUE(!symbolic_angle &&
                    symbolic_angle.error().code == lamina::CasErrc::NumericFailure,
                "checked vector_angle rejects symbolic components");

    auto two_as_expr = SymbolicExpr::add(one, one);
    auto expression_angle = lamina::vector_angle_checked(
        {two_as_expr, zero}, {SymbolicExpr::number(2), zero});
    EXPECT_TRUE(expression_angle.has_value(),
                "checked vector_angle accepts finite numeric expressions");
    if (expression_angle) {
        EXPECT_TRUE(std::abs(expression_angle.value()) < 1e-9,
                    "checked vector_angle evaluates expression components");
    }

    lamina::CancellationToken token;
    token.cancel();
    lamina::ComputationContext context({}, token);
    auto cancelled = lamina::vector_dot_checked({one, zero}, {zero, one}, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == lamina::CasErrc::Cancelled,
                "checked vector_dot observes cancelled context");
}

void test_line_plane_intersection() {
    TEST_CASE("Line-Plane Intersection: Known intersection point");
    {
        // Line: point (0,0,0), direction (1,1,1)
        // Plane: x + y + z = 3 (normal (1,1,1), d=3)
        // Intersection: t such that (1+1+1)*t = 3 => t=1 => point (1,1,1)
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);
        auto three = SymbolicExpr::number(3);

        lamina::LineSymbolic line;
        line.point = {zero, zero, zero};
        line.direction = {one, one, one};

        lamina::PlaneSymbolic plane;
        plane.normal = {one, one, one};
        plane.d = three;

        auto result = lamina::line_plane_intersection(line, plane);
        EXPECT_TRUE(result.size() == 3, "intersection has 3 coordinates");
        for (int i = 0; i < 3; ++i) {
            auto s = result[i] ? result[i]->simplify() : nullptr;
            std::string str = s ? s->to_string() : "null";
            std::cout << "  intersection[" << i << "] = " << str << std::endl;
            EXPECT_EQ_EXPR_STR(s, "1", "intersection coord " + std::to_string(i) + " is 1");
        }
    }
}

void test_point_plane_distance() {
    TEST_CASE("Point-Plane Distance: Known distance");
    {
        // Point (1,2,3), Plane: x + y + z = 0 (normal (1,1,1), d=0)
        // Distance = |1+2+3 - 0| / sqrt(1+1+1) = 6/sqrt(3) = 2*sqrt(3)
        auto one = SymbolicExpr::number(1);
        auto two = SymbolicExpr::number(2);
        auto three = SymbolicExpr::number(3);
        auto zero = SymbolicExpr::number(0);

        std::vector<std::shared_ptr<SymbolicExpr>> point = {one, two, three};

        lamina::PlaneSymbolic plane;
        plane.normal = {one, one, one};
        plane.d = zero;

        auto result = lamina::point_plane_distance(point, plane);
        std::string s = result ? result->to_string() : "null";
        std::cout << "  distance = " << s << std::endl;
        // 6/sqrt(3) = 2*sqrt(3) ≈ 3.4641
        // Use numeric evaluation to verify
        auto val = test_numeric_eval(result ? result->simplify() : nullptr);
        double expected = 6.0 / std::sqrt(3.0);
        std::cout << "  numeric = " << (val ? std::to_string(*val) : "non-numeric") << " (expected " << expected << ")" << std::endl;
        if (val) {
            EXPECT_TRUE(std::abs(*val - expected) < 1e-6, "point-plane distance is 6/sqrt(3)");
        } else {
            // If numeric eval fails, check structural content
            EXPECT_CONTAINS(s, {"6"}, "distance expression contains 6");
        }
    }
}

void test_skew_lines_distance() {
    TEST_CASE("Skew Lines Distance: Analytically known distance");
    {
        // Line 1: point (0,0,0), direction (1,0,0) — the x-axis
        // Line 2: point (0,1,0), direction (0,0,1) — parallel to z-axis, offset by 1 in y
        // These are skew lines. Distance = |(a2-a1) . (d1 x d2)| / |d1 x d2|
        // d1 x d2 = (1,0,0) x (0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        // |d1 x d2| = 1
        // (a2-a1) = (0,1,0)
        // (a2-a1) . (d1 x d2) = (0,1,0).(0,-1,0) = -1
        // Distance = |-1| / 1 = 1
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);

        lamina::LineSymbolic l1;
        l1.point = {zero, zero, zero};
        l1.direction = {one, zero, zero};

        lamina::LineSymbolic l2;
        l2.point = {zero, one, zero};
        l2.direction = {zero, zero, one};

        auto result = lamina::skew_lines_distance(l1, l2);
        std::string s = result ? result->to_string() : "null";
        std::cout << "  skew distance = " << s << std::endl;
        auto simplified = result ? result->simplify() : nullptr;
        std::string ss = simplified ? simplified->to_string() : "null";
        std::cout << "  simplified = " << ss << std::endl;
        // The distance should be 1
        auto val = test_numeric_eval(simplified);
        if (val) {
            std::cout << "  numeric = " << *val << std::endl;
            EXPECT_TRUE(std::abs(*val - 1.0) < 1e-9, "skew lines distance is 1");
        } else {
            // Structural check: result should simplify to 1
            EXPECT_EQ_EXPR_STR(simplified, "1", "skew lines distance simplifies to 1");
        }
    }
}

void test_line_plane_checked_contracts() {
    TEST_CASE("Line/Plane Vector Geometry Checked API Contracts");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    lamina::LineSymbolic line{{zero, zero, zero}, {one, one, one}};
    lamina::PlaneSymbolic plane{{one, one, one}, three};
    auto intersection = lamina::line_plane_intersection_checked(line, plane);
    EXPECT_TRUE(intersection.has_value(),
                "checked line_plane_intersection succeeds");
    if (intersection) {
        EXPECT_TRUE(intersection.value().size() == 3,
                    "checked intersection has three coordinates");
        EXPECT_EQ_EXPR_STR(intersection.value()[0]->simplify(), "1",
                           "checked intersection x = 1");
    }

    lamina::LineSymbolic parallel{{zero, zero, zero}, {one, zero, zero}};
    lamina::PlaneSymbolic z_plane{{zero, zero, one}, one};
    auto no_unique = lamina::line_plane_intersection_checked(parallel, z_plane);
    EXPECT_TRUE(!no_unique &&
                    no_unique.error().code == lamina::CasErrc::DomainError,
                "checked line_plane_intersection rejects parallel line-plane");
    EXPECT_TRUE(lamina::line_plane_intersection(parallel, z_plane).empty(),
                "legacy line_plane_intersection unwraps no-unique result to empty");

    lamina::LineSymbolic zero_direction{{zero, zero, zero}, {zero, zero, zero}};
    auto bad_line = lamina::line_plane_intersection_checked(zero_direction, plane);
    EXPECT_TRUE(!bad_line &&
                    bad_line.error().code == lamina::CasErrc::DomainError,
                "checked line_plane_intersection rejects zero direction");

    lamina::PlaneSymbolic zero_plane{{zero, zero, zero}, one};
    auto bad_distance = lamina::point_plane_distance_checked({one, two, three}, zero_plane);
    EXPECT_TRUE(!bad_distance &&
                    bad_distance.error().code == lamina::CasErrc::DomainError,
                "checked point_plane_distance rejects zero plane normal");
    EXPECT_TRUE(lamina::point_plane_distance({one, two, three}, zero_plane) == nullptr,
                "legacy point_plane_distance unwraps zero normal to nullptr");

    lamina::LineSymbolic l1{{zero, zero, zero}, {one, zero, zero}};
    lamina::LineSymbolic l2{{zero, one, zero}, {zero, zero, one}};
    auto skew = lamina::skew_lines_distance_checked(l1, l2);
    EXPECT_TRUE(skew.has_value(), "checked skew_lines_distance succeeds");
    if (skew) {
        auto value = test_numeric_eval(skew.value()->simplify());
        EXPECT_TRUE(value.has_value() && std::abs(*value - 1.0) < 1e-9,
                    "checked skew distance is 1");
    }

    lamina::LineSymbolic parallel_l2{{zero, one, zero}, {one, zero, zero}};
    auto parallel_distance = lamina::skew_lines_distance_checked(l1, parallel_l2);
    EXPECT_TRUE(!parallel_distance &&
                    parallel_distance.error().code == lamina::CasErrc::DomainError,
                "checked skew_lines_distance rejects parallel directions");

    auto constructed_line = lamina::line_from_two_points_checked(
        {zero, zero, zero}, {one, two, three});
    EXPECT_TRUE(constructed_line.has_value(), "checked line_from_two_points succeeds");
    if (constructed_line) {
        EXPECT_EQ_EXPR_STR(constructed_line.value().direction[2]->simplify(), "3",
                           "checked constructed line z direction = 3");
    }

    auto identical_points = lamina::line_from_two_points_checked(
        {one, one, one}, {one, one, one});
    EXPECT_TRUE(!identical_points &&
                    identical_points.error().code == lamina::CasErrc::DomainError,
                "checked line_from_two_points rejects identical points");

    auto checked_plane = lamina::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {zero, one, zero});
    EXPECT_TRUE(checked_plane.has_value(), "checked plane_from_three_points succeeds");
    if (checked_plane) {
        EXPECT_EQ_EXPR_STR(checked_plane.value().normal[2]->simplify(), "1",
                           "checked plane normal z = 1");
    }

    auto collinear_plane = lamina::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {two, zero, zero});
    EXPECT_TRUE(!collinear_plane &&
                    collinear_plane.error().code == lamina::CasErrc::DomainError,
                "checked plane_from_three_points rejects collinear points");

    auto angle = lamina::dihedral_angle_checked(
        lamina::PlaneSymbolic{{one, zero, zero}, zero},
        lamina::PlaneSymbolic{{zero, one, zero}, zero});
    EXPECT_TRUE(angle.has_value(), "checked dihedral_angle succeeds");
    if (angle) {
        auto value = test_numeric_eval(angle.value()->simplify());
        EXPECT_TRUE(value.has_value() && std::abs(*value - LMMC_CONST_PI / 2.0) < 1e-6,
                    "checked dihedral angle is pi/2");
    }

    lamina::CancellationToken token;
    token.cancel();
    lamina::ComputationContext context({}, token);
    auto cancelled = lamina::point_plane_distance_checked({one, two, three}, plane, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == lamina::CasErrc::Cancelled,
                "checked point_plane_distance observes cancellation");
}

void test_geometry_extensions() {
    using lamina::SurfaceSymbolic;
    auto N = [](int n){ return SymbolicExpr::number(n); };

    TEST_CASE("line_from_two_points: direction = p2 - p1");
    {
        std::vector<std::shared_ptr<SymbolicExpr>> p1 = {N(0), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p2 = {N(1), N(2), N(3)};
        auto line = lamina::line_from_two_points(p1, p2);
        EXPECT_EQ_EXPR_STR(line.direction[0]->simplify(), "1", "dir.x = 1");
        EXPECT_EQ_EXPR_STR(line.direction[1]->simplify(), "2", "dir.y = 2");
        EXPECT_EQ_EXPR_STR(line.direction[2]->simplify(), "3", "dir.z = 3");
    }

    TEST_CASE("plane_from_three_points: xy-plane normal is (0,0,c)");
    {
        std::vector<std::shared_ptr<SymbolicExpr>> p1 = {N(0), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p2 = {N(1), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p3 = {N(0), N(1), N(0)};
        auto plane = lamina::plane_from_three_points(p1, p2, p3);
        // normal = (p2-p1)x(p3-p1) = (1,0,0)x(0,1,0) = (0,0,1)
        EXPECT_EQ_EXPR_STR(plane.normal[0]->simplify(), "0", "n.x = 0");
        EXPECT_EQ_EXPR_STR(plane.normal[1]->simplify(), "0", "n.y = 0");
        EXPECT_EQ_EXPR_STR(plane.normal[2]->simplify(), "1", "n.z = 1");
    }

    TEST_CASE("classify_quadric: unit sphere");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        // x^2 + y^2 + z^2 - 1 = 0
        auto F = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
            SymbolicExpr::add(SymbolicExpr::multiply(z, z), N(-1)));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};
        std::string c = lamina::classify_quadric(surf);
        EXPECT_TRUE(c == "sphere", "x^2+y^2+z^2=1 classified as sphere");

        auto checked = lamina::classify_quadric_checked(surf);
        EXPECT_TRUE(checked.has_value(), "checked classify_quadric succeeds");
        if (checked) {
            EXPECT_TRUE(checked.value() == "sphere",
                        "checked classify_quadric reports sphere");
        }
    }

    TEST_CASE("classify_quadric_checked: explicit failures and unknowns");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto a = SymbolicExpr::variable("a");

        SurfaceSymbolic invalid{nullptr, {"x", "y", "z"}};
        auto invalid_result = lamina::classify_quadric_checked(invalid);
        EXPECT_TRUE(!invalid_result &&
                        invalid_result.error().code == lamina::CasErrc::InvalidArgument,
                    "checked classify_quadric rejects null surface equation");

        auto linear = SymbolicExpr::add(x, y);
        SurfaceSymbolic non_quadric{linear, {"x", "y", "z"}};
        auto unknown = lamina::classify_quadric_checked(non_quadric);
        EXPECT_TRUE(unknown.has_value() && unknown.value() == "unknown",
                    "checked classify_quadric returns unknown for supported non-quadric");

        auto symbolic_coeff = SymbolicExpr::add(
            SymbolicExpr::multiply(a, SymbolicExpr::multiply(x, x)),
            SymbolicExpr::add(SymbolicExpr::multiply(y, y),
                              SymbolicExpr::multiply(z, z)));
        SurfaceSymbolic symbolic{symbolic_coeff, {"x", "y", "z"}};
        auto unsupported = lamina::classify_quadric_checked(symbolic);
        EXPECT_TRUE(!unsupported &&
                        unsupported.error().code == lamina::CasErrc::Inconclusive,
                    "checked classify_quadric rejects unproved symbolic coefficients");
        EXPECT_TRUE(lamina::classify_quadric(symbolic) == "unknown",
                    "legacy classify_quadric unwraps checked failure to unknown");

        lamina::CancellationToken token;
        token.cancel();
        lamina::ComputationContext context({}, token);
        auto cancelled = lamina::classify_quadric_checked(non_quadric, context);
        EXPECT_TRUE(!cancelled &&
                        cancelled.error().code == lamina::CasErrc::Cancelled,
                    "checked classify_quadric observes cancellation");
    }

    TEST_CASE("surface_normal / tangent_plane of sphere at (1,0,0)");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto F = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
            SymbolicExpr::add(SymbolicExpr::multiply(z, z), N(-1)));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};
        std::vector<std::shared_ptr<SymbolicExpr>> pt = {N(1), N(0), N(0)};
        auto n = lamina::surface_normal(surf, pt);
        // grad = (2x,2y,2z) at (1,0,0) = (2,0,0); normalized = (1,0,0)
        EXPECT_EQ_EXPR_STR(n[0]->simplify(), "1", "unit normal x = 1");
        EXPECT_EQ_EXPR_STR(n[1]->simplify(), "0", "unit normal y = 0");

        auto checked_normal = lamina::surface_normal_checked(surf, pt);
        EXPECT_TRUE(checked_normal.has_value(), "checked surface normal succeeds");
        if (checked_normal) {
            EXPECT_EQ_EXPR_STR(checked_normal.value()[0]->simplify(), "1",
                               "checked unit normal x = 1");
            EXPECT_EQ_EXPR_STR(checked_normal.value()[1]->simplify(), "0",
                               "checked unit normal y = 0");
        }

        auto checked_plane = lamina::tangent_plane_checked(surf, pt);
        EXPECT_TRUE(checked_plane.has_value(), "checked tangent plane succeeds");
        if (checked_plane) {
            EXPECT_EQ_EXPR_STR(checked_plane.value().normal[0]->simplify(), "2",
                               "checked tangent plane normal x = 2");
            EXPECT_EQ_EXPR_STR(checked_plane.value().d->simplify(), "2",
                               "checked tangent plane d = 2");
        }
    }

    TEST_CASE("surface_normal_checked / tangent_plane_checked: explicit failures");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto zero = N(0);

        auto singular_F = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y)),
            SymbolicExpr::multiply(z, z));
        SurfaceSymbolic singular{singular_F, {"x", "y", "z"}};
        std::vector<std::shared_ptr<SymbolicExpr>> origin = {zero, zero, zero};
        auto singular_normal = lamina::surface_normal_checked(singular, origin);
        EXPECT_TRUE(!singular_normal.has_value(),
                    "checked surface normal rejects singular point");
        EXPECT_TRUE(singular_normal.error().code == lamina::CasErrc::DomainError,
                    "checked surface normal reports DomainError at singular point");

        auto singular_plane = lamina::tangent_plane_checked(singular, origin);
        EXPECT_TRUE(!singular_plane.has_value(),
                    "checked tangent plane rejects singular point");
        EXPECT_TRUE(singular_plane.error().code == lamina::CasErrc::DomainError,
                    "checked tangent plane reports DomainError at singular point");

        auto unsupported_F = SymbolicExpr::eq(x, zero);
        SurfaceSymbolic unsupported{unsupported_F, {"x", "y", "z"}};
        auto unsupported_normal = lamina::surface_normal_checked(unsupported, origin);
        EXPECT_TRUE(!unsupported_normal.has_value(),
                    "checked surface normal rejects unsupported derivatives");
        EXPECT_TRUE(unsupported_normal.error().code == lamina::CasErrc::Inconclusive,
                    "checked surface normal reports Inconclusive for unsupported derivatives");

        lamina::CancellationToken token;
        token.cancel();
        lamina::ComputationContext context({}, token);
        auto cancelled = lamina::tangent_plane_checked(singular, origin, context);
        EXPECT_TRUE(!cancelled.has_value(),
                    "checked tangent plane observes cancellation");
        EXPECT_TRUE(cancelled.error().code == lamina::CasErrc::Cancelled,
                    "checked tangent plane reports Cancelled");
    }

    TEST_CASE("dihedral_angle: perpendicular planes = pi/2");
    {
        lamina::PlaneSymbolic p1{{N(1), N(0), N(0)}, N(0)};
        lamina::PlaneSymbolic p2{{N(0), N(1), N(0)}, N(0)};
        auto ang = lamina::dihedral_angle(p1, p2);
        auto v = test_numeric_eval(ang ? ang->simplify() : nullptr);
        EXPECT_TRUE(v.has_value() && std::abs(*v - LMMC_CONST_PI/2.0) < 1e-6,
            "dihedral angle of perpendicular planes is pi/2");
    }
}

int main() {
    try {
        test_vector_dot();
        test_vector_cross();
        test_vector_angle();
        test_vector_checked_contracts();
        test_line_plane_intersection();
        test_point_plane_distance();
        test_skew_lines_distance();
        test_line_plane_checked_contracts();
        test_geometry_extensions();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
