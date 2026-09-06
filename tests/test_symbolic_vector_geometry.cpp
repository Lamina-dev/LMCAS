#include "test_common.hpp"
#include "symbolic_vector_geometry.hpp"
#include <cfloat>
#include <cmath>

using namespace LMCAS;

void test_vector_dot() {
    TEST_CASE("Vector Dot: Orthogonal vectors (expect 0)");
    {
        // i . j = 0
        auto one = SymbolicExpr::number(1);
        auto zero = SymbolicExpr::number(0);
        std::vector<std::shared_ptr<SymbolicExpr>> a = {one, zero, zero};
        std::vector<std::shared_ptr<SymbolicExpr>> b = {zero, one, zero};
        auto result = LMCAS::vector_dot(a, b);
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
        auto result = LMCAS::vector_dot(
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
        auto result = LMCAS::vector_dot(
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
        auto result = LMCAS::vector_cross(
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
        auto result = LMCAS::vector_cross(
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
        auto result = LMCAS::vector_cross(
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
        double angle = LMCAS::vector_angle_checked({one, zero, zero}, {zero, one, zero}).value();
        std::cout << "  angle(i, j) = " << angle << " (expected " << M_PI / 2.0 << ")" << std::endl;
        EXPECT_TRUE(std::abs(angle - M_PI / 2.0) < 1e-9, "orthogonal angle is pi/2");
    }

    TEST_CASE("Vector Angle: Parallel vectors (0)");
    {
        // angle between (1,2,3) and (2,4,6) = 0
        double angle = LMCAS::vector_angle_checked({SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)}, {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}).value();
        std::cout << "  angle((1,2,3), (2,4,6)) = " << angle << " (expected 0)" << std::endl;
        EXPECT_TRUE(std::abs(angle - 0.0) < 1e-9, "parallel angle is 0");
    }

    TEST_CASE("Vector Angle: Extreme finite scale");
    {
        const auto large = SymbolicExpr::number(1.0e200);
        const auto negative_large = SymbolicExpr::number(-1.0e200);
        auto angle = LMCAS::vector_angle_checked(
            {large, large}, {large, negative_large});
        EXPECT_TRUE(angle.has_value(),
                    "extreme-scale finite vectors have a defined angle");
        if (angle) {
            EXPECT_TRUE(std::isfinite(angle.value()),
                        "extreme-scale vector angle remains finite");
            EXPECT_TRUE(std::abs(angle.value() - M_PI / 2.0) < 1e-12,
                        "extreme-scale orthogonal vectors retain pi/2 angle");
        }
    }
}

void test_vector_checked_contracts() {
    TEST_CASE("Vector Geometry Checked API Contracts");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto dot = LMCAS::vector_dot_checked({one, zero}, {zero, one});
    EXPECT_TRUE(dot.has_value(), "checked vector_dot succeeds");
    if (dot) {
        EXPECT_EQ_EXPR_STR(dot.value()->simplify(), "0",
                           "checked dot product of orthogonal vectors is 0");
    }

    auto bad_dot = LMCAS::vector_dot_checked({one}, {one, zero});
    EXPECT_TRUE(!bad_dot &&
                    bad_dot.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked vector_dot rejects dimension mismatch");
    bool legacy_dot_threw = false;
    try {
        (void)LMCAS::vector_dot({one}, {one, zero});
    } catch (const std::invalid_argument&) {
        legacy_dot_threw = true;
    }
    EXPECT_TRUE(legacy_dot_threw,
                "legacy vector_dot preserves invalid_argument for dimension mismatch");

    auto bad_component = LMCAS::vector_cross_checked({one, nullptr, zero},
                                                     {zero, one, zero});
    EXPECT_TRUE(!bad_component &&
                    bad_component.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked vector_cross rejects null component");

    auto bad_cross_dim = LMCAS::vector_cross_checked({one, zero},
                                                     {zero, one});
    EXPECT_TRUE(!bad_cross_dim &&
                    bad_cross_dim.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked vector_cross rejects non-3D vectors");

    auto angle = LMCAS::vector_angle_checked({one, zero}, {zero, one});
    EXPECT_TRUE(angle.has_value(), "checked vector_angle succeeds");
    if (angle) {
        EXPECT_TRUE(std::abs(angle.value() - M_PI / 2.0) < 1e-9,
                    "checked vector_angle of orthogonal vectors is pi/2");
    }

    auto zero_angle = LMCAS::vector_angle_checked({zero, zero}, {one, zero});
    EXPECT_TRUE(!zero_angle &&
                    zero_angle.error().code == LMCAS::CasErrc::DomainError,
                "checked vector_angle rejects zero-length vectors");

    auto x = SymbolicExpr::variable("x");
    auto symbolic_angle = LMCAS::vector_angle_checked({x, zero}, {one, zero});
    EXPECT_TRUE(!symbolic_angle &&
                    symbolic_angle.error().code == LMCAS::CasErrc::NumericFailure,
                "checked vector_angle rejects symbolic components");

    auto two_as_expr = SymbolicExpr::add(one, one);
    auto expression_angle = LMCAS::vector_angle_checked(
        {two_as_expr, zero}, {SymbolicExpr::number(2), zero});
    EXPECT_TRUE(expression_angle.has_value(),
                "checked vector_angle accepts finite numeric expressions");
    if (expression_angle) {
        EXPECT_TRUE(std::abs(expression_angle.value()) < 1e-9,
                    "checked vector_angle evaluates expression components");
    }

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::vector_dot_checked({one, zero}, {zero, one}, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == LMCAS::CasErrc::Cancelled,
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

        LMCAS::LineSymbolic line;
        line.point = {zero, zero, zero};
        line.direction = {one, one, one};

        LMCAS::PlaneSymbolic plane;
        plane.normal = {one, one, one};
        plane.d = three;

        auto result = LMCAS::line_plane_intersection_checked(line, plane).value();
        EXPECT_TRUE(result.size() == 3, "intersection has 3 coordinates");
        for (int i = 0; i < 3; ++i) {
            auto s = result[i] ? result[i]->simplify() : nullptr;
            std::string str = s ? s->to_string() : "null";
            std::cout << "  intersection[" << i << "] = " << str << std::endl;
            EXPECT_EQ_EXPR_STR(s, "1", "intersection coord " + std::to_string(i) + " is 1");
        }
    }

    TEST_CASE("Line-plane intersection tolerates tiny finite coefficient scales");
    {
        auto zero = SymbolicExpr::number(0);
        auto tiny = SymbolicExpr::number(1.0e-200);
        LMCAS::LineSymbolic line{
            {zero, zero, zero}, {tiny, zero, zero}};
        LMCAS::PlaneSymbolic plane{{tiny, zero, zero}, tiny};

        auto result = LMCAS::line_plane_intersection_checked(line, plane);
        EXPECT_TRUE(result.has_value(),
                    "nonparallel tiny directions retain a unique intersection");
        if (result) {
            auto x = test_numeric_eval(result.value()[0]);
            EXPECT_TRUE(x.has_value() && std::abs(*x - 1.0) < 1e-12,
                        "line and plane scales cancel from the intersection");
        }
    }

    TEST_CASE("Line-plane intersection scales cancelling extreme point dots");
    {
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
        EXPECT_TRUE(result.has_value(),
                    "a finite on-plane point survives cancelling dot products");
        if (result) {
            auto x = test_numeric_eval(result.value()[0]);
            auto y = test_numeric_eval(result.value()[1]);
            EXPECT_TRUE(
                x.has_value() && y.has_value() &&
                    *x == std::numeric_limits<double>::max() &&
                    *y == -std::numeric_limits<double>::max(),
                "scaled intersection preserves the extreme on-plane point");
        }
    }

    TEST_CASE("Line-plane intersection scales cancelling coordinate updates");
    {
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
        EXPECT_TRUE(result.has_value(),
                    "finite intersection survives an overflowing displacement");
        if (result) {
            auto x = test_numeric_eval(result.value()[0]);
            EXPECT_TRUE(
                x.has_value() &&
                    *x == std::numeric_limits<double>::max(),
                "scaled coordinate update preserves the finite endpoint");
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

        LMCAS::PlaneSymbolic plane;
        plane.normal = {one, one, one};
        plane.d = zero;

        auto result = LMCAS::point_plane_distance_checked(point, plane).value();
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

    TEST_CASE("Point-plane distance tolerates huge finite normal scale");
    {
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);
        auto huge = SymbolicExpr::number(1e200);
        LMCAS::PlaneSymbolic plane{{huge, zero, zero}, zero};

        auto result = LMCAS::point_plane_distance_checked(
            {one, zero, zero}, plane);
        EXPECT_TRUE(result.has_value(),
                    "finite huge plane normal has a representable distance");
        if (result) {
            auto value = test_numeric_eval(result.value());
            EXPECT_TRUE(value.has_value() && std::abs(*value - 1.0) < 1e-12,
                        "huge normal scale cancels from point-plane distance");
        }
    }

    TEST_CASE("Point-plane distance scales cancelling extreme dot products");
    {
        auto zero = SymbolicExpr::number(0);
        auto two = SymbolicExpr::number(2);
        auto maximum = SymbolicExpr::number(
            std::numeric_limits<double>::max());
        auto negative_maximum = SymbolicExpr::number(
            -std::numeric_limits<double>::max());
        LMCAS::PlaneSymbolic plane{{two, two, zero}, zero};

        auto result = LMCAS::point_plane_distance_checked(
            {maximum, negative_maximum, zero}, plane);
        EXPECT_TRUE(result.has_value(),
                    "cancelling extreme dot products define a finite distance");
        if (result) {
            auto value = test_numeric_eval(result.value());
            EXPECT_TRUE(value.has_value() && *value == 0.0,
                        "scaled point-normal dot product preserves cancellation");
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

        LMCAS::LineSymbolic l1;
        l1.point = {zero, zero, zero};
        l1.direction = {one, zero, zero};

        LMCAS::LineSymbolic l2;
        l2.point = {zero, one, zero};
        l2.direction = {zero, zero, one};

        auto result = LMCAS::skew_lines_distance_checked(l1, l2).value();
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

    TEST_CASE("Skew-line distance tolerates huge finite direction scales");
    {
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);
        auto huge = SymbolicExpr::number(1e200);
        LMCAS::LineSymbolic l1{
            {zero, zero, zero}, {huge, zero, zero}};
        LMCAS::LineSymbolic l2{
            {zero, one, zero}, {zero, zero, huge}};

        auto result = LMCAS::skew_lines_distance_checked(l1, l2);
        EXPECT_TRUE(result.has_value(),
                    "finite huge directions define a representable skew distance");
        if (result) {
            auto value = test_numeric_eval(result.value());
            EXPECT_TRUE(value.has_value() && std::abs(*value - 1.0) < 1e-12,
                        "direction scales cancel from skew-line distance");
        }
    }

    TEST_CASE("Skew-line distance scales extreme finite point offsets");
    {
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
        EXPECT_TRUE(result.has_value(),
                    "extreme finite point offsets have a representable skew distance");
        if (result) {
            auto value = test_numeric_eval(result.value());
            const double expected =
                (std::numeric_limits<double>::max() * epsilon) * 2.0;
            EXPECT_TRUE(
                value.has_value() && std::isfinite(*value) &&
                    std::abs(*value / expected - 1.0) < 1e-12,
                "point-offset scaling preserves skew-line distance");
        }
    }

}

void test_line_plane_checked_contracts() {
    TEST_CASE("Line/Plane Vector Geometry Checked API Contracts");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    LMCAS::LineSymbolic line{{zero, zero, zero}, {one, one, one}};
    LMCAS::PlaneSymbolic plane{{one, one, one}, three};
    auto intersection = LMCAS::line_plane_intersection_checked(line, plane);
    EXPECT_TRUE(intersection.has_value(),
                "checked line_plane_intersection succeeds");
    if (intersection) {
        EXPECT_TRUE(intersection.value().size() == 3,
                    "checked intersection has three coordinates");
        EXPECT_EQ_EXPR_STR(intersection.value()[0]->simplify(), "1",
                           "checked intersection x = 1");
    }

    LMCAS::LineSymbolic parallel{{zero, zero, zero}, {one, zero, zero}};
    LMCAS::PlaneSymbolic z_plane{{zero, zero, one}, one};
    auto no_unique = LMCAS::line_plane_intersection_checked(parallel, z_plane);
    EXPECT_TRUE(!no_unique &&
                    no_unique.error().code == LMCAS::CasErrc::DomainError,
                "checked line_plane_intersection rejects parallel line-plane");

    LMCAS::LineSymbolic zero_direction{{zero, zero, zero}, {zero, zero, zero}};
    auto bad_line = LMCAS::line_plane_intersection_checked(zero_direction, plane);
    EXPECT_TRUE(!bad_line &&
                    bad_line.error().code == LMCAS::CasErrc::DomainError,
                "checked line_plane_intersection rejects zero direction");

    LMCAS::PlaneSymbolic zero_plane{{zero, zero, zero}, one};
    auto bad_distance = LMCAS::point_plane_distance_checked({one, two, three}, zero_plane);
    EXPECT_TRUE(!bad_distance &&
                    bad_distance.error().code == LMCAS::CasErrc::DomainError,
                "checked point_plane_distance rejects zero plane normal");

    LMCAS::LineSymbolic l1{{zero, zero, zero}, {one, zero, zero}};
    LMCAS::LineSymbolic l2{{zero, one, zero}, {zero, zero, one}};
    auto skew = LMCAS::skew_lines_distance_checked(l1, l2);
    EXPECT_TRUE(skew.has_value(), "checked skew_lines_distance succeeds");
    if (skew) {
        auto value = test_numeric_eval(skew.value()->simplify());
        EXPECT_TRUE(value.has_value() && std::abs(*value - 1.0) < 1e-9,
                    "checked skew distance is 1");
    }

    LMCAS::LineSymbolic parallel_l2{{zero, one, zero}, {one, zero, zero}};
    auto parallel_distance = LMCAS::skew_lines_distance_checked(l1, parallel_l2);
    EXPECT_TRUE(!parallel_distance &&
                    parallel_distance.error().code == LMCAS::CasErrc::DomainError,
                "checked skew_lines_distance rejects parallel directions");

    auto constructed_line = LMCAS::line_from_two_points_checked(
        {zero, zero, zero}, {one, two, three});
    EXPECT_TRUE(constructed_line.has_value(), "checked line_from_two_points succeeds");
    if (constructed_line) {
        EXPECT_EQ_EXPR_STR(constructed_line.value().direction[2]->simplify(), "3",
                           "checked constructed line z direction = 3");
    }

    {
        auto negative_max = SymbolicExpr::number(
            -std::numeric_limits<double>::max());
        auto positive_max = SymbolicExpr::number(
            std::numeric_limits<double>::max());
        auto extreme_line = LMCAS::line_from_two_points_checked(
            {negative_max, zero, zero}, {positive_max, zero, zero});
        EXPECT_TRUE(extreme_line.has_value(),
                    "opposite extreme finite points define a line");
        if (extreme_line) {
            auto dx = test_numeric_eval(extreme_line.value().direction[0]);
            EXPECT_TRUE(dx.has_value() && std::isfinite(*dx) && *dx > 0.0,
                        "line direction is scaled before point subtraction overflows");
        }
    }

    auto identical_points = LMCAS::line_from_two_points_checked(
        {one, one, one}, {one, one, one});
    EXPECT_TRUE(!identical_points &&
                    identical_points.error().code == LMCAS::CasErrc::DomainError,
                "checked line_from_two_points rejects identical points");

    auto checked_plane = LMCAS::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {zero, one, zero});
    EXPECT_TRUE(checked_plane.has_value(), "checked plane_from_three_points succeeds");
    if (checked_plane) {
        EXPECT_EQ_EXPR_STR(checked_plane.value().normal[2]->simplify(), "1",
                           "checked plane normal z = 1");
    }

    {
        auto huge = SymbolicExpr::number(1.0e200);
        auto huge_plane = LMCAS::plane_from_three_points_checked(
            {zero, zero, zero}, {huge, zero, zero}, {zero, huge, zero});
        EXPECT_TRUE(huge_plane.has_value(),
                    "finite huge edges define a representable plane direction");
        if (huge_plane) {
            auto nx = test_numeric_eval(huge_plane.value().normal[0]);
            auto ny = test_numeric_eval(huge_plane.value().normal[1]);
            auto nz = test_numeric_eval(huge_plane.value().normal[2]);
            EXPECT_TRUE(nx.has_value() && ny.has_value() && nz.has_value() &&
                            *nx == 0.0 && *ny == 0.0 && *nz > 0.0 &&
                            std::isfinite(*nz),
                        "plane construction scales edge vectors before crossing");
        }
    }

    {
        auto negative_max = SymbolicExpr::number(
            -std::numeric_limits<double>::max());
        auto positive_max = SymbolicExpr::number(
            std::numeric_limits<double>::max());
        auto extreme_plane = LMCAS::plane_from_three_points_checked(
            {negative_max, zero, zero},
            {positive_max, zero, zero},
            {negative_max, positive_max, zero});
        EXPECT_TRUE(extreme_plane.has_value(),
                    "extreme finite point differences define a plane");
        if (extreme_plane) {
            auto nz = test_numeric_eval(extreme_plane.value().normal[2]);
            EXPECT_TRUE(nz.has_value() && std::isfinite(*nz) && *nz > 0.0,
                        "plane edges are formed after common point scaling");
        }
    }

    auto collinear_plane = LMCAS::plane_from_three_points_checked(
        {zero, zero, zero}, {one, zero, zero}, {two, zero, zero});
    EXPECT_TRUE(!collinear_plane &&
                    collinear_plane.error().code == LMCAS::CasErrc::DomainError,
                "checked plane_from_three_points rejects collinear points");

    auto angle = LMCAS::dihedral_angle_checked(
        LMCAS::PlaneSymbolic{{one, zero, zero}, zero},
        LMCAS::PlaneSymbolic{{zero, one, zero}, zero});
    EXPECT_TRUE(angle.has_value(), "checked dihedral_angle succeeds");
    if (angle) {
        auto value = test_numeric_eval(angle.value()->simplify());
        EXPECT_TRUE(value.has_value() && std::abs(*value - LMMC_CONST_PI / 2.0) < 1e-6,
                    "checked dihedral angle is pi/2");
    }

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::point_plane_distance_checked({one, two, three}, plane, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == LMCAS::CasErrc::Cancelled,
                "checked point_plane_distance observes cancellation");
}

void test_geometry_extensions() {
    using LMCAS::SurfaceSymbolic;
    auto N = [](int n){ return SymbolicExpr::number(n); };

    TEST_CASE("line_from_two_points: direction = p2 - p1");
    {
        std::vector<std::shared_ptr<SymbolicExpr>> p1 = {N(0), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p2 = {N(1), N(2), N(3)};
        auto line = LMCAS::line_from_two_points_checked(p1, p2).value();
        EXPECT_EQ_EXPR_STR(line.direction[0]->simplify(), "1", "dir.x = 1");
        EXPECT_EQ_EXPR_STR(line.direction[1]->simplify(), "2", "dir.y = 2");
        EXPECT_EQ_EXPR_STR(line.direction[2]->simplify(), "3", "dir.z = 3");
    }

    TEST_CASE("plane_from_three_points: xy-plane normal is (0,0,c)");
    {
        std::vector<std::shared_ptr<SymbolicExpr>> p1 = {N(0), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p2 = {N(1), N(0), N(0)};
        std::vector<std::shared_ptr<SymbolicExpr>> p3 = {N(0), N(1), N(0)};
        auto plane = LMCAS::plane_from_three_points_checked(p1, p2, p3).value();
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
        std::string c = LMCAS::classify_quadric_checked(surf).value();
        EXPECT_TRUE(c == "sphere", "x^2+y^2+z^2=1 classified as sphere");

        auto checked = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(checked.has_value(), "checked classify_quadric succeeds");
        if (checked) {
            EXPECT_TRUE(checked.value() == "sphere",
                        "checked classify_quadric reports sphere");
        }
    }

    TEST_CASE("quadric classification is invariant under equation scaling");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto quadratic = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(y, y)),
            SymbolicExpr::add(
                SymbolicExpr::multiply(N(-1), SymbolicExpr::multiply(z, z)),
                N(-1)));
        auto F = SymbolicExpr::multiply(
            SymbolicExpr::number(1e-12), quadratic);
        SurfaceSymbolic surf{F, {"x", "y", "z"}};

        auto classification = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(classification.has_value() &&
                        classification.value() == "hyperboloid",
                    "nonzero equation scaling preserves hyperboloid classification");
    }

    TEST_CASE("quadric classification follows the unsquared axis");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto F = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(z, z)),
            SymbolicExpr::multiply(N(-1), y));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};

        auto classification = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(classification.has_value() &&
                        classification.value() == "paraboloid",
                    "x^2+z^2-y=0 is a paraboloid");
    }

    TEST_CASE("quadric null-axis detection uses componentwise error bounds");
    {
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
        EXPECT_TRUE(
            classification && classification.value() == "paraboloid",
            "a large range-space term does not hide the null-axis term");
    }

    TEST_CASE("translated cone classification completes the square");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto shifted_x_square = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::multiply(x, x),
                SymbolicExpr::multiply(N(-2), x)),
            N(1));
        auto F = SymbolicExpr::add(
            SymbolicExpr::add(shifted_x_square,
                              SymbolicExpr::multiply(y, y)),
            SymbolicExpr::multiply(N(-1), SymbolicExpr::multiply(z, z)));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};

        auto classification = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(classification.has_value() &&
                        classification.value() == "cone",
                    "(x-1)^2+y^2-z^2=0 is a translated cone");
    }

    TEST_CASE("quadric classifier reduces mixed terms to principal axes");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto xy = SymbolicExpr::multiply(x, y);
        auto rotated_ellipsoid = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(N(2), SymbolicExpr::multiply(x, x)),
                    SymbolicExpr::multiply(N(2), xy)),
                SymbolicExpr::multiply(N(2), SymbolicExpr::multiply(y, y))),
            SymbolicExpr::add(SymbolicExpr::multiply(z, z), N(-1)));
        SurfaceSymbolic ellipsoid{
            rotated_ellipsoid, {"x", "y", "z"}};

        auto ellipsoid_classification =
            LMCAS::classify_quadric_checked(ellipsoid);
        EXPECT_TRUE(
            ellipsoid_classification &&
                ellipsoid_classification.value() == "ellipsoid",
            "positive-definite mixed quadric is classified after rotation");

        auto rotated_paraboloid = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(x, x),
                    SymbolicExpr::multiply(N(2), xy)),
                SymbolicExpr::multiply(y, y)),
            SymbolicExpr::add(
                SymbolicExpr::multiply(z, z),
                SymbolicExpr::add(
                    SymbolicExpr::multiply(N(-1), x), y)));
        SurfaceSymbolic paraboloid{
            rotated_paraboloid, {"x", "y", "z"}};

        auto paraboloid_classification =
            LMCAS::classify_quadric_checked(paraboloid);
        EXPECT_TRUE(
            paraboloid_classification &&
                paraboloid_classification.value() == "paraboloid",
            "mixed rank-two quadric follows its rotated null axis");
    }

    TEST_CASE("rank-one quadric identifies a parabolic cylinder");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto F = SymbolicExpr::add(
            SymbolicExpr::multiply(x, x),
            SymbolicExpr::multiply(N(-1), y));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};

        auto classification = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(
            classification && classification.value() == "cylinder",
            "x^2-y=0 is a parabolic cylinder along the z-axis");
    }

    TEST_CASE("quadric classifier does not erase unresolved eigenvalues");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto F = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::multiply(x, x),
                SymbolicExpr::multiply(
                    SymbolicExpr::number(1e-16),
                    SymbolicExpr::multiply(y, y))),
            SymbolicExpr::add(SymbolicExpr::multiply(z, z), N(-1)));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};

        auto classification = LMCAS::classify_quadric_checked(surf);
        EXPECT_TRUE(
            !classification &&
                classification.error().code == LMCAS::CasErrc::Inconclusive,
            "an eigenvalue inside the backward-error bound is unresolved");
    }

    TEST_CASE("quadric classifier rejects empty and degenerate real loci");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto squares = SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                              SymbolicExpr::multiply(y, y)),
            SymbolicExpr::multiply(z, z));

        SurfaceSymbolic empty_sphere{
            SymbolicExpr::add(squares, N(1)), {"x", "y", "z"}};
        auto empty = LMCAS::classify_quadric_checked(empty_sphere);
        EXPECT_TRUE(!empty &&
                        empty.error().code == LMCAS::CasErrc::Inconclusive,
                    "x^2+y^2+z^2+1=0 is not reported as a real sphere");

        auto point_equation = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                                  SymbolicExpr::multiply(N(-2), x)),
                N(1)),
            SymbolicExpr::add(SymbolicExpr::multiply(y, y),
                              SymbolicExpr::multiply(z, z)));
        SurfaceSymbolic point_surface{
            point_equation, {"x", "y", "z"}};
        auto point = LMCAS::classify_quadric_checked(point_surface);
        EXPECT_TRUE(!point &&
                        point.error().code == LMCAS::CasErrc::Inconclusive,
                    "(x-1)^2+y^2+z^2=0 is a point, not a sphere");

        SurfaceSymbolic empty_cylinder{
            SymbolicExpr::add(
                SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                                  SymbolicExpr::multiply(y, y)),
                N(1)),
            {"x", "y", "z"}};
        auto empty_rank_two =
            LMCAS::classify_quadric_checked(empty_cylinder);
        EXPECT_TRUE(!empty_rank_two &&
                        empty_rank_two.error().code ==
                            LMCAS::CasErrc::Inconclusive,
                    "x^2+y^2+1=0 is not reported as a real cylinder");

        SurfaceSymbolic real_cylinder{
            SymbolicExpr::add(
                SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                                  SymbolicExpr::multiply(y, y)),
                N(-1)),
            {"x", "y", "z"}};
        auto real_rank_two =
            LMCAS::classify_quadric_checked(real_cylinder);
        EXPECT_TRUE(real_rank_two &&
                        real_rank_two.value() == "cylinder",
                    "x^2+y^2-1=0 remains a real cylinder");

        SurfaceSymbolic plane_pair{
            SymbolicExpr::add(
                SymbolicExpr::multiply(x, x),
                SymbolicExpr::multiply(
                    N(-1), SymbolicExpr::multiply(y, y))),
            {"x", "y", "z"}};
        auto degenerate_rank_two =
            LMCAS::classify_quadric_checked(plane_pair);
        EXPECT_TRUE(!degenerate_rank_two &&
                        degenerate_rank_two.error().code ==
                            LMCAS::CasErrc::Inconclusive,
                    "x^2-y^2=0 is a plane pair, not a cylinder");
    }

    TEST_CASE("classify_quadric_checked: explicit failures and unknowns");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto a = SymbolicExpr::variable("a");

        SurfaceSymbolic invalid{nullptr, {"x", "y", "z"}};
        auto invalid_result = LMCAS::classify_quadric_checked(invalid);
        EXPECT_TRUE(!invalid_result &&
                        invalid_result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked classify_quadric rejects null surface equation");

        auto linear = SymbolicExpr::add(x, y);
        SurfaceSymbolic non_quadric{linear, {"x", "y", "z"}};
        auto unknown = LMCAS::classify_quadric_checked(non_quadric);
        EXPECT_TRUE(!unknown &&
                        unknown.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked classify_quadric reports an unclassified surface as Inconclusive");

        auto symbolic_coeff = SymbolicExpr::add(
            SymbolicExpr::multiply(a, SymbolicExpr::multiply(x, x)),
            SymbolicExpr::add(SymbolicExpr::multiply(y, y),
                              SymbolicExpr::multiply(z, z)));
        SurfaceSymbolic symbolic{symbolic_coeff, {"x", "y", "z"}};
        auto unsupported = LMCAS::classify_quadric_checked(symbolic);
        EXPECT_TRUE(!unsupported &&
                        unsupported.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked classify_quadric rejects unproved symbolic coefficients");

        LMCAS::CancellationToken token;
        token.cancel();
        LMCAS::ComputationContext context({}, token);
        auto cancelled = LMCAS::classify_quadric_checked(non_quadric, context);
        EXPECT_TRUE(!cancelled &&
                        cancelled.error().code == LMCAS::CasErrc::Cancelled,
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
        auto n = LMCAS::surface_normal_checked(surf, pt).value();
        // grad = (2x,2y,2z) at (1,0,0) = (2,0,0); normalized = (1,0,0)
        EXPECT_EQ_EXPR_STR(n[0]->simplify(), "1", "unit normal x = 1");
        EXPECT_EQ_EXPR_STR(n[1]->simplify(), "0", "unit normal y = 0");

        auto checked_normal = LMCAS::surface_normal_checked(surf, pt);
        EXPECT_TRUE(checked_normal.has_value(), "checked surface normal succeeds");
        if (checked_normal) {
            EXPECT_EQ_EXPR_STR(checked_normal.value()[0]->simplify(), "1",
                               "checked unit normal x = 1");
            EXPECT_EQ_EXPR_STR(checked_normal.value()[1]->simplify(), "0",
                               "checked unit normal y = 0");
        }

        auto checked_plane = LMCAS::tangent_plane_checked(surf, pt);
        EXPECT_TRUE(checked_plane.has_value(), "checked tangent plane succeeds");
        if (checked_plane) {
            EXPECT_EQ_EXPR_STR(checked_plane.value().normal[0]->simplify(), "2",
                               "checked tangent plane normal x = 2");
            EXPECT_EQ_EXPR_STR(checked_plane.value().d->simplify(), "2",
                               "checked tangent plane d = 2");
        }
    }

    TEST_CASE("surface normal and tangent plane tolerate huge finite gradients");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto huge = SymbolicExpr::number(1e200);
        auto F = SymbolicExpr::add(
            SymbolicExpr::multiply(huge, x),
            SymbolicExpr::multiply(huge, y));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};
        std::vector<std::shared_ptr<SymbolicExpr>> origin = {N(0), N(0), N(0)};

        auto normal = LMCAS::surface_normal_checked(surf, origin);
        EXPECT_TRUE(normal.has_value(),
                    "finite huge gradient has a representable unit normal");
        if (normal) {
            auto nx = test_numeric_eval(normal.value()[0]);
            auto ny = test_numeric_eval(normal.value()[1]);
            EXPECT_TRUE(nx.has_value() && ny.has_value() &&
                            std::abs(*nx - std::sqrt(0.5)) < 1e-12 &&
                            std::abs(*ny - std::sqrt(0.5)) < 1e-12,
                        "huge gradient normal is scaled before normalization");
        }

        auto plane = LMCAS::tangent_plane_checked(surf, origin);
        EXPECT_TRUE(plane.has_value(),
                    "finite huge gradient defines a tangent plane");
        if (plane) {
            auto nx = test_numeric_eval(plane.value().normal[0]);
            auto ny = test_numeric_eval(plane.value().normal[1]);
            EXPECT_TRUE(nx.has_value() && ny.has_value() &&
                            *nx == 1e200 && *ny == 1e200,
                        "tangent plane preserves finite huge coefficients");
        }
    }

    TEST_CASE("tangent plane scales cancelling extreme point dots");
    {
        const double maximum = std::numeric_limits<double>::max();
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto z = SymbolicExpr::variable("z");
        auto F = SymbolicExpr::add(
            SymbolicExpr::multiply(N(2), x),
            SymbolicExpr::multiply(N(2), y));
        SurfaceSymbolic surf{F, {"x", "y", "z"}};
        std::vector<std::shared_ptr<SymbolicExpr>> point = {
            SymbolicExpr::number(maximum),
            SymbolicExpr::number(-maximum),
            N(0)};

        auto plane = LMCAS::tangent_plane_checked(surf, point);
        EXPECT_TRUE(plane.has_value(),
                    "finite tangent plane survives a cancelling gradient dot");
        if (plane) {
            auto d = test_numeric_eval(plane.value().d);
            EXPECT_TRUE(d.has_value() && *d == 0.0,
                        "cancelling extreme tangent-plane constant is zero");
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
        auto singular_normal = LMCAS::surface_normal_checked(singular, origin);
        EXPECT_TRUE(!singular_normal.has_value(),
                    "checked surface normal rejects singular point");
        EXPECT_TRUE(singular_normal.error().code == LMCAS::CasErrc::DomainError,
                    "checked surface normal reports DomainError at singular point");

        auto singular_plane = LMCAS::tangent_plane_checked(singular, origin);
        EXPECT_TRUE(!singular_plane.has_value(),
                    "checked tangent plane rejects singular point");
        EXPECT_TRUE(singular_plane.error().code == LMCAS::CasErrc::DomainError,
                    "checked tangent plane reports DomainError at singular point");

        auto unsupported_F = SymbolicExpr::eq(x, zero);
        SurfaceSymbolic unsupported{unsupported_F, {"x", "y", "z"}};
        auto unsupported_normal = LMCAS::surface_normal_checked(unsupported, origin);
        EXPECT_TRUE(!unsupported_normal.has_value(),
                    "checked surface normal rejects unsupported derivatives");
        EXPECT_TRUE(unsupported_normal.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked surface normal reports Inconclusive for unsupported derivatives");

        LMCAS::CancellationToken token;
        token.cancel();
        LMCAS::ComputationContext context({}, token);
        auto cancelled = LMCAS::tangent_plane_checked(singular, origin, context);
        EXPECT_TRUE(!cancelled.has_value(),
                    "checked tangent plane observes cancellation");
        EXPECT_TRUE(cancelled.error().code == LMCAS::CasErrc::Cancelled,
                    "checked tangent plane reports Cancelled");
    }

    TEST_CASE("dihedral_angle: perpendicular planes = pi/2");
    {
        LMCAS::PlaneSymbolic p1{{N(1), N(0), N(0)}, N(0)};
        LMCAS::PlaneSymbolic p2{{N(0), N(1), N(0)}, N(0)};
        auto ang = LMCAS::dihedral_angle_checked(p1, p2).value();
        auto v = test_numeric_eval(ang ? ang->simplify() : nullptr);
        EXPECT_TRUE(v.has_value() && std::abs(*v - LMMC_CONST_PI/2.0) < 1e-6,
            "dihedral angle of perpendicular planes is pi/2");
    }

    TEST_CASE("dihedral angle tolerates huge finite plane normals");
    {
        auto huge = SymbolicExpr::number(1e200);
        LMCAS::PlaneSymbolic p1{{huge, N(0), N(0)}, N(0)};
        LMCAS::PlaneSymbolic p2{{N(0), huge, N(0)}, N(0)};
        auto angle = LMCAS::dihedral_angle_checked(p1, p2);
        EXPECT_TRUE(angle.has_value(),
                    "finite huge plane normals remain valid");
        if (angle) {
            auto value = test_numeric_eval(angle.value());
            EXPECT_TRUE(value.has_value() &&
                            std::abs(*value - LMMC_CONST_PI / 2.0) < 1e-12,
                        "huge perpendicular normals retain pi/2 angle");
        }
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
