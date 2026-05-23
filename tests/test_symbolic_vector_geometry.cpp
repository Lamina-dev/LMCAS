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

int main() {
    try {
        test_vector_dot();
        test_vector_cross();
        test_vector_angle();
        test_line_plane_intersection();
        test_point_plane_distance();
        test_skew_lines_distance();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
