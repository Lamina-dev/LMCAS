#include "test_common.hpp"
#include "symbolic_implicit_diff.hpp"

void test_implicit_diff_circle() {
    TEST_CASE("Implicit Diff: Circle x^2 + y^2 - r^2 = 0 => dy/dx = -x/y");
    {
        // F(x,y) = x^2 + y^2 - r^2
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto r = SymbolicExpr::variable("r");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto r2 = SymbolicExpr::power(r, SymbolicExpr::number(2));
        auto neg_r2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), r2);
        auto F = SymbolicExpr::add(SymbolicExpr::add(x2, y2), neg_r2);

        auto result = lamina::implicit_diff(F, "x", "y");
        std::string s = result ? result->to_string() : "null";
        // dy/dx = -F_x / F_y = -(2x) / (2y) = -x/y
        // The result should contain x and y, and represent -x/y (possibly unsimplified as (-1*2*x)/(2*y))
        EXPECT_CONTAINS(s, {"x", "y"}, "circle dy/dx contains x and y");
    }
}

void test_implicit_diff_ellipse() {
    TEST_CASE("Implicit Diff: Ellipse x^2/a^2 + y^2/b^2 - 1 = 0");
    {
        // F(x,y) = x^2/a^2 + y^2/b^2 - 1
        // dy/dx = -F_x / F_y = -(2x/a^2) / (2y/b^2) = -(b^2 * x) / (a^2 * y)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto a = SymbolicExpr::variable("a");
        auto b = SymbolicExpr::variable("b");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto a2 = SymbolicExpr::power(a, SymbolicExpr::number(2));
        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto term1 = SymbolicExpr::divide(x2, a2);
        auto term2 = SymbolicExpr::divide(y2, b2);
        auto neg_one = SymbolicExpr::number(-1);
        auto F = SymbolicExpr::add(SymbolicExpr::add(term1, term2), neg_one);

        auto result = lamina::implicit_diff(F, "x", "y");
        std::string s = result ? result->to_string() : "null";
        // Result should contain x, y, a, b representing -(b^2*x)/(a^2*y)
        EXPECT_CONTAINS(s, {"x", "y"}, "ellipse dy/dx contains x and y");
        EXPECT_CONTAINS(s, {"a", "b"}, "ellipse dy/dx contains a and b");
    }
}

void test_implicit_diff_polynomial() {
    TEST_CASE("Implicit Diff: Polynomial x^3 + y^3 - 3xy = 0");
    {
        // F(x,y) = x^3 + y^3 - 3xy
        // F_x = 3x^2 - 3y
        // F_y = 3y^2 - 3x
        // dy/dx = -(3x^2 - 3y) / (3y^2 - 3x) = -(x^2 - y) / (y^2 - x)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto y3 = SymbolicExpr::power(y, SymbolicExpr::number(3));
        auto three = SymbolicExpr::number(3);
        auto xy = SymbolicExpr::multiply(x, y);
        auto three_xy = SymbolicExpr::multiply(three, xy);
        auto neg_three_xy = SymbolicExpr::multiply(SymbolicExpr::number(-1), three_xy);
        auto F = SymbolicExpr::add(SymbolicExpr::add(x3, y3), neg_three_xy);

        auto result = lamina::implicit_diff(F, "x", "y");
        std::string s = result ? result->to_string() : "null";
        // Result should contain x and y terms
        EXPECT_CONTAINS(s, {"x", "y"}, "polynomial dy/dx contains x and y");
    }
}

void test_implicit_diff_transcendental() {
    TEST_CASE("Implicit Diff: Transcendental sin(x) + y^2 = 0");
    {
        // F(x,y) = sin(x) + y^2
        // F_x = cos(x)
        // F_y = 2y
        // dy/dx = -cos(x) / (2y)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto sin_x = SymbolicExpr::sin(x);
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto F = SymbolicExpr::add(sin_x, y2);

        auto result = lamina::implicit_diff(F, "x", "y");
        std::string s = result ? result->to_string() : "null";
        // Result should contain cos (from derivative of sin(x)) and y
        EXPECT_CONTAINS(s, {"cos", "y"}, "transcendental dy/dx contains cos and y");
    }
}

int main() {
    try {
        test_implicit_diff_circle();
        test_implicit_diff_ellipse();
        test_implicit_diff_polynomial();
        test_implicit_diff_transcendental();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
