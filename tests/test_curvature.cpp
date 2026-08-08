#include "test_common.hpp"
#include "calculus_utils.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace lamina;

int main()
{
    auto x = SymbolicExpr::variable("x");
    auto t = SymbolicExpr::variable("t");


    TEST_CASE("curvature: y = x^2 at x=0 should be 2");
    {
        // f(x) = x^2, f'=2x, f''=2
        // κ = |2| / (1 + (2*0)^2)^(3/2) = 2 / 1 = 2
        auto f = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto kappa = curvature(f, "x");
        EXPECT_TRUE(kappa != nullptr, "curvature(x^2) should not be null");

        if (kappa) {
            // Evaluate at x=0
            auto at_zero = kappa->substitute("x", SymbolicExpr::number(0));
            auto simplified = at_zero->simplify();
            auto val = simplified ? simplified : at_zero;
            std::cout << "  curvature(x^2) at x=0 = " << val->to_string() << std::endl;
            double num_val = val->to_numeric();
            EXPECT_NEAR(num_val, 2.0, 1e-9, "curvature of x^2 at x=0 is 2");
        }
    }

    TEST_CASE("curvature: y = x (straight line) should be 0");
    {
        // f(x) = x, f'=1, f''=0
        // κ = |0| / (1 + 1)^(3/2) = 0
        auto f = x;
        auto kappa = curvature(f, "x");
        EXPECT_TRUE(kappa != nullptr, "curvature(x) should not be null");

        if (kappa) {
            auto at_one = kappa->substitute("x", SymbolicExpr::number(1));
            auto simplified = at_one->simplify();
            auto val = simplified ? simplified : at_one;
            std::cout << "  curvature(x) at x=1 = " << val->to_string() << std::endl;
            double num_val = val->to_numeric();
            EXPECT_NEAR(num_val, 0.0, 1e-9, "curvature of straight line is 0");
        }
    }

    TEST_CASE("curvature_parametric: circle (cos(t), sin(t)) should be 1");
    {
        // x(t) = cos(t), y(t) = sin(t)
        // x' = -sin(t), x'' = -cos(t)
        // y' = cos(t), y'' = -sin(t)
        // |x'y'' - y'x''| = |(-sin)(-sin) - (cos)(-cos)| = |sin^2 + cos^2| = 1
        // (x'^2 + y'^2)^(3/2) = (sin^2 + cos^2)^(3/2) = 1
        // κ = 1
        auto x_t = SymbolicExpr::cos(t);
        auto y_t = SymbolicExpr::sin(t);
        auto kappa = curvature_parametric(x_t, y_t, "t");
        EXPECT_TRUE(kappa != nullptr, "curvature_parametric(cos,sin) should not be null");

        if (kappa) {
            // Evaluate at t = 0
            auto at_zero = kappa->substitute("t", SymbolicExpr::number(0));
            auto simplified = at_zero->simplify();
            auto val = simplified ? simplified : at_zero;
            std::cout << "  curvature_parametric(cos,sin) at t=0 = " << val->to_string() << std::endl;
            double num_val = val->to_numeric();
            EXPECT_NEAR(num_val, 1.0, 1e-9, "curvature of unit circle is 1");
        }
    }

    TEST_CASE("curvature_parametric: circle radius 2 should be 1/2");
    {
        // x(t) = 2*cos(t), y(t) = 2*sin(t)
        // κ = 1/R = 1/2
        auto x_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::cos(t));
        auto y_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::sin(t));
        auto kappa = curvature_parametric(x_t, y_t, "t");
        EXPECT_TRUE(kappa != nullptr, "curvature_parametric(2cos,2sin) should not be null");

        if (kappa) {
            auto at_zero = kappa->substitute("t", SymbolicExpr::number(0));
            auto simplified = at_zero->simplify();
            auto val = simplified ? simplified : at_zero;
            std::cout << "  curvature_parametric(2cos,2sin) at t=0 = " << val->to_string() << std::endl;
            double num_val = val->to_numeric();
            EXPECT_NEAR(num_val, 0.5, 1e-9, "curvature of circle radius 2 is 1/2");
        }
    }


    TEST_CASE("surface_area_revolution_x: sphere from y=sqrt(r^2-x^2)");
    {
        // Revolving y = sqrt(1-x^2) around x-axis from -1 to 1 gives sphere of radius 1
        // Surface area = 4π
        // f(x) = sqrt(1-x^2), f'(x) = -x/sqrt(1-x^2)
        // 1 + f'^2 = 1 + x^2/(1-x^2) = 1/(1-x^2)
        // sqrt(1+f'^2) = 1/sqrt(1-x^2)
        // |f| * sqrt(1+f'^2) = sqrt(1-x^2) * 1/sqrt(1-x^2) = 1
        // S = 2π * ∫_{-1}^{1} 1 dx = 2π * 2 = 4π
        auto one_minus_x2 = SymbolicExpr::add(
            SymbolicExpr::number(1),
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                   SymbolicExpr::power(x, SymbolicExpr::number(2))));
        auto f = SymbolicExpr::sqrt(one_minus_x2);
        auto a = SymbolicExpr::number(-1);
        auto b = SymbolicExpr::number(1);

        auto sa = surface_area_revolution_x(f, "x", a, b);
        EXPECT_TRUE(sa != nullptr, "surface_area_revolution_x(sqrt(1-x^2)) should not be null");
        std::cout << "  surface_area_revolution_x(sqrt(1-x^2), -1, 1) = " << sa->to_string() << std::endl;
        double num_val = sa->to_numeric();
        // Integrated via numerical fallback (Simpson), so use a looser tolerance.
        EXPECT_NEAR(num_val, 4.0 * M_PI, 0.05, "surface area of unit sphere is 4*pi");
    }

    TEST_CASE("surface_area_revolution_x: cone from y=x, [0,1]");
    {
        // Revolving y = x around x-axis from 0 to 1
        // f(x) = x, f'(x) = 1
        // S = 2π ∫₀¹ |x| * sqrt(1+1) dx = 2π * sqrt(2) * ∫₀¹ x dx
        //   = 2π * sqrt(2) * 1/2 = π*sqrt(2)
        auto f = x;
        auto a = SymbolicExpr::number(0);
        auto b = SymbolicExpr::number(1);

        auto sa = surface_area_revolution_x(f, "x", a, b);
        EXPECT_TRUE(sa != nullptr, "surface_area_revolution_x(x, 0, 1) should not be null");
        std::cout << "  surface_area_revolution_x(x, 0, 1) = " << sa->to_string() << std::endl;
        double num_val = sa->to_numeric();
        double expected = M_PI * std::sqrt(2.0);
        EXPECT_NEAR(num_val, expected, 0.01,
                    "surface area of cone y=x from 0 to 1 is pi*sqrt(2)");
    }

    TEST_CASE("surface_area_revolution_y: cone from y=x, [0,1]");
    {
        // Revolving y = x around y-axis from 0 to 1
        // f(x) = x, f'(x) = 1
        // S = 2π ∫₀¹ |x| * sqrt(1+1) dx = 2π * sqrt(2) * ∫₀¹ x dx
        //   = 2π * sqrt(2) * 1/2 = π*sqrt(2)
        auto f = x;
        auto a = SymbolicExpr::number(0);
        auto b = SymbolicExpr::number(1);

        auto sa = surface_area_revolution_y(f, "x", a, b);
        EXPECT_TRUE(sa != nullptr, "surface_area_revolution_y(x, 0, 1) should not be null");
        std::cout << "  surface_area_revolution_y(x, 0, 1) = " << sa->to_string() << std::endl;
        double num_val = sa->to_numeric();
        double expected = M_PI * std::sqrt(2.0);
        EXPECT_NEAR(num_val, expected, 0.01,
                    "surface area of cone y=x around y-axis from 0 to 1 is pi*sqrt(2)");
    }

    return TEST_REPORT();
}
