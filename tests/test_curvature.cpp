#include "test_common.hpp"
#include "calculus_utils.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace LMCAS;

TEST(Curvature, ParabolicCurvature) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto kappa = curvature_checked(f, "x").value();
    EXPECT_TRUE((kappa != nullptr)) << "curvature(x^2) should not be null";

    if (kappa) {
        auto at_zero = kappa->substitute("x", SymbolicExpr::number(0));
        auto simplified = at_zero->simplify();
        auto val = simplified ? simplified : at_zero;
        double num_val = val->to_numeric();
        EXPECT_NEAR(num_val, 2.0, 1e-9) << "curvature of x^2 at x=0 is 2";
    }
}

TEST(Curvature, LineCurvature) {
    auto x = SymbolicExpr::variable("x");

    auto f = x;
    auto kappa = curvature_checked(f, "x").value();
    EXPECT_TRUE((kappa != nullptr)) << "curvature(x) should not be null";

    if (kappa) {
        auto at_one = kappa->substitute("x", SymbolicExpr::number(1));
        auto simplified = at_one->simplify();
        auto val = simplified ? simplified : at_one;
        double num_val = val->to_numeric();
        EXPECT_NEAR(num_val, 0.0, 1e-9) << "curvature of straight line is 0";
    }
}

TEST(Curvature, UnitCircleCurvature) {
    auto t = SymbolicExpr::variable("t");

    /**
     * @brief 单位圆 x(t)=cos(t)、y(t)=sin(t) 的曲率为 1。
     * x'=-sin(t)，x''=-cos(t)，y'=cos(t)，y''=-sin(t)。
     * |x'y''-y'x''|=|sin^2(t)+cos^2(t)|=1，
     * (x'^2+y'^2)^(3/2)=(sin^2(t)+cos^2(t))^(3/2)=1，故 κ=1。
     */
    auto x_t = SymbolicExpr::cos(t);
    auto y_t = SymbolicExpr::sin(t);
    auto kappa = curvature_parametric_checked(x_t, y_t, "t").value();
    EXPECT_TRUE((kappa != nullptr)) << "curvature_parametric(cos,sin) should not be null";

    if (kappa) {
        auto at_zero = kappa->substitute("t", SymbolicExpr::number(0));
        auto simplified = at_zero->simplify();
        auto val = simplified ? simplified : at_zero;
        double num_val = val->to_numeric();
        EXPECT_NEAR(num_val, 1.0, 1e-9) << "curvature of unit circle is 1";
    }
}

TEST(Curvature, ScaledCircleCurvature) {
    auto t = SymbolicExpr::variable("t");

    auto x_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::cos(t));
    auto y_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::sin(t));
    auto kappa = curvature_parametric_checked(x_t, y_t, "t").value();
    EXPECT_TRUE((kappa != nullptr)) << "curvature_parametric(2cos,2sin) should not be null";

    if (kappa) {
        auto at_zero = kappa->substitute("t", SymbolicExpr::number(0));
        auto simplified = at_zero->simplify();
        auto val = simplified ? simplified : at_zero;
        double num_val = val->to_numeric();
        EXPECT_NEAR(num_val, 0.5, 1e-9) << "curvature of circle radius 2 is 1/2";
    }
}

TEST(Curvature, SphereSurfaceDomain) {
    auto x = SymbolicExpr::variable("x");

    /**
     * @brief y=sqrt(1-x^2) 在 [-1,1] 绕 x 轴旋转，得到面积为 4π 的单位球面。
     * f'=-x/sqrt(1-x^2)，1+f'^2=1+x^2/(1-x^2)=1/(1-x^2)，
     * sqrt(1+f'^2)=1/sqrt(1-x^2)，|f|*sqrt(1+f'^2)=1。
     * S=2π*∫_{-1}^{1}1 dx=2π*2=4π。
     */
    auto one_minus_x2 = SymbolicExpr::add(
        SymbolicExpr::number(1),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(x, SymbolicExpr::number(2))));
    auto f = SymbolicExpr::sqrt(one_minus_x2);
    auto a = SymbolicExpr::number(-1);
    auto b = SymbolicExpr::number(1);

    auto surface =
        surface_area_revolution_x_checked(f, "x", a, b);
    EXPECT_TRUE((!surface &&
                 surface.error().code == LMCAS::CasErrc::Inconclusive))
        << "sphere surface requiring implicit numeric fallback is Inconclusive";
}

TEST(Curvature, ConeSurfaceX) {
    auto x = SymbolicExpr::variable("x");

    /**
     * @brief y=x 在 [0,1] 绕 x 轴旋转的曲面面积为 π*sqrt(2)。
     * f(x)=x，f'=1，S=2π∫₀¹|x|*sqrt(1+1) dx
     * =2π*sqrt(2)*∫₀¹x dx=2π*sqrt(2)*1/2=π*sqrt(2)。
     */
    auto f = x;
    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto sa = surface_area_revolution_x_checked(f, "x", a, b).value();
    EXPECT_TRUE((sa != nullptr)) << "surface_area_revolution_x(x, 0, 1) should not be null";
    double num_val = sa->to_numeric();
    double expected = M_PI * std::sqrt(2.0);
    {
        const double actual_value = (num_val);
        const double expected_value = (expected);
        const double tolerance = (0.01);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(Curvature, ConeSurfaceY) {
    auto x = SymbolicExpr::variable("x");

    /**
     * @brief y=x 在 [0,1] 绕 y 轴旋转的曲面面积为 π*sqrt(2)。
     * f(x)=x，f'=1，S=2π∫₀¹|x|*sqrt(1+1) dx
     * =2π*sqrt(2)*∫₀¹x dx=2π*sqrt(2)*1/2=π*sqrt(2)。
     */
    auto f = x;
    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto sa = surface_area_revolution_y_checked(f, "x", a, b).value();
    EXPECT_TRUE((sa != nullptr)) << "surface_area_revolution_y(x, 0, 1) should not be null";
    double num_val = sa->to_numeric();
    double expected = M_PI * std::sqrt(2.0);
    {
        const double actual_value = (num_val);
        const double expected_value = (expected);
        const double tolerance = (0.01);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}
