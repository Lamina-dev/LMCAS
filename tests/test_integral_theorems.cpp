/**
 * @file test_integral_theorems.cpp
 * @brief 积分定理单元测试：格林定理、散度定理（高斯定理）、斯托克斯定理。
 */

#include "test_common.hpp"
#include "vector_calculus.hpp"
#include "internal/symbolic_ast.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>

using namespace LMCAS;

TEST(IntegralTheorems, GreensTheoremConstantField) {
    // ∮ P dx + Q dy = ∬ (∂Q/∂x - ∂P/∂y) dA
    // P = 0, Q = x => integrand = ∂x/∂x - 0 = 1
    // ∬ 1 dA over [0,1]x[0,1] = 1
    auto P = SymbolicExpr::number(0);
    auto Q = SymbolicExpr::variable("x");

    auto x_lo = SymbolicExpr::number(0);
    auto x_hi = SymbolicExpr::number(1);
    auto y_lo = SymbolicExpr::number(0);
    auto y_hi = SymbolicExpr::number(1);

    auto result = greens_theorem(P, Q, {"x", "y"},
                                 {x_lo, x_hi}, {y_lo, y_hi});

    EXPECT_TRUE((result != nullptr)) << "Green's theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-10) << "∬ 1 dA over [0,1]x[0,1] = 1";
        }
    }
}

TEST(IntegralTheorems, GreensTheoremLinearField) {
    // P = y, Q = x => ∂Q/∂x - ∂P/∂y = 1 - 1 = 0
    auto P = SymbolicExpr::variable("y");
    auto Q = SymbolicExpr::variable("x");

    auto x_lo = SymbolicExpr::number(0);
    auto x_hi = SymbolicExpr::number(2);
    auto y_lo = SymbolicExpr::number(0);
    auto y_hi = SymbolicExpr::number(3);

    auto result = greens_theorem(P, Q, {"x", "y"},
                                 {x_lo, x_hi}, {y_lo, y_hi});

    EXPECT_TRUE((result != nullptr)) << "Green's theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10) << "∬(∂x/∂x - ∂y/∂y) dA = ∬ 0 dA = 0";
        }
    }
}

TEST(IntegralTheorems, GreensTheoremQuadratic) {
    // P = -y^2, Q = x^2
    // ∂Q/∂x = 2x, ∂P/∂y = -2y
    // integrand = 2x - (-2y) = 2x + 2y
    // ∬(2x + 2y) dA over [0,1]x[0,1]
    // = ∫₀¹ ∫₀¹ (2x + 2y) dy dx
    // = ∫₀¹ [2xy + y²]₀¹ dx = ∫₀¹ (2x + 1) dx = [x² + x]₀¹ = 2
    auto y = SymbolicExpr::variable("y");
    auto x = SymbolicExpr::variable("x");

    auto P = SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                    SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto Q = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto x_lo = SymbolicExpr::number(0);
    auto x_hi = SymbolicExpr::number(1);
    auto y_lo = SymbolicExpr::number(0);
    auto y_hi = SymbolicExpr::number(1);

    auto result = greens_theorem(P, Q, {"x", "y"},
                                 {x_lo, x_hi}, {y_lo, y_hi});

    EXPECT_TRUE((result != nullptr)) << "Green's theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 2.0, 1e-10) << "∬(2x+2y) dA over [0,1]x[0,1] = 2";
        }
    }
}

TEST(IntegralTheorems, GreensAreaUnitSquare) {
    // For a unit circle: r(t) = (cos(t), sin(t)), t ∈ [0, 2π]
    // A = (1/2) ∮ (x dy - y dx)
    // = (1/2) ∫₀²π (cos(t)·cos(t) - sin(t)·(-sin(t))) dt
    // = (1/2) ∫₀²π (cos²t + sin²t) dt = (1/2) · 2π = π
    auto t = SymbolicExpr::variable("t");
    auto cos_t = SymbolicExpr::cos(t);
    auto sin_t = SymbolicExpr::sin(t);

    VectorField circle_param = {cos_t, sin_t};

    auto a = SymbolicExpr::number(0);
    auto pi2 = SymbolicExpr::number(2.0 * 3.14159265358979323846);

    auto area = greens_theorem_area(circle_param, "t", a, pi2);

    EXPECT_TRUE((area != nullptr)) << "Green's area result is not null";
    if (area) {
        auto val = test_numeric_eval(area);
        ASSERT_TRUE((val.has_value())) << "area is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.14159265358979323846, 0.001) << "area of unit circle = pi";
        }
    }
}

TEST(IntegralTheorems, GreensAreaEllipse) {
    // Ellipse with semi-axes a=2, b=3
    // Area = π·a·b = 6π
    auto t = SymbolicExpr::variable("t");
    auto x_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::cos(t));
    auto y_t = SymbolicExpr::multiply(SymbolicExpr::number(3), SymbolicExpr::sin(t));

    VectorField ellipse_param = {x_t, y_t};

    auto a = SymbolicExpr::number(0);
    auto pi2 = SymbolicExpr::number(2.0 * 3.14159265358979323846);

    auto area = greens_theorem_area(ellipse_param, "t", a, pi2);

    EXPECT_TRUE((area != nullptr)) << "Green's area result is not null";
    if (area) {
        auto val = test_numeric_eval(area);
        ASSERT_TRUE((val.has_value())) << "area is numeric";
        if (val.has_value()) {
            {
                const double actual_value = (*val);
                const double expected_value = (6.0 * 3.14159265358979323846);
                const double tolerance = (0.01);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}
TEST(IntegralTheorems, GreensTheoremInvalidInput) {
    auto P = SymbolicExpr::number(0);
    auto Q = SymbolicExpr::variable("x");
    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    bool threw = false;
    try {
        // Wrong number of vars
        greens_theorem(P, Q, {"x"}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for wrong number of vars";

    threw = false;
    try {
        // Null P
        greens_theorem(nullptr, Q, {"x", "y"}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for null P";
}
