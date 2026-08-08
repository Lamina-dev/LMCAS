/**
 * @file test_integral_theorems.cpp
 * @brief 积分定理单元测试：格林定理、散度定理（高斯定理）、斯托克斯定理。
 */

#include "test_common.hpp"
#include "vector_calculus.hpp"
#include "symbolic_ast.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>

using namespace lamina;


static void test_greens_theorem_constant_field()
{
    TEST_CASE("Green's theorem: P=0, Q=x => area of rectangle [0,1]x[0,1]");

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

    EXPECT_TRUE(result != nullptr, "Green's theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-10,
                "∬ 1 dA over [0,1]x[0,1] = 1");
        }
    }
}

static void test_greens_theorem_linear_field()
{
    TEST_CASE("Green's theorem: P=y, Q=x => ∬(1-1)dA = 0");

    // P = y, Q = x => ∂Q/∂x - ∂P/∂y = 1 - 1 = 0
    auto P = SymbolicExpr::variable("y");
    auto Q = SymbolicExpr::variable("x");

    auto x_lo = SymbolicExpr::number(0);
    auto x_hi = SymbolicExpr::number(2);
    auto y_lo = SymbolicExpr::number(0);
    auto y_hi = SymbolicExpr::number(3);

    auto result = greens_theorem(P, Q, {"x", "y"},
        {x_lo, x_hi}, {y_lo, y_hi});

    EXPECT_TRUE(result != nullptr, "Green's theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10,
                "∬(∂x/∂x - ∂y/∂y) dA = ∬ 0 dA = 0");
        }
    }
}

static void test_greens_theorem_quadratic()
{
    TEST_CASE("Green's theorem: P=-y^2, Q=x^2 => ∬(2x+2y)dA over [0,1]x[0,1]");

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

    EXPECT_TRUE(result != nullptr, "Green's theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 2.0, 1e-10,
                "∬(2x+2y) dA over [0,1]x[0,1] = 2");
        }
    }
}


static void test_greens_area_unit_square()
{
    TEST_CASE("Green's area: unit square parametrized as 4 segments");

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

    EXPECT_TRUE(area != nullptr, "Green's area result is not null");
    if (area) {
        auto val = test_numeric_eval(area);
        EXPECT_TRUE(val.has_value(), "area is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.14159265358979323846, 0.001,
                "area of unit circle = pi");
        }
    }
}

static void test_greens_area_ellipse()
{
    TEST_CASE("Green's area: ellipse r(t) = (2cos(t), 3sin(t))");

    // Ellipse with semi-axes a=2, b=3
    // Area = π·a·b = 6π
    auto t = SymbolicExpr::variable("t");
    auto x_t = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::cos(t));
    auto y_t = SymbolicExpr::multiply(SymbolicExpr::number(3), SymbolicExpr::sin(t));

    VectorField ellipse_param = {x_t, y_t};

    auto a = SymbolicExpr::number(0);
    auto pi2 = SymbolicExpr::number(2.0 * 3.14159265358979323846);

    auto area = greens_theorem_area(ellipse_param, "t", a, pi2);

    EXPECT_TRUE(area != nullptr, "Green's area result is not null");
    if (area) {
        auto val = test_numeric_eval(area);
        EXPECT_TRUE(val.has_value(), "area is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 6.0 * 3.14159265358979323846, 0.01,
                "area of ellipse (a=2, b=3) = 6*pi");
        }
    }
}


static void test_divergence_theorem_constant_field()
{
    TEST_CASE("Divergence theorem: F=(1,0,0) over unit cube => ∭ 0 dV = 0");

    // F = (1, 0, 0), ∇·F = 0
    // ∭ 0 dV = 0
    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)
    };

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
        {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "divergence theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10,
                "∭ div(1,0,0) dV = 0");
        }
    }
}

static void test_divergence_theorem_linear_field()
{
    TEST_CASE("Divergence theorem: F=(x,y,z) over unit cube => ∭ 3 dV = 3");

    // F = (x, y, z), ∇·F = 1 + 1 + 1 = 3
    // ∭ 3 dV over [0,1]³ = 3
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {x, y, z};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
        {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "divergence theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-10,
                "∭ div(x,y,z) dV over [0,1]^3 = 3");
        }
    }
}

static void test_divergence_theorem_quadratic_field()
{
    TEST_CASE("Divergence theorem: F=(x^2, y^2, z^2) over [0,1]^3");

    // F = (x², y², z²), ∇·F = 2x + 2y + 2z
    // ∭(2x + 2y + 2z) dV over [0,1]³
    // = 2·∫₀¹∫₀¹∫₀¹ x dz dy dx + 2·∫₀¹∫₀¹∫₀¹ y dz dy dx + 2·∫₀¹∫₀¹∫₀¹ z dz dy dx
    // = 2·(1/2) + 2·(1/2) + 2·(1/2) = 3
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)),
        SymbolicExpr::power(z, SymbolicExpr::number(2))
    };

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
        {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "divergence theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-10,
                "∭ div(x^2,y^2,z^2) dV over [0,1]^3 = 3");
        }
    }
}


static void test_stokes_theorem_constant_curl()
{
    TEST_CASE("Stokes' theorem: F=(y,-x,0), curl=(0,0,-2), flat surface z=0");

    // F = (y, -x, 0)
    // curl(F) = (0, 0, -2)
    // Surface: z=0 plane, parametrized as r(u,v) = (u, v, 0), u,v ∈ [0,1]
    // r_u × r_v = (0, 0, 1)
    // ∬ curl(F)·(r_u × r_v) du dv = ∬ (0,0,-2)·(0,0,1) du dv = ∬ -2 du dv = -2
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)
    };

    // Surface parametrization: r(u,v) = (u, v, 0)
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
        surface_param, "u", "v",
        {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "Stokes' theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, -2.0, 1e-10,
                "∬ curl(y,-x,0)·dS over [0,1]^2 = -2");
        }
    }
}

static void test_stokes_theorem_zero_curl()
{
    TEST_CASE("Stokes' theorem: F=grad(x^2+y^2+z^2), curl=0");

    // F = grad(x²+y²+z²) = (2x, 2y, 2z)
    // curl(grad(f)) = 0 for any f
    // ∬ 0·dS = 0
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y),
        SymbolicExpr::multiply(SymbolicExpr::number(2), z)
    };

    // Surface parametrization: r(u,v) = (u, v, 0)
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
        surface_param, "u", "v",
        {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "Stokes' theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10,
                "∬ curl(grad(f))·dS = 0");
        } else {
            EXPECT_TRUE(result->is_zero(),
                "∬ curl(grad(f))·dS = 0 (symbolic)");
        }
    }
}

static void test_stokes_theorem_linear_field()
{
    TEST_CASE("Stokes' theorem: F=(0,0,x), curl=(0,-1,0), flat surface z=0");

    // F = (0, 0, x)
    // curl(F) = (∂x/∂y - 0, 0 - ∂x/∂z, 0 - 0) = (0, -1, 0)
    // Wait: curl(F) = (∂F₃/∂y - ∂F₂/∂z, ∂F₁/∂z - ∂F₃/∂x, ∂F₂/∂x - ∂F₁/∂y)
    //              = (∂x/∂y - 0, 0 - ∂x/∂x, 0 - 0) = (0, -1, 0)
    // Surface: r(u,v) = (u, v, 0), u,v ∈ [0,1]
    // r_u × r_v = (0, 0, 1)
    // ∬ (0,-1,0)·(0,0,1) du dv = ∬ 0 du dv = 0
    auto x = SymbolicExpr::variable("x");

    VectorField F = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        x
    };

    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
        surface_param, "u", "v",
        {lo, hi}, {lo, hi});

    EXPECT_TRUE(result != nullptr, "Stokes' theorem result is not null");
    if (result) {
        auto val = test_numeric_eval(result);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10,
                "∬ curl(0,0,x)·(0,0,1) du dv = 0");
        } else {
            EXPECT_TRUE(result->is_zero(),
                "∬ curl(0,0,x)·(0,0,1) du dv = 0 (symbolic)");
        }
    }
}

static void test_integral_theorems_checked_contracts()
{
    TEST_CASE("Integral theorem checked APIs: explicit errors and context");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto green_ok = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "y"}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(green_ok.has_value(), "checked Green's theorem succeeds");
    if (green_ok) {
        auto val = test_numeric_eval(green_ok.value());
        EXPECT_TRUE(val.has_value(), "checked Green's theorem result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-10,
                        "checked Green's theorem integral equals 1");
        }
    }

    VectorField circle = {SymbolicExpr::cos(u_var), SymbolicExpr::sin(u_var)};
    auto area_ok = greens_theorem_area_checked(
        circle, "u", lo, SymbolicExpr::number(2.0 * 3.14159265358979323846));
    EXPECT_TRUE(area_ok.has_value(), "checked Green's area formula succeeds");

    VectorField div_field = {x, y, z};
    auto div_ok = divergence_theorem_checked(
        div_field, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(div_ok.has_value(), "checked divergence theorem succeeds");
    if (div_ok) {
        auto val = test_numeric_eval(div_ok.value());
        EXPECT_TRUE(val.has_value(), "checked divergence theorem result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-10,
                        "checked divergence theorem integral equals 3");
        }
    }

    VectorField stokes_field = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)
    };
    VectorField surface = {u_var, v_var, SymbolicExpr::number(0)};
    auto stokes_ok = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, surface, "u", "v", {lo, hi}, {lo, hi});
    EXPECT_TRUE(stokes_ok.has_value(), "checked Stokes theorem succeeds");
    if (stokes_ok) {
        auto val = test_numeric_eval(stokes_ok.value());
        EXPECT_TRUE(val.has_value(), "checked Stokes theorem result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, -2.0, 1e-10,
                        "checked Stokes theorem integral equals -2");
        }
    }

    std::shared_ptr<SymbolicExpr> null_root;
    auto null_green = greens_theorem_checked(
        null_root, x, {"x", "y"}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(!null_green.has_value(),
                "checked Green's theorem rejects null components");
    EXPECT_TRUE(null_green.error().code == CasErrc::InvalidArgument,
                "checked Green's theorem reports InvalidArgument for null input");

    auto duplicate_vars = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "x"}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(!duplicate_vars.has_value(),
                "checked Green's theorem rejects duplicate variables");
    EXPECT_TRUE(duplicate_vars.error().code == CasErrc::InvalidArgument,
                "checked Green's theorem reports InvalidArgument for duplicate variables");

    auto null_bound = divergence_theorem_checked(
        div_field, {"x", "y", "z"}, {lo, nullptr}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(!null_bound.has_value(),
                "checked divergence theorem rejects null bounds");
    EXPECT_TRUE(null_bound.error().code == CasErrc::InvalidArgument,
                "checked divergence theorem reports InvalidArgument for null bounds");

    auto bad_area_dim = greens_theorem_area_checked(
        {u_var, v_var, SymbolicExpr::number(0)}, "u", lo, hi);
    EXPECT_TRUE(!bad_area_dim.has_value(),
                "checked Green's area rejects non-2D parametrization");
    EXPECT_TRUE(bad_area_dim.error().code == CasErrc::InvalidArgument,
                "checked Green's area reports InvalidArgument for non-2D parametrization");

    VectorField bad_stokes_surface = {u_var, null_root, SymbolicExpr::number(0)};
    auto null_surface = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, bad_stokes_surface, "u", "v",
        {lo, hi}, {lo, hi});
    EXPECT_TRUE(!null_surface.has_value(),
                "checked Stokes theorem rejects null parametrization");
    EXPECT_TRUE(null_surface.error().code == CasErrc::InvalidArgument,
                "checked Stokes theorem reports InvalidArgument for null parametrization");

    lamina::CancellationToken cancellation;
    lamina::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "y"}, {lo, hi}, {lo, hi},
        cancelled_context);
    EXPECT_TRUE(!cancelled.has_value(), "checked Green's theorem observes cancellation");
    EXPECT_TRUE(cancelled.error().code == CasErrc::Cancelled,
                "checked Green's theorem reports Cancelled");

    lamina::ResourceLimits limits;
    limits.max_steps = 1;
    lamina::ComputationContext limited_context(limits);
    auto limited = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, surface, "u", "v",
        {lo, hi}, {lo, hi}, limited_context);
    EXPECT_TRUE(!limited.has_value(),
                "checked Stokes theorem observes exhausted step budget");
    EXPECT_TRUE(limited.error().code == CasErrc::ResourceLimit,
                "checked Stokes theorem reports ResourceLimit");

    auto unsupported_dx = SymbolicExpr::eq(x, lo);
    auto unsupported_du = SymbolicExpr::eq(u_var, lo);

    auto green_unsupported = greens_theorem_checked(
        SymbolicExpr::number(0), unsupported_dx, {"x", "y"},
        {lo, hi}, {lo, hi});
    EXPECT_TRUE(!green_unsupported.has_value(),
                "checked Green's theorem rejects unsupported derivatives");
    EXPECT_TRUE(green_unsupported.error().code == CasErrc::Inconclusive,
                "checked Green's theorem reports Inconclusive for unsupported derivatives");

    auto area_unsupported = greens_theorem_area_checked(
        {unsupported_du, v_var}, "u", lo, hi);
    EXPECT_TRUE(!area_unsupported.has_value(),
                "checked Green's area rejects unsupported path derivatives");
    EXPECT_TRUE(area_unsupported.error().code == CasErrc::Inconclusive,
                "checked Green's area reports Inconclusive for unsupported derivatives");

    VectorField unsupported_div_field = {unsupported_dx, y, z};
    auto div_unsupported = divergence_theorem_checked(
        unsupported_div_field, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    EXPECT_TRUE(!div_unsupported.has_value(),
                "checked divergence theorem rejects unsupported derivatives");
    EXPECT_TRUE(div_unsupported.error().code == CasErrc::Inconclusive,
                "checked divergence theorem reports Inconclusive for unsupported derivatives");

    VectorField unsupported_stokes_field = {unsupported_dx, y, z};
    auto stokes_unsupported = stokes_theorem_checked(
        unsupported_stokes_field, {"x", "y", "z"}, surface, "u", "v",
        {lo, hi}, {lo, hi});
    EXPECT_TRUE(!stokes_unsupported.has_value(),
                "checked Stokes theorem rejects unsupported derivatives");
    EXPECT_TRUE(stokes_unsupported.error().code == CasErrc::Inconclusive,
                "checked Stokes theorem reports Inconclusive for unsupported derivatives");
}


static void test_greens_theorem_invalid_input()
{
    TEST_CASE("Green's theorem: invalid inputs throw");

    auto P = SymbolicExpr::number(0);
    auto Q = SymbolicExpr::variable("x");
    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    bool threw = false;
    try {
        // Wrong number of vars
        greens_theorem(P, Q, {"x"}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for wrong number of vars");

    threw = false;
    try {
        // Null P
        greens_theorem(nullptr, Q, {"x", "y"}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for null P");
}

static void test_divergence_theorem_invalid_input()
{
    TEST_CASE("Divergence theorem: invalid inputs throw");

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    bool threw = false;
    try {
        // Wrong dimension F
        VectorField F2 = {SymbolicExpr::number(1), SymbolicExpr::number(0)};
        divergence_theorem(F2, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for 2D vector field");

    threw = false;
    try {
        // Wrong number of vars
        VectorField F3 = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};
        divergence_theorem(F3, {"x", "y"}, {lo, hi}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for wrong number of vars");
}

static void test_stokes_theorem_invalid_input()
{
    TEST_CASE("Stokes' theorem: invalid inputs throw");

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    bool threw = false;
    try {
        // Wrong dimension F
        VectorField F2 = {SymbolicExpr::number(1), SymbolicExpr::number(0)};
        stokes_theorem(F2, {"x", "y", "z"}, surface_param, "u", "v",
            {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for 2D vector field");

    threw = false;
    try {
        // Wrong dimension parametrization
        VectorField F3 = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};
        VectorField bad_param = {u_var, v_var};
        stokes_theorem(F3, {"x", "y", "z"}, bad_param, "u", "v",
            {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "throws for 2D parametrization");
}

int main()
{
    // Green's theorem tests
    test_greens_theorem_constant_field();
    test_greens_theorem_linear_field();
    test_greens_theorem_quadratic();

    // Green's theorem area formula tests
    test_greens_area_unit_square();
    test_greens_area_ellipse();

    // Divergence theorem tests
    test_divergence_theorem_constant_field();
    test_divergence_theorem_linear_field();
    test_divergence_theorem_quadratic_field();

    // Stokes' theorem tests
    test_stokes_theorem_constant_curl();
    test_stokes_theorem_zero_curl();
    test_stokes_theorem_linear_field();
    test_integral_theorems_checked_contracts();

    // Input validation tests
    test_greens_theorem_invalid_input();
    test_divergence_theorem_invalid_input();
    test_stokes_theorem_invalid_input();

    return TEST_REPORT();
}
