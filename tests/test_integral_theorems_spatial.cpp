/**
 * @file test_integral_theorems_spatial.cpp
 * @brief 散度定理与 Stokes 定理契约。
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

TEST(IntegralTheoremsSpatial, DivergenceTheoremConstantField) {
    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
                                     {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "divergence theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10) << "∭ div(1,0,0) dV = 0";
        }
    }
}

TEST(IntegralTheoremsSpatial, DivergenceTheoremLinearField) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {x, y, z};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
                                     {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "divergence theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-10) << "∭ div(x,y,z) dV over [0,1]^3 = 3";
        }
    }
}

TEST(IntegralTheoremsSpatial, DivergenceTheoremQuadraticField) {
    /**
     * @brief 向量场 F = (x², y², z²) 的散度为 2x + 2y + 2z。
     * 在 [0,1]³ 上，∭(2x + 2y + 2z) dV = 2·(1/2) + 2·(1/2) + 2·(1/2) = 3。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)),
        SymbolicExpr::power(z, SymbolicExpr::number(2))};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = divergence_theorem(F, {"x", "y", "z"},
                                     {lo, hi}, {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "divergence theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-10) << "∭ div(x^2,y^2,z^2) dV over [0,1]^3 = 3";
        }
    }
}

TEST(IntegralTheoremsSpatial, StokesTheoremConstantCurl) {
    /**
     * @brief 向量场 F = (y, -x, 0) 的旋度为 (0, 0, -2)。
     * z=0 曲面取 r(u,v) = (u, v, 0)，u,v ∈ [0,1]，r_u × r_v = (0, 0, 1)。
     * 曲面积分 ∬ curl(F)·(r_u × r_v) du dv = ∬ -2 du dv = -2。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)};

    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
                                 surface_param, "u", "v",
                                 {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "Stokes' theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "result is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, -2.0, 1e-10) << "∬ curl(y,-x,0)·dS over [0,1]^2 = -2";
        }
    }
}

TEST(IntegralTheoremsSpatial, StokesTheoremUsesCallerCoordinateNames) {
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto u = SymbolicExpr::variable("u");
    auto v = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    VectorField field{
        zero, zero, SymbolicExpr::multiply(a, b)};
    VectorField surface{zero, u, v};

    auto result = stokes_theorem(
        field, {"a", "b", "c"}, surface, "u", "v",
        {zero, one}, {zero, one});
    ASSERT_NE(result, nullptr);
    auto value = test_numeric_eval(result);
    ASSERT_TRUE(value.has_value());
    ASSERT_TRUE(std::isfinite(*value));
    EXPECT_NEAR(*value, 0.0, 1e-12);
}

TEST(IntegralTheoremsSpatial, StokesTheoremZeroCurl) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y),
        SymbolicExpr::multiply(SymbolicExpr::number(2), z)};

    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
                                 surface_param, "u", "v",
                                 {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "Stokes' theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10) << "∬ curl(grad(f))·dS = 0";
        } else {
            EXPECT_TRUE((result->is_zero())) << "∬ curl(grad(f))·dS = 0 (symbolic)";
        }
    }
}

TEST(IntegralTheoremsSpatial, StokesTheoremLinearField) {
    /**
     * @brief 向量场 F = (0, 0, x) 的旋度为 (∂x/∂y, -∂x/∂x, 0) = (0, -1, 0)。
     * 曲面取 r(u,v) = (u, v, 0)，u,v ∈ [0,1]，r_u × r_v = (0, 0, 1)。
     * 曲面积分 ∬ (0,-1,0)·(0,0,1) du dv = 0。
     */
    auto x = SymbolicExpr::variable("x");

    VectorField F = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        x};

    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    auto result = stokes_theorem(F, {"x", "y", "z"},
                                 surface_param, "u", "v",
                                 {lo, hi}, {lo, hi});

    EXPECT_TRUE((result != nullptr)) << "Stokes' theorem result is not null";
    if (result) {
        auto val = test_numeric_eval(result);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10) << "∬ curl(0,0,x)·(0,0,1) du dv = 0";
        } else {
            EXPECT_TRUE((result->is_zero())) << "∬ curl(0,0,x)·(0,0,1) du dv = 0 (symbolic)";
        }
    }
}
TEST(IntegralTheoremsSpatial, DivergenceTheoremInvalidInput) {
    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);

    bool threw = false;
    try {
        VectorField F2 = {SymbolicExpr::number(1), SymbolicExpr::number(0)};
        divergence_theorem(F2, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for 2D vector field";

    threw = false;
    try {
        VectorField F3 = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};
        divergence_theorem(F3, {"x", "y"}, {lo, hi}, {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for wrong number of vars";
}

TEST(IntegralTheoremsSpatial, StokesTheoremInvalidInput) {
    auto lo = SymbolicExpr::number(0);
    auto hi = SymbolicExpr::number(1);
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    VectorField surface_param = {u_var, v_var, SymbolicExpr::number(0)};

    bool threw = false;
    try {
        VectorField F2 = {SymbolicExpr::number(1), SymbolicExpr::number(0)};
        stokes_theorem(F2, {"x", "y", "z"}, surface_param, "u", "v",
                       {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for 2D vector field";

    threw = false;
    try {
        VectorField F3 = {SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(0)};
        VectorField bad_param = {u_var, v_var};
        stokes_theorem(F3, {"x", "y", "z"}, bad_param, "u", "v",
                       {lo, hi}, {lo, hi});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "throws for 2D parametrization";
}
