/**
 * @file test_integral_theorems_checked.cpp
 * @brief 带检查积分定理的错误、预算与支持域测试。
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

class LmcasIntegralTheoremsChecked : public ::testing::Test {
  protected:
    std::shared_ptr<SymbolicExpr> x = SymbolicExpr::variable("x");
    std::shared_ptr<SymbolicExpr> y = SymbolicExpr::variable("y");
    std::shared_ptr<SymbolicExpr> z = SymbolicExpr::variable("z");
    std::shared_ptr<SymbolicExpr> u_var = SymbolicExpr::variable("u");
    std::shared_ptr<SymbolicExpr> v_var = SymbolicExpr::variable("v");
    std::shared_ptr<SymbolicExpr> lo = SymbolicExpr::number(0);
    std::shared_ptr<SymbolicExpr> hi = SymbolicExpr::number(1);
    VectorField div_field = {x, y, z};
    VectorField stokes_field = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)};
    VectorField surface = {u_var, v_var, SymbolicExpr::number(0)};
    std::shared_ptr<SymbolicExpr> null_root;
    std::shared_ptr<SymbolicExpr> unsupported_dx = SymbolicExpr::eq(x, lo);
    std::shared_ptr<SymbolicExpr> unsupported_du = SymbolicExpr::eq(u_var, lo);
};

TEST_F(LmcasIntegralTheoremsChecked, GreenSuccess) {
    auto green_ok = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "y"}, {lo, hi}, {lo, hi});
    ASSERT_TRUE((green_ok.has_value())) << "checked Green's theorem succeeds";
    if (green_ok) {
        auto val = test_numeric_eval(green_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked Green's theorem result is numeric";
        if (val.has_value()) {
            {
                const double actual_value = *val;
                const double expected_value = 1.0;
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, 1e-10);
            }
        }
    }
}

TEST_F(LmcasIntegralTheoremsChecked, GreenAreaSuccess) {
    VectorField circle = {SymbolicExpr::cos(u_var), SymbolicExpr::sin(u_var)};
    auto area_ok = greens_theorem_area_checked(
        circle, "u", lo, SymbolicExpr::number(2.0 * 3.14159265358979323846));
    EXPECT_TRUE((area_ok.has_value())) << "checked Green's area formula succeeds";
}

TEST_F(LmcasIntegralTheoremsChecked, DivergenceSuccess) {
    auto div_ok = divergence_theorem_checked(
        div_field, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    ASSERT_TRUE((div_ok.has_value())) << "checked divergence theorem succeeds";
    if (div_ok) {
        auto val = test_numeric_eval(div_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked divergence theorem result is numeric";
        if (val.has_value()) {
            {
                const double actual_value = *val;
                const double expected_value = 3.0;
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, 1e-10);
            }
        }
    }
}

TEST_F(LmcasIntegralTheoremsChecked, StokesSuccess) {
    auto stokes_ok = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, surface, "u", "v", {lo, hi}, {lo, hi});
    ASSERT_TRUE((stokes_ok.has_value())) << "checked Stokes theorem succeeds";
    if (stokes_ok) {
        auto val = test_numeric_eval(stokes_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked Stokes theorem result is numeric";
        if (val.has_value()) {
            {
                const double actual_value = *val;
                const double expected_value = -2.0;
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, 1e-10);
            }
        }
    }
}

TEST_F(LmcasIntegralTheoremsChecked, StokesUsesCallerCoordinateNames) {
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    VectorField field{
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::multiply(a, b)};
    VectorField caller_surface{
        SymbolicExpr::number(0), u_var, v_var};

    auto result = stokes_theorem_checked(
        field, {"a", "b", "c"}, caller_surface, "u", "v",
        {lo, hi}, {lo, hi});
    ASSERT_TRUE(result.has_value()) << result.error().message;
    auto value = test_numeric_eval(result.value());
    ASSERT_TRUE(value.has_value());
    ASSERT_TRUE(std::isfinite(*value));
    EXPECT_NEAR(*value, 0.0, 1e-12);
}

TEST_F(LmcasIntegralTheoremsChecked, GreenInvalidInputs) {
    auto null_green = greens_theorem_checked(
        null_root, x, {"x", "y"}, {lo, hi}, {lo, hi});
    EXPECT_TRUE((!null_green.has_value())) << "checked Green's theorem rejects null components";
    EXPECT_TRUE((null_green.error().code == CasErrc::InvalidArgument)) << "checked Green's theorem reports InvalidArgument for null input";

    auto duplicate_vars = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "x"}, {lo, hi}, {lo, hi});
    EXPECT_TRUE((!duplicate_vars.has_value())) << "checked Green's theorem rejects duplicate variables";
    EXPECT_TRUE((duplicate_vars.error().code == CasErrc::InvalidArgument)) << "checked Green's theorem reports InvalidArgument for duplicate variables";
}

TEST_F(LmcasIntegralTheoremsChecked, DivergenceInvalidBounds) {
    auto null_bound = divergence_theorem_checked(
        div_field, {"x", "y", "z"}, {lo, nullptr}, {lo, hi}, {lo, hi});
    EXPECT_TRUE((!null_bound.has_value())) << "checked divergence theorem rejects null bounds";
    EXPECT_TRUE((null_bound.error().code == CasErrc::InvalidArgument)) << "checked divergence theorem reports InvalidArgument for null bounds";
}

TEST_F(LmcasIntegralTheoremsChecked, AreaInvalidDimension) {
    auto bad_area_dim = greens_theorem_area_checked(
        {u_var, v_var, SymbolicExpr::number(0)}, "u", lo, hi);
    EXPECT_TRUE((!bad_area_dim.has_value())) << "checked Green's area rejects non-2D parametrization";
    EXPECT_TRUE((bad_area_dim.error().code == CasErrc::InvalidArgument)) << "checked Green's area reports InvalidArgument for non-2D parametrization";
}

TEST_F(LmcasIntegralTheoremsChecked, StokesNullSurface) {
    VectorField bad_stokes_surface = {u_var, null_root, SymbolicExpr::number(0)};
    auto null_surface = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, bad_stokes_surface, "u", "v",
        {lo, hi}, {lo, hi});
    EXPECT_TRUE((!null_surface.has_value())) << "checked Stokes theorem rejects null parametrization";
    EXPECT_TRUE((null_surface.error().code == CasErrc::InvalidArgument)) << "checked Stokes theorem reports InvalidArgument for null parametrization";
}

TEST_F(LmcasIntegralTheoremsChecked, GreenCancellation) {
    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = greens_theorem_checked(
        SymbolicExpr::number(0), x, {"x", "y"}, {lo, hi}, {lo, hi},
        cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked Green's theorem observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked Green's theorem reports Cancelled";
}

TEST_F(LmcasIntegralTheoremsChecked, StokesBudget) {
    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = stokes_theorem_checked(
        stokes_field, {"x", "y", "z"}, surface, "u", "v",
        {lo, hi}, {lo, hi}, limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked Stokes theorem observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked Stokes theorem reports ResourceLimit";
}

TEST_F(LmcasIntegralTheoremsChecked, GreenUnsupportedDerivatives) {
    auto green_unsupported = greens_theorem_checked(
        SymbolicExpr::number(0), unsupported_dx, {"x", "y"},
        {lo, hi}, {lo, hi});
    EXPECT_TRUE((!green_unsupported.has_value())) << "checked Green's theorem rejects unsupported derivatives";
    EXPECT_TRUE((green_unsupported.error().code == CasErrc::Inconclusive)) << "checked Green's theorem reports Inconclusive for unsupported derivatives";

    auto area_unsupported = greens_theorem_area_checked(
        {unsupported_du, v_var}, "u", lo, hi);
    EXPECT_TRUE((!area_unsupported.has_value())) << "checked Green's area rejects unsupported path derivatives";
    EXPECT_TRUE((area_unsupported.error().code == CasErrc::Inconclusive)) << "checked Green's area reports Inconclusive for unsupported derivatives";
}

TEST_F(LmcasIntegralTheoremsChecked, SpatialUnsupportedDerivatives) {
    VectorField unsupported_div_field = {unsupported_dx, y, z};
    auto div_unsupported = divergence_theorem_checked(
        unsupported_div_field, {"x", "y", "z"}, {lo, hi}, {lo, hi}, {lo, hi});
    EXPECT_TRUE((!div_unsupported.has_value())) << "checked divergence theorem rejects unsupported derivatives";
    EXPECT_TRUE((div_unsupported.error().code == CasErrc::Inconclusive)) << "checked divergence theorem reports Inconclusive for unsupported derivatives";

    VectorField unsupported_stokes_field = {unsupported_dx, y, z};
    auto stokes_unsupported = stokes_theorem_checked(
        unsupported_stokes_field, {"x", "y", "z"}, surface, "u", "v",
        {lo, hi}, {lo, hi});
    EXPECT_TRUE((!stokes_unsupported.has_value())) << "checked Stokes theorem rejects unsupported derivatives";
    EXPECT_TRUE((stokes_unsupported.error().code == CasErrc::Inconclusive)) << "checked Stokes theorem reports Inconclusive for unsupported derivatives";
}
