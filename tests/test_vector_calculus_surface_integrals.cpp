#include "test_common.hpp"
#include "vector_calculus.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

TEST(VectorCalculusSurfaceIntegrals, SurfaceIntegralScalarPlane) {
    /**
     * @brief 单位平面方块的面积积分为 1。
     * r(u,v) = (u, v, 0)，u,v ∈ [0,1]；
     * r_u = (1, 0, 0)，r_v = (0, 1, 0)；
     * r_u × r_v = (0, 0, 1)，|r_u × r_v| = 1，∬ 1 * 1 du dv = 1。
     */
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto f = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_scalar(f, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((result != nullptr)) << "surface_integral_scalar result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "surface area is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "area of unit square = 1";
        }
    }
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceIntegralVectorFlux) {
    /**
     * @brief F = (0, 0, 1) 穿过单位平面方块的通量为 1。
     * r(u,v) = (u, v, 0)，u,v ∈ [0,1]，r_u × r_v = (0, 0, 1)；
     * F · (r_u × r_v) = 1，∬ 1 du dv = 1。
     */
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField F = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1)};

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_vector(F, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((result != nullptr)) << "surface_integral_vector result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "flux is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "flux of (0,0,1) through unit square = 1";
        }
    }
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceIntegralVectorZeroFlux) {
    /**
     * @brief 切向场 F = (1, 0, 0) 穿过平面方块的通量为 0。
     * r(u,v) = (u, v, 0)，r_u × r_v = (0, 0, 1)；
     * F · (r_u × r_v) = 0，∬ 0 du dv = 0。
     */
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)};

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_vector(F, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((result != nullptr)) << "zero flux result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "zero flux is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-6) << "flux of tangent field through surface = 0";
        }
    }
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedValue) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto scalar_ok = surface_integral_scalar_checked(
        SymbolicExpr::number(1), plane, "u", "v", zero, one, zero, one);
    ASSERT_TRUE((scalar_ok.has_value())) << "checked scalar surface integral succeeds";
    if (scalar_ok) {
        auto val = test_numeric_eval(scalar_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked scalar surface integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "checked scalar surface integral equals 1";
        }
    }
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceVectorCheckedValue) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    VectorField flux_field = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1)};

    auto vector_ok = surface_integral_vector_checked(
        flux_field, plane, "u", "v", zero, one, zero, one);
    ASSERT_TRUE((vector_ok.has_value())) << "checked vector surface integral succeeds";
    if (vector_ok) {
        auto val = test_numeric_eval(vector_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked vector surface integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "checked vector surface integral equals 1";
        }
    }
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedNullField) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto null_scalar = surface_integral_scalar_checked(
        nullptr, plane, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!null_scalar.has_value())) << "checked scalar surface integral rejects null scalar field";
    EXPECT_TRUE((null_scalar.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for null scalar";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedNullPath) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    std::shared_ptr<SymbolicExpr> null_root;
    VectorField bad_param = {u_var, null_root, SymbolicExpr::number(0)};
    auto null_param = surface_integral_scalar_checked(
        SymbolicExpr::number(1), bad_param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!null_param.has_value())) << "checked scalar surface integral rejects null parametrization";
    EXPECT_TRUE((null_param.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for null parametrization";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedDimension) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField bad_dim = {x, y};
    auto dim_error = surface_integral_scalar_checked(
        SymbolicExpr::number(1), bad_dim, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!dim_error.has_value())) << "checked scalar surface integral rejects non-3D parametrization";
    EXPECT_TRUE((dim_error.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for non-3D parametrization";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedParameterNames) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto empty_parameter = surface_integral_scalar_checked(
        SymbolicExpr::number(1), plane, "", "v", zero, one, zero, one);
    EXPECT_TRUE((!empty_parameter.has_value())) << "checked scalar surface integral rejects empty parameter name";
    EXPECT_TRUE((empty_parameter.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for empty parameter";

    auto same_parameter = surface_integral_scalar_checked(
        SymbolicExpr::number(1), plane, "u", "u", zero, one, zero, one);
    EXPECT_TRUE((!same_parameter.has_value())) << "checked scalar surface integral rejects duplicate parameter names";
    EXPECT_TRUE((same_parameter.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for duplicate parameters";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedBound) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto null_bound = surface_integral_scalar_checked(
        SymbolicExpr::number(1), plane, "u", "v", zero, nullptr, zero, one);
    EXPECT_TRUE((!null_bound.has_value())) << "checked scalar surface integral rejects null bounds";
    EXPECT_TRUE((null_bound.error().code == CasErrc::InvalidArgument)) << "checked scalar surface integral reports InvalidArgument for null bounds";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceVectorCheckedNullField) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    std::shared_ptr<SymbolicExpr> null_root;

    VectorField bad_field = {SymbolicExpr::number(0), null_root, SymbolicExpr::number(1)};
    auto null_field_component = surface_integral_vector_checked(
        bad_field, plane, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!null_field_component.has_value())) << "checked vector surface integral rejects null field components";
    EXPECT_TRUE((null_field_component.error().code == CasErrc::InvalidArgument)) << "checked vector surface integral reports InvalidArgument for null field";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceVectorCheckedDimension) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    auto field_dim_error = surface_integral_vector_checked(
        {SymbolicExpr::number(0), SymbolicExpr::number(1)},
        plane, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!field_dim_error.has_value())) << "checked vector surface integral rejects non-3D vector field";
    EXPECT_TRUE((field_dim_error.error().code == CasErrc::InvalidArgument)) << "checked vector surface integral reports InvalidArgument for non-3D field";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceScalarCheckedCancellation) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = surface_integral_scalar_checked(
        SymbolicExpr::number(1), plane, "u", "v", zero, one, zero, one,
        cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked scalar surface integral observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked scalar surface integral reports Cancelled";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceVectorCheckedBudget) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    VectorField flux_field = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1)};

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = surface_integral_vector_checked(
        flux_field, plane, "u", "v", zero, one, zero, one, limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked vector surface integral observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked vector surface integral reports ResourceLimit";
}

TEST(VectorCalculusSurfaceIntegrals, SurfaceCheckedUnsupportedDerivatives) {
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField plane = {u_var, v_var, SymbolicExpr::number(0)};

    VectorField flux_field = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1)};

    auto unsupported_derivative = SymbolicExpr::eq(u_var, zero);
    VectorField unsupported_surface = {unsupported_derivative, v_var, zero};
    auto unsupported_scalar = surface_integral_scalar_checked(
        SymbolicExpr::number(1), unsupported_surface, "u", "v",
        zero, one, zero, one);
    EXPECT_TRUE((!unsupported_scalar.has_value())) << "checked scalar surface integral rejects unsupported surface derivatives";
    EXPECT_TRUE((unsupported_scalar.error().code == CasErrc::Inconclusive)) << "checked scalar surface integral reports Inconclusive for unsupported derivatives";

    auto unsupported_vector = surface_integral_vector_checked(
        flux_field, unsupported_surface, "u", "v", zero, one, zero, one);
    EXPECT_TRUE((!unsupported_vector.has_value())) << "checked vector surface integral rejects unsupported surface derivatives";
    EXPECT_TRUE((unsupported_vector.error().code == CasErrc::Inconclusive)) << "checked vector surface integral reports Inconclusive for unsupported derivatives";
}
