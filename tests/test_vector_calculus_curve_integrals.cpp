#include "test_common.hpp"
#include "vector_calculus.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

TEST(VectorCalculusCurveIntegrals, CurveIntegralScalarLine) {
    /**
     * @brief 线段 r(t) = (t, 2t)，t 属于 [0, 1]。
     * 弧长元 |r'(t)| = sqrt(1 + 4) = sqrt(5)，
     * 故 integral_0^1 1 * sqrt(5) dt = sqrt(5)。
     */
    auto t_var = SymbolicExpr::variable("t");
    VectorField param = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    auto f = SymbolicExpr::number(1);
    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_scalar(f, param, "t", a, b);
    EXPECT_TRUE((result != nullptr)) << "curve_integral_scalar result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        double expected = std::sqrt(5.0);
        ASSERT_TRUE((val.has_value())) << "arc length result is numeric";
        if (val.has_value()) {
            {
                const double actual_value = (*val);
                const double expected_value = (expected);
                const double tolerance = (1e-6);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralScalarCircle) {
    /**
     * @brief 单位圆 r(t) = (cos(t), sin(t))，t 属于 [0, 2*pi]。
     * 弧长元 |r'(t)| = sqrt(sin^2(t) + cos^2(t)) = 1，
     * 故 integral_0^{2pi} 1 * 1 dt = 2*pi。
     */
    auto t_var = SymbolicExpr::variable("t");
    VectorField param = {
        SymbolicExpr::cos(t_var),
        SymbolicExpr::sin(t_var)};

    auto f = SymbolicExpr::number(1);
    auto a = SymbolicExpr::number(0);
    auto pi2 = SymbolicExpr::number(2.0 * 3.14159265358979323846);

    auto result = curve_integral_scalar(f, param, "t", a, pi2);
    EXPECT_TRUE((result != nullptr)) << "circle arc length result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        double expected = 2.0 * 3.14159265358979323846;
        ASSERT_TRUE((val.has_value())) << "circle arc length is numeric";
        if (val.has_value()) {
            {
                const double actual_value = (*val);
                const double expected_value = (expected);
                const double tolerance = (1e-3);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralVectorConservative) {
    /**
     * @brief 梯度场 F = (2x, 2y) = grad(x^2 + y^2) 沿 r(t) = (t, t) 积分。
     * t 属于 [0, 1]，F(r(t)) dot r'(t) = (2t, 2t) dot (1, 1) = 4t，
     * 故 integral_0^1 4t dt = 2。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)};

    VectorField param = {t_var, t_var};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE((result != nullptr)) << "curve_integral_vector result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "conservative field integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 2.0, 1e-6) << "work of grad(x^2+y^2) along (t,t) = 2";
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralVectorWork) {
    /**
     * @brief 向量场 F = (y, -x) 沿 r(t) = (t, 0) 积分，t 属于 [0, 1]。
     * F(r(t)) = (0, -t)，r'(t) = (1, 0)，
     * 故 integral_0^1 (0, -t) dot (1, 0) dt = integral_0^1 0 dt = 0。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)};

    VectorField param = {t_var, SymbolicExpr::number(0)};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE((result != nullptr)) << "work integral result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "work integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-6) << "work of (y,-x) along x-axis = 0";
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralVector3d) {
    /**
     * @brief 常向量场 F = (1, 0, 0) 沿 r(t) = (t, 0, 0) 积分，t 属于 [0, 3]。
     * 积分为 integral_0^3 (1,0,0) dot (1,0,0) dt = 3。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)};

    VectorField param = {t_var, SymbolicExpr::number(0), SymbolicExpr::number(0)};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(3);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE((result != nullptr)) << "3D curve integral result is not null";

    if (result) {
        auto val = test_numeric_eval(result);
        ASSERT_TRUE((val.has_value())) << "3D curve integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-6) << "work of (1,0,0) along x-axis [0,3] = 3";
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralNumericFallbackFailureIsInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto w = SymbolicExpr::variable("w");
    auto t = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto unsupported = SymbolicExpr::sin(SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), w));
    VectorField line = {t, SymbolicExpr::number(0)};

    auto legacy = curve_integral_scalar(
        unsupported, line, "t", zero, one);
    EXPECT_TRUE((legacy == nullptr)) << "legacy curve integral returns null when fallback cannot evaluate samples";

    auto checked = curve_integral_scalar_checked(
        unsupported, line, "t", zero, one);
    EXPECT_TRUE((!checked.has_value())) << "checked curve integral rejects unsupported numeric fallback samples";
    if (!checked.has_value()) {
        EXPECT_TRUE((checked.error().code == CasErrc::Inconclusive)) << "checked curve integral reports Inconclusive for unsupported samples";
    }
}

TEST(VectorCalculusCurveIntegrals, CurveIntegralCheckedRejectsImplicitNumericFallback) {
    auto x = SymbolicExpr::variable("x");
    auto t = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto fresnel_like = SymbolicExpr::sin(
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    VectorField line = {t, SymbolicExpr::number(0)};

    auto legacy = curve_integral_scalar(fresnel_like, line, "t", zero, one);
    EXPECT_TRUE((legacy != nullptr)) << "legacy curve integral may use numeric fallback for non-elementary integrals";

    auto checked = curve_integral_scalar_checked(fresnel_like, line, "t", zero, one);
    EXPECT_TRUE((!checked.has_value())) << "checked curve integral rejects implicit numeric fallback";
    if (!checked.has_value()) {
        EXPECT_TRUE((checked.error().code == CasErrc::Inconclusive)) << "checked curve integral reports Inconclusive when exact integral is unsupported";
    }
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedValue) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField line = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    auto scalar_ok = curve_integral_scalar_checked(
        SymbolicExpr::number(1), line, "t", zero, one);
    ASSERT_TRUE((scalar_ok.has_value())) << "checked scalar curve integral succeeds";
    if (scalar_ok) {
        auto val = test_numeric_eval(scalar_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked scalar curve integral is numeric";
        if (val.has_value()) {
            {
                const double actual_value = (*val);
                const double expected_value = (std::sqrt(5.0));
                const double tolerance = (1e-6);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveVectorCheckedValue) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField field = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)};
    VectorField diagonal = {t_var, t_var};

    auto vector_ok = curve_integral_vector_checked(field, diagonal, "t", zero, one);
    ASSERT_TRUE((vector_ok.has_value())) << "checked vector curve integral succeeds";
    if (vector_ok) {
        auto val = test_numeric_eval(vector_ok.value());
        ASSERT_TRUE((val.has_value())) << "checked vector curve integral is numeric";
        if (val.has_value()) {
            EXPECT_NEAR(*val, 2.0, 1e-6) << "checked vector curve integral equals 2";
        }
    }
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedNullField) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField line = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    auto null_scalar = curve_integral_scalar_checked(nullptr, line, "t", zero, one);
    EXPECT_TRUE((!null_scalar.has_value())) << "checked scalar curve integral rejects null scalar field";
    EXPECT_TRUE((null_scalar.error().code == CasErrc::InvalidArgument)) << "checked scalar curve integral reports InvalidArgument for null scalar";
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedNullPath) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    std::shared_ptr<SymbolicExpr> null_root;
    VectorField bad_param = {t_var, null_root};
    auto null_param = curve_integral_scalar_checked(
        SymbolicExpr::number(1), bad_param, "t", zero, one);
    EXPECT_TRUE((!null_param.has_value())) << "checked scalar curve integral rejects null parametrization";
    EXPECT_TRUE((null_param.error().code == CasErrc::InvalidArgument)) << "checked scalar curve integral reports InvalidArgument for null parametrization";
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedDimension) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField one_dim = {t_var};
    auto bad_dim = curve_integral_scalar_checked(
        SymbolicExpr::number(1), one_dim, "t", zero, one);
    EXPECT_TRUE((!bad_dim.has_value())) << "checked scalar curve integral rejects unsupported dimension";
    EXPECT_TRUE((bad_dim.error().code == CasErrc::InvalidArgument)) << "checked scalar curve integral reports InvalidArgument for unsupported dimension";
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedParameter) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField line = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    auto empty_param_name = curve_integral_scalar_checked(
        SymbolicExpr::number(1), line, "", zero, one);
    EXPECT_TRUE((!empty_param_name.has_value())) << "checked scalar curve integral rejects empty parameter name";
    EXPECT_TRUE((empty_param_name.error().code == CasErrc::InvalidArgument)) << "checked scalar curve integral reports InvalidArgument for empty parameter";
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedBound) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField line = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    auto null_bound = curve_integral_scalar_checked(
        SymbolicExpr::number(1), line, "t", nullptr, one);
    EXPECT_TRUE((!null_bound.has_value())) << "checked scalar curve integral rejects null bounds";
    EXPECT_TRUE((null_bound.error().code == CasErrc::InvalidArgument)) << "checked scalar curve integral reports InvalidArgument for null bounds";
}

TEST(VectorCalculusCurveIntegrals, CurveVectorCheckedNullField) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    std::shared_ptr<SymbolicExpr> null_root;
    VectorField diagonal = {t_var, t_var};

    VectorField bad_field = {x, null_root};
    auto null_field_component = curve_integral_vector_checked(
        bad_field, diagonal, "t", zero, one);
    EXPECT_TRUE((!null_field_component.has_value())) << "checked vector curve integral rejects null field components";
    EXPECT_TRUE((null_field_component.error().code == CasErrc::InvalidArgument)) << "checked vector curve integral reports InvalidArgument for null field";
}

TEST(VectorCalculusCurveIntegrals, CurveVectorCheckedDimension) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField field = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)};
    VectorField diagonal = {t_var, t_var};

    auto mismatch = curve_integral_vector_checked(field, {t_var, t_var, t_var},
                                                  "t", zero, one);
    EXPECT_TRUE((!mismatch.has_value())) << "checked vector curve integral rejects dimension mismatch";
    EXPECT_TRUE((mismatch.error().code == CasErrc::InvalidArgument)) << "checked vector curve integral reports InvalidArgument for dimension mismatch";
}

TEST(VectorCalculusCurveIntegrals, CurveScalarCheckedCancellation) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    VectorField line = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)};

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = curve_integral_scalar_checked(
        SymbolicExpr::number(1), line, "t", zero, one, cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked scalar curve integral observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked scalar curve integral reports Cancelled";
}

TEST(VectorCalculusCurveIntegrals, CurveVectorCheckedBudget) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField field = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)};
    VectorField diagonal = {t_var, t_var};

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = curve_integral_vector_checked(field, diagonal, "t", zero, one,
                                                 limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked vector curve integral observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked vector curve integral reports ResourceLimit";
}

TEST(VectorCalculusCurveIntegrals, CurveCheckedUnsupportedDerivatives) {
    auto t_var = SymbolicExpr::variable("t");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField field = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)};
    VectorField diagonal = {t_var, t_var};

    auto unsupported_derivative = SymbolicExpr::eq(t_var, zero);
    VectorField unsupported_path = {unsupported_derivative, t_var};
    auto unsupported_scalar = curve_integral_scalar_checked(
        SymbolicExpr::number(1), unsupported_path, "t", zero, one);
    EXPECT_TRUE((!unsupported_scalar.has_value())) << "checked scalar curve integral rejects unsupported path derivatives";
    EXPECT_TRUE((unsupported_scalar.error().code == CasErrc::Inconclusive)) << "checked scalar curve integral reports Inconclusive for unsupported derivatives";

    auto unsupported_vector = curve_integral_vector_checked(
        field, unsupported_path, "t", zero, one);
    EXPECT_TRUE((!unsupported_vector.has_value())) << "checked vector curve integral rejects unsupported path derivatives";
    EXPECT_TRUE((unsupported_vector.error().code == CasErrc::Inconclusive)) << "checked vector curve integral reports Inconclusive for unsupported derivatives";
}
