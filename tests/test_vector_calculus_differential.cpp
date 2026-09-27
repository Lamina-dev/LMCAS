#include "test_common.hpp"
#include "expr.hpp"
#include "vector_calculus.hpp"
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace LMCAS;

TEST(VectorCalculusDifferential, GradientBasic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::power(z, SymbolicExpr::number(2)));

    auto grad = gradient(f, {"x", "y", "z"});

    EXPECT_TRUE((grad.size() == 3)) << "gradient has 3 components";

    auto expected_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x)->simplify();
    EXPECT_TRUE(test_expression_text(grad[0], expected_x)) << "df/dx = 2x";

    auto expected_y = SymbolicExpr::multiply(SymbolicExpr::number(2), y)->simplify();
    EXPECT_TRUE(test_expression_text(grad[1], expected_y)) << "df/dy = 2y";

    auto expected_z = SymbolicExpr::multiply(SymbolicExpr::number(2), z)->simplify();
    EXPECT_TRUE(test_expression_text(grad[2], expected_z)) << "df/dz = 2z";
}

TEST(VectorCalculusDifferential, GradientMixed) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::add(
        SymbolicExpr::multiply(x, y),
        SymbolicExpr::multiply(y, z));

    auto grad = gradient(f, {"x", "y", "z"});

    EXPECT_TRUE((grad.size() == 3)) << "gradient has 3 components";

    EXPECT_TRUE(test_expression_text(grad[0], y)) << "df/dx = y";

    auto expected_y = SymbolicExpr::add(x, z)->simplify();
    EXPECT_TRUE(test_expression_text(grad[1], expected_y)) << "df/dy = x + z";

    EXPECT_TRUE(test_expression_text(grad[2], y)) << "df/dz = y";
}

TEST(VectorCalculusDifferential, DivergenceBasic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)),
        SymbolicExpr::power(z, SymbolicExpr::number(2))};

    auto div = divergence(F, {"x", "y", "z"});

    auto expected = SymbolicExpr::add(
                        SymbolicExpr::add(
                            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
                            SymbolicExpr::multiply(SymbolicExpr::number(2), y)),
                        SymbolicExpr::multiply(SymbolicExpr::number(2), z))
                        ->simplify();

    EXPECT_TRUE(test_expression_text(div, expected)) << "div(x^2, y^2, z^2) = 2x + 2y + 2z";
}

TEST(VectorCalculusDifferential, DivergenceConstantField) {
    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(2),
        SymbolicExpr::number(3)};

    auto div = divergence(F, {"x", "y", "z"});

    EXPECT_TRUE((div->is_zero())) << "div(constant field) = 0";
}

TEST(VectorCalculusDifferential, Curl3d) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)};

    auto c = curl(F, {"x", "y", "z"});

    EXPECT_TRUE((c.size() == 3)) << "curl has 3 components";

    EXPECT_TRUE((c[0]->is_zero())) << "curl_x = 0";

    EXPECT_TRUE((c[1]->is_zero())) << "curl_y = 0";

    auto val = test_numeric_eval(c[2]);
    EXPECT_TRUE((val.has_value() && std::abs(*val - (-2.0)) < 1e-10)) << "curl_z = -2";
}

TEST(VectorCalculusDifferential, Curl2d) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)};

    auto c = curl(F, {"x", "y"});

    EXPECT_TRUE((c.size() == 1)) << "2D curl returns single element";

    auto val = test_numeric_eval(c[0]);
    EXPECT_TRUE((val.has_value() && std::abs(*val - (-2.0)) < 1e-10)) << "2D scalar curl = -2";
}

TEST(VectorCalculusDifferential, CurlGradIsZero) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::multiply(SymbolicExpr::multiply(x, y), z);

    auto grad = gradient(f, {"x", "y", "z"});
    auto c = curl(grad, {"x", "y", "z"});

    EXPECT_TRUE((c.size() == 3)) << "curl(grad) has 3 components";

    for (int i = 0; i < 3; ++i) {
        auto val = test_numeric_eval(c[i]);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10) << "curl(grad(f))_" + std::to_string(i) + " = 0 (numeric)";
        } else {
            auto test_expr = c[i]->substitute("x", SymbolicExpr::number(2));
            test_expr = test_expr->substitute("y", SymbolicExpr::number(3));
            test_expr = test_expr->substitute("z", SymbolicExpr::number(5));
            test_expr = test_expr->simplify();
            EXPECT_TRUE((test_expr->is_zero())) << "curl(grad(f))_" + std::to_string(i) + " = 0 (at point)";
        }
    }
}

TEST(VectorCalculusDifferential, LaplacianBasic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::power(z, SymbolicExpr::number(2)));

    auto lap = laplacian(f, {"x", "y", "z"});

    auto val = test_numeric_eval(lap);
    ASSERT_TRUE((val.has_value())) << "laplacian is numeric";
    if (val.has_value()) {
        EXPECT_NEAR(*val, 6.0, 1e-10) << "laplacian(x^2+y^2+z^2) = 6";
    }
}

TEST(VectorCalculusDifferential, LaplacianHarmonic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto lap = laplacian(f, {"x", "y"});

    auto val = test_numeric_eval(lap);
    EXPECT_TRUE((val.has_value() && std::abs(*val) < 1e-10)) << "laplacian(x^2-y^2) = 0 (harmonic)";
}

TEST(VectorCalculusDifferential, LaplacianEqualsDivGrad) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::power(y, SymbolicExpr::number(3)));

    auto lap = laplacian(f, {"x", "y"});
    auto grad_f = gradient(f, {"x", "y"});
    auto div_grad = divergence(grad_f, {"x", "y"});

    auto lap_at = lap->substitute("x", SymbolicExpr::number(2));
    lap_at = lap_at->substitute("y", SymbolicExpr::number(3));
    lap_at = lap_at->simplify();

    auto dg_at = div_grad->substitute("x", SymbolicExpr::number(2));
    dg_at = dg_at->substitute("y", SymbolicExpr::number(3));
    dg_at = dg_at->simplify();

    auto v1 = test_numeric_eval(lap_at);
    auto v2 = test_numeric_eval(dg_at);
    EXPECT_TRUE((v1.has_value() && v2.has_value() && std::abs(*v1 - *v2) < 1e-10)) << "laplacian(f) = div(grad(f)) at (2,3)";
}

TEST(VectorCalculusDifferential, DirectionalDerivativeAxis) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    VectorField dir = {SymbolicExpr::number(1), SymbolicExpr::number(0)};

    auto dd = directional_derivative(f, {"x", "y"}, dir);
    EXPECT_TRUE((dd != nullptr)) << "directional derivative is not null";

    auto dd_at = dd->substitute("x", SymbolicExpr::number(3));
    dd_at = dd_at->substitute("y", SymbolicExpr::number(5));
    dd_at = dd_at->simplify();
    auto val = test_numeric_eval(dd_at);
    EXPECT_TRUE((val.has_value() && std::abs(*val - 6.0) < 1e-10)) << "D_(1,0) f at (3,5) = 2*3 = 6";
}

TEST(VectorCalculusDifferential, DirectionalDerivativeNormalization) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(x, y);

    VectorField dir = {SymbolicExpr::number(3), SymbolicExpr::number(4)};

    auto dd = directional_derivative(f, {"x", "y"}, dir);
    EXPECT_TRUE((dd != nullptr)) << "directional derivative is not null";

    auto val = test_numeric_eval(dd);
    if (val.has_value()) {
        EXPECT_NEAR(*val, 1.4, 1e-10) << "D_(3,4) (x+y) = 7/5 = 1.4";
    } else {
        auto dd_s = dd->simplify();
        auto val2 = test_numeric_eval(dd_s);
        EXPECT_TRUE((val2.has_value() && std::abs(*val2 - 1.4) < 1e-10)) << "D_(3,4) (x+y) = 7/5 after simplify";
    }
}

TEST(VectorCalculusDifferential, DirectionalDerivativeZeroVector) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(2));

    VectorField dir = {SymbolicExpr::number(0)};

    auto dd = directional_derivative(f, {"x"}, dir);

    EXPECT_TRUE((dd == nullptr)) << "zero direction returns nullptr";
}

TEST(VectorCalculusDifferential, DirectionalDerivativeHigherOrder) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));

    VectorField dir = {SymbolicExpr::number(1)};

    auto dd2 = directional_derivative(f, {"x"}, dir, 2);
    EXPECT_TRUE((dd2 != nullptr)) << "second order derivative is not null";

    auto dd2_at = dd2->substitute("x", SymbolicExpr::number(4));
    dd2_at = dd2_at->simplify();
    auto val = test_numeric_eval(dd2_at);
    EXPECT_TRUE((val.has_value() && std::abs(*val - 24.0) < 1e-10)) << "D^2_(1) x^3 at x=4 = 6*4 = 24";
}

TEST(VectorCalculusDifferential, DivCurlIsZero) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::multiply(x, y),
        SymbolicExpr::multiply(y, z),
        SymbolicExpr::multiply(z, x)};

    auto c = curl(F, {"x", "y", "z"});
    auto div_curl = divergence(c, {"x", "y", "z"});

    auto val = test_numeric_eval(div_curl);
    if (val.has_value()) {
        EXPECT_NEAR(*val, 0.0, 1e-10) << "div(curl(F)) = 0 (numeric)";
    } else {
        auto test_expr = div_curl->substitute("x", SymbolicExpr::number(2));
        test_expr = test_expr->substitute("y", SymbolicExpr::number(3));
        test_expr = test_expr->substitute("z", SymbolicExpr::number(5));
        test_expr = test_expr->simplify();
        auto v = test_numeric_eval(test_expr);
        EXPECT_TRUE((v.has_value() && std::abs(*v) < 1e-10)) << "div(curl(F)) = 0 at (2,3,5)";
    }
}

TEST(VectorCalculusDifferential, GradientCheckedShape) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto grad = gradient_checked(f, {"x", "y"});
    ASSERT_TRUE((grad.has_value())) << "checked gradient succeeds";
    if (grad) {
        EXPECT_TRUE((grad.value().size() == 2)) << "checked gradient returns two components";
    }
}

TEST(VectorCalculusDifferential, DifferentialNullInputs) {
    auto x = SymbolicExpr::variable("x");

    auto null_grad = gradient_checked(nullptr, {"x"});
    EXPECT_TRUE((!null_grad.has_value())) << "checked gradient rejects null expression";
    EXPECT_TRUE((null_grad.error().code == CasErrc::InvalidArgument)) << "checked gradient reports InvalidArgument";

    std::shared_ptr<SymbolicExpr> null_root;
    VectorField bad_field = {x, null_root};
    auto bad_div = divergence_checked(bad_field, {"x", "y"});
    EXPECT_TRUE((!bad_div.has_value())) << "checked divergence rejects null component";
    EXPECT_TRUE((bad_div.error().code == CasErrc::InvalidArgument)) << "checked divergence reports InvalidArgument for null component";
}

TEST(VectorCalculusDifferential, CurlCheckedDimensions) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField two_d = {x, y};
    auto bad_curl = curl_checked(two_d, {"x"});
    EXPECT_TRUE((!bad_curl.has_value())) << "checked curl rejects dimension mismatch";
    EXPECT_TRUE((bad_curl.error().code == CasErrc::InvalidArgument)) << "checked curl reports InvalidArgument for dimension mismatch";
}

TEST(VectorCalculusDifferential, LaplacianCheckedVariables) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto empty_lap = laplacian_checked(f, {});
    EXPECT_TRUE((!empty_lap.has_value())) << "checked laplacian rejects empty variables";
    EXPECT_TRUE((empty_lap.error().code == CasErrc::InvalidArgument)) << "checked laplacian reports InvalidArgument for empty variables";
}

TEST(VectorCalculusDifferential, DirectionalCheckedZero) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    VectorField zero_dir = {SymbolicExpr::number(0), SymbolicExpr::number(0)};
    auto zero_direction = directional_derivative_checked(f, {"x", "y"}, zero_dir);
    EXPECT_TRUE((!zero_direction.has_value())) << "checked directional derivative rejects zero direction";
    EXPECT_TRUE((zero_direction.error().code == CasErrc::DomainError)) << "checked directional derivative reports DomainError for zero direction";
    EXPECT_TRUE((directional_derivative(f, {"x", "y"}, zero_dir) == nullptr)) << "legacy directional derivative unwraps zero direction to nullptr";
}

TEST(VectorCalculusDifferential, GradientCheckedCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = gradient_checked(f, {"x"}, cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked gradient observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked gradient reports Cancelled";
}

TEST(VectorCalculusDifferential, LaplacianCheckedBudget) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    LMCAS::ResourceLimits limits;
    limits.max_steps = 0;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = laplacian_checked(f, {"x"}, limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked laplacian observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked laplacian reports ResourceLimit";
}

TEST(VectorCalculusDifferential, UnsupportedDerivativeErrors) {
    const auto floor_result = LMCAS::floor(SymbolicExpr::variable("x"));
    ASSERT_TRUE(floor_result.has_value());
    const auto& floor_x = floor_result.value();
    const auto zero = SymbolicExpr::number(0);
    ComputationContext context;
    auto grad = gradient_checked(floor_x, {"x"}, context);
    ASSERT_FALSE(grad.has_value());
    EXPECT_EQ(grad.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(grad.error().operation, "gradient");

    auto div = divergence_checked({floor_x}, {"x"});
    ASSERT_FALSE(div.has_value());
    EXPECT_EQ(div.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(div.error().operation, "divergence");

    auto curl_2d = curl_checked({zero, floor_x}, {"x", "y"});
    ASSERT_FALSE(curl_2d.has_value());
    EXPECT_EQ(curl_2d.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(curl_2d.error().operation, "curl");

    const auto floor_y_result = LMCAS::floor(SymbolicExpr::variable("y"));
    ASSERT_TRUE(floor_y_result.has_value());
    const auto& floor_y = floor_y_result.value();
    auto curl_3d = curl_checked({zero, zero, floor_y}, {"x", "y", "z"});
    ASSERT_FALSE(curl_3d.has_value());
    EXPECT_EQ(curl_3d.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(curl_3d.error().operation, "curl");

    auto lap = laplacian_checked(floor_x, {"x"});
    ASSERT_FALSE(lap.has_value());
    EXPECT_EQ(lap.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(lap.error().operation, "laplacian");

    auto directional = directional_derivative_checked(
        floor_x, {"x"}, {SymbolicExpr::number(1)}, 1);
    ASSERT_FALSE(directional.has_value());
    EXPECT_EQ(directional.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(directional.error().operation, "directional_derivative");
    EXPECT_THROW(gradient(floor_x, {"x"}), std::runtime_error);
}

TEST(VectorCalculusDifferential, UnsupportedSecondDerivative) {
    const auto floor_result = LMCAS::floor(SymbolicExpr::variable("x"));
    ASSERT_TRUE(floor_result.has_value());
    const auto& floor_x = floor_result.value();
    const auto integral = SymbolicExpr::make_integral(floor_x, "x");
    auto first = gradient_checked(integral, {"x"});
    ASSERT_TRUE(first.has_value());
    ASSERT_EQ(first.value().size(), 1u);
    EXPECT_TRUE(test_same_expression(first.value()[0], floor_x));

    auto lap = laplacian_checked(integral, {"x"});
    ASSERT_FALSE(lap.has_value());
    EXPECT_EQ(lap.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(lap.error().operation, "laplacian");

    auto directional = directional_derivative_checked(
        integral, {"x"}, {SymbolicExpr::number(1)}, 2);
    ASSERT_FALSE(directional.has_value());
    EXPECT_EQ(directional.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(directional.error().operation, "directional_derivative");
    EXPECT_THROW(laplacian(integral, {"x"}), std::runtime_error);
}
