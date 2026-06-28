/**
 * @file test_vector_calculus.cpp
 * @brief 向量微积分模块单元测试：梯度、散度、旋度、拉普拉斯算子、方向导数。
 */

#include "test_common.hpp"
#include "vector_calculus.hpp"
#include "symbolic_ast.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <memory>

using namespace lamina;

static void test_gradient_basic()
{
    TEST_CASE("gradient: f = x^2 + y^2 + z^2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    // f = x^2 + y^2 + z^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::power(z, SymbolicExpr::number(2)));

    auto grad = gradient(f, {"x", "y", "z"});

    EXPECT_TRUE(grad.size() == 3, "gradient has 3 components");

    // ∂f/∂x = 2x
    auto expected_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x)->simplify();
    EXPECT_EQ_EXPR(grad[0], expected_x, "df/dx = 2x");

    // ∂f/∂y = 2y
    auto expected_y = SymbolicExpr::multiply(SymbolicExpr::number(2), y)->simplify();
    EXPECT_EQ_EXPR(grad[1], expected_y, "df/dy = 2y");

    // ∂f/∂z = 2z
    auto expected_z = SymbolicExpr::multiply(SymbolicExpr::number(2), z)->simplify();
    EXPECT_EQ_EXPR(grad[2], expected_z, "df/dz = 2z");
}

static void test_gradient_mixed()
{
    TEST_CASE("gradient: f = x*y + y*z");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    // f = x*y + y*z
    auto f = SymbolicExpr::add(
        SymbolicExpr::multiply(x, y),
        SymbolicExpr::multiply(y, z));

    auto grad = gradient(f, {"x", "y", "z"});

    EXPECT_TRUE(grad.size() == 3, "gradient has 3 components");

    // ∂f/∂x = y
    EXPECT_EQ_EXPR(grad[0], y, "df/dx = y");

    // ∂f/∂y = x + z
    auto expected_y = SymbolicExpr::add(x, z)->simplify();
    EXPECT_EQ_EXPR(grad[1], expected_y, "df/dy = x + z");

    // ∂f/∂z = y
    EXPECT_EQ_EXPR(grad[2], y, "df/dz = y");
}

static void test_divergence_basic()
{
    TEST_CASE("divergence: F = (x^2, y^2, z^2)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)),
        SymbolicExpr::power(z, SymbolicExpr::number(2))
    };

    auto div = divergence(F, {"x", "y", "z"});

    // ∂(x²)/∂x + ∂(y²)/∂y + ∂(z²)/∂z = 2x + 2y + 2z
    auto expected = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::multiply(SymbolicExpr::number(2), y)),
        SymbolicExpr::multiply(SymbolicExpr::number(2), z))->simplify();

    EXPECT_EQ_EXPR(div, expected, "div(x^2, y^2, z^2) = 2x + 2y + 2z");
}

static void test_divergence_constant_field()
{
    TEST_CASE("divergence: constant field F = (1, 2, 3)");

    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(2),
        SymbolicExpr::number(3)
    };

    auto div = divergence(F, {"x", "y", "z"});

    // Divergence of constant field is 0
    EXPECT_TRUE(div->is_zero(), "div(constant field) = 0");
}

static void test_curl_3d()
{
    TEST_CASE("curl 3D: F = (y, -x, 0)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // F = (y, -x, 0)
    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
        SymbolicExpr::number(0)
    };

    auto c = curl(F, {"x", "y", "z"});

    EXPECT_TRUE(c.size() == 3, "curl has 3 components");

    // curl_x = ∂(0)/∂y - ∂(-x)/∂z = 0 - 0 = 0
    EXPECT_TRUE(c[0]->is_zero(), "curl_x = 0");

    // curl_y = ∂(y)/∂z - ∂(0)/∂x = 0 - 0 = 0
    EXPECT_TRUE(c[1]->is_zero(), "curl_y = 0");

    // curl_z = ∂(-x)/∂x - ∂(y)/∂y = -1 - 1 = -2
    // Check numerically
    auto val = test_numeric_eval(c[2]);
    EXPECT_TRUE(val.has_value() && std::abs(*val - (-2.0)) < 1e-10, "curl_z = -2");
}

static void test_curl_2d()
{
    TEST_CASE("curl 2D: F = (y, -x)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    };

    auto c = curl(F, {"x", "y"});

    EXPECT_TRUE(c.size() == 1, "2D curl returns single element");

    // scalar curl = ∂(-x)/∂x - ∂(y)/∂y = -1 - 1 = -2
    auto val = test_numeric_eval(c[0]);
    EXPECT_TRUE(val.has_value() && std::abs(*val - (-2.0)) < 1e-10, "2D scalar curl = -2");
}

static void test_curl_grad_is_zero()
{
    TEST_CASE("curl(grad(f)) = 0 for f = x*y*z");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    // f = x*y*z
    auto f = SymbolicExpr::multiply(SymbolicExpr::multiply(x, y), z);

    auto grad = gradient(f, {"x", "y", "z"});
    auto c = curl(grad, {"x", "y", "z"});

    EXPECT_TRUE(c.size() == 3, "curl(grad) has 3 components");

    // Each component should simplify to zero
    for (int i = 0; i < 3; ++i) {
        auto val = test_numeric_eval(c[i]);
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-10,
                "curl(grad(f))_" + std::to_string(i) + " = 0 (numeric)");
        } else {
            // Contains variables — try substituting specific values
            auto test_expr = c[i]->substitute("x", SymbolicExpr::number(2));
            test_expr = test_expr->substitute("y", SymbolicExpr::number(3));
            test_expr = test_expr->substitute("z", SymbolicExpr::number(5));
            test_expr = test_expr->simplify();
            EXPECT_TRUE(test_expr->is_zero(),
                "curl(grad(f))_" + std::to_string(i) + " = 0 (at point)");
        }
    }
}

static void test_laplacian_basic()
{
    TEST_CASE("laplacian: f = x^2 + y^2 + z^2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    // f = x^2 + y^2 + z^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::power(z, SymbolicExpr::number(2)));

    auto lap = laplacian(f, {"x", "y", "z"});

    // ∂²(x²)/∂x² + ∂²(y²)/∂y² + ∂²(z²)/∂z² = 2 + 2 + 2 = 6
    auto val = test_numeric_eval(lap);
    EXPECT_TRUE(val.has_value(), "laplacian is numeric");
    if (val.has_value()) {
        EXPECT_NEAR(*val, 6.0, 1e-10, "laplacian(x^2+y^2+z^2) = 6");
    }
}

static void test_laplacian_harmonic()
{
    TEST_CASE("laplacian: harmonic function f = x^2 - y^2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 - y^2 (harmonic in 2D)
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto lap = laplacian(f, {"x", "y"});

    // ∂²(x²)/∂x² + ∂²(-y²)/∂y² = 2 + (-2) = 0
    auto val = test_numeric_eval(lap);
    EXPECT_TRUE(val.has_value() && std::abs(*val) < 1e-10,
        "laplacian(x^2-y^2) = 0 (harmonic)");
}

static void test_laplacian_equals_div_grad()
{
    TEST_CASE("laplacian = div(grad(f)) for f = x^3 + y^3");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^3 + y^3
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::power(y, SymbolicExpr::number(3)));

    auto lap = laplacian(f, {"x", "y"});
    auto grad_f = gradient(f, {"x", "y"});
    auto div_grad = divergence(grad_f, {"x", "y"});

    // Both should give 6x + 6y — compare at a specific point
    auto lap_at = lap->substitute("x", SymbolicExpr::number(2));
    lap_at = lap_at->substitute("y", SymbolicExpr::number(3));
    lap_at = lap_at->simplify();

    auto dg_at = div_grad->substitute("x", SymbolicExpr::number(2));
    dg_at = dg_at->substitute("y", SymbolicExpr::number(3));
    dg_at = dg_at->simplify();

    auto v1 = test_numeric_eval(lap_at);
    auto v2 = test_numeric_eval(dg_at);
    EXPECT_TRUE(v1.has_value() && v2.has_value() && std::abs(*v1 - *v2) < 1e-10,
        "laplacian(f) = div(grad(f)) at (2,3)");
}

static void test_directional_derivative_axis()
{
    TEST_CASE("directional derivative: f = x^2 + y^2, dir = (1, 0)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 + y^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    VectorField dir = {SymbolicExpr::number(1), SymbolicExpr::number(0)};

    auto dd = directional_derivative(f, {"x", "y"}, dir);
    EXPECT_TRUE(dd != nullptr, "directional derivative is not null");

    // Direction (1,0) is unit. D_u f = ∂f/∂x * 1 + ∂f/∂y * 0 = 2x
    // Evaluate at x=3, y=5: should be 6
    auto dd_at = dd->substitute("x", SymbolicExpr::number(3));
    dd_at = dd_at->substitute("y", SymbolicExpr::number(5));
    dd_at = dd_at->simplify();
    auto val = test_numeric_eval(dd_at);
    EXPECT_TRUE(val.has_value() && std::abs(*val - 6.0) < 1e-10,
        "D_(1,0) f at (3,5) = 2*3 = 6");
}

static void test_directional_derivative_normalization()
{
    TEST_CASE("directional derivative: f = x + y, dir = (3, 4)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x + y
    auto f = SymbolicExpr::add(x, y);

    // direction = (3, 4), |dir| = 5, unit = (3/5, 4/5)
    VectorField dir = {SymbolicExpr::number(3), SymbolicExpr::number(4)};

    auto dd = directional_derivative(f, {"x", "y"}, dir);
    EXPECT_TRUE(dd != nullptr, "directional derivative is not null");

    // D_u f = 1 * (3/5) + 1 * (4/5) = 7/5 = 1.4
    auto val = test_numeric_eval(dd);
    if (val.has_value()) {
        EXPECT_NEAR(*val, 1.4, 1e-10, "D_(3,4) (x+y) = 7/5 = 1.4");
    } else {
        // Try substituting (shouldn't have free vars for constant gradient)
        auto dd_s = dd->simplify();
        auto val2 = test_numeric_eval(dd_s);
        EXPECT_TRUE(val2.has_value() && std::abs(*val2 - 1.4) < 1e-10,
            "D_(3,4) (x+y) = 7/5 after simplify");
    }
}

static void test_directional_derivative_zero_vector()
{
    TEST_CASE("directional derivative: zero direction vector returns nullptr");

    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(2));

    VectorField dir = {SymbolicExpr::number(0)};

    auto dd = directional_derivative(f, {"x"}, dir);

    EXPECT_TRUE(dd == nullptr, "zero direction returns nullptr");
}

static void test_directional_derivative_higher_order()
{
    TEST_CASE("directional derivative: order 2, f = x^3, dir = (1)");

    auto x = SymbolicExpr::variable("x");

    // f = x^3
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));

    VectorField dir = {SymbolicExpr::number(1)};

    // First order: 3x^2, second order: 6x
    auto dd2 = directional_derivative(f, {"x"}, dir, 2);
    EXPECT_TRUE(dd2 != nullptr, "second order derivative is not null");

    // Evaluate at x=4: should be 6*4 = 24
    auto dd2_at = dd2->substitute("x", SymbolicExpr::number(4));
    dd2_at = dd2_at->simplify();
    auto val = test_numeric_eval(dd2_at);
    EXPECT_TRUE(val.has_value() && std::abs(*val - 24.0) < 1e-10,
        "D^2_(1) x^3 at x=4 = 6*4 = 24");
}

static void test_div_curl_is_zero()
{
    TEST_CASE("div(curl(F)) = 0 for F = (x*y, y*z, z*x)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    VectorField F = {
        SymbolicExpr::multiply(x, y),
        SymbolicExpr::multiply(y, z),
        SymbolicExpr::multiply(z, x)
    };

    auto c = curl(F, {"x", "y", "z"});
    auto div_curl = divergence(c, {"x", "y", "z"});

    // div(curl(F)) should be 0 — evaluate at a point
    auto val = test_numeric_eval(div_curl);
    if (val.has_value()) {
        EXPECT_NEAR(*val, 0.0, 1e-10, "div(curl(F)) = 0 (numeric)");
    } else {
        auto test_expr = div_curl->substitute("x", SymbolicExpr::number(2));
        test_expr = test_expr->substitute("y", SymbolicExpr::number(3));
        test_expr = test_expr->substitute("z", SymbolicExpr::number(5));
        test_expr = test_expr->simplify();
        auto v = test_numeric_eval(test_expr);
        EXPECT_TRUE(v.has_value() && std::abs(*v) < 1e-10,
            "div(curl(F)) = 0 at (2,3,5)");
    }
}

// ============================================================
// 雅可比矩阵测试
// ============================================================

static std::shared_ptr<SymbolicExpr> get_mat_entry(
    const std::shared_ptr<SymbolicExpr>& mat, size_t r, size_t c)
{
    if (!mat || !mat->root) return nullptr;
    auto mn = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    if (!mn) return nullptr;
    auto node = mn->get(r, c);
    if (!node) return SymbolicExpr::number(0);
    return std::make_shared<SymbolicExpr>(node);
}

static size_t get_mat_rows(const std::shared_ptr<SymbolicExpr>& mat)
{
    if (!mat || !mat->root) return 0;
    auto mn = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    return mn ? mn->rows : 0;
}

static size_t get_mat_cols(const std::shared_ptr<SymbolicExpr>& mat)
{
    if (!mat || !mat->root) return 0;
    auto mn = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    return mn ? mn->cols : 0;
}

static void test_jacobian_square()
{
    TEST_CASE("Jacobian: square 2x2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f1 = x^2 + y^2, f2 = x*y
    auto f1 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto f2 = SymbolicExpr::multiply(x, y);

    auto J = jacobian({f1, f2}, {"x", "y"});

    EXPECT_TRUE(J != nullptr, "Jacobian is not null");
    EXPECT_TRUE(get_mat_rows(J) == 2, "Jacobian has 2 rows");
    EXPECT_TRUE(get_mat_cols(J) == 2, "Jacobian has 2 columns");

    // J[0][0] = 2x, J[0][1] = 2y at (3,5)
    auto j00 = get_mat_entry(J, 0, 0);
    auto j00_at = j00->substitute("x", SymbolicExpr::number(3));
    j00_at = j00_at->simplify();
    auto v00 = test_numeric_eval(j00_at);
    EXPECT_TRUE(v00.has_value() && std::abs(*v00 - 6.0) < 1e-10,
        "J[0][0] = 2x => 6 at x=3");

    auto j01 = get_mat_entry(J, 0, 1);
    auto j01_at = j01->substitute("y", SymbolicExpr::number(5));
    j01_at = j01_at->simplify();
    auto v01 = test_numeric_eval(j01_at);
    EXPECT_TRUE(v01.has_value() && std::abs(*v01 - 10.0) < 1e-10,
        "J[0][1] = 2y => 10 at y=5");

    // J[1][0] = y, J[1][1] = x
    EXPECT_EQ_EXPR(get_mat_entry(J, 1, 0), y, "J[1][0] = y");
    EXPECT_EQ_EXPR(get_mat_entry(J, 1, 1), x, "J[1][1] = x");
}

static void test_jacobian_non_square()
{
    TEST_CASE("Jacobian: non-square 3x2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f1 = x + y, f2 = x*y, f3 = x^2
    auto f1 = SymbolicExpr::add(x, y);
    auto f2 = SymbolicExpr::multiply(x, y);
    auto f3 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto J = jacobian({f1, f2, f3}, {"x", "y"});

    EXPECT_TRUE(get_mat_rows(J) == 3, "3x2 Jacobian has 3 rows");
    EXPECT_TRUE(get_mat_cols(J) == 2, "3x2 Jacobian has 2 columns");

    // J[0][0] = 1, J[0][1] = 1
    auto j00 = get_mat_entry(J, 0, 0);
    auto v00 = test_numeric_eval(j00);
    EXPECT_TRUE(v00.has_value() && std::abs(*v00 - 1.0) < 1e-10,
        "J[0][0] = 1");

    auto j01 = get_mat_entry(J, 0, 1);
    auto v01 = test_numeric_eval(j01);
    EXPECT_TRUE(v01.has_value() && std::abs(*v01 - 1.0) < 1e-10,
        "J[0][1] = 1");

    // J[1][0] = y, J[1][1] = x
    EXPECT_EQ_EXPR(get_mat_entry(J, 1, 0), y, "J[1][0] = y");
    EXPECT_EQ_EXPR(get_mat_entry(J, 1, 1), x, "J[1][1] = x");

    // J[2][0] = 2x at x=4 => 8
    auto j20 = get_mat_entry(J, 2, 0);
    auto j20_at = j20->substitute("x", SymbolicExpr::number(4));
    j20_at = j20_at->simplify();
    auto v20 = test_numeric_eval(j20_at);
    EXPECT_TRUE(v20.has_value() && std::abs(*v20 - 8.0) < 1e-10,
        "J[2][0] = 2x => 8 at x=4");

    // J[2][1] = 0
    auto j21 = get_mat_entry(J, 2, 1);
    EXPECT_TRUE(j21->is_zero(), "J[2][1] = 0");
}

static void test_jacobian_wide()
{
    TEST_CASE("Jacobian: wide 1x3");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    // f = x + 2y + 3z
    auto f = SymbolicExpr::add(x,
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), y),
            SymbolicExpr::multiply(SymbolicExpr::number(3), z)));

    auto J = jacobian({f}, {"x", "y", "z"});

    EXPECT_TRUE(get_mat_rows(J) == 1, "1x3 Jacobian has 1 row");
    EXPECT_TRUE(get_mat_cols(J) == 3, "1x3 Jacobian has 3 columns");

    auto v0 = test_numeric_eval(get_mat_entry(J, 0, 0));
    auto v1 = test_numeric_eval(get_mat_entry(J, 0, 1));
    auto v2 = test_numeric_eval(get_mat_entry(J, 0, 2));
    EXPECT_TRUE(v0.has_value() && std::abs(*v0 - 1.0) < 1e-10,
        "J[0][0] = 1");
    EXPECT_TRUE(v1.has_value() && std::abs(*v1 - 2.0) < 1e-10,
        "J[0][1] = 2");
    EXPECT_TRUE(v2.has_value() && std::abs(*v2 - 3.0) < 1e-10,
        "J[0][2] = 3");
}

// ============================================================
// 海森矩阵测试
// ============================================================

static void test_hessian_quadratic()
{
    TEST_CASE("Hessian: quadratic f = x^2 + 2xy + 3y^2");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 + 2*x*y + 3*y^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2),
                SymbolicExpr::multiply(x, y)),
            SymbolicExpr::multiply(SymbolicExpr::number(3),
                SymbolicExpr::power(y, SymbolicExpr::number(2)))));

    auto H = hessian(f, {"x", "y"});

    EXPECT_TRUE(H != nullptr, "Hessian is not null");
    EXPECT_TRUE(get_mat_rows(H) == 2, "Hessian has 2 rows");
    EXPECT_TRUE(get_mat_cols(H) == 2, "Hessian has 2 columns");

    // H[0][0] = 2, H[0][1] = 2, H[1][0] = 2, H[1][1] = 6
    auto v00 = test_numeric_eval(get_mat_entry(H, 0, 0));
    auto v01 = test_numeric_eval(get_mat_entry(H, 0, 1));
    auto v10 = test_numeric_eval(get_mat_entry(H, 1, 0));
    auto v11 = test_numeric_eval(get_mat_entry(H, 1, 1));

    EXPECT_TRUE(v00.has_value() && std::abs(*v00 - 2.0) < 1e-10,
        "H[0][0] = 2");
    EXPECT_TRUE(v01.has_value() && std::abs(*v01 - 2.0) < 1e-10,
        "H[0][1] = 2");
    EXPECT_TRUE(v10.has_value() && std::abs(*v10 - 2.0) < 1e-10,
        "H[1][0] = 2 (symmetric)");
    EXPECT_TRUE(v11.has_value() && std::abs(*v11 - 6.0) < 1e-10,
        "H[1][1] = 6");
}

static void test_hessian_3d()
{
    TEST_CASE("Hessian: 3D f = x*y*z");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::multiply(
        SymbolicExpr::multiply(x, y), z);

    auto H = hessian(f, {"x", "y", "z"});

    EXPECT_TRUE(get_mat_rows(H) == 3, "3D Hessian has 3 rows");
    EXPECT_TRUE(get_mat_cols(H) == 3, "3D Hessian has 3 columns");

    // Diagonal: all 0 (no x^2, y^2, z^2 terms)
    auto d00 = test_numeric_eval(get_mat_entry(H, 0, 0));
    auto d11 = test_numeric_eval(get_mat_entry(H, 1, 1));
    auto d22 = test_numeric_eval(get_mat_entry(H, 2, 2));
    EXPECT_TRUE(d00.has_value() && std::abs(*d00) < 1e-10,
        "H[0][0] = 0");
    EXPECT_TRUE(d11.has_value() && std::abs(*d11) < 1e-10,
        "H[1][1] = 0");
    EXPECT_TRUE(d22.has_value() && std::abs(*d22) < 1e-10,
        "H[2][2] = 0");

    // Off-diagonal: H[0][1] = z, H[0][2] = y, H[1][2] = x
    // Evaluate at (2, 3, 5)
    auto h01 = get_mat_entry(H, 0, 1);
    auto h01_at = h01->substitute("x", SymbolicExpr::number(2));
    h01_at = h01_at->substitute("y", SymbolicExpr::number(3));
    h01_at = h01_at->substitute("z", SymbolicExpr::number(5));
    h01_at = h01_at->simplify();
    auto v01 = test_numeric_eval(h01_at);
    EXPECT_TRUE(v01.has_value() && std::abs(*v01 - 5.0) < 1e-10,
        "H[0][1] = z => 5 at z=5");

    auto h02 = get_mat_entry(H, 0, 2);
    auto h02_at = h02->substitute("x", SymbolicExpr::number(2));
    h02_at = h02_at->substitute("y", SymbolicExpr::number(3));
    h02_at = h02_at->substitute("z", SymbolicExpr::number(5));
    h02_at = h02_at->simplify();
    auto v02 = test_numeric_eval(h02_at);
    EXPECT_TRUE(v02.has_value() && std::abs(*v02 - 3.0) < 1e-10,
        "H[0][2] = y => 3 at y=3");

    auto h12 = get_mat_entry(H, 1, 2);
    auto h12_at = h12->substitute("x", SymbolicExpr::number(2));
    h12_at = h12_at->substitute("y", SymbolicExpr::number(3));
    h12_at = h12_at->substitute("z", SymbolicExpr::number(5));
    h12_at = h12_at->simplify();
    auto v12 = test_numeric_eval(h12_at);
    EXPECT_TRUE(v12.has_value() && std::abs(*v12 - 2.0) < 1e-10,
        "H[1][2] = x => 2 at x=2");

    // Symmetry: H[1][0] = H[0][1], H[2][0] = H[0][2], H[2][1] = H[1][2]
    auto h10 = get_mat_entry(H, 1, 0);
    auto h10_at = h10->substitute("z", SymbolicExpr::number(5));
    h10_at = h10_at->simplify();
    auto v10 = test_numeric_eval(h10_at);
    EXPECT_TRUE(v10.has_value() && std::abs(*v10 - 5.0) < 1e-10,
        "H[1][0] = z => 5 (symmetric)");
}

static void test_hessian_single_var()
{
    TEST_CASE("Hessian: single variable f = x^3");

    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));

    auto H = hessian(f, {"x"});

    EXPECT_TRUE(get_mat_rows(H) == 1, "1x1 Hessian has 1 row");
    EXPECT_TRUE(get_mat_cols(H) == 1, "1x1 Hessian has 1 column");

    // H[0][0] = 6x at x=2 => 12
    auto h00 = get_mat_entry(H, 0, 0);
    auto h00_at = h00->substitute("x", SymbolicExpr::number(2));
    h00_at = h00_at->simplify();
    auto v = test_numeric_eval(h00_at);
    EXPECT_TRUE(v.has_value() && std::abs(*v - 12.0) < 1e-10,
        "H[0][0] = 6x => 12 at x=2");
}

// ============================================================
// 曲线积分测试 (Requirements 48)
// ============================================================

static void test_curve_integral_scalar_line()
{
    TEST_CASE("curve_integral_scalar: f=1 along line segment (arc length)");

    // Parametrization: r(t) = (t, 2t) for t in [0, 1]
    // |r'(t)| = sqrt(1 + 4) = sqrt(5)
    // integral_0^1 1 * sqrt(5) dt = sqrt(5)
    auto t_var = SymbolicExpr::variable("t");
    VectorField param = {
        t_var,
        SymbolicExpr::multiply(SymbolicExpr::number(2), t_var)
    };

    auto f = SymbolicExpr::number(1);
    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_scalar(f, param, "t", a, b);
    EXPECT_TRUE(result != nullptr, "curve_integral_scalar result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        double expected = std::sqrt(5.0);
        EXPECT_TRUE(val.has_value(), "arc length result is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, expected, 1e-6,
                "arc length of r(t)=(t,2t) on [0,1] = sqrt(5)");
        }
    }
}

static void test_curve_integral_scalar_circle()
{
    TEST_CASE("curve_integral_scalar: f=1 along unit circle (circumference)");

    // Parametrization: r(t) = (cos(t), sin(t)) for t in [0, 2*pi]
    // |r'(t)| = sqrt(sin^2(t) + cos^2(t)) = 1
    // integral_0^{2pi} 1 * 1 dt = 2*pi
    auto t_var = SymbolicExpr::variable("t");
    VectorField param = {
        SymbolicExpr::cos(t_var),
        SymbolicExpr::sin(t_var)
    };

    auto f = SymbolicExpr::number(1);
    auto a = SymbolicExpr::number(0);
    auto pi2 = SymbolicExpr::number(2.0 * 3.14159265358979323846);

    auto result = curve_integral_scalar(f, param, "t", a, pi2);
    EXPECT_TRUE(result != nullptr, "circle arc length result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        double expected = 2.0 * 3.14159265358979323846;
        EXPECT_TRUE(val.has_value(), "circle arc length is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, expected, 1e-3,
                "circumference of unit circle = 2*pi");
        }
    }
}

static void test_curve_integral_vector_conservative()
{
    TEST_CASE("curve_integral_vector: F=grad(x^2+y^2) along path");

    // F = (2x, 2y) = grad(x^2 + y^2)
    // Parametrization: r(t) = (t, t) for t in [0, 1]
    // integral_0^1 F(r(t)) dot r'(t) dt = integral_0^1 (2t, 2t) dot (1, 1) dt
    //   = integral_0^1 4t dt = 2
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::multiply(SymbolicExpr::number(2), y)
    };

    VectorField param = {t_var, t_var};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE(result != nullptr, "curve_integral_vector result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "conservative field integral is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 2.0, 1e-6,
                "work of grad(x^2+y^2) along (t,t) = 2");
        }
    }
}

static void test_curve_integral_vector_work()
{
    TEST_CASE("curve_integral_vector: work done by F=(y, -x) along x-axis");

    // F = (y, -x)
    // Parametrization: r(t) = (t, 0) for t in [0, 1]
    // F(r(t)) = (0, -t), r'(t) = (1, 0)
    // integral_0^1 (0, -t) dot (1, 0) dt = integral_0^1 0 dt = 0
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        y,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    };

    VectorField param = {t_var, SymbolicExpr::number(0)};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE(result != nullptr, "work integral result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "work integral is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-6,
                "work of (y,-x) along x-axis = 0");
        }
    }
}

static void test_curve_integral_vector_3d()
{
    TEST_CASE("curve_integral_vector: 3D F=(1,0,0) along x-axis");

    // F = (1, 0, 0), r(t) = (t, 0, 0) for t in [0, 3]
    // integral_0^3 (1,0,0) dot (1,0,0) dt = 3
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto t_var = SymbolicExpr::variable("t");

    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)
    };

    VectorField param = {t_var, SymbolicExpr::number(0), SymbolicExpr::number(0)};

    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(3);

    auto result = curve_integral_vector(F, param, "t", a, b);
    EXPECT_TRUE(result != nullptr, "3D curve integral result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "3D curve integral is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 3.0, 1e-6,
                "work of (1,0,0) along x-axis [0,3] = 3");
        }
    }
}

// ============================================================
// 曲面积分测试 (Requirements 49)
// ============================================================

static void test_surface_integral_scalar_plane()
{
    TEST_CASE("surface_integral_scalar: f=1 on flat square (area)");

    // Parametrization: r(u,v) = (u, v, 0) for u,v in [0,1]
    // r_u = (1, 0, 0), r_v = (0, 1, 0)
    // r_u x r_v = (0, 0, 1), |r_u x r_v| = 1
    // double integral of 1 * 1 du dv = 1
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto f = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_scalar(f, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE(result != nullptr, "surface_integral_scalar result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "surface area is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6,
                "area of unit square = 1");
        }
    }
}

static void test_surface_integral_vector_flux()
{
    TEST_CASE("surface_integral_vector: flux of F=(0,0,1) through flat square");

    // F = (0, 0, 1), r(u,v) = (u, v, 0) for u,v in [0,1]
    // r_u x r_v = (0, 0, 1)
    // F dot (r_u x r_v) = 1
    // double integral of 1 du dv = 1
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField F = {
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1)
    };

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_vector(F, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE(result != nullptr, "surface_integral_vector result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "flux is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 1.0, 1e-6,
                "flux of (0,0,1) through unit square = 1");
        }
    }
}

static void test_surface_integral_vector_zero_flux()
{
    TEST_CASE("surface_integral_vector: zero flux of tangent field");

    // F = (1, 0, 0) (tangent to the surface), r(u,v) = (u, v, 0)
    // r_u x r_v = (0, 0, 1)
    // F dot (r_u x r_v) = 0
    // double integral of 0 du dv = 0
    auto u_var = SymbolicExpr::variable("u");
    auto v_var = SymbolicExpr::variable("v");

    VectorField F = {
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(0)
    };

    VectorField param = {u_var, v_var, SymbolicExpr::number(0)};

    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);

    auto result = surface_integral_vector(F, param, "u", "v", zero, one, zero, one);
    EXPECT_TRUE(result != nullptr, "zero flux result is not null");

    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "zero flux is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.0, 1e-6,
                "flux of tangent field through surface = 0");
        }
    }
}

// ============================================================
// 多元极值测试
// ============================================================

static void test_find_extrema_quadratic_min()
{
    TEST_CASE("find_extrema: f = x^2 + y^2 (minimum at origin)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 + y^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE(!extrema.empty(), "found at least one critical point");
    if (!extrema.empty()) {
        EXPECT_TRUE(extrema[0].classification == "minimum",
            "x^2+y^2 has minimum at origin, got: " + extrema[0].classification);

        // Check the point is (0, 0)
        auto x_val = extrema[0].point.at("x");
        auto y_val = extrema[0].point.at("y");
        EXPECT_TRUE(x_val && x_val->is_zero(), "critical point x = 0");
        EXPECT_TRUE(y_val && y_val->is_zero(), "critical point y = 0");
    }
}

static void test_find_extrema_quadratic_max()
{
    TEST_CASE("find_extrema: f = -x^2 - y^2 (maximum at origin)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = -x^2 - y^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(x, SymbolicExpr::number(2))),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE(!extrema.empty(), "found at least one critical point");
    if (!extrema.empty()) {
        EXPECT_TRUE(extrema[0].classification == "maximum",
            "-x^2-y^2 has maximum at origin, got: " + extrema[0].classification);
    }
}

static void test_find_extrema_saddle()
{
    TEST_CASE("find_extrema: f = x^2 - y^2 (saddle at origin)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 - y^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE(!extrema.empty(), "found at least one critical point");
    if (!extrema.empty()) {
        EXPECT_TRUE(extrema[0].classification == "saddle",
            "x^2-y^2 has saddle at origin, got: " + extrema[0].classification);
    }
}

static void test_find_extrema_degenerate()
{
    TEST_CASE("find_extrema: f = x^2 (degenerate in y direction)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // f = x^2 (Hessian is [[2, 0], [0, 0]] — singular)
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto extrema = find_extrema(f, {"x", "y"});

    // The system df/dx = 2x = 0, df/dy = 0 has solution x=0, y=free
    if (!extrema.empty()) {
        EXPECT_TRUE(extrema[0].classification == "degenerate",
            "x^2 in 2D has degenerate critical point, got: " + extrema[0].classification);
    }
}

static void test_find_extrema_single_var()
{
    TEST_CASE("find_extrema: f = x^2 - 2x + 1 (single variable min at x=1)");

    auto x = SymbolicExpr::variable("x");

    // f = x^2 - 2x + 1 = (x-1)^2
    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));

    auto extrema = find_extrema(f, {"x"});

    EXPECT_TRUE(!extrema.empty(), "found critical point for (x-1)^2");
    if (!extrema.empty()) {
        EXPECT_TRUE(extrema[0].classification == "minimum",
            "(x-1)^2 has minimum, got: " + extrema[0].classification);

        auto x_val = extrema[0].point.at("x");
        if (x_val) {
            auto val = x_val->simplify();
            auto num = test_numeric_eval(val);
            EXPECT_TRUE(num.has_value() && std::abs(*num - 1.0) < 1e-10,
                "critical point at x=1");
        }
    }
}

// ============================================================
// 拉格朗日乘数法测试
// ============================================================

static void test_lagrange_basic()
{
    TEST_CASE("lagrange_multipliers: max x+y subject to x^2+y^2=1");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // Maximize f = x + y subject to g = x^2 + y^2 - 1 = 0
    auto f = SymbolicExpr::add(x, y);
    auto g = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::number(-1));

    auto solutions = lagrange_multipliers(f, {g}, {"x", "y"});

    EXPECT_TRUE(!solutions.empty(), "found at least one critical point");

    // The solutions should be (1/sqrt(2), 1/sqrt(2)) and (-1/sqrt(2), -1/sqrt(2))
    bool found_max = false;
    for (const auto& sol : solutions) {
        if (sol.count("x") && sol.count("y")) {
            auto x_val = sol.at("x")->simplify();
            auto y_val = sol.at("y")->simplify();
            auto xn = test_numeric_eval(x_val);
            auto yn = test_numeric_eval(y_val);
            if (xn.has_value() && yn.has_value()) {
                // Check if x approx equals y (both solutions have x = y)
                if (std::abs(*xn - *yn) < 1e-8) {
                    found_max = true;
                }
            }
        }
    }
    EXPECT_TRUE(found_max, "found solution where x = y");
}

static void test_lagrange_linear_constraint()
{
    TEST_CASE("lagrange_multipliers: min x^2+y^2 subject to x+y=1");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // Minimize f = x^2 + y^2 subject to g = x + y - 1 = 0
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto g = SymbolicExpr::add(
        SymbolicExpr::add(x, y),
        SymbolicExpr::number(-1));

    auto solutions = lagrange_multipliers(f, {g}, {"x", "y"});

    EXPECT_TRUE(!solutions.empty(), "found critical point");

    // Solution should be x = y = 1/2
    if (!solutions.empty()) {
        auto x_val = solutions[0].at("x")->simplify();
        auto y_val = solutions[0].at("y")->simplify();
        auto xn = test_numeric_eval(x_val);
        auto yn = test_numeric_eval(y_val);
        if (xn.has_value() && yn.has_value()) {
            EXPECT_NEAR(*xn, 0.5, 1e-8, "x = 1/2");
            EXPECT_NEAR(*yn, 0.5, 1e-8, "y = 1/2");
        }
    }
}


int main()
{
    test_gradient_basic();
    test_gradient_mixed();
    test_divergence_basic();
    test_divergence_constant_field();
    test_curl_3d();
    test_curl_2d();
    test_curl_grad_is_zero();
    test_laplacian_basic();
    test_laplacian_harmonic();
    test_laplacian_equals_div_grad();
    test_directional_derivative_axis();
    test_directional_derivative_normalization();
    test_directional_derivative_zero_vector();
    test_directional_derivative_higher_order();
    test_div_curl_is_zero();

    // Jacobian tests
    test_jacobian_square();
    test_jacobian_non_square();
    test_jacobian_wide();

    // Hessian tests
    test_hessian_quadratic();
    test_hessian_3d();
    test_hessian_single_var();

    // Curve integral tests
    test_curve_integral_scalar_line();
    test_curve_integral_scalar_circle();
    test_curve_integral_vector_conservative();
    test_curve_integral_vector_work();
    test_curve_integral_vector_3d();

    // Surface integral tests
    test_surface_integral_scalar_plane();
    test_surface_integral_vector_flux();
    test_surface_integral_vector_zero_flux();

    // Extrema tests
    test_find_extrema_quadratic_min();
    test_find_extrema_quadratic_max();
    test_find_extrema_saddle();
    test_find_extrema_degenerate();
    test_find_extrema_single_var();

    // Lagrange multiplier tests
    test_lagrange_basic();
    test_lagrange_linear_constraint();

    return TEST_REPORT();
}
