/**
 * @file test_prop_vector_calculus.cpp
 * @brief Property-based tests for vector calculus module.
 *
 * Properties tested:
 *   - Property 6: Directional derivative equals gradient dot unit direction
 *   - Property 7: Gradient components equal partial derivatives
 *   - Property 8: Jacobian entry correctness
 *   - Property 9: Hessian symmetry and correctness
 *   - Property 22: curl(grad(f)) = 0
 *   - Property 23: div(curl(F)) = 0
 *   - Property 24: Laplacian equals divergence of gradient
 *
 * **Validates: Requirements 8.1-11.3, 45.1-47.3, 88.1-88.3**
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "vector_calculus.hpp"
#include "symbolic_ast.hpp"

#include <cmath>
#include <algorithm>

using namespace lamina;
using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string& name) { return SE::variable(name); }

// ============================================================
// Helpers: matrix entry access
// ============================================================

static std::shared_ptr<SymbolicExpr> get_mat_entry(
    const std::shared_ptr<SymbolicExpr>& mat, size_t r, size_t c)
{
    if (!mat || !mat->root) return nullptr;
    auto mn = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    if (!mn) return nullptr;
    auto node = mn->get(r, c);
    if (!node) return SE::number(0);
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

// ============================================================
// Generators: random polynomial expressions in multiple variables
// ============================================================

namespace {

/**
 * @brief Generate a random polynomial in variables {x, y, z} with
 *        integer coefficients in [-3, 3] and total degree in [1, 3].
 *
 * Generates monomials like: a*x^i*y^j*z^k where i+j+k <= max_degree.
 */
std::shared_ptr<SymbolicExpr> gen_poly_3d(int max_degree = 3) {
    auto x = var("x");
    auto y = var("y");
    auto z = var("z");

    std::shared_ptr<SymbolicExpr> result = nullptr;
    int num_terms = rc::gen::inRange(2, 5);

    for (int t = 0; t < num_terms; ++t) {
        int coeff = rc::gen::inRange(-3, 3);
        if (coeff == 0) coeff = 1;

        int ix = rc::gen::inRange(0, max_degree);
        int iy = rc::gen::inRange(0, max_degree - ix);
        int iz = rc::gen::inRange(0, max_degree - ix - iy);

        std::shared_ptr<SymbolicExpr> term = num(coeff);

        if (ix > 0) {
            auto xp = (ix == 1) ? x : SE::power(x, num(ix));
            term = SE::multiply(term, xp);
        }
        if (iy > 0) {
            auto yp = (iy == 1) ? y : SE::power(y, num(iy));
            term = SE::multiply(term, yp);
        }
        if (iz > 0) {
            auto zp = (iz == 1) ? z : SE::power(z, num(iz));
            term = SE::multiply(term, zp);
        }

        if (!result) {
            result = term;
        } else {
            result = SE::add(result, term);
        }
    }

    return result ? result : num(1);
}

/**
 * @brief Generate a random 3D vector field with polynomial components.
 */
VectorField gen_vector_field_3d(int max_degree = 2) {
    return {gen_poly_3d(max_degree), gen_poly_3d(max_degree), gen_poly_3d(max_degree)};
}

/**
 * @brief Generate a random non-zero direction vector with integer components in [-3, 3].
 *        Guarantees at least one non-zero component.
 */
VectorField gen_direction_3d() {
    int d0 = rc::gen::inRange(-3, 3);
    int d1 = rc::gen::inRange(-3, 3);
    int d2 = rc::gen::inRange(-3, 3);

    // Ensure not all zero
    if (d0 == 0 && d1 == 0 && d2 == 0) {
        d0 = 1;
    }

    return {num(d0), num(d1), num(d2)};
}

/**
 * @brief Evaluate a symbolic expression at a specific 3D point (x, y, z).
 *        Returns std::nullopt if evaluation fails.
 */
std::optional<double> eval_at_point(const std::shared_ptr<SymbolicExpr>& expr,
                                     double px, double py, double pz)
{
    if (!expr) return std::nullopt;
    auto e = expr->substitute("x", SE::number(px));
    e = e->substitute("y", SE::number(py));
    e = e->substitute("z", SE::number(pz));
    e = e->simplify();
    return test_numeric_eval(e);
}

} // anonymous namespace

// ============================================================
// Property 6: Directional derivative equals gradient dot unit direction
// **Validates: Requirements 8.1, 8.2, 88.1, 88.2**
// ============================================================

static void test_property_6_directional_derivative() {
    TEST_CASE("Property 6: Directional derivative equals gradient dot unit direction");

    rc::check("D_u(f) == grad(f) . unit(u) for random polynomials", []() {
        auto f = gen_poly_3d(2);
        auto dir = gen_direction_3d();
        std::vector<std::string> vars = {"x", "y", "z"};

        auto dd = directional_derivative(f, vars, dir);
        RC_ASSERT(dd != nullptr);

        auto grad_f = gradient(f, vars);
        RC_ASSERT(grad_f.size() == 3);

        // Compute |dir| for normalization
        double d0 = 0, d1 = 0, d2 = 0;
        auto v0 = test_numeric_eval(dir[0]);
        auto v1 = test_numeric_eval(dir[1]);
        auto v2 = test_numeric_eval(dir[2]);
        RC_ASSERT(v0.has_value() && v1.has_value() && v2.has_value());
        d0 = *v0; d1 = *v1; d2 = *v2;
        double mag = std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
        RC_ASSERT(mag > 1e-10);
        double u0 = d0 / mag, u1 = d1 / mag, u2 = d2 / mag;

        // Evaluate at a random point
        double px = 1.0 + rc::gen::inRange(0, 3);
        double py = 1.0 + rc::gen::inRange(0, 3);
        double pz = 1.0 + rc::gen::inRange(0, 3);

        auto dd_val = eval_at_point(dd, px, py, pz);
        auto g0_val = eval_at_point(grad_f[0], px, py, pz);
        auto g1_val = eval_at_point(grad_f[1], px, py, pz);
        auto g2_val = eval_at_point(grad_f[2], px, py, pz);

        if (dd_val && g0_val && g1_val && g2_val) {
            double expected = (*g0_val) * u0 + (*g1_val) * u1 + (*g2_val) * u2;
            if (std::isfinite(*dd_val) && std::isfinite(expected)) {
                double diff = std::abs(*dd_val - expected);
                double scale = std::max(1.0, std::abs(expected));
                RC_ASSERT(diff / scale < 1e-4);
            }
        }
    });
}

// ============================================================
// Property 7: Gradient components equal partial derivatives
// **Validates: Requirements 9.1, 9.2, 9.3**
// ============================================================

static void test_property_7_gradient_components() {
    TEST_CASE("Property 7: Gradient components equal partial derivatives");

    rc::check("grad(f)[i] == df/dxi for random polynomials", []() {
        auto f = gen_poly_3d(3);
        std::vector<std::string> vars = {"x", "y", "z"};

        auto grad_f = gradient(f, vars);
        RC_ASSERT(grad_f.size() == 3);

        // Evaluate at a random point
        double px = 1.0 + rc::gen::inRange(0, 4);
        double py = 1.0 + rc::gen::inRange(0, 4);
        double pz = 1.0 + rc::gen::inRange(0, 4);

        for (int i = 0; i < 3; ++i) {
            auto partial = f->differentiate(vars[i]);
            RC_ASSERT(partial != nullptr);

            auto grad_val = eval_at_point(grad_f[i], px, py, pz);
            auto partial_val = eval_at_point(partial, px, py, pz);

            if (grad_val && partial_val) {
                if (std::isfinite(*grad_val) && std::isfinite(*partial_val)) {
                    double diff = std::abs(*grad_val - *partial_val);
                    double scale = std::max(1.0, std::abs(*partial_val));
                    RC_ASSERT(diff / scale < 1e-6);
                }
            }
        }
    });
}

// ============================================================
// Property 8: Jacobian entry correctness
// **Validates: Requirements 10.1, 10.2, 10.3**
// ============================================================

static void test_property_8_jacobian_entries() {
    TEST_CASE("Property 8: Jacobian entry correctness");

    rc::check("J[i][j] == d(fi)/d(xj) for random polynomial vector functions", []() {
        std::vector<std::string> vars = {"x", "y", "z"};

        // Generate 2-3 random polynomial functions
        int m = rc::gen::inRange(2, 3);
        std::vector<std::shared_ptr<SymbolicExpr>> funcs;
        for (int i = 0; i < m; ++i) {
            funcs.push_back(gen_poly_3d(2));
        }

        auto J = jacobian(funcs, vars);
        RC_ASSERT(J != nullptr);
        RC_ASSERT(get_mat_rows(J) == (size_t)m);
        RC_ASSERT(get_mat_cols(J) == 3);

        // Evaluate at a random point
        double px = 1.0 + rc::gen::inRange(0, 3);
        double py = 1.0 + rc::gen::inRange(0, 3);
        double pz = 1.0 + rc::gen::inRange(0, 3);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto entry = get_mat_entry(J, i, j);
                auto expected = funcs[i]->differentiate(vars[j]);

                auto entry_val = eval_at_point(entry, px, py, pz);
                auto expected_val = eval_at_point(expected, px, py, pz);

                if (entry_val && expected_val) {
                    if (std::isfinite(*entry_val) && std::isfinite(*expected_val)) {
                        double diff = std::abs(*entry_val - *expected_val);
                        double scale = std::max(1.0, std::abs(*expected_val));
                        RC_ASSERT(diff / scale < 1e-6);
                    }
                }
            }
        }
    });
}

// ============================================================
// Property 9: Hessian symmetry and correctness
// **Validates: Requirements 11.1, 11.2, 11.3**
// ============================================================

static void test_property_9_hessian_symmetry() {
    TEST_CASE("Property 9: Hessian symmetry and correctness");

    rc::check("H[i][j] == H[j][i] and H[i][j] == d2f/(dxi dxj) for random polynomials", []() {
        auto f = gen_poly_3d(3);
        std::vector<std::string> vars = {"x", "y", "z"};

        auto H = hessian(f, vars);
        RC_ASSERT(H != nullptr);
        RC_ASSERT(get_mat_rows(H) == 3);
        RC_ASSERT(get_mat_cols(H) == 3);

        // Evaluate at a random point
        double px = 1.0 + rc::gen::inRange(0, 3);
        double py = 1.0 + rc::gen::inRange(0, 3);
        double pz = 1.0 + rc::gen::inRange(0, 3);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto hij = get_mat_entry(H, i, j);
                auto hij_val = eval_at_point(hij, px, py, pz);

                // Check correctness: H[i][j] == d2f/(dxi dxj)
                auto df_dxi = f->differentiate(vars[i]);
                auto d2f_dxidxj = df_dxi->differentiate(vars[j]);
                auto expected_val = eval_at_point(d2f_dxidxj, px, py, pz);

                if (hij_val && expected_val) {
                    if (std::isfinite(*hij_val) && std::isfinite(*expected_val)) {
                        double diff = std::abs(*hij_val - *expected_val);
                        double scale = std::max(1.0, std::abs(*expected_val));
                        RC_ASSERT(diff / scale < 1e-6);
                    }
                }

                // Check symmetry: H[i][j] == H[j][i]
                if (i != j) {
                    auto hji = get_mat_entry(H, j, i);
                    auto hji_val = eval_at_point(hji, px, py, pz);

                    if (hij_val && hji_val) {
                        if (std::isfinite(*hij_val) && std::isfinite(*hji_val)) {
                            double diff = std::abs(*hij_val - *hji_val);
                            RC_ASSERT(diff < 1e-6);
                        }
                    }
                }
            }
        }
    });
}

// ============================================================
// Property 22: curl(grad(f)) = 0
// **Validates: Requirements 45.1, 45.2, 45.3, 46.1, 46.2**
// ============================================================

static void test_property_22_curl_grad_zero() {
    TEST_CASE("Property 22: curl(grad(f)) = 0");

    rc::check("curl(grad(f)) == 0 for random polynomial scalar fields", []() {
        auto f = gen_poly_3d(3);
        std::vector<std::string> vars = {"x", "y", "z"};

        auto grad_f = gradient(f, vars);
        RC_ASSERT(grad_f.size() == 3);

        auto curl_grad = curl(grad_f, vars);
        RC_ASSERT(curl_grad.size() == 3);

        // Evaluate each component at a random point — should be zero
        double px = 1.0 + rc::gen::inRange(0, 4);
        double py = 1.0 + rc::gen::inRange(0, 4);
        double pz = 1.0 + rc::gen::inRange(0, 4);

        for (int i = 0; i < 3; ++i) {
            auto val = eval_at_point(curl_grad[i], px, py, pz);
            if (val) {
                RC_ASSERT(std::abs(*val) < 1e-6);
            }
        }
    });
}

// ============================================================
// Property 23: div(curl(F)) = 0
// **Validates: Requirements 45.1, 45.3, 46.1, 46.3**
// ============================================================

static void test_property_23_div_curl_zero() {
    TEST_CASE("Property 23: div(curl(F)) = 0");

    rc::check("div(curl(F)) == 0 for random polynomial vector fields", []() {
        auto F = gen_vector_field_3d(2);
        std::vector<std::string> vars = {"x", "y", "z"};

        auto curl_F = curl(F, vars);
        RC_ASSERT(curl_F.size() == 3);

        auto div_curl = divergence(curl_F, vars);
        RC_ASSERT(div_curl != nullptr);

        // Evaluate at a random point — should be zero
        double px = 1.0 + rc::gen::inRange(0, 4);
        double py = 1.0 + rc::gen::inRange(0, 4);
        double pz = 1.0 + rc::gen::inRange(0, 4);

        auto val = eval_at_point(div_curl, px, py, pz);
        if (val) {
            RC_ASSERT(std::abs(*val) < 1e-6);
        }
    });
}

// ============================================================
// Property 24: Laplacian equals divergence of gradient
// **Validates: Requirements 47.1, 47.2, 47.3, 88.1, 88.3**
// ============================================================

static void test_property_24_laplacian_div_grad() {
    TEST_CASE("Property 24: Laplacian equals divergence of gradient");

    rc::check("laplacian(f) == div(grad(f)) for random polynomial scalar fields", []() {
        auto f = gen_poly_3d(3);
        std::vector<std::string> vars = {"x", "y", "z"};

        auto lap = laplacian(f, vars);
        RC_ASSERT(lap != nullptr);

        auto grad_f = gradient(f, vars);
        RC_ASSERT(grad_f.size() == 3);

        auto div_grad = divergence(grad_f, vars);
        RC_ASSERT(div_grad != nullptr);

        // Evaluate both at a random point — should be equal
        double px = 1.0 + rc::gen::inRange(0, 4);
        double py = 1.0 + rc::gen::inRange(0, 4);
        double pz = 1.0 + rc::gen::inRange(0, 4);

        auto lap_val = eval_at_point(lap, px, py, pz);
        auto dg_val = eval_at_point(div_grad, px, py, pz);

        if (lap_val && dg_val) {
            if (std::isfinite(*lap_val) && std::isfinite(*dg_val)) {
                double diff = std::abs(*lap_val - *dg_val);
                double scale = std::max(1.0, std::abs(*lap_val));
                RC_ASSERT(diff / scale < 1e-6);
            }
        }
    });
}

// ============================================================
// Property 20: Vector projection perpendicularity
// ============================================================

static void test_property_20_projection_perpendicular() {
    TEST_CASE("Property 20: Vector projection perpendicularity");

    rc::check("dot(a - proj(a, b), b) == 0 for random vectors", []() {
        auto a = gen_vector_field_3d(0);
        auto b = gen_direction_3d();
        std::vector<std::string> vars = {"x", "y", "z"};

        auto proj = vector_project(a, b);
        RC_ASSERT(proj.size() == 3);

        VectorField diff = {
            SE::add(a[0], SE::multiply(num(-1), proj[0])),
            SE::add(a[1], SE::multiply(num(-1), proj[1])),
            SE::add(a[2], SE::multiply(num(-1), proj[2]))
        };

        auto d = dot_product(diff, b);
        auto d_val = test_numeric_eval(d->simplify());
        if (d_val) {
            RC_ASSERT(std::abs(*d_val) < 1e-6);
        }
    });
}

// ============================================================
// Property 21: Mixed product equals determinant
// ============================================================

static void test_property_21_mixed_product_det() {
    TEST_CASE("Property 21: Mixed product equals determinant");

    rc::check("dot(a, cross(b, c)) == det([a, b, c]) for random vectors", []() {
        auto a = gen_direction_3d();
        auto b = gen_direction_3d();
        auto c = gen_direction_3d();

        auto cr = cross_product(b, c);
        auto dp = dot_product(a, cr);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_elems = {
            {a[0], a[1], a[2]},
            {b[0], b[1], b[2]},
            {c[0], c[1], c[2]}
        };
        auto mat = SE::matrix(mat_elems);
        auto d = SE::determinant(mat);

        auto dp_val = test_numeric_eval(dp->simplify());
        auto d_val = test_numeric_eval(d->simplify());

        if (dp_val && d_val) {
            RC_ASSERT(std::abs(*dp_val - *d_val) < 1e-6);
        }
    });
}

// ============================================================
// Main
// ============================================================

int main() {
    test_property_6_directional_derivative();
    test_property_7_gradient_components();
    test_property_8_jacobian_entries();
    test_property_9_hessian_symmetry();
    test_property_20_projection_perpendicular();
    test_property_21_mixed_product_det();
    test_property_22_curl_grad_zero();
    test_property_23_div_curl_zero();
    test_property_24_laplacian_div_grad();

    return TEST_REPORT();
}
