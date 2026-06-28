#include "numerical_integration.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"
#include <cmath>

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static bool close(const std::shared_ptr<SymbolicExpr>& e, double expected, double tol = 1e-6) {
    if (!e) return false;
    auto s = e->simplify();
    if (s->is_number()) return std::abs(s->to_numeric() - expected) < tol;
    // Fall back to recursive numeric evaluation (handles sqrt, etc.).
    auto v = test_numeric_eval(s);
    if (v.has_value()) return std::abs(*v - expected) < tol;
    return false;
}

int main() {
    auto x = SymbolicExpr::variable("x");

    // ∫₀¹ x dx = 1/2 (Simpson exact for linear)
    {
        auto r = quadrature_simpson(x, "x", num(0), num(1), 10);
        EXPECT_TRUE(close(r, 0.5), "Simpson ∫₀¹ x dx = 0.5");
    }

    // ∫₀¹ x^2 dx = 1/3 (Simpson exact for quadratic)
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = quadrature_simpson(f, "x", num(0), num(1), 10);
        EXPECT_TRUE(close(r, 1.0/3.0), "Simpson ∫₀¹ x² dx = 1/3");
    }

    // ∫₀¹ x^3 dx = 1/4 via Gauss-Legendre (n=2 exact up to degree 3)
    {
        auto f = SymbolicExpr::power(x, num(3));
        auto r = quadrature_gaussian(f, "x", num(0), num(1), 2);
        EXPECT_TRUE(close(r, 0.25, 1e-6), "Gauss ∫₀¹ x³ dx = 1/4");
    }

    // adaptive_simpson on ∫₀¹ x^2 dx = 1/3
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = adaptive_simpson(f, "x", num(0), num(1), 1e-10);
        EXPECT_TRUE(close(r, 1.0/3.0, 1e-8), "adaptive Simpson ∫₀¹ x² dx = 1/3");
    }

    // adaptive_simpson on a transcendental: ∫₀^π sin(x) dx = 2
    {
        auto f = SymbolicExpr::sin(x);
        auto pi = SymbolicExpr::number(LMMC_CONST_PI);
        auto r = adaptive_simpson(f, "x", num(0), pi, 1e-10);
        EXPECT_TRUE(close(r, 2.0, 1e-6), "adaptive Simpson ∫₀^π sin x dx = 2");
    }

    // numerical_integrate wrapper
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = numerical_integrate(f, "x", num(0), num(1), 100);
        EXPECT_TRUE(close(r, 1.0/3.0, 1e-6), "numerical_integrate ∫₀¹ x² dx = 1/3");
    }

    return TEST_REPORT();
}
