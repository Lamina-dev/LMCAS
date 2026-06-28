#include "differential_geometry.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

int main() {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    std::vector<std::string> vars = {"x", "y", "z"};

    // ---- Lie derivative of f = x^2 along X = (1,0,0) is 2x ----
    {
        auto f = SymbolicExpr::multiply(x, x);
        std::vector<std::shared_ptr<SymbolicExpr>> X = {num(1), num(0), num(0)};
        auto L = lie_derivative(f, X, vars, 1);
        auto expected = SymbolicExpr::multiply(num(2), x);
        EXPECT_EQ_EXPR(L->simplify(), expected->simplify(), "L_(1,0,0)(x^2) = 2x");
    }

    // ---- Lie derivative of f = x*y along X = (y, x, 0) = y*y + x*x ----
    {
        auto f = SymbolicExpr::multiply(x, y);
        std::vector<std::shared_ptr<SymbolicExpr>> X = {y, x, num(0)};
        auto L = lie_derivative(f, X, vars, 1);
        // d/dx(xy)=y, d/dy(xy)=x => y*y + x*x
        auto expected = SymbolicExpr::add(SymbolicExpr::multiply(y, y), SymbolicExpr::multiply(x, x));
        EXPECT_EQ_EXPR(L->simplify(), expected->simplify(), "L_(y,x,0)(xy) = y^2 + x^2");
    }

    // ---- exterior derivative of 0-form f = x^2 + y^2 + z^2 equals gradient (2x,2y,2z) ----
    {
        auto f = SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                 SymbolicExpr::add(SymbolicExpr::multiply(y, y), SymbolicExpr::multiply(z, z)));
        auto d = exterior_derivative({f}, 0, vars);
        EXPECT_TRUE(d.size() == 3, "d(0-form) has 3 components");
        EXPECT_EQ_EXPR(d[0]->simplify(), SymbolicExpr::multiply(num(2), x)->simplify(), "d f / dx = 2x");
        EXPECT_EQ_EXPR(d[1]->simplify(), SymbolicExpr::multiply(num(2), y)->simplify(), "d f / dy = 2y");
        EXPECT_EQ_EXPR(d[2]->simplify(), SymbolicExpr::multiply(num(2), z)->simplify(), "d f / dz = 2z");
    }

    // ---- d^2 = 0: exterior derivative of (d of a 0-form) is the zero 2-form ----
    {
        auto f = SymbolicExpr::multiply(x, y); // smooth scalar
        auto d1 = exterior_derivative({f}, 0, vars); // 1-form coeffs
        auto d2 = exterior_derivative(d1, 1, vars);   // 2-form coeffs, should all be 0
        EXPECT_TRUE(d2.size() == 3, "2-form has 3 components in 3D");
        bool all_zero = true;
        for (auto& c : d2) if (!(c->simplify()->root && c->simplify()->root->is_zero())) all_zero = false;
        EXPECT_TRUE(all_zero, "d^2(f) = 0");
    }

    // ---- existing Christoffel smoke check still works ----
    {
        // metric_inverse expects a matrix; use a 2x2 identity metric.
        auto g = SymbolicExpr::matrix({{num(1), num(0)}, {num(0), num(1)}});
        auto mi = metric_inverse(g);
        EXPECT_TRUE(mi != nullptr, "metric_inverse returns non-null");
    }

    return TEST_REPORT();
}
