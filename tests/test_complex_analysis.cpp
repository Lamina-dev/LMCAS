#include "complex_analysis.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> cmplx(int a, int b) {
    return std::make_shared<SymbolicExpr>(
        SymbolicFactory::create_complex(num(a)->root, num(b)->root));
}

int main() {
    // ---- residue of 1/(z-a) at z=a is 1 ----
    {
        auto z = SymbolicExpr::variable("z");
        auto a = num(2);
        auto denom = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), a));
        auto f = SymbolicExpr::divide(num(1), denom);
        auto res = residue(f, "z", a, 1);
        EXPECT_TRUE(res != nullptr, "residue not null");
        EXPECT_EQ_EXPR(res->simplify(), num(1), "Res(1/(z-2), z=2) = 1");
    }

    // ---- real_part / imag_part of (3 + 4i) ----
    {
        auto z = cmplx(3, 4);
        EXPECT_EQ_EXPR(real_part(z), num(3), "Re(3+4i) = 3");
        EXPECT_EQ_EXPR(imag_part(z), num(4), "Im(3+4i) = 4");
    }

    // ---- real/imag of (2+i)*(1+i) = 2 + 2i + i + i^2 = 1 + 3i ----
    {
        auto prod = SymbolicExpr::multiply(cmplx(2,1), cmplx(1,1));
        EXPECT_EQ_EXPR(real_part(prod), num(1), "Re((2+i)(1+i)) = 1");
        EXPECT_EQ_EXPR(imag_part(prod), num(3), "Im((2+i)(1+i)) = 3");
    }

    // ---- conjugate of (3+4i) = 3-4i ----
    {
        auto z = cmplx(3, 4);
        auto c = conjugate(z);
        EXPECT_EQ_EXPR(real_part(c), num(3), "Re(conj(3+4i)) = 3");
        EXPECT_EQ_EXPR(imag_part(c)->simplify(), num(-4), "Im(conj(3+4i)) = -4");
    }

    // ---- i^2 = -1 ----
    {
        auto i = cmplx(0, 1);
        auto i2 = SymbolicExpr::multiply(i, i)->simplify();
        EXPECT_EQ_EXPR(real_part(i2), num(-1), "Re(i^2) = -1");
        EXPECT_EQ_EXPR(imag_part(i2)->simplify(), num(0), "Im(i^2) = 0");
    }

    // ---- is_analytic: f(z) = z^2 is analytic ----
    {
        auto z = SymbolicExpr::variable("z");
        auto f = SymbolicExpr::multiply(z, z);
        EXPECT_TRUE(is_analytic(f, "z"), "z^2 is analytic");
    }

    return TEST_REPORT();
}
