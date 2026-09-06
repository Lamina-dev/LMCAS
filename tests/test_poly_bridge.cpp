#include "poly_utils.hpp"
#include "polynomial.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"
#include <iostream>

using namespace LMCAS;

void test_bridge() {
    std::cout << "Testing Symbolic <-> Polynomial Bridge..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto one = SymbolicExpr::number(1);

    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, two),
        SymbolicExpr::add(
            SymbolicExpr::multiply(two, x),
            one
        )
    );

    std::cout << "Symbolic Expr: " << expr->to_string() << std::endl;

    Polynomial<BigInt> poly = symbolic_to_poly<BigInt>(expr, "x");

    std::cout << "Polynomial coeffs: ";
    for (const auto& c : poly.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    EXPECT_TRUE(poly.degree() == 2, "symbolic_to_poly preserves degree");
    EXPECT_TRUE(poly.coeffs.size() >= 3, "symbolic_to_poly returns three coefficients");
    if (poly.coeffs.size() >= 3) {
        EXPECT_TRUE(poly.coeffs[0] == BigInt(1), "constant coefficient is 1");
        EXPECT_TRUE(poly.coeffs[1] == BigInt(2), "linear coefficient is 2");
        EXPECT_TRUE(poly.coeffs[2] == BigInt(1), "quadratic coefficient is 1");
    }

    Polynomial<BigInt> poly2({BigInt(1), BigInt(1)}, "x");

    auto gcd = Polynomial<BigInt>::gcd(poly, poly2);
    std::cout << "GCD(P, x+1): ";
    for (const auto& c : gcd.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    EXPECT_TRUE(gcd.degree() == 1, "Polynomial::gcd finds shared linear factor");

    auto res_expr = poly_to_symbolic(gcd);
    std::cout << "Result Symbolic: " << res_expr->to_string() << std::endl;

    std::cout << "\nTesting integrated SymbolicExpr::poly_gcd..." << std::endl;
    auto sp2 = SymbolicExpr::variable("x");
    sp2 = SymbolicExpr::add(sp2, SymbolicExpr::number(1));

    ComputationContext gcd_context;
    auto sgcd = symbolic_polynomial_gcd(
        *expr, *sp2, gcd_context).value();

    std::cout << "SymbolicExpr::poly_gcd result: " << sgcd->to_string() << std::endl;

    EXPECT_TRUE(!(sgcd->is_one() && !expr->is_one() && !sp2->is_one()),
                "SymbolicExpr::poly_gcd does not silently return 1 for expressions sharing x+1");

    std::cout << "Bridge Test Passed." << std::endl;
}

int main() {
    try {
        test_bridge();
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
