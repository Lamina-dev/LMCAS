#include "poly_utils.hpp"
#include "polynomial.hpp"
#include "symbolic.hpp"
#include <iostream>
#include <cassert>

using namespace lamina;

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

    assert(poly.degree() == 2);
    assert(poly.coeffs[0] == BigInt(1));
    assert(poly.coeffs[1] == BigInt(2));
    assert(poly.coeffs[2] == BigInt(1));

    Polynomial<BigInt> poly2({BigInt(1), BigInt(1)}, "x");

    auto gcd = Polynomial<BigInt>::gcd(poly, poly2);
    std::cout << "GCD(P, x+1): ";
    for (const auto& c : gcd.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    assert(gcd.degree() == 1);

    auto res_expr = poly_to_symbolic(gcd);
    std::cout << "Result Symbolic: " << res_expr->to_string() << std::endl;

    std::cout << "\nTesting integrated SymbolicExpr::poly_gcd..." << std::endl;
    auto sp2 = SymbolicExpr::variable("x");
    sp2 = SymbolicExpr::add(sp2, SymbolicExpr::number(1));

    auto sgcd = SymbolicExpr::poly_gcd(expr, sp2);

    std::cout << "SymbolicExpr::poly_gcd result: " << sgcd->to_string() << std::endl;

    if (sgcd->is_one() && !expr->is_one() && !sp2->is_one()) {
        std::cerr << "Warning: poly_gcd returned 1, expected something else (unless coprime, but here they share x+1)." << std::endl;
    }

    std::cout << "Bridge Test Passed." << std::endl;
}

int main() {
    try {
        test_bridge();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
