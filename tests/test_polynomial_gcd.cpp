#include "polynomial.hpp"
#include "bigint.hpp"
#include <iostream>
#include <vector>

using namespace lamina;

void test_gcd_primitive() {
    std::cout << "Testing Polynomial GCD with BigInt..." << std::endl;

    std::vector<BigInt> c1 = {BigInt("-2"), BigInt("1"), BigInt("1")};
    std::vector<BigInt> c2 = {BigInt("-3"), BigInt("2"), BigInt("1")};

    Polynomial<BigInt> p1(c1, "x");
    Polynomial<BigInt> p2(c2, "x");

    Polynomial<BigInt> g = Polynomial<BigInt>::gcd(p1, p2);

    std::cout << "GCD((x-1)(x+2), (x-1)(x+3)) coefficients: ";
    for (const auto& c : g.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    std::vector<BigInt> ca = {BigInt("10"), BigInt("21"), BigInt("12"), BigInt("1")};
    std::vector<BigInt> cb = {BigInt("20"), BigInt("41"), BigInt("22"), BigInt("1")};

    Polynomial<BigInt> pa(ca, "x");
    Polynomial<BigInt> pb(cb, "x");

    Polynomial<BigInt> gab = Polynomial<BigInt>::gcd(pa, pb);

    std::cout << "GCD coefficients (Expected 1 2 1): ";
    for (const auto& c : gab.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    std::vector<BigInt> c3 = {BigInt("-2"), BigInt("2")};
    std::vector<BigInt> c4 = {BigInt("-6"), BigInt("6")};

    Polynomial<BigInt> p3(c3, "x");
    Polynomial<BigInt> p4(c4, "x");

    Polynomial<BigInt> g34 = Polynomial<BigInt>::gcd(p3, p4);

    std::cout << "GCD(2(x-1), 6(x-1)) coefficients (Expected -2 2): ";
    for (const auto& c : g34.coeffs) std::cout << c.ToString() << " ";
    std::cout << std::endl;

    std::cout << "\nTesting Polynomial GCD with Rational..." << std::endl;

    std::vector<Rational> r1 = {Rational(-1), Rational(0), Rational(1)};
    std::vector<Rational> r2 = {Rational(1), Rational(2), Rational(1)};

    Polynomial<Rational> pr1(r1, "x");
    Polynomial<Rational> pr2(r2, "x");

    Polynomial<Rational> gr = Polynomial<Rational>::gcd(pr1, pr2);

    std::cout << "GCD(x^2-1, x^2+2x+1) coefficients (Expected 1 1): ";
    for (const auto& c : gr.coeffs) std::cout << c.to_string() << " ";
    std::cout << std::endl;

    if (gr.degree() == 1 && gr.coeffs[0].to_double() == 1.0 && gr.coeffs[1].to_double() == 1.0) {
        std::cout << "Rational GCD Test Passed." << std::endl;
    } else {
        std::cout << "Rational GCD Test Failed." << std::endl;
    }

}

int main() {
    try {
        test_gcd_primitive();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
