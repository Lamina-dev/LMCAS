#include "polynomial.hpp"
#include "bigint.hpp"
#include "poly_utils.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"
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
    EXPECT_TRUE(g.degree() == 1, "BigInt GCD has degree 1 for shared x-1 factor");
    EXPECT_TRUE(gab.degree() == 2, "BigInt GCD has degree 2 for shared quadratic factor");
    EXPECT_TRUE(g34.degree() == 1, "BigInt GCD preserves shared linear primitive factor");
    EXPECT_TRUE(gr.degree() == 1 &&
                    gr.coeffs.size() >= 2 &&
                    gr.coeffs[0].to_string() == "1" &&
                    gr.coeffs[1].to_string() == "1",
                "Rational GCD finds x + 1");

}

void test_symbolic_polynomial_gcd() {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto one = SymbolicExpr::number(1);
    auto common = SymbolicExpr::add(x, y);
    auto lhs = SymbolicExpr::multiply(
        common, SymbolicExpr::add(x, one))->expand();
    auto rhs = SymbolicExpr::multiply(
        common, SymbolicExpr::add(y, one))->expand();

    ComputationContext context;
    auto result = symbolic_polynomial_gcd(*lhs, *rhs, context);
    EXPECT_TRUE(result.has_value(),
                "symbolic GCD accepts exact multivariate polynomials");
    if (result) {
        auto difference = SymbolicExpr::add(
            result.value(), SymbolicExpr::multiply(common, SymbolicExpr::number(-1)));
        EXPECT_TRUE(difference->expand()->simplify()->is_zero(),
                    "multivariate symbolic GCD recovers x + y");
    }

    auto quadratic_common = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(x, y),
                          SymbolicExpr::power(y, SymbolicExpr::number(2))));
    auto harder_lhs = SymbolicExpr::multiply(
        quadratic_common, SymbolicExpr::add(x, SymbolicExpr::number(2)))->expand();
    auto harder_rhs = SymbolicExpr::multiply(
        quadratic_common, SymbolicExpr::add(y, SymbolicExpr::number(3)))->expand();
    ComputationContext harder_context;
    auto harder_result = symbolic_polynomial_gcd(
        *harder_lhs, *harder_rhs, harder_context);
    EXPECT_TRUE(harder_result.has_value(),
                "symbolic GCD recovers a non-linear multivariate factor");
    if (harder_result) {
        auto difference = SymbolicExpr::add(
            harder_result.value(),
            SymbolicExpr::multiply(quadratic_common, SymbolicExpr::number(-1)));
        EXPECT_TRUE(difference->expand()->simplify()->is_zero(),
                    "non-linear multivariate symbolic GCD is maximal");
    }

    auto half = SymbolicExpr::number(Rational(1, 2));
    auto rational_common = SymbolicExpr::add(x, half);
    auto rational_lhs = SymbolicExpr::multiply(
        rational_common, SymbolicExpr::add(x, SymbolicExpr::number(2)));
    auto rational_rhs = SymbolicExpr::multiply(
        rational_common, SymbolicExpr::add(x, SymbolicExpr::number(3)));
    ComputationContext rational_context;
    auto rational_result = symbolic_polynomial_gcd(
        *rational_lhs, *rational_rhs, rational_context);
    EXPECT_TRUE(rational_result.has_value(),
                "symbolic GCD accepts exact rational coefficients");
    if (rational_result) {
        auto difference = SymbolicExpr::add(
            rational_result.value(),
            SymbolicExpr::multiply(rational_common, SymbolicExpr::number(-1)));
        EXPECT_TRUE(difference->expand()->simplify()->is_zero(),
                    "rational symbolic GCD is monic");
    }

    auto sine = SymbolicExpr::sin(x);
    ComputationContext unsupported_context;
    auto unsupported = symbolic_polynomial_gcd(*sine, *lhs, unsupported_context);
    EXPECT_TRUE(!unsupported &&
                    unsupported.error().code == CasErrc::UnsupportedExpression,
                "symbolic GCD rejects non-polynomial expressions");

    auto approximate = SymbolicExpr::add(
        x, SymbolicExpr::number(static_cast<lmmc_real_t>(0.5)));
    ComputationContext approximate_context;
    auto approximate_result = symbolic_polynomial_gcd(
        *approximate, *lhs, approximate_context);
    EXPECT_TRUE(!approximate_result &&
                    approximate_result.error().code ==
                        CasErrc::UnsupportedExpression,
                "symbolic GCD rejects approximate coefficients");

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited = symbolic_polynomial_gcd(*lhs, *rhs, limited_context);
    EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                "symbolic GCD observes the computation step budget");
}

int main() {
    try {
        test_gcd_primitive();
        test_symbolic_polynomial_gcd();
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
