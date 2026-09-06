#include "../include/symbolic.hpp"
#include <iostream>
#include <cassert>
#include "test_common.hpp"

using namespace LMCAS;

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto integrated = expr->integrate("x");
    EXPECT_TRUE(integrated != nullptr, "integral of x^2 exists");
    std::cout << "int(x^2) = " << integrated->to_string() << "\n";

    auto point = SymbolicExpr::number(2);
    auto lim = LMCAS::limit_expression_checked(expr, "x", point).value();
    EXPECT_TRUE(lim != nullptr && std::abs(lim->to_numeric() - 4.0) < 1e-12,
                "limit x^2 at 2 equals 4");
    std::cout << "limit(x^2, x->2) = " << lim->to_string() << "\n";

    auto five = SymbolicExpr::number(5);
    auto sum = SymbolicExpr::add(five, x);
    auto simple_sum = sum->simplify();
    std::cout << "5 + x simplified -> " << simple_sum->to_string() << "\n";

    std::cout << "--- LHopital Check ---\n";
    auto four = SymbolicExpr::number(4);
    auto two = SymbolicExpr::number(2);
    auto neg_four = SymbolicExpr::multiply(four, SymbolicExpr::number(-1));
    auto neg_two = SymbolicExpr::multiply(two, SymbolicExpr::number(-1));

    auto num = SymbolicExpr::add(expr, neg_four);
    auto den = SymbolicExpr::add(x, neg_two);

    std::cout << "Num: " << num->to_string() << "\n";
    std::cout << "Den: " << den->to_string() << "\n";

    auto dN = num->differentiate("x");
    auto dD = den->differentiate("x");
    std::cout << "d(Num)/dx: " << dN->to_string() << "\n";
    std::cout << "d(Den)/dx: " << dD->to_string() << "\n";

    auto minus_one = SymbolicExpr::number(-1);
    auto dD_inv = SymbolicExpr::power(dD, minus_one);
    std::cout << "dD^(-1): " << dD_inv->to_string() << "\n";

    auto ratio = SymbolicExpr::multiply(dN, dD_inv);
    std::cout << "Ratio dN/dD: " << ratio->to_string() << "\n";

    auto lim_ratio = LMCAS::limit_expression_checked(ratio, "x", two).value();
    std::cout << "Limit of Ratio: " << lim_ratio->to_string() << "\n";

    auto rational = SymbolicExpr::multiply(num, SymbolicExpr::power(den, minus_one));
    std::cout << "Rational: " << rational->to_string() << "\n";
    auto lim_rational = LMCAS::limit_expression_checked(rational, "x", two).value();
    EXPECT_TRUE(lim_rational != nullptr && std::abs(lim_rational->to_numeric() - 4.0) < 1e-12,
                "removable singularity limit equals 4");
    std::cout << "Result: " << lim_rational->to_string() << "\n";

    return TEST_REPORT();
}
