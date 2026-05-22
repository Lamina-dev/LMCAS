#include "symbolic.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>

template<typename T>
std::shared_ptr<SymbolicExpr> num(T n) {
    return SymbolicExpr::number(n);
}

std::shared_ptr<SymbolicExpr> var(const std::string& name) {
    return SymbolicExpr::variable(name);
}

void test_fraction_arithmetic() {
    std::cout << "--- Fraction Arithmetic ---" << std::endl;

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto sum = SymbolicExpr::add(half, third)->simplify();

    std::cout << "1/2 + 1/3 = " << sum->to_string() << std::endl;
}

void test_negative_power() {
    std::cout << "--- Negative Power ---" << std::endl;

    auto base = num(2);
    auto exponent = num(-2);
    auto res = SymbolicExpr::power(base, exponent)->simplify();

    std::cout << "pow(2, -2) = " << res->to_string() << std::endl;

    auto base3 = num(3);
    auto exp_neg1 = num(-1);
    auto res2 = SymbolicExpr::power(base3, exp_neg1)->simplify();
    std::cout << "pow(3, -1) = " << res2->to_string() << std::endl;
}

void test_fraction_mixed_with_var() {
    std::cout << "--- Mixed with Var ---" << std::endl;

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto x = var("x");

    auto term1 = SymbolicExpr::multiply(half, x);
    auto term2 = SymbolicExpr::multiply(third, x);
    auto res = SymbolicExpr::add(term1, term2)->simplify();

    std::cout << "(1/2)x + (1/3)x = " << res->to_string() << std::endl;
}

void test_rational_simplification() {
    std::cout << "--- Rational Simplification ---" << std::endl;

    auto six = num(6);
    auto twelve = num(12);
    auto res = SymbolicExpr::divide(six, twelve)->simplify();
    std::cout << "6/12 = " << res->to_string() << std::endl;

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto sq = SymbolicExpr::power(half, num(2))->simplify();
    std::cout << "(1/2)^2 = " << sq->to_string() << std::endl;
}

void test_fraction_matrix() {
    std::cout << "--- Fraction Matrix ---" << std::endl;

    auto m11 = SymbolicExpr::divide(num(1), num(2));
    auto m12 = SymbolicExpr::divide(num(1), num(3));
    auto m21 = SymbolicExpr::divide(num(1), num(4));
    auto m22 = SymbolicExpr::divide(num(1), num(5));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> elements = {
        {m11, m12},
        {m21, m22}
    };

    auto mat = SymbolicExpr::matrix(elements);
    auto det = SymbolicExpr::determinant(mat)->simplify();

    std::cout << "det([[1/2, 1/3], [1/4, 1/5]]) = " << det->to_string() << std::endl;

    auto i1 = num(1); auto i2 = num(2);
    auto i3 = num(3); auto i4 = num(4);
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> inv_elems = {{i1, i2}, {i3, i4}};
    auto mat_inv = SymbolicExpr::matrix(inv_elems);
    auto inv = SymbolicExpr::inverse(mat_inv)->simplify();

    std::cout << "inv([[1, 2], [3, 4]]) = " << inv->to_string() << std::endl;
}

int main() {
    try {
        test_fraction_arithmetic();
        test_negative_power();
        test_fraction_mixed_with_var();
        test_rational_simplification();
        test_fraction_matrix();
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
