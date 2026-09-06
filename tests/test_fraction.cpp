#include "test_common.hpp"
#include "symbolic.hpp"
#include "symbolic_matrix.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>

using namespace LMCAS;

template<typename T>
std::shared_ptr<SymbolicExpr> num(T n) {
    return SymbolicExpr::number(n);
}

std::shared_ptr<SymbolicExpr> var(const std::string& name) {
    return SymbolicExpr::variable(name);
}

void test_fraction_arithmetic() {
    TEST_CASE("Fraction arithmetic");

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto sum = SymbolicExpr::add(half, third)->simplify();

    EXPECT_EQ_EXPR_STR(sum, "5/6", "1/2 + 1/3 = 5/6");
}

void test_negative_power() {
    TEST_CASE("Negative powers");

    auto base = num(2);
    auto exponent = num(-2);
    auto res = SymbolicExpr::power(base, exponent)->simplify();

    EXPECT_EQ_EXPR_STR(res, "1/4", "2^-2 = 1/4");

    auto base3 = num(3);
    auto exp_neg1 = num(-1);
    auto res2 = SymbolicExpr::power(base3, exp_neg1)->simplify();
    EXPECT_EQ_EXPR_STR(res2, "1/3", "3^-1 = 1/3");
}

void test_fraction_mixed_with_var() {
    TEST_CASE("Fraction coefficients with variables");

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto x = var("x");

    auto term1 = SymbolicExpr::multiply(half, x);
    auto term2 = SymbolicExpr::multiply(third, x);
    auto res = SymbolicExpr::add(term1, term2)->simplify();

    EXPECT_EQ_EXPR_STR(res, "(5/6)*x", "(1/2)x + (1/3)x = (5/6)x");
}

void test_rational_simplification() {
    TEST_CASE("Rational simplification");

    auto six = num(6);
    auto twelve = num(12);
    auto res = SymbolicExpr::divide(six, twelve)->simplify();
    EXPECT_EQ_EXPR_STR(res, "1/2", "6/12 simplifies to 1/2");

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto sq = SymbolicExpr::power(half, num(2))->simplify();
    EXPECT_EQ_EXPR_STR(sq, "1/4", "(1/2)^2 = 1/4");
}

void test_fraction_matrix() {
    TEST_CASE("Fraction matrix operations");

    auto m11 = SymbolicExpr::divide(num(1), num(2));
    auto m12 = SymbolicExpr::divide(num(1), num(3));
    auto m21 = SymbolicExpr::divide(num(1), num(4));
    auto m22 = SymbolicExpr::divide(num(1), num(5));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> elements = {
        {m11, m12},
        {m21, m22}
    };

    auto mat = SymbolicExpr::matrix(elements);
    auto det = LMCAS::matrix_determinant_checked(mat).value()->simplify();

    EXPECT_EQ_EXPR_STR(det, "1/60", "det([[1/2, 1/3], [1/4, 1/5]]) = 1/60");

    auto i1 = num(1); auto i2 = num(2);
    auto i3 = num(3); auto i4 = num(4);
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> inv_elems = {{i1, i2}, {i3, i4}};
    auto mat_inv = SymbolicExpr::matrix(inv_elems);
    auto inv = LMCAS::matrix_inverse_checked(mat_inv).value()->simplify();

    EXPECT_EQ_EXPR_STR(inv, "[[-2, 1], [3/2, -1/2]]", "inverse of [[1, 2], [3, 4]] is exact");
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
    return TEST_REPORT();
}
