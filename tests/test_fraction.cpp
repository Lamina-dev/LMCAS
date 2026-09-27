#include "test_common.hpp"
#include "symbolic.hpp"
#include "symbolic_matrix.hpp"
#include <cmath>
#include <memory>

using namespace LMCAS;

template <typename T>
std::shared_ptr<SymbolicExpr> num(T n) {
    return SymbolicExpr::number(n);
}

std::shared_ptr<SymbolicExpr> var(const std::string &name) {
    return SymbolicExpr::variable(name);
}

TEST(Fraction, FractionArithmetic) {
    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto sum = SymbolicExpr::add(half, third)->simplify();

    EXPECT_TRUE(test_proved_equivalent(
        sum, SymbolicExpr::number(Rational(5, 6))));
}

TEST(Fraction, NegativePower) {
    auto base = num(2);
    auto exponent = num(-2);
    auto res = SymbolicExpr::power(base, exponent)->simplify();

    EXPECT_TRUE(test_proved_equivalent(
        res, SymbolicExpr::number(Rational(1, 4))));

    auto res2 = SymbolicExpr::power(num(3), num(-1))->simplify();
    EXPECT_TRUE(test_proved_equivalent(
        res2, SymbolicExpr::number(Rational(1, 3))));
}

TEST(Fraction, FractionMixedWithVar) {
    auto half = SymbolicExpr::divide(num(1), num(2));
    auto third = SymbolicExpr::divide(num(1), num(3));
    auto x = var("x");

    auto term1 = SymbolicExpr::multiply(half, x);
    auto term2 = SymbolicExpr::multiply(third, x);
    auto res = SymbolicExpr::add(term1, term2)->simplify();

    EXPECT_TRUE(test_proved_equivalent(
        res,
        SymbolicExpr::multiply(
            SymbolicExpr::number(Rational(5, 6)), x)));
}

TEST(Fraction, RationalSimplification) {
    auto six = num(6);
    auto twelve = num(12);
    auto res = SymbolicExpr::divide(six, twelve)->simplify();
    EXPECT_TRUE(test_proved_equivalent(
        res, SymbolicExpr::number(Rational(1, 2))));

    auto half = SymbolicExpr::divide(num(1), num(2));
    auto sq = SymbolicExpr::power(half, num(2))->simplify();
    EXPECT_TRUE(test_proved_equivalent(
        sq, SymbolicExpr::number(Rational(1, 4))));
}

TEST(Fraction, FractionMatrix) {
    auto m11 = SymbolicExpr::divide(num(1), num(2));
    auto m12 = SymbolicExpr::divide(num(1), num(3));
    auto m21 = SymbolicExpr::divide(num(1), num(4));
    auto m22 = SymbolicExpr::divide(num(1), num(5));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> elements = {
        {m11, m12},
        {m21, m22}};
    auto mat = SymbolicExpr::matrix(elements);

    auto determinant = LMCAS::matrix_determinant_checked(mat);
    ASSERT_TRUE(determinant) << determinant.error().message;
    auto det = determinant.value()->simplify();
    EXPECT_TRUE(test_proved_equivalent(
        det, SymbolicExpr::number(Rational(1, 60))));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> inv_elems = {
        {num(1), num(2)}, {num(3), num(4)}};
    auto mat_inv = SymbolicExpr::matrix(inv_elems);
    auto inverse = LMCAS::matrix_inverse_checked(mat_inv);
    ASSERT_TRUE(inverse) << inverse.error().message;
    auto simplified = inverse.value()->simplify();
    ASSERT_NE(simplified, nullptr);
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(simplified));
    ASSERT_NE(matrix, nullptr);
    ASSERT_EQ(matrix->rows(), 2U);
    ASSERT_EQ(matrix->cols(), 2U);
    const std::vector<Rational> expected = {
        Rational(-2), Rational(1),
        Rational(3, 2), Rational(-1, 2)};
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t column = 0; column < 2; ++column) {
            auto element =
                detail::make_expression_ptr(matrix->get(row, column));
            EXPECT_TRUE(test_proved_equivalent(
                element, SymbolicExpr::number(
                             expected[row * 2 + column])));
        }
    }
}
