/**
 * @file test_limits_at_infinity.cpp
 * @brief 测试极限在无穷处的增强功能：有理函数次数比较、增长速率比较、
 *        复合函数、负无穷代换。
 */
#include "test_common.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto inf = SymbolicExpr::infinity(1);
    auto neg_inf = SymbolicExpr::infinity(-1);

    // =========================================================================
    // Requirement 5.1: Rational function degree comparison at infinity
    // =========================================================================
    TEST_CASE("Rational function: deg(P) < deg(Q) -> 0");
    {
        // lim(x→∞) x / x^2 = lim(x→∞) 1/x = 0
        auto num = x;
        auto den = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->inf) x/x^2 = 0");
    }

    TEST_CASE("Rational function: deg(P) = deg(Q) -> ratio of leading coefficients");
    {
        // lim(x→∞) (3x^2 + x) / (2x^2 + 1) = 3/2
        auto three = SymbolicExpr::number(3);
        auto two = SymbolicExpr::number(2);
        auto one = SymbolicExpr::number(1);
        auto x2 = SymbolicExpr::power(x, two);
        auto num = SymbolicExpr::add(SymbolicExpr::multiply(three, x2), x);
        auto den = SymbolicExpr::add(SymbolicExpr::multiply(two, x2), one);
        auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "3/2", "lim(x->inf) (3x^2+x)/(2x^2+1) = 3/2");
    }

    TEST_CASE("Rational function: deg(P) > deg(Q) -> infinity");
    {
        // lim(x→∞) x^3 / x = lim(x→∞) x^2 = ∞
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto expr = SymbolicExpr::multiply(x3, SymbolicExpr::power(x, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        // Should be infinity
        EXPECT_TRUE(result != nullptr, "lim(x->inf) x^3/x is not null");
        if (result) {
            auto str = result->to_string();
            EXPECT_TRUE(str.find("inf") != std::string::npos || str.find("Inf") != std::string::npos || str.find("∞") != std::string::npos,
                "lim(x->inf) x^3/x = infinity");
        }
    }

    TEST_CASE("Rational function: equal degree with negative leading coeff");
    {
        // lim(x→∞) (-2x^2) / (x^2) = -2
        auto neg_two = SymbolicExpr::number(-2);
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto num = SymbolicExpr::multiply(neg_two, x2);
        auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(x2, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "-2", "lim(x->inf) -2x^2/x^2 = -2");
    }

    // =========================================================================
    // Requirement 5.2: Growth-rate comparison
    // =========================================================================
    TEST_CASE("Growth rate: exp(x) dominates polynomial");
    {
        // lim(x→∞) x^2 / exp(x) = 0
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto exp_x = SymbolicExpr::exp(x);
        auto expr = SymbolicExpr::multiply(x2, SymbolicExpr::power(exp_x, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->inf) x^2/exp(x) = 0");
    }

    TEST_CASE("Growth rate: polynomial dominates logarithmic");
    {
        // lim(x→∞) ln(x) / x = 0
        auto ln_x = SymbolicExpr::ln(x);
        auto expr = SymbolicExpr::multiply(ln_x, SymbolicExpr::power(x, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->inf) ln(x)/x = 0");
    }

    TEST_CASE("Growth rate: exp(x) dominates x^n for large n");
    {
        // lim(x→∞) x^10 / exp(x) = 0
        auto x10 = SymbolicExpr::power(x, SymbolicExpr::number(10));
        auto exp_x = SymbolicExpr::exp(x);
        auto expr = SymbolicExpr::multiply(x10, SymbolicExpr::power(exp_x, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->inf) x^10/exp(x) = 0");
    }

    // =========================================================================
    // Requirement 5.3: Composed functions
    // =========================================================================
    TEST_CASE("Composed function: exp(-x) as x->inf");
    {
        // lim(x→∞) exp(-x) = 0
        auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
        auto expr = SymbolicExpr::exp(neg_x);
        auto result = expr->limit("x", inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->inf) exp(-x) = 0");
    }

    // =========================================================================
    // Requirement 5.4: Limits at negative infinity
    // =========================================================================
    TEST_CASE("Negative infinity: lim(x->-inf) x^2 / (x^2 + 1) = 1");
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto den = SymbolicExpr::add(x2, SymbolicExpr::number(1));
        auto expr = SymbolicExpr::multiply(x2, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
        auto result = expr->limit("x", neg_inf);
        EXPECT_EQ_EXPR_STR(result, "1", "lim(x->-inf) x^2/(x^2+1) = 1");
    }

    TEST_CASE("Negative infinity: lim(x->-inf) exp(x) = 0");
    {
        auto expr = SymbolicExpr::exp(x);
        auto result = expr->limit("x", neg_inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->-inf) exp(x) = 0");
    }

    TEST_CASE("Negative infinity: lim(x->-inf) 1/x = 0");
    {
        auto expr = SymbolicExpr::power(x, SymbolicExpr::number(-1));
        auto result = expr->limit("x", neg_inf);
        EXPECT_EQ_EXPR_STR(result, "0", "lim(x->-inf) 1/x = 0");
    }

    return TEST_REPORT();
}
