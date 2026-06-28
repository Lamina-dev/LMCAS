#include "test_common.hpp"

int main() {
    try {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto one = SymbolicExpr::number(1);
        auto two = SymbolicExpr::number(2);
        auto three = SymbolicExpr::number(3);

        TEST_CASE("Cancel (x^2 - 1) / (x - 1) = x + 1");
        {
            // (x^2 - 1) / (x - 1)
            auto x2 = SymbolicExpr::power(x, two);
            auto numerator = SymbolicExpr::add(x2, SymbolicExpr::number(-1));
            auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 验证：代入 x=5 应得 6
            auto val = result->substitute("x", SymbolicExpr::number(5))->simplify();
            EXPECT_EQ_EXPR_STR(val, "6", "(x^2-1)/(x-1) at x=5 = 6");

            // 验证结果不含分母（无负指数）
            auto result_str = result->to_string();
            bool no_negative_power = result_str.find("^(-1)") == std::string::npos;
            EXPECT_TRUE(no_negative_power, "(x^2-1)/(x-1) cancels to polynomial");
        }

        TEST_CASE("Cancel (x^2 + 2x + 1) / (x + 1) = x + 1");
        {
            // (x+1)^2 / (x+1) = x+1
            auto x2 = SymbolicExpr::power(x, two);
            auto numerator = SymbolicExpr::add(x2, SymbolicExpr::add(
                SymbolicExpr::multiply(two, x), one));
            auto denominator = SymbolicExpr::add(x, one);
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=3 应得 4
            auto val = result->substitute("x", SymbolicExpr::number(3))->simplify();
            EXPECT_EQ_EXPR_STR(val, "4", "(x^2+2x+1)/(x+1) at x=3 = 4");
        }

        TEST_CASE("Cancel (x^3 - x) / (x^2 - 1) = x");
        {
            // x(x^2-1) / (x^2-1) = x
            auto x3 = SymbolicExpr::power(x, three);
            auto numerator = SymbolicExpr::add(x3, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
            auto denominator = SymbolicExpr::add(SymbolicExpr::power(x, two), SymbolicExpr::number(-1));
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=7 应得 7
            auto val = result->substitute("x", SymbolicExpr::number(7))->simplify();
            EXPECT_EQ_EXPR_STR(val, "7", "(x^3-x)/(x^2-1) at x=7 = 7");
        }

        TEST_CASE("Cancel with no common factor: (x + 1) / (x + 2)");
        {
            auto numerator = SymbolicExpr::add(x, one);
            auto denominator = SymbolicExpr::add(x, two);
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=3 应得 4/5 = 0.8，验证数值正确性
            // 由于无公因式，结果仍为分式形式
            auto val_num = result->substitute("x", SymbolicExpr::number(3))->simplify();
            std::cout << "  Value at x=3: " << val_num->to_string() << std::endl;
            // (3+1)/(3+2) = 4/5
            auto expected = SymbolicExpr::divide(SymbolicExpr::number(4), SymbolicExpr::number(5))->simplify();
            auto diff = SymbolicExpr::add(val_num, SymbolicExpr::multiply(SymbolicExpr::number(-1), expected))->simplify();
            EXPECT_TRUE(diff->is_zero(), "(x+1)/(x+2) at x=3 = 4/5");
        }

        TEST_CASE("Cancel constant: 6/3 = 2");
        {
            auto expr = SymbolicExpr::divide(SymbolicExpr::number(6), three);
            auto result = expr->cancel();
            std::cout << "  Result: " << result->to_string() << std::endl;
            EXPECT_EQ_EXPR_STR(result, "2", "6/3 = 2");
        }

        TEST_CASE("Cancel multivariate: (x*y - y) / (x - 1) = y");
        {
            // (xy - y) / (x - 1) = y(x-1)/(x-1) = y
            auto xy = SymbolicExpr::multiply(x, y);
            auto numerator = SymbolicExpr::add(xy, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
            auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=3, y=5 应得 5
            auto val = result->substitute("x", SymbolicExpr::number(3))
                             ->substitute("y", SymbolicExpr::number(5))->simplify();
            EXPECT_EQ_EXPR_STR(val, "5", "(xy-y)/(x-1) at x=3,y=5 = 5");
        }

        TEST_CASE("Cancel polynomial / polynomial: already simplified");
        {
            // x / 1 = x (no denominator)
            auto result = x->cancel();
            std::cout << "  Result: " << result->to_string() << std::endl;
            EXPECT_EQ_EXPR_STR(result, "x", "x cancels to x");
        }

        TEST_CASE("Cancel (2x^2 + 2x) / (2x) = x + 1");
        {
            auto x2 = SymbolicExpr::power(x, two);
            auto numerator = SymbolicExpr::add(
                SymbolicExpr::multiply(two, x2),
                SymbolicExpr::multiply(two, x));
            auto denominator = SymbolicExpr::multiply(two, x);
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=4 应得 5
            auto val = result->substitute("x", SymbolicExpr::number(4))->simplify();
            EXPECT_EQ_EXPR_STR(val, "5", "(2x^2+2x)/(2x) at x=4 = 5");
        }

        TEST_CASE("Factor cubic: x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)");
        {
            // x^3 - 6x^2 + 11x - 6
            auto x3 = SymbolicExpr::power(x, three);
            auto x2 = SymbolicExpr::power(x, two);
            auto expr = SymbolicExpr::add(x3,
                SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-6), x2),
                SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(11), x),
                    SymbolicExpr::number(-6))));

            auto factored = expr->factor();
            std::cout << "  Factored: " << factored->to_string() << std::endl;

            // 验证：代入 x=1, x=2, x=3 应得 0
            auto v1 = factored->substitute("x", SymbolicExpr::number(1))->simplify();
            auto v2 = factored->substitute("x", SymbolicExpr::number(2))->simplify();
            auto v3 = factored->substitute("x", SymbolicExpr::number(3))->simplify();
            EXPECT_TRUE(v1->is_zero(), "x^3-6x^2+11x-6 at x=1 = 0");
            EXPECT_TRUE(v2->is_zero(), "x^3-6x^2+11x-6 at x=2 = 0");
            EXPECT_TRUE(v3->is_zero(), "x^3-6x^2+11x-6 at x=3 = 0");

            // 验证是乘积形式
            auto s = factored->to_string();
            bool is_product = s.find("*") != std::string::npos;
            EXPECT_TRUE(is_product, "cubic factored into product form");
        }

        TEST_CASE("Cancel higher degree: (x^3 - x) / (x^2 + x) = x - 1");
        {
            // (x^3 - x) / (x^2 + x) = x(x^2-1) / x(x+1) = (x-1)(x+1)/(x+1) = x-1
            auto x3 = SymbolicExpr::power(x, three);
            auto x2 = SymbolicExpr::power(x, two);
            auto numerator = SymbolicExpr::add(x3, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
            auto denominator = SymbolicExpr::add(x2, x);
            auto expr = SymbolicExpr::divide(numerator, denominator);

            auto result = expr->cancel();
            std::cout << "  Input:  " << expr->simplify()->to_string() << std::endl;
            std::cout << "  Result: " << result->to_string() << std::endl;

            // 代入 x=5 应得 4
            auto val = result->substitute("x", SymbolicExpr::number(5))->simplify();
            EXPECT_EQ_EXPR_STR(val, "4", "(x^3-x)/(x^2+x) at x=5 = 4");
        }

    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception" << std::endl;
        g_failures++;
    }

    return TEST_REPORT();
}
