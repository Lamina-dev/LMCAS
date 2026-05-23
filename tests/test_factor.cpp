#include "test_common.hpp"

bool str_contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

int main() {
    try {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto two = SymbolicExpr::number(2);

        TEST_CASE("Factor Common Term");
        {
            auto t1 = SymbolicExpr::multiply(two, x);
            auto t2 = SymbolicExpr::multiply(two, y);
            auto expr1 = SymbolicExpr::add(t1, t2);

            auto factored1 = expr1->factor();
            std::cout << "  Factored: " << factored1->to_string() << std::endl;

            std::string s = factored1->to_string();

            bool ok = str_contains(s, "2*(") || str_contains(s, "2 (");
            EXPECT_TRUE(ok, "Factor 2x+2y");
        }

        TEST_CASE("Factor Quadratic Simple");
        {

            auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
            auto x5 = SymbolicExpr::multiply(SymbolicExpr::number(5), x);
            auto n6 = SymbolicExpr::number(6);

            auto expr2 = SymbolicExpr::add(x2, SymbolicExpr::add(x5, n6));

            auto factored2 = expr2->factor();
            std::cout << "  Factored: " << factored2->to_string() << std::endl;

            std::string s = factored2->to_string();
            bool ok = str_contains(s, "x + 2") && str_contains(s, "x + 3");
            EXPECT_TRUE(ok, "Factor x^2+5x+6");
        }

        TEST_CASE("Factor Difference of Squares");
        {

            auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
            auto n4 = SymbolicExpr::number(-4);
            auto expr3 = SymbolicExpr::add(x2, n4);

            auto factored3 = expr3->factor();
            std::cout << "  Factored: " << factored3->to_string() << std::endl;

            std::string s = factored3->to_string();
            bool ok = str_contains(s, "x + 2") && (str_contains(s, "x - 2") || str_contains(s, "x + -2"));
            EXPECT_TRUE(ok, "Factor x^2-4");
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
