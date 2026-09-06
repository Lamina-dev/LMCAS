#include "test_common.hpp"
#include "poly_utils.hpp"

using namespace LMCAS;

void test_symbolic_to_poly() {
    TEST_CASE("symbolic_to_poly: linear 2x+1");
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)
        );
        auto poly = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        // coeffs[0] = 1, coeffs[1] = 2
        EXPECT_TRUE(poly.degree() == 1, "2x+1 degree is 1");
        EXPECT_TRUE(poly.coeffs.size() == 2, "2x+1 has 2 coefficients");
        EXPECT_TRUE(poly.coeffs[0] == BigInt(1), "2x+1 constant term is 1");
        EXPECT_TRUE(poly.coeffs[1] == BigInt(2), "2x+1 linear term is 2");
    }

    TEST_CASE("symbolic_to_poly: quadratic x^2-3x+2");
    {
        auto x = SymbolicExpr::variable("x");
        // x^2 - 3x + 2
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto neg3x = SymbolicExpr::multiply(SymbolicExpr::number(-3), SymbolicExpr::variable("x"));
        auto two = SymbolicExpr::number(2);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(x2, neg3x), two);
        auto poly = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        // coeffs[0] = 2, coeffs[1] = -3, coeffs[2] = 1
        EXPECT_TRUE(poly.degree() == 2, "x^2-3x+2 degree is 2");
        EXPECT_TRUE(poly.coeffs[0] == BigInt(2), "x^2-3x+2 constant term is 2");
        EXPECT_TRUE(poly.coeffs[1] == BigInt(-3), "x^2-3x+2 linear term is -3");
        EXPECT_TRUE(poly.coeffs[2] == BigInt(1), "x^2-3x+2 quadratic term is 1");
    }

    TEST_CASE("symbolic_to_poly: cubic x^3+x");
    {
        auto x = SymbolicExpr::variable("x");
        // x^3 + x
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto expr = SymbolicExpr::add(x3, SymbolicExpr::variable("x"));
        auto poly = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        // coeffs[0] = 0, coeffs[1] = 1, coeffs[2] = 0, coeffs[3] = 1
        EXPECT_TRUE(poly.degree() == 3, "x^3+x degree is 3");
        EXPECT_TRUE(poly.coeffs[0] == BigInt(0), "x^3+x constant term is 0");
        EXPECT_TRUE(poly.coeffs[1] == BigInt(1), "x^3+x linear term is 1");
        EXPECT_TRUE(poly.coeffs[2] == BigInt(0), "x^3+x quadratic term is 0");
        EXPECT_TRUE(poly.coeffs[3] == BigInt(1), "x^3+x cubic term is 1");
    }
}

void test_poly_to_symbolic() {
    TEST_CASE("poly_to_symbolic: linear polynomial");
    {
        // 2x + 1 => coeffs = {1, 2}
        LMCAS::Polynomial<BigInt> poly({BigInt(1), BigInt(2)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        std::string s = expr ? expr->to_string() : "null";
        // Should contain "2" and "x"
        EXPECT_CONTAINS(s, {"2", "x"}, "poly_to_symbolic linear contains 2 and x");
    }

    TEST_CASE("poly_to_symbolic: quadratic polynomial");
    {
        // x^2 - 3x + 2 => coeffs = {2, -3, 1}
        LMCAS::Polynomial<BigInt> poly({BigInt(2), BigInt(-3), BigInt(1)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        std::string s = expr ? expr->to_string() : "null";
        // Should contain x^2 or x and structural elements
        EXPECT_CONTAINS(s, {"x", "2"}, "poly_to_symbolic quadratic contains x and 2");
    }

    TEST_CASE("poly_to_symbolic: constant polynomial");
    {
        // 5 => coeffs = {5}
        LMCAS::Polynomial<BigInt> poly({BigInt(5)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        std::string s = expr ? expr->to_string() : "null";
        EXPECT_EQ_STR(s, "5", "poly_to_symbolic constant is 5");
    }
}

void test_depends_on_var() {
    TEST_CASE("depends_on_var: expression containing variable");
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(1)
        );
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE(result, "3x+1 depends on x");
    }

    TEST_CASE("depends_on_var: expression not containing variable");
    {
        auto y = SymbolicExpr::variable("y");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), y),
            SymbolicExpr::number(5)
        );
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE(!result, "2y+5 does not depend on x");
    }

    TEST_CASE("depends_on_var: constant expression");
    {
        auto expr = SymbolicExpr::number(42);
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE(!result, "constant 42 does not depend on x");
    }

    TEST_CASE("depends_on_var: nested expression containing variable");
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::power(
            SymbolicExpr::add(x, SymbolicExpr::number(1)),
            SymbolicExpr::number(2)
        );
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE(result, "(x+1)^2 depends on x");
    }
}


int main() {
    try {
        test_symbolic_to_poly();
        test_poly_to_symbolic();
        test_depends_on_var();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
