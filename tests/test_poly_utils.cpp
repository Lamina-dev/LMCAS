#include "test_common.hpp"
#include "poly_utils.hpp"

void test_symbolic_to_poly() {
    TEST_CASE("symbolic_to_poly: linear 2x+1");
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)
        );
        auto poly = lamina::symbolic_to_poly<BigInt>(expr, "x");
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
        auto poly = lamina::symbolic_to_poly<BigInt>(expr, "x");
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
        auto poly = lamina::symbolic_to_poly<BigInt>(expr, "x");
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
        lamina::Polynomial<BigInt> poly({BigInt(1), BigInt(2)}, "x");
        auto expr = lamina::poly_to_symbolic<BigInt>(poly);
        std::string s = expr ? expr->to_string() : "null";
        // Should contain "2" and "x"
        EXPECT_CONTAINS(s, {"2", "x"}, "poly_to_symbolic linear contains 2 and x");
    }

    TEST_CASE("poly_to_symbolic: quadratic polynomial");
    {
        // x^2 - 3x + 2 => coeffs = {2, -3, 1}
        lamina::Polynomial<BigInt> poly({BigInt(2), BigInt(-3), BigInt(1)}, "x");
        auto expr = lamina::poly_to_symbolic<BigInt>(poly);
        std::string s = expr ? expr->to_string() : "null";
        // Should contain x^2 or x and structural elements
        EXPECT_CONTAINS(s, {"x", "2"}, "poly_to_symbolic quadratic contains x and 2");
    }

    TEST_CASE("poly_to_symbolic: constant polynomial");
    {
        // 5 => coeffs = {5}
        lamina::Polynomial<BigInt> poly({BigInt(5)}, "x");
        auto expr = lamina::poly_to_symbolic<BigInt>(poly);
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
        bool result = lamina::contains(*expr, "x");
        EXPECT_TRUE(result, "3x+1 depends on x");
    }

    TEST_CASE("depends_on_var: expression not containing variable");
    {
        auto y = SymbolicExpr::variable("y");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), y),
            SymbolicExpr::number(5)
        );
        bool result = lamina::contains(*expr, "x");
        EXPECT_TRUE(!result, "2y+5 does not depend on x");
    }

    TEST_CASE("depends_on_var: constant expression");
    {
        auto expr = SymbolicExpr::number(42);
        bool result = lamina::contains(*expr, "x");
        EXPECT_TRUE(!result, "constant 42 does not depend on x");
    }

    TEST_CASE("depends_on_var: nested expression containing variable");
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::power(
            SymbolicExpr::add(x, SymbolicExpr::number(1)),
            SymbolicExpr::number(2)
        );
        bool result = lamina::contains(*expr, "x");
        EXPECT_TRUE(result, "(x+1)^2 depends on x");
    }
}

void test_gaussian_eliminate() {
    TEST_CASE("gaussian_eliminate: 2x3 augmented matrix");
    {
        // System: x + 2y = 5, 3x + 4y = 11
        // Augmented matrix:
        // [1 2 5]
        // [3 4 11]
        // After elimination, should get row-echelon form with pivots at columns 0 and 1
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(5)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4), SymbolicExpr::number(11)}
        };

        std::vector<size_t> pivot_col_for_row;
        int sign;
        lamina::gaussian_eliminate(A, 2, 3, pivot_col_for_row, sign);

        // Verify pivot positions: row 0 has pivot at col 0, row 1 has pivot at col 1
        EXPECT_TRUE(pivot_col_for_row[0] == 0, "gaussian_eliminate pivot row 0 at col 0");
        EXPECT_TRUE(pivot_col_for_row[1] == 1, "gaussian_eliminate pivot row 1 at col 1");

        // After full elimination (reduced row echelon), A[0][0] should be 1 and A[1][1] should be 1
        std::string a00 = A[0][0] ? A[0][0]->simplify()->to_string() : "null";
        std::string a11 = A[1][1] ? A[1][1]->simplify()->to_string() : "null";
        EXPECT_CONTAINS(a00, {"1"}, "gaussian_eliminate A[0][0] is 1 (pivot)");
        EXPECT_CONTAINS(a11, {"1"}, "gaussian_eliminate A[1][1] is 1 (pivot)");

        // Solution: x=1, y=2 => A[0][2]=1, A[1][2]=2
        std::string a02 = A[0][2] ? A[0][2]->simplify()->to_string() : "null";
        std::string a12 = A[1][2] ? A[1][2]->simplify()->to_string() : "null";
        EXPECT_EQ_STR(a02, "1", "gaussian_eliminate solution x=1");
        EXPECT_EQ_STR(a12, "2", "gaussian_eliminate solution y=2");
    }
}

void test_gaussian_eliminate_singular() {
    TEST_CASE("gaussian_eliminate: singular matrix (rank deficient)");
    {
        // System: x + 2y = 3, 2x + 4y = 6 (second row is 2x first)
        // Augmented matrix:
        // [1 2 3]
        // [2 4 6]
        // After elimination, second row should become all zeros (rank 1)
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4), SymbolicExpr::number(6)}
        };

        std::vector<size_t> pivot_col_for_row;
        int sign;
        lamina::gaussian_eliminate(A, 2, 3, pivot_col_for_row, sign);

        // Row 0 should have a pivot, row 1 should have no pivot (SIZE_MAX)
        EXPECT_TRUE(pivot_col_for_row[0] == 0, "singular matrix pivot row 0 at col 0");
        EXPECT_TRUE(pivot_col_for_row[1] == (size_t)-1, "singular matrix row 1 has no pivot (rank deficient)");

        // Second row should be all zeros after elimination
        std::string a10 = A[1][0] ? A[1][0]->simplify()->to_string() : "null";
        std::string a11 = A[1][1] ? A[1][1]->simplify()->to_string() : "null";
        std::string a12 = A[1][2] ? A[1][2]->simplify()->to_string() : "null";
        EXPECT_EQ_STR(a10, "0", "singular matrix A[1][0] is 0");
        EXPECT_EQ_STR(a11, "0", "singular matrix A[1][1] is 0");
        EXPECT_EQ_STR(a12, "0", "singular matrix A[1][2] is 0");
    }
}

int main() {
    try {
        test_symbolic_to_poly();
        test_poly_to_symbolic();
        test_depends_on_var();
        test_gaussian_eliminate();
        test_gaussian_eliminate_singular();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
