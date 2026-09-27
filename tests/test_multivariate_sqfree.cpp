#include "test_multivariate_sqfree_support.hpp"

TEST(MultivariateSqfree, AlreadySquareFreePolynomialX2Y2) {
    // x^2 - y^2 = (x+y)(x-y), already square-free
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), // x^2
        make_term({0, 2}, Rational(-1)) // -y^2
    };
    MultiPoly f(terms, vars);

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "x^2 - y^2: product matches and components are square-free";

    // Should have single component (already square-free)
    EXPECT_TRUE((decomp.components.size() == 1)) << "x^2 - y^2: single component (already square-free)";
}

TEST(MultivariateSqfree, PolynomialWithSquaredFactorXY2XY) {
    // f = (x+y)^2 * (x-y) = (x^2 + 2xy + y^2)(x - y)
    //   = x^3 + 2x^2*y + x*y^2 - x^2*y - 2x*y^2 - y^3
    //   = x^3 + x^2*y - x*y^2 - y^3
    std::vector<std::string> vars = {"x", "y"};

    // Build (x+y)
    std::vector<MultiPoly::Term> xpy_terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(1))  // y
    };
    MultiPoly xpy(xpy_terms, vars);

    // Build (x-y)
    std::vector<MultiPoly::Term> xmy_terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(-1)) // -y
    };
    MultiPoly xmy(xmy_terms, vars);

    // f = (x+y)^2 * (x-y)
    MultiPoly f = xpy * xpy * xmy;

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "(x+y)^2*(x-y): product matches and components are square-free";
}

TEST(MultivariateSqfree, PolynomialWithCubedFactorX13) {
    // f = (x+1)^3 = x^3 + 3x^2 + 3x + 1
    std::vector<std::string> vars = {"x", "y"};

    // Build (x+1) in {x, y} variable set
    std::vector<MultiPoly::Term> xp1_terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 0}, Rational(1))  // 1
    };
    MultiPoly xp1(xp1_terms, vars);

    // f = (x+1)^3
    MultiPoly f = xp1 * xp1 * xp1;

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "(x+1)^3: product matches and components are square-free";
}

TEST(MultivariateSqfree, UnivariateWithRepeatedRootsX12X1) {
    // f = (x-1)^2 * (x+1) = (x^2 - 2x + 1)(x + 1)
    //   = x^3 - x^2 - x + 1
    std::vector<std::string> vars = {"x"};

    // Build (x-1)
    std::vector<MultiPoly::Term> xm1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly xm1(xm1_terms, vars);

    // Build (x+1)
    std::vector<MultiPoly::Term> xp1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // 1
    };
    MultiPoly xp1(xp1_terms, vars);

    // f = (x-1)^2 * (x+1)
    MultiPoly f = xm1 * xm1 * xp1;

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "(x-1)^2*(x+1): product matches and components are square-free";
}

TEST(MultivariateSqfree, ConstantPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly f(Rational(42), vars);

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    // For constant, decomposition should return the constant itself
    EXPECT_TRUE((!decomp.components.empty())) << "constant 42: decomposition is non-empty";
    EXPECT_TRUE((decomp.components[0].is_constant())) << "constant 42: first component is constant";
}

TEST(MultivariateSqfree, LinearPolynomialX2y3) {
    // Linear in x → already square-free
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(2)), // 2y
        make_term({0, 0}, Rational(3))  // 3
    };
    MultiPoly f(terms, vars);

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "x + 2y + 3: product matches and components are square-free";
    EXPECT_TRUE((decomp.components.size() == 1)) << "x + 2y + 3: single component (linear, already square-free)";
}

TEST(MultivariateSqfree, PolynomialWithSquaredAndCubedFactorsX12X13) {
    // f = (x+1)^2 * (x-1)^3
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> xp1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // 1
    };
    MultiPoly xp1(xp1_terms, vars);

    std::vector<MultiPoly::Term> xm1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly xm1(xm1_terms, vars);

    // f = (x+1)^2 * (x-1)^3
    MultiPoly f = xp1 * xp1 * xm1 * xm1 * xm1;

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "(x+1)^2*(x-1)^3: product matches and components are square-free";
}

TEST(MultivariateSqfree, MultivariateWithSquaredFactorXY12XY) {
    // f = (x+y+1)^2 * (x-y)
    std::vector<std::string> vars = {"x", "y"};

    std::vector<MultiPoly::Term> xpy1_terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(1)), // y
        make_term({0, 0}, Rational(1))  // 1
    };
    MultiPoly xpy1(xpy1_terms, vars);

    std::vector<MultiPoly::Term> xmy_terms = {
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(-1)) // -y
    };
    MultiPoly xmy(xmy_terms, vars);

    // f = (x+y+1)^2 * (x-y)
    MultiPoly f = xpy1 * xpy1 * xmy;

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((verify_sqfree_decomp(f, decomp, "x"))) << "(x+y+1)^2*(x-y): product matches and components are square-free";
}

TEST(MultivariateSqfree, ZeroPolynomial) {
    MultiPoly f; // zero polynomial

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    // Zero polynomial: decomposition should contain zero
    EXPECT_TRUE((!decomp.components.empty())) << "zero poly: decomposition is non-empty";
    EXPECT_TRUE((decomp.components[0].is_zero())) << "zero poly: first component is zero";
}
