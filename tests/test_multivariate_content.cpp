#include "test_multivariate_support.hpp"

TEST(MultivariateContent, MultivariateContentZeroPolynomialReturnsZero) {
    MultiPoly zero;
    MultiPoly content = multivariate_content(zero, "x");
    EXPECT_TRUE((content.is_zero())) << "content of zero poly is zero";
}

TEST(MultivariateContent, MultivariateContentConstantPolynomial) {
    // poly = 6, main_var = "x"
    // 视为 x 的 0 次多项式，系数为 6（常数）
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(6), vars);
    MultiPoly content = multivariate_content(poly, "x");
    EXPECT_EQ(content, poly);
}

TEST(MultivariateContent, MultivariateContentPolyDoesNotContainMainVar) {
    // poly = y^2 + y, main_var = "x"
    // poly 不含 x，容度 = poly 本身
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({0, 2}, Rational(1)), // y^2
        make_term({0, 1}, Rational(1))  // y
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");
    EXPECT_TRUE((content == poly)) << "content when main_var absent equals poly itself";
}

TEST(MultivariateContent, MultivariateContentX2XContent1AllCoeffsConstant) {
    // poly = x^2 + x, main_var = "x"
    // coeff(x^2) = 1, coeff(x^1) = 1 → content = 1
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), // x^2
        make_term({1, 0}, Rational(1))  // x
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");
    MultiPoly expected(Rational(1), {"y"});
    EXPECT_EQ(content, expected);
}

TEST(MultivariateContent, MultivariateContentUnivariatePolyX3X2XContent1) {
    // All coefficients are 1 (constants), content = gcd(1, 1, 1) = 1
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        make_term({3}, Rational(1)),
        make_term({2}, Rational(1)),
        make_term({1}, Rational(1))};
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");
    MultiPoly expected(Rational(1), {});
    EXPECT_EQ(content, expected);
}

TEST(MultivariateContent, MultivariateContentSingleTermX3Y2ContentY2) {
    // poly = x^3*y^2, main_var = "x"
    // 只有一个系数多项式：coeff(x^3) = y^2 → content = y^2
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({3, 2}, Rational(1)) // x^3 * y^2
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");

    std::vector<std::string> remaining_vars = {"y"};
    std::vector<MultiPoly::Term> expected_terms = {
        make_term({2}, Rational(1)) // y^2
    };
    MultiPoly expected(expected_terms, remaining_vars);

    EXPECT_EQ(content, expected);
}

TEST(MultivariateContent, MultivariateContentMainVarNotInVariableList) {
    // poly = x + 1, main_var = "z" (not present)
    // 容度 = poly 本身
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // 1
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "z");
    EXPECT_TRUE((content == poly)) << "content with absent main_var returns poly";
}

TEST(MultivariateContent, MultivariateContentMultipleTermsSameXDegreeGroupedCorrectly) {
    // poly = x^2*y + x^2*z, main_var = "x"
    // coeff(x^2) = y + z → single coefficient, content = y + z
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1, 0}, Rational(1)), // x^2 * y
        make_term({2, 0, 1}, Rational(1))  // x^2 * z
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");

    // 只有一个系数多项式 (y + z)，容度 = y + z
    std::vector<std::string> remaining_vars = {"y", "z"};
    std::vector<MultiPoly::Term> expected_terms = {
        make_term({1, 0}, Rational(1)), // y
        make_term({0, 1}, Rational(1))  // z
    };
    MultiPoly expected(expected_terms, remaining_vars);

    EXPECT_EQ(content, expected);
}

TEST(MultivariateContent, MultivariateContentX2YXYContentYRequiresGcd) {
    // poly = x^2*y + x*y, main_var = "x"
    // coeff(x^2) = y, coeff(x^1) = y → content = gcd(y, y) = y
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), // x^2 * y
        make_term({1, 1}, Rational(1))  // x * y
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");

    std::vector<std::string> remaining_vars = {"y"};
    std::vector<MultiPoly::Term> expected_terms = {
        make_term({1}, Rational(1)) // y
    };
    MultiPoly expected(expected_terms, remaining_vars);

    EXPECT_EQ(content, expected);
}

TEST(MultivariateContent, MultivariateContent2X2Y4XY2Content2yRequiresGcd) {
    // coeff(x^2) = 2y, coeff(x^1) = 4y^2
    // content = gcd(2y, 4y^2) = 2y
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(2)), // 2*x^2*y
        make_term({1, 2}, Rational(4))  // 4*x*y^2
    };
    MultiPoly poly(terms, vars);
    MultiPoly content = multivariate_content(poly, "x");

    std::vector<std::string> remaining_vars = {"y"};
    std::vector<MultiPoly::Term> expected_terms = {
        make_term({1}, Rational(2)) // 2*y
    };
    MultiPoly expected(expected_terms, remaining_vars);

    EXPECT_EQ(content, expected);
}
