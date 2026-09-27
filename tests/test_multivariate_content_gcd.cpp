#include "test_multivariate_support.hpp"

TEST(MultivariateContentGcd, MultivariateGcdGcdX2YXY2XYXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms_a = {
        make_term({2, 1}, Rational(1)), /**< 项 x^2*y。 */
        make_term({1, 2}, Rational(1))  /**< 项 x*y^2。 */
    };
    MultiPoly a(terms_a, vars);

    std::vector<MultiPoly::Term> terms_b = {
        make_term({1, 1}, Rational(1)) /**< 项 x*y。 */
    };
    MultiPoly b(terms_b, vars);

    MultiPoly g = multivariate_gcd(a, b);
    MultiPoly g_prim = g.make_primitive();

    MultiPoly expected(terms_b, vars);
    MultiPoly expected_prim = expected.make_primitive();

    EXPECT_EQ(g_prim, expected_prim);
    EXPECT_EQ(a.exact_div(g) * g, a);
    EXPECT_EQ(b.exact_div(g) * g, b);
}

TEST(MultivariateContentGcd, MultivariateGcdCoprimePolynomialsReturn1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms_a = {
        make_term({1}, Rational(1)), /**< 项 x。 */
        make_term({0}, Rational(1))  /**< 常数项 1。 */
    };
    MultiPoly a(terms_a, vars);

    std::vector<MultiPoly::Term> terms_b = {
        make_term({1}, Rational(1)), /**< 项 x。 */
        make_term({0}, Rational(2))  /**< 常数项 2。 */
    };
    MultiPoly b(terms_b, vars);

    MultiPoly g = multivariate_gcd(a, b);
    MultiPoly expected(Rational(1), vars);
    EXPECT_EQ(g.make_primitive(), expected);
    EXPECT_EQ(a.exact_div(g) * g, a);
    EXPECT_EQ(b.exact_div(g) * g, b);
}

TEST(MultivariateContentGcd, MultivariateGcdGcdX2Y2XYXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms_a = {
        make_term({2, 0}, Rational(1)), /**< 项 x^2。 */
        make_term({0, 2}, Rational(-1)) /**< 项 -y^2。 */
    };
    MultiPoly a(terms_a, vars);

    std::vector<MultiPoly::Term> terms_b = {
        make_term({1, 0}, Rational(1)), /**< 项 x。 */
        make_term({0, 1}, Rational(1))  /**< 项 y。 */
    };
    MultiPoly b(terms_b, vars);

    MultiPoly g = multivariate_gcd(a, b);
    MultiPoly g_prim = g.make_primitive();

    MultiPoly expected_prim = b.make_primitive();

    EXPECT_EQ(g_prim, expected_prim);
    EXPECT_EQ(a.exact_div(g) * g, a);
    EXPECT_EQ(b.exact_div(g) * g, b);
}

TEST(MultivariateContentGcd, MultivariateGcdGcd0FFUpToScalar) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly zero;

    std::vector<MultiPoly::Term> terms_f = {
        make_term({2, 1}, Rational(3)), /**< 项 3*x^2*y。 */
        make_term({1, 0}, Rational(6))  /**< 项 6*x。 */
    };
    MultiPoly f(terms_f, vars);

    MultiPoly g = multivariate_gcd(zero, f);
    MultiPoly g_prim = g.make_primitive();
    MultiPoly f_prim = f.make_primitive();

    EXPECT_EQ(g_prim, f_prim);
    EXPECT_EQ(f.exact_div(g) * g, f);
}

TEST(MultivariateContentGcd, MultivariateGcdGcdFFFUpToScalar) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms_f = {
        make_term({2}, Rational(1)), /**< 项 x^2。 */
        make_term({1}, Rational(3)), /**< 项 3*x。 */
        make_term({0}, Rational(2))  /**< 常数项 2。 */
    };
    MultiPoly f(terms_f, vars);

    MultiPoly g = multivariate_gcd(f, f);
    MultiPoly g_prim = g.make_primitive();
    MultiPoly f_prim = f.make_primitive();

    EXPECT_EQ(g_prim, f_prim);
    EXPECT_EQ(f.exact_div(g) * g, f);
}

TEST(MultivariateContentGcd, MultivariateGcdGcdOfTwoConstantsNumericGcd) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a(Rational(12), vars);
    MultiPoly b(Rational(8), vars);

    MultiPoly g = multivariate_gcd(a, b);
    MultiPoly expected(Rational(4), vars);
    EXPECT_EQ(g, expected);
    EXPECT_EQ(a.exact_div(g) * g, a);
    EXPECT_EQ(b.exact_div(g) * g, b);
}
