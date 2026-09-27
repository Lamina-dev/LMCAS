/**
 * @file test_multipoly_division_degree.cpp
 * @brief MultiPoly 精确除法与次数契约测试。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyDivisionDegree, ExactDivBasic) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0}, Rational(1)},  /**< x^2 项。 */
        {{0, 2}, Rational(-1)}, /**< -y^2 项。 */
    };
    MultiPoly f(f_terms, vars);

    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)}, /**< x 项。 */
        {{0, 1}, Rational(1)}, /**< y 项。 */
    };
    MultiPoly g(g_terms, vars);

    MultiPoly q = f.exact_div(g);

    std::vector<MultiPoly::Term> expected_terms = {
        {{1, 0}, Rational(1)},  /**< x 项。 */
        {{0, 1}, Rational(-1)}, /**< -y 项。 */
    };
    MultiPoly expected(expected_terms, vars);

    EXPECT_TRUE((q == expected)) << "exact_div: (x^2 - y^2) / (x + y) == x - y";
}

TEST(MultipolyDivisionDegree, ExactDivProductRoundTrip) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(2)}, /**< 2x 项。 */
        {{0, 1}, Rational(3)}, /**< 3y 项。 */
    };
    MultiPoly g(g_terms, vars);

    std::vector<MultiPoly::Term> h_terms = {
        {{1, 0}, Rational(1)},  /**< x 项。 */
        {{0, 1}, Rational(-1)}, /**< -y 项。 */
    };
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;

    MultiPoly q1 = f.exact_div(g);
    EXPECT_TRUE((q1 == h)) << "exact_div: (g*h) / g == h";

    MultiPoly q2 = f.exact_div(h);
    EXPECT_TRUE((q2 == g)) << "exact_div: (g*h) / h == g";
}

TEST(MultipolyDivisionDegree, ExactDivByConstant) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0}, Rational(6)}, /**< 6x^2 项。 */
        {{1, 1}, Rational(4)}, /**< 4xy 项。 */
    };
    MultiPoly f(f_terms, vars);

    MultiPoly divisor(Rational(2), vars);

    MultiPoly q = f.exact_div(divisor);

    std::vector<MultiPoly::Term> expected_terms = {
        {{2, 0}, Rational(3)}, /**< 3x^2 项。 */
        {{1, 1}, Rational(2)}, /**< 2xy 项。 */
    };
    MultiPoly expected(expected_terms, vars);

    EXPECT_TRUE((q == expected)) << "exact_div by constant: (6x^2 + 4xy) / 2 == 3x^2 + 2xy";
}

TEST(MultipolyDivisionDegree, ExactDivZeroDividend) {
    std::vector<std::string> vars = {"x", "y"};
    const MultiPoly zero;
    const MultiPoly other_ring_zero(std::vector<MultiPoly::Term>{}, {"z"});

    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)}, /**< x 项。 */
    };
    MultiPoly g(g_terms, vars);

    for (const auto *dividend : {&zero, &other_ring_zero}) {
        const MultiPoly q = dividend->exact_div(g);
        EXPECT_TRUE((q.is_zero() && q.variables() == vars)) << "exact_div: 0 / g is zero in the divisor ring";
    }
}

TEST(MultipolyDivisionDegree, ExactDivThrowsOnZeroDivisor) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 0}, Rational(1)}, /**< x 项。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly zero(std::vector<MultiPoly::Term>{}, vars);

    bool threw = false;
    try {
        f.exact_div(zero);
    } catch (const std::runtime_error &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "exact_div throws on zero divisor";
}

TEST(MultipolyDivisionDegree, ExactDivThrowsOnNonExactDivision) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0}, Rational(1)}, /**< x^2 项。 */
        {{0, 0}, Rational(1)}, /**< 1 项。 */
    };
    MultiPoly f(f_terms, vars);

    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)}, /**< x 项。 */
        {{0, 1}, Rational(1)}, /**< y 项。 */
    };
    MultiPoly g(g_terms, vars);

    bool threw = false;
    try {
        f.exact_div(g);
    } catch (const std::runtime_error &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "exact_div throws on non-exact division";
}

TEST(MultipolyDivisionDegree, ExactDivRejectsIncompatibleRings) {
    const MultiPoly x({{{1}, Rational(1)}}, {"x"});
    const MultiPoly y({{{1}, Rational(1)}}, {"y"});
    const MultiPoly xy({{{1, 0}, Rational(1)}}, {"x", "y"});
    const MultiPoly yx({{{1, 0}, Rational(1)}}, {"y", "x"});
    const MultiPoly one_x(Rational(1), {"x"});
    const MultiPoly one_y(Rational(1), {"y"});
    for (const auto &operands : {std::pair{x, y}, std::pair{xy, yx},
                                 std::pair{x, xy}, std::pair{one_x, one_y}}) {
        bool rejected = false;
        try {
            (void)operands.first.exact_div(operands.second);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        EXPECT_TRUE((rejected)) << "nonzero cross-ring exact division raises invalid_argument";
    }

    bool zero_rejected = false;
    try {
        (void)x.exact_div(MultiPoly(std::vector<MultiPoly::Term>{}, {"y"}));
    } catch (const std::runtime_error &) {
        zero_rejected = true;
    }
    EXPECT_TRUE((zero_rejected)) << "zero division remains runtime_error even with different variables";
}

TEST(MultipolyDivisionDegree, ExactDivRespectsTheDividendOrder) {
    const std::vector<std::string> vars{"x", "y"};
    const std::vector<MultiPoly::Term> numerator{
        {{2, 0}, Rational(1)}, {{0, 4}, Rational(-1)}};
    const std::vector<MultiPoly::Term> denominator{
        {{1, 0}, Rational(1)}, {{0, 2}, Rational(1)}};
    const std::vector<MultiPoly::Term> quotient{
        {{1, 0}, Rational(1)}, {{0, 2}, Rational(-1)}};
    for (const auto dividend_order : {MonomialOrderType::Lex, MonomialOrderType::GrevLex}) {
        for (const auto divisor_order : {MonomialOrderType::Lex, MonomialOrderType::GrevLex}) {
            const MultiPoly dividend(numerator, vars, dividend_order);
            const MultiPoly divisor(denominator, vars, divisor_order);
            const MultiPoly expected(quotient, vars, dividend_order);
            const auto result = dividend.exact_div(divisor);
            EXPECT_TRUE((result == expected && result.variables() == vars)) << "(x^2-y^4)/(x+y^2) is x-y^2 for either ordering";
            EXPECT_TRUE((result.terms() == expected.terms())) << "the quotient retains the dividend monomial ordering";
            EXPECT_TRUE((divisor.terms() == MultiPoly(denominator, vars, divisor_order).terms())) << "exact division leaves the differently ordered divisor unchanged";
        }
    }
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripSimpleMonomialsXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> g_terms = {{{1, 0}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{0, 1}, Rational(1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (x*y)/x == y";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (x*y)/y == x";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripBinomialsX1X1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(1)}, {{0}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(1)}, {{0}, Rational(-1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (x^2-1)/(x+1) == x-1";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (x^2-1)/(x-1) == x+1";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripBivariateXYXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> g_terms = {{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{1, 0}, Rational(1)}, {{0, 1}, Rational(-1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (x^2-y^2)/(x+y) == x-y";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (x^2-y^2)/(x-y) == x+y";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripWithCoefficients2x13x2) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(2)}, {{0}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(3)}, {{0}, Rational(2)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (6x^2+7x+2)/(2x+1) == 3x+2";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (6x^2+7x+2)/(3x+2) == 2x+1";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripHigherDegreeX2X1X1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> g_terms = {
        {{2}, Rational(1)}, {{1}, Rational(1)}, {{0}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(1)}, {{0}, Rational(-1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (x^3-1)/(x^2+x+1) == x-1";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (x^3-1)/(x-1) == x^2+x+1";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripTrivariateXYZXY) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(1)}, {{0, 0, 1}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {
        {{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(-1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: ((x+y+z)*(x-y))/(x+y+z) == x-y";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: ((x+y+z)*(x-y))/(x-y) == x+y+z";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripMonomialTimesPolynomialXyX2Y) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> g_terms = {{{1, 1}, Rational(1)}};
    std::vector<MultiPoly::Term> h_terms = {{{2, 0}, Rational(1)}, {{0, 1}, Rational(1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: (x^3*y+x*y^2)/(xy) == x^2+y";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: (x^3*y+x*y^2)/(x^2+y) == xy";
}

TEST(MultipolyDivisionDegree, ExactDivRoundTripRationalCoefficientsX2133x1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(1, 2)}, {{0}, Rational(1, 3)}};
    std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(3)}, {{0}, Rational(-1)}};
    MultiPoly g(g_terms, vars);
    MultiPoly h(h_terms, vars);

    MultiPoly f = g * h;
    EXPECT_TRUE((f.exact_div(g) == h)) << "exact_div: f/(x/2+1/3) == 3x-1";
    EXPECT_TRUE((f.exact_div(h) == g)) << "exact_div: f/(3x-1) == x/2+1/3";
}

TEST(MultipolyDivisionDegree, DegreeOfZeroPolynomialIs1) {
    MultiPoly zero;
    EXPECT_TRUE((zero.total_degree() == -1)) << "zero poly total_degree is -1";
    EXPECT_TRUE((zero.degree("x") == -1)) << "zero poly degree(x) is -1";
    EXPECT_TRUE((zero.degree("y") == -1)) << "zero poly degree(y) is -1";
}

TEST(MultipolyDivisionDegree, DegreeOfConstantPolynomialIs0) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly c(Rational(42), vars);
    EXPECT_TRUE((c.total_degree() == 0)) << "constant poly total_degree is 0";
    EXPECT_TRUE((c.degree("x") == 0)) << "constant poly degree(x) is 0";
    EXPECT_TRUE((c.degree("y") == 0)) << "constant poly degree(y) is 0";
}

TEST(MultipolyDivisionDegree, DegreeOfUnivariatePolynomial) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{4}, Rational(5)},
        {{2}, Rational(3)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 4)) << "5x^4+3x^2+1 total_degree is 4";
    EXPECT_TRUE((p.degree("x") == 4)) << "5x^4+3x^2+1 degree(x) is 4";
}

TEST(MultipolyDivisionDegree, DegreeOfBivariatePolynomial) {
    /**
     * @brief 双变量多项式 p = 2x^3*y^2 + x*y^5 + 3x^2 的次数。
     * 各项总次数为 3+2=5、1+5=6、2+0=2。
     * total_degree = 6、degree(y) = 5，均由 x*y^5 决定；
     * degree(x) = 3，由 x^3*y^2 决定。
     */
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{3, 2}, Rational(2)}, /**< 2x^3*y^2 项。 */
        {{1, 5}, Rational(1)}, /**< x*y^5 项。 */
        {{2, 0}, Rational(3)}, /**< 3x^2 项。 */
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 6)) << "2x^3y^2+xy^5+3x^2 total_degree is 6";
    EXPECT_TRUE((p.degree("x") == 3)) << "2x^3y^2+xy^5+3x^2 degree(x) is 3";
    EXPECT_TRUE((p.degree("y") == 5)) << "2x^3y^2+xy^5+3x^2 degree(y) is 5";
}

TEST(MultipolyDivisionDegree, DegreeOfTrivariatePolynomial) {
    /**
     * @brief 三变量多项式 p = x^2*y*z^3 + x*y^4*z + z^7 的次数。
     * 各项总次数为 2+1+3=6、1+4+1=6、0+0+7=7。
     * total_degree = degree(z) = 7，均由 z^7 决定；
     * degree(x) = 2 由 x^2*y*z^3 决定，degree(y) = 4 由 x*y^4*z 决定。
     */
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 1, 3}, Rational(1)}, /**< x^2*y*z^3 项。 */
        {{1, 4, 1}, Rational(1)}, /**< x*y^4*z 项。 */
        {{0, 0, 7}, Rational(1)}, /**< z^7 项。 */
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 7)) << "x^2yz^3+xy^4z+z^7 total_degree is 7";
    EXPECT_TRUE((p.degree("x") == 2)) << "x^2yz^3+xy^4z+z^7 degree(x) is 2";
    EXPECT_TRUE((p.degree("y") == 4)) << "x^2yz^3+xy^4z+z^7 degree(y) is 4";
    EXPECT_TRUE((p.degree("z") == 7)) << "x^2yz^3+xy^4z+z^7 degree(z) is 7";
}

TEST(MultipolyDivisionDegree, DegreeWithHighDegreeSingleTerm) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{10, 8, 6}, Rational(7)},
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 24)) << "x^10*y^8*z^6 total_degree is 24";
    EXPECT_TRUE((p.degree("x") == 10)) << "x^10*y^8*z^6 degree(x) is 10";
    EXPECT_TRUE((p.degree("y") == 8)) << "x^10*y^8*z^6 degree(y) is 8";
    EXPECT_TRUE((p.degree("z") == 6)) << "x^10*y^8*z^6 degree(z) is 6";
}

TEST(MultipolyDivisionDegree, DegreeForVariableNotAppearingInPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{3, 0}, Rational(1)}, /**< x^3 项。 */
        {{0, 0}, Rational(1)}, /**< 1 项。 */
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 3)) << "x^3+1 total_degree is 3";
    EXPECT_TRUE((p.degree("x") == 3)) << "x^3+1 degree(x) is 3";
    EXPECT_TRUE((p.degree("y") == 0)) << "x^3+1 degree(y) is 0 (y not used)";
}

TEST(MultipolyDivisionDegree, DegreeConsistencyWithArithmetic) {
    /**
     * @brief 乘积满足总次数等式 deg(f*g) = deg(f) + deg(g)。
     * f = x^2 + y、g = x*y + 1，总次数均为 2；
     * f*g = x^3*y + x^2 + x*y^2 + y，总次数为 4。
     */
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0}, Rational(1)}, /**< x^2 项。 */
        {{0, 1}, Rational(1)}, /**< y 项。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 1}, Rational(1)}, /**< xy 项。 */
        {{0, 0}, Rational(1)}, /**< 1 项。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);
    MultiPoly fg = f * g;

    int deg_f = f.total_degree();
    int deg_g = g.total_degree();
    int deg_fg = fg.total_degree();

    EXPECT_TRUE((deg_f == 2)) << "f=x^2+y total_degree is 2";
    EXPECT_TRUE((deg_g == 2)) << "g=xy+1 total_degree is 2";
    EXPECT_TRUE((deg_fg == deg_f + deg_g)) << "deg(f*g) == deg(f) + deg(g)";
}

TEST(MultipolyDivisionDegree, DegreeOfHomogeneousPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{3, 0}, Rational(1)}, /**< x^3 项。 */
        {{2, 1}, Rational(1)}, /**< x^2*y 项。 */
        {{1, 2}, Rational(1)}, /**< x*y^2 项。 */
        {{0, 3}, Rational(1)}, /**< y^3 项。 */
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.total_degree() == 3)) << "homogeneous poly total_degree is 3";
    EXPECT_TRUE((p.degree("x") == 3)) << "homogeneous poly degree(x) is 3";
    EXPECT_TRUE((p.degree("y") == 3)) << "homogeneous poly degree(y) is 3";
    EXPECT_TRUE((p.is_homogeneous())) << "polynomial is homogeneous";
}

} // namespace
