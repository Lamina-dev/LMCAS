/**
 * @file test_multipoly_construction.cpp
 * @brief 多元多项式构造契约测试。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyConstruction, ZeroConstructor) {
    MultiPoly zero;
    EXPECT_TRUE((zero.is_zero())) << "default constructor creates zero polynomial";
    EXPECT_TRUE((zero.terms().empty())) << "zero polynomial has no terms";
}

TEST(MultipolyConstruction, ConstantConstructorNonZero) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly c(Rational(5), vars);
    EXPECT_FALSE((c.is_zero())) << "constant 5 is not zero";
    EXPECT_TRUE((c.num_terms() == 1)) << "constant has one term";
    auto &terms = c.terms();
    EXPECT_TRUE((terms[0].first == Monomial({0, 0}))) << "constant monomial is [0,0]";
    EXPECT_TRUE((terms[0].second == Rational(5))) << "constant coefficient is 5";
}

TEST(MultipolyConstruction, ConstantConstructorZero) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly c(Rational(0), vars);
    EXPECT_TRUE((c.is_zero())) << "constant 0 creates zero polynomial";
    EXPECT_TRUE((c.terms().empty())) << "zero constant has no terms";
}

TEST(MultipolyConstruction, FromTermsConstructorWithNormalization) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 1}, Rational(3)},
        {{1, 1}, Rational(2)}, /**< 与前项同类的 2xy。 */
        {{2, 0}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.num_terms() == 2)) << "merged like terms: 2 terms remain";
    auto &t = p.terms();
    EXPECT_TRUE((t[0].first == Monomial({2, 0}))) << "first term is x^2";
    EXPECT_TRUE((t[0].second == Rational(1))) << "x^2 coefficient is 1";
    EXPECT_TRUE((t[1].first == Monomial({1, 1}))) << "second term is xy";
    EXPECT_TRUE((t[1].second == Rational(5))) << "xy coefficient is 5 (merged)";
}

TEST(MultipolyConstruction, NormalizeRemovesZeroCoefficientTerms) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 0}, Rational(3)},
        {{1, 0}, Rational(-3)}, /**< 与 3x 抵消。 */
        {{0, 1}, Rational(2)},
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.num_terms() == 1)) << "cancelled terms removed: 1 term remains";
    auto &t = p.terms();
    EXPECT_TRUE((t[0].first == Monomial({0, 1}))) << "remaining term is y";
    EXPECT_TRUE((t[0].second == Rational(2))) << "y coefficient is 2";
}

TEST(MultipolyConstruction, NormalizeSortsByMonomialOrderGrevlex) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{0, 0, 1}, Rational(1)}, /**< 一次项 z。 */
        {{1, 1, 1}, Rational(2)}, /**< 三次项 xyz。 */
        {{2, 0, 0}, Rational(3)}, /**< 二次项 x^2。 */
        {{0, 2, 0}, Rational(4)}, /**< 二次项 y^2。 */
    };
    MultiPoly p(terms, vars);

    auto &t = p.terms();
    EXPECT_TRUE((t[0].first == Monomial({1, 1, 1}))) << "first term is xyz (highest degree)";
    EXPECT_TRUE((t[1].first == Monomial({2, 0, 0}))) << "second term is x^2";
    EXPECT_TRUE((t[2].first == Monomial({0, 2, 0}))) << "third term is y^2";
    EXPECT_TRUE((t[3].first == Monomial({0, 0, 1}))) << "fourth term is z";
}

TEST(MultipolyConstruction, MonomialPaddingToVariableCount) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(7)}, /**< 仅给出 x 指数，补齐为 [1,0,0]。 */
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.num_terms() == 1)) << "one term after padding";
    EXPECT_TRUE((p.terms()[0].first.size() == 3)) << "monomial padded to 3 components";
    EXPECT_TRUE((p.terms()[0].first == Monomial({1, 0, 0}))) << "padded monomial is [1,0,0]";
}

TEST(MultipolyConstruction, LeadingCoeffBasicBivariate) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 1}, Rational(3)},
        {{2, 0}, Rational(2)},
        {{1, 1}, Rational(1)},
        {{0, 0}, Rational(5)},
    };
    MultiPoly p(terms, vars);

    MultiPoly lc = p.leading_coeff("x");
    EXPECT_TRUE((lc.num_terms() == 2)) << "leading_coeff(x) has 2 terms";
    EXPECT_TRUE((lc.variables().size() == 1)) << "leading_coeff(x) has 1 variable";
    EXPECT_TRUE((lc.variables()[0] == "y")) << "remaining variable is y";
    auto &lc_t = lc.terms();
    EXPECT_TRUE((lc_t[0].first == Monomial({1}))) << "first term monomial is [1] (y)";
    EXPECT_TRUE((lc_t[0].second == Rational(3))) << "first term coeff is 3";
    EXPECT_TRUE((lc_t[1].first == Monomial({0}))) << "second term monomial is [0] (constant)";
    EXPECT_TRUE((lc_t[1].second == Rational(2))) << "second term coeff is 2";
}

TEST(MultipolyConstruction, LeadingCoeffWithRespectToY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 3}, Rational(1)},
        {{1, 3}, Rational(2)},
        {{0, 2}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    MultiPoly lc = p.leading_coeff("y");
    EXPECT_TRUE((lc.num_terms() == 2)) << "leading_coeff(y) has 2 terms";
    EXPECT_TRUE((lc.variables().size() == 1)) << "leading_coeff(y) has 1 variable";
    EXPECT_TRUE((lc.variables()[0] == "x")) << "remaining variable is x";
    auto &lc_t = lc.terms();
    EXPECT_TRUE((lc_t[0].first == Monomial({2}))) << "first term is x^2";
    EXPECT_TRUE((lc_t[0].second == Rational(1))) << "x^2 coeff is 1";
    EXPECT_TRUE((lc_t[1].first == Monomial({1}))) << "second term is x";
    EXPECT_TRUE((lc_t[1].second == Rational(2))) << "x coeff is 2";
}

TEST(MultipolyConstruction, LeadingCoeffOfConstantPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly c(Rational(7), vars);

    MultiPoly lc = c.leading_coeff("x");
    EXPECT_TRUE((lc.num_terms() == 1)) << "leading_coeff of constant has 1 term";
    EXPECT_TRUE((lc.terms()[0].second == Rational(7))) << "leading_coeff of constant is 7";
}

TEST(MultipolyConstruction, LeadingCoeffOfZeroPolynomial) {
    MultiPoly zero;
    MultiPoly lc = zero.leading_coeff("x");
    EXPECT_TRUE((lc.is_zero())) << "leading_coeff of zero is zero";
}

TEST(MultipolyConstruction, LeadingCoeffWithVariableNotInPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 0}, Rational(3)},
        {{0, 0}, Rational(2)},
    };
    MultiPoly p(terms, vars);

    MultiPoly lc = p.leading_coeff("z");
    EXPECT_TRUE((lc == p)) << "leading_coeff of unknown var returns polynomial itself";
}

TEST(MultipolyConstruction, LeadingCoeffTrivariate) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 1, 1}, Rational(2)},
        {{2, 0, 2}, Rational(1)},
        {{1, 1, 0}, Rational(3)},
    };
    MultiPoly p(terms, vars);

    MultiPoly lc = p.leading_coeff("x");
    EXPECT_TRUE((lc.num_terms() == 2)) << "leading_coeff(x) trivariate has 2 terms";
    EXPECT_TRUE((lc.variables().size() == 2)) << "remaining vars are y, z";
    EXPECT_TRUE((lc.variables()[0] == "y")) << "first remaining var is y";
    EXPECT_TRUE((lc.variables()[1] == "z")) << "second remaining var is z";
}

TEST(MultipolyConstruction, ToStringOfZeroPolynomial) {
    MultiPoly zero;
    EXPECT_TRUE((zero.to_string() == "0")) << "zero poly to_string is '0'";
}

TEST(MultipolyConstruction, ToStringOfConstant) {
    std::vector<std::string> vars = {"x"};
    MultiPoly c(Rational(42), vars);
    EXPECT_TRUE((c.to_string() == "42")) << "constant 42 to_string is '42'";
}

TEST(MultipolyConstruction, ToStringBasicPolynomial) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(3)},
        {{1}, Rational(2)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    std::string s = p.to_string();
    EXPECT_TRUE((s == "3*x^2 + 2*x + 1")) << "to_string of 3x^2+2x+1";
}

TEST(MultipolyConstruction, ToStringWithCoefficient1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)},
        {{1}, Rational(1)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    std::string s = p.to_string();
    EXPECT_TRUE((s == "x^2 + x + 1")) << "to_string omits coefficient 1";
}

TEST(MultipolyConstruction, ToStringWithNegativeCoefficients) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)},
        {{1}, Rational(-2)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    std::string s = p.to_string();
    EXPECT_TRUE((s == "x^2 - 2*x + 1")) << "to_string handles negative coefficients";
}

TEST(MultipolyConstruction, ToStringMultivariate) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 1}, Rational(2)},
        {{1, 2}, Rational(3)},
    };
    MultiPoly p(terms, vars);
    std::string s = p.to_string();
    EXPECT_TRUE((s == "2*x^2*y + 3*x*y^2")) << "to_string multivariate";
}

TEST(MultipolyConstruction, ToStringNegativeLeadingCoefficient) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(-1)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    std::string s = p.to_string();
    EXPECT_TRUE((s == "-x + 1")) << "to_string with -1 leading coefficient";
}

} // namespace
