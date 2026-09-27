/**
 * @file test_multipoly_content.cpp
 * @brief MultiPoly 的数值内容与本原部分契约。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyContent, NumericContentOfZeroPolynomial) {
    MultiPoly zero;
    Rational content = zero.numeric_content();
    EXPECT_TRUE((content == Rational(0))) << "content of zero poly is 0";
}

TEST(MultipolyContent, NumericContentOfIntegerPolynomial) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(6)},
        {{1}, Rational(4)},
        {{0}, Rational(2)},
    };
    MultiPoly p(terms, vars);
    Rational content = p.numeric_content();
    EXPECT_TRUE((content == Rational(2))) << "content of 6x^2+4x+2 is 2";
}

TEST(MultipolyContent, NumericContentOfRationalPolynomial) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(1, 2)},
        {{0}, Rational(1, 3)},
    };
    MultiPoly p(terms, vars);
    Rational content = p.numeric_content();
    EXPECT_TRUE((content == Rational(1, 6))) << "content of (1/2)x+(1/3) is 1/6";
}

TEST(MultipolyContent, NumericContentAlreadyPrimitive) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)},
        {{1}, Rational(1)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    Rational content = p.numeric_content();
    EXPECT_TRUE((content == Rational(1))) << "content of x^2+x+1 is 1";
}

TEST(MultipolyContent, NumericContentWithNegativeCoefficients) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(-6)},
        {{0}, Rational(9)},
    };
    MultiPoly p(terms, vars);
    Rational content = p.numeric_content();
    EXPECT_TRUE((content == Rational(3))) << "content of -6x+9 is 3";
}

TEST(MultipolyContent, MakePrimitiveOfZeroPolynomial) {
    MultiPoly zero;
    MultiPoly prim = zero.make_primitive();
    EXPECT_TRUE((prim.is_zero())) << "primitive of zero is zero";
}

TEST(MultipolyContent, MakePrimitiveBasicInteger) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(6)},
        {{1}, Rational(4)},
        {{0}, Rational(2)},
    };
    MultiPoly p(terms, vars);
    MultiPoly prim = p.make_primitive();

    EXPECT_TRUE((prim.num_terms() == 3)) << "primitive has 3 terms";
    EXPECT_TRUE((prim.terms()[0].second == Rational(3))) << "leading coeff is 3";
    EXPECT_TRUE((prim.terms()[1].second == Rational(2))) << "middle coeff is 2";
    EXPECT_TRUE((prim.terms()[2].second == Rational(1))) << "constant is 1";
}

TEST(MultipolyContent, MakePrimitiveEnsuresPositiveLeadingCoefficient) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(-6)},
        {{0}, Rational(9)},
    };
    MultiPoly p(terms, vars);
    MultiPoly prim = p.make_primitive();

    EXPECT_TRUE((prim.terms()[0].second > Rational(0))) << "leading coefficient is positive";
    EXPECT_TRUE((prim.terms()[0].second == Rational(2))) << "leading coeff is 2";
    EXPECT_TRUE((prim.terms()[1].second == Rational(-3))) << "constant is -3";
}

TEST(MultipolyContent, MakePrimitiveOfAlreadyPrimitivePolynomial) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)},
        {{1}, Rational(1)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);
    MultiPoly prim = p.make_primitive();
    EXPECT_TRUE((prim == p)) << "already primitive polynomial unchanged";
}

TEST(MultipolyContent, MakePrimitiveWithRationalCoefficients) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{1}, Rational(1, 2)},
        {{0}, Rational(1, 3)},
    };
    MultiPoly p(terms, vars);
    MultiPoly prim = p.make_primitive();

    EXPECT_TRUE((prim.terms()[0].second == Rational(3))) << "leading coeff is 3";
    EXPECT_TRUE((prim.terms()[1].second == Rational(2))) << "constant is 2";
}

} // namespace
