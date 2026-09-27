/**
 * @file test_multipoly_construction_univariate.cpp
 * @brief MultiPoly 与一元多项式的构造和转换约定。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyConstructionUnivariate, ToUnivariateBasic) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(3)},
        {{1}, Rational(2)},
        {{0}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    Polynomial<Rational> uni = p.to_univariate();
    EXPECT_TRUE((uni.degree() == 2)) << "to_univariate degree is 2";
    EXPECT_TRUE((uni.coeffs[0] == Rational(1))) << "constant term is 1";
    EXPECT_TRUE((uni.coeffs[1] == Rational(2))) << "x^1 coefficient is 2";
    EXPECT_TRUE((uni.coeffs[2] == Rational(3))) << "x^2 coefficient is 3";
    EXPECT_TRUE((uni.variable_name == "x")) << "variable name is x";
}

TEST(MultipolyConstructionUnivariate, ToUnivariateWithMultiVarListButUnivariate) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{0, 3}, Rational(5)},
        {{0, 1}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    EXPECT_TRUE((p.is_univariate())) << "polynomial is univariate";
    Polynomial<Rational> uni = p.to_univariate();
    EXPECT_TRUE((uni.degree() == 3)) << "to_univariate degree is 3";
    EXPECT_TRUE((uni.coeffs[0] == Rational(0))) << "constant term is 0";
    EXPECT_TRUE((uni.coeffs[1] == Rational(1))) << "y^1 coefficient is 1";
    EXPECT_TRUE((uni.coeffs[2] == Rational(0))) << "y^2 coefficient is 0";
    EXPECT_TRUE((uni.coeffs[3] == Rational(5))) << "y^3 coefficient is 5";
    EXPECT_TRUE((uni.variable_name == "y")) << "variable name is y";
}

TEST(MultipolyConstructionUnivariate, ToUnivariateZeroPolynomial) {
    std::vector<std::string> vars = {"x"};
    MultiPoly zero(std::vector<MultiPoly::Term>{}, vars);
    Polynomial<Rational> uni = zero.to_univariate();
    EXPECT_TRUE((uni.is_zero())) << "zero MultiPoly converts to zero Polynomial";
}

TEST(MultipolyConstructionUnivariate, ToUnivariateConstantPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly c(Rational(7), vars);
    Polynomial<Rational> uni = c.to_univariate();
    EXPECT_TRUE((uni.degree() == 0)) << "constant converts to degree-0 polynomial";
    EXPECT_TRUE((uni.coeffs[0] == Rational(7))) << "constant value is 7";
}

TEST(MultipolyConstructionUnivariate, ToUnivariateThrowsForMultivariate) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 0}, Rational(1)},
        {{0, 1}, Rational(1)},
    };
    MultiPoly p(terms, vars);

    bool threw = false;
    try {
        p.to_univariate();
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "to_univariate throws for multivariate polynomial";
}

TEST(MultipolyConstructionUnivariate, FromUnivariateBasic) {
    std::vector<Rational> coeffs = {Rational(4), Rational(1), Rational(0), Rational(2)};
    Polynomial<Rational> uni(coeffs, "x");

    MultiPoly mp = MultiPoly::from_univariate(uni, "x");
    EXPECT_TRUE((mp.num_terms() == 3)) << "from_univariate has 3 non-zero terms";
    EXPECT_TRUE((mp.variables().size() == 1)) << "from_univariate has 1 variable";
    EXPECT_TRUE((mp.variables()[0] == "x")) << "variable is x";
    EXPECT_TRUE((mp.total_degree() == 3)) << "total degree is 3";
}

TEST(MultipolyConstructionUnivariate, FromUnivariateZeroPolynomial) {
    Polynomial<Rational> uni("t");
    MultiPoly mp = MultiPoly::from_univariate(uni, "t");
    EXPECT_TRUE((mp.is_zero())) << "from_univariate of zero is zero";
    EXPECT_TRUE((mp.variables().size() == 1)) << "zero poly still has variable list";
    EXPECT_TRUE((mp.variables()[0] == "t")) << "variable name preserved";
}

TEST(MultipolyConstructionUnivariate, FromUnivariateConstant) {
    Polynomial<Rational> uni(Rational(9), "z");
    MultiPoly mp = MultiPoly::from_univariate(uni, "z");
    EXPECT_TRUE((mp.is_constant())) << "from_univariate of constant is constant";
    EXPECT_TRUE((mp.num_terms() == 1)) << "constant has 1 term";
    EXPECT_TRUE((mp.terms()[0].second == Rational(9))) << "coefficient is 9";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTrip) {
    std::vector<Rational> coeffs = {Rational(7), Rational(0), Rational(-2), Rational(0), Rational(3)};
    Polynomial<Rational> original(coeffs, "x");

    MultiPoly mp = MultiPoly::from_univariate(original, "x");
    Polynomial<Rational> recovered = mp.to_univariate();

    EXPECT_TRUE((recovered == original)) << "round-trip from_univariate -> to_univariate preserves polynomial";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripMultipolyFirst) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{3}, Rational(2)},
        {{1}, Rational(5)},
    };
    MultiPoly original(terms, vars);

    Polynomial<Rational> uni = original.to_univariate();
    MultiPoly recovered = MultiPoly::from_univariate(uni, "x");

    EXPECT_TRUE((recovered == original)) << "round-trip to_univariate -> from_univariate preserves polynomial";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripZeroPolynomial) {
    Polynomial<Rational> p("x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves zero polynomial";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripConstantPolynomial) {
    Polynomial<Rational> p(Rational(42), "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves constant polynomial 42";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripLinearPolynomial) {
    std::vector<Rational> coeffs = {Rational(7), Rational(3)};
    Polynomial<Rational> p(coeffs, "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves linear polynomial 3x + 7";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripQuadraticPolynomial) {
    std::vector<Rational> coeffs = {Rational(1), Rational(-5), Rational(2)};
    Polynomial<Rational> p(coeffs, "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves quadratic 2x^2 - 5x + 1";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripHighDegreePolynomial) {
    std::vector<Rational> coeffs = {Rational(9), Rational(-2), Rational(0), Rational(0),
                                    Rational(3), Rational(0), Rational(0), Rational(1)};
    Polynomial<Rational> p(coeffs, "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves high-degree x^7 + 3x^4 - 2x + 9";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripRationalCoefficients) {
    std::vector<Rational> coeffs = {Rational(-5, 7), Rational(2, 3), Rational(0), Rational(1, 2)};
    Polynomial<Rational> p(coeffs, "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves rational coefficients (1/2)x^3 + (2/3)x - 5/7";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripMonomialSingleTerm) {
    std::vector<Rational> coeffs = {Rational(0), Rational(0), Rational(0), Rational(0), Rational(0), Rational(4)};
    Polynomial<Rational> p(coeffs, "x");
    MultiPoly mp = MultiPoly::from_univariate(p, "x");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves monomial 4x^5";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripDifferentVariableName) {
    std::vector<Rational> coeffs = {Rational(1), Rational(1), Rational(1)};
    Polynomial<Rational> p(coeffs, "t");
    MultiPoly mp = MultiPoly::from_univariate(p, "t");
    Polynomial<Rational> recovered = mp.to_univariate();
    EXPECT_TRUE((recovered == p)) << "round-trip preserves variable name 't' in t^2 + t + 1";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripReverseMultipolyToUnivariateAndBack) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{3}, Rational(5)},
        {{2}, Rational(-1)},
        {{0}, Rational(2)},
    };
    MultiPoly original(terms, vars);

    Polynomial<Rational> uni = original.to_univariate();
    MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
    EXPECT_TRUE((recovered == original)) << "reverse round-trip preserves 5x^3 - x^2 + 2";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripReverseConstantMultipoly) {
    std::vector<std::string> vars = {"x"};
    MultiPoly original(Rational(13), vars);

    Polynomial<Rational> uni = original.to_univariate();
    MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
    EXPECT_TRUE((recovered == original)) << "reverse round-trip preserves constant MultiPoly 13";
}

TEST(MultipolyConstructionUnivariate, UnivariateRoundTripReverseZeroMultipoly) {
    std::vector<std::string> vars = {"x"};
    MultiPoly original(std::vector<MultiPoly::Term>{}, vars);

    Polynomial<Rational> uni = original.to_univariate();
    MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
    EXPECT_TRUE((recovered == original)) << "reverse round-trip preserves zero MultiPoly";
}

} // namespace
