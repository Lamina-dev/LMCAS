#include "test_multivariate_factor_support.hpp"

TEST(MultivariateFactorDivision, TrialDivisionSingleFactorDividesExactly) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)),
        make_term({0, 2}, Rational(-1))
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "factor product equals original for x^2 - y^2";
}

TEST(MultivariateFactorDivision, TrialDivisionCommonMonomialExtractionThenFactor) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)),
        make_term({1, 2}, Rational(1))
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "factor product equals original for x^2*y + x*y^2";

    EXPECT_TRUE((result.factors.size() >= 2)) << "x^2*y + x*y^2 has at least 2 factors";
}

TEST(MultivariateFactorDivision, TrialDivisionIrreduciblePolynomialReturnsItself) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)),
        make_term({0, 2}, Rational(1)),
        make_term({0, 0}, Rational(1))
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    EXPECT_TRUE((result.factors.size() == 1)) << "x^2 + y^2 + 1 is irreducible (single factor)";

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "factor product equals original for irreducible poly";
}

TEST(MultivariateFactorDivision, TrialDivisionConstantPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(42), vars);

    MultiFactorResult result = checked_factor_multivariate(poly);
    EXPECT_TRUE((result.factors.empty())) << "constant polynomial has no factors";
    EXPECT_TRUE((result.constant == Rational(42))) << "constant is 42";
}

TEST(MultivariateFactorDivision, TrialDivisionZeroPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(0), vars);

    MultiFactorResult result = checked_factor_multivariate(poly);
    EXPECT_TRUE((result.factors.empty())) << "zero polynomial has no factors";
    EXPECT_TRUE((result.constant == Rational(0))) << "constant is 0";
}

TEST(MultivariateFactorDivision, FactorCombinationProductOfTwoLinearFactors) {
    std::vector<std::string> vars = {"x", "y"};

    std::vector<MultiPoly::Term> t1 = {
        make_term({1, 0}, Rational(1)),
        make_term({0, 1}, Rational(1)),
        make_term({0, 0}, Rational(1))
    };
    MultiPoly f1(t1, vars);

    std::vector<MultiPoly::Term> t2 = {
        make_term({1, 0}, Rational(1)),
        make_term({0, 1}, Rational(-1)),
        make_term({0, 0}, Rational(2))
    };
    MultiPoly f2(t2, vars);

    MultiPoly poly = f1 * f2;

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "factor product equals original for (x+y+1)(x-y+2)";
    EXPECT_TRUE((result.factors.size() == 2)) << "both nonconstant linear factors are recovered";
    bool found_first = false;
    bool found_second = false;
    for (const auto &factor : result.factors) {
        found_first = found_first || factor == f1 || factor == -f1;
        found_second = found_second || factor == f2 || factor == -f2;
    }
    EXPECT_TRUE((found_first && found_second)) << "factorization retains x+y+1 and x-y+2 up to units";
}

TEST(MultivariateFactorDivision, FactorCombinationNumericContentExtraction) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(6)),
        make_term({1, 2}, Rational(-3))
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "factor product equals original for 6x^2*y - 3x*y^2";
}
