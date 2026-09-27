#include "test_multivariate_factor_support.hpp"

static MultiPoly reconstruct_factorization(const MultiFactorResult &result,
                                           const std::vector<std::string> &vars) {
    MultiPoly reconstructed(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            reconstructed = reconstructed * result.factors[i];
        }
    }
    return reconstructed;
}

static MultiPoly homogeneous_bivariate_linear(int x, int y, bool fallback_y,
                                              const std::vector<std::string> &vars) {
    std::vector<MultiPoly::Term> terms;
    if (x != 0) {
        terms.push_back({Monomial({1, 0}), Rational(x)});
    }
    if (y != 0) {
        terms.push_back({Monomial({0, 1}), Rational(y)});
    }
    if (terms.empty()) {
        terms.push_back({fallback_y ? Monomial({0, 1}) : Monomial({1, 0}), Rational(1)});
    }
    return MultiPoly(terms, vars);
}

TEST(MultivariateFactorProperties, ConstantHenselCorrectionKeepsTheUnivariateFactorDomain) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly({make_term({1, 1}, Rational(-2)),
                    make_term({0, 2}, Rational(1)),
                    make_term({1, 0}, Rational(6)),
                    make_term({0, 1}, Rational(-5)),
                    make_term({0, 0}, Rational(6))},
                   vars);

    MultiFactorResult result = checked_factor_multivariate(poly);
    MultiPoly reconstructed = reconstruct_factorization(result, vars);
    EXPECT_TRUE((reconstructed == poly)) << "constant Hensel correction reconstructs -2xy+y^2+6x-5y+6";
    const MultiPoly first({make_term({0, 1}, Rational(1)),
                           make_term({0, 0}, Rational(-3))},
                          vars);
    const MultiPoly second({make_term({1, 0}, Rational(-2)),
                            make_term({0, 1}, Rational(1)),
                            make_term({0, 0}, Rational(-2))},
                           vars);
    EXPECT_TRUE((result.factors.size() == 2)) << "both linear factors are recovered";
    bool found_first = false;
    bool found_second = false;
    for (const auto &factor : result.factors) {
        found_first = found_first || factor == first || factor == -first;
        found_second = found_second || factor == second || factor == -second;
    }
    EXPECT_TRUE((found_first && found_second)) << "factorization retains y-3 and y-2-2x up to units";
}

TEST(MultivariateFactorProperties, LinearRecombinationPreservesMainVariableUnit) {
    const std::vector<std::string> vars = {"lift", "spare", "main"};
    const MultiPoly first({make_term({0, 0, 1}, Rational(2)),
                           make_term({1, 0, 0}, Rational(1)),
                           make_term({0, 1, 0}, Rational(1)),
                           make_term({0, 0, 0}, Rational(1))},
                          vars);
    const MultiPoly second({make_term({0, 0, 1}, Rational(1)),
                            make_term({0, 0, 0}, Rational(-3))},
                           vars);
    const MultiPoly poly = first * second;
    const auto result = checked_factor_multivariate(poly);
    EXPECT_TRUE((reconstruct_factorization(result, vars) == poly)) << "recombined factors and unit exactly reconstruct the ordered ring polynomial";
    EXPECT_TRUE((result.factors.size() == 2)) << "both primitive linear factors are recovered";
    bool found_first = false;
    bool found_second = false;
    for (const auto &factor : result.factors) {
        found_first = found_first || factor == first || factor == -first;
        found_second = found_second || factor == second || factor == -second;
    }
    EXPECT_TRUE((found_first && found_second)) << "monic interpolation does not discard a rational main-variable unit";
}

TEST(MultivariateFactorProperties, FactorizationProductCorrectnessRandomFactorableBivariate) {
    EXPECT_TRUE(rc::check("For random products of linear factors, factorization reconstructs original", []() {
        std::vector<std::string> vars = {"x", "y"};
        int num_factors = 2 + (*rc::gen::inRange(0, (2) + 1));

        MultiPoly product(Rational(1), vars);
        for (int i = 0; i < num_factors; ++i) {
            int a = (*rc::gen::inRange(-3, (3) + 1));
            int b = (*rc::gen::inRange(-3, (3) + 1));
            int c = (*rc::gen::inRange(-3, (3) + 1));
            if (a == 0 && b == 0) {
                a = 1;
            }

            std::vector<MultiPoly::Term> terms;
            if (a != 0)
                terms.push_back({Monomial({1, 0}), Rational(a)});
            if (b != 0)
                terms.push_back({Monomial({0, 1}), Rational(b)});
            if (c != 0)
                terms.push_back({Monomial({0, 0}), Rational(c)});
            if (terms.empty())
                terms.push_back({Monomial({1, 0}), Rational(1)});

            MultiPoly factor(terms, vars);
            product = product * factor;
        }

        MultiFactorResult result = checked_factor_multivariate(product);

        MultiPoly reconstructed = reconstruct_factorization(result, vars);

        RC_ASSERT(reconstructed == product);
    }));
}

TEST(MultivariateFactorProperties, EachFactorIsPrimitiveWithPositiveLeadingCoefficient) {
    EXPECT_TRUE(rc::check("Factors from factorization are primitive with positive leading coefficient", []() {
        std::vector<std::string> vars = {"x", "y"};
        int a = (*rc::gen::inRange(1, (4) + 1)); /**< 正首项系数。 */
        int b = (*rc::gen::inRange(-3, (3) + 1));
        int c = (*rc::gen::inRange(-3, (3) + 1));
        int d = (*rc::gen::inRange(-3, (3) + 1));
        int e = (*rc::gen::inRange(-3, (3) + 1));
        int f = (*rc::gen::inRange(-3, (3) + 1));

        std::vector<MultiPoly::Term> terms;
        terms.push_back({Monomial({2, 0}), Rational(a)});
        if (b != 0)
            terms.push_back({Monomial({1, 1}), Rational(b)});
        if (c != 0)
            terms.push_back({Monomial({0, 2}), Rational(c)});
        if (d != 0)
            terms.push_back({Monomial({1, 0}), Rational(d)});
        if (e != 0)
            terms.push_back({Monomial({0, 1}), Rational(e)});
        if (f != 0)
            terms.push_back({Monomial({0, 0}), Rational(f)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = checked_factor_multivariate(poly);

        for (size_t i = 0; i < result.factors.size(); ++i) {
            Rational content = result.factors[i].numeric_content();
            RC_ASSERT(content == Rational(1));
        }

        MultiPoly reconstructed = reconstruct_factorization(result, vars);
        RC_ASSERT(reconstructed == poly);
    }));
}

TEST(MultivariateFactorProperties, FactorizationProductCorrectnessRandomMonomialLinear) {
    EXPECT_TRUE(rc::check("For monomial * linear factor products, factorization reconstructs original", []() {
        std::vector<std::string> vars = {"x", "y"};
        int exp_x = (*rc::gen::inRange(0, (3) + 1));
        int exp_y = (*rc::gen::inRange(0, (3) + 1));
        if (exp_x == 0 && exp_y == 0) {
            exp_x = 1;
        }

        std::vector<MultiPoly::Term> mono_terms = {
            {Monomial({exp_x, exp_y}), Rational(1)}};
        MultiPoly monomial_factor(mono_terms, vars);

        int a = (*rc::gen::inRange(-3, (3) + 1));
        int b = (*rc::gen::inRange(-3, (3) + 1));
        int c = (*rc::gen::inRange(-3, (3) + 1));
        if (a == 0 && b == 0) {
            a = 1;
        }

        std::vector<MultiPoly::Term> lin_terms;
        if (a != 0)
            lin_terms.push_back({Monomial({1, 0}), Rational(a)});
        if (b != 0)
            lin_terms.push_back({Monomial({0, 1}), Rational(b)});
        if (c != 0)
            lin_terms.push_back({Monomial({0, 0}), Rational(c)});
        if (lin_terms.empty())
            lin_terms.push_back({Monomial({1, 0}), Rational(1)});
        MultiPoly linear_factor(lin_terms, vars);

        int scalar = (*rc::gen::inRange(1, (5) + 1));
        MultiPoly poly = monomial_factor * linear_factor * Rational(scalar);

        MultiFactorResult result = checked_factor_multivariate(poly);

        MultiPoly reconstructed = reconstruct_factorization(result, vars);
        RC_ASSERT(reconstructed == poly);
    }));
}

TEST(MultivariateFactorProperties, LinearPolynomialsAreIrreducibleRandomBivariateLinear) {
    EXPECT_TRUE(rc::check("Linear bivariate polynomial returns itself as sole factor", []() {
        std::vector<std::string> vars = {"x", "y"};
        int a = (*rc::gen::inRange(-5, (5) + 1));
        int b = (*rc::gen::inRange(-5, (5) + 1));
        int c = (*rc::gen::inRange(-5, (5) + 1));
        if (a == 0 && b == 0) {
            a = 1;
        }

        std::vector<MultiPoly::Term> terms;
        if (a != 0)
            terms.push_back({Monomial({1, 0}), Rational(a)});
        if (b != 0)
            terms.push_back({Monomial({0, 1}), Rational(b)});
        if (c != 0)
            terms.push_back({Monomial({0, 0}), Rational(c)});
        if (terms.empty())
            terms.push_back({Monomial({1, 0}), Rational(1)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = checked_factor_multivariate(poly);

        RC_ASSERT(result.factors.size() == 1);
        RC_ASSERT(result.multiplicities[0] == 1);

        MultiPoly reconstructed(Rational(result.constant), vars);
        reconstructed = reconstructed * result.factors[0];
        RC_ASSERT(reconstructed == poly);
    }));
}

TEST(MultivariateFactorProperties, LinearPolynomialsAreIrreducibleRandomTrivariateLinear) {
    EXPECT_TRUE(rc::check("Linear trivariate polynomial returns itself as sole factor", []() {
        std::vector<std::string> vars = {"x", "y", "z"};
        int a = (*rc::gen::inRange(-4, (4) + 1));
        int b = (*rc::gen::inRange(-4, (4) + 1));
        int c = (*rc::gen::inRange(-4, (4) + 1));
        int d = (*rc::gen::inRange(-4, (4) + 1));
        if (a == 0 && b == 0 && c == 0) {
            a = 1;
        }

        std::vector<MultiPoly::Term> terms;
        if (a != 0)
            terms.push_back({Monomial({1, 0, 0}), Rational(a)});
        if (b != 0)
            terms.push_back({Monomial({0, 1, 0}), Rational(b)});
        if (c != 0)
            terms.push_back({Monomial({0, 0, 1}), Rational(c)});
        if (d != 0)
            terms.push_back({Monomial({0, 0, 0}), Rational(d)});
        if (terms.empty())
            terms.push_back({Monomial({1, 0, 0}), Rational(1)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = checked_factor_multivariate(poly);

        RC_ASSERT(result.factors.size() == 1);
        RC_ASSERT(result.multiplicities[0] == 1);

        MultiPoly reconstructed(Rational(result.constant), vars);
        reconstructed = reconstructed * result.factors[0];
        RC_ASSERT(reconstructed == poly);
    }));
}

TEST(MultivariateFactorProperties, DifferenceOfSquaresFactorizationRandomA2B2) {
    EXPECT_TRUE(rc::check("a^2 - b^2 factors into product containing (a+b) and (a-b)", []() {
        std::vector<std::string> vars = {"x", "y"};
        int c1 = (*rc::gen::inRange(-3, (3) + 1));
        int c2 = (*rc::gen::inRange(-3, (3) + 1));
        if (c1 == 0 && c2 == 0) {
            c1 = 1;
        }

        int c3 = (*rc::gen::inRange(-3, (3) + 1));
        int c4 = (*rc::gen::inRange(-3, (3) + 1));
        if (c3 == 0 && c4 == 0) {
            c4 = 1;
        }

        if (c1 == c3 && c2 == c4) {
            c3 = c3 + 1;
        }
        if (c1 == -c3 && c2 == -c4) {
            c4 = c4 + 1;
        }

        MultiPoly a_poly = homogeneous_bivariate_linear(c1, c2, false, vars);
        MultiPoly b_poly = homogeneous_bivariate_linear(c3, c4, true, vars);

        MultiPoly a_sq = a_poly * a_poly;
        MultiPoly b_sq = b_poly * b_poly;
        MultiPoly diff = a_sq - b_sq;

        if (diff.is_zero()) {
            return;
        }

        MultiFactorResult result = checked_factor_multivariate(diff);

        MultiPoly reconstructed = reconstruct_factorization(result, vars);
        RC_ASSERT(reconstructed == diff);

        if (diff.num_terms() == 2) {
            int total_factor_count = 0;
            for (size_t i = 0; i < result.multiplicities.size(); ++i) {
                total_factor_count += result.multiplicities[i];
            }
            RC_ASSERT(total_factor_count >= 2);
        }
    }));
}

TEST(MultivariateFactorProperties, DifferenceOfSquaresWithSingleVariables) {
    EXPECT_TRUE(rc::check("x_i^2 - x_j^2 factors into (x_i + x_j)(x_i - x_j)", []() {
        std::vector<std::string> vars = {"x", "y", "z"};
        int i = (*rc::gen::inRange(0, (2) + 1));
        int j = (*rc::gen::inRange(0, (2) + 1));
        if (i == j) {
            j = (i + 1) % 3;
        }

        std::vector<int> exp_pos(3, 0);
        exp_pos[i] = 2;
        std::vector<int> exp_neg(3, 0);
        exp_neg[j] = 2;

        std::vector<MultiPoly::Term> terms = {
            {Monomial(exp_pos.begin(), exp_pos.end()), Rational(1)},
            {Monomial(exp_neg.begin(), exp_neg.end()), Rational(-1)}};
        MultiPoly poly(terms, vars);

        MultiFactorResult result = checked_factor_multivariate(poly);

        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t k = 0; k < result.factors.size(); ++k) {
            for (int m = 0; m < result.multiplicities[k]; ++m) {
                reconstructed = reconstructed * result.factors[k];
            }
        }
        RC_ASSERT(reconstructed == poly);

        RC_ASSERT(result.factors.size() == 2);
        RC_ASSERT(result.multiplicities[0] == 1);
        RC_ASSERT(result.multiplicities[1] == 1);
    }));
}
