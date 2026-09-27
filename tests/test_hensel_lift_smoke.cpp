/**
 * @file test_hensel_lift_smoke.cpp
 * @brief 多元 Hensel 提升冒烟测试。
 */
#include "test_multivariate_support.hpp"

using namespace LMCAS;

TEST(HenselLiftSmoke, HenselLiftXY1X1FromFactorsAtY0) {
    std::vector<std::string> vars = {"x", "y"};

    // poly = (x+y+1)(x-1) = x^2 + xy - y - 1
    std::vector<MultiPoly::Term> poly_terms = {
        make_term({2, 0}, Rational(1)),  // x^2
        make_term({1, 1}, Rational(1)),  // xy
        make_term({0, 1}, Rational(-1)), // -y
        make_term({0, 0}, Rational(-1))  // -1
    };
    MultiPoly poly(poly_terms, vars);

    // univariate factors at y=0: x+1, x-1
    Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");  // x + 1
    Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x"); // x - 1

    std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 1);

    EXPECT_TRUE((lifted.size() == 2)) << "lifted has 2 factors";

    // Check product equals poly
    if (lifted.size() == 2) {
        MultiPoly product = lifted[0] * lifted[1];
        EXPECT_TRUE((product == poly)) << "product of lifted factors == poly";
    }
}

TEST(HenselLiftSmoke, HenselLiftX1X1NoYDependenceDegreeBound2) {
    std::vector<std::string> vars = {"x", "y"};

    // poly = x^2 - 1 (no y terms)
    std::vector<MultiPoly::Term> poly_terms = {
        make_term({2, 0}, Rational(1)), // x^2
        make_term({0, 0}, Rational(-1)) // -1
    };
    MultiPoly poly(poly_terms, vars);

    // univariate factors at y=0: x+1, x-1
    Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");  // x + 1
    Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x"); // x - 1

    std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

    EXPECT_TRUE((lifted.size() == 2)) << "lifted has 2 factors";

    if (lifted.size() == 2) {
        MultiPoly product = lifted[0] * lifted[1];
        EXPECT_TRUE((product == poly)) << "product of lifted factors == poly (no y)";
    }
}

TEST(HenselLiftSmoke, HenselLiftSingleFactorReturnsPolyItself) {
    std::vector<std::string> vars = {"x", "y"};

    std::vector<MultiPoly::Term> poly_terms = {
        make_term({2, 1}, Rational(1)), // x^2*y
        make_term({1, 0}, Rational(1))  // x
    };
    MultiPoly poly(poly_terms, vars);

    Polynomial<Rational> f1({Rational(0), Rational(1)}, "x"); // x

    std::vector<Polynomial<Rational>> uni_factors = {f1};

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

    EXPECT_TRUE((lifted.size() == 1)) << "single factor: 1 result";
    if (lifted.size() == 1) {
        EXPECT_TRUE((lifted[0] == poly)) << "single factor returns poly itself";
    }
}

TEST(HenselLiftSmoke, HenselLiftEmptyFactorsReturnsEmpty) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(1), vars);
    std::vector<Polynomial<Rational>> uni_factors;

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);
    EXPECT_TRUE((lifted.empty())) << "empty factors returns empty";
}

TEST(HenselLiftSmoke, HenselLiftX2y1XY1FromFactorsAtY0) {
    std::vector<std::string> vars = {"x", "y"};

    // poly = (x+2y+1)(x-y+1) = x^2 + xy + 2x - 2y^2 - y + 1
    // At y=0: (x+1)(x+1) = x^2 + 2x + 1 — but this is not square-free!
    // Let's use a different example.
    // poly = (x+y)(x-y) = x^2 - y^2
    // At y=0: x^2 = x*x — also not coprime.
    // Better: poly = (x+y+1)(x-y-1) = x^2 - y^2 - 2y - 1 + x*0
    // Wait: (x+y+1)(x-y-1) = x^2 - xy - x + xy - y^2 - y + x - y - 1
    //     = x^2 - y^2 - 2y - 1
    // At y=0: x^2 - 1 = (x+1)(x-1)
    std::vector<MultiPoly::Term> poly_terms = {
        make_term({2, 0}, Rational(1)),  // x^2
        make_term({0, 2}, Rational(-1)), // -y^2
        make_term({0, 1}, Rational(-2)), // -2y
        make_term({0, 0}, Rational(-1))  // -1
    };
    MultiPoly poly(poly_terms, vars);

    // univariate factors at y=0: x+1, x-1
    Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");  // x + 1
    Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x"); // x - 1

    std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

    EXPECT_TRUE((lifted.size() == 2)) << "lifted has 2 factors";

    if (lifted.size() == 2) {
        MultiPoly product = lifted[0] * lifted[1];
        EXPECT_TRUE((product == poly)) << "product of lifted factors == poly (degree 2 in y)";
    }
}

TEST(HenselLiftSmoke, HenselLiftNonZeroEvaluationPointXY1X1AtY1) {
    std::vector<std::string> vars = {"x", "y"};

    // poly = (x+y+1)(x-1) = x^2 + xy - y - 1
    std::vector<MultiPoly::Term> poly_terms = {
        make_term({2, 0}, Rational(1)),  // x^2
        make_term({1, 1}, Rational(1)),  // xy
        make_term({0, 1}, Rational(-1)), // -y
        make_term({0, 0}, Rational(-1))  // -1
    };
    MultiPoly poly(poly_terms, vars);

    // At y=1: x^2 + x - 1 - 1 = x^2 + x - 2 = (x+2)(x-1)
    Polynomial<Rational> f1({Rational(2), Rational(1)}, "x");  // x + 2
    Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x"); // x - 1

    std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

    auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(1), 1);

    EXPECT_TRUE((lifted.size() == 2)) << "non-zero eval: lifted has 2 factors";

    if (lifted.size() == 2) {
        MultiPoly product = lifted[0] * lifted[1];
        EXPECT_TRUE((product == poly)) << "non-zero eval: product == poly";
    }
}

TEST(HenselLiftSmoke, HenselLiftRestoresNamedCoordinates) {
    const std::vector<std::string> vars = {"lift", "unused", "main"};
    const MultiPoly poly({make_term({0, 0, 2}, Rational(1)),
                          make_term({1, 0, 1}, Rational(2)),
                          make_term({1, 0, 0}, Rational(-2)),
                          make_term({0, 0, 0}, Rational(-1))},
                         vars);
    const std::vector<Polynomial<Rational>> base = {
        Polynomial<Rational>({Rational(3), Rational(1)}, "main"),
        Polynomial<Rational>({Rational(-1), Rational(1)}, "main")};
    const auto lifted = multivariate_hensel_lift(poly, base, "lift", Rational(1), 1);
    EXPECT_TRUE((lifted.size() == 2)) << "two lifted factors retain the base factorization";
    if (lifted.size() != 2) {
        return;
    }
    EXPECT_TRUE((lifted[0] * lifted[1] == poly)) << "lifting reconstructs all coefficients in the original ordered ring";
    for (size_t i = 0; i < lifted.size(); ++i) {
        EXPECT_TRUE((lifted[i].variables() == vars)) << "lifting retains the unused coordinate";
        EXPECT_TRUE((lifted[i].eval("lift", Rational(1)).to_univariate() == base[i])) << "each lifted factor specializes to its original named base factor";
    }
}

TEST(HenselLiftSmoke, DiophantineConstantCofactorUsesActiveRing) {
    const std::vector<std::string> vars = {"lift", "main"};
    const std::vector<MultiPoly> cofactors = {
        MultiPoly(Rational(2), vars),
        MultiPoly({make_term({0, 1}, Rational(1)),
                   make_term({0, 0}, Rational(1))},
                  vars)};
    const MultiPoly target({make_term({0, 2}, Rational(1)),
                            make_term({0, 0}, Rational(3))},
                           vars);
    const auto corrections = multivariate_diophantine(
        cofactors, target, "lift", Rational(0), 1);
    EXPECT_TRUE((corrections.size() == 2)) << "both cofactor corrections are returned";
    if (corrections.size() != 2) {
        return;
    }
    EXPECT_TRUE((corrections[0] * cofactors[0] + corrections[1] * cofactors[1] == target)) << "corrections satisfy the exact Bezout equation in the original ring";
}

TEST(HenselLiftSmoke, HenselLiftRequiresTheSpecializationUnit) {
    const std::vector<std::string> vars = {"lift", "main"};
    const MultiPoly first({make_term({0, 1}, Rational(1)), make_term({0, 0}, Rational(-3))}, vars);
    const MultiPoly second({make_term({1, 0}, Rational(2)), make_term({0, 1}, Rational(-1)),
                            make_term({0, 0}, Rational(2))},
                           vars);
    const MultiPoly poly = first * second;
    std::vector<Polynomial<Rational>> base = {
        Polynomial<Rational>({Rational(-3), Rational(1)}, "main"),
        Polynomial<Rational>({Rational(-2), Rational(1)}, "main")};
    EXPECT_TRUE((multivariate_hensel_lift(poly, base, "lift", Rational(0), 1).empty())) << "a missing unit cannot be repaired by positive-order Hensel corrections";
    base[1] = Polynomial<Rational>({Rational(2), Rational(-1)}, "main");
    const auto lifted = multivariate_hensel_lift(poly, base, "lift", Rational(0), 1);
    EXPECT_TRUE((lifted.size() == 2)) << "the signed specialization lifts both factors";
    if (lifted.size() != 2) {
        return;
    }
    EXPECT_TRUE((lifted[0] == first && lifted[1] == second)) << "exact named-coordinate corrections recover the signed factors";
    EXPECT_TRUE((lifted[0] * lifted[1] == poly)) << "the lifted product equals the original polynomial exactly";
}

TEST(HenselLiftSmoke, DiophantineRejectsNondivisibleTarget) {
    const std::vector<std::string> vars = {"lift", "main"};
    const MultiPoly cofactor({make_term({0, 1}, Rational(1)), make_term({0, 0}, Rational(1))}, vars);
    const auto corrections = multivariate_diophantine(
        {cofactor}, MultiPoly(Rational(1), vars), "lift", Rational(0), 1);
    EXPECT_TRUE((corrections.empty())) << "a nonzero division remainder cannot produce a correction";
}

TEST(HenselLiftSmoke, SignedConstantFactorsEmbedInEmptyRing) {
    const std::vector<std::string> vars;
    const MultiPoly poly(Rational(-6), vars);
    const std::vector<Polynomial<Rational>> base = {
        Polynomial<Rational>({Rational(-2)}, "main"),
        Polynomial<Rational>({Rational(3)}, "main")};

    const auto lifted = multivariate_hensel_lift(
        poly, base, "lift", Rational(0), 1);

    ASSERT_EQ(lifted.size(), base.size());
    EXPECT_TRUE(lifted[0] == MultiPoly(Rational(-2), vars));
    EXPECT_TRUE(lifted[1] == MultiPoly(Rational(3), vars));
    EXPECT_TRUE(lifted[0] * lifted[1] == poly);
    for (const auto& factor : lifted) {
        EXPECT_EQ(factor.variables(), vars);
    }
}

TEST(HenselLiftSmoke, ActiveVariableMissingFromTargetRingIsRejected) {
    const MultiPoly poly(Rational(1), std::vector<std::string>{});
    const std::vector<Polynomial<Rational>> base = {
        Polynomial<Rational>({Rational(0), Rational(1)}, "main"),
        Polynomial<Rational>({Rational(1)}, "main")};

    EXPECT_THROW(
        multivariate_hensel_lift(poly, base, "lift", Rational(0), 1),
        std::invalid_argument);
}
