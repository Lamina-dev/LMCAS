/**
 * @file test_hensel_lift_smoke.cpp
 * @brief 多元 Hensel 提升冒烟测试。
 */
#include "test_common.hpp"
#include "multivariate_factor.hpp"

using namespace LMCAS;

static MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff) {
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}

int main() {
    TEST_CASE("Hensel lift: (x+y+1)(x-1) from factors at y=0");
    {
        std::vector<std::string> vars = {"x", "y"};

        // poly = (x+y+1)(x-1) = x^2 + xy - y - 1
        std::vector<MultiPoly::Term> poly_terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({1, 1}, Rational(1)),   // xy
            make_term({0, 1}, Rational(-1)),  // -y
            make_term({0, 0}, Rational(-1))   // -1
        };
        MultiPoly poly(poly_terms, vars);

        // univariate factors at y=0: x+1, x-1
        Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");   // x + 1
        Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x");  // x - 1

        std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 1);

        EXPECT_TRUE(lifted.size() == 2, "lifted has 2 factors");

        // Check product equals poly
        if (lifted.size() == 2) {
            MultiPoly product = lifted[0] * lifted[1];
            EXPECT_TRUE(product == poly, "product of lifted factors == poly");
        }
    }

    TEST_CASE("Hensel lift: (x+1)(x-1) no y dependence, degree_bound=2");
    {
        std::vector<std::string> vars = {"x", "y"};

        // poly = x^2 - 1 (no y terms)
        std::vector<MultiPoly::Term> poly_terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 0}, Rational(-1))   // -1
        };
        MultiPoly poly(poly_terms, vars);

        // univariate factors at y=0: x+1, x-1
        Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");   // x + 1
        Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x");  // x - 1

        std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

        EXPECT_TRUE(lifted.size() == 2, "lifted has 2 factors");

        if (lifted.size() == 2) {
            MultiPoly product = lifted[0] * lifted[1];
            EXPECT_TRUE(product == poly, "product of lifted factors == poly (no y)");
        }
    }

    TEST_CASE("Hensel lift: single factor returns poly itself");
    {
        std::vector<std::string> vars = {"x", "y"};

        std::vector<MultiPoly::Term> poly_terms = {
            make_term({2, 1}, Rational(1)),   // x^2*y
            make_term({1, 0}, Rational(1))    // x
        };
        MultiPoly poly(poly_terms, vars);

        Polynomial<Rational> f1({Rational(0), Rational(1)}, "x");  // x

        std::vector<Polynomial<Rational>> uni_factors = {f1};

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

        EXPECT_TRUE(lifted.size() == 1, "single factor: 1 result");
        if (lifted.size() == 1) {
            EXPECT_TRUE(lifted[0] == poly, "single factor returns poly itself");
        }
    }

    TEST_CASE("Hensel lift: empty factors returns empty");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(1), vars);
        std::vector<Polynomial<Rational>> uni_factors;

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);
        EXPECT_TRUE(lifted.empty(), "empty factors returns empty");
    }

    TEST_CASE("Hensel lift: (x+2y+1)(x-y+1) from factors at y=0");
    {
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
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1)),  // -y^2
            make_term({0, 1}, Rational(-2)),  // -2y
            make_term({0, 0}, Rational(-1))   // -1
        };
        MultiPoly poly(poly_terms, vars);

        // univariate factors at y=0: x+1, x-1
        Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");   // x + 1
        Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x");  // x - 1

        std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(0), 2);

        EXPECT_TRUE(lifted.size() == 2, "lifted has 2 factors");

        if (lifted.size() == 2) {
            MultiPoly product = lifted[0] * lifted[1];
            EXPECT_TRUE(product == poly, "product of lifted factors == poly (degree 2 in y)");
        }
    }

    TEST_CASE("Hensel lift: non-zero evaluation point (x+y+1)(x-1) at y=1");
    {
        std::vector<std::string> vars = {"x", "y"};

        // poly = (x+y+1)(x-1) = x^2 + xy - y - 1
        std::vector<MultiPoly::Term> poly_terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({1, 1}, Rational(1)),   // xy
            make_term({0, 1}, Rational(-1)),  // -y
            make_term({0, 0}, Rational(-1))   // -1
        };
        MultiPoly poly(poly_terms, vars);

        // At y=1: x^2 + x - 1 - 1 = x^2 + x - 2 = (x+2)(x-1)
        Polynomial<Rational> f1({Rational(2), Rational(1)}, "x");   // x + 2
        Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x");  // x - 1

        std::vector<Polynomial<Rational>> uni_factors = {f1, f2};

        auto lifted = multivariate_hensel_lift(poly, uni_factors, "y", Rational(1), 1);

        EXPECT_TRUE(lifted.size() == 2, "non-zero eval: lifted has 2 factors");

        if (lifted.size() == 2) {
            MultiPoly product = lifted[0] * lifted[1];
            EXPECT_TRUE(product == poly, "non-zero eval: product == poly");
        }
    }

    return TEST_REPORT();
}
