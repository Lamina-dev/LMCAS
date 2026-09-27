#include "test_multivariate_diophantine_support.hpp"

TEST(MultivariateDiophantineProperties, DiophantineSolverCorrectnessRandomCoprimePairs) {
    EXPECT_TRUE(rc::check("For random coprime linear factor pairs, s1*f1 + s2*f2 == target", []() {
        std::vector<std::string> vars = {"x"};
        int a = (*rc::gen::inRange(-10, (10) + 1));
        int b = (*rc::gen::inRange(-10, (10) + 1));
        if (a == b) {
            b = a + 1;
        }

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}};
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}};
        MultiPoly f2(f2_terms, vars);

        int t = (*rc::gen::inRange(-10, (10) + 1));
        if (t == 0) {
            t = 1;
        }
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    }));
}

TEST(MultivariateDiophantineProperties, DiophantineSolverCorrectnessRandomLinearTarget) {
    EXPECT_TRUE(rc::check("For coprime linear factors with linear target, s1*f1 + s2*f2 == target", []() {
        std::vector<std::string> vars = {"x"};
        int a = (*rc::gen::inRange(-5, (5) + 1));
        int b = (*rc::gen::inRange(-5, (5) + 1));
        if (a == b) {
            b = a + 1;
        }

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}};
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}};
        MultiPoly f2(f2_terms, vars);

        int c1 = (*rc::gen::inRange(-5, (5) + 1));
        int c0 = (*rc::gen::inRange(-5, (5) + 1));
        if (c1 == 0 && c0 == 0) {
            c0 = 1;
        }

        std::vector<MultiPoly::Term> target_terms;
        if (c1 != 0)
            target_terms.push_back({Monomial({1}), Rational(c1)});
        if (c0 != 0)
            target_terms.push_back({Monomial({0}), Rational(c0)});
        if (target_terms.empty())
            target_terms.push_back({Monomial({0}), Rational(1)});
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    }));
}

TEST(MultivariateDiophantineProperties, DiophantineSolverCorrectnessThreeCoprimeFactors) {
    EXPECT_TRUE(rc::check("For three coprime linear factors, s1*f1 + s2*f2 + s3*f3 == target", []() {
        std::vector<std::string> vars = {"x"};
        int a = (*rc::gen::inRange(-5, (5) + 1));
        int b = a + 1 + (*rc::gen::inRange(0, (3) + 1));
        int c = b + 1 + (*rc::gen::inRange(0, (3) + 1));

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}};
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}};
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> f3_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(c)}};
        MultiPoly f3(f3_terms, vars);

        int t = (*rc::gen::inRange(1, (5) + 1));
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2, f3};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 3);

        MultiPoly sum = solution[0] * f1 + solution[1] * f2 + solution[2] * f3;
        RC_ASSERT(sum == target);
    }));
}

TEST(MultivariateDiophantineProperties, DiophantineSolverCorrectnessQuadraticAndLinearCoprime) {
    EXPECT_TRUE(rc::check("For coprime quadratic+linear factors, s1*f1 + s2*f2 == target", []() {
        std::vector<std::string> vars = {"x"};
        int a = (*rc::gen::inRange(1, (5) + 1)); /**< a > 0 使 x^2+a 的有理根集合为空。 */
        int b = (*rc::gen::inRange(-5, (5) + 1));

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({2}), Rational(1)},
            {Monomial({0}), Rational(a)}};
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}};
        MultiPoly f2(f2_terms, vars);

        int t = (*rc::gen::inRange(1, (5) + 1));
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    }));
}

TEST(MultivariateDiophantineProperties, DiophantineSolverCorrectnessDegreeConstraint) {
    EXPECT_TRUE(rc::check("Solution components satisfy degree constraints: deg(si) < deg(product/fi)", []() {
        std::vector<std::string> vars = {"x"};
        int a = (*rc::gen::inRange(-10, (10) + 1));
        int b = (*rc::gen::inRange(-10, (10) + 1));
        if (a == b) {
            b = a + 1;
        }

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}};
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}};
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        RC_ASSERT(solution[0].degree("x") < 1);
        RC_ASSERT(solution[1].degree("x") < 1);

        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    }));
}
