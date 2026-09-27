#include "test_multivariate_diophantine_support.hpp"

TEST(MultivariateDiophantineSolutions, UnitTwoFactorF1XF2X1Target1KnownSolutionS11S21) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "two-factor x,x+1: solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "(-1)*x + (1)*(x+1) == 1";
}

TEST(MultivariateDiophantineSolutions, UnitTwoFactorF1X1F2X1TargetXNonTrivialTarget) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({1}, Rational(1))};
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "non-trivial target x: solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x+1) + s2*(x-1) == x (non-trivial target)";
}

TEST(MultivariateDiophantineSolutions, UnitSingleFactorF1X1TargetX1SolutionS11) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 1)) << "single factor: solution has 1 component";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x+1) == x+1, so s1 = 1";

    MultiPoly expected_s1(Rational(1), vars);
    EXPECT_TRUE((solution[0] == expected_s1)) << "single factor: s1 == 1";
}

TEST(MultivariateDiophantineSolutions, UnitThreeFactorF1XF2X1F3X1Target1VerifySum) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> f3_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly f3(f3_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2, f3};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 3)) << "three-factor unit: solution has 3 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*x + s2*(x+1) + s3*(x-1) == 1 (unit test)";

    for (size_t i = 0; i < solution.size(); ++i) {
        int deg_si = solution[i].degree("x");
        EXPECT_TRUE((deg_si < 2)) << "three-factor unit: deg(s_i) < 2";
    }
}

TEST(MultivariateDiophantineSolutions, UnitTwoFactorF1X1F2X1Target1KnownSolutionS112S212) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "unit two-factor (x+1,x-1): 2 solutions";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "unit two-factor (x+1,x-1): s1*(x+1) + s2*(x-1) == 1";

    EXPECT_TRUE((solution[0].is_constant())) << "s1 is constant for (x+1,x-1) target=1";
    EXPECT_TRUE((solution[1].is_constant())) << "s2 is constant for (x+1,x-1) target=1";

    MultiPoly expected_s1(Rational(1, 2), vars);
    MultiPoly expected_s2(Rational(-1, 2), vars);
    EXPECT_TRUE((solution[0] == expected_s1)) << "s1 == 1/2";
    EXPECT_TRUE((solution[1] == expected_s2)) << "s2 == -1/2";
}

TEST(MultivariateDiophantineSolutions, UnitTwoFactorF1XF2X1Target1KnownSolutionS11S21Contract) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "unit two-factor (x,x+1): 2 solutions";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "unit two-factor (x,x+1): s1*x + s2*(x+1) == 1";

    EXPECT_TRUE((solution[0].is_constant())) << "s1 is constant for (x,x+1) target=1";
    EXPECT_TRUE((solution[1].is_constant())) << "s2 is constant for (x,x+1) target=1";

    MultiPoly expected_s1(Rational(-1), vars);
    MultiPoly expected_s2(Rational(1), vars);
    EXPECT_TRUE((solution[0] == expected_s1)) << "s1 == -1";
    EXPECT_TRUE((solution[1] == expected_s2)) << "s2 == 1";
}

TEST(MultivariateDiophantineSolutions, UnitThreeFactorF1XF2X1F3X1Target1) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> f3_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly f3(f3_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2, f3};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 3)) << "unit three-factor (x,x+1,x-1): 3 solutions";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "unit three-factor: s1*x + s2*(x+1) + s3*(x-1) == 1";

    /**
     * @brief 验证三次乘积 product = x(x+1)(x-1) = x^3-x 对应的解次数界。
     * deg(s_1) < deg((x+1)(x-1)) = 2；
     * deg(s_2) < deg(x(x-1)) = 2；deg(s₃) < deg(x(x+1)) = 2。
     */
    EXPECT_TRUE((solution[0].degree("x") < 2)) << "three-factor: deg(s1) < 2";
    EXPECT_TRUE((solution[1].degree("x") < 2)) << "three-factor: deg(s2) < 2";
    EXPECT_TRUE((solution[2].degree("x") < 2)) << "three-factor: deg(s3) < 2";
}

TEST(MultivariateDiophantineSolutions, UnitTwoFactorF1X1F2X1TargetXNonTrivialTargetContract) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({1}, Rational(1))};
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "unit non-trivial target: 2 solutions";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "unit non-trivial target: s1*(x+1) + s2*(x-1) == x";

    EXPECT_TRUE((solution[0].degree("x") < 1)) << "non-trivial target: deg(s1) < 1";
    EXPECT_TRUE((solution[1].degree("x") < 1)) << "non-trivial target: deg(s2) < 1";

    MultiPoly expected_s1(Rational(1, 2), vars);
    MultiPoly expected_s2(Rational(1, 2), vars);
    EXPECT_TRUE((solution[0] == expected_s1)) << "non-trivial target: s1 == 1/2";
    EXPECT_TRUE((solution[1] == expected_s2)) << "non-trivial target: s2 == 1/2";
}

TEST(MultivariateDiophantineSolutions, UnitDegreeBoundTruncationSolutionsRespectDegreeBound) {
    /**
     * @brief 对 f_1 = x^2+1、f_2 = x+1、target = x^2 验证 degree_bound = 3 的截断。
     * product = (x^2+1)(x+1) = x^3+x^2+x+1，次数为 3，解截断到次数小于 3。
     * deg(s_1) < deg(product/f_1) = deg(x+1) = 1；
     * deg(s_2) < deg(product/f_2) = deg(x^2+1) = 2。
     */
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({2}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({2}, Rational(1))};
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 3);

    EXPECT_TRUE((solution.size() == 2)) << "degree bound: 2 solutions";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "degree bound: s1*(x^2+1) + s2*(x+1) == x^2";

    int deg_s1 = solution[0].degree("x");
    int deg_s2 = solution[1].degree("x");
    EXPECT_TRUE((deg_s1 < 1)) << "degree bound: deg(s1) < deg(product/f1) = 1";
    EXPECT_TRUE((deg_s2 < 2)) << "degree bound: deg(s2) < deg(product/f2) = 2";

    std::vector<MultiPoly> solution_tight = multivariate_diophantine(
        factors, target, "x", Rational(0), 1);
    EXPECT_TRUE((solution_tight.size() == 2)) << "tight degree bound: 2 solutions returned";
    EXPECT_TRUE((solution_tight[0].degree("x") < 1)) << "tight degree bound: s1 is constant";
    EXPECT_TRUE((solution_tight[1].degree("x") < 1)) << "tight degree bound: s2 is constant";
}

TEST(MultivariateDiophantineSolutions, ConstantFactorsInEmptyRing) {
    const std::vector<std::string> vars;
    const std::vector<MultiPoly> factors = {
        MultiPoly(Rational(2), vars),
        MultiPoly(Rational(3), vars)};

    for (int value : {1, 0}) {
        SCOPED_TRACE(value);
        const MultiPoly target(Rational(value), vars);
        const auto solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 1);

        ASSERT_EQ(solution.size(), factors.size());
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target));
        for (const auto& coefficient : solution) {
            EXPECT_EQ(coefficient.variables(), vars);
        }
    }
}
