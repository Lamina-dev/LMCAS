#include "test_multivariate_diophantine_support.hpp"

TEST(MultivariateDiophantine, TwoCoprimeLinearFactorsTarget1) {
    // f_1 = x+1, f_2 = x-1, target = 1
    // 求解 s_1*(x+1) + s_2*(x-1) = 1
    // 已知解:s_1 = 1/2, s_2 = -1/2(因为 (1/2)(x+1) + (-1/2)(x-1) = 1)
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x+1) + s2*(x-1) == 1";
}

TEST(MultivariateDiophantine, TwoCoprimeLinearFactorsTargetX) {
    // f_1 = x+1, f_2 = x-1, target = x
    // 求解 s_1*(x+1) + s_2*(x-1) = x
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({1}, Rational(1)) // x
    };
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x+1) + s2*(x-1) == x";
}

TEST(MultivariateDiophantine, ThreeCoprimeFactorsTarget1) {
    // f_1 = x, f_2 = x+1, f₃ = x-1, target = 1
    // 求解 s_1*x + s_2*(x+1) + s₃*(x-1) = 1
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)) // x
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> f3_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly f3(f3_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2, f3};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 3)) << "solution has 3 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*x + s2*(x+1) + s3*(x-1) == 1";
}

TEST(MultivariateDiophantine, FactorsWithHigherDegreeTarget1) {
    // f_1 = x^2+1, f_2 = x+1, target = 1
    // 求解 s_1*(x^2+1) + s_2*(x+1) = 1
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({2}, Rational(1)), // x^2
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x^2+1) + s2*(x+1) == 1";
}

TEST(MultivariateDiophantine, NonTrivialTargetTwoCoprimeFactors) {
    // f_1 = x+1, f_2 = x+2, target = x+3
    // 求解 s_1*(x+1) + s_2*(x+2) = x+3
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(2))  // +2
    };
    MultiPoly f2(f2_terms, vars);

    std::vector<MultiPoly::Term> target_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(3))  // +3
    };
    MultiPoly target(target_terms, vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x+1) + s2*(x+2) == x+3";
}

TEST(MultivariateDiophantine, Target0GivesAllZeroSolution) {
    // f_1 = x+1, f_2 = x-1, target = 0
    // 求解 s_1*(x+1) + s_2*(x-1) = 0
    // 平凡解:s_1 = 0, s_2 = 0
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)), // x
        make_term({0}, Rational(-1)) // -1
    };
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(0), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "zero target: solution has 2 components";

    // All solution components should be zero
    bool all_zero = true;
    for (const auto &s : solution) {
        if (!s.is_zero()) {
            all_zero = false;
            break;
        }
    }
    EXPECT_TRUE((all_zero)) << "zero target: all solution components are zero";

    // Also verify via the general correctness check
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "zero target: s1*(x+1) + s2*(x-1) == 0";
}

TEST(MultivariateDiophantine, TwoQuadraticCoprimeFactorsTarget1) {
    // f_1 = x^2+1, f_2 = x^2+x+1, target = 1
    // These are coprime since gcd(x^2+1, x^2+x+1) = 1
    // 求解 s_1*(x^2+1) + s_2*(x^2+x+1) = 1
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({2}, Rational(1)), // x^2
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({2}, Rational(1)), // x^2
        make_term({1}, Rational(1)), // +x
        make_term({0}, Rational(1))  // +1
    };
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "two quadratic: solution has 2 components";
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "s1*(x^2+1) + s2*(x^2+x+1) == 1";
}

TEST(MultivariateDiophantine, DiophantineDegreeConstraintsDegSIDegProductFI) {
    // f_1 = x+1, f_2 = x-1, product = (x+1)(x-1) = x^2-1
    // deg(s_1) < deg(product/f_1) = deg(x-1) = 1 -> s_1 is constant
    // deg(s_2) < deg(product/f_2) = deg(x+1) = 1 -> s_2 is constant
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

    EXPECT_TRUE((solution.size() == 2)) << "degree constraint: 2 solutions";

    // product/f_1 = f_2, deg(f_2) = 1 -> deg(s_1) < 1
    int deg_s1 = solution[0].degree("x");
    EXPECT_TRUE((deg_s1 < 1)) << "deg(s1) < deg(product/f1) = 1";

    // product/f_2 = f_1, deg(f_1) = 1 -> deg(s_2) < 1
    int deg_s2 = solution[1].degree("x");
    EXPECT_TRUE((deg_s2 < 1)) << "deg(s2) < deg(product/f2) = 1";

    // Verify correctness still holds
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "degree constraint: sum still equals target";
}

TEST(MultivariateDiophantine, DiophantineDegreeConstraintsHigherDegreeFactors) {
    // f_1 = x^2+1, f_2 = x+1, product = (x^2+1)(x+1)
    // deg(s_1) < deg(product/f_1) = deg(x+1) = 1 -> s_1 is constant
    // deg(s_2) < deg(product/f_2) = deg(x^2+1) = 2 -> deg(s_2) <= 1
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> f1_terms = {
        make_term({2}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f1(f1_terms, vars);

    std::vector<MultiPoly::Term> f2_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly f2(f2_terms, vars);

    MultiPoly target(Rational(1), vars);

    std::vector<MultiPoly> factors = {f1, f2};
    std::vector<MultiPoly> solution = multivariate_diophantine(
        factors, target, "x", Rational(0), 10);

    EXPECT_TRUE((solution.size() == 2)) << "higher deg constraint: 2 solutions";

    // deg(s_1) < deg(product/f_1) = deg(x+1) = 1
    int deg_s1 = solution[0].degree("x");
    EXPECT_TRUE((deg_s1 < 1)) << "deg(s1) < deg(x+1) = 1";

    // deg(s_2) < deg(product/f_2) = deg(x^2+1) = 2
    int deg_s2 = solution[1].degree("x");
    EXPECT_TRUE((deg_s2 < 2)) << "deg(s2) < deg(x^2+1) = 2";

    // Verify correctness
    EXPECT_TRUE((verify_diophantine_solution(factors, solution, target))) << "higher deg constraint: sum equals target";
}
