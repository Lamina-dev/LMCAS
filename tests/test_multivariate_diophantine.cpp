
#include "test_common.hpp"
#include "multivariate_factor.hpp"
#include "rapidcheck/rapidcheck.h"

using namespace lamina;

/**
 * @brief 辅助函数：构造 MultiPoly 项
 */
static MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff)
{
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}

/**
 * @brief 辅助函数：验证丢番图方程解的正确性
 *
 * 检查 s₁*f₁ + ... + sᵣ*fᵣ == target（截断到 degree < degree_bound in var）
 */
static bool verify_diophantine_solution(
    const std::vector<MultiPoly>& factors,
    const std::vector<MultiPoly>& solution,
    const MultiPoly& target,
    const std::string& var,
    int degree_bound)
{
    if (solution.size() != factors.size()) return false;

    // 计算 sum = s₁*f₁ + s₂*f₂ + ... + sᵣ*fᵣ
    const auto& vars = factors[0].variables();
    MultiPoly sum(Rational(0), vars);
    for (size_t i = 0; i < factors.size(); ++i) {
        sum = sum + solution[i] * factors[i];
    }

    /// 按 var 将 sum 截断到 degree < degree_bound。
    /// 当前用例的 degree_bound 覆盖全部项，因此直接比较 sum 与 target。
    return sum == target;
}

int main()
{

    TEST_CASE("Two coprime linear factors, target = 1");
    {
        // f₁ = x+1, f₂ = x-1, target = 1
        // 求解 s₁*(x+1) + s₂*(x-1) = 1
        // 已知解：s₁ = 1/2, s₂ = -1/2（因为 (1/2)(x+1) + (-1/2)(x-1) = 1）
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x+1) + s2*(x-1) == 1");
    }

    TEST_CASE("Two coprime linear factors, target = x");
    {
        // f₁ = x+1, f₂ = x-1, target = x
        // 求解 s₁*(x+1) + s₂*(x-1) = x
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x+1) + s2*(x-1) == x");
    }

    TEST_CASE("Three coprime factors, target = 1");
    {
        // f₁ = x, f₂ = x+1, f₃ = x-1, target = 1
        // 求解 s₁*x + s₂*(x+1) + s₃*(x-1) = 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> f3_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f3(f3_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2, f3};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 3, "solution has 3 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*x + s2*(x+1) + s3*(x-1) == 1");
    }

    TEST_CASE("Factors with higher degree, target = 1");
    {
        // f₁ = x²+1, f₂ = x+1, target = 1
        // 求解 s₁*(x²+1) + s₂*(x+1) = 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({2}, Rational(1)),  // x²
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x^2+1) + s2*(x+1) == 1");
    }

    TEST_CASE("Non-trivial target, two coprime factors");
    {
        // f₁ = x+1, f₂ = x+2, target = x+3
        // 求解 s₁*(x+1) + s₂*(x+2) = x+3
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(2))   // +2
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(3))   // +3
        };
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x+1) + s2*(x+2) == x+3");
    }

    TEST_CASE("Target = 0 gives all-zero solution");
    {
        // f₁ = x+1, f₂ = x-1, target = 0
        // 求解 s₁*(x+1) + s₂*(x-1) = 0
        // 平凡解：s₁ = 0, s₂ = 0
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(0), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "zero target: solution has 2 components");

        // All solution components should be zero
        bool all_zero = true;
        for (const auto& s : solution) {
            if (!s.is_zero()) {
                all_zero = false;
                break;
            }
        }
        EXPECT_TRUE(all_zero, "zero target: all solution components are zero");

        // Also verify via the general correctness check
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "zero target: s1*(x+1) + s2*(x-1) == 0");
    }

    TEST_CASE("Two quadratic coprime factors, target = 1");
    {
        // f₁ = x²+1, f₂ = x²+x+1, target = 1
        // These are coprime since gcd(x²+1, x²+x+1) = 1
        // 求解 s₁*(x²+1) + s₂*(x²+x+1) = 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({2}, Rational(1)),  // x²
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({2}, Rational(1)),  // x²
            make_term({1}, Rational(1)),  // +x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "two quadratic: solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x^2+1) + s2*(x^2+x+1) == 1");
    }


    TEST_CASE("Diophantine degree constraints: deg(s_i) < deg(product/f_i)");
    {
        // f₁ = x+1, f₂ = x-1, product = (x+1)(x-1) = x²-1
        // deg(s₁) < deg(product/f₁) = deg(x-1) = 1 → s₁ is constant
        // deg(s₂) < deg(product/f₂) = deg(x+1) = 1 → s₂ is constant
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),
            make_term({0}, Rational(1))
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),
            make_term({0}, Rational(-1))
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "degree constraint: 2 solutions");

        // product/f₁ = f₂, deg(f₂) = 1 → deg(s₁) < 1
        int deg_s1 = solution[0].degree("x");
        EXPECT_TRUE(deg_s1 < 1, "deg(s1) < deg(product/f1) = 1");

        // product/f₂ = f₁, deg(f₁) = 1 → deg(s₂) < 1
        int deg_s2 = solution[1].degree("x");
        EXPECT_TRUE(deg_s2 < 1, "deg(s2) < deg(product/f2) = 1");

        // Verify correctness still holds
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "degree constraint: sum still equals target");
    }

    TEST_CASE("Diophantine degree constraints: higher degree factors");
    {
        // f₁ = x²+1, f₂ = x+1, product = (x²+1)(x+1)
        // deg(s₁) < deg(product/f₁) = deg(x+1) = 1 → s₁ is constant
        // deg(s₂) < deg(product/f₂) = deg(x²+1) = 2 → deg(s₂) ≤ 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({2}, Rational(1)),
            make_term({0}, Rational(1))
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),
            make_term({0}, Rational(1))
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "higher deg constraint: 2 solutions");

        // deg(s₁) < deg(product/f₁) = deg(x+1) = 1
        int deg_s1 = solution[0].degree("x");
        EXPECT_TRUE(deg_s1 < 1, "deg(s1) < deg(x+1) = 1");

        // deg(s₂) < deg(product/f₂) = deg(x²+1) = 2
        int deg_s2 = solution[1].degree("x");
        EXPECT_TRUE(deg_s2 < 2, "deg(s2) < deg(x^2+1) = 2");

        // Verify correctness
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "higher deg constraint: sum equals target");
    }


    TEST_CASE("Unit: Two-factor f1=x, f2=x+1, target=1, known solution s1=-1, s2=1");
    {
        // f₁ = x, f₂ = x+1, target = 1
        // Known: s₁ = -1, s₂ = 1 (since (-1)(x) + (1)(x+1) = -x + x + 1 = 1)
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "two-factor x,x+1: solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "(-1)*x + (1)*(x+1) == 1");
    }

    TEST_CASE("Unit: Two-factor f1=x+1, f2=x-1, target=x (non-trivial target)");
    {
        // f₁ = x+1, f₂ = x-1, target = x
        // Verify s₁*(x+1) + s₂*(x-1) = x
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "non-trivial target x: solution has 2 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x+1) + s2*(x-1) == x (non-trivial target)");
    }

    TEST_CASE("Unit: Single factor f1=x+1, target=x+1, solution s1=1");
    {
        // f₁ = x+1, target = x+1
        // Trivially: s₁ = 1 (since 1*(x+1) = x+1)
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 1, "single factor: solution has 1 component");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*(x+1) == x+1, so s1 = 1");

        // Additionally verify s₁ is exactly the constant 1
        MultiPoly expected_s1(Rational(1), vars);
        EXPECT_TRUE(solution[0] == expected_s1, "single factor: s1 == 1");
    }

    TEST_CASE("Unit: Three-factor f1=x, f2=x+1, f3=x-1, target=1, verify sum");
    {
        // f₁ = x, f₂ = x+1, f₃ = x-1, target = 1
        // Verify s₁*x + s₂*(x+1) + s₃*(x-1) = 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> f3_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f3(f3_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2, f3};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 3, "three-factor unit: solution has 3 components");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "s1*x + s2*(x+1) + s3*(x-1) == 1 (unit test)");

        // Verify degree constraints: each s_i should be constant (degree 0)
        // since product = x*(x+1)*(x-1) = x³-x, deg(product/f_i) ≤ 2
        // and deg(s_i) < deg(product/f_i)
        for (size_t i = 0; i < solution.size(); ++i) {
            int deg_si = solution[i].degree("x");
            EXPECT_TRUE(deg_si < 2, "three-factor unit: deg(s_i) < 2");
        }
    }


    TEST_CASE("Unit: Two-factor f1=x+1, f2=x-1, target=1 — known solution s1=1/2, s2=-1/2");
    {
        // f₁ = x+1, f₂ = x-1, target = 1
        // Expected: s₁ = 1/2, s₂ = -1/2 (since (1/2)(x+1) + (-1/2)(x-1) = 1)
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "unit two-factor (x+1,x-1): 2 solutions");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "unit two-factor (x+1,x-1): s1*(x+1) + s2*(x-1) == 1");

        // Verify specific expected values: s₁ = 1/2, s₂ = -1/2
        // Both should be constants
        EXPECT_TRUE(solution[0].is_constant(), "s1 is constant for (x+1,x-1) target=1");
        EXPECT_TRUE(solution[1].is_constant(), "s2 is constant for (x+1,x-1) target=1");

        // Check s₁ = 1/2 by evaluating at x=0 (since it's constant, eval doesn't matter)
        MultiPoly expected_s1(Rational(1, 2), vars);
        MultiPoly expected_s2(Rational(-1, 2), vars);
        EXPECT_TRUE(solution[0] == expected_s1, "s1 == 1/2");
        EXPECT_TRUE(solution[1] == expected_s2, "s2 == -1/2");
    }

    TEST_CASE("Unit: Two-factor f1=x, f2=x+1, target=1 — known solution s1=-1, s2=1");
    {
        // f₁ = x, f₂ = x+1, target = 1
        // Expected: s₁ = -1, s₂ = 1 (since -1*x + 1*(x+1) = -x + x + 1 = 1)
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "unit two-factor (x,x+1): 2 solutions");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "unit two-factor (x,x+1): s1*x + s2*(x+1) == 1");

        // Verify specific expected values: s₁ = -1, s₂ = 1
        EXPECT_TRUE(solution[0].is_constant(), "s1 is constant for (x,x+1) target=1");
        EXPECT_TRUE(solution[1].is_constant(), "s2 is constant for (x,x+1) target=1");

        MultiPoly expected_s1(Rational(-1), vars);
        MultiPoly expected_s2(Rational(1), vars);
        EXPECT_TRUE(solution[0] == expected_s1, "s1 == -1");
        EXPECT_TRUE(solution[1] == expected_s2, "s2 == 1");
    }

    TEST_CASE("Unit: Three-factor f1=x, f2=x+1, f3=x-1, target=1");
    {
        // f₁ = x, f₂ = x+1, f₃ = x-1, target = 1
        // Verify s₁*x + s₂*(x+1) + s₃*(x-1) == 1
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> f3_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f3(f3_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2, f3};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 3, "unit three-factor (x,x+1,x-1): 3 solutions");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "unit three-factor: s1*x + s2*(x+1) + s3*(x-1) == 1");

        // Degree constraints: product = x(x+1)(x-1) = x³-x, degree 3
        // deg(s₁) < deg(product/f₁) = deg((x+1)(x-1)) = 2
        // deg(s₂) < deg(product/f₂) = deg(x(x-1)) = 2
        // deg(s₃) < deg(product/f₃) = deg(x(x+1)) = 2
        EXPECT_TRUE(solution[0].degree("x") < 2, "three-factor: deg(s1) < 2");
        EXPECT_TRUE(solution[1].degree("x") < 2, "three-factor: deg(s2) < 2");
        EXPECT_TRUE(solution[2].degree("x") < 2, "three-factor: deg(s3) < 2");
    }

    TEST_CASE("Unit: Two-factor f1=x+1, f2=x-1, target=x — non-trivial target");
    {
        // f₁ = x+1, f₂ = x-1, target = x
        // Verify s₁*(x+1) + s₂*(x-1) == x
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({1}, Rational(1))   // x
        };
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        EXPECT_TRUE(solution.size() == 2, "unit non-trivial target: 2 solutions");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 10),
                    "unit non-trivial target: s1*(x+1) + s2*(x-1) == x");

        // Both solutions should be constants (deg < deg(product/fi) = 1)
        EXPECT_TRUE(solution[0].degree("x") < 1, "non-trivial target: deg(s1) < 1");
        EXPECT_TRUE(solution[1].degree("x") < 1, "non-trivial target: deg(s2) < 1");

        // Verify: s₁ = 1/2, s₂ = 1/2 (since (1/2)(x+1) + (1/2)(x-1) = x)
        MultiPoly expected_s1(Rational(1, 2), vars);
        MultiPoly expected_s2(Rational(1, 2), vars);
        EXPECT_TRUE(solution[0] == expected_s1, "non-trivial target: s1 == 1/2");
        EXPECT_TRUE(solution[1] == expected_s2, "non-trivial target: s2 == 1/2");
    }

    TEST_CASE("Unit: Degree bound truncation — solutions respect degree_bound");
    {
        // f₁ = x²+1, f₂ = x+1, target = x² (higher degree target)
        // With degree_bound = 3, solutions should be truncated to degree < 3
        // product = (x²+1)(x+1) = x³+x²+x+1, degree 3
        // deg(s₁) < deg(product/f₁) = deg(x+1) = 1
        // deg(s₂) < deg(product/f₂) = deg(x²+1) = 2
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> f1_terms = {
            make_term({2}, Rational(1)),  // x²
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // +1
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> target_terms = {
            make_term({2}, Rational(1))   // x²
        };
        MultiPoly target(target_terms, vars);

        // Use degree_bound = 3 (sufficient for this problem)
        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 3);

        EXPECT_TRUE(solution.size() == 2, "degree bound: 2 solutions");
        EXPECT_TRUE(verify_diophantine_solution(factors, solution, target, "x", 3),
                    "degree bound: s1*(x^2+1) + s2*(x+1) == x^2");

        // Verify degree constraints are respected
        int deg_s1 = solution[0].degree("x");
        int deg_s2 = solution[1].degree("x");
        EXPECT_TRUE(deg_s1 < 1, "degree bound: deg(s1) < deg(product/f1) = 1");
        EXPECT_TRUE(deg_s2 < 2, "degree bound: deg(s2) < deg(product/f2) = 2");

        // With a very small degree_bound = 1, the truncation should still produce
        // valid results within the truncated space (or at least not crash)
        std::vector<MultiPoly> solution_tight = multivariate_diophantine(
            factors, target, "x", Rational(0), 1);
        EXPECT_TRUE(solution_tight.size() == 2, "tight degree bound: 2 solutions returned");
        // With degree_bound=1, solutions are truncated to degree < 1 (constants only)
        EXPECT_TRUE(solution_tight[0].degree("x") < 1,
                    "tight degree bound: s1 is constant");
        EXPECT_TRUE(solution_tight[1].degree("x") < 1,
                    "tight degree bound: s2 is constant");
    }


    TEST_CASE("Diophantine solver correctness (random coprime pairs)");
    rc::check("For random coprime linear factor pairs, s1*f1 + s2*f2 == target", []() {
        // Generate two coprime linear factors: f1 = x + a, f2 = x + b with a != b
        std::vector<std::string> vars = {"x"};
        int a = rc::gen::inRange(-10, 10);
        int b = rc::gen::inRange(-10, 10);
        // Ensure a != b so factors are coprime
        if (a == b) b = a + 1;

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}
        };
        MultiPoly f2(f2_terms, vars);

        // Generate a random constant target
        int t = rc::gen::inRange(-10, 10);
        if (t == 0) t = 1;
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        // Verify: s1*f1 + s2*f2 == target
        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    });

    TEST_CASE("Diophantine solver correctness (random linear target)");
    rc::check("For coprime linear factors with linear target, s1*f1 + s2*f2 == target", []() {
        // Generate two coprime linear factors: f1 = x + a, f2 = x + b with a != b
        std::vector<std::string> vars = {"x"};
        int a = rc::gen::inRange(-5, 5);
        int b = rc::gen::inRange(-5, 5);
        if (a == b) b = a + 1;

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}
        };
        MultiPoly f2(f2_terms, vars);

        // Generate a random linear target: c1*x + c0
        int c1 = rc::gen::inRange(-5, 5);
        int c0 = rc::gen::inRange(-5, 5);
        if (c1 == 0 && c0 == 0) c0 = 1;

        std::vector<MultiPoly::Term> target_terms;
        if (c1 != 0) target_terms.push_back({Monomial({1}), Rational(c1)});
        if (c0 != 0) target_terms.push_back({Monomial({0}), Rational(c0)});
        if (target_terms.empty()) target_terms.push_back({Monomial({0}), Rational(1)});
        MultiPoly target(target_terms, vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        // Verify: s1*f1 + s2*f2 == target
        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    });

    TEST_CASE("Diophantine solver correctness (three coprime factors)");
    rc::check("For three coprime linear factors, s1*f1 + s2*f2 + s3*f3 == target", []() {
        // Generate three coprime linear factors: fi = x + ai with distinct ai
        std::vector<std::string> vars = {"x"};
        int a = rc::gen::inRange(-5, 5);
        int b = a + 1 + rc::gen::inRange(0, 3);
        int c = b + 1 + rc::gen::inRange(0, 3);

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}
        };
        MultiPoly f2(f2_terms, vars);

        std::vector<MultiPoly::Term> f3_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(c)}
        };
        MultiPoly f3(f3_terms, vars);

        // Target is a constant
        int t = rc::gen::inRange(1, 5);
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2, f3};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 3);

        // Verify: s1*f1 + s2*f2 + s3*f3 == target
        MultiPoly sum = solution[0] * f1 + solution[1] * f2 + solution[2] * f3;
        RC_ASSERT(sum == target);
    });

    TEST_CASE("Diophantine solver correctness (quadratic and linear coprime)");
    rc::check("For coprime quadratic+linear factors, s1*f1 + s2*f2 == target", []() {
        // f1 = x^2 + a (quadratic), f2 = x + b (linear), coprime when b^2 + a != 0
        std::vector<std::string> vars = {"x"};
        int a = rc::gen::inRange(1, 5);  // positive ensures x^2+a has no rational roots
        int b = rc::gen::inRange(-5, 5);

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({2}), Rational(1)},
            {Monomial({0}), Rational(a)}
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}
        };
        MultiPoly f2(f2_terms, vars);

        // Target is a constant
        int t = rc::gen::inRange(1, 5);
        MultiPoly target(Rational(t), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        // Verify: s1*f1 + s2*f2 == target
        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    });

    TEST_CASE("Diophantine solver correctness (degree constraint)");
    rc::check("Solution components satisfy degree constraints: deg(si) < deg(product/fi)", []() {
        // f1 = x + a, f2 = x + b with a != b
        std::vector<std::string> vars = {"x"};
        int a = rc::gen::inRange(-10, 10);
        int b = rc::gen::inRange(-10, 10);
        if (a == b) b = a + 1;

        std::vector<MultiPoly::Term> f1_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(a)}
        };
        MultiPoly f1(f1_terms, vars);

        std::vector<MultiPoly::Term> f2_terms = {
            {Monomial({1}), Rational(1)},
            {Monomial({0}), Rational(b)}
        };
        MultiPoly f2(f2_terms, vars);

        MultiPoly target(Rational(1), vars);

        std::vector<MultiPoly> factors = {f1, f2};
        std::vector<MultiPoly> solution = multivariate_diophantine(
            factors, target, "x", Rational(0), 10);

        RC_ASSERT(solution.size() == 2);

        // deg(s1) < deg(product/f1) = deg(f2) = 1 → s1 is constant
        RC_ASSERT(solution[0].degree("x") < 1);
        // deg(s2) < deg(product/f2) = deg(f1) = 1 → s2 is constant
        RC_ASSERT(solution[1].degree("x") < 1);

        // Also verify correctness
        MultiPoly sum = solution[0] * f1 + solution[1] * f2;
        RC_ASSERT(sum == target);
    });

    return TEST_REPORT();
}
