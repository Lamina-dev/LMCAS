/**
 * @file test_multivariate_sqfree.cpp
 * @brief 多元无平方因子分解的属性测试与单元测试。
 *
 * Property 7: Square-free decomposition correctness
 * - product of components^multiplicity equals original (up to constant)
 * - each component is square-free: gcd(fᵢ, ∂fᵢ/∂main_var) is constant
 *
 * **Validates: Requirements 3.1, 3.2**
 */

#include "test_common.hpp"
#include "multivariate_factor.hpp"

using namespace lamina;

/**
 * @brief 辅助函数：构造 MultiPoly 项
 */
static MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff)
{
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}

/**
 * @brief 辅助函数：计算多元多项式关于指定变量的形式导数
 *
 * 用于测试中验证各分量是否无平方。
 */
static MultiPoly test_formal_derivative(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return poly;

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) {
        return MultiPoly(Rational(0), vars);
    }

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp == 0) continue;

        Rational new_coeff = term.second * Rational(exp);
        Monomial new_mono = mono;
        new_mono[var_idx] = exp - 1;
        result_terms.emplace_back(std::move(new_mono), std::move(new_coeff));
    }

    if (result_terms.empty()) return MultiPoly(Rational(0), vars);
    return MultiPoly(std::move(result_terms), vars);
}

/**
 * @brief 辅助函数：检查多项式是否无平方
 *
 * 无平方条件：gcd(f, ∂f/∂main_var) 为常数。
 */
static bool is_square_free(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero() || poly.is_constant()) return true;

    MultiPoly deriv = test_formal_derivative(poly, main_var);
    if (deriv.is_zero()) return true;

    MultiPoly g = multivariate_gcd(poly, deriv);
    return g.is_constant();
}

/**
 * @brief 辅助函数：计算 base^exp（MultiPoly 幂运算）
 */
static MultiPoly poly_pow(const MultiPoly& base, int exp)
{
    if (exp == 0) {
        return MultiPoly(Rational(1), base.variables());
    }
    MultiPoly result = base;
    for (int i = 1; i < exp; ++i) {
        result = result * base;
    }
    return result;
}

/**
 * @brief 辅助函数：验证无平方因子分解的两个核心性质
 *
 * 1. f₁ * f₂² * f₃³ * ... == f（至多差常数倍）
 * 2. 每个 fᵢ 无平方
 *
 * @return true 如果两个性质都满足
 */
static bool verify_sqfree_decomp(const MultiPoly& original,
                                 const SquareFreeDecomp& decomp,
                                 const std::string& main_var)
{
    const auto& comps = decomp.components;
    if (comps.empty()) return original.is_zero();

    // 性质 1：重构乘积 = f₁ * f₂² * f₃³ * ...
    MultiPoly product = poly_pow(comps[0], 1);
    for (size_t i = 1; i < comps.size(); ++i) {
        product = product * poly_pow(comps[i], static_cast<int>(i + 1));
    }

    // 比较乘积与原多项式（至多差常数倍）
    // 将两者本原化后比较
    if (original.is_zero() && product.is_zero()) return true;
    if (original.is_zero() || product.is_zero()) return false;

    MultiPoly orig_prim = original.make_primitive();
    MultiPoly prod_prim = product.make_primitive();

    // 检查是否相等或相差符号
    bool product_matches = (orig_prim == prod_prim) ||
                           (orig_prim == (prod_prim * Rational(-1)));

    if (!product_matches) return false;

    // 性质 2：每个非常数分量无平方
    for (size_t i = 0; i < comps.size(); ++i) {
        if (!comps[i].is_constant() && !is_square_free(comps[i], main_var)) {
            return false;
        }
    }

    return true;
}

int main()
{
    // ================================================================
    // Property 7: Square-free decomposition correctness
    // **Validates: Requirements 3.1, 3.2**
    // ================================================================

    TEST_CASE("Property 7: already square-free polynomial (x^2 - y^2)");
    {
        // x^2 - y^2 = (x+y)(x-y), already square-free
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1))   // -y^2
        };
        MultiPoly f(terms, vars);

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "x^2 - y^2: product matches and components are square-free");

        // Should have single component (already square-free)
        EXPECT_TRUE(decomp.components.size() == 1,
                    "x^2 - y^2: single component (already square-free)");
    }

    TEST_CASE("Property 7: polynomial with squared factor (x+y)^2*(x-y)");
    {
        // f = (x+y)^2 * (x-y) = (x^2 + 2xy + y^2)(x - y)
        //   = x^3 + 2x^2*y + x*y^2 - x^2*y - 2x*y^2 - y^3
        //   = x^3 + x^2*y - x*y^2 - y^3
        std::vector<std::string> vars = {"x", "y"};

        // Build (x+y)
        std::vector<MultiPoly::Term> xpy_terms = {
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 1}, Rational(1))   // y
        };
        MultiPoly xpy(xpy_terms, vars);

        // Build (x-y)
        std::vector<MultiPoly::Term> xmy_terms = {
            make_term({1, 0}, Rational(1)),   // x
            make_term({0, 1}, Rational(-1))   // -y
        };
        MultiPoly xmy(xmy_terms, vars);

        // f = (x+y)^2 * (x-y)
        MultiPoly f = xpy * xpy * xmy;

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "(x+y)^2*(x-y): product matches and components are square-free");
    }

    TEST_CASE("Property 7: polynomial with cubed factor (x+1)^3");
    {
        // f = (x+1)^3 = x^3 + 3x^2 + 3x + 1
        std::vector<std::string> vars = {"x", "y"};

        // Build (x+1) in {x, y} variable set
        std::vector<MultiPoly::Term> xp1_terms = {
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 0}, Rational(1))   // 1
        };
        MultiPoly xp1(xp1_terms, vars);

        // f = (x+1)^3
        MultiPoly f = xp1 * xp1 * xp1;

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "(x+1)^3: product matches and components are square-free");
    }

    TEST_CASE("Property 7: univariate with repeated roots (x-1)^2*(x+1)");
    {
        // f = (x-1)^2 * (x+1) = (x^2 - 2x + 1)(x + 1)
        //   = x^3 - x^2 - x + 1
        std::vector<std::string> vars = {"x"};

        // Build (x-1)
        std::vector<MultiPoly::Term> xm1_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly xm1(xm1_terms, vars);

        // Build (x+1)
        std::vector<MultiPoly::Term> xp1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // 1
        };
        MultiPoly xp1(xp1_terms, vars);

        // f = (x-1)^2 * (x+1)
        MultiPoly f = xm1 * xm1 * xp1;

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "(x-1)^2*(x+1): product matches and components are square-free");
    }

    TEST_CASE("Property 7: constant polynomial");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly f(Rational(42), vars);

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        // For constant, decomposition should return the constant itself
        EXPECT_TRUE(!decomp.components.empty(),
                    "constant 42: decomposition is non-empty");
        EXPECT_TRUE(decomp.components[0].is_constant(),
                    "constant 42: first component is constant");
    }

    TEST_CASE("Property 7: linear polynomial (x + 2y + 3)");
    {
        // Linear in x → already square-free
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 1}, Rational(2)),  // 2y
            make_term({0, 0}, Rational(3))   // 3
        };
        MultiPoly f(terms, vars);

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "x + 2y + 3: product matches and components are square-free");
        EXPECT_TRUE(decomp.components.size() == 1,
                    "x + 2y + 3: single component (linear, already square-free)");
    }

    TEST_CASE("Property 7: polynomial with squared and cubed factors (x+1)^2*(x-1)^3");
    {
        // f = (x+1)^2 * (x-1)^3
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> xp1_terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // 1
        };
        MultiPoly xp1(xp1_terms, vars);

        std::vector<MultiPoly::Term> xm1_terms = {
            make_term({1}, Rational(1)),   // x
            make_term({0}, Rational(-1))   // -1
        };
        MultiPoly xm1(xm1_terms, vars);

        // f = (x+1)^2 * (x-1)^3
        MultiPoly f = xp1 * xp1 * xm1 * xm1 * xm1;

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "(x+1)^2*(x-1)^3: product matches and components are square-free");
    }

    TEST_CASE("Property 7: multivariate with squared factor (x+y+1)^2*(x-y)");
    {
        // f = (x+y+1)^2 * (x-y)
        std::vector<std::string> vars = {"x", "y"};

        std::vector<MultiPoly::Term> xpy1_terms = {
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 1}, Rational(1)),  // y
            make_term({0, 0}, Rational(1))   // 1
        };
        MultiPoly xpy1(xpy1_terms, vars);

        std::vector<MultiPoly::Term> xmy_terms = {
            make_term({1, 0}, Rational(1)),   // x
            make_term({0, 1}, Rational(-1))   // -y
        };
        MultiPoly xmy(xmy_terms, vars);

        // f = (x+y+1)^2 * (x-y)
        MultiPoly f = xpy1 * xpy1 * xmy;

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(verify_sqfree_decomp(f, decomp, "x"),
                    "(x+y+1)^2*(x-y): product matches and components are square-free");
    }

    TEST_CASE("Property 7: zero polynomial");
    {
        MultiPoly f;  // zero polynomial

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        // Zero polynomial: decomposition should contain zero
        EXPECT_TRUE(!decomp.components.empty(),
                    "zero poly: decomposition is non-empty");
        EXPECT_TRUE(decomp.components[0].is_zero(),
                    "zero poly: first component is zero");
    }

    // ================================================================
    // Unit tests: verify specific component assignments
    // **Validates: Requirements 3.1, 3.2**
    // ================================================================

    TEST_CASE("unit: (x+y)^2*(x-y) component[0]=(x-y), component[1]=(x+y)");
    {
        std::vector<std::string> vars = {"x", "y"};

        std::vector<MultiPoly::Term> xpy_terms = {
            make_term({1, 0}, Rational(1)),
            make_term({0, 1}, Rational(1))
        };
        MultiPoly xpy(xpy_terms, vars);

        std::vector<MultiPoly::Term> xmy_terms = {
            make_term({1, 0}, Rational(1)),
            make_term({0, 1}, Rational(-1))
        };
        MultiPoly xmy(xmy_terms, vars);

        MultiPoly f = xpy * xpy * xmy;
        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(decomp.components.size() >= 2,
                    "(x+y)^2*(x-y): at least 2 components");

        // component[0] (multiplicity 1) == (x-y) up to constant
        MultiPoly comp0_prim = decomp.components[0].make_primitive();
        MultiPoly xmy_prim = xmy.make_primitive();
        EXPECT_TRUE(comp0_prim == xmy_prim || comp0_prim == (-xmy).make_primitive(),
                    "(x+y)^2*(x-y): component[0] == (x-y) up to sign");

        // component[1] (multiplicity 2) == (x+y) up to constant
        MultiPoly comp1_prim = decomp.components[1].make_primitive();
        MultiPoly xpy_prim = xpy.make_primitive();
        EXPECT_TRUE(comp1_prim == xpy_prim || comp1_prim == (-xpy).make_primitive(),
                    "(x+y)^2*(x-y): component[1] == (x+y) up to sign");

        // Verify reconstruction
        MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
        MultiPoly f_prim = f.make_primitive();
        MultiPoly recon_prim = reconstructed.make_primitive();
        EXPECT_TRUE(recon_prim == f_prim || recon_prim == (-f).make_primitive(),
                    "(x+y)^2*(x-y): reconstruction matches original");
    }

    TEST_CASE("unit: already square-free (x^2-y^2) returns single component");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),
            make_term({0, 2}, Rational(-1))
        };
        MultiPoly f(terms, vars);

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(decomp.components.size() == 1,
                    "x^2-y^2: exactly 1 component (already square-free)");

        MultiPoly comp_prim = decomp.components[0].make_primitive();
        MultiPoly f_prim = f.make_primitive();
        EXPECT_TRUE(comp_prim == f_prim || comp_prim == (-f).make_primitive(),
                    "x^2-y^2: component[0] == f up to sign");
    }

    TEST_CASE("unit: univariate (x-1)^2*(x+1) component[0]=(x+1), component[1]=(x-1)");
    {
        std::vector<std::string> vars = {"x"};

        std::vector<MultiPoly::Term> xm1_terms = {
            make_term({1}, Rational(1)),
            make_term({0}, Rational(-1))
        };
        MultiPoly xm1(xm1_terms, vars);

        std::vector<MultiPoly::Term> xp1_terms = {
            make_term({1}, Rational(1)),
            make_term({0}, Rational(1))
        };
        MultiPoly xp1(xp1_terms, vars);

        MultiPoly f = xm1 * xm1 * xp1;
        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(decomp.components.size() >= 2,
                    "(x-1)^2*(x+1): at least 2 components");

        // component[0] (multiplicity 1) == (x+1) up to constant
        MultiPoly comp0_prim = decomp.components[0].make_primitive();
        MultiPoly xp1_prim = xp1.make_primitive();
        EXPECT_TRUE(comp0_prim == xp1_prim || comp0_prim == (-xp1).make_primitive(),
                    "(x-1)^2*(x+1): component[0] == (x+1) up to sign");

        // component[1] (multiplicity 2) == (x-1) up to constant
        MultiPoly comp1_prim = decomp.components[1].make_primitive();
        MultiPoly xm1_prim = xm1.make_primitive();
        EXPECT_TRUE(comp1_prim == xm1_prim || comp1_prim == (-xm1).make_primitive(),
                    "(x-1)^2*(x+1): component[1] == (x-1) up to sign");

        // Verify reconstruction
        MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
        MultiPoly f_prim = f.make_primitive();
        MultiPoly recon_prim = reconstructed.make_primitive();
        EXPECT_TRUE(recon_prim == f_prim || recon_prim == (-f).make_primitive(),
                    "(x-1)^2*(x+1): reconstruction matches original");
    }

    TEST_CASE("unit: pure square (x+y)^2 has constant component[0], (x+y) as component[1]");
    {
        std::vector<std::string> vars = {"x", "y"};

        std::vector<MultiPoly::Term> xpy_terms = {
            make_term({1, 0}, Rational(1)),
            make_term({0, 1}, Rational(1))
        };
        MultiPoly xpy(xpy_terms, vars);

        MultiPoly f = xpy * xpy;
        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(decomp.components.size() >= 2,
                    "(x+y)^2: at least 2 components");

        // component[0] (multiplicity 1) should be constant (trivial)
        EXPECT_TRUE(decomp.components[0].is_constant(),
                    "(x+y)^2: component[0] is constant");

        // component[1] (multiplicity 2) == (x+y) up to constant
        MultiPoly comp1_prim = decomp.components[1].make_primitive();
        MultiPoly xpy_prim = xpy.make_primitive();
        EXPECT_TRUE(comp1_prim == xpy_prim || comp1_prim == (-xpy).make_primitive(),
                    "(x+y)^2: component[1] == (x+y) up to sign");

        // Verify reconstruction: 1 * (x+y)^2 == f
        MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
        MultiPoly f_prim = f.make_primitive();
        MultiPoly recon_prim = reconstructed.make_primitive();
        EXPECT_TRUE(recon_prim == f_prim || recon_prim == (-f).make_primitive(),
                    "(x+y)^2: reconstruction matches original");
    }

    TEST_CASE("unit: zero polynomial decomposition");
    {
        MultiPoly zero;
        SquareFreeDecomp decomp = square_free_decompose(zero, "x");

        EXPECT_TRUE(decomp.components.size() == 1,
                    "zero: exactly 1 component");
        EXPECT_TRUE(decomp.components[0].is_zero(),
                    "zero: component is zero polynomial");
    }

    TEST_CASE("unit: constant polynomial decomposition");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly f(Rational(42), vars);

        SquareFreeDecomp decomp = square_free_decompose(f, "x");

        EXPECT_TRUE(decomp.components.size() == 1,
                    "constant 42: exactly 1 component");
        EXPECT_TRUE(decomp.components[0].is_constant(),
                    "constant 42: component is constant");
    }

    return TEST_REPORT();
}
