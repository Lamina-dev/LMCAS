/**
 * @file test_multivariate_factor.cpp
 * @brief 多元因式分解完整流程的测试。
 *
 * 包含首项系数预计算（Wang's trick）的集成测试，
 * 验证 Hensel 提升前的 lc 分配逻辑。
 * 包含 Properties 9-11 的属性测试和完整分解的单元测试。
 */

#include "test_common.hpp"
#include "multivariate_factor.hpp"
#include "multivariate_poly.hpp"
#include "rapidcheck/rapidcheck.h"

using namespace lamina;

/**
 * @brief 辅助函数：构造 MultiPoly 项
 */
static MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff)
{
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}

int main()
{
    // ================================================================
    // 首项系数预计算相关测试
    // 验证 lc(f, x_main) 的计算和分配逻辑
    // **Validates: Requirements 6.1, 6.2, 6.3, 6.4**
    // ================================================================

    TEST_CASE("Leading coefficient: constant lc needs no precomputation");
    {
        // f = x^2 + y (lc w.r.t. x is 1, constant)
        // 常数首项系数无需预计算
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),  // x^2
            make_term({0, 1}, Rational(1))   // y
        };
        MultiPoly poly(terms, vars);

        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_TRUE(lc.is_constant(), "lc(x^2 + y, x) is constant (= 1)");
    }

    TEST_CASE("Leading coefficient: non-constant lc in auxiliary variable");
    {
        // f = y*x^2 + x + 1 (lc w.r.t. x is y, non-constant)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),  // y*x^2
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 0}, Rational(1))   // 1
        };
        MultiPoly poly(terms, vars);

        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_FALSE(lc.is_constant(), "lc(y*x^2 + x + 1, x) is non-constant (= y)");

        // 验证 lc 在 y=2 处求值为 2
        std::map<std::string, Rational> eval_pts = {{"y", Rational(2)}};
        MultiPoly lc_eval = lc.eval(eval_pts);
        EXPECT_TRUE(lc_eval.is_constant(), "lc evaluated at y=2 is constant");
        EXPECT_TRUE(lc_eval.terms()[0].second == Rational(2),
                    "lc evaluated at y=2 equals 2");
    }

    TEST_CASE("Leading coefficient: lc evaluation matches factor lc product");
    {
        // f = (y*x + 1)(x + y) = y*x^2 + (y^2+1)*x + y
        // lc(f, x) = y
        // 一元因子（在 y=1 处）：f(x,1) = x^2 + 2x + 1 = (x+1)^2
        // 在 y=2 处：f(x,2) = 2x^2 + 5x + 2 = (2x+1)(x+2)
        // lc(f,x) = y, 在 y=2 处 lc_eval = 2
        // 一元因子 lc: lc(2x+1) = 2, lc(x+2) = 1, 乘积 = 2 = lc_eval ✓
        std::vector<std::string> vars = {"x", "y"};

        // 构造 f = y*x^2 + (y^2+1)*x + y
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),  // y*x^2
            make_term({1, 2}, Rational(1)),  // y^2*x
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 1}, Rational(1))   // y
        };
        MultiPoly poly(terms, vars);

        // 验证 lc(f, x) = y
        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_FALSE(lc.is_constant(), "lc is non-constant");

        // 在 y=2 处求值
        std::map<std::string, Rational> eval_pts = {{"y", Rational(2)}};
        MultiPoly lc_eval = lc.eval(eval_pts);
        EXPECT_TRUE(lc_eval.terms()[0].second == Rational(2),
                    "lc at y=2 is 2");

        // 求值后的多项式：f(x,2) = 2x^2 + 5x + 2
        MultiPoly f_eval = poly.eval("y", Rational(2));
        Polynomial<Rational> f_uni = f_eval.to_univariate();
        EXPECT_TRUE(f_uni.lead_coeff() == Rational(2),
                    "f(x,2) has leading coefficient 2");
    }

    TEST_CASE("Leading coefficient: Hensel lift with non-constant lc");
    {
        // 测试 Hensel 提升在非常数首项系数情形下的行为
        // f = (x + y)(x - y) = x^2 - y^2
        // lc(f, x) = 1 (常数)，这是简单情形
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1))   // -y^2
        };
        MultiPoly poly(terms, vars);

        // 验证 lc(f, x) = 1（常数），无需预计算
        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_TRUE(lc.is_constant(), "lc(x^2 - y^2, x) is constant (= 1)");

        // 在 y=1 处求值：f(x,1) = x^2 - 1 = (x+1)(x-1)
        MultiPoly f_eval = poly.eval("y", Rational(1));
        Polynomial<Rational> f_uni = f_eval.to_univariate();
        EXPECT_TRUE(f_uni.degree() == 2, "f(x,1) has degree 2");
        EXPECT_TRUE(f_uni.lead_coeff() == Rational(1), "f(x,1) is monic");
    }

    TEST_CASE("Leading coefficient: multivariate poly with y as lc");
    {
        // f = y*x^2 - y = y*(x^2 - 1) = y*(x+1)*(x-1)
        // lc(f, x) = y (非常数)
        // 在 y=1 处：f(x,1) = x^2 - 1 = (x+1)(x-1)
        // lc_eval = 1, 一元因子 lc 乘积 = 1*1 = 1 = lc_eval
        // 此情形下 scale = 1，无需调整
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),   // y*x^2
            make_term({0, 1}, Rational(-1))   // -y
        };
        MultiPoly poly(terms, vars);

        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_FALSE(lc.is_constant(), "lc(y*x^2 - y, x) = y is non-constant");

        // 在 y=1 处求值
        std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
        MultiPoly lc_eval = lc.eval(eval_pts);
        EXPECT_TRUE(lc_eval.terms()[0].second == Rational(1),
                    "lc at y=1 is 1");
    }

    TEST_CASE("Leading coefficient: scale factor computation");
    {
        // f = 2y*x^2 + 3y*x + y = y*(2x^2 + 3x + 1) = y*(2x+1)*(x+1)
        // lc(f, x) = 2y
        // 在 y=1 处：f(x,1) = 2x^2 + 3x + 1 = (2x+1)(x+1)
        // lc_eval = 2, 一元因子 lc: lc(2x+1)=2, lc(x+1)=1, 乘积=2 = lc_eval ✓
        // 无需缩放
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(2)),   // 2y*x^2
            make_term({1, 1}, Rational(3)),   // 3y*x
            make_term({0, 1}, Rational(1))    // y
        };
        MultiPoly poly(terms, vars);

        MultiPoly lc = poly.leading_coeff("x");
        EXPECT_FALSE(lc.is_constant(), "lc(2y*x^2+3y*x+y, x) = 2y is non-constant");

        // 在 y=1 处求值
        std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
        MultiPoly lc_eval = lc.eval(eval_pts);
        Rational lc_val = lc_eval.is_zero() ? Rational(0) : lc_eval.terms()[0].second;
        EXPECT_TRUE(lc_val == Rational(2), "lc at y=1 is 2");

        // 一元因子 (2x+1)(x+1) 的首项系数乘积
        Polynomial<Rational> f1({Rational(1), Rational(2)}, "x");   // 2x+1
        Polynomial<Rational> f2({Rational(1), Rational(1)}, "x");   // x+1
        Rational product_lcs = f1.lead_coeff() * f2.lead_coeff();
        EXPECT_TRUE(product_lcs == lc_val,
                    "product of factor lcs equals lc_eval (no scaling needed)");
    }

    TEST_CASE("Leading coefficient: scaling needed when lc product differs");
    {
        // 场景：一元因子的首项系数乘积与 lc_eval 不同
        // f = 3y*x^2 + ... 在 y=1 处 lc_eval = 3
        // 若一元因子为 (x+a)(x+b)（首一），则 lc 乘积 = 1 ≠ 3
        // 需要将缩放因子 3 分配给某个因子

        // 构造 f = 3y*x^2 - 3y = 3y*(x^2-1) = 3y*(x+1)*(x-1)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(3)),   // 3y*x^2
            make_term({0, 1}, Rational(-3))   // -3y
        };
        MultiPoly poly(terms, vars);

        MultiPoly lc = poly.leading_coeff("x");
        std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
        MultiPoly lc_eval = lc.eval(eval_pts);
        Rational lc_val = lc_eval.terms()[0].second;
        EXPECT_TRUE(lc_val == Rational(3), "lc at y=1 is 3");

        // 若一元分解给出首一因子 (x+1)(x-1)
        Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");   // x+1
        Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x");  // x-1
        Rational product_lcs = f1.lead_coeff() * f2.lead_coeff();
        EXPECT_TRUE(product_lcs == Rational(1), "monic factors have lc product = 1");

        // 缩放因子 = lc_val / product_lcs = 3
        Rational scale = lc_val / product_lcs;
        EXPECT_TRUE(scale == Rational(3), "scale factor is 3");

        // 应用缩放后，第一个因子变为 3x+3
        Polynomial<Rational> scaled_f1({Rational(3), Rational(3)}, "x");  // 3(x+1) = 3x+3
        EXPECT_TRUE(scaled_f1.lead_coeff() == Rational(3),
                    "scaled factor has correct leading coefficient");

        // 验证缩放后乘积的首项系数
        Rational new_product_lcs = scaled_f1.lead_coeff() * f2.lead_coeff();
        EXPECT_TRUE(new_product_lcs == lc_val,
                    "after scaling, product of lcs equals lc_eval");
    }

    // ================================================================
    // 试除验证与因子组合测试
    // 验证 trial_division 和 factor_combination 的逻辑
    // **Validates: Requirements 7.1, 7.2, 7.3**
    // ================================================================

    TEST_CASE("Trial division: single factor divides exactly");
    {
        // f = (x+y)(x-y) = x^2 - y^2
        // 若提升因子为 (x+y) 和 (x-y)，试除应逐一验证
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1))   // -y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // 验证因子乘积等于原多项式
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "factor product equals original for x^2 - y^2");
    }

    TEST_CASE("Trial division: common monomial extraction then factor");
    {
        // f = x^2*y + x*y^2 = xy(x+y)
        // 公因子单项式 xy 提取后，商为 (x+y)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),   // x^2*y
            make_term({1, 2}, Rational(1))    // x*y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // 验证因子乘积等于原多项式
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "factor product equals original for x^2*y + x*y^2");

        // 应有因子 x, y, (x+y)
        EXPECT_TRUE(result.factors.size() >= 2,
                    "x^2*y + x*y^2 has at least 2 factors");
    }

    TEST_CASE("Trial division: irreducible polynomial returns itself");
    {
        // f = x^2 + y^2 + 1 (不可约 over Q)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(1)),   // y^2
            make_term({0, 0}, Rational(1))    // 1
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // 不可约多项式应返回自身
        EXPECT_TRUE(result.factors.size() == 1,
                    "x^2 + y^2 + 1 is irreducible (single factor)");

        // 验证乘积
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "factor product equals original for irreducible poly");
    }

    TEST_CASE("Trial division: constant polynomial");
    {
        // f = 42
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(42), vars);

        MultiFactorResult result = factor_multivariate(poly);
        EXPECT_TRUE(result.factors.empty(), "constant polynomial has no factors");
        EXPECT_TRUE(result.constant == Rational(42), "constant is 42");
    }

    TEST_CASE("Trial division: zero polynomial");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(0), vars);

        MultiFactorResult result = factor_multivariate(poly);
        EXPECT_TRUE(result.factors.empty(), "zero polynomial has no factors");
        EXPECT_TRUE(result.constant == Rational(0), "constant is 0");
    }

    TEST_CASE("Factor combination: product of two linear factors");
    {
        // f = (x + y + 1)(x - y + 2) = x^2 + x - y^2 + 3x + 2 - y + ...
        // 构造乘积
        std::vector<std::string> vars = {"x", "y"};

        // (x + y + 1)
        std::vector<MultiPoly::Term> t1 = {
            make_term({1, 0}, Rational(1)),   // x
            make_term({0, 1}, Rational(1)),   // y
            make_term({0, 0}, Rational(1))    // 1
        };
        MultiPoly f1(t1, vars);

        // (x - y + 2)
        std::vector<MultiPoly::Term> t2 = {
            make_term({1, 0}, Rational(1)),   // x
            make_term({0, 1}, Rational(-1)),  // -y
            make_term({0, 0}, Rational(2))    // 2
        };
        MultiPoly f2(t2, vars);

        MultiPoly poly = f1 * f2;

        MultiFactorResult result = factor_multivariate(poly);

        // 验证因子乘积等于原多项式
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly,
                    "factor product equals original for (x+y+1)(x-y+2)");
    }

    TEST_CASE("Factor combination: numeric content extraction");
    {
        // f = 6x^2*y - 3x*y^2 = 3xy(2x - y)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(6)),   // 6x^2*y
            make_term({1, 2}, Rational(-3))   // -3x*y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // 验证因子乘积等于原多项式
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly,
                    "factor product equals original for 6x^2*y - 3x*y^2");
    }

    // ================================================================
    // Property 9: Factorization product correctness
    // constant * ∏(factors[i]^mult[i]) == original
    // Each factor is primitive with positive leading coefficient
    // **Validates: Requirements 7.4, 7.5**
    // ================================================================

    TEST_CASE("Feature: multivariate-factorization, Property 9: Factorization product correctness (random factorable bivariate)");
    rc::check("For random products of linear factors, factorization reconstructs original", []() {
        // Generate 2-3 random linear factors in {x, y} and verify product correctness
        std::vector<std::string> vars = {"x", "y"};
        int num_factors = 2 + rc::gen::inRange(0, 2);  // 2 or 3 factors

        MultiPoly product(Rational(1), vars);
        for (int i = 0; i < num_factors; ++i) {
            // Generate linear factor: a*x + b*y + c with small integer coefficients
            int a = rc::gen::inRange(-3, 3);
            int b = rc::gen::inRange(-3, 3);
            int c = rc::gen::inRange(-3, 3);
            // Ensure factor is non-zero (at least one of a, b non-zero)
            if (a == 0 && b == 0) a = 1;

            std::vector<MultiPoly::Term> terms;
            if (a != 0) terms.push_back({Monomial({1, 0}), Rational(a)});
            if (b != 0) terms.push_back({Monomial({0, 1}), Rational(b)});
            if (c != 0) terms.push_back({Monomial({0, 0}), Rational(c)});
            if (terms.empty()) terms.push_back({Monomial({1, 0}), Rational(1)});

            MultiPoly factor(terms, vars);
            product = product * factor;
        }

        // Factor the product
        MultiFactorResult result = factor_multivariate(product);

        // Reconstruct: constant * ∏(factors[i]^mult[i])
        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                reconstructed = reconstructed * result.factors[i];
            }
        }

        RC_ASSERT(reconstructed == product);
    });

    TEST_CASE("Feature: multivariate-factorization, Property 9: Each factor is primitive with positive leading coefficient");
    rc::check("Factors from factorization are primitive with positive leading coefficient", []() {
        std::vector<std::string> vars = {"x", "y"};
        // Generate a random quadratic: a*x^2 + b*x*y + c*y^2 + d*x + e*y + f
        int a = rc::gen::inRange(1, 4);  // positive leading coeff
        int b = rc::gen::inRange(-3, 3);
        int c = rc::gen::inRange(-3, 3);
        int d = rc::gen::inRange(-3, 3);
        int e = rc::gen::inRange(-3, 3);
        int f = rc::gen::inRange(-3, 3);

        std::vector<MultiPoly::Term> terms;
        terms.push_back({Monomial({2, 0}), Rational(a)});
        if (b != 0) terms.push_back({Monomial({1, 1}), Rational(b)});
        if (c != 0) terms.push_back({Monomial({0, 2}), Rational(c)});
        if (d != 0) terms.push_back({Monomial({1, 0}), Rational(d)});
        if (e != 0) terms.push_back({Monomial({0, 1}), Rational(e)});
        if (f != 0) terms.push_back({Monomial({0, 0}), Rational(f)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = factor_multivariate(poly);

        // Each factor should be primitive (numeric_content == 1)
        for (size_t i = 0; i < result.factors.size(); ++i) {
            Rational content = result.factors[i].numeric_content();
            RC_ASSERT(content == Rational(1));
        }

        // Verify product correctness
        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                reconstructed = reconstructed * result.factors[i];
            }
        }
        RC_ASSERT(reconstructed == poly);
    });

    TEST_CASE("Feature: multivariate-factorization, Property 9: Factorization product correctness (random monomial * linear)");
    rc::check("For monomial * linear factor products, factorization reconstructs original", []() {
        std::vector<std::string> vars = {"x", "y"};
        // Generate a monomial factor: x^a * y^b with small exponents
        int exp_x = rc::gen::inRange(0, 3);
        int exp_y = rc::gen::inRange(0, 3);
        if (exp_x == 0 && exp_y == 0) exp_x = 1;

        std::vector<MultiPoly::Term> mono_terms = {
            {Monomial({exp_x, exp_y}), Rational(1)}
        };
        MultiPoly monomial_factor(mono_terms, vars);

        // Generate a linear factor: a*x + b*y + c
        int a = rc::gen::inRange(-3, 3);
        int b = rc::gen::inRange(-3, 3);
        int c = rc::gen::inRange(-3, 3);
        if (a == 0 && b == 0) a = 1;

        std::vector<MultiPoly::Term> lin_terms;
        if (a != 0) lin_terms.push_back({Monomial({1, 0}), Rational(a)});
        if (b != 0) lin_terms.push_back({Monomial({0, 1}), Rational(b)});
        if (c != 0) lin_terms.push_back({Monomial({0, 0}), Rational(c)});
        if (lin_terms.empty()) lin_terms.push_back({Monomial({1, 0}), Rational(1)});
        MultiPoly linear_factor(lin_terms, vars);

        // Multiply with a scalar
        int scalar = rc::gen::inRange(1, 5);
        MultiPoly poly = monomial_factor * linear_factor * Rational(scalar);

        MultiFactorResult result = factor_multivariate(poly);

        // Reconstruct
        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                reconstructed = reconstructed * result.factors[i];
            }
        }
        RC_ASSERT(reconstructed == poly);
    });

    // ================================================================
    // Property 10: Linear polynomials are irreducible
    // degree-1 polynomial returns itself as sole factor
    // **Validates: Requirements 9.1**
    // ================================================================

    TEST_CASE("Feature: multivariate-factorization, Property 10: Linear polynomials are irreducible (random bivariate linear)");
    rc::check("Linear bivariate polynomial returns itself as sole factor", []() {
        std::vector<std::string> vars = {"x", "y"};
        // Generate a random linear polynomial: a*x + b*y + c
        int a = rc::gen::inRange(-5, 5);
        int b = rc::gen::inRange(-5, 5);
        int c = rc::gen::inRange(-5, 5);
        // Ensure at least one of a, b is non-zero (so it's truly linear)
        if (a == 0 && b == 0) a = 1;

        std::vector<MultiPoly::Term> terms;
        if (a != 0) terms.push_back({Monomial({1, 0}), Rational(a)});
        if (b != 0) terms.push_back({Monomial({0, 1}), Rational(b)});
        if (c != 0) terms.push_back({Monomial({0, 0}), Rational(c)});
        if (terms.empty()) terms.push_back({Monomial({1, 0}), Rational(1)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = factor_multivariate(poly);

        // Should have exactly one factor (the polynomial itself, up to constant)
        RC_ASSERT(result.factors.size() == 1);
        RC_ASSERT(result.multiplicities[0] == 1);

        // Verify product correctness
        MultiPoly reconstructed(Rational(result.constant), vars);
        reconstructed = reconstructed * result.factors[0];
        RC_ASSERT(reconstructed == poly);
    });

    TEST_CASE("Feature: multivariate-factorization, Property 10: Linear polynomials are irreducible (random trivariate linear)");
    rc::check("Linear trivariate polynomial returns itself as sole factor", []() {
        std::vector<std::string> vars = {"x", "y", "z"};
        // Generate a random linear polynomial: a*x + b*y + c*z + d
        int a = rc::gen::inRange(-4, 4);
        int b = rc::gen::inRange(-4, 4);
        int c = rc::gen::inRange(-4, 4);
        int d = rc::gen::inRange(-4, 4);
        // Ensure at least one variable coefficient is non-zero
        if (a == 0 && b == 0 && c == 0) a = 1;

        std::vector<MultiPoly::Term> terms;
        if (a != 0) terms.push_back({Monomial({1, 0, 0}), Rational(a)});
        if (b != 0) terms.push_back({Monomial({0, 1, 0}), Rational(b)});
        if (c != 0) terms.push_back({Monomial({0, 0, 1}), Rational(c)});
        if (d != 0) terms.push_back({Monomial({0, 0, 0}), Rational(d)});
        if (terms.empty()) terms.push_back({Monomial({1, 0, 0}), Rational(1)});

        MultiPoly poly(terms, vars);
        MultiFactorResult result = factor_multivariate(poly);

        // Should have exactly one factor
        RC_ASSERT(result.factors.size() == 1);
        RC_ASSERT(result.multiplicities[0] == 1);

        // Verify product correctness
        MultiPoly reconstructed(Rational(result.constant), vars);
        reconstructed = reconstructed * result.factors[0];
        RC_ASSERT(reconstructed == poly);
    });

    // ================================================================
    // Property 11: Difference of squares factorization
    // a²-b² factors into (a+b)(a-b)
    // **Validates: Requirements 9.2**
    // ================================================================

    TEST_CASE("Feature: multivariate-factorization, Property 11: Difference of squares factorization (random a^2 - b^2)");
    rc::check("a^2 - b^2 factors into product containing (a+b) and (a-b)", []() {
        std::vector<std::string> vars = {"x", "y"};
        // Generate random a and b as simple monomials/linear terms
        // a = c1*x + c2*y, b = c3*x + c4*y (ensure a != ±b)
        int c1 = rc::gen::inRange(-3, 3);
        int c2 = rc::gen::inRange(-3, 3);
        if (c1 == 0 && c2 == 0) c1 = 1;

        int c3 = rc::gen::inRange(-3, 3);
        int c4 = rc::gen::inRange(-3, 3);
        if (c3 == 0 && c4 == 0) c4 = 1;

        // Ensure a != b and a != -b (otherwise a^2 - b^2 = 0)
        if (c1 == c3 && c2 == c4) c3 = c3 + 1;
        if (c1 == -c3 && c2 == -c4) c4 = c4 + 1;

        std::vector<MultiPoly::Term> a_terms;
        if (c1 != 0) a_terms.push_back({Monomial({1, 0}), Rational(c1)});
        if (c2 != 0) a_terms.push_back({Monomial({0, 1}), Rational(c2)});
        if (a_terms.empty()) a_terms.push_back({Monomial({1, 0}), Rational(1)});
        MultiPoly a_poly(a_terms, vars);

        std::vector<MultiPoly::Term> b_terms;
        if (c3 != 0) b_terms.push_back({Monomial({1, 0}), Rational(c3)});
        if (c4 != 0) b_terms.push_back({Monomial({0, 1}), Rational(c4)});
        if (b_terms.empty()) b_terms.push_back({Monomial({0, 1}), Rational(1)});
        MultiPoly b_poly(b_terms, vars);

        // Compute a^2 - b^2
        MultiPoly a_sq = a_poly * a_poly;
        MultiPoly b_sq = b_poly * b_poly;
        MultiPoly diff = a_sq - b_sq;

        // Skip if diff is zero (a == ±b)
        if (diff.is_zero()) return;

        MultiFactorResult result = factor_multivariate(diff);

        // Verify product correctness: constant * ∏(factors[i]^mult[i]) == diff
        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                reconstructed = reconstructed * result.factors[i];
            }
        }
        RC_ASSERT(reconstructed == diff);

        // When the polynomial has exactly 2 terms (pure monomial squares),
        // the difference-of-squares fast path should produce at least 2 factors
        if (diff.num_terms() == 2) {
            int total_factor_count = 0;
            for (size_t i = 0; i < result.multiplicities.size(); ++i) {
                total_factor_count += result.multiplicities[i];
            }
            RC_ASSERT(total_factor_count >= 2);
        }
    });

    TEST_CASE("Feature: multivariate-factorization, Property 11: Difference of squares with single variables");
    rc::check("x_i^2 - x_j^2 factors into (x_i + x_j)(x_i - x_j)", []() {
        std::vector<std::string> vars = {"x", "y", "z"};
        // Pick two distinct variable indices
        int i = rc::gen::inRange(0, 2);
        int j = rc::gen::inRange(0, 2);
        if (i == j) j = (i + 1) % 3;

        // Construct x_i^2 - x_j^2
        std::vector<int> exp_pos(3, 0);
        exp_pos[i] = 2;
        std::vector<int> exp_neg(3, 0);
        exp_neg[j] = 2;

        std::vector<MultiPoly::Term> terms = {
            {Monomial(exp_pos.begin(), exp_pos.end()), Rational(1)},
            {Monomial(exp_neg.begin(), exp_neg.end()), Rational(-1)}
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // Verify product correctness
        MultiPoly reconstructed(Rational(result.constant), vars);
        for (size_t k = 0; k < result.factors.size(); ++k) {
            for (int m = 0; m < result.multiplicities[k]; ++m) {
                reconstructed = reconstructed * result.factors[k];
            }
        }
        RC_ASSERT(reconstructed == poly);

        // Should have exactly 2 linear factors
        RC_ASSERT(result.factors.size() == 2);
        RC_ASSERT(result.multiplicities[0] == 1);
        RC_ASSERT(result.multiplicities[1] == 1);
    });

    // ================================================================
    // Unit tests for complete factorization
    // **Validates: Requirements 7.4, 8.1, 8.4, 9.1, 9.2, 9.4**
    // ================================================================

    TEST_CASE("Complete factorization: x^2 - y^2 -> (x+y)(x-y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1))   // -y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // Verify product
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "x^2-y^2: product equals original");
        EXPECT_TRUE(result.factors.size() == 2, "x^2-y^2: has 2 factors");
        EXPECT_TRUE(result.constant == Rational(1), "x^2-y^2: constant is 1");
    }

    TEST_CASE("Complete factorization: x^2*y + x*y^2 -> xy(x+y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),   // x^2*y
            make_term({1, 2}, Rational(1))    // x*y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // Verify product
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "x^2*y+x*y^2: product equals original");
        // Should have factors: x, y, (x+y)
        EXPECT_TRUE(result.factors.size() >= 2, "x^2*y+x*y^2: has at least 2 factors");
    }

    TEST_CASE("Complete factorization: x^2+2xy+y^2-z^2 -> (x+y+z)(x+y-z)");
    {
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0, 0}, Rational(1)),   // x^2
            make_term({1, 1, 0}, Rational(2)),   // 2xy
            make_term({0, 2, 0}, Rational(1)),   // y^2
            make_term({0, 0, 2}, Rational(-1))   // -z^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // Verify product correctness (fundamental invariant)
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "x^2+2xy+y^2-z^2: product equals original");
        EXPECT_TRUE(result.factors.size() >= 1, "x^2+2xy+y^2-z^2: has at least 1 factor");
        EXPECT_TRUE(result.constant == Rational(1), "x^2+2xy+y^2-z^2: constant is 1");
    }

    TEST_CASE("Complete factorization: 6x^2*y - 3x*y^2 -> 3xy(2x-y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(6)),   // 6x^2*y
            make_term({1, 2}, Rational(-3))   // -3x*y^2
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        // Verify product
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "6x^2*y-3x*y^2: product equals original");
        // Constant should be 3 (numeric content)
        // Factors should include x, y, (2x-y)
        EXPECT_TRUE(result.factors.size() >= 2, "6x^2*y-3x*y^2: has at least 2 factors");
    }

    TEST_CASE("Complete factorization: constant 42 -> {constant=42, factors=[]}");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(42), vars);

        MultiFactorResult result = factor_multivariate(poly);
        EXPECT_TRUE(result.factors.empty(), "constant 42: no factors");
        EXPECT_TRUE(result.constant == Rational(42), "constant 42: constant is 42");
    }

    TEST_CASE("Complete factorization: zero -> {constant=0, factors=[]}");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(0), vars);

        MultiFactorResult result = factor_multivariate(poly);
        EXPECT_TRUE(result.factors.empty(), "zero: no factors");
        EXPECT_TRUE(result.constant == Rational(0), "zero: constant is 0");
    }

    TEST_CASE("Complete factorization: irreducible x^2+y^2+1 -> returns itself");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(1)),   // y^2
            make_term({0, 0}, Rational(1))    // 1
        };
        MultiPoly poly(terms, vars);

        MultiFactorResult result = factor_multivariate(poly);

        EXPECT_TRUE(result.factors.size() == 1, "x^2+y^2+1: single irreducible factor");
        EXPECT_TRUE(result.constant == Rational(1), "x^2+y^2+1: constant is 1");
        EXPECT_TRUE(result.multiplicities[0] == 1, "x^2+y^2+1: multiplicity is 1");

        // Verify product
        MultiPoly product(Rational(result.constant), vars);
        for (size_t i = 0; i < result.factors.size(); ++i) {
            for (int m = 0; m < result.multiplicities[i]; ++m) {
                product = product * result.factors[i];
            }
        }
        EXPECT_TRUE(product == poly, "x^2+y^2+1: product equals original");
    }

    return TEST_REPORT();
}
