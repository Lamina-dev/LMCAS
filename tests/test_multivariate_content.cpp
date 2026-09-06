/**
 * @file test_multivariate_content.cpp
 * @brief 多元容度/本原部分和多元 GCD 的测试。
 */

#include "test_common.hpp"
#include "multivariate_factor.hpp"

using namespace LMCAS;

/**
 * @brief 辅助函数：构造 MultiPoly 项
 */
static MultiPoly::Term make_term(const std::vector<int>& exponents, const Rational& coeff)
{
    return {Monomial(exponents.begin(), exponents.end()), coeff};
}

/**
 * @brief 将 content（辅助变量多项式）嵌入到完整变量集中
 *
 * multivariate_content 返回的多项式不含 main_var，需要将其嵌入到
 * 包含 main_var 的完整变量集中才能与原多项式做乘法比较。
 */
static MultiPoly embed_content(const MultiPoly& content,
                               const std::string& main_var,
                               const std::vector<std::string>& full_vars)
{
    // 找到 main_var 在 full_vars 中的位置
    int var_idx = -1;
    for (size_t i = 0; i < full_vars.size(); ++i) {
        if (full_vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) return content;  // main_var 不在 full_vars 中

    // 将 content 的每个项的单项式扩展到 full_vars 维度
    std::vector<MultiPoly::Term> new_terms;
    for (const auto& term : content.terms()) {
        Monomial full_mono(full_vars.size(), 0);
        size_t ri = 0;
        for (size_t i = 0; i < full_vars.size(); ++i) {
            if (static_cast<int>(i) == var_idx) {
                full_mono[i] = 0;  // main_var 指数为 0
            } else {
                if (ri < term.first.size()) full_mono[i] = term.first[ri];
                ++ri;
            }
        }
        new_terms.emplace_back(std::move(full_mono), term.second);
    }
    return MultiPoly(std::move(new_terms), full_vars);
}

int main()
{

    TEST_CASE("multivariate_content: zero polynomial returns zero");
    {
        MultiPoly zero;
        MultiPoly content = multivariate_content(zero, "x");
        EXPECT_TRUE(content.is_zero(), "content of zero poly is zero");
    }

    TEST_CASE("multivariate_content: constant polynomial");
    {
        // poly = 6, main_var = "x"
        // 视为 x 的 0 次多项式，系数为 6（常数）
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly poly(Rational(6), vars);
        MultiPoly content = multivariate_content(poly, "x");
        EXPECT_TRUE(content.is_constant(), "content of constant is constant");
    }

    TEST_CASE("multivariate_content: poly does not contain main_var");
    {
        // poly = y^2 + y, main_var = "x"
        // poly 不含 x，容度 = poly 本身
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({0, 2}, Rational(1)),  // y^2
            make_term({0, 1}, Rational(1))   // y
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");
        EXPECT_TRUE(content == poly, "content when main_var absent equals poly itself");
    }

    TEST_CASE("multivariate_content: x^2 + x (content = 1, all coeffs constant)");
    {
        // poly = x^2 + x, main_var = "x"
        // coeff(x^2) = 1, coeff(x^1) = 1 → content = 1
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),  // x^2
            make_term({1, 0}, Rational(1))   // x
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");
        EXPECT_TRUE(content.is_constant(), "content of x^2 + x w.r.t. x is constant");
    }

    TEST_CASE("multivariate_content: univariate poly x^3 + x^2 + x (content = 1)");
    {
        // All coefficients are 1 (constants), content = gcd(1, 1, 1) = 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            make_term({3}, Rational(1)),
            make_term({2}, Rational(1)),
            make_term({1}, Rational(1))
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");
        EXPECT_TRUE(content.is_constant(), "content of x^3+x^2+x is constant 1");
    }

    TEST_CASE("multivariate_content: single term x^3*y^2 (content = y^2)");
    {
        // poly = x^3*y^2, main_var = "x"
        // 只有一个系数多项式：coeff(x^3) = y^2 → content = y^2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({3, 2}, Rational(1))  // x^3 * y^2
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");

        std::vector<std::string> remaining_vars = {"y"};
        std::vector<MultiPoly::Term> expected_terms = {
            make_term({2}, Rational(1))  // y^2
        };
        MultiPoly expected(expected_terms, remaining_vars);

        EXPECT_EQ_STR(content.to_string(), expected.to_string(),
                      "content of x^3*y^2 w.r.t. x is y^2");
    }

    TEST_CASE("multivariate_content: main_var not in variable list");
    {
        // poly = x + 1, main_var = "z" (not present)
        // 容度 = poly 本身
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // 1
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "z");
        EXPECT_TRUE(content == poly, "content with absent main_var returns poly");
    }

    TEST_CASE("multivariate_content: multiple terms same x-degree grouped correctly");
    {
        // poly = x^2*y + x^2*z, main_var = "x"
        // coeff(x^2) = y + z → single coefficient, content = y + z
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1, 0}, Rational(1)),  // x^2 * y
            make_term({2, 0, 1}, Rational(1))   // x^2 * z
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");

        // 只有一个系数多项式 (y + z)，容度 = y + z
        std::vector<std::string> remaining_vars = {"y", "z"};
        std::vector<MultiPoly::Term> expected_terms = {
            make_term({1, 0}, Rational(1)),  // y
            make_term({0, 1}, Rational(1))   // z
        };
        MultiPoly expected(expected_terms, remaining_vars);

        EXPECT_EQ_STR(content.to_string(), expected.to_string(),
                      "content of x^2*y + x^2*z w.r.t. x is y + z");
    }


    TEST_CASE("multivariate_content: x^2*y + x*y (content = y) [requires GCD]");
    {
        // poly = x^2*y + x*y, main_var = "x"
        // coeff(x^2) = y, coeff(x^1) = y → content = gcd(y, y) = y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),  // x^2 * y
            make_term({1, 1}, Rational(1))   // x * y
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");

        std::vector<std::string> remaining_vars = {"y"};
        std::vector<MultiPoly::Term> expected_terms = {
            make_term({1}, Rational(1))  // y
        };
        MultiPoly expected(expected_terms, remaining_vars);

        // 注意：此测试在 multivariate_gcd 实现前会失败
        EXPECT_EQ_STR(content.to_string(), expected.to_string(),
                      "content of x^2*y + x*y w.r.t. x is y");
    }

    TEST_CASE("multivariate_content: 2*x^2*y + 4*x*y^2 (content = 2y) [requires GCD]");
    {
        // coeff(x^2) = 2y, coeff(x^1) = 4y^2
        // content = gcd(2y, 4y^2) = 2y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(2)),  // 2*x^2*y
            make_term({1, 2}, Rational(4))   // 4*x*y^2
        };
        MultiPoly poly(terms, vars);
        MultiPoly content = multivariate_content(poly, "x");

        std::vector<std::string> remaining_vars = {"y"};
        std::vector<MultiPoly::Term> expected_terms = {
            make_term({1}, Rational(2))  // 2*y
        };
        MultiPoly expected(expected_terms, remaining_vars);

        // 注意：此测试在 multivariate_gcd 实现前会失败
        EXPECT_EQ_STR(content.to_string(), expected.to_string(),
                      "content of 2x^2*y + 4x*y^2 w.r.t. x is 2y");
    }


    TEST_CASE("multivariate_gcd: gcd(x^2*y + x*y^2, x*y) == x*y");
    {
        // gcd(x^2*y + x*y^2, x*y) should be x*y (up to scalar)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms_a = {
            make_term({2, 1}, Rational(1)),  // x^2*y
            make_term({1, 2}, Rational(1))   // x*y^2
        };
        MultiPoly a(terms_a, vars);

        std::vector<MultiPoly::Term> terms_b = {
            make_term({1, 1}, Rational(1))   // x*y
        };
        MultiPoly b(terms_b, vars);

        MultiPoly g = multivariate_gcd(a, b);
        MultiPoly g_prim = g.make_primitive();

        // Expected: x*y (primitive)
        MultiPoly expected(terms_b, vars);
        MultiPoly expected_prim = expected.make_primitive();

        EXPECT_EQ_STR(g_prim.to_string(), expected_prim.to_string(),
                      "gcd(x^2*y + x*y^2, x*y) == x*y");
    }

    TEST_CASE("multivariate_gcd: coprime polynomials return 1");
    {
        // gcd(x+1, x+2) should be 1 (coprime)
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms_a = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(1))   // 1
        };
        MultiPoly a(terms_a, vars);

        std::vector<MultiPoly::Term> terms_b = {
            make_term({1}, Rational(1)),  // x
            make_term({0}, Rational(2))   // 2
        };
        MultiPoly b(terms_b, vars);

        MultiPoly g = multivariate_gcd(a, b);
        EXPECT_TRUE(g.is_constant(), "gcd(x+1, x+2) is constant (coprime)");
    }

    TEST_CASE("multivariate_gcd: gcd(x^2 - y^2, x + y) == x + y");
    {
        // x^2 - y^2 = (x+y)(x-y), so gcd(x^2-y^2, x+y) = x+y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms_a = {
            make_term({2, 0}, Rational(1)),   // x^2
            make_term({0, 2}, Rational(-1))   // -y^2
        };
        MultiPoly a(terms_a, vars);

        std::vector<MultiPoly::Term> terms_b = {
            make_term({1, 0}, Rational(1)),   // x
            make_term({0, 1}, Rational(1))    // y
        };
        MultiPoly b(terms_b, vars);

        MultiPoly g = multivariate_gcd(a, b);
        MultiPoly g_prim = g.make_primitive();

        MultiPoly expected_prim = b.make_primitive();

        EXPECT_EQ_STR(g_prim.to_string(), expected_prim.to_string(),
                      "gcd(x^2 - y^2, x + y) == x + y");
    }

    TEST_CASE("multivariate_gcd: gcd(0, f) == f (up to scalar)");
    {
        // gcd(0, f) should be f.make_primitive()
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly zero;

        std::vector<MultiPoly::Term> terms_f = {
            make_term({2, 1}, Rational(3)),  // 3*x^2*y
            make_term({1, 0}, Rational(6))   // 6*x
        };
        MultiPoly f(terms_f, vars);

        MultiPoly g = multivariate_gcd(zero, f);
        MultiPoly g_prim = g.make_primitive();
        MultiPoly f_prim = f.make_primitive();

        EXPECT_EQ_STR(g_prim.to_string(), f_prim.to_string(),
                      "gcd(0, f) == f (up to scalar)");
    }

    TEST_CASE("multivariate_gcd: gcd(f, f) == f (up to scalar)");
    {
        // gcd(f, f) should be f (up to scalar multiple)
        // Use a univariate polynomial to test this basic property
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms_f = {
            make_term({2}, Rational(1)),  // x^2
            make_term({1}, Rational(3)),  // 3*x
            make_term({0}, Rational(2))   // 2
        };
        MultiPoly f(terms_f, vars);

        MultiPoly g = multivariate_gcd(f, f);
        MultiPoly g_prim = g.make_primitive();
        MultiPoly f_prim = f.make_primitive();

        EXPECT_EQ_STR(g_prim.to_string(), f_prim.to_string(),
                      "gcd(f, f) == f (up to scalar)");
    }

    TEST_CASE("multivariate_gcd: gcd of two constants == numeric GCD");
    {
        // gcd(12, 8) should be 4
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly a(Rational(12), vars);
        MultiPoly b(Rational(8), vars);

        MultiPoly g = multivariate_gcd(a, b);
        EXPECT_TRUE(g.is_constant(), "gcd of constants is constant");

        MultiPoly expected(Rational(4), vars);
        EXPECT_EQ_STR(g.to_string(), expected.to_string(),
                      "gcd(12, 8) == 4");
    }


    TEST_CASE("trivial content (x^2 + x + 1, content = 1)");
    {
        // poly = x^2 + x + 1, main_var = "x"
        // All coefficients are constants (1, 1, 1), content = 1
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(1)),  // x^2
            make_term({1, 0}, Rational(1)),  // x
            make_term({0, 0}, Rational(1))   // 1
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for x^2+x+1");
        EXPECT_TRUE(content.is_constant(), "content of x^2+x+1 is trivial (constant)");
    }

    TEST_CASE("non-trivial content (x^2*y + x*y, content = y)");
    {
        // poly = x^2*y + x*y, main_var = "x"
        // coeff(x^2) = y, coeff(x^1) = y → content = gcd(y, y) = y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(1)),  // x^2 * y
            make_term({1, 1}, Rational(1))   // x * y
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for x^2*y + x*y");

        // Verify content divides each coefficient
        // coeff(x^2) = y, coeff(x^1) = y; content = y divides both
        std::vector<std::string> rem_vars = {"y"};
        std::vector<MultiPoly::Term> coeff_x2 = {make_term({1}, Rational(1))};  // y
        std::vector<MultiPoly::Term> coeff_x1 = {make_term({1}, Rational(1))};  // y
        MultiPoly c_x2(coeff_x2, rem_vars);
        MultiPoly c_x1(coeff_x1, rem_vars);
        bool divides_x2 = true, divides_x1 = true;
        try { c_x2.exact_div(content); } catch (...) { divides_x2 = false; }
        try { c_x1.exact_div(content); } catch (...) { divides_x1 = false; }
        EXPECT_TRUE(divides_x2, "content divides coeff(x^2) for x^2*y + x*y");
        EXPECT_TRUE(divides_x1, "content divides coeff(x^1) for x^2*y + x*y");
    }

    TEST_CASE("numeric content (6*x^2 + 4*x, content = 2)");
    {
        // poly = 6*x^2 + 4*x, main_var = "x"
        // coeff(x^2) = 6, coeff(x^1) = 4 → content = gcd(6, 4) = 2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 0}, Rational(6)),  // 6*x^2
            make_term({1, 0}, Rational(4))   // 4*x
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for 6x^2 + 4x");
        EXPECT_TRUE(content.is_constant(), "content of 6x^2+4x is constant");
    }

    TEST_CASE("trivariate polynomial (x^2*y*z + x*y*z^2, content = y*z)");
    {
        // poly = x^2*y*z + x*y*z^2, main_var = "x"
        // coeff(x^2) = y*z, coeff(x^1) = y*z^2 → content = gcd(y*z, y*z^2) = y*z
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1, 1}, Rational(1)),  // x^2 * y * z
            make_term({1, 1, 2}, Rational(1))   // x * y * z^2
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for trivariate poly");

        // Verify content divides each coefficient
        std::vector<std::string> rem_vars = {"y", "z"};
        std::vector<MultiPoly::Term> coeff_x2_terms = {make_term({1, 1}, Rational(1))};  // y*z
        std::vector<MultiPoly::Term> coeff_x1_terms = {make_term({1, 2}, Rational(1))};  // y*z^2
        MultiPoly c_x2(coeff_x2_terms, rem_vars);
        MultiPoly c_x1(coeff_x1_terms, rem_vars);
        bool div2 = true, div1 = true;
        try { c_x2.exact_div(content); } catch (...) { div2 = false; }
        try { c_x1.exact_div(content); } catch (...) { div1 = false; }
        EXPECT_TRUE(div2, "content divides coeff(x^2) for trivariate poly");
        EXPECT_TRUE(div1, "content divides coeff(x^1) for trivariate poly");
    }

    TEST_CASE("mixed numeric and polynomial content (2*x^2*y + 4*x*y^2, content = 2y)");
    {
        // poly = 2*x^2*y + 4*x*y^2, main_var = "x"
        // coeff(x^2) = 2y, coeff(x^1) = 4y^2 → content = gcd(2y, 4y^2) = 2y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            make_term({2, 1}, Rational(2)),  // 2*x^2*y
            make_term({1, 2}, Rational(4))   // 4*x*y^2
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for 2x^2*y + 4x*y^2");

        // Verify content divides each coefficient
        std::vector<std::string> rem_vars = {"y"};
        std::vector<MultiPoly::Term> coeff_x2_terms = {make_term({1}, Rational(2))};  // 2y
        std::vector<MultiPoly::Term> coeff_x1_terms = {make_term({2}, Rational(4))};  // 4y^2
        MultiPoly c_x2(coeff_x2_terms, rem_vars);
        MultiPoly c_x1(coeff_x1_terms, rem_vars);
        bool div2 = true, div1 = true;
        try { c_x2.exact_div(content); } catch (...) { div2 = false; }
        try { c_x1.exact_div(content); } catch (...) { div1 = false; }
        EXPECT_TRUE(div2, "content divides coeff(x^2) for 2x^2*y + 4x*y^2");
        EXPECT_TRUE(div1, "content divides coeff(x^1) for 2x^2*y + 4x*y^2");
    }

    TEST_CASE("single term polynomial (x^3*y^2*z, content = y^2*z)");
    {
        // poly = x^3*y^2*z, main_var = "x"
        // Only one coefficient: coeff(x^3) = y^2*z → content = y^2*z
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            make_term({3, 2, 1}, Rational(1))  // x^3 * y^2 * z
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for single term x^3*y^2*z");
    }

    TEST_CASE("constant polynomial (content = poly itself)");
    {
        // poly = 7, main_var = "x"
        // Constant polynomial has degree 0 in x, content = poly itself
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly f(Rational(7), vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for constant 7");
    }

    TEST_CASE("three-term trivariate (x^3*y + x^2*y*z + x*y*z^2, content = y)");
    {
        // poly = x^3*y + x^2*y*z + x*y*z^2, main_var = "x"
        // coeff(x^3) = y, coeff(x^2) = y*z, coeff(x^1) = y*z^2
        // content = gcd(y, y*z, y*z^2) = y
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            make_term({3, 1, 0}, Rational(1)),  // x^3 * y
            make_term({2, 1, 1}, Rational(1)),  // x^2 * y * z
            make_term({1, 1, 2}, Rational(1))   // x * y * z^2
        };
        MultiPoly f(terms, vars);
        MultiPoly content = multivariate_content(f, "x");
        MultiPoly pp = multivariate_primitive_part(f, "x");
        MultiPoly content_full = embed_content(content, "x", vars);
        MultiPoly product = content_full * pp;
        EXPECT_TRUE(product == f, "content * primitive_part == f for x^3*y + x^2*y*z + x*y*z^2");

        // Verify content divides each coefficient
        std::vector<std::string> rem_vars = {"y", "z"};
        std::vector<MultiPoly::Term> c3_terms = {make_term({1, 0}, Rational(1))};  // y
        std::vector<MultiPoly::Term> c2_terms = {make_term({1, 1}, Rational(1))};  // y*z
        std::vector<MultiPoly::Term> c1_terms = {make_term({1, 2}, Rational(1))};  // y*z^2
        MultiPoly c3(c3_terms, rem_vars);
        MultiPoly c2(c2_terms, rem_vars);
        MultiPoly c1(c1_terms, rem_vars);
        bool d3 = true, d2 = true, d1 = true;
        try { c3.exact_div(content); } catch (...) { d3 = false; }
        try { c2.exact_div(content); } catch (...) { d2 = false; }
        try { c1.exact_div(content); } catch (...) { d1 = false; }
        EXPECT_TRUE(d3, "content divides coeff(x^3) for trivariate 3-term");
        EXPECT_TRUE(d2, "content divides coeff(x^2) for trivariate 3-term");
        EXPECT_TRUE(d1, "content divides coeff(x^1) for trivariate 3-term");
    }

    return TEST_REPORT();
}
