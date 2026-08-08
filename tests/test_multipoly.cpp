/**
 * @file test_multipoly.cpp
 * @brief MultiPoly 类基本操作的单元测试和属性测试。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace lamina;

int main()
{
    TEST_CASE("MultiPoly: zero constructor");
    {
        MultiPoly zero;
        EXPECT_TRUE(zero.is_zero(), "default constructor creates zero polynomial");
        EXPECT_TRUE(zero.terms().empty(), "zero polynomial has no terms");
    }

    TEST_CASE("MultiPoly: constant constructor (non-zero)");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly c(Rational(5), vars);
        EXPECT_FALSE(c.is_zero(), "constant 5 is not zero");
        EXPECT_TRUE(c.num_terms() == 1, "constant has one term");
        // The monomial should be all zeros
        auto& terms = c.terms();
        EXPECT_TRUE(terms[0].first == Monomial({0, 0}), "constant monomial is [0,0]");
        EXPECT_TRUE(terms[0].second == Rational(5), "constant coefficient is 5");
    }

    TEST_CASE("MultiPoly: constant constructor (zero)");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly c(Rational(0), vars);
        EXPECT_TRUE(c.is_zero(), "constant 0 creates zero polynomial");
        EXPECT_TRUE(c.terms().empty(), "zero constant has no terms");
    }

    TEST_CASE("MultiPoly: from terms constructor with normalization");
    {
        std::vector<std::string> vars = {"x", "y"};
        // Create terms: 3xy + 2xy + x^2 (should merge 3xy+2xy = 5xy)
        std::vector<MultiPoly::Term> terms = {
            {{1, 1}, Rational(3)},  // 3xy
            {{1, 1}, Rational(2)},  // 2xy (duplicate monomial)
            {{2, 0}, Rational(1)},  // x^2
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.num_terms() == 2, "merged like terms: 2 terms remain");
        // After normalization with GrevLex: x^2 (degree 2) > xy (degree 2, but
        // GrevLex compares last variable reversed)
        auto& t = p.terms();
        // x^2 has total degree 2, xy has total degree 2
        // GrevLex: compare last index reversed. x^2=[2,0], xy=[1,1]
        // total degrees equal (2==2), compare from last: y-exponent: 0 < 1, so x^2 > xy in GrevLex
        EXPECT_TRUE(t[0].first == Monomial({2, 0}), "first term is x^2");
        EXPECT_TRUE(t[0].second == Rational(1), "x^2 coefficient is 1");
        EXPECT_TRUE(t[1].first == Monomial({1, 1}), "second term is xy");
        EXPECT_TRUE(t[1].second == Rational(5), "xy coefficient is 5 (merged)");
    }

    TEST_CASE("MultiPoly: normalize removes zero-coefficient terms");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 0}, Rational(3)},   // 3x
            {{1, 0}, Rational(-3)},  // -3x (cancels)
            {{0, 1}, Rational(2)},   // 2y
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.num_terms() == 1, "cancelled terms removed: 1 term remains");
        auto& t = p.terms();
        EXPECT_TRUE(t[0].first == Monomial({0, 1}), "remaining term is y");
        EXPECT_TRUE(t[0].second == Rational(2), "y coefficient is 2");
    }

    TEST_CASE("MultiPoly: normalize sorts by monomial order (GrevLex)");
    {
        std::vector<std::string> vars = {"x", "y", "z"};
        // Terms in random order
        std::vector<MultiPoly::Term> terms = {
            {{0, 0, 1}, Rational(1)},  // z (degree 1)
            {{1, 1, 1}, Rational(2)},  // xyz (degree 3)
            {{2, 0, 0}, Rational(3)},  // x^2 (degree 2)
            {{0, 2, 0}, Rational(4)},  // y^2 (degree 2)
        };
        MultiPoly p(terms, vars);

        auto& t = p.terms();
        // GrevLex: highest total degree first
        // xyz (deg 3) > x^2 (deg 2) and y^2 (deg 2) > z (deg 1)
        EXPECT_TRUE(t[0].first == Monomial({1, 1, 1}), "first term is xyz (highest degree)");
        // Among degree-2 terms: x^2=[2,0,0] vs y^2=[0,2,0]
        // GrevLex: same total degree, compare from last reversed: z-exp: 0==0, y-exp: 0<2, so x^2 > y^2
        EXPECT_TRUE(t[1].first == Monomial({2, 0, 0}), "second term is x^2");
        EXPECT_TRUE(t[2].first == Monomial({0, 2, 0}), "third term is y^2");
        EXPECT_TRUE(t[3].first == Monomial({0, 0, 1}), "fourth term is z");
    }

    TEST_CASE("MultiPoly: monomial padding to variable count");
    {
        std::vector<std::string> vars = {"x", "y", "z"};
        // Provide a monomial shorter than vars count
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(7)},  // only x exponent given, should pad to [1,0,0]
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.num_terms() == 1, "one term after padding");
        EXPECT_TRUE(p.terms()[0].first.size() == 3, "monomial padded to 3 components");
        EXPECT_TRUE(p.terms()[0].first == Monomial({1, 0, 0}), "padded monomial is [1,0,0]");
    }

    TEST_CASE("MultiPoly: leading_coeff basic bivariate");
    {
        // p = 3x^2*y + 2x^2 + x*y + 5
        // Viewed as univariate in x: lc(p, x) = 3y + 2 (coefficient of x^2)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 1}, Rational(3)},  // 3x^2*y
            {{2, 0}, Rational(2)},  // 2x^2
            {{1, 1}, Rational(1)},  // xy
            {{0, 0}, Rational(5)},  // 5
        };
        MultiPoly p(terms, vars);

        MultiPoly lc = p.leading_coeff("x");
        // lc should be 3y + 2 in variable y
        EXPECT_TRUE(lc.num_terms() == 2, "leading_coeff(x) has 2 terms");
        EXPECT_TRUE(lc.variables().size() == 1, "leading_coeff(x) has 1 variable");
        EXPECT_TRUE(lc.variables()[0] == "y", "remaining variable is y");
        // Check terms: 3y + 2
        auto& lc_t = lc.terms();
        EXPECT_TRUE(lc_t[0].first == Monomial({1}), "first term monomial is [1] (y)");
        EXPECT_TRUE(lc_t[0].second == Rational(3), "first term coeff is 3");
        EXPECT_TRUE(lc_t[1].first == Monomial({0}), "second term monomial is [0] (constant)");
        EXPECT_TRUE(lc_t[1].second == Rational(2), "second term coeff is 2");
    }

    TEST_CASE("MultiPoly: leading_coeff with respect to y");
    {
        // p = x^2*y^3 + 2x*y^3 + y^2
        // Viewed as univariate in y: lc(p, y) = x^2 + 2x (coefficient of y^3)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 3}, Rational(1)},  // x^2*y^3
            {{1, 3}, Rational(2)},  // 2x*y^3
            {{0, 2}, Rational(1)},  // y^2
        };
        MultiPoly p(terms, vars);

        MultiPoly lc = p.leading_coeff("y");
        // lc should be x^2 + 2x in variable x
        EXPECT_TRUE(lc.num_terms() == 2, "leading_coeff(y) has 2 terms");
        EXPECT_TRUE(lc.variables().size() == 1, "leading_coeff(y) has 1 variable");
        EXPECT_TRUE(lc.variables()[0] == "x", "remaining variable is x");
        auto& lc_t = lc.terms();
        EXPECT_TRUE(lc_t[0].first == Monomial({2}), "first term is x^2");
        EXPECT_TRUE(lc_t[0].second == Rational(1), "x^2 coeff is 1");
        EXPECT_TRUE(lc_t[1].first == Monomial({1}), "second term is x");
        EXPECT_TRUE(lc_t[1].second == Rational(2), "x coeff is 2");
    }

    TEST_CASE("MultiPoly: leading_coeff of constant polynomial");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly c(Rational(7), vars);

        MultiPoly lc = c.leading_coeff("x");
        // Constant poly: degree in x is 0, so lc is the constant itself in y
        EXPECT_TRUE(lc.num_terms() == 1, "leading_coeff of constant has 1 term");
        EXPECT_TRUE(lc.terms()[0].second == Rational(7), "leading_coeff of constant is 7");
    }

    TEST_CASE("MultiPoly: leading_coeff of zero polynomial");
    {
        MultiPoly zero;
        MultiPoly lc = zero.leading_coeff("x");
        EXPECT_TRUE(lc.is_zero(), "leading_coeff of zero is zero");
    }

    TEST_CASE("MultiPoly: leading_coeff with variable not in polynomial");
    {
        // p = 3x + 2 in vars {x, y}
        // leading_coeff("z") where z is not in vars: returns p itself
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 0}, Rational(3)},  // 3x
            {{0, 0}, Rational(2)},  // 2
        };
        MultiPoly p(terms, vars);

        MultiPoly lc = p.leading_coeff("z");
        EXPECT_TRUE(lc == p, "leading_coeff of unknown var returns polynomial itself");
    }

    TEST_CASE("MultiPoly: leading_coeff trivariate");
    {
        // p = 2x^2*y*z + x^2*z^2 + 3x*y
        // Viewed as univariate in x: lc(p, x) = 2yz + z^2 (coefficient of x^2)
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 1, 1}, Rational(2)},  // 2x^2*y*z
            {{2, 0, 2}, Rational(1)},  // x^2*z^2
            {{1, 1, 0}, Rational(3)},  // 3x*y
        };
        MultiPoly p(terms, vars);

        MultiPoly lc = p.leading_coeff("x");
        // lc should be 2yz + z^2 in variables {y, z}
        EXPECT_TRUE(lc.num_terms() == 2, "leading_coeff(x) trivariate has 2 terms");
        EXPECT_TRUE(lc.variables().size() == 2, "remaining vars are y, z");
        EXPECT_TRUE(lc.variables()[0] == "y", "first remaining var is y");
        EXPECT_TRUE(lc.variables()[1] == "z", "second remaining var is z");
    }

    TEST_CASE("MultiPoly: to_univariate basic");
    {
        // p = 3x^2 + 2x + 1 in variable x
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(3)},  // 3x^2
            {{1}, Rational(2)},  // 2x
            {{0}, Rational(1)},  // 1
        };
        MultiPoly p(terms, vars);

        Polynomial<Rational> uni = p.to_univariate();
        EXPECT_TRUE(uni.degree() == 2, "to_univariate degree is 2");
        EXPECT_TRUE(uni.coeffs[0] == Rational(1), "constant term is 1");
        EXPECT_TRUE(uni.coeffs[1] == Rational(2), "x^1 coefficient is 2");
        EXPECT_TRUE(uni.coeffs[2] == Rational(3), "x^2 coefficient is 3");
        EXPECT_TRUE(uni.variable_name == "x", "variable name is x");
    }

    TEST_CASE("MultiPoly: to_univariate with multi-var list but univariate");
    {
        // p = 5y^3 + y in vars {x, y} — only y has non-zero exponents
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{0, 3}, Rational(5)},  // 5y^3
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.is_univariate(), "polynomial is univariate");
        Polynomial<Rational> uni = p.to_univariate();
        EXPECT_TRUE(uni.degree() == 3, "to_univariate degree is 3");
        EXPECT_TRUE(uni.coeffs[0] == Rational(0), "constant term is 0");
        EXPECT_TRUE(uni.coeffs[1] == Rational(1), "y^1 coefficient is 1");
        EXPECT_TRUE(uni.coeffs[2] == Rational(0), "y^2 coefficient is 0");
        EXPECT_TRUE(uni.coeffs[3] == Rational(5), "y^3 coefficient is 5");
        EXPECT_TRUE(uni.variable_name == "y", "variable name is y");
    }

    TEST_CASE("MultiPoly: to_univariate zero polynomial");
    {
        std::vector<std::string> vars = {"x"};
        MultiPoly zero(std::vector<MultiPoly::Term>{}, vars);
        Polynomial<Rational> uni = zero.to_univariate();
        EXPECT_TRUE(uni.is_zero(), "zero MultiPoly converts to zero Polynomial");
    }

    TEST_CASE("MultiPoly: to_univariate constant polynomial");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly c(Rational(7), vars);
        Polynomial<Rational> uni = c.to_univariate();
        EXPECT_TRUE(uni.degree() == 0, "constant converts to degree-0 polynomial");
        EXPECT_TRUE(uni.coeffs[0] == Rational(7), "constant value is 7");
    }

    TEST_CASE("MultiPoly: to_univariate throws for multivariate");
    {
        // p = x + y — not univariate
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly p(terms, vars);

        bool threw = false;
        try {
            p.to_univariate();
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "to_univariate throws for multivariate polynomial");
    }

    TEST_CASE("MultiPoly: from_univariate basic");
    {
        // Create Polynomial<Rational>: 2x^3 + x + 4
        std::vector<Rational> coeffs = {Rational(4), Rational(1), Rational(0), Rational(2)};
        Polynomial<Rational> uni(coeffs, "x");

        MultiPoly mp = MultiPoly::from_univariate(uni, "x");
        EXPECT_TRUE(mp.num_terms() == 3, "from_univariate has 3 non-zero terms");
        EXPECT_TRUE(mp.variables().size() == 1, "from_univariate has 1 variable");
        EXPECT_TRUE(mp.variables()[0] == "x", "variable is x");
        EXPECT_TRUE(mp.total_degree() == 3, "total degree is 3");
    }

    TEST_CASE("MultiPoly: from_univariate zero polynomial");
    {
        Polynomial<Rational> uni("t");  // zero polynomial
        MultiPoly mp = MultiPoly::from_univariate(uni, "t");
        EXPECT_TRUE(mp.is_zero(), "from_univariate of zero is zero");
        EXPECT_TRUE(mp.variables().size() == 1, "zero poly still has variable list");
        EXPECT_TRUE(mp.variables()[0] == "t", "variable name preserved");
    }

    TEST_CASE("MultiPoly: from_univariate constant");
    {
        Polynomial<Rational> uni(Rational(9), "z");
        MultiPoly mp = MultiPoly::from_univariate(uni, "z");
        EXPECT_TRUE(mp.is_constant(), "from_univariate of constant is constant");
        EXPECT_TRUE(mp.num_terms() == 1, "constant has 1 term");
        EXPECT_TRUE(mp.terms()[0].second == Rational(9), "coefficient is 9");
    }

    TEST_CASE("MultiPoly: univariate round-trip");
    {
        // Create a Polynomial<Rational>: 3x^4 - 2x^2 + 7
        std::vector<Rational> coeffs = {Rational(7), Rational(0), Rational(-2), Rational(0), Rational(3)};
        Polynomial<Rational> original(coeffs, "x");

        MultiPoly mp = MultiPoly::from_univariate(original, "x");
        Polynomial<Rational> recovered = mp.to_univariate();

        EXPECT_TRUE(recovered == original, "round-trip from_univariate -> to_univariate preserves polynomial");
    }

    TEST_CASE("MultiPoly: univariate round-trip (MultiPoly first)");
    {
        // Start with MultiPoly: 2x^3 + 5x
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{3}, Rational(2)},  // 2x^3
            {{1}, Rational(5)},  // 5x
        };
        MultiPoly original(terms, vars);

        Polynomial<Rational> uni = original.to_univariate();
        MultiPoly recovered = MultiPoly::from_univariate(uni, "x");

        EXPECT_TRUE(recovered == original, "round-trip to_univariate -> from_univariate preserves polynomial");
    }

    TEST_CASE("MultiPoly: exact_div basic");
    {
        // f = x^2 - y^2, g = x + y => f / g = x - y
        std::vector<std::string> vars = {"x", "y"};
        // f = x^2 - y^2
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0}, Rational(1)},   // x^2
            {{0, 2}, Rational(-1)},  // -y^2
        };
        MultiPoly f(f_terms, vars);

        // g = x + y
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly g(g_terms, vars);

        MultiPoly q = f.exact_div(g);

        // Expected: x - y
        std::vector<MultiPoly::Term> expected_terms = {
            {{1, 0}, Rational(1)},   // x
            {{0, 1}, Rational(-1)},  // -y
        };
        MultiPoly expected(expected_terms, vars);

        EXPECT_TRUE(q == expected, "exact_div: (x^2 - y^2) / (x + y) == x - y");
    }

    TEST_CASE("MultiPoly: exact_div product round-trip");
    {
        // g = 2x + 3y, h = x - y => f = g*h = 2x^2 + xy - 3y^2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(2)},  // 2x
            {{0, 1}, Rational(3)},  // 3y
        };
        MultiPoly g(g_terms, vars);

        std::vector<MultiPoly::Term> h_terms = {
            {{1, 0}, Rational(1)},   // x
            {{0, 1}, Rational(-1)},  // -y
        };
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;

        MultiPoly q1 = f.exact_div(g);
        EXPECT_TRUE(q1 == h, "exact_div: (g*h) / g == h");

        MultiPoly q2 = f.exact_div(h);
        EXPECT_TRUE(q2 == g, "exact_div: (g*h) / h == g");
    }

    TEST_CASE("MultiPoly: exact_div by constant");
    {
        std::vector<std::string> vars = {"x", "y"};
        // f = 6x^2 + 4xy
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0}, Rational(6)},  // 6x^2
            {{1, 1}, Rational(4)},  // 4xy
        };
        MultiPoly f(f_terms, vars);

        // divisor = 2
        MultiPoly divisor(Rational(2), vars);

        MultiPoly q = f.exact_div(divisor);

        // Expected: 3x^2 + 2xy
        std::vector<MultiPoly::Term> expected_terms = {
            {{2, 0}, Rational(3)},  // 3x^2
            {{1, 1}, Rational(2)},  // 2xy
        };
        MultiPoly expected(expected_terms, vars);

        EXPECT_TRUE(q == expected, "exact_div by constant: (6x^2 + 4xy) / 2 == 3x^2 + 2xy");
    }

    TEST_CASE("MultiPoly: exact_div zero dividend");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly zero;
        zero = MultiPoly(std::vector<MultiPoly::Term>{}, vars);

        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},  // x
        };
        MultiPoly g(g_terms, vars);

        MultiPoly q = zero.exact_div(g);
        EXPECT_TRUE(q.is_zero(), "exact_div: 0 / g == 0");
    }

    TEST_CASE("MultiPoly: exact_div throws on zero divisor");
    {
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 0}, Rational(1)},  // x
        };
        MultiPoly f(f_terms, vars);
        MultiPoly zero(std::vector<MultiPoly::Term>{}, vars);

        bool threw = false;
        try {
            f.exact_div(zero);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "exact_div throws on zero divisor");
    }

    TEST_CASE("MultiPoly: exact_div throws on non-exact division");
    {
        std::vector<std::string> vars = {"x", "y"};
        // f = x^2 + 1, g = x + y (does not divide exactly)
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0}, Rational(1)},  // x^2
            {{0, 0}, Rational(1)},  // 1
        };
        MultiPoly f(f_terms, vars);

        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly g(g_terms, vars);

        bool threw = false;
        try {
            f.exact_div(g);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "exact_div throws on non-exact division");
    }


    TEST_CASE("MultiPoly: eval single variable basic");
    {
        // p = 3x^2*y + 2x*y + y  in vars {x, y}
        // eval(x=2): 3*(4)*y + 2*(2)*y + y = 12y + 4y + y = 17y
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 1}, Rational(3)},  // 3x^2*y
            {{1, 1}, Rational(2)},  // 2x*y
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("x", Rational(2));
        EXPECT_TRUE(result.num_vars() == 1, "eval reduces variable count by 1");
        EXPECT_TRUE(result.variables()[0] == "y", "remaining variable is y");
        EXPECT_TRUE(result.num_terms() == 1, "result has 1 term (merged)");
        EXPECT_TRUE(result.terms()[0].second == Rational(17), "coefficient is 17");
        EXPECT_TRUE(result.terms()[0].first == Monomial({1}), "monomial is y^1");
    }

    TEST_CASE("MultiPoly: eval single variable with val=0");
    {
        // p = x^2 + 3x + 5  in vars {x}
        // eval(x=0): 0 + 0 + 5 = 5
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},  // x^2
            {{1}, Rational(3)},  // 3x
            {{0}, Rational(5)},  // 5
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("x", Rational(0));
        EXPECT_TRUE(result.is_constant(), "eval at 0 gives constant");
        EXPECT_TRUE(result.num_terms() == 1, "result has 1 term");
        EXPECT_TRUE(result.terms()[0].second == Rational(5), "constant is 5");
    }

    TEST_CASE("MultiPoly: eval single variable with val=1");
    {
        // p = 2x^3 + x^2 - x + 4  in vars {x}
        // eval(x=1): 2 + 1 - 1 + 4 = 6
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{3}, Rational(2)},   // 2x^3
            {{2}, Rational(1)},   // x^2
            {{1}, Rational(-1)},  // -x
            {{0}, Rational(4)},   // 4
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("x", Rational(1));
        EXPECT_TRUE(result.is_constant(), "eval at 1 gives constant");
        EXPECT_TRUE(result.terms()[0].second == Rational(6), "value is 6");
    }

    TEST_CASE("MultiPoly: eval variable not in polynomial");
    {
        // p = x + y in vars {x, y}
        // eval(z=5): polynomial unchanged
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("z", Rational(5));
        EXPECT_TRUE(result == p, "eval of unknown variable returns polynomial unchanged");
    }

    TEST_CASE("MultiPoly: eval reduces to zero");
    {
        // p = x^2 - 4  in vars {x}
        // eval(x=2): 4 - 4 = 0
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},   // x^2
            {{0}, Rational(-4)},  // -4
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("x", Rational(2));
        EXPECT_TRUE(result.is_zero(), "eval at root gives zero");
    }

    TEST_CASE("MultiPoly: eval with rational value");
    {
        // p = 4x^2 + 2x + 1  in vars {x}
        // eval(x=1/2): 4*(1/4) + 2*(1/2) + 1 = 1 + 1 + 1 = 3
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(4)},  // 4x^2
            {{1}, Rational(2)},  // 2x
            {{0}, Rational(1)},  // 1
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("x", Rational(1, 2));
        EXPECT_TRUE(result.is_constant(), "eval gives constant");
        EXPECT_TRUE(result.terms()[0].second == Rational(3), "value is 3");
    }

    TEST_CASE("MultiPoly: eval trivariate partial substitution");
    {
        // p = x*y*z + x*z + y  in vars {x, y, z}
        // eval(y=3): x*3*z + x*z + 3 = 3xz + xz + 3 = 4xz + 3
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 1, 1}, Rational(1)},  // xyz
            {{1, 0, 1}, Rational(1)},  // xz
            {{0, 1, 0}, Rational(1)},  // y
        };
        MultiPoly p(terms, vars);

        MultiPoly result = p.eval("y", Rational(3));
        EXPECT_TRUE(result.num_vars() == 2, "result has 2 variables");
        EXPECT_TRUE(result.variables()[0] == "x", "first var is x");
        EXPECT_TRUE(result.variables()[1] == "z", "second var is z");
        EXPECT_TRUE(result.num_terms() == 2, "result has 2 terms");
        // 4xz term
        EXPECT_TRUE(result.degree("x") == 1, "degree in x is 1");
        EXPECT_TRUE(result.degree("z") == 1, "degree in z is 1");
    }

    TEST_CASE("MultiPoly: eval multi-variable substitution");
    {
        // p = x^2 + y^2 + z^2  in vars {x, y, z}
        // eval({x=1, y=2, z=3}): 1 + 4 + 9 = 14
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 0, 0}, Rational(1)},  // x^2
            {{0, 2, 0}, Rational(1)},  // y^2
            {{0, 0, 2}, Rational(1)},  // z^2
        };
        MultiPoly p(terms, vars);

        std::map<std::string, Rational> sub = {
            {"x", Rational(1)},
            {"y", Rational(2)},
            {"z", Rational(3)}
        };
        MultiPoly result = p.eval(sub);
        EXPECT_TRUE(result.is_constant(), "full substitution gives constant");
        EXPECT_TRUE(result.terms()[0].second == Rational(14), "value is 14");
    }

    TEST_CASE("MultiPoly: eval multi-variable partial substitution");
    {
        // p = 2x*y + 3y*z + x  in vars {x, y, z}
        // eval({x=1, z=2}): 2*1*y + 3*y*2 + 1 = 2y + 6y + 1 = 8y + 1
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{1, 1, 0}, Rational(2)},  // 2xy
            {{0, 1, 1}, Rational(3)},  // 3yz
            {{1, 0, 0}, Rational(1)},  // x
        };
        MultiPoly p(terms, vars);

        std::map<std::string, Rational> sub = {
            {"x", Rational(1)},
            {"z", Rational(2)}
        };
        MultiPoly result = p.eval(sub);
        EXPECT_TRUE(result.num_vars() == 1, "result has 1 variable");
        EXPECT_TRUE(result.variables()[0] == "y", "remaining variable is y");
        EXPECT_TRUE(result.num_terms() == 2, "result has 2 terms (8y + 1)");
    }

    TEST_CASE("MultiPoly: eval zero polynomial");
    {
        MultiPoly zero;
        MultiPoly result = zero.eval("x", Rational(5));
        EXPECT_TRUE(result.is_zero(), "eval of zero polynomial is zero");
    }

    TEST_CASE("MultiPoly: eval is ring homomorphism (additive)");
    {
        // Verify eval(f + g, x=a) == eval(f, x=a) + eval(g, x=a)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 1}, Rational(1)},  // x^2*y
            {{0, 0}, Rational(3)},  // 3
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 1}, Rational(2)},  // 2xy
            {{0, 1}, Rational(-1)}, // -y
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(2);
        MultiPoly sum_then_eval = (f + g).eval("x", a);
        MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum, "eval is additive homomorphism");
    }

    TEST_CASE("MultiPoly: eval is ring homomorphism (multiplicative)");
    {
        // Verify eval(f * g, x=a) == eval(f, x=a) * eval(g, x=a)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(-1)}, // -y
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(3);
        MultiPoly prod_then_eval = (f * g).eval("x", a);
        MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod, "eval is multiplicative homomorphism");
    }


    TEST_CASE("MultiPoly: numeric_content of zero polynomial");
    {
        MultiPoly zero;
        Rational content = zero.numeric_content();
        EXPECT_TRUE(content == Rational(0), "content of zero poly is 0");
    }

    TEST_CASE("MultiPoly: numeric_content of integer polynomial");
    {
        // 6x^2 + 4x + 2 → content = 2
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(6)},
            {{1}, Rational(4)},
            {{0}, Rational(2)},
        };
        MultiPoly p(terms, vars);
        Rational content = p.numeric_content();
        EXPECT_TRUE(content == Rational(2), "content of 6x^2+4x+2 is 2");
    }

    TEST_CASE("MultiPoly: numeric_content of rational polynomial");
    {
        // (1/2)x + (1/3) → content = gcd(1/2, 1/3) = 1/6
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(1, 2)},
            {{0}, Rational(1, 3)},
        };
        MultiPoly p(terms, vars);
        Rational content = p.numeric_content();
        EXPECT_TRUE(content == Rational(1, 6), "content of (1/2)x+(1/3) is 1/6");
    }

    TEST_CASE("MultiPoly: numeric_content already primitive");
    {
        // x^2 + x + 1 → content = 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},
            {{1}, Rational(1)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        Rational content = p.numeric_content();
        EXPECT_TRUE(content == Rational(1), "content of x^2+x+1 is 1");
    }

    TEST_CASE("MultiPoly: numeric_content with negative coefficients");
    {
        // -6x + 9 → content = 3 (always positive)
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(-6)},
            {{0}, Rational(9)},
        };
        MultiPoly p(terms, vars);
        Rational content = p.numeric_content();
        EXPECT_TRUE(content == Rational(3), "content of -6x+9 is 3");
    }


    TEST_CASE("MultiPoly: make_primitive of zero polynomial");
    {
        MultiPoly zero;
        MultiPoly prim = zero.make_primitive();
        EXPECT_TRUE(prim.is_zero(), "primitive of zero is zero");
    }

    TEST_CASE("MultiPoly: make_primitive basic integer");
    {
        // 6x^2 + 4x + 2 → primitive = 3x^2 + 2x + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(6)},
            {{1}, Rational(4)},
            {{0}, Rational(2)},
        };
        MultiPoly p(terms, vars);
        MultiPoly prim = p.make_primitive();

        EXPECT_TRUE(prim.num_terms() == 3, "primitive has 3 terms");
        EXPECT_TRUE(prim.terms()[0].second == Rational(3), "leading coeff is 3");
        EXPECT_TRUE(prim.terms()[1].second == Rational(2), "middle coeff is 2");
        EXPECT_TRUE(prim.terms()[2].second == Rational(1), "constant is 1");
    }

    TEST_CASE("MultiPoly: make_primitive ensures positive leading coefficient");
    {
        // -6x + 9 → content=3, divide: -2x + 3, leading is negative → negate: 2x - 3
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(-6)},
            {{0}, Rational(9)},
        };
        MultiPoly p(terms, vars);
        MultiPoly prim = p.make_primitive();

        EXPECT_TRUE(prim.terms()[0].second > Rational(0), "leading coefficient is positive");
        EXPECT_TRUE(prim.terms()[0].second == Rational(2), "leading coeff is 2");
        EXPECT_TRUE(prim.terms()[1].second == Rational(-3), "constant is -3");
    }

    TEST_CASE("MultiPoly: make_primitive of already primitive polynomial");
    {
        // x^2 + x + 1 → already primitive, positive leading coeff
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},
            {{1}, Rational(1)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        MultiPoly prim = p.make_primitive();
        EXPECT_TRUE(prim == p, "already primitive polynomial unchanged");
    }

    TEST_CASE("MultiPoly: make_primitive with rational coefficients");
    {
        // (1/2)x + (1/3) → content=1/6, primitive = 3x + 2
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(1, 2)},
            {{0}, Rational(1, 3)},
        };
        MultiPoly p(terms, vars);
        MultiPoly prim = p.make_primitive();

        EXPECT_TRUE(prim.terms()[0].second == Rational(3), "leading coeff is 3");
        EXPECT_TRUE(prim.terms()[1].second == Rational(2), "constant is 2");
    }


    TEST_CASE("MultiPoly: to_string of zero polynomial");
    {
        MultiPoly zero;
        EXPECT_TRUE(zero.to_string() == "0", "zero poly to_string is '0'");
    }

    TEST_CASE("MultiPoly: to_string of constant");
    {
        std::vector<std::string> vars = {"x"};
        MultiPoly c(Rational(42), vars);
        EXPECT_TRUE(c.to_string() == "42", "constant 42 to_string is '42'");
    }

    TEST_CASE("MultiPoly: to_string basic polynomial");
    {
        // 3x^2 + 2x + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(3)},
            {{1}, Rational(2)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        std::string s = p.to_string();
        EXPECT_TRUE(s == "3*x^2 + 2*x + 1", "to_string of 3x^2+2x+1");
    }

    TEST_CASE("MultiPoly: to_string with coefficient 1");
    {
        // x^2 + x + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},
            {{1}, Rational(1)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        std::string s = p.to_string();
        EXPECT_TRUE(s == "x^2 + x + 1", "to_string omits coefficient 1");
    }

    TEST_CASE("MultiPoly: to_string with negative coefficients");
    {
        // x^2 - 2x + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{2}, Rational(1)},
            {{1}, Rational(-2)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        std::string s = p.to_string();
        EXPECT_TRUE(s == "x^2 - 2*x + 1", "to_string handles negative coefficients");
    }

    TEST_CASE("MultiPoly: to_string multivariate");
    {
        // 2x^2*y + 3x*y^2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 1}, Rational(2)},
            {{1, 2}, Rational(3)},
        };
        MultiPoly p(terms, vars);
        std::string s = p.to_string();
        EXPECT_TRUE(s == "2*x^2*y + 3*x*y^2", "to_string multivariate");
    }

    TEST_CASE("MultiPoly: to_string negative leading coefficient");
    {
        // -x + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{1}, Rational(-1)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);
        std::string s = p.to_string();
        EXPECT_TRUE(s == "-x + 1", "to_string with -1 leading coefficient");
    }


    TEST_CASE("Property 2: Additive homomorphism — bivariate pair 1");
    {
        // f = x^2*y + 3, g = 2x*y - y^2, eval x=2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 1}, Rational(1)},  // x^2*y
            {{0, 0}, Rational(3)},  // 3
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 1}, Rational(2)},   // 2xy
            {{0, 2}, Rational(-1)},  // -y^2
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(2);
        MultiPoly sum_then_eval = (f + g).eval("x", a);
        MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum,
            "P2 additive: eval(f+g, x=2) == eval(f,x=2)+eval(g,x=2) [bivariate 1]");
    }

    TEST_CASE("Property 2: Additive homomorphism — bivariate pair 2");
    {
        // f = 5x^3 + x*y, g = -x^3 + 2y^2 + 1, eval y=-1
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{3, 0}, Rational(5)},  // 5x^3
            {{1, 1}, Rational(1)},  // xy
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{3, 0}, Rational(-1)},  // -x^3
            {{0, 2}, Rational(2)},   // 2y^2
            {{0, 0}, Rational(1)},   // 1
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(-1);
        MultiPoly sum_then_eval = (f + g).eval("y", a);
        MultiPoly eval_then_sum = f.eval("y", a) + g.eval("y", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum,
            "P2 additive: eval(f+g, y=-1) == eval(f,y=-1)+eval(g,y=-1) [bivariate 2]");
    }

    TEST_CASE("Property 2: Additive homomorphism — bivariate pair 3 (rational point)");
    {
        // f = 4x^2*y^2, g = 6x*y + 2, eval x=1/2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 2}, Rational(4)},  // 4x^2*y^2
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 1}, Rational(6)},  // 6xy
            {{0, 0}, Rational(2)},  // 2
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(1, 2);
        MultiPoly sum_then_eval = (f + g).eval("x", a);
        MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum,
            "P2 additive: eval(f+g, x=1/2) == eval(f,x=1/2)+eval(g,x=1/2) [bivariate 3]");
    }

    TEST_CASE("Property 2: Additive homomorphism — trivariate pair 4");
    {
        // f = x*y*z + z^2, g = 2x*z - y, eval z=3
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 1, 1}, Rational(1)},  // xyz
            {{0, 0, 2}, Rational(1)},  // z^2
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0, 1}, Rational(2)},   // 2xz
            {{0, 1, 0}, Rational(-1)},  // -y
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(3);
        MultiPoly sum_then_eval = (f + g).eval("z", a);
        MultiPoly eval_then_sum = f.eval("z", a) + g.eval("z", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum,
            "P2 additive: eval(f+g, z=3) == eval(f,z=3)+eval(g,z=3) [trivariate 4]");
    }

    TEST_CASE("Property 2: Additive homomorphism — trivariate pair 5 (eval at 0)");
    {
        // f = x^2 + y^2 + z^2, g = -x^2 + 2y*z + 7, eval x=0
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0, 0}, Rational(1)},  // x^2
            {{0, 2, 0}, Rational(1)},  // y^2
            {{0, 0, 2}, Rational(1)},  // z^2
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{2, 0, 0}, Rational(-1)},  // -x^2
            {{0, 1, 1}, Rational(2)},   // 2yz
            {{0, 0, 0}, Rational(7)},   // 7
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(0);
        MultiPoly sum_then_eval = (f + g).eval("x", a);
        MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
        EXPECT_TRUE(sum_then_eval == eval_then_sum,
            "P2 additive: eval(f+g, x=0) == eval(f,x=0)+eval(g,x=0) [trivariate 5]");
    }

    TEST_CASE("Property 2: Multiplicative homomorphism — bivariate pair 1");
    {
        // f = x + y, g = x - y, eval x=3
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 1}, Rational(1)},  // y
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},   // x
            {{0, 1}, Rational(-1)},  // -y
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(3);
        MultiPoly prod_then_eval = (f * g).eval("x", a);
        MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod,
            "P2 multiplicative: eval(f*g, x=3) == eval(f,x=3)*eval(g,x=3) [bivariate 1]");
    }

    TEST_CASE("Property 2: Multiplicative homomorphism — bivariate pair 2");
    {
        // f = 2x^2 + y, g = x + 3y^2, eval y=2
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0}, Rational(2)},  // 2x^2
            {{0, 1}, Rational(1)},  // y
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},  // x
            {{0, 2}, Rational(3)},  // 3y^2
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(2);
        MultiPoly prod_then_eval = (f * g).eval("y", a);
        MultiPoly eval_then_prod = f.eval("y", a) * g.eval("y", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod,
            "P2 multiplicative: eval(f*g, y=2) == eval(f,y=2)*eval(g,y=2) [bivariate 2]");
    }

    TEST_CASE("Property 2: Multiplicative homomorphism — bivariate pair 3 (rational point)");
    {
        // f = x*y + 1, g = x - y + 2, eval x=2/3
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 1}, Rational(1)},  // xy
            {{0, 0}, Rational(1)},  // 1
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0}, Rational(1)},   // x
            {{0, 1}, Rational(-1)},  // -y
            {{0, 0}, Rational(2)},   // 2
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(2, 3);
        MultiPoly prod_then_eval = (f * g).eval("x", a);
        MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod,
            "P2 multiplicative: eval(f*g, x=2/3) == eval(f,x=2/3)*eval(g,x=2/3) [bivariate 3]");
    }

    TEST_CASE("Property 2: Multiplicative homomorphism — trivariate pair 4");
    {
        // f = x + y + z, g = x*y - z, eval y=-2
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 0, 0}, Rational(1)},  // x
            {{0, 1, 0}, Rational(1)},  // y
            {{0, 0, 1}, Rational(1)},  // z
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 1, 0}, Rational(1)},   // xy
            {{0, 0, 1}, Rational(-1)},  // -z
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(-2);
        MultiPoly prod_then_eval = (f * g).eval("y", a);
        MultiPoly eval_then_prod = f.eval("y", a) * g.eval("y", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod,
            "P2 multiplicative: eval(f*g, y=-2) == eval(f,y=-2)*eval(g,y=-2) [trivariate 4]");
    }

    TEST_CASE("Property 2: Multiplicative homomorphism — trivariate pair 5 (eval at 0)");
    {
        // f = x*z + y^2 + 1, g = 3z - x, eval z=0
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> f_terms = {
            {{1, 0, 1}, Rational(1)},  // xz
            {{0, 2, 0}, Rational(1)},  // y^2
            {{0, 0, 0}, Rational(1)},  // 1
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{0, 0, 1}, Rational(3)},   // 3z
            {{1, 0, 0}, Rational(-1)},  // -x
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);

        Rational a(0);
        MultiPoly prod_then_eval = (f * g).eval("z", a);
        MultiPoly eval_then_prod = f.eval("z", a) * g.eval("z", a);
        EXPECT_TRUE(prod_then_eval == eval_then_prod,
            "P2 multiplicative: eval(f*g, z=0) == eval(f,z=0)*eval(g,z=0) [trivariate 5]");
    }


    TEST_CASE("Property 3: Univariate round-trip — zero polynomial");
    {
        // Zero polynomial: from_univariate(zero).to_univariate() == zero
        Polynomial<Rational> p("x");  // zero polynomial
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves zero polynomial");
    }

    TEST_CASE("Property 3: Univariate round-trip — constant polynomial");
    {
        // Constant: p = 42
        Polynomial<Rational> p(Rational(42), "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves constant polynomial 42");
    }

    TEST_CASE("Property 3: Univariate round-trip — linear polynomial");
    {
        // Linear: p = 3x + 7
        std::vector<Rational> coeffs = {Rational(7), Rational(3)};
        Polynomial<Rational> p(coeffs, "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves linear polynomial 3x + 7");
    }

    TEST_CASE("Property 3: Univariate round-trip — quadratic polynomial");
    {
        // Quadratic: p = 2x^2 - 5x + 1
        std::vector<Rational> coeffs = {Rational(1), Rational(-5), Rational(2)};
        Polynomial<Rational> p(coeffs, "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves quadratic 2x^2 - 5x + 1");
    }

    TEST_CASE("Property 3: Univariate round-trip — high-degree polynomial");
    {
        // High-degree: p = x^7 + 3x^4 - 2x + 9
        std::vector<Rational> coeffs = {Rational(9), Rational(-2), Rational(0), Rational(0),
                                         Rational(3), Rational(0), Rational(0), Rational(1)};
        Polynomial<Rational> p(coeffs, "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves high-degree x^7 + 3x^4 - 2x + 9");
    }

    TEST_CASE("Property 3: Univariate round-trip — rational coefficients");
    {
        // Rational coefficients: p = (1/2)x^3 + (2/3)x - (5/7)
        std::vector<Rational> coeffs = {Rational(-5, 7), Rational(2, 3), Rational(0), Rational(1, 2)};
        Polynomial<Rational> p(coeffs, "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves rational coefficients (1/2)x^3 + (2/3)x - 5/7");
    }

    TEST_CASE("Property 3: Univariate round-trip — monomial (single term)");
    {
        // Monomial: p = 4x^5
        std::vector<Rational> coeffs = {Rational(0), Rational(0), Rational(0), Rational(0), Rational(0), Rational(4)};
        Polynomial<Rational> p(coeffs, "x");
        MultiPoly mp = MultiPoly::from_univariate(p, "x");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves monomial 4x^5");
    }

    TEST_CASE("Property 3: Univariate round-trip — different variable name");
    {
        // Different variable: p = t^2 + t + 1
        std::vector<Rational> coeffs = {Rational(1), Rational(1), Rational(1)};
        Polynomial<Rational> p(coeffs, "t");
        MultiPoly mp = MultiPoly::from_univariate(p, "t");
        Polynomial<Rational> recovered = mp.to_univariate();
        EXPECT_TRUE(recovered == p, "round-trip preserves variable name 't' in t^2 + t + 1");
    }

    TEST_CASE("Property 3: Univariate round-trip reverse — MultiPoly to univariate and back");
    {
        // Start from MultiPoly: 5x^3 - x^2 + 2
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{3}, Rational(5)},   // 5x^3
            {{2}, Rational(-1)},  // -x^2
            {{0}, Rational(2)},   // 2
        };
        MultiPoly original(terms, vars);

        Polynomial<Rational> uni = original.to_univariate();
        MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
        EXPECT_TRUE(recovered == original, "reverse round-trip preserves 5x^3 - x^2 + 2");
    }

    TEST_CASE("Property 3: Univariate round-trip reverse — constant MultiPoly");
    {
        // Constant MultiPoly: 13
        std::vector<std::string> vars = {"x"};
        MultiPoly original(Rational(13), vars);

        Polynomial<Rational> uni = original.to_univariate();
        MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
        EXPECT_TRUE(recovered == original, "reverse round-trip preserves constant MultiPoly 13");
    }

    TEST_CASE("Property 3: Univariate round-trip reverse — zero MultiPoly");
    {
        // Zero MultiPoly
        std::vector<std::string> vars = {"x"};
        MultiPoly original(std::vector<MultiPoly::Term>{}, vars);

        Polynomial<Rational> uni = original.to_univariate();
        MultiPoly recovered = MultiPoly::from_univariate(uni, "x");
        EXPECT_TRUE(recovered == original, "reverse round-trip preserves zero MultiPoly");
    }


    TEST_CASE("Property 5: exact_div round-trip — simple monomials (x, y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        // g = x, h = y → f = xy
        std::vector<MultiPoly::Term> g_terms = {{{1, 0}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{0, 1}, Rational(1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (x*y)/x == y");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (x*y)/y == x");
    }

    TEST_CASE("Property 5: exact_div round-trip — binomials (x+1, x-1)");
    {
        std::vector<std::string> vars = {"x"};
        // g = x+1, h = x-1 → f = x²-1
        std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(1)}, {{0}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(1)}, {{0}, Rational(-1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (x^2-1)/(x+1) == x-1");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (x^2-1)/(x-1) == x+1");
    }

    TEST_CASE("Property 5: exact_div round-trip — bivariate (x+y, x-y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        // g = x+y, h = x-y → f = x²-y²
        std::vector<MultiPoly::Term> g_terms = {{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{1, 0}, Rational(1)}, {{0, 1}, Rational(-1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (x^2-y^2)/(x+y) == x-y");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (x^2-y^2)/(x-y) == x+y");
    }

    TEST_CASE("Property 5: exact_div round-trip — with coefficients (2x+1, 3x+2)");
    {
        std::vector<std::string> vars = {"x"};
        // g = 2x+1, h = 3x+2 → f = 6x²+7x+2
        std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(2)}, {{0}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(3)}, {{0}, Rational(2)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (6x^2+7x+2)/(2x+1) == 3x+2");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (6x^2+7x+2)/(3x+2) == 2x+1");
    }

    TEST_CASE("Property 5: exact_div round-trip — higher degree (x^2+x+1, x-1)");
    {
        std::vector<std::string> vars = {"x"};
        // g = x^2+x+1, h = x-1 → f = x³-1
        std::vector<MultiPoly::Term> g_terms = {
            {{2}, Rational(1)}, {{1}, Rational(1)}, {{0}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(1)}, {{0}, Rational(-1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (x^3-1)/(x^2+x+1) == x-1");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (x^3-1)/(x-1) == x^2+x+1");
    }

    TEST_CASE("Property 5: exact_div round-trip — trivariate (x+y+z, x-y)");
    {
        std::vector<std::string> vars = {"x", "y", "z"};
        // g = x+y+z, h = x-y
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(1)}, {{0, 0, 1}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {
            {{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(-1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: ((x+y+z)*(x-y))/(x+y+z) == x-y");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: ((x+y+z)*(x-y))/(x-y) == x+y+z");
    }

    TEST_CASE("Property 5: exact_div round-trip — monomial times polynomial (xy, x^2+y)");
    {
        std::vector<std::string> vars = {"x", "y"};
        // g = xy, h = x^2+y → f = x^3*y + x*y^2
        std::vector<MultiPoly::Term> g_terms = {{{1, 1}, Rational(1)}};
        std::vector<MultiPoly::Term> h_terms = {{{2, 0}, Rational(1)}, {{0, 1}, Rational(1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: (x^3*y+x*y^2)/(xy) == x^2+y");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: (x^3*y+x*y^2)/(x^2+y) == xy");
    }

    TEST_CASE("Property 5: exact_div round-trip — rational coefficients (x/2+1/3, 3x-1)");
    {
        std::vector<std::string> vars = {"x"};
        // g = (1/2)x + 1/3, h = 3x - 1
        std::vector<MultiPoly::Term> g_terms = {{{1}, Rational(1, 2)}, {{0}, Rational(1, 3)}};
        std::vector<MultiPoly::Term> h_terms = {{{1}, Rational(3)}, {{0}, Rational(-1)}};
        MultiPoly g(g_terms, vars);
        MultiPoly h(h_terms, vars);

        MultiPoly f = g * h;
        EXPECT_TRUE(f.exact_div(g) == h, "exact_div: f/(x/2+1/3) == 3x-1");
        EXPECT_TRUE(f.exact_div(h) == g, "exact_div: f/(3x-1) == x/2+1/3");
    }


    TEST_CASE("Property 4: degree of zero polynomial is -1");
    {
        MultiPoly zero;
        EXPECT_TRUE(zero.total_degree() == -1, "zero poly total_degree is -1");
        EXPECT_TRUE(zero.degree("x") == -1, "zero poly degree(x) is -1");
        EXPECT_TRUE(zero.degree("y") == -1, "zero poly degree(y) is -1");
    }

    TEST_CASE("Property 4: degree of constant polynomial is 0");
    {
        std::vector<std::string> vars = {"x", "y"};
        MultiPoly c(Rational(42), vars);
        EXPECT_TRUE(c.total_degree() == 0, "constant poly total_degree is 0");
        EXPECT_TRUE(c.degree("x") == 0, "constant poly degree(x) is 0");
        EXPECT_TRUE(c.degree("y") == 0, "constant poly degree(y) is 0");
    }

    TEST_CASE("Property 4: degree of univariate polynomial");
    {
        // p = 5x^4 + 3x^2 + 1
        std::vector<std::string> vars = {"x"};
        std::vector<MultiPoly::Term> terms = {
            {{4}, Rational(5)},
            {{2}, Rational(3)},
            {{0}, Rational(1)},
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 4, "5x^4+3x^2+1 total_degree is 4");
        EXPECT_TRUE(p.degree("x") == 4, "5x^4+3x^2+1 degree(x) is 4");
    }

    TEST_CASE("Property 4: degree of bivariate polynomial");
    {
        // p = 2x^3*y^2 + x*y^5 + 3x^2
        // term degrees: 3+2=5, 1+5=6, 2+0=2
        // total_degree = 6 (from x*y^5)
        // degree(x) = 3 (from x^3*y^2)
        // degree(y) = 5 (from x*y^5)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{3, 2}, Rational(2)},  // 2x^3*y^2
            {{1, 5}, Rational(1)},  // x*y^5
            {{2, 0}, Rational(3)},  // 3x^2
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 6, "2x^3y^2+xy^5+3x^2 total_degree is 6");
        EXPECT_TRUE(p.degree("x") == 3, "2x^3y^2+xy^5+3x^2 degree(x) is 3");
        EXPECT_TRUE(p.degree("y") == 5, "2x^3y^2+xy^5+3x^2 degree(y) is 5");
    }

    TEST_CASE("Property 4: degree of trivariate polynomial");
    {
        // p = x^2*y*z^3 + x*y^4*z + z^7
        // term degrees: 2+1+3=6, 1+4+1=6, 0+0+7=7
        // total_degree = 7 (from z^7)
        // degree(x) = 2 (from x^2*y*z^3)
        // degree(y) = 4 (from x*y^4*z)
        // degree(z) = 7 (from z^7)
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{2, 1, 3}, Rational(1)},  // x^2*y*z^3
            {{1, 4, 1}, Rational(1)},  // x*y^4*z
            {{0, 0, 7}, Rational(1)},  // z^7
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 7, "x^2yz^3+xy^4z+z^7 total_degree is 7");
        EXPECT_TRUE(p.degree("x") == 2, "x^2yz^3+xy^4z+z^7 degree(x) is 2");
        EXPECT_TRUE(p.degree("y") == 4, "x^2yz^3+xy^4z+z^7 degree(y) is 4");
        EXPECT_TRUE(p.degree("z") == 7, "x^2yz^3+xy^4z+z^7 degree(z) is 7");
    }

    TEST_CASE("Property 4: degree with high-degree single term");
    {
        // p = x^10*y^8*z^6 (single term)
        // total_degree = 10+8+6 = 24
        std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<MultiPoly::Term> terms = {
            {{10, 8, 6}, Rational(7)},
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 24, "x^10*y^8*z^6 total_degree is 24");
        EXPECT_TRUE(p.degree("x") == 10, "x^10*y^8*z^6 degree(x) is 10");
        EXPECT_TRUE(p.degree("y") == 8, "x^10*y^8*z^6 degree(y) is 8");
        EXPECT_TRUE(p.degree("z") == 6, "x^10*y^8*z^6 degree(z) is 6");
    }

    TEST_CASE("Property 4: degree for variable not appearing in polynomial");
    {
        // p = x^3 + 1 in vars {x, y}
        // degree(y) should be 0 (y appears with exponent 0 in all terms)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{3, 0}, Rational(1)},  // x^3
            {{0, 0}, Rational(1)},  // 1
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 3, "x^3+1 total_degree is 3");
        EXPECT_TRUE(p.degree("x") == 3, "x^3+1 degree(x) is 3");
        EXPECT_TRUE(p.degree("y") == 0, "x^3+1 degree(y) is 0 (y not used)");
    }

    TEST_CASE("Property 4: degree consistency with arithmetic");
    {
        // Verify: deg(f*g) == deg(f) + deg(g) for total degree
        // f = x^2 + y, g = x*y + 1
        // f*g = x^3*y + x^2 + x*y^2 + y
        // total_degree(f) = 2, total_degree(g) = 2, total_degree(f*g) should be <= 4
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> f_terms = {
            {{2, 0}, Rational(1)},  // x^2
            {{0, 1}, Rational(1)},  // y
        };
        std::vector<MultiPoly::Term> g_terms = {
            {{1, 1}, Rational(1)},  // xy
            {{0, 0}, Rational(1)},  // 1
        };
        MultiPoly f(f_terms, vars);
        MultiPoly g(g_terms, vars);
        MultiPoly fg = f * g;

        int deg_f = f.total_degree();
        int deg_g = g.total_degree();
        int deg_fg = fg.total_degree();

        EXPECT_TRUE(deg_f == 2, "f=x^2+y total_degree is 2");
        EXPECT_TRUE(deg_g == 2, "g=xy+1 total_degree is 2");
        // For polynomials over a field, deg(f*g) = deg(f) + deg(g)
        EXPECT_TRUE(deg_fg == deg_f + deg_g, "deg(f*g) == deg(f) + deg(g)");
    }

    TEST_CASE("Property 4: degree of homogeneous polynomial");
    {
        // p = x^3 + x^2*y + x*y^2 + y^3 (all terms have total degree 3)
        std::vector<std::string> vars = {"x", "y"};
        std::vector<MultiPoly::Term> terms = {
            {{3, 0}, Rational(1)},  // x^3
            {{2, 1}, Rational(1)},  // x^2*y
            {{1, 2}, Rational(1)},  // x*y^2
            {{0, 3}, Rational(1)},  // y^3
        };
        MultiPoly p(terms, vars);

        EXPECT_TRUE(p.total_degree() == 3, "homogeneous poly total_degree is 3");
        EXPECT_TRUE(p.degree("x") == 3, "homogeneous poly degree(x) is 3");
        EXPECT_TRUE(p.degree("y") == 3, "homogeneous poly degree(y) is 3");
        EXPECT_TRUE(p.is_homogeneous(), "polynomial is homogeneous");
    }


    TEST_CASE("Property: addition commutativity (a + b == b + a)");
    {
        std::vector<std::string> vars = {"x", "y"};

        // Case 1: simple linear polynomials
        {
            MultiPoly a({{{1, 0}, Rational(3)}, {{0, 1}, Rational(2)}}, vars);
            MultiPoly b({{{1, 0}, Rational(-1)}, {{0, 0}, Rational(5)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: (3x+2y) + (-x+5) == (-x+5) + (3x+2y)");
        }
        // Case 2: zero polynomial
        {
            MultiPoly a({{{2, 1}, Rational(1)}, {{0, 0}, Rational(7)}}, vars);
            MultiPoly zero;
            EXPECT_TRUE(a + zero == zero + a, "commutativity: a + 0 == 0 + a");
        }
        // Case 3: same polynomial (a + a == a + a trivially, but tests path)
        {
            MultiPoly a({{{1, 1}, Rational(4)}, {{2, 0}, Rational(-3)}}, vars);
            EXPECT_TRUE(a + a == a + a, "commutativity: a + a == a + a");
        }
        // Case 4: polynomials with many terms
        {
            MultiPoly a({{{3, 0}, Rational(1)}, {{2, 1}, Rational(-2)}, {{1, 2}, Rational(3)}, {{0, 3}, Rational(-4)}}, vars);
            MultiPoly b({{{0, 0}, Rational(1)}, {{1, 0}, Rational(2)}, {{0, 1}, Rational(3)}, {{1, 1}, Rational(4)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: dense cubic + dense linear");
        }
        // Case 5: polynomials with rational coefficients
        {
            MultiPoly a({{{1, 0}, Rational(1, 3)}, {{0, 1}, Rational(2, 7)}}, vars);
            MultiPoly b({{{1, 0}, Rational(5, 6)}, {{0, 0}, Rational(-1, 2)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: rational coefficients");
        }
        // Case 6: trivariate
        {
            std::vector<std::string> vars3 = {"x", "y", "z"};
            MultiPoly a({{{1, 1, 1}, Rational(2)}, {{2, 0, 0}, Rational(1)}}, vars3);
            MultiPoly b({{{0, 0, 2}, Rational(3)}, {{1, 1, 0}, Rational(-1)}}, vars3);
            EXPECT_TRUE(a + b == b + a, "commutativity: trivariate");
        }
        // Case 7: cancellation scenario
        {
            MultiPoly a({{{2, 0}, Rational(5)}, {{1, 1}, Rational(3)}}, vars);
            MultiPoly b({{{2, 0}, Rational(-5)}, {{0, 2}, Rational(1)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: with cancellation");
        }
        // Case 8: constant polynomials
        {
            MultiPoly a(Rational(42), vars);
            MultiPoly b(Rational(-17), vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: constants");
        }
        // Case 9: high degree
        {
            MultiPoly a({{{5, 0}, Rational(1)}, {{0, 5}, Rational(1)}}, vars);
            MultiPoly b({{{3, 2}, Rational(2)}, {{2, 3}, Rational(-2)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: high degree");
        }
        // Case 10: negative coefficients
        {
            MultiPoly a({{{1, 0}, Rational(-7)}, {{0, 1}, Rational(-3)}, {{0, 0}, Rational(-1)}}, vars);
            MultiPoly b({{{1, 0}, Rational(7)}, {{0, 1}, Rational(3)}, {{0, 0}, Rational(1)}}, vars);
            EXPECT_TRUE(a + b == b + a, "commutativity: negatives (sum is zero)");
        }
    }

    TEST_CASE("Property: multiplication associativity ((a*b)*c == a*(b*c))");
    {
        std::vector<std::string> vars = {"x", "y"};

        // Case 1: simple linear polynomials
        {
            MultiPoly a({{{1, 0}, Rational(1)}, {{0, 0}, Rational(1)}}, vars);  // x + 1
            MultiPoly b({{{0, 1}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); // y - 1
            MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  // x + y
            EXPECT_TRUE((a * b) * c == a * (b * c), "associativity: (x+1)(y-1)(x+y)");
        }
        // Case 2: one factor is constant
        {
            MultiPoly a({{{1, 0}, Rational(2)}, {{0, 1}, Rational(3)}}, vars);  // 2x + 3y
            MultiPoly b(Rational(5), vars);                                      // 5
            MultiPoly c({{{1, 1}, Rational(1)}, {{0, 0}, Rational(-2)}}, vars); // xy - 2
            EXPECT_TRUE((a * b) * c == a * (b * c), "associativity: constant factor");
        }
        // Case 3: one factor is 1 (identity)
        {
            MultiPoly a({{{2, 0}, Rational(1)}, {{0, 2}, Rational(-1)}}, vars); // x^2 - y^2
            MultiPoly one(Rational(1), vars);
            MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  // x + y
            EXPECT_TRUE((a * one) * c == a * (one * c), "associativity: identity element");
        }
        // Case 4: quadratic factors
        {
            MultiPoly a({{{2, 0}, Rational(1)}, {{0, 0}, Rational(1)}}, vars);  // x^2 + 1
            MultiPoly b({{{0, 2}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); // y^2 - 1
            MultiPoly c({{{1, 1}, Rational(1)}}, vars);                          // xy
            EXPECT_TRUE((a * b) * c == a * (b * c), "associativity: quadratic factors");
        }
        // Case 5: rational coefficients
        {
            MultiPoly a({{{1, 0}, Rational(1, 2)}, {{0, 0}, Rational(1, 3)}}, vars);
            MultiPoly b({{{0, 1}, Rational(2, 3)}, {{0, 0}, Rational(3, 4)}}, vars);
            MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(-1)}}, vars);
            EXPECT_TRUE((a * b) * c == a * (b * c), "associativity: rational coefficients");
        }
        // Case 6: trivariate
        {
            std::vector<std::string> vars3 = {"x", "y", "z"};
            MultiPoly a({{{1, 0, 0}, Rational(1)}, {{0, 0, 1}, Rational(1)}}, vars3); // x + z
            MultiPoly b({{{0, 1, 0}, Rational(1)}, {{0, 0, 0}, Rational(2)}}, vars3); // y + 2
            MultiPoly c({{{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(-1)}}, vars3); // x - y
            EXPECT_TRUE((a * b) * c == a * (b * c), "associativity: trivariate");
        }
        // Case 7: zero factor
        {
            MultiPoly a({{{1, 0}, Rational(3)}, {{0, 1}, Rational(2)}}, vars);
            MultiPoly zero;
            MultiPoly c({{{1, 0}, Rational(1)}}, vars);
            EXPECT_TRUE((a * zero) * c == a * (zero * c), "associativity: zero factor");
        }
    }

    TEST_CASE("Property: distributive law (a * (b + c) == a*b + a*c)");
    {
        std::vector<std::string> vars = {"x", "y"};

        // Case 1: simple linear
        {
            MultiPoly a({{{1, 0}, Rational(2)}, {{0, 0}, Rational(1)}}, vars);  // 2x + 1
            MultiPoly b({{{0, 1}, Rational(1)}, {{0, 0}, Rational(3)}}, vars);  // y + 3
            MultiPoly c({{{1, 0}, Rational(-1)}, {{0, 1}, Rational(2)}}, vars); // -x + 2y
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: (2x+1)*((y+3)+(-x+2y))");
        }
        // Case 2: a is monomial
        {
            MultiPoly a({{{1, 1}, Rational(3)}}, vars);                          // 3xy
            MultiPoly b({{{2, 0}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); // x^2 - 1
            MultiPoly c({{{0, 2}, Rational(2)}, {{1, 0}, Rational(1)}}, vars);  // 2y^2 + x
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: monomial * sum");
        }
        // Case 3: b + c cancels to zero
        {
            MultiPoly a({{{1, 0}, Rational(5)}, {{0, 1}, Rational(7)}}, vars);
            MultiPoly b({{{2, 1}, Rational(3)}, {{1, 0}, Rational(1)}}, vars);
            MultiPoly c = -b;
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: b + c = 0");
        }
        // Case 4: a is constant
        {
            MultiPoly a(Rational(7), vars);
            MultiPoly b({{{1, 0}, Rational(1)}, {{0, 1}, Rational(-2)}}, vars);
            MultiPoly c({{{2, 0}, Rational(3)}, {{0, 0}, Rational(4)}}, vars);
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: constant * sum");
        }
        // Case 5: all quadratic
        {
            MultiPoly a({{{2, 0}, Rational(1)}, {{1, 1}, Rational(-1)}, {{0, 0}, Rational(2)}}, vars);
            MultiPoly b({{{0, 2}, Rational(1)}, {{1, 0}, Rational(3)}}, vars);
            MultiPoly c({{{2, 0}, Rational(-1)}, {{0, 1}, Rational(4)}}, vars);
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: quadratic * sum of quadratics");
        }
        // Case 6: rational coefficients
        {
            MultiPoly a({{{1, 0}, Rational(1, 2)}, {{0, 0}, Rational(3, 4)}}, vars);
            MultiPoly b({{{0, 1}, Rational(2, 3)}, {{0, 0}, Rational(1, 5)}}, vars);
            MultiPoly c({{{1, 0}, Rational(4, 7)}, {{0, 1}, Rational(-1, 3)}}, vars);
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: rational coefficients");
        }
        // Case 7: trivariate
        {
            std::vector<std::string> vars3 = {"x", "y", "z"};
            MultiPoly a({{{1, 0, 1}, Rational(1)}, {{0, 1, 0}, Rational(2)}}, vars3);
            MultiPoly b({{{1, 1, 0}, Rational(1)}, {{0, 0, 1}, Rational(-3)}}, vars3);
            MultiPoly c({{{0, 0, 2}, Rational(2)}, {{1, 0, 0}, Rational(1)}}, vars3);
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: trivariate");
        }
        // Case 8: right distributive (b + c) * a == b*a + c*a
        {
            MultiPoly a({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  // x + y
            MultiPoly b({{{2, 0}, Rational(1)}}, vars);                          // x^2
            MultiPoly c({{{0, 2}, Rational(1)}}, vars);                          // y^2
            EXPECT_TRUE((b + c) * a == b * a + c * a, "right distributive: (x^2+y^2)*(x+y)");
        }
        // Case 9: high degree
        {
            MultiPoly a({{{3, 0}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); // x^3 - 1
            MultiPoly b({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  // x + y
            MultiPoly c({{{0, 1}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); // y - 1
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: cubic * sum");
        }
        // Case 10: single-term polynomials
        {
            MultiPoly a({{{1, 0}, Rational(1)}}, vars);  // x
            MultiPoly b({{{0, 1}, Rational(1)}}, vars);  // y
            MultiPoly c({{{0, 0}, Rational(1)}}, vars);  // 1
            EXPECT_TRUE(a * (b + c) == a * b + a * c, "distributive: x*(y+1) == xy + x");
        }
    }

    return TEST_REPORT();
}
