/**
 * @file test_multipoly_evaluation.cpp
 * @brief MultiPoly 求值契约。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyEvaluation, EvalSingleVariableBasic) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 1}, Rational(3)}, /**< 单项式 3x^2*y。 */
        {{1, 1}, Rational(2)}, /**< 单项式 2x*y。 */
        {{0, 1}, Rational(1)}, /**< 单项式 y。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("x", Rational(2));
    EXPECT_TRUE((result.num_vars() == 1)) << "eval reduces variable count by 1";
    EXPECT_TRUE((result.variables()[0] == "y")) << "remaining variable is y";
    EXPECT_TRUE((result.num_terms() == 1)) << "result has 1 term (merged)";
    EXPECT_TRUE((result.terms()[0].second == Rational(17))) << "coefficient is 17";
    EXPECT_TRUE((result.terms()[0].first == Monomial({1}))) << "monomial is y^1";
}

TEST(MultipolyEvaluation, EvalSingleVariableWithVal0) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)}, /**< 单项式 x^2。 */
        {{1}, Rational(3)}, /**< 单项式 3x。 */
        {{0}, Rational(5)}, /**< 单项式 5。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("x", Rational(0));
    EXPECT_TRUE((result.is_constant())) << "eval at 0 gives constant";
    EXPECT_TRUE((result.num_terms() == 1)) << "result has 1 term";
    EXPECT_TRUE((result.terms()[0].second == Rational(5))) << "constant is 5";
}

TEST(MultipolyEvaluation, EvalSingleVariableWithVal1) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{3}, Rational(2)},  /**< 单项式 2x^3。 */
        {{2}, Rational(1)},  /**< 单项式 x^2。 */
        {{1}, Rational(-1)}, /**< 单项式 -x。 */
        {{0}, Rational(4)},  /**< 单项式 4。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("x", Rational(1));
    EXPECT_TRUE((result.is_constant())) << "eval at 1 gives constant";
    EXPECT_TRUE((result.terms()[0].second == Rational(6))) << "value is 6";
}

TEST(MultipolyEvaluation, EvalVariableNotInPolynomial) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 0}, Rational(1)}, /**< 单项式 x。 */
        {{0, 1}, Rational(1)}, /**< 单项式 y。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("z", Rational(5));
    EXPECT_TRUE((result == p)) << "eval of unknown variable returns polynomial unchanged";
}

TEST(MultipolyEvaluation, EvalReducesToZero) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(1)},  /**< 单项式 x^2。 */
        {{0}, Rational(-4)}, /**< 单项式 -4。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("x", Rational(2));
    EXPECT_TRUE((result.is_zero())) << "eval at root gives zero";
}

TEST(MultipolyEvaluation, EvalWithRationalValue) {
    std::vector<std::string> vars = {"x"};
    std::vector<MultiPoly::Term> terms = {
        {{2}, Rational(4)}, /**< 单项式 4x^2。 */
        {{1}, Rational(2)}, /**< 单项式 2x。 */
        {{0}, Rational(1)}, /**< 单项式 1。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("x", Rational(1, 2));
    EXPECT_TRUE((result.is_constant())) << "eval gives constant";
    EXPECT_TRUE((result.terms()[0].second == Rational(3))) << "value is 3";
}

TEST(MultipolyEvaluation, EvalTrivariatePartialSubstitution) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 1, 1}, Rational(1)}, /**< 单项式 xyz。 */
        {{1, 0, 1}, Rational(1)}, /**< 单项式 xz。 */
        {{0, 1, 0}, Rational(1)}, /**< 单项式 y。 */
    };
    MultiPoly p(terms, vars);

    MultiPoly result = p.eval("y", Rational(3));
    EXPECT_TRUE((result.num_vars() == 2)) << "result has 2 variables";
    EXPECT_TRUE((result.variables()[0] == "x")) << "first var is x";
    EXPECT_TRUE((result.variables()[1] == "z")) << "second var is z";
    EXPECT_TRUE((result.num_terms() == 2)) << "result has 2 terms";
    EXPECT_TRUE((result.degree("x") == 1)) << "degree in x is 1";
    EXPECT_TRUE((result.degree("z") == 1)) << "degree in z is 1";
}

TEST(MultipolyEvaluation, EvalMultiVariableSubstitution) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{2, 0, 0}, Rational(1)}, /**< 单项式 x^2。 */
        {{0, 2, 0}, Rational(1)}, /**< 单项式 y^2。 */
        {{0, 0, 2}, Rational(1)}, /**< 单项式 z^2。 */
    };
    MultiPoly p(terms, vars);

    std::map<std::string, Rational> sub = {
        {"x", Rational(1)},
        {"y", Rational(2)},
        {"z", Rational(3)}};
    MultiPoly result = p.eval(sub);
    EXPECT_TRUE((result.is_constant())) << "full substitution gives constant";
    EXPECT_TRUE((result.terms()[0].second == Rational(14))) << "value is 14";
}

TEST(MultipolyEvaluation, EvalMultiVariablePartialSubstitution) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        {{1, 1, 0}, Rational(2)}, /**< 单项式 2xy。 */
        {{0, 1, 1}, Rational(3)}, /**< 单项式 3yz。 */
        {{1, 0, 0}, Rational(1)}, /**< 单项式 x。 */
    };
    MultiPoly p(terms, vars);

    std::map<std::string, Rational> sub = {
        {"x", Rational(1)},
        {"z", Rational(2)}};
    MultiPoly result = p.eval(sub);
    EXPECT_TRUE((result.num_vars() == 1)) << "result has 1 variable";
    EXPECT_TRUE((result.variables()[0] == "y")) << "remaining variable is y";
    EXPECT_TRUE((result.num_terms() == 2)) << "result has 2 terms (8y + 1)";
}

TEST(MultipolyEvaluation, EvalZeroPolynomial) {
    MultiPoly zero;
    MultiPoly result = zero.eval("x", Rational(5));
    EXPECT_TRUE((result.is_zero())) << "eval of zero polynomial is zero";
}

TEST(MultipolyEvaluation, EvalIsRingHomomorphismAdditive) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 1}, Rational(1)}, /**< 单项式 x^2*y。 */
        {{0, 0}, Rational(3)}, /**< 单项式 3。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 1}, Rational(2)},  /**< 单项式 2xy。 */
        {{0, 1}, Rational(-1)}, /**< 单项式 -y。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(2);
    MultiPoly sum_then_eval = (f + g).eval("x", a);
    MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "eval is additive homomorphism";
}

TEST(MultipolyEvaluation, EvalIsRingHomomorphismMultiplicative) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 0}, Rational(1)}, /**< 单项式 x。 */
        {{0, 1}, Rational(1)}, /**< 单项式 y。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)},  /**< 单项式 x。 */
        {{0, 1}, Rational(-1)}, /**< 单项式 -y。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(3);
    MultiPoly prod_then_eval = (f * g).eval("x", a);
    MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "eval is multiplicative homomorphism";
}

TEST(MultipolyEvaluation, AdditiveHomomorphismBivariateX) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 1}, Rational(1)}, /**< 单项式 x^2*y。 */
        {{0, 0}, Rational(3)}, /**< 单项式 3。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 1}, Rational(2)},  /**< 单项式 2xy。 */
        {{0, 2}, Rational(-1)}, /**< 单项式 -y^2。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(2);
    MultiPoly sum_then_eval = (f + g).eval("x", a);
    MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "P2 additive: eval(f+g, x=2) == eval(f,x=2)+eval(g,x=2) [bivariate 1]";
}

TEST(MultipolyEvaluation, AdditiveHomomorphismBivariateYNegative) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{3, 0}, Rational(5)}, /**< 单项式 5x^3。 */
        {{1, 1}, Rational(1)}, /**< 单项式 xy。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{3, 0}, Rational(-1)}, /**< 单项式 -x^3。 */
        {{0, 2}, Rational(2)},  /**< 单项式 2y^2。 */
        {{0, 0}, Rational(1)},  /**< 单项式 1。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(-1);
    MultiPoly sum_then_eval = (f + g).eval("y", a);
    MultiPoly eval_then_sum = f.eval("y", a) + g.eval("y", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "P2 additive: eval(f+g, y=-1) == eval(f,y=-1)+eval(g,y=-1) [bivariate 2]";
}

TEST(MultipolyEvaluation, AdditiveHomomorphismBivariateRational) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 2}, Rational(4)}, /**< 单项式 4x^2*y^2。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 1}, Rational(6)}, /**< 单项式 6xy。 */
        {{0, 0}, Rational(2)}, /**< 单项式 2。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(1, 2);
    MultiPoly sum_then_eval = (f + g).eval("x", a);
    MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "P2 additive: eval(f+g, x=1/2) == eval(f,x=1/2)+eval(g,x=1/2) [bivariate 3]";
}

TEST(MultipolyEvaluation, AdditiveHomomorphismTrivariateZ) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 1, 1}, Rational(1)}, /**< 单项式 xyz。 */
        {{0, 0, 2}, Rational(1)}, /**< 单项式 z^2。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0, 1}, Rational(2)},  /**< 单项式 2xz。 */
        {{0, 1, 0}, Rational(-1)}, /**< 单项式 -y。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(3);
    MultiPoly sum_then_eval = (f + g).eval("z", a);
    MultiPoly eval_then_sum = f.eval("z", a) + g.eval("z", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "P2 additive: eval(f+g, z=3) == eval(f,z=3)+eval(g,z=3) [trivariate 4]";
}

TEST(MultipolyEvaluation, AdditiveHomomorphismTrivariateZero) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0, 0}, Rational(1)}, /**< 单项式 x^2。 */
        {{0, 2, 0}, Rational(1)}, /**< 单项式 y^2。 */
        {{0, 0, 2}, Rational(1)}, /**< 单项式 z^2。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{2, 0, 0}, Rational(-1)}, /**< 单项式 -x^2。 */
        {{0, 1, 1}, Rational(2)},  /**< 单项式 2yz。 */
        {{0, 0, 0}, Rational(7)},  /**< 单项式 7。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(0);
    MultiPoly sum_then_eval = (f + g).eval("x", a);
    MultiPoly eval_then_sum = f.eval("x", a) + g.eval("x", a);
    EXPECT_TRUE((sum_then_eval == eval_then_sum)) << "P2 additive: eval(f+g, x=0) == eval(f,x=0)+eval(g,x=0) [trivariate 5]";
}

TEST(MultipolyEvaluation, MultiplicativeHomomorphismBivariateX) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 0}, Rational(1)}, /**< 单项式 x。 */
        {{0, 1}, Rational(1)}, /**< 单项式 y。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)},  /**< 单项式 x。 */
        {{0, 1}, Rational(-1)}, /**< 单项式 -y。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(3);
    MultiPoly prod_then_eval = (f * g).eval("x", a);
    MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "P2 multiplicative: eval(f*g, x=3) == eval(f,x=3)*eval(g,x=3) [bivariate 1]";
}

TEST(MultipolyEvaluation, MultiplicativeHomomorphismBivariateY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{2, 0}, Rational(2)}, /**< 单项式 2x^2。 */
        {{0, 1}, Rational(1)}, /**< 单项式 y。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)}, /**< 单项式 x。 */
        {{0, 2}, Rational(3)}, /**< 单项式 3y^2。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(2);
    MultiPoly prod_then_eval = (f * g).eval("y", a);
    MultiPoly eval_then_prod = f.eval("y", a) * g.eval("y", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "P2 multiplicative: eval(f*g, y=2) == eval(f,y=2)*eval(g,y=2) [bivariate 2]";
}

TEST(MultipolyEvaluation, MultiplicativeHomomorphismBivariateRational) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 1}, Rational(1)}, /**< 单项式 xy。 */
        {{0, 0}, Rational(1)}, /**< 单项式 1。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 0}, Rational(1)},  /**< 单项式 x。 */
        {{0, 1}, Rational(-1)}, /**< 单项式 -y。 */
        {{0, 0}, Rational(2)},  /**< 单项式 2。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(2, 3);
    MultiPoly prod_then_eval = (f * g).eval("x", a);
    MultiPoly eval_then_prod = f.eval("x", a) * g.eval("x", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "P2 multiplicative: eval(f*g, x=2/3) == eval(f,x=2/3)*eval(g,x=2/3) [bivariate 3]";
}

TEST(MultipolyEvaluation, MultiplicativeHomomorphismTrivariateYNegative) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 0, 0}, Rational(1)}, /**< 单项式 x。 */
        {{0, 1, 0}, Rational(1)}, /**< 单项式 y。 */
        {{0, 0, 1}, Rational(1)}, /**< 单项式 z。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{1, 1, 0}, Rational(1)},  /**< 单项式 xy。 */
        {{0, 0, 1}, Rational(-1)}, /**< 单项式 -z。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(-2);
    MultiPoly prod_then_eval = (f * g).eval("y", a);
    MultiPoly eval_then_prod = f.eval("y", a) * g.eval("y", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "P2 multiplicative: eval(f*g, y=-2) == eval(f,y=-2)*eval(g,y=-2) [trivariate 4]";
}

TEST(MultipolyEvaluation, MultiplicativeHomomorphismTrivariateZero) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> f_terms = {
        {{1, 0, 1}, Rational(1)}, /**< 单项式 xz。 */
        {{0, 2, 0}, Rational(1)}, /**< 单项式 y^2。 */
        {{0, 0, 0}, Rational(1)}, /**< 单项式 1。 */
    };
    std::vector<MultiPoly::Term> g_terms = {
        {{0, 0, 1}, Rational(3)},  /**< 单项式 3z。 */
        {{1, 0, 0}, Rational(-1)}, /**< 单项式 -x。 */
    };
    MultiPoly f(f_terms, vars);
    MultiPoly g(g_terms, vars);

    Rational a(0);
    MultiPoly prod_then_eval = (f * g).eval("z", a);
    MultiPoly eval_then_prod = f.eval("z", a) * g.eval("z", a);
    EXPECT_TRUE((prod_then_eval == eval_then_prod)) << "P2 multiplicative: eval(f*g, z=0) == eval(f,z=0)*eval(g,z=0) [trivariate 5]";
}

} // namespace
