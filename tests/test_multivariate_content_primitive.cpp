#include "test_multivariate_support.hpp"

namespace {

MultiPoly embed_content(
    const MultiPoly &content,
    const std::string &main_variable,
    const std::vector<std::string> &full_variables) {
    std::size_t main_index = full_variables.size();
    for (std::size_t index = 0;
         index < full_variables.size(); ++index) {
        if (full_variables[index] == main_variable) {
            main_index = index;
            break;
        }
    }
    if (main_index == full_variables.size()) {
        return content;
    }

    std::vector<MultiPoly::Term> embedded;
    embedded.reserve(content.terms().size());
    for (const auto &[monomial, coefficient] : content.terms()) {
        Monomial full_monomial(full_variables.size(), 0);
        std::size_t source_index = 0;
        for (std::size_t index = 0;
             index < full_variables.size(); ++index) {
            if (index != main_index) {
                if (source_index < monomial.size()) {
                    full_monomial[index] = monomial[source_index];
                }
                ++source_index;
            }
        }
        embedded.emplace_back(
            std::move(full_monomial), coefficient);
    }
    return MultiPoly(std::move(embedded), full_variables);
}

} // namespace

TEST(MultivariateContentPrimitive, TrivialContentX2X1Content1) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), /**< 项 x^2。 */
        make_term({1, 0}, Rational(1)), /**< 项 x。 */
        make_term({0, 0}, Rational(1))  /**< 常数项 1。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for x^2+x+1";
    EXPECT_TRUE((content.is_constant())) << "content of x^2+x+1 is trivial (constant)";
}

TEST(MultivariateContentPrimitive, NonTrivialContentX2YXYContentY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), /**< 项 x^2*y。 */
        make_term({1, 1}, Rational(1))  /**< 项 x*y。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for x^2*y + x*y";

    std::vector<std::string> rem_vars = {"y"};
    std::vector<MultiPoly::Term> coeff_x2 = {make_term({1}, Rational(1))}; /**< 系数 y。 */
    std::vector<MultiPoly::Term> coeff_x1 = {make_term({1}, Rational(1))}; /**< 系数 y。 */
    MultiPoly c_x2(coeff_x2, rem_vars);
    MultiPoly c_x1(coeff_x1, rem_vars);
    bool divides_x2 = true, divides_x1 = true;
    try {
        c_x2.exact_div(content);
    } catch (...) {
        divides_x2 = false;
    }
    try {
        c_x1.exact_div(content);
    } catch (...) {
        divides_x1 = false;
    }
    EXPECT_TRUE((divides_x2)) << "content divides coeff(x^2) for x^2*y + x*y";
    EXPECT_TRUE((divides_x1)) << "content divides coeff(x^1) for x^2*y + x*y";
}

TEST(MultivariateContentPrimitive, NumericContent6X24XContent2) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(6)), /**< 项 6*x^2。 */
        make_term({1, 0}, Rational(4))  /**< 项 4*x。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for 6x^2 + 4x";
    EXPECT_TRUE((content.is_constant())) << "content of 6x^2+4x is constant";
}

TEST(MultivariateContentPrimitive, TrivariatePolynomialX2YZXYZ2ContentYZ) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1, 1}, Rational(1)), /**< 项 x^2*y*z。 */
        make_term({1, 1, 2}, Rational(1))  /**< 项 x*y*z^2。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for trivariate poly";

    std::vector<std::string> rem_vars = {"y", "z"};
    std::vector<MultiPoly::Term> coeff_x2_terms = {make_term({1, 1}, Rational(1))}; /**< 系数 y*z。 */
    std::vector<MultiPoly::Term> coeff_x1_terms = {make_term({1, 2}, Rational(1))}; /**< 系数 y*z^2。 */
    MultiPoly c_x2(coeff_x2_terms, rem_vars);
    MultiPoly c_x1(coeff_x1_terms, rem_vars);
    bool div2 = true, div1 = true;
    try {
        c_x2.exact_div(content);
    } catch (...) {
        div2 = false;
    }
    try {
        c_x1.exact_div(content);
    } catch (...) {
        div1 = false;
    }
    EXPECT_TRUE((div2)) << "content divides coeff(x^2) for trivariate poly";
    EXPECT_TRUE((div1)) << "content divides coeff(x^1) for trivariate poly";
}

TEST(MultivariateContentPrimitive, MixedNumericAndPolynomialContent2X2Y4XY2Content2y) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(2)), /**< 项 2*x^2*y。 */
        make_term({1, 2}, Rational(4))  /**< 项 4*x*y^2。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for 2x^2*y + 4x*y^2";

    std::vector<std::string> rem_vars = {"y"};
    std::vector<MultiPoly::Term> coeff_x2_terms = {make_term({1}, Rational(2))}; /**< 系数 2y。 */
    std::vector<MultiPoly::Term> coeff_x1_terms = {make_term({2}, Rational(4))}; /**< 系数 4y^2。 */
    MultiPoly c_x2(coeff_x2_terms, rem_vars);
    MultiPoly c_x1(coeff_x1_terms, rem_vars);
    bool div2 = true, div1 = true;
    try {
        c_x2.exact_div(content);
    } catch (...) {
        div2 = false;
    }
    try {
        c_x1.exact_div(content);
    } catch (...) {
        div1 = false;
    }
    EXPECT_TRUE((div2)) << "content divides coeff(x^2) for 2x^2*y + 4x*y^2";
    EXPECT_TRUE((div1)) << "content divides coeff(x^1) for 2x^2*y + 4x*y^2";
}

TEST(MultivariateContentPrimitive, SingleTermPolynomialX3Y2ZContentY2Z) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        make_term({3, 2, 1}, Rational(1)) /**< 项 x^3*y^2*z。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for single term x^3*y^2*z";
}

TEST(MultivariateContentPrimitive, ConstantPolynomialContentPolyItself) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly f(Rational(7), vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for constant 7";
}

TEST(MultivariateContentPrimitive, ThreeTermTrivariateX3YX2YZXYZ2ContentY) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        make_term({3, 1, 0}, Rational(1)), /**< 项 x^3*y。 */
        make_term({2, 1, 1}, Rational(1)), /**< 项 x^2*y*z。 */
        make_term({1, 1, 2}, Rational(1))  /**< 项 x*y*z^2。 */
    };
    MultiPoly f(terms, vars);
    MultiPoly content = multivariate_content(f, "x");
    MultiPoly pp = multivariate_primitive_part(f, "x");
    MultiPoly content_full = embed_content(content, "x", vars);
    MultiPoly product = content_full * pp;
    EXPECT_TRUE((product == f)) << "content * primitive_part == f for x^3*y + x^2*y*z + x*y*z^2";

    std::vector<std::string> rem_vars = {"y", "z"};
    std::vector<MultiPoly::Term> c3_terms = {make_term({1, 0}, Rational(1))}; /**< 系数 y。 */
    std::vector<MultiPoly::Term> c2_terms = {make_term({1, 1}, Rational(1))}; /**< 系数 y*z。 */
    std::vector<MultiPoly::Term> c1_terms = {make_term({1, 2}, Rational(1))}; /**< 系数 y*z^2。 */
    MultiPoly c3(c3_terms, rem_vars);
    MultiPoly c2(c2_terms, rem_vars);
    MultiPoly c1(c1_terms, rem_vars);
    bool d3 = true, d2 = true, d1 = true;
    try {
        c3.exact_div(content);
    } catch (...) {
        d3 = false;
    }
    try {
        c2.exact_div(content);
    } catch (...) {
        d2 = false;
    }
    try {
        c1.exact_div(content);
    } catch (...) {
        d1 = false;
    }
    EXPECT_TRUE((d3)) << "content divides coeff(x^3) for trivariate 3-term";
    EXPECT_TRUE((d2)) << "content divides coeff(x^2) for trivariate 3-term";
    EXPECT_TRUE((d1)) << "content divides coeff(x^1) for trivariate 3-term";
}
