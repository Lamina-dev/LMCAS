#include "test_common.hpp"
#include "symbolic_ode_engine.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include <limits>

using namespace LMCAS;

TEST(OdeEngineClassification, BernoulliExponents) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto linear_term = SymbolicExpr::multiply(
        SymbolicExpr::number(-2), y);
    std::shared_ptr<SymbolicExpr> detected_p;
    std::shared_ptr<SymbolicExpr> detected_q;
    int detected_n = 0;

    auto quartic_rhs = SymbolicExpr::add(
        linear_term,
        SymbolicExpr::multiply(
            x, SymbolicExpr::power(y, SymbolicExpr::number(4))));
    EXPECT_TRUE((is_bernoulli_ode(
        quartic_rhs, "x", "y", detected_p, detected_q, detected_n)))
        << "quartic Bernoulli equation is recognized";
    EXPECT_TRUE((detected_n == 4)) << "quartic Bernoulli exponent is preserved";
    EXPECT_TRUE(test_expression_text(detected_p, SymbolicExpr::number(2))) << "quartic Bernoulli P is extracted";
    EXPECT_TRUE(test_expression_text(detected_q, x)) << "quartic Bernoulli Q is extracted";

    auto inverse_rhs = SymbolicExpr::add(
        linear_term,
        SymbolicExpr::multiply(
            x, SymbolicExpr::power(y, SymbolicExpr::number(-1))));
    EXPECT_TRUE((is_bernoulli_ode(
        inverse_rhs, "x", "y", detected_p, detected_q, detected_n)))
        << "negative-exponent Bernoulli equation is recognized";
    EXPECT_TRUE((detected_n == -1)) << "negative Bernoulli exponent is preserved";
    EXPECT_TRUE(test_expression_text(detected_p, SymbolicExpr::number(2))) << "negative-exponent Bernoulli P is extracted";
    EXPECT_TRUE(test_expression_text(detected_q, x)) << "negative-exponent Bernoulli Q is extracted";
}

TEST(OdeEngineClassification, ClassifyRatio) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(y, x);

    auto cls = classify_first_order_ode(rhs, "x", "y");
    EXPECT_TRUE((cls.type == ODEType::Separable ||
                 cls.type == ODEType::Linear1 ||
                 cls.type == ODEType::Homogeneous))
        << "y/x classified as separable, linear, or homogeneous";
}

TEST(OdeEngineClassification, ClassifyNonlinearOrder) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto cls = classify_first_order_ode(rhs, "x", "y");
    EXPECT_TRUE((cls.order == 1)) << "classified as first order";
}

TEST(OdeEngineClassification, ClassifyUnrepresentable) {
    std::shared_ptr<SymbolicExpr> null_rhs;
    auto first_order = classify_first_order_ode(null_rhs, "x", "y");
    EXPECT_TRUE((first_order.type == ODEType::Unknown)) << "null first-order RHS is not classified as separable";

    auto parameter = SymbolicExpr::variable("a");
    std::vector<std::shared_ptr<SymbolicExpr>> coeffs{
        parameter, SymbolicExpr::number(1), SymbolicExpr::number(1)};
    auto higher_order =
        classify_higher_order_ode(coeffs, SymbolicExpr::number(0), "x", "y");
    EXPECT_TRUE((higher_order.type == ODEType::Unknown)) << "symbolic constants are not represented as numeric zero coefficients";
    EXPECT_TRUE((higher_order.const_coeffs.empty())) << "unsupported constant coefficients do not expose fabricated values";
}

TEST(OdeEngineClassification, ClassificationProofs) {
    std::shared_ptr<SymbolicExpr> P;
    std::shared_ptr<SymbolicExpr> Q;
    std::shared_ptr<SymbolicExpr> null_rhs;
    EXPECT_FALSE((is_linear_first_order(null_rhs, "x", "y", P, Q))) << "null RHS is not a linear ODE";
    EXPECT_FALSE((is_separable(null_rhs, "x", "y"))) << "null RHS is not a separable ODE";

    auto y = SymbolicExpr::variable("y");
    auto nonlinear = SymbolicExpr::power(
        y, SymbolicExpr::number(2));
    P = SymbolicExpr::number(9);
    Q = SymbolicExpr::number(9);
    EXPECT_FALSE((is_linear_first_order(nonlinear, "x", "y", P, Q))) << "nonlinear RHS is not a linear ODE";
    EXPECT_TRUE((!P && !Q)) << "failed linear classification clears extracted coefficients";

    auto one = SymbolicExpr::number(1);
    EXPECT_FALSE((is_separable(one, "x", "x"))) << "separable classifier rejects duplicate variable names";
    EXPECT_FALSE((is_homogeneous_ode(one, "", "y"))) << "homogeneous classifier rejects an empty variable name";
    int bernoulli_n = 7;
    EXPECT_FALSE((is_bernoulli_ode(
        nonlinear, "y", "y", P, Q, bernoulli_n)))
        << "Bernoulli classifier rejects duplicate variable names";
    EXPECT_TRUE((!P && !Q && bernoulli_n == 0)) << "failed Bernoulli classification clears extracted outputs";
    EXPECT_FALSE((is_exact_ode(one, one, "x", "x"))) << "exact classifier rejects duplicate variable names";
    EXPECT_FALSE((is_constant_coefficient({}, "x"))) << "an empty coefficient list is not a constant-coefficient ODE";

    auto x = SymbolicExpr::variable("x");
    auto rhs = SymbolicExpr::multiply(
        SymbolicExpr::multiply(
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(2), x),
                SymbolicExpr::number(-1)),
            SymbolicExpr::add(x, SymbolicExpr::number(-1))),
        SymbolicExpr::multiply(
            SymbolicExpr::add(x, SymbolicExpr::number(-2)),
            SymbolicExpr::add(x, SymbolicExpr::number(-4))));
    EXPECT_FALSE((is_homogeneous_ode(rhs, "x", "y"))) << "zeros at fixed sample points do not prove homogeneity";
}

TEST(OdeEngineClassification, EulerClassification) {
    auto x = SymbolicExpr::variable("x");
    auto ratio = SymbolicExpr::add(
        SymbolicExpr::multiply(
            SymbolicExpr::add(x, SymbolicExpr::number(-1)),
            SymbolicExpr::add(x, SymbolicExpr::number(-2))),
        SymbolicExpr::number(1));
    auto sampled_coefficient =
        SymbolicExpr::multiply(x, ratio)->simplify();
    std::vector<double> constants{42.0};
    EXPECT_FALSE((is_euler_equation(
        {sampled_coefficient, SymbolicExpr::number(1)},
        "x", constants)))
        << "matching at two sample points does not prove an Euler coefficient";
    EXPECT_TRUE((constants.empty())) << "failed Euler classification clears extracted constants";
    EXPECT_FALSE((is_euler_equation({}, "x", constants))) << "an empty coefficient list is not an Euler equation";
    EXPECT_FALSE((is_euler_equation(
        {SymbolicExpr::number(0), SymbolicExpr::number(1)},
        "x", constants)))
        << "a zero leading coefficient does not define the stated order";
}
