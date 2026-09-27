#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

TEST(Cancel, CancelX21X1X1) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);

    // (x^2 - 1) / (x - 1)
    auto x2 = SymbolicExpr::power(x, two);
    auto numerator = SymbolicExpr::add(x2, SymbolicExpr::number(-1));
    auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 验证：代入 x=5 应得 6
    auto val = result->substitute("x", SymbolicExpr::number(5))->simplify();
    EXPECT_TRUE(test_expression_text((val), ("6"))) << "(x^2-1)/(x-1) at x=5 = 6";

    // 验证结果不含分母（无负指数）
    auto result_str = result->to_string();
    bool no_negative_power = result_str.find("^(-1)") == std::string::npos;
    EXPECT_TRUE((no_negative_power)) << "(x^2-1)/(x-1) cancels to polynomial";
}

TEST(Cancel, CancelX22x1X1X1) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);

    // (x+1)^2 / (x+1) = x+1
    auto x2 = SymbolicExpr::power(x, two);
    auto numerator = SymbolicExpr::add(x2, SymbolicExpr::add(
                                               SymbolicExpr::multiply(two, x), one));
    auto denominator = SymbolicExpr::add(x, one);
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=3 应得 4
    auto val = result->substitute("x", SymbolicExpr::number(3))->simplify();
    EXPECT_TRUE(test_expression_text((val), ("4"))) << "(x^2+2x+1)/(x+1) at x=3 = 4";
}

TEST(Cancel, CancelX3XX21X) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    // x(x^2-1) / (x^2-1) = x
    auto x3 = SymbolicExpr::power(x, three);
    auto numerator = SymbolicExpr::add(x3, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    auto denominator = SymbolicExpr::add(SymbolicExpr::power(x, two), SymbolicExpr::number(-1));
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=7 应得 7
    auto val = result->substitute("x", SymbolicExpr::number(7))->simplify();
    EXPECT_TRUE(test_expression_text((val), ("7"))) << "(x^3-x)/(x^2-1) at x=7 = 7";
}

TEST(Cancel, CancelWithNoCommonFactorX1X2) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);

    auto numerator = SymbolicExpr::add(x, one);
    auto denominator = SymbolicExpr::add(x, two);
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=3 应得 4/5 = 0.8，验证数值正确性
    // 由于无公因式，结果仍为分式形式
    auto val_num = result->substitute("x", SymbolicExpr::number(3))->simplify();
    // (3+1)/(3+2) = 4/5
    auto expected = SymbolicExpr::divide(SymbolicExpr::number(4), SymbolicExpr::number(5))->simplify();
    auto diff = SymbolicExpr::add(val_num, SymbolicExpr::multiply(SymbolicExpr::number(-1), expected))->simplify();
    EXPECT_TRUE((diff->is_zero())) << "(x+1)/(x+2) at x=3 = 4/5";
}

TEST(Cancel, CancelConstant632) {
    auto three = SymbolicExpr::number(3);

    auto expr = SymbolicExpr::divide(SymbolicExpr::number(6), three);
    auto result = expr->cancel();
    EXPECT_TRUE(test_expression_text((result), ("2"))) << "6/3 = 2";
}

TEST(Cancel, CancelMultivariateXYYX1Y) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    // (xy - y) / (x - 1) = y(x-1)/(x-1) = y
    auto xy = SymbolicExpr::multiply(x, y);
    auto numerator = SymbolicExpr::add(xy, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
    auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=3, y=5 应得 5
    auto val = result->substitute("x", SymbolicExpr::number(3))
                   ->substitute("y", SymbolicExpr::number(5))
                   ->simplify();
    EXPECT_TRUE(test_expression_text((val), ("5"))) << "(xy-y)/(x-1) at x=3,y=5 = 5";
}

TEST(Cancel, CancelPolynomialPolynomialAlreadySimplified) {
    auto x = SymbolicExpr::variable("x");

    // x / 1 = x (no denominator)
    auto result = x->cancel();
    EXPECT_TRUE(test_expression_text((result), ("x"))) << "x cancels to x";
}

TEST(Cancel, Cancel2x22x2xX1) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);

    auto x2 = SymbolicExpr::power(x, two);
    auto numerator = SymbolicExpr::add(
        SymbolicExpr::multiply(two, x2),
        SymbolicExpr::multiply(two, x));
    auto denominator = SymbolicExpr::multiply(two, x);
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=4 应得 5
    auto val = result->substitute("x", SymbolicExpr::number(4))->simplify();
    EXPECT_TRUE(test_expression_text((val), ("5"))) << "(2x^2+2x)/(2x) at x=4 = 5";
}

TEST(Cancel, FactorCubicX36x211x6X1X2X3) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    // x^3 - 6x^2 + 11x - 6
    auto x3 = SymbolicExpr::power(x, three);
    auto x2 = SymbolicExpr::power(x, two);
    auto expr = SymbolicExpr::add(x3,
                                  SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-6), x2),
                                                    SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(11), x),
                                                                      SymbolicExpr::number(-6))));

    auto factored = expr->factor_checked().value();

    // 验证：代入 x=1, x=2, x=3 应得 0
    auto v1 = factored->substitute("x", SymbolicExpr::number(1))->simplify();
    auto v2 = factored->substitute("x", SymbolicExpr::number(2))->simplify();
    auto v3 = factored->substitute("x", SymbolicExpr::number(3))->simplify();
    EXPECT_TRUE((v1->is_zero())) << "x^3-6x^2+11x-6 at x=1 = 0";
    EXPECT_TRUE((v2->is_zero())) << "x^3-6x^2+11x-6 at x=2 = 0";
    EXPECT_TRUE((v3->is_zero())) << "x^3-6x^2+11x-6 at x=3 = 0";

    // 验证是乘积形式
    auto s = factored->to_string();
    bool is_product = s.find("*") != std::string::npos;
    EXPECT_TRUE((is_product)) << "cubic factored into product form";
}

TEST(Cancel, CancelHigherDegreeX3XX2XX1) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);

    // (x^3 - x) / (x^2 + x) = x(x^2-1) / x(x+1) = (x-1)(x+1)/(x+1) = x-1
    auto x3 = SymbolicExpr::power(x, three);
    auto x2 = SymbolicExpr::power(x, two);
    auto numerator = SymbolicExpr::add(x3, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    auto denominator = SymbolicExpr::add(x2, x);
    auto expr = SymbolicExpr::divide(numerator, denominator);

    auto result = expr->cancel();

    // 代入 x=5 应得 4
    auto val = result->substitute("x", SymbolicExpr::number(5))->simplify();
    EXPECT_TRUE(test_expression_text((val), ("4"))) << "(x^3-x)/(x^2+x) at x=5 = 4";
}

TEST(Cancel, CancelPreservesUnsupportedNumerators) {
    auto x = SymbolicExpr::variable("x");
    auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto radical = SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 2)));
    auto sine = SymbolicExpr::sin(x);
    auto radical_fraction = SymbolicExpr::divide(radical, denominator);
    auto sine_fraction = SymbolicExpr::divide(sine, denominator);
    auto sum = SymbolicExpr::add(radical_fraction, sine_fraction);
    auto nested = SymbolicExpr::divide(
        SymbolicExpr::multiply(radical,
                               SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(2))),
        denominator);
    const std::vector<std::pair<ExprPtr, double>> cases{
        {radical_fraction, 0.5}, {sine_fraction, std::sin(1.0) / 2.0}, {sum, (1.0 + std::sin(1.0)) / 2.0}, {nested, 1.5}};
    for (const auto &sample : cases) {
        auto cancelled = sample.first->cancel();
        ASSERT_TRUE(cancelled != nullptr) << "cancel returns an expression";
        auto evaluated = evaluate_numeric(*cancelled, {{"x", 1.0}});
        ASSERT_TRUE(evaluated.has_value())
            << "nonpolynomial fraction evaluates at a regular point";
        const double actual_value = evaluated.value().value;
        const double expected_value = sample.second;
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, 1e-12);
        auto at_pole = evaluate_numeric(*cancelled, {{"x", -1.0}});
        EXPECT_TRUE((!at_pole && at_pole.error().code == CasErrc::DomainError)) << "cancel does not turn an original singularity into zero";
        ComputationContext context;
        auto residual = check_zero_residual(sample.first, context);
        EXPECT_TRUE((residual && !std::holds_alternative<ProvedZeroResidual>(residual.value()))) << "failed polynomial conversion is not a zero residual proof";
    }
}

TEST(Cancel, CombinesTermsAlreadyUsingTheCommonDenominator) {
    auto x = SymbolicExpr::variable("x");
    auto reciprocal = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto expression = SymbolicExpr::add(
        SymbolicExpr::add(
            reciprocal,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x)),
        SymbolicExpr::divide(
            SymbolicExpr::power(x, SymbolicExpr::number(2)), x));

    auto cancelled = expression->cancel();

    ASSERT_TRUE(cancelled);
    EXPECT_TRUE(test_same_expression(cancelled, reciprocal->simplify()));
}
