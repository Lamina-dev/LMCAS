#include "test_multivariate_factor_support.hpp"
#include <algorithm>

TEST(MultivariateFactorComplete, CompleteFactorizationX2Y2XYXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), /**< 单项式 x^2。 */
        make_term({0, 2}, Rational(-1)) /**< 单项式 -y^2。 */
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "x^2-y^2: product equals original";
    EXPECT_TRUE((result.factors.size() == 2)) << "x^2-y^2: has 2 factors";
    EXPECT_TRUE((result.constant == Rational(1))) << "x^2-y^2: constant is 1";
}

TEST(MultivariateFactorComplete, CompleteFactorizationX2YXY2XyXY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), /**< 单项式 x^2*y。 */
        make_term({1, 2}, Rational(1))  /**< 单项式 x*y^2。 */
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "x^2*y+x*y^2: product equals original";
    EXPECT_TRUE((result.factors.size() >= 2)) << "x^2*y+x*y^2: has at least 2 factors";
}

TEST(MultivariateFactorComplete, CompleteFactorizationX22xyY2Z2XYZXYZ) {
    std::vector<std::string> vars = {"x", "y", "z"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0, 0}, Rational(1)), /**< 单项式 x^2。 */
        make_term({1, 1, 0}, Rational(2)), /**< 单项式 2xy。 */
        make_term({0, 2, 0}, Rational(1)), /**< 单项式 y^2。 */
        make_term({0, 0, 2}, Rational(-1)) /**< 单项式 -z^2。 */
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "x^2+2xy+y^2-z^2: product equals original";
    EXPECT_TRUE((result.factors.size() >= 1)) << "x^2+2xy+y^2-z^2: has at least 1 factor";
    EXPECT_TRUE((result.constant == Rational(1))) << "x^2+2xy+y^2-z^2: constant is 1";
}

TEST(MultivariateFactorComplete, CompleteFactorization6x2Y3xY23xy2xY) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(6)), /**< 单项式 6x^2*y。 */
        make_term({1, 2}, Rational(-3)) /**< 单项式 -3x*y^2。 */
    };
    MultiPoly poly(terms, vars);

    MultiFactorResult result = checked_factor_multivariate(poly);

    MultiPoly product(Rational(result.constant), vars);
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int m = 0; m < result.multiplicities[i]; ++m) {
            product = product * result.factors[i];
        }
    }
    EXPECT_TRUE((product == poly)) << "6x^2*y-3x*y^2: product equals original";
    EXPECT_TRUE((result.factors.size() >= 2)) << "6x^2*y-3x*y^2: has at least 2 factors";
}

TEST(MultivariateFactorComplete, CompleteFactorizationConstant42Constant42Factors) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(42), vars);

    MultiFactorResult result = checked_factor_multivariate(poly);
    EXPECT_TRUE((result.factors.empty())) << "constant 42: no factors";
    EXPECT_TRUE((result.constant == Rational(42))) << "constant 42: constant is 42";
}

TEST(MultivariateFactorComplete, CompleteFactorizationZeroConstant0Factors) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly(Rational(0), vars);

    MultiFactorResult result = checked_factor_multivariate(poly);
    EXPECT_TRUE((result.factors.empty())) << "zero: no factors";
    EXPECT_TRUE((result.constant == Rational(0))) << "zero: constant is 0";
}

TEST(MultivariateFactorComplete, CheckedFactorizationX2Y21) {
    const std::vector<std::string> vars = {"x", "y"};
    MultiPoly poly({make_term({2, 0}, Rational(1)),
                    make_term({0, 2}, Rational(1)),
                    make_term({0, 0}, Rational(1))}, vars);

    auto checked = factor_multivariate_checked(poly);
    ASSERT_TRUE((checked.has_value()));
    EXPECT_TRUE((checked.value().completeness == Completeness::Inconclusive));
    ASSERT_TRUE((checked.value().value.factors.size() == 1));
    EXPECT_TRUE((checked.value().value.factors[0] == poly));
    EXPECT_TRUE((checked.value().value.multiplicities[0] == 1));
    EXPECT_TRUE((checked.value().value.constant == Rational(1)));
}

TEST(MultivariateFactorComplete, CheckedEmbeddedUnivariateFactor) {
    const std::vector<std::string> vars = {"x", "y"};
    const MultiPoly poly({make_term({2, 0}, Rational(1)),
                          make_term({0, 0}, Rational(1))}, vars);
    auto checked = factor_multivariate_checked(poly);
    ASSERT_TRUE((checked.has_value()));
    EXPECT_TRUE((checked.value().completeness == Completeness::Complete));
    ASSERT_TRUE((checked.value().value.factors.size() == 1));
    EXPECT_TRUE((checked.value().value.factors[0] == poly));
}

TEST(MultivariateFactorComplete, InconclusiveFactorizationPreservesOriginal) {
    const std::vector<std::string> vars = {"x", "y"};
    const MultiPoly residual({make_term({2, 0}, Rational(1)),
                              make_term({0, 2}, Rational(1)),
                              make_term({0, 0}, Rational(1))}, vars);
    const MultiPoly x({make_term({1, 0}, Rational(1))}, vars);
    const MultiPoly poly = x * residual;

    auto checked = factor_multivariate_checked(poly);
    ASSERT_TRUE(checked.has_value());
    EXPECT_EQ(checked.value().completeness, Completeness::Inconclusive);
    const auto& result = checked.value().value;
    ASSERT_EQ(result.factors.size(), 1u);
    EXPECT_TRUE(result.factors[0] * result.constant == poly);
    EXPECT_EQ(result.multiplicities[0], 1);
}

TEST(MultivariateFactorComplete, CheckedFactorization) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly difference({make_term({2, 0}, Rational(1)),
                          make_term({0, 2}, Rational(-1))},
                         vars);
    auto complete = factor_multivariate_checked(difference);
    ASSERT_TRUE((complete.has_value())) << "平方差分解应成功";
    EXPECT_TRUE((complete.value().completeness == Completeness::Complete));
    const auto& factors = complete.value().value;
    ASSERT_TRUE((factors.factors.size() == 2));
    MultiPoly product(factors.constant, vars);
    for (size_t i = 0; i < factors.factors.size(); ++i) {
        for (int m = 0; m < factors.multiplicities[i]; ++m) {
            product = product * factors.factors[i];
        }
    }
    EXPECT_TRUE((product == difference));

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);
    auto limited =
        factor_multivariate_checked(difference, context);
    EXPECT_FALSE((limited.has_value())) << "耗尽预算应返回错误";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "预算错误应报告 ResourceLimit";

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled =
        factor_multivariate_checked(difference, cancelled_context);
    EXPECT_FALSE((cancelled.has_value())) << "取消的受检分解不得降级为未决成功";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "取消应通过 CasError 显式传播";
}
