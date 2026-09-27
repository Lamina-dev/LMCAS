#include "limit_result.hpp"
/**
 * @file test_taylor_fallback.cpp
 * @brief Taylor 展开回退策略测试.
 *
 * 验证 L'Hôpital 规则达到有限迭代边界后,
 * Taylor 级数展开继续完成极限计算.
 *
 * 覆盖需求: 3.1, 3.2, 3.3, 3.4
 */
#include "test_common.hpp"

using namespace LMCAS;

TEST(TaylorFallback, CubicSineRemainder) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto three = SymbolicExpr::number(3);
    auto neg_one = SymbolicExpr::number(-1);

    auto sin_x = SymbolicExpr::sin(x);
    auto neg_x = SymbolicExpr::multiply(x, neg_one);
    auto num = SymbolicExpr::add(sin_x, neg_x);
    auto den = SymbolicExpr::power(x, three);
    auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE((lim != nullptr)) << "limit((sin(x)-x)/x^3, x->0) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            {
                const double actual_value = (*val);
                const double expected_value = (-1.0 / 6.0);
                const double tolerance = (1e-6);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        } else {
            ADD_FAILURE() << "limit((sin(x)-x)/x^3, x->0) should be numeric -1/6; actual: "
                          << lim->to_string();
        }
    }
}

TEST(TaylorFallback, QuadraticCosineRemainder) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto neg_one = SymbolicExpr::number(-1);

    auto cos_x = SymbolicExpr::cos(x);
    auto neg_cos = SymbolicExpr::multiply(cos_x, neg_one);
    auto num = SymbolicExpr::add(one, neg_cos);
    auto den = SymbolicExpr::power(x, two);
    auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE((lim != nullptr)) << "limit((1-cos(x))/x^2, x->0) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            EXPECT_NEAR(*val, 0.5, 1e-6) << "limit((1-cos(x))/x^2, x->0) = 1/2";
        } else {
            ADD_FAILURE() << "limit((1-cos(x))/x^2, x->0) should be numeric 1/2; actual: "
                          << lim->to_string();
        }
    }
}

TEST(TaylorFallback, QuadraticExponentialRemainder) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto neg_one = SymbolicExpr::number(-1);

    auto exp_x = SymbolicExpr::exp(x);
    auto neg_1 = SymbolicExpr::multiply(one, neg_one);
    auto neg_x = SymbolicExpr::multiply(x, neg_one);
    auto num = SymbolicExpr::add(SymbolicExpr::add(exp_x, neg_1), neg_x);
    auto den = SymbolicExpr::power(x, two);
    auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE((lim != nullptr)) << "limit((e^x-1-x)/x^2, x->0) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            EXPECT_NEAR(*val, 0.5, 1e-6) << "limit((e^x-1-x)/x^2, x->0) = 1/2";
        } else {
            ADD_FAILURE() << "limit((e^x-1-x)/x^2, x->0) should be numeric 1/2; actual: "
                          << lim->to_string();
        }
    }
}

TEST(TaylorFallback, SineLinearTerm) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::multiply(sin_x, SymbolicExpr::power(x, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE((lim != nullptr)) << "limit(sin(x)/x, x->0) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "limit(sin(x)/x, x->0) = 1";
        } else {
            EXPECT_TRUE(test_expression_text((lim), ("1"))) << "limit(sin(x)/x, x->0) = 1";
        }
    }
}

TEST(TaylorFallback, TangentLinearTerm) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    auto tan_x = SymbolicExpr::tan(x);
    auto expr = SymbolicExpr::multiply(tan_x, SymbolicExpr::power(x, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE((lim != nullptr)) << "limit(tan(x)/x, x->0) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            EXPECT_NEAR(*val, 1.0, 1e-6) << "limit(tan(x)/x, x->0) = 1";
        } else {
            EXPECT_TRUE(test_expression_text((lim), ("1"))) << "limit(tan(x)/x, x->0) = 1";
        }
    }
}

TEST(TaylorFallback, ReciprocalSineAtInfinity) {
    auto x = SymbolicExpr::variable("x");
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto expr = SymbolicExpr::multiply(sin_inv_x, SymbolicExpr::power(inv_x, neg_one));
    auto result = LMCAS::limit_expression_checked(expr, "x", inf);
    ASSERT_TRUE(result.has_value()) << result.error().message;

    auto lim = result.value();
    ASSERT_TRUE(lim) << "limit(sin(1/x)/(1/x), x->inf) is not null";
    auto value = test_numeric_eval(lim);
    ASSERT_TRUE(value.has_value()) << "limit should be numeric; actual: " << lim->to_string();
    EXPECT_NEAR(*value, 1.0, 1e-6);
}

TEST(TaylorFallback, RationalLeadingTerms) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    auto x_sq = SymbolicExpr::power(x, two);
    auto num = SymbolicExpr::add(x_sq, x);
    auto den = SymbolicExpr::add(
        SymbolicExpr::multiply(two, x_sq),
        three);
    auto expr = SymbolicExpr::multiply(num, SymbolicExpr::power(den, neg_one));
    auto lim = LMCAS::limit_expression_checked(expr, "x", inf).value();
    EXPECT_TRUE((lim != nullptr)) << "limit((x^2+x)/(2x^2+3), x->inf) is not null";
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            EXPECT_NEAR(*val, 0.5, 1e-6) << "limit((x^2+x)/(2x^2+3), x->inf) = 1/2";
        } else {
            ADD_FAILURE() << "limit((x^2+x)/(2x^2+3), x->inf) should be numeric 1/2; actual: "
                          << lim->to_string();
        }
    }
}
