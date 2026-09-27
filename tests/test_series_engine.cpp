/**
 * @file test_series_engine.cpp
 * @brief 级数引擎单元测试：收敛半径、收敛判定、幂级数运算、傅里叶级数、洛朗级数、符号求和与乘积。
 */

#include "test_common.hpp"
#include "series_engine.hpp"

using namespace LMCAS;

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }

TEST(SeriesEngine, PowerSeriesAddBasic) {
    // a(x) = 1 + 2x + 3x²
    Coeffs a = {num(1), num(2), num(3)};
    // b(x) = 4 + 5x + 6x²
    Coeffs b = {num(4), num(5), num(6)};

    auto result = LMCAS::power_series_add(a, b);

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(5))) << "c[0] = 1+4 = 5";
    EXPECT_TRUE(test_same_expression(result[1], num(7))) << "c[1] = 2+5 = 7";
    EXPECT_TRUE(test_same_expression(result[2], num(9))) << "c[2] = 3+6 = 9";
}

TEST(SeriesEngine, PowerSeriesAddDifferentLengths) {
    // a(x) = 1 + 2x
    Coeffs a = {num(1), num(2)};
    // b(x) = 3 + 4x + 5x² + 6x³
    Coeffs b = {num(3), num(4), num(5), num(6)};

    auto result = LMCAS::power_series_add(a, b);

    ASSERT_EQ(result.size(), 4u) << "Result has 4 coefficients (max length)";
    EXPECT_TRUE(test_same_expression(result[0], num(4))) << "c[0] = 1+3 = 4";
    EXPECT_TRUE(test_same_expression(result[1], num(6))) << "c[1] = 2+4 = 6";
    EXPECT_TRUE(test_same_expression(result[2], num(5))) << "c[2] = 0+5 = 5";
    EXPECT_TRUE(test_same_expression(result[3], num(6))) << "c[3] = 0+6 = 6";
}

TEST(SeriesEngine, PowerSeriesAddEmpty) {
    Coeffs a = {};
    Coeffs b = {num(1), num(2)};

    auto result = LMCAS::power_series_add(a, b);
    ASSERT_EQ(result.size(), 2u) << "Result has 2 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 0+1 = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(2))) << "c[1] = 0+2 = 2";
}

TEST(SeriesEngine, PowerSeriesMultiplyBasic) {
    // a(x) = 1 + x
    Coeffs a = {num(1), num(1)};
    // b(x) = 1 + x
    Coeffs b = {num(1), num(1)};

    auto checked_result = LMCAS::power_series_multiply_checked(a, b, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 1*1 = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(2))) << "c[1] = 1*1 + 1*1 = 2";
    EXPECT_TRUE(test_same_expression(result[2], num(1))) << "c[2] = 1*1 = 1";
}

TEST(SeriesEngine, PowerSeriesMultiplyTruncation) {
    // a(x) = 1 + x + x²
    Coeffs a = {num(1), num(1), num(1)};
    // b(x) = 1 + x + x²
    Coeffs b = {num(1), num(1), num(1)};

    // Full product: 1 + 2x + 3x² + 2x³ + x⁴
    // Truncate to order 3: 1 + 2x + 3x²
    auto checked_result = LMCAS::power_series_multiply_checked(a, b, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients (truncated)";
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(2))) << "c[1] = 2";
    EXPECT_TRUE(test_same_expression(result[2], num(3))) << "c[2] = 3";
}

TEST(SeriesEngine, PowerSeriesMultiplyZero) {
    Coeffs a = {num(1), num(2), num(3)};
    Coeffs b = {num(0)};

    auto checked_result = LMCAS::power_series_multiply_checked(a, b, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(0))) << "c[0] = 0";
    EXPECT_TRUE(test_same_expression(result[1], num(0))) << "c[1] = 0";
    EXPECT_TRUE(test_same_expression(result[2], num(0))) << "c[2] = 0";
}

TEST(SeriesEngine, PowerSeriesMultiplyConstants) {
    // a(x) = 3
    Coeffs a = {num(3)};
    // b(x) = 1 + 2x + 4x²
    Coeffs b = {num(1), num(2), num(4)};

    auto checked_result = LMCAS::power_series_multiply_checked(a, b, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(3))) << "c[0] = 3*1 = 3";
    EXPECT_TRUE(test_same_expression(result[1], num(6))) << "c[1] = 3*2 = 6";
    EXPECT_TRUE(test_same_expression(result[2], num(12))) << "c[2] = 3*4 = 12";
}

TEST(SeriesEngine, PowerSeriesComposeBasic) {
    // f(x) = 1 + x
    Coeffs f = {num(1), num(1)};
    // g(x) = x (i.e., g(0)=0, g[1]=1)
    Coeffs g = {num(0), num(1)};

    auto checked_result = LMCAS::power_series_compose_checked(f, g, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    // f(g(x)) = f(x) = 1 + x
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(1))) << "c[1] = 1";
    EXPECT_TRUE(test_same_expression(result[2], num(0))) << "c[2] = 0";
}

TEST(SeriesEngine, PowerSeriesComposeQuadratic) {
    // f(x) = 1 + x + x²
    Coeffs f = {num(1), num(1), num(1)};
    // g(x) = 2x (g(0)=0)
    Coeffs g = {num(0), num(2)};

    auto checked_result = LMCAS::power_series_compose_checked(f, g, 3);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    // f(2x) = 1 + 2x + 4x²
    ASSERT_EQ(result.size(), 3u) << "Result has 3 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(2))) << "c[1] = 2";
    EXPECT_TRUE(test_same_expression(result[2], num(4))) << "c[2] = 4";
}

TEST(SeriesEngine, PowerSeriesComposeG0Nonzero) {
    Coeffs f = {num(1), num(1)};
    // g(x) = 1 + x (g(0) = 1 ≠ 0)
    Coeffs g = {num(1), num(1)};

    auto result = LMCAS::power_series_compose_checked(f, g, 3);
    EXPECT_TRUE((!result && result.error().code == LMCAS::CasErrc::DomainError)) << "g(0) != 0 is a checked domain error";
}

TEST(SeriesEngine, PowerSeriesComposeExpLike) {
    // f(x) = 1 + x + x²/2 + x³/6 (exp(x) truncated to order 4)
    // Using Rational for 1/2 and 1/6
    Coeffs f = {num(1), num(1),
                SymbolicExpr::number(Rational(1, 2)),
                SymbolicExpr::number(Rational(1, 6))};
    // g(x) = 2x
    Coeffs g = {num(0), num(2)};

    auto checked_result = LMCAS::power_series_compose_checked(f, g, 4);
    ASSERT_TRUE(checked_result.has_value());
    const auto &result = checked_result.value();

    // exp(2x) = 1 + 2x + 2x² + (4/3)x³
    ASSERT_EQ(result.size(), 4u) << "Result has 4 coefficients";
    EXPECT_TRUE(test_same_expression(result[0], num(1))) << "c[0] = 1";
    EXPECT_TRUE(test_same_expression(result[1], num(2))) << "c[1] = 2";
    EXPECT_TRUE(test_same_expression(result[2], num(2))) << "c[2] = (1/2)*4 = 2";
    EXPECT_TRUE(test_same_expression(result[3], SymbolicExpr::number(Rational(4, 3)))) << "c[3] = (1/6)*8 = 4/3";
}

TEST(SeriesEngine, PowerSeriesMultiplyOrderZero) {
    Coeffs a = {num(1), num(2)};
    Coeffs b = {num(3), num(4)};

    auto result = LMCAS::power_series_multiply_checked(a, b, 0);
    EXPECT_TRUE((!result && result.error().code == LMCAS::CasErrc::InvalidArgument)) << "order zero is a checked invalid argument";
}

TEST(SeriesEngine, PowerSeriesMultiplyCheckedContracts) {
    Coeffs a = {num(1), num(1)};
    Coeffs b = {num(1), num(1)};
    auto checked_product = LMCAS::power_series_multiply_checked(a, b, 3);
    ASSERT_TRUE(checked_product.has_value()) << "checked multiply succeeds";
    ASSERT_EQ(checked_product.value().size(), 3u) << "checked multiply has requested order";
    ASSERT_TRUE(checked_product.value()[1]);
    EXPECT_EQ(checked_product.value()[1]->to_string(), num(2)->to_string())
        << "checked multiply coefficient c[1] = 2";

    auto bad_order = LMCAS::power_series_multiply_checked(a, b, 0);
    EXPECT_TRUE((!bad_order &&
                 bad_order.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked multiply rejects non-positive order";

    Coeffs with_null = {num(1), nullptr};
    auto null_coeff = LMCAS::power_series_multiply_checked(with_null, b, 2);
    EXPECT_TRUE((!null_coeff &&
                 null_coeff.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked multiply rejects null coefficients";
}

TEST(SeriesEngine, PowerSeriesComposeCheckedContracts) {
    Coeffs f = {num(1), num(1)};
    Coeffs g = {num(0), num(2)};
    auto checked_compose = LMCAS::power_series_compose_checked(f, g, 3);
    ASSERT_TRUE(checked_compose.has_value()) << "checked compose succeeds";
    ASSERT_GE(checked_compose.value().size(), 2u);
    ASSERT_TRUE(checked_compose.value()[0]);
    ASSERT_TRUE(checked_compose.value()[1]);
    EXPECT_EQ(checked_compose.value()[0]->to_string(), num(1)->to_string())
        << "checked compose coefficient c[0] = 1";
    EXPECT_EQ(checked_compose.value()[1]->to_string(), num(2)->to_string())
        << "checked compose coefficient c[1] = 2";

    Coeffs bad_g = {num(1), num(1)};
    auto nonzero_g0 = LMCAS::power_series_compose_checked(f, bad_g, 3);
    EXPECT_TRUE((!nonzero_g0 &&
                 nonzero_g0.error().code == LMCAS::CasErrc::DomainError))
        << "checked compose reports g(0) != 0 as domain error";
}

TEST(SeriesEngine, PowerSeriesMultiplyCancelled) {
    Coeffs a = {num(1), num(1)};
    Coeffs b = {num(1), num(1)};

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::power_series_multiply_checked(a, b, 3, context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked multiply observes cancelled context";
}
