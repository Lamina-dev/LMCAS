/**
 * @file test_calculus_utils_inverse.cpp
 * @brief 反函数分支与反函数求导约定。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"
#include <string>
#include <cmath>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string &name) { return SE::variable(name); }

TEST(CalculusUtilsInverse, LinearInverse) {
    auto x = var("x");
    auto expression = SE::add(SE::multiply(num(2), x), num(1));
    auto result = LMCAS::inverse_function_checked(
        expression, "x", num(5));
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_EQ(result.value().size(), 1u);
    auto value = test_numeric_eval(result.value().front());
    ASSERT_TRUE(value.has_value());
    EXPECT_NEAR(*value, 2.0, 1e-9);
}

TEST(CalculusUtilsInverse, QuadraticInverseReturnsBothBranchesOnce) {
    auto x = var("x");
    auto expression = SE::power(x, num(2));
    auto result = LMCAS::inverse_function_checked(
        expression, "x", num(4));
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_EQ(result.value().size(), 2u);

    bool found_negative = false;
    bool found_positive = false;
    for (const auto &root : result.value()) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        ASSERT_TRUE(std::isfinite(*value));
        if (std::abs(*value + 2.0) <= 1e-9) {
            EXPECT_FALSE(found_negative) << "negative branch is unique";
            found_negative = true;
        } else if (std::abs(*value - 2.0) <= 1e-9) {
            EXPECT_FALSE(found_positive) << "positive branch is unique";
            found_positive = true;
        } else {
            ADD_FAILURE() << "unexpected quadratic inverse root " << *value;
        }
    }
    EXPECT_TRUE(found_negative);
    EXPECT_TRUE(found_positive);
}

static bool inverse_has_error(
    const LMCAS::SymbolicExprVectorResult &result, LMCAS::CasErrc code) {
    return !result && result.error().code == code;
}

TEST(CalculusUtilsInverse, InverseSolutionSets) {
    auto x = var("x");
    auto linear = SE::add(SE::multiply(num(2), x), num(1));

    auto finite = LMCAS::inverse_function_checked(linear, "x", num(5));
    EXPECT_TRUE((finite && finite.value().size() == 1)) << "checked inverse function returns finite exact candidates";

    auto empty = LMCAS::inverse_function_checked(num(1), "x", num(2));
    EXPECT_TRUE((empty && empty.value().empty())) << "checked inverse function preserves a proven empty set";

    auto null_input = LMCAS::inverse_function_checked(nullptr, "x", num(1));
    EXPECT_TRUE((inverse_has_error(null_input, LMCAS::CasErrc::InvalidArgument))) << "checked inverse function rejects null expression";

    auto null_target = LMCAS::inverse_function_checked(x, "x", nullptr);
    EXPECT_TRUE((inverse_has_error(null_target, LMCAS::CasErrc::InvalidArgument))) << "checked inverse function rejects null target";

    auto trigonometric = LMCAS::inverse_function_checked(
        SE::sin(x), "x", num(1));
    EXPECT_TRUE((!trigonometric && trigonometric.error().code == LMCAS::CasErrc::Inconclusive)) << "finite inverse projection does not discard infinitely many sine preimages";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::inverse_function_checked(
        linear, "x", num(5), limited_context);
    EXPECT_TRUE((inverse_has_error(limited, LMCAS::CasErrc::ResourceLimit))) << "checked inverse function observes exhausted step budget";
}

TEST(CalculusUtilsInverse, LinearInverseDerivative) {
    auto x = var("x");
    auto expression = SE::add(SE::multiply(num(2), x), num(1));
    auto result = LMCAS::inverse_derivative_checked(
        expression, "x", num(5));
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_NE(result.value(), nullptr);
    auto value = test_numeric_eval(result.value());
    ASSERT_TRUE(value.has_value());
    ASSERT_TRUE(std::isfinite(*value));
    EXPECT_NEAR(*value, 0.5, 1e-9);
}

TEST(CalculusUtilsInverse, MultibranchInverseDerivative) {
    auto x = var("x");
    auto f = SE::power(x, num(2));
    auto result = LMCAS::inverse_derivative_checked(f, "x", num(4));
    EXPECT_TRUE((!result && result.error().code == LMCAS::CasErrc::Inconclusive)) << "checked inverse derivative rejects multiple branches";
}

TEST(CalculusUtilsInverse, CubicInverseDerivative) {
    /**
     * @brief 验证 f(x) = x^3 的反函数导数 (f^{-1})'(8) = 1/f'(2) = 1/12。
     * f'(x) = 3x^2，f^{-1}(8) = 2。
     */
    auto x = var("x");
    auto f = SE::power(x, num(3));
    auto checked = LMCAS::inverse_derivative_checked(f, "x", num(8));
    ASSERT_TRUE((checked.has_value())) << (checked ? "cubic inverse derivative succeeds" : checked.error().message);
    if (!checked)
        return;
    auto result = std::move(checked.value());
    EXPECT_TRUE((result != nullptr)) << "inverse_derivative(x^3, 8) non-null";
    if (result) {
        auto val = test_numeric_eval(result);
        if (val) {
            {
                const double actual_value = (*val);
                const double expected_value = (1.0 / 12.0);
                const double tolerance = (1e-9);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}

TEST(CalculusUtilsInverse, CheckedInverseDerivative) {
    auto x = var("x");
    auto linear = SE::add(SE::multiply(num(2), x), num(1));

    auto checked = LMCAS::inverse_derivative_checked(linear, "x", num(5));
    ASSERT_TRUE((checked.has_value())) << "checked inverse derivative succeeds for a unique linear inverse";
    if (checked) {
        auto val = test_numeric_eval(checked.value());
        EXPECT_TRUE((val.has_value() && std::abs(*val - 0.5) < 1e-9)) << "checked inverse derivative of 2x+1 is 1/2";
    }

    auto cubic_derivative = LMCAS::inverse_derivative_checked(
        SE::power(x, num(3)), "x", num(8));
    ASSERT_TRUE((cubic_derivative.has_value())) << "checked inverse derivative accepts the unique real cubic branch";
    if (cubic_derivative) {
        auto val = test_numeric_eval(cubic_derivative.value());
        EXPECT_TRUE((val.has_value() && std::abs(*val - (1.0 / 12.0)) < 1e-9)) << "checked inverse derivative of x^3 at 8 is 1/12";
    }

    auto no_preimage = LMCAS::inverse_derivative_checked(num(1), "x", num(2));
    EXPECT_TRUE((!no_preimage &&
                 no_preimage.error().code == LMCAS::CasErrc::DomainError))
        << "checked inverse derivative reports DomainError for no inverse point";

    auto multibranch = LMCAS::inverse_derivative_checked(
        SE::power(x, num(2)), "x", num(4));
    EXPECT_TRUE((!multibranch &&
                 multibranch.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked inverse derivative rejects non-unique inverse branches";

    auto zero_derivative = LMCAS::inverse_derivative_checked(
        SE::power(x, num(3)), "x", num(0));
    EXPECT_TRUE((!zero_derivative &&
                 zero_derivative.error().code == LMCAS::CasErrc::DomainError))
        << "checked inverse derivative reports DomainError when f' is zero";

    auto null_input = LMCAS::inverse_derivative_checked(nullptr, "x", num(1));
    EXPECT_TRUE((!null_input &&
                 null_input.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked inverse derivative rejects null expression";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::inverse_derivative_checked(
        linear, "x", num(5), limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked inverse derivative observes exhausted step budget";
}
