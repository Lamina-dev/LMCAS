#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprMath, MathInvalidInputs) {
    auto null_sin = LMCAS::sin(nullptr);
    auto null_pow = LMCAS::pow(SymbolicExpr::number(2), nullptr);
    auto null_clamp = LMCAS::clamp(SymbolicExpr::number(1), nullptr,
                                   SymbolicExpr::number(2));
    EXPECT_TRUE((!null_sin &&
                 null_sin.error().code == LMCAS::CasErrc::InvalidArgument))
        << "std.math Expr unary wrappers reject null input";
    EXPECT_TRUE((!null_pow &&
                 null_pow.error().code == LMCAS::CasErrc::InvalidArgument))
        << "std.math Expr binary wrappers reject null input";
    EXPECT_TRUE((!null_clamp &&
                 null_clamp.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "std.math Expr clamp rejects null input";
}

TEST(ExprMath, MathBudget) {
    LMCAS::ResourceLimits exhausted_math_limits;
    exhausted_math_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_math_context(exhausted_math_limits);
    auto exhausted_sin =
        LMCAS::sin(SymbolicExpr::number(1), exhausted_math_context);
    EXPECT_TRUE((!exhausted_sin &&
                 exhausted_sin.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "std.math Expr wrappers observe the computation budget";
}

TEST(ExprMath, TrigonometricValues) {
    auto pi_constant = LMCAS::pi();
    auto quarter_pi = SymbolicExpr::divide(pi_constant.value(),
                                           SymbolicExpr::number(4));
    auto math_sin_quarter_pi = LMCAS::sin(quarter_pi);
    auto math_cos_zero = LMCAS::cos(SymbolicExpr::number(0));
    auto math_tan_zero = LMCAS::tan(SymbolicExpr::number(0));
    auto sin_quarter_pi_value =
        math_sin_quarter_pi
            ? LMCAS::evalf(*math_sin_quarter_pi.value())
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "sin construction failed", "test");
    auto cos_zero_value =
        math_cos_zero ? LMCAS::evalf(*math_cos_zero.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "cos construction failed", "test");
    auto tan_zero_value =
        math_tan_zero ? LMCAS::evalf(*math_tan_zero.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "tan construction failed", "test");
    {
        const double actual_value = (sin_quarter_pi_value.value().value);
        const double expected_value = (std::sqrt(0.5));
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    EXPECT_NEAR(cos_zero_value.value().value, 1.0, 0.0) << "cos(0) evaluates explicitly";
    EXPECT_NEAR(tan_zero_value.value().value, 0.0, 0.0) << "tan(0) evaluates explicitly";
}

TEST(ExprMath, RootsAndPowers) {
    auto math_sqrt_four = LMCAS::sqrt(SymbolicExpr::number(4));
    auto math_pow_two_three =
        LMCAS::pow(SymbolicExpr::number(2), SymbolicExpr::number(3));
    auto sqrt_four_value =
        math_sqrt_four ? LMCAS::evalf(*math_sqrt_four.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "sqrt construction failed", "test");
    auto pow_two_three_value =
        math_pow_two_three
            ? LMCAS::evalf(*math_pow_two_three.value())
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "pow construction failed", "test");
    EXPECT_NEAR(sqrt_four_value.value().value, 2.0, 0.0) << "sqrt(4) evaluates explicitly";
    EXPECT_NEAR(pow_two_three_value.value().value, 8.0, 0.0) << "pow(2, 3) evaluates explicitly";
}

TEST(ExprMath, InverseTrigonometricValues) {
    auto math_asin_half = LMCAS::asin(SymbolicExpr::number(0.5));
    auto math_acos_one = LMCAS::acos(SymbolicExpr::number(1));
    auto math_atan_one = LMCAS::atan(SymbolicExpr::number(1));
    auto asin_half_value =
        math_asin_half ? LMCAS::evalf(*math_asin_half.value())
                       : LMCAS::Result<LMCAS::ApproxReal>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "asin construction failed", "test");
    auto acos_one_value =
        math_acos_one ? LMCAS::evalf(*math_acos_one.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "acos construction failed", "test");
    auto atan_one_value =
        math_atan_one ? LMCAS::evalf(*math_atan_one.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "atan construction failed", "test");
    {
        const double actual_value = (asin_half_value.value().value);
        const double expected_value = (std::asin(0.5));
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    EXPECT_NEAR(acos_one_value.value().value, 0.0, 0.0) << "acos(1) evaluates explicitly";
    {
        const double actual_value = (atan_one_value.value().value);
        const double expected_value = (std::atan(1.0));
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(ExprMath, ExponentialAndLogarithmicValues) {
    auto e_constant = LMCAS::e();
    auto math_exp_zero = LMCAS::exp(SymbolicExpr::number(0));
    auto math_log_e = LMCAS::log(e_constant.value());
    auto math_log10_hundred = LMCAS::log10(SymbolicExpr::number(100));
    auto exp_zero_value =
        math_exp_zero ? LMCAS::evalf(*math_exp_zero.value())
                      : LMCAS::Result<LMCAS::ApproxReal>::failure(
                            LMCAS::CasErrc::InternalInvariant,
                            "exp construction failed", "test");
    auto log_e_value =
        math_log_e ? LMCAS::evalf(*math_log_e.value())
                   : LMCAS::Result<LMCAS::ApproxReal>::failure(
                         LMCAS::CasErrc::InternalInvariant,
                         "log construction failed", "test");
    auto log10_hundred_value =
        math_log10_hundred
            ? LMCAS::evalf(*math_log10_hundred.value())
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "log10 construction failed", "test");
    EXPECT_NEAR(exp_zero_value.value().value, 1.0, 0.0) << "exp(0) evaluates explicitly";
    EXPECT_NEAR(log_e_value.value().value, 1.0, 1e-12) << "log(e) evaluates explicitly";
    EXPECT_NEAR(log10_hundred_value.value().value, 2.0, 1e-12) << "log10(100) evaluates explicitly";
}

TEST(ExprMath, RoundingValues) {
    auto math_floor = LMCAS::floor(SymbolicExpr::number(2.75));
    auto math_ceil = LMCAS::ceil(SymbolicExpr::number(2.25));
    auto math_round = LMCAS::round(SymbolicExpr::number(-2.5));
    auto floor_value =
        math_floor ? LMCAS::evalf(*math_floor.value())
                   : LMCAS::Result<LMCAS::ApproxReal>::failure(
                         LMCAS::CasErrc::InternalInvariant,
                         "floor construction failed", "test");
    auto ceil_value =
        math_ceil ? LMCAS::evalf(*math_ceil.value())
                  : LMCAS::Result<LMCAS::ApproxReal>::failure(
                        LMCAS::CasErrc::InternalInvariant,
                        "ceil construction failed", "test");
    auto round_value =
        math_round ? LMCAS::evalf(*math_round.value())
                   : LMCAS::Result<LMCAS::ApproxReal>::failure(
                         LMCAS::CasErrc::InternalInvariant,
                         "round construction failed", "test");
    EXPECT_NEAR(floor_value.value().value, 2.0, 0.0) << "floor(2.75) evaluates explicitly";
    EXPECT_NEAR(ceil_value.value().value, 3.0, 0.0) << "ceil(2.25) evaluates explicitly";
    EXPECT_NEAR(round_value.value().value, -3.0, 0.0) << "round(-2.5) evaluates explicitly";
}

TEST(ExprMath, ClampingValues) {
    auto math_clamp = LMCAS::clamp(SymbolicExpr::number(7),
                                   SymbolicExpr::number(0),
                                   SymbolicExpr::number(5));
    auto clamp_value =
        math_clamp ? LMCAS::evalf(*math_clamp.value())
                   : LMCAS::Result<LMCAS::ApproxReal>::failure(
                         LMCAS::CasErrc::InternalInvariant,
                         "clamp construction failed", "test");
    EXPECT_NEAR(clamp_value.value().value, 5.0, 0.0) << "clamp(7, 0, 5) evaluates explicitly";
}

} // namespace
