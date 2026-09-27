#include "limit_result.hpp"
#include "integration.hpp"
#include "internal/symbolic_ast.hpp"
#include "test_common.hpp"
#include <cmath>

using namespace LMCAS;

TEST(LmcasSingularity, OneSidedLimit) {
    auto x = SymbolicExpr::variable("x");

    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto zero = SymbolicExpr::number(0);

    auto right = limit_checked(expr, "x", zero, LimitDirection::FromAbove);
    auto left = limit_checked(expr, "x", zero, LimitDirection::FromBelow);
    ASSERT_TRUE(right.has_value()) << "right reciprocal limit completes";
    ASSERT_TRUE(left.has_value()) << "left reciprocal limit completes";
    EXPECT_TRUE(std::holds_alternative<PositiveInfinityLimit>(
        right.value().value))
        << "right reciprocal limit is positive infinity";
    EXPECT_TRUE(std::holds_alternative<NegativeInfinityLimit>(
        left.value().value))
        << "left reciprocal limit is negative infinity";
}

TEST(LmcasSingularity, ImproperIntegralSingularity) {
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");

    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(-1));

    auto lower = SymbolicExpr::number(-1);
    auto upper = SymbolicExpr::number(1);

    auto res = integrator.integrate_def(*expr, "x", *lower, *upper);
    EXPECT_FALSE((res.has_value())) << ("opposite logarithmic divergences have no ordinary integral");
    if (!res) {
        EXPECT_TRUE((res.error().code == CasErrc::DomainError)) << ("ordinary nonexistence is a domain error, not an unevaluated principal value");
    }
}
