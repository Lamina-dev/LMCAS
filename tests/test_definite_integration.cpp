#include "integration.hpp"
#include "expr.hpp"
#include "numeric_evaluation.hpp"
#include "test_common.hpp"
#include "internal/symbolic_ast.hpp"

#include <cmath>
#include <string>

using namespace LMCAS;

namespace {
using Expr = std::shared_ptr<SymbolicExpr>;

Expr integer(int value) { return SymbolicExpr::number(value); }
Expr power(const Expr &base, int exponent) { return SymbolicExpr::power(base, integer(exponent)); }
Expr add_expr(const Expr &a, const Expr &b) { return SymbolicExpr::add(a, b); }
Expr multiply(const Expr &a, const Expr &b) { return SymbolicExpr::multiply(a, b); }
Expr inverse(const Expr &value) { return power(value, -1); }

void expect_finite(const Expr &expression, const Expr &lower, const Expr &upper,
                   double expected, const char *message) {
    Integrator integrator;
    auto result = integrator.integrate_def(*expression, "x", *lower, *upper);
    ASSERT_TRUE((result.has_value())) << (message);
    if (!result)
        return;
    auto number = evaluate_numeric(result.value());
    ASSERT_TRUE((number.has_value())) << (message);
    if (number) {
        {
            const auto actual_value = number.value().value;
            const auto expected_value = expected;
            EXPECT_TRUE(std::isfinite(actual_value)) << (message);
            EXPECT_NEAR(actual_value, expected_value, 1e-9) << (message);
        }
    }
}

void expect_infinity(const Expr &expression, const Expr &lower, const Expr &upper,
                     int sign, const char *message) {
    Integrator integrator;
    auto result = integrator.integrate_def(*expression, "x", *lower, *upper);
    ASSERT_TRUE((result.has_value())) << (message);
    if (!result)
        return;
    const auto expected = SymbolicExpr::infinity(sign);
    EXPECT_TRUE((result.value().compare(expected) == 0)) << (message);
}

void expect_error(const Expr &expression, const Expr &lower, const Expr &upper,
                  CasErrc code, const char *message) {
    Integrator integrator;
    auto result = integrator.integrate_def(*expression, "x", *lower, *upper);
    EXPECT_FALSE((result.has_value())) << (message);
    if (!result) {
        EXPECT_TRUE((result.error().code == code)) << (message);
    }
}

void expect_integral_node(const SymbolicExpr &actual, const Expr &expression,
                          const Expr &lower, const Expr &upper, const char *message) {
    auto expected = parse_expr("Integral(" + expression->to_string() + ",x," +
                               lower->to_string() + "," + upper->to_string() + ")");
    ASSERT_TRUE((expected.has_value())) << ("public parser constructs the complete integral obligation");
    if (expected) {
        EXPECT_TRUE((actual.compare(expected.value()) == 0)) << (message);
    }
}

void expect_unresolved(const Expr &expression, const Expr &lower, const Expr &upper,
                       const char *message) {
    Integrator integrator;
    auto result = integrator.integrate_def(*expression, "x", *lower, *upper);
    ASSERT_TRUE((result.has_value())) << (message);
    if (!result)
        return;
    expect_integral_node(result.value(), expression, lower, upper, message);
}

class SuppliedPrimitive final : public BuiltInIntegrationStrategy {
  public:
    explicit SuppliedPrimitive(Expr expression) : expression_(std::move(expression)) {}
    std::string name() const override { return "SuppliedPrimitive"; }
    Result<Expr> try_integrate_raw(const SymbolicExpr &, const std::string &,
                                   Integrator &, ComputationContext &, int) override { return expression_; }

  private:
    Expr expression_;
};

TEST(LmcasDefiniteIntegration, PrimitiveCertificates) {
    const auto x = SymbolicExpr::variable("x");
    const auto expression = power(x, 2);
    const auto nested_integral = SymbolicExpr::make_integral(SymbolicExpr::sin(power(x, 2)), "x");
    for (const auto &candidate : {integer(0), add_expr(x, nested_integral)}) {
        Integrator integrator;
        auto installed = integrator.add_strategy(std::make_unique<SuppliedPrimitive>(candidate), 0);
        ASSERT_TRUE((installed.has_value())) << ("install a primitive candidate through the strategy API");
        if (!installed)
            continue;
        auto result = integrator.integrate_def(*expression, "x", *integer(0), *integer(1));
        ASSERT_TRUE((result.has_value())) << ("uncertified primitive leaves an unresolved integral");
        if (!result)
            continue;
        expect_integral_node(result.value(), expression, integer(0), integer(1),
                             "uncertified primitive preserves the entire original integral and bounds");
    }
}

TEST(LmcasDefiniteIntegration, AbsoluteValueSegments) {
    const auto x = SymbolicExpr::variable("x");
    const auto absolute_result = LMCAS::abs(x);
    ASSERT_TRUE((absolute_result.has_value())) << ("construct an absolute-value expression");
    if (!absolute_result)
        return;
    const auto absolute = absolute_result.value();
    expect_finite(absolute, integer(-1), integer(1), 1.0,
                  "absolute-value branch primitive is differentiated only on certified segments");
    const auto shifted = LMCAS::abs(add_expr(x, integer(-2)));
    ASSERT_TRUE((shifted.has_value())) << ("construct a shifted absolute value");
    if (shifted)
        expect_finite(shifted.value(), integer(1), integer(3), 1.0,
                      "the branch cut follows the argument root rather than the origin");
}

TEST(LmcasDefiniteIntegration, PiecewiseSegments) {
    const auto x = SymbolicExpr::variable("x");
    const auto negative = detail::make_node<RelationalNode>(
        detail::node(x), detail::node(integer(0)), RelationOp::LT);
    const auto make_piecewise = [&](const Expr &left, const Expr &right) {
        return detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{{detail::node(left), negative}},
            right ? detail::node(right) : nullptr));
    };
    auto linear = make_piecewise(x, multiply(integer(2), x));
    expect_finite(linear, integer(-1), integer(1), 0.5,
                  "piecewise branches use independent primitives on both open segments");
    expect_finite(linear, integer(1), integer(-1), -0.5,
                  "piecewise integral reverses orientation after combining segments");
    expect_finite(make_piecewise(integer(1), integer(3)), integer(-1), integer(1), 4.0,
                  "a finite jump does not invalidate an ordinary definite integral");
    expect_finite(make_piecewise(integer(-1), SymbolicExpr::ln(x)), integer(-2), integer(-1), -1.0,
                  "the undefined inactive branch does not contaminate the selected segment");
    expect_error(make_piecewise(integer(1), nullptr), integer(-1), integer(1), CasErrc::DomainError,
                 "a missing branch on an open segment is not a complete integral");
    expect_infinity(make_piecewise(inverse(x), multiply(integer(-1), inverse(x))),
                    integer(-1), integer(1), -1,
                    "piecewise pole contributions combine without principal-value cancellation");
}

TEST(LmcasDefiniteIntegration, RegularIntegrals) {
    const auto x = SymbolicExpr::variable("x");
    expect_finite(power(x, 2), integer(0), integer(1), 1.0 / 3.0, "polynomial integral");
    expect_finite(SymbolicExpr::sin(x), integer(0), SymbolicExpr::number(std::acos(-1.0)),
                  2.0, "sine integral with binary64 endpoint");
    expect_finite(inverse(add_expr(integer(1), power(x, 2))), integer(-1), integer(1),
                  std::acos(-1.0) / 2, "arctangent primitive");
    expect_finite(SymbolicExpr::exp(x), integer(0), integer(1), std::exp(1.0) - 1,
                  "exponential continuous primitive");
    const auto a = SymbolicExpr::variable("a");
    const auto b = SymbolicExpr::variable("b");
    Integrator integrator;
    auto result = integrator.integrate_def(*x, "x", *a, *b);
    ASSERT_TRUE((result.has_value())) << ("polynomial symbolic bounds remain supported");
    if (result) {
        auto value = evaluate_numeric(result.value(), {{"a", 2.0}, {"b", 4.0}});
        ASSERT_TRUE((value.has_value())) << ("symbolic endpoint expression is numerically evaluable");
        if (value) {
            const auto actual_value = value.value().value;
            const auto expected_value = 6.0;
            EXPECT_TRUE(std::isfinite(actual_value)) << ("integral x[a,b]");
            EXPECT_NEAR(actual_value, expected_value, 1e-9) << ("integral x[a,b]");
        }
    }
}

TEST(LmcasDefiniteIntegration, OrdinaryPoles) {
    const auto x = SymbolicExpr::variable("x");
    const auto square_pole = power(x, -2);
    expect_infinity(square_pole, integer(-1), integer(1), 1, "double pole is not endpoint subtraction");
    const auto parsed_pole = parse_expr("x^-2");
    ASSERT_TRUE((parsed_pole.has_value())) << ("raw signed-power fixture parses");
    if (parsed_pole) {
        expect_infinity(parsed_pole.value(), integer(-1), integer(1), 1,
                        "raw arithmetic exponent preserves both real pole neighborhoods");
    }
    expect_infinity(multiply(integer(2), square_pole), integer(-1), integer(1), 1,
                    "positive pole coefficient");
    expect_infinity(multiply(integer(-1), square_pole), integer(-1), integer(1), -1,
                    "negative pole coefficient");
    expect_error(inverse(x), integer(-1), integer(1), CasErrc::DomainError,
                 "opposite divergences do not form a principal value");
    expect_infinity(square_pole, integer(0), integer(1), 1, "endpoint double pole");
    expect_infinity(square_pole, integer(1), integer(-1), -1, "reversal negates divergent integral");
    expect_finite(square_pole, integer(1), integer(2), 0.5, "positive regular interval");
    expect_finite(square_pole, integer(-2), integer(-1), 0.5, "negative regular interval");
    expect_finite(square_pole, integer(2), integer(1), -0.5, "reversed regular interval");
    const auto shifted = add_expr(x, integer(-2));
    expect_infinity(power(shifted, -2), integer(1), integer(3), 1, "translated double pole");
    expect_infinity(multiply(integer(2), power(shifted, -4)), integer(1), integer(3), 1,
                    "translated fourth-order pole");
    expect_error(multiply(integer(-2), power(shifted, -3)), integer(1), integer(3),
                 CasErrc::DomainError, "translated odd-order pole");
    const auto irrational = add_expr(power(x, 2), integer(-2));
    expect_infinity(power(irrational, -2), integer(-2), integer(2), 1,
                    "both irrational double poles are isolated exactly");
    expect_error(inverse(irrational), integer(-2), integer(2), CasErrc::DomainError,
                 "irrational simple poles do not cancel");
    const auto rational_cut = multiply(add_expr(multiply(integer(3), x), integer(-1)), irrational);
    expect_infinity(power(rational_cut, -2), integer(0), integer(1), 1,
                    "a nondyadic rational pole is recovered among algebraic polynomial cuts");
    const auto high_degree = add_expr(power(x, 5), integer(-2));
    expect_infinity(power(high_degree, -2), integer(0), integer(2), 1,
                    "arbitrary-degree algebraic double pole needs no closed primitive");
    expect_error(inverse(high_degree), integer(0), integer(2), CasErrc::DomainError,
                 "arbitrary-degree odd pole needs no closed primitive");
    const auto mixed = multiply(power(x, -2), inverse(add_expr(x, integer(-2))));
    expect_error(mixed, integer(-1), integer(3), CasErrc::DomainError,
                 "all singular contributions are combined before returning infinity");
    const std::string huge = "1" + std::string(400, '0');
    expect_error(inverse(x), SymbolicExpr::number(BigInt("-" + huge)),
                 SymbolicExpr::number(BigInt(huge)), CasErrc::DomainError,
                 "exact bounds outside double range still straddle a simple pole");
}

TEST(LmcasDefiniteIntegration, HolesAndRealBranches) {
    const auto x = SymbolicExpr::variable("x");
    const auto numerator = add_expr(power(x, 2), integer(-1));
    const auto denominator = add_expr(x, integer(-1));
    const auto hole = multiply(numerator, inverse(denominator));
    expect_finite(hole, integer(0), integer(2), 4.0, "removable hole joins finite one-sided contributions");
    expect_error(hole, integer(1), integer(1), CasErrc::DomainError,
                 "zero length does not fill the original hole");
    expect_finite(inverse(x), integer(1), integer(2), std::log(2.0), "positive logarithmic primitive");
    expect_finite(inverse(x), integer(-2), integer(-1), -std::log(2.0),
                  "negative logarithmic primitive is proved on its real interval");
    expect_finite(inverse(add_expr(x, integer(-3))), integer(1), integer(2), -std::log(2.0),
                  "translated negative logarithmic primitive");
    expect_finite(SymbolicExpr::ln(x), integer(1), integer(2), 2 * std::log(2.0) - 1,
                  "restricted logarithmic integrand uses interval-authorized derivative cancellation");
    expect_error(SymbolicExpr::ln(x), integer(-2), integer(-1), CasErrc::DomainError,
                 "nonreal open interval is rejected");
    const auto zero_power = power(x, 0);
    expect_finite(zero_power, integer(-1), integer(1), 2.0, "zero exponent retains its removable domain boundary");
    expect_error(zero_power, integer(0), integer(0), CasErrc::DomainError,
                 "zero exponent does not swallow an illegal singleton");
    const auto nested_inverse = inverse(inverse(x));
    expect_finite(nested_inverse, integer(-1), integer(1), 0.0,
                  "nested inverse retains the original hole while integrating its reduced value");
}

TEST(LmcasDefiniteIntegration, InfiniteAndIntegrableBoundaries) {
    const auto x = SymbolicExpr::variable("x");
    expect_finite(SymbolicExpr::power(x, SymbolicExpr::number(Rational(BigInt(-1), BigInt(2)))),
                  integer(0), integer(1), 2.0, "integrable fractional-power endpoint");
    expect_finite(power(x, -2), integer(1), SymbolicExpr::infinity(1), 1.0,
                  "convergent positive rational tail");
    expect_finite(power(x, -2), SymbolicExpr::infinity(-1), integer(-1), 1.0,
                  "convergent negative rational tail");
    expect_finite(inverse(add_expr(integer(1), power(x, 2))), SymbolicExpr::infinity(-1),
                  SymbolicExpr::infinity(1), std::acos(-1.0), "two convergent rational tails");
    expect_infinity(inverse(x), integer(1), SymbolicExpr::infinity(1), 1,
                    "logarithmically divergent positive tail");
    expect_infinity(inverse(x), SymbolicExpr::infinity(-1), integer(-1), -1,
                    "logarithmically divergent negative tail");
    expect_error(x, SymbolicExpr::infinity(-1), SymbolicExpr::infinity(1), CasErrc::DomainError,
                 "opposite polynomial tails are not a principal value");
    expect_infinity(power(x, 2), SymbolicExpr::infinity(-1), SymbolicExpr::infinity(1), 1,
                    "same-sign polynomial tails");
    expect_error(SymbolicExpr::sin(x), integer(0), SymbolicExpr::infinity(1), CasErrc::DomainError,
                 "oscillatory primitive has no ordinary improper integral");
}

TEST(LmcasDefiniteIntegration, UncertifiedAndErrors) {
    const auto x = SymbolicExpr::variable("x");
    const auto a = SymbolicExpr::variable("a");
    const auto b = SymbolicExpr::variable("b");
    expect_unresolved(inverse(add_expr(x, a)), integer(-1), integer(1),
                      "symbolic pole position is not guessed");
    expect_unresolved(inverse(x), a, b, "unordered restricted symbolic endpoints preserve whole integral");
    expect_unresolved(inverse(SymbolicExpr::sin(x)), integer(-1), integer(1),
                      "unisolated transcendental singularity set preserves whole integral");
    expect_unresolved(add_expr(x, SymbolicExpr::sin(power(x, 2))), integer(0), integer(1),
                      "uncertified non-elementary primitive preserves the whole integral");
    expect_finite(x, integer(2), integer(2), 0.0, "legal zero-length interval");
    expect_error(x, x, integer(1), CasErrc::InvalidArgument, "variable-dependent endpoint");
    Integrator integrator;
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    auto budget = integrator.integrate_def_checked(*power(x, -2), "x", *integer(-1), *integer(1), exhausted);
    EXPECT_FALSE((budget.has_value())) << ("definite integration preserves resource exhaustion");
    if (!budget) {
        EXPECT_TRUE((budget.error().code == CasErrc::ResourceLimit)) << ("resource error code");
    }
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto cancellation = integrator.integrate_def_checked(*x, "x", *integer(0), *integer(1), cancelled);
    EXPECT_FALSE((cancellation.has_value())) << ("definite integration preserves cancellation");
    if (!cancellation) {
        EXPECT_TRUE((cancellation.error().code == CasErrc::Cancelled)) << ("cancellation error code");
    }
}
} // namespace
