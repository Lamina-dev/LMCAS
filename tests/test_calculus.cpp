#include "test_common.hpp"
#include "symbolic.hpp"
#include "limit_result.hpp"
#include "expr.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"
#include "assumption_context.hpp"
#include <optional>
#include <stdexcept>

using namespace LMCAS;

TEST(Calculus, Atan2Constant) {
    auto x = SymbolicExpr::variable("x");

    auto angle = SymbolicExpr::atan2(x, SymbolicExpr::number(1));
    auto derivative = angle->differentiate("x");
    auto value = evaluate_numeric(*derivative, {{"x", 0.0}});
    EXPECT_TRUE((value && value.value().is_finite())) << "atan2(x, 1) derivative is finite at zero";
    if (value) {
        EXPECT_NEAR(value.value().value, 1.0, 1e-12) << "d atan2(x, 1)/dx at zero equals one";
    }
}

TEST(Calculus, Atan2PartialDerivatives) {
    auto x = SymbolicExpr::variable("x");

    auto z = SymbolicExpr::variable("z");
    auto y_argument = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), z);
    auto x_argument = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto angle = SymbolicExpr::atan2(y_argument, x_argument);
    auto dx = LMCAS::differentiate(angle, "x");
    auto dz = LMCAS::differentiate(angle, "z");
    EXPECT_TRUE((dx && dz)) << "atan2 supports partial derivatives";
    if (dx && dz) {
        const NumericBindings point{{"x", 2.0}, {"z", 3.0}};
        auto dx_value = evaluate_numeric(*dx.value(), point);
        auto dz_value = evaluate_numeric(*dz.value(), point);
        EXPECT_TRUE((dx_value && dx_value.value().is_finite() &&
                     dz_value && dz_value.value().is_finite()))
            << "atan2 partial derivatives evaluate at a regular point";
        if (dx_value && dz_value) {
            {
                const double actual_value = (dx_value.value().value);
                const double expected_value = (5.0 / 58.0);
                const double tolerance = (1e-12);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
            {
                const double actual_value = (dz_value.value().value);
                const double expected_value = (3.0 / 58.0);
                const double tolerance = (1e-12);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
}

TEST(Calculus, UnsupportedDerivatives) {
    auto x = SymbolicExpr::variable("x");

    auto logarithm = SymbolicExpr::log(x, SymbolicExpr::number(10));
    auto log_derivative = LMCAS::differentiate(logarithm, "x");
    EXPECT_TRUE((!log_derivative &&
                 log_derivative.error().code == CasErrc::UnsupportedExpression))
        << "unsupported non-unary derivative does not become zero";

    auto floor_expression = LMCAS::floor(x);
    ASSERT_TRUE((floor_expression.has_value())) << "symbolic floor is constructible";
    if (floor_expression) {
        auto floor_derivative = LMCAS::differentiate(floor_expression.value(), "x");
        EXPECT_TRUE((!floor_derivative &&
                     floor_derivative.error().code == CasErrc::UnsupportedExpression))
            << "unsupported unary derivative reports UnsupportedExpression";
        bool threw = false;
        try {
            (void)floor_expression.value()->differentiate("x");
        } catch (const std::runtime_error &) {
            threw = true;
        }
        EXPECT_TRUE((threw)) << "unchecked unsupported differentiation throws";
    }
}

TEST(Calculus, FiniteLimit) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto lim = LMCAS::limit_checked(f, "x", SymbolicExpr::number(2));
    const auto *finite = lim ? std::get_if<FiniteLimit>(&lim.value().value) : nullptr;
    EXPECT_TRUE((finite && finite->value->compare(SymbolicExpr::number(3)) == 0)) << "limit(x+1, x->2) is exactly three";
}

TEST(Calculus, IndeterminateLimit) {
    auto x = SymbolicExpr::variable("x");

    auto num = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
    auto den = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto f = SymbolicExpr::divide(num, den);

    auto lim = LMCAS::limit_checked(f, "x", SymbolicExpr::number(1));
    const auto *finite = lim ? std::get_if<FiniteLimit>(&lim.value().value) : nullptr;
    EXPECT_TRUE((finite && finite->value->compare(SymbolicExpr::number(2)) == 0)) << "removable rational singularity has the exact finite limit two";
}

TEST(Calculus, LimitOutcomes) {
    auto x = SymbolicExpr::variable("x");

    auto finite = LMCAS::limit_checked(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        "x", SymbolicExpr::number(2));
    const auto *finite_value = finite
                                   ? std::get_if<LMCAS::FiniteLimit>(&finite.value().value)
                                   : nullptr;
    EXPECT_TRUE((finite_value && finite_value->value->compare(SymbolicExpr::number(3)) == 0)) << "checked finite limit returns a proved finite outcome";
    auto oscillatory = LMCAS::limit_checked(
        SymbolicExpr::sin(x), "x", SymbolicExpr::infinity());
    EXPECT_TRUE((oscillatory &&
                 std::holds_alternative<LMCAS::LimitDoesNotExist>(
                     oscillatory.value().value)))
        << "sin(x) at infinity is DoesNotExist rather than Inconclusive";
}

static void check_primitive(const std::shared_ptr<SymbolicExpr> &primitive,
                            const std::shared_ptr<SymbolicExpr> &integrand,
                            double expected_increment) {
    ASSERT_TRUE(primitive) << "primitive is available";
    auto derivative = LMCAS::differentiate(primitive, "x");
    ASSERT_TRUE(derivative.has_value()) << "primitive derivative is available";
    EXPECT_TRUE(test_proved_equivalent(derivative.value(), integrand))
        << "primitive derivative equals the integrand";
    auto lower = evaluate_numeric(*primitive, {{"x", 1.0}});
    auto upper = evaluate_numeric(*primitive, {{"x", 2.0}});
    EXPECT_TRUE((lower && upper && lower.value().is_finite() && upper.value().is_finite())) << "primitive evaluates at both regular endpoints";
    if (lower && upper) {
        const double actual_value = (upper.value().value - lower.value().value);
        const double expected_value = (expected_increment);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(Calculus, PolynomialIntegrals) {
    auto x = SymbolicExpr::variable("x");

    auto integ = x->integrate("x");
    check_primitive(integ, x, 1.5);

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto integ2 = x2->integrate("x");
    check_primitive(integ2, x2, 7.0 / 3.0);
}

TEST(Calculus, StaticIntegralReturnsCertifiedPrimitive) {
    auto x = SymbolicExpr::variable("x");
    auto integrand = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto primitive = SymbolicExpr::integral(integrand, "x");
    ASSERT_TRUE(primitive);

    auto derivative = primitive->differentiate("x");
    ASSERT_TRUE(derivative);
    EXPECT_TRUE(test_proved_equivalent(derivative, integrand));
}

TEST(Calculus, SinOverXExpressionLimitReturnsOne) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::divide(SymbolicExpr::sin(x), x);
    auto result = limit_expression_checked(
        expression, "x", SymbolicExpr::number(0));
    ASSERT_TRUE(result);
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->compare(SymbolicExpr::number(1)), 0);
}

TEST(Calculus, LogarithmicIntegral) {
    auto x = SymbolicExpr::variable("x");

    auto x_inv = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto integ = x_inv->integrate("x");
    check_primitive(integ, x_inv, std::log(2.0));
}

TEST(Calculus, IntegralLinearity) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto integ = f->integrate("x");

    check_primitive(integ, f, 2.5);
}

static void check_lambert_origin_derivatives(const std::shared_ptr<SymbolicExpr> &x) {
    for (const int scale : {1, 2}) {
        auto argument = scale == 1 ? x : SymbolicExpr::multiply(SymbolicExpr::number(scale), x);
        auto derivative = LMCAS::differentiate(SymbolicExpr::lambertw(argument), "x");
        ASSERT_TRUE((derivative.has_value())) << "Lambert derivative is constructible";
        if (!derivative) {
            continue;
        }
        auto direct = evaluate_numeric(*derivative.value(), {{"x", 0.0}});
        auto at_zero = derivative.value()->substitute("x", SymbolicExpr::number(0));
        auto substituted = evaluate_numeric(*at_zero);
        EXPECT_TRUE((direct && direct.value().is_finite() && direct.value().value == scale)) << "direct evaluation at zero includes the inner derivative";
        EXPECT_TRUE((substituted && substituted.value().is_finite() &&
                     substituted.value().value == scale))
            << "substitution at zero preserves the removable derivative value";
        auto exact = at_zero->simplify();
        EXPECT_TRUE((exact && detail::node(exact)->equals(
                                  *detail::node(SymbolicExpr::number(scale)))))
            << "zero-point derivative simplifies to its exact integer value";
    }
}

static void check_lambert_regular_and_branch_derivatives(const std::shared_ptr<SymbolicExpr> &x) {
    auto derivative = LMCAS::differentiate(SymbolicExpr::lambertw(x), "x");
    if (!derivative) {
        return;
    }
    const double e = std::exp(1.0);
    const auto regular = evaluate_numeric(*derivative.value(), {{"x", e}});
    EXPECT_TRUE((regular && regular.value().is_finite() &&
                 std::abs(regular.value().value - 1.0 / (2.0 * e)) <= 1e-12))
        << "W(e)=1 gives the independent regular derivative 1/(2e)";
    const auto branch = evaluate_numeric(*derivative.value(), {{"x", -1.0 / e}});
    EXPECT_TRUE((!branch && branch.error().code == CasErrc::DomainError)) << "the genuine W=-1 branch-point singularity is not removed";
}

TEST(Calculus, LambertDerivativeAtOrigin) {
    auto x = SymbolicExpr::variable("x");

    check_lambert_origin_derivatives(x);
    check_lambert_regular_and_branch_derivatives(x);
}

TEST(Calculus, OneSidedAndFinitePayloadContracts) {
    struct FiniteCase {
        const char *source;
        const char *point;
        LimitDirection direction;
        Rational expected;
    };
    const FiniteCase finite_cases[] = {
        {"abs(x)/x", "0", LimitDirection::FromBelow, Rational(-1)},
        {"abs(x)/x", "0", LimitDirection::FromAbove, Rational(1)},
        {"abs(x-2)/(x-2)", "2", LimitDirection::FromBelow, Rational(-1)},
        {"abs(x-2)/(x-2)", "2", LimitDirection::FromAbove, Rational(1)},
        {"sin(x)/x", "0", LimitDirection::Both, Rational(1)},
        {"x*sin(1/x)", "0", LimitDirection::Both, Rational(0)},
        {"(1-cos(x))/x^2", "0", LimitDirection::Both, Rational(1, 2)},
        {"(sin(x)-x)/x^3", "0", LimitDirection::Both, Rational(-1, 6)},
        {"(x-2)^12/(x-2)^10", "2", LimitDirection::Both, Rational(0)},
        {"x^0", "0", LimitDirection::Both, Rational(1)},
        {"1^(1/x)", "0", LimitDirection::Both, Rational(1)},
        {"1/x-1/x", "0", LimitDirection::Both, Rational(0)},
        {"1", "0", LimitDirection::Both, Rational(1)}};
    for (const auto &sample : finite_cases) {
        auto expression = parse_expr(sample.source);
        auto point = parse_expr(sample.point);
        EXPECT_TRUE((expression && point)) << "limit fixture parses";
        if (!expression || !point) {
            continue;
        }
        auto value = limit_checked(expression.value(), "x", point.value(), sample.direction);
        auto finite = value ? std::get_if<FiniteLimit>(&value.value().value) : nullptr;
        EXPECT_TRUE((finite != nullptr)) << sample.source;
        if (!finite) {
            continue;
        }
        ComputationContext context;
        auto proof = check_equivalent(finite->value, SymbolicExpr::number(sample.expected), context);
        EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << sample.source;
        auto numeric = evaluate_numeric(*finite->value);
        EXPECT_TRUE((numeric && numeric.value().is_finite())) << "finite payload is actually evaluable";
    }
}

TEST(Calculus, SignedInfiniteLimitContracts) {
    struct InfiniteCase {
        const char *source;
        const char *point;
        LimitDirection direction;
        int sign;
    };
    const InfiniteCase infinite_cases[] = {
        {"1/x", "0", LimitDirection::FromBelow, -1},
        {"1/x", "0", LimitDirection::FromAbove, 1},
        {"2/x^2", "0", LimitDirection::FromBelow, 1},
        {"2*x^-2", "0", LimitDirection::FromBelow, 1},
        {"2/x^2", "0", LimitDirection::FromAbove, 1},
        {"2/x^2", "0", LimitDirection::Both, 1},
        {"2/(x-3)^4", "3", LimitDirection::FromBelow, 1},
        {"-2/(x-3)^3", "3", LimitDirection::FromBelow, 1},
        {"-2/(x-3)^3", "3", LimitDirection::FromAbove, -1},
        {"2/(x-a)^4", "a", LimitDirection::Both, 1}};
    for (const auto &sample : infinite_cases) {
        auto expression = parse_expr(sample.source);
        auto point = parse_expr(sample.point);
        EXPECT_TRUE((expression && point)) << "infinite fixture parses";
        if (!expression || !point) {
            continue;
        }
        auto result = limit_checked(expression.value(), "x", point.value(), sample.direction);
        EXPECT_TRUE((result && (sample.sign > 0
                                    ? std::holds_alternative<PositiveInfinityLimit>(result.value().value)
                                    : std::holds_alternative<NegativeInfinityLimit>(result.value().value))))
            << sample.source;
    }
}

TEST(Calculus, LimitNonexistenceContracts) {
    for (const char *source : {"1/x", "abs(x)/x", "sin(1/x)"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "nonexistence fixture parses";
        if (!expression) {
            continue;
        }
        auto result = limit_checked(expression.value(), "x", SymbolicExpr::number(0));
        EXPECT_TRUE((result && std::holds_alternative<LimitDoesNotExist>(result.value().value))) << source;
    }
}

TEST(Calculus, ExtendedRealLimitProducts) {
    auto infinity = SymbolicExpr::infinity();
    for (int sign : {-1, 1}) {
        auto expression = SymbolicExpr::multiply(SymbolicExpr::number(2 * sign), infinity);
        auto result = limit_checked(expression, "x", SymbolicExpr::number(0));
        EXPECT_TRUE((result && (sign > 0
                                    ? std::holds_alternative<PositiveInfinityLimit>(result.value().value)
                                    : std::holds_alternative<NegativeInfinityLimit>(result.value().value))))
            << "extended-real multiplication retains its sign";
    }
    auto zero_infinity = SymbolicExpr::multiply(SymbolicExpr::number(0), infinity);
    auto result = limit_checked(zero_infinity, "x", SymbolicExpr::number(0));
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "zero times infinity is not a finite payload";
}

TEST(Calculus, FinitePayloadDomainContracts) {
    for (const char *source : {"0^-1", "ln(-1)"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "domain fixture parses";
        if (!expression) {
            continue;
        }
        auto invalid = limit_checked(expression.value(), "x", SymbolicExpr::number(0));
        EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::Inconclusive)) << "undefined constant cannot be certified finite";
    }
    auto enormous = SymbolicExpr::number(BigInt("1" + std::string(400, '0')));
    auto enormous_limit = limit_checked(enormous, "x", SymbolicExpr::number(0));
    auto finite_enormous = enormous_limit ? std::get_if<FiniteLimit>(&enormous_limit.value().value) : nullptr;
    EXPECT_TRUE((finite_enormous && finite_enormous->value->compare(enormous) == 0)) << "exact values beyond double range are still finite";
}

TEST(Calculus, LimitComputationErrors) {
    auto source = parse_expr("(x^2-1)/(x-1)");
    ASSERT_TRUE((source.has_value())) << "resource fixture parses";
    if (!source) {
        return;
    }
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext exhausted(limits);
    auto limited = limit_checked(source.value(), "x", SymbolicExpr::number(1), LimitDirection::Both, exhausted);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "nested one-sided work shares the caller budget";
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = limit_checked(source.value(), "x", SymbolicExpr::number(1), LimitDirection::Both, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "cancellation is not converted to an unknown value";
}

TEST(Calculus, ParameterDependentSideComparison) {
    auto a = SymbolicExpr::variable("a");
    auto zero = SymbolicExpr::number(0);
    auto condition = detail::make_node<RelationalNode>(
        detail::node(SymbolicExpr::variable("x")), detail::node(zero), RelationOp::LT);
    auto expression = detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{detail::node(a), condition}}, detail::node(zero)));
    for (const auto sign : {std::optional<Sign>{}, std::optional<Sign>{Sign::Zero},
                            std::optional<Sign>{Sign::NonZero}}) {
        auto assumptions = std::make_shared<AssumptionContext>();
        EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
        if (sign) {
            EXPECT_TRUE((assumptions->assume_sign_checked("a", *sign).has_value())) << "parameter sign is attached";
        }
        ComputationContext context;
        EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "use the supplied parameter facts";
        auto result = limit_checked(expression, "x", zero, LimitDirection::Both, context);
        if (!sign) {
            EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "piecewise(a, x<0, 0) is undecided when a is merely real";
        } else if (*sign == Sign::Zero) {
            auto finite = result ? std::get_if<FiniteLimit>(&result.value().value) : nullptr;
            EXPECT_TRUE((finite != nullptr)) << "equal sides under a=0 have a finite limit";
            if (finite) {
                auto value = evaluate_numeric(*finite->value, {{"a", 0.0}});
                EXPECT_TRUE((value && value.value().is_finite() && value.value().value == 0.0)) << "the finite limit under a=0 evaluates to zero";
            }
        } else {
            EXPECT_TRUE((result && std::holds_alternative<LimitDoesNotExist>(result.value().value))) << "a nonzero real parameter proves unequal sides";
        }
    }
}
