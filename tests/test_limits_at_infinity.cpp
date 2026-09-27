#include "limit_result.hpp"
#include "test_common.hpp"
#include "expr.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

TEST(LimitsAtInfinity, FiniteTailLimits) {
    struct Sample {
        const char *expression;
        int tail;
        Rational expected;
    };
    const Sample samples[] = {
        {"x/x^2", 1, Rational(0)},
        {"(3*x^2+x)/(2*x^2+1)", 1, Rational(3, 2)},
        {"-2*x^2/x^2", 1, Rational(-2)},
        {"x^2/exp(x)", 1, Rational(0)},
        {"x^10/exp(x)", 1, Rational(0)},
        {"ln(x)/x", 1, Rational(0)},
        {"exp(-x)", 1, Rational(0)},
        {"x^2/(x^2+1)", -1, Rational(1)},
        {"exp(x)", -1, Rational(0)},
        {"1/x", -1, Rational(0)},
        {"sin(1/x)", 1, Rational(0)},
        {"sin(1/x)", -1, Rational(0)},
        {"cos(1/x)", 1, Rational(1)},
        {"cos(1/x)", -1, Rational(1)},
        {"x*exp(-ln(x))", 1, Rational(1)},
        {"x*exp(-x)", 1, Rational(0)},
        {"x^10*exp(-x)", 1, Rational(0)},
        {"x*exp(-sqrt(x))", 1, Rational(0)}};
    for (const auto &sample : samples) {
        auto expression = parse_expr(sample.expression);
        ASSERT_TRUE((expression.has_value())) << "tail fixture parses";
        if (!expression)
            continue;
        auto result = limit_checked(expression.value(), "x", SymbolicExpr::infinity(sample.tail));
        const auto *finite = result ? std::get_if<FiniteLimit>(&result.value().value) : nullptr;
        EXPECT_TRUE((finite != nullptr)) << sample.expression;
        if (!finite)
            continue;
        ComputationContext context;
        auto proof = check_equivalent(finite->value, SymbolicExpr::number(sample.expected), context);
        EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << sample.expression;
        auto numeric = evaluate_numeric(*finite->value);
        EXPECT_TRUE((numeric && numeric.value().is_finite())) << "finite tail payload evaluates finitely";
    }
}

TEST(LimitsAtInfinity, InfiniteAndOscillatoryTails) {
    auto polynomial = parse_expr("x^3/x");
    ASSERT_TRUE((polynomial.has_value())) << "polynomial fixture parses";
    if (polynomial) {
        auto value = limit_checked(polynomial.value(), "x", SymbolicExpr::infinity());
        EXPECT_TRUE((value && std::holds_alternative<PositiveInfinityLimit>(value.value().value))) << "positive quadratic tail diverges positively";
    }
    for (int sign : {-1, 1}) {
        for (const char *source : {"sin(x)", "cos(x)", "sin(x^3)"}) {
            auto expression = parse_expr(source);
            ASSERT_TRUE((expression.has_value())) << "oscillation fixture parses";
            if (!expression)
                continue;
            auto value = limit_checked(expression.value(), "x", SymbolicExpr::infinity(sign));
            EXPECT_TRUE((value && std::holds_alternative<LimitDoesNotExist>(value.value().value))) << "continuous argument escaping to infinity produces oscillation";
        }
    }
    auto invalid = parse_expr("ln(-x)");
    ASSERT_TRUE((invalid.has_value())) << "invalid-tail fixture parses";
    if (invalid) {
        auto value = limit_checked(invalid.value(), "x", SymbolicExpr::infinity());
        EXPECT_TRUE((!value && value.error().code == CasErrc::Inconclusive)) << "logarithm of a negative real tail is not a finite or infinite real limit";
    }
}
