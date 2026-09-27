#include "test_common.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include <cmath>

using namespace LMCAS;

TEST(PolyUtils, IntegerCoefficientsRequireExactValues) {
    for (double value : {std::nextafter(1.0, 2.0), std::nextafter(1.0, 0.0)}) {
        auto expression = SymbolicExpr::number(value);
        auto integer = extract_coeff_value<BigInt>(expression);
        ASSERT_FALSE(integer);
        EXPECT_EQ(
            integer.error().code, CasErrc::UnsupportedExpression);
        auto polynomial = symbolic_to_poly<BigInt>(expression, "x");
        ASSERT_FALSE(polynomial);
        EXPECT_EQ(
            polynomial.error().code,
            CasErrc::UnsupportedExpression);
        auto rational = extract_coeff_value<Rational>(expression);
        ASSERT_TRUE(rational);
        EXPECT_EQ(rational.value(), Rational::from_double(value));
    }
    auto one = extract_coeff_value<BigInt>(SymbolicExpr::number(1.0));
    ASSERT_TRUE(one);
    EXPECT_EQ(one.value(), BigInt(1));
}

TEST(PolyUtils, SymbolicToPoly) {
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1));
        auto converted = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        ASSERT_TRUE(converted) << converted.error().message;
        const auto &poly = converted.value();
        EXPECT_EQ(
            poly.coeffs,
            (std::vector<BigInt>{BigInt(1), BigInt(2)}));
        EXPECT_TRUE(test_proved_equivalent(
            poly_to_symbolic(poly), expr));
    }

    {
        auto x = SymbolicExpr::variable("x");
        // x^2 - 3x + 2
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto neg3x = SymbolicExpr::multiply(SymbolicExpr::number(-3), SymbolicExpr::variable("x"));
        auto two = SymbolicExpr::number(2);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(x2, neg3x), two);
        auto converted = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        ASSERT_TRUE(converted) << converted.error().message;
        const auto &poly = converted.value();
        EXPECT_EQ(
            poly.coeffs,
            (std::vector<BigInt>{
                BigInt(2), BigInt(-3), BigInt(1)}));
        EXPECT_TRUE(test_proved_equivalent(
            poly_to_symbolic(poly), expr));
    }

    {
        auto x = SymbolicExpr::variable("x");
        // x^3 + x
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto expr = SymbolicExpr::add(x3, SymbolicExpr::variable("x"));
        auto converted = LMCAS::symbolic_to_poly<BigInt>(expr, "x");
        ASSERT_TRUE(converted) << converted.error().message;
        const auto &poly = converted.value();
        EXPECT_EQ(
            poly.coeffs,
            (std::vector<BigInt>{
                BigInt(0), BigInt(1), BigInt(0), BigInt(1)}));
        EXPECT_TRUE(test_proved_equivalent(
            poly_to_symbolic(poly), expr));
    }
}

TEST(PolyUtils, PolyToSymbolic) {
    {
        LMCAS::Polynomial<BigInt> poly(
            {BigInt(1), BigInt(2)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        auto expected = SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::number(2),
                SymbolicExpr::variable("x")),
            SymbolicExpr::number(1));
        ASSERT_NE(expr, nullptr);
        EXPECT_TRUE(test_proved_equivalent(expr, expected));
    }

    {
        LMCAS::Polynomial<BigInt> poly(
            {BigInt(2), BigInt(-3), BigInt(1)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        auto x = SymbolicExpr::variable("x");
        auto expected = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::add(
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-3), x),
                SymbolicExpr::number(2)));
        ASSERT_NE(expr, nullptr);
        EXPECT_TRUE(test_proved_equivalent(expr, expected));
    }

    {
        LMCAS::Polynomial<BigInt> poly({BigInt(5)}, "x");
        auto expr = LMCAS::poly_to_symbolic<BigInt>(poly);
        ASSERT_NE(expr, nullptr);
        EXPECT_TRUE(test_proved_equivalent(
            expr, SymbolicExpr::number(5)));
    }
}

TEST(PolyUtils, DependsOnVar) {
    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(1));
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE((result)) << "3x+1 depends on x";
    }

    {
        auto y = SymbolicExpr::variable("y");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), y),
            SymbolicExpr::number(5));
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE((!result)) << "2y+5 does not depend on x";
    }

    {
        auto expr = SymbolicExpr::number(42);
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE((!result)) << "constant 42 does not depend on x";
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::power(
            SymbolicExpr::add(x, SymbolicExpr::number(1)),
            SymbolicExpr::number(2));
        bool result = LMCAS::contains(*expr, "x");
        EXPECT_TRUE((result)) << "(x+1)^2 depends on x";
    }
}

TEST(PolyUtils, UnsupportedSubtreesAreAtomic) {
    auto x = SymbolicExpr::variable("x");
    const auto remainder = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
    const auto huge = SymbolicExpr::power(
        x, SymbolicExpr::number(BigInt("18446744073709551616")));
    const auto fractional = SymbolicExpr::power(x, SymbolicExpr::number(Rational(3, 2)));
    auto oversized = symbolic_to_poly<Rational>(SymbolicExpr::add(huge, remainder), "x");
    EXPECT_TRUE((!oversized && oversized.error().code == CasErrc::ResourceLimit)) << "oversized power rejects the whole sum as a resource limit";
    auto unsupported = symbolic_to_poly<Rational>(
        SymbolicExpr::multiply(SymbolicExpr::add(fractional, x), remainder), "x");
    EXPECT_TRUE((!unsupported && unsupported.error().code == CasErrc::UnsupportedExpression)) << "fractional nested factor rejects the whole product explicitly";
    const auto supported = symbolic_to_poly<Rational>(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(Rational(3))), remainder),
        "x");
    EXPECT_TRUE((supported && supported.value().degree() == 3 &&
                 supported.value().eval(Rational(2)) == Rational(11)))
        << "supported integer powers still convert with all sibling terms";
}

static void check_unsupported_polynomial_forms(const ExprPtr &x, const ExprPtr &a) {
    for (const auto &unsupported : {
             SymbolicExpr::sin(x), SymbolicExpr::power(x, SymbolicExpr::number(-1)),
             SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 2))),
             SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), a)}) {
        auto converted = symbolic_to_poly<Rational>(unsupported, "x");
        ASSERT_FALSE(converted);
        EXPECT_EQ(
            converted.error().code,
            CasErrc::UnsupportedExpression);
}
}

static void check_symbolic_coefficient_conversion(const ExprPtr &x, const ExprPtr &a) {
    auto parameter_poly = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), a);
    auto symbolic =
        symbolic_to_poly<SymbolicPolyCoeff>(parameter_poly, "x");
    ASSERT_TRUE(symbolic) << symbolic.error().message;
    EXPECT_EQ(symbolic.value().degree(), 2);
    const auto reconstructed = poly_to_symbolic(symbolic.value());
    const auto evaluated = test_numeric_eval(
        reconstructed->substitute("a", SymbolicExpr::number(7))
            ->substitute("x", SymbolicExpr::number(3)));
    ASSERT_TRUE(evaluated.has_value());
    EXPECT_TRUE(std::isfinite(*evaluated));
    EXPECT_EQ(*evaluated, 16.0);
}

static void check_unsupported_rational_coefficients(const ExprPtr &a) {
    for (const auto &coefficient : {a, SymbolicExpr::sin(a),
                                    SymbolicExpr::power(SymbolicExpr::number(2), SymbolicExpr::number(Rational(1, 2)))}) {
        auto value = extract_coeff_value<Rational>(coefficient);
        ASSERT_FALSE(value);
        EXPECT_EQ(
            value.error().code, CasErrc::UnsupportedExpression);
    }
}

static void check_exact_coefficient_and_power_conversion(const ExprPtr &x) {
    auto binary =
        extract_coeff_value<Rational>(SymbolicExpr::number(0.1));
    ASSERT_TRUE(binary) << binary.error().message;
    EXPECT_EQ(binary.value(), Rational::from_double(0.1));
    EXPECT_NE(binary.value(), Rational(1, 10));
    auto accepted = symbolic_to_poly<Rational>(
        SymbolicExpr::power(x, SymbolicExpr::number(999)), "x");
    ASSERT_TRUE(accepted) << accepted.error().message;
    EXPECT_EQ(accepted.value().degree(), 999);
    EXPECT_EQ(accepted.value().eval(Rational(1)), Rational(1));
}

static void check_conversion_failures_and_zero_powers(const ExprPtr &x) {
    auto rejected = symbolic_to_poly<Rational>(
        SymbolicExpr::power(x, SymbolicExpr::number(1000)), "x");
    ASSERT_FALSE(rejected);
    EXPECT_EQ(rejected.error().code, CasErrc::ResourceLimit);
    auto null_input = symbolic_to_poly<Rational>(nullptr, "x");
    auto empty_variable = symbolic_to_poly<Rational>(x, "");
    ASSERT_FALSE(null_input);
    EXPECT_EQ(null_input.error().code, CasErrc::InvalidArgument);
    ASSERT_FALSE(empty_variable);
    EXPECT_EQ(
        empty_variable.error().code, CasErrc::InvalidArgument);
    auto unknown_zero_power = symbolic_to_poly<Rational>(
        SymbolicExpr::power(x, SymbolicExpr::number(0)), "x");
    auto defined_zero_power = symbolic_to_poly<Rational>(
        SymbolicExpr::power(
            SymbolicExpr::number(2), SymbolicExpr::number(0)),
        "x");
    auto undefined_zero_power = symbolic_to_poly<Rational>(
        SymbolicExpr::power(
            SymbolicExpr::number(0), SymbolicExpr::number(0)),
        "x");
    EXPECT_FALSE(unknown_zero_power);
    EXPECT_FALSE(undefined_zero_power);
    ASSERT_TRUE(defined_zero_power)
        << defined_zero_power.error().message;
    EXPECT_EQ(
        defined_zero_power.value().eval(Rational(0)), Rational(1));
}

TEST(PolyUtils, ConversionResultContract) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto zero = symbolic_to_poly<Rational>(
        SymbolicExpr::number(0), "x");
    ASSERT_TRUE(zero) << zero.error().message;
    EXPECT_TRUE(zero.value().is_zero());
    check_unsupported_polynomial_forms(x, a);
    check_symbolic_coefficient_conversion(x, a);
    check_unsupported_rational_coefficients(a);
    check_exact_coefficient_and_power_conversion(x);
    check_conversion_failures_and_zero_powers(x);
}

TEST(PolyUtils, AffineRecognition) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    struct Case {
        ExprPtr expression;
        ExprPtr slope;
        ExprPtr offset;
    };
    for (const auto &sample : std::vector<Case>{
             {SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x),
                                SymbolicExpr::number(3)),
              SymbolicExpr::number(2), SymbolicExpr::number(3)},
             {SymbolicExpr::add(SymbolicExpr::divide(x, SymbolicExpr::number(2)), b),
              SymbolicExpr::number(Rational(1, 2)), b},
             {SymbolicExpr::add(SymbolicExpr::multiply(a, x), b), a, b},
             {SymbolicExpr::add(
                  SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::number(1)),
                                      SymbolicExpr::number(2)),
                  SymbolicExpr::multiply(SymbolicExpr::number(-1),
                      SymbolicExpr::power(x, SymbolicExpr::number(2)))),
              SymbolicExpr::number(2), SymbolicExpr::number(1)},
             {SymbolicExpr::number(0), SymbolicExpr::number(0), SymbolicExpr::number(0)},
             {b, SymbolicExpr::number(0), b}}) {
        ComputationContext context;
        auto affine = detail::recognize_affine(*sample.expression, "x", context);
        EXPECT_TRUE((affine && affine.value())) << "affine expression matches";
        if (affine && affine.value()) {
            EXPECT_TRUE((detail::node(affine.value()->slope->simplify())->equals(*detail::node(sample.slope)))) << "affine slope is preserved";
            EXPECT_TRUE((detail::node(affine.value()->offset->simplify())->equals(*detail::node(sample.offset)))) << "affine offset is preserved";
        }
    }
    for (const auto &expression : {SymbolicExpr::multiply(x, x), SymbolicExpr::sin(x)}) {
        ComputationContext context;
        auto affine = detail::recognize_affine(*expression, "x", context);
        EXPECT_TRUE((affine && !affine.value())) << "non-affine expression is a successful non-match";
    }
    ComputationContext context;
    auto invalid = detail::recognize_affine(*x, "", context);
    EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::InvalidArgument)) << "empty variable is invalid";
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted_context(limits);
    auto exhausted = detail::recognize_affine(*x, "x", exhausted_context);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "affine recognition preserves the caller's exhausted budget";
}

TEST(PolyUtils, AffineRecognitionSharesConversionLimits) {
    auto x = SymbolicExpr::variable("x");
    auto power = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(1)), SymbolicExpr::number(8));
    ComputationContext ordinary_context;
    auto ordinary = detail::recognize_affine(*power, "x", ordinary_context);
    ASSERT_TRUE(ordinary);
    EXPECT_FALSE(ordinary.value());

    struct LimitCase {
        const char* name;
        std::size_t ResourceLimits::* limit;
        std::size_t value;
    };
    for (const auto& sample : {
             LimitCase{"steps", &ResourceLimits::max_steps, 1},
             LimitCase{"nodes", &ResourceLimits::max_ast_nodes, 0},
             LimitCase{"depth", &ResourceLimits::max_recursion_depth, 0},
             LimitCase{"terms", &ResourceLimits::max_expansion_terms, 0},
             LimitCase{"expanded nodes", &ResourceLimits::max_ast_nodes, 8},
             LimitCase{"expanded terms", &ResourceLimits::max_expansion_terms, 8}}) {
        SCOPED_TRACE(sample.name);
        ResourceLimits limits;
        limits.*sample.limit = sample.value;
        ComputationContext context(limits);
        auto result = detail::recognize_affine(*power, "x", context);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().code, CasErrc::ResourceLimit);
    }

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = detail::recognize_affine(*power, "x", cancelled_context);
    ASSERT_FALSE(cancelled);
    EXPECT_EQ(cancelled.error().code, CasErrc::Cancelled);
}

TEST(PolyUtils, NestedPolynomialPowerChecksFinalCoefficientCount) {
    auto x = SymbolicExpr::variable("x");
    auto base = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto inner = detail::make_node<PowerNode>(
        detail::node(base), detail::node(SymbolicExpr::number(20)));
    auto nested = detail::make_expression_ptr(detail::make_node<PowerNode>(
        inner, detail::node(SymbolicExpr::number(20))));
    ResourceLimits limits;
    limits.max_expansion_terms = 100;
    ComputationContext context(limits);
    auto result = detail::recognize_affine(*nested, "x", context);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::ResourceLimit);
}

TEST(PolyUtils, PolynomialCoefficientsShareAggregateNodeLimit) {
    auto expression = SymbolicExpr::power(
        SymbolicExpr::add(SymbolicExpr::variable("a"),
            SymbolicExpr::multiply(SymbolicExpr::variable("b"), SymbolicExpr::variable("x"))),
        SymbolicExpr::number(2));
    ResourceLimits limits;
    limits.max_ast_nodes = 9;
    ComputationContext context(limits);
    auto result = detail::symbolic_to_poly_checked(*expression, "x", context);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::ResourceLimit);

    ComputationContext ordinary_context;
    auto ordinary = detail::symbolic_to_poly_checked(*expression, "x", ordinary_context);
    ASSERT_TRUE(ordinary);
    EXPECT_TRUE(test_proved_equivalent(poly_to_symbolic(ordinary.value()), expression));
}
