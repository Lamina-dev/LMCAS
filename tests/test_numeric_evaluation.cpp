#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "expr.hpp"
#include "internal/symbolic_ast.hpp"
#include <limits>

TEST(NumericEvaluation, TestEvaluatorPolicies) {
    auto cube_root = detail::make_expression_ptr(detail::make_node<PowerNode>(
        detail::node(SymbolicExpr::number(-8)),
        detail::node(SymbolicExpr::number(Rational(1, 3)))));
    EXPECT_FALSE(test_numeric_eval(cube_root).has_value());
    EXPECT_NEAR(test_numeric_value(cube_root), -2.0, 1e-12);

    auto scaled = detail::make_expression_ptr(detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::node(cube_root), detail::node(SymbolicExpr::number(2))}));
    EXPECT_FALSE(test_numeric_eval(scaled).has_value());
    EXPECT_NEAR(test_numeric_value(scaled), -4.0, 1e-12);

    auto reciprocal_zero = detail::make_expression_ptr(detail::make_node<PowerNode>(
        detail::node(SymbolicExpr::number(0)),
        detail::node(SymbolicExpr::number(-1))));
    EXPECT_FALSE(test_numeric_eval(reciprocal_zero).has_value());
    EXPECT_TRUE(std::isinf(test_numeric_value(reciprocal_zero)));

    auto variable = SymbolicExpr::variable("x");
    EXPECT_FALSE(test_numeric_eval(variable).has_value());
    EXPECT_TRUE(std::isnan(test_numeric_value(variable)));
    EXPECT_FALSE(test_numeric_eval(nullptr).has_value());
    EXPECT_DOUBLE_EQ(test_numeric_value(nullptr), 0.0);

    auto sinh = detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sinh,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::node(SymbolicExpr::number(1))}));
    auto finite_sinh = test_numeric_eval(sinh);
    ASSERT_TRUE(finite_sinh.has_value());
    EXPECT_NEAR(*finite_sinh, std::sinh(1.0), 1e-12);
    EXPECT_TRUE(std::isnan(test_numeric_value(sinh)));
}

TEST(NumericEvaluation, NumericBindingAndBudgets) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(x, SymbolicExpr::number(2));
    auto bound = evaluate_numeric(*expression, {{"x", 5.0}});
    ASSERT_TRUE((bound.has_value())) << "bound expression evaluates";
    if (bound) {
        EXPECT_NEAR(bound.value().value, 7.0, 1e-12) << "x + 2 at x=5";
    }

    auto unbound = evaluate_numeric(*expression);
    EXPECT_FALSE((unbound.has_value())) << "unbound x does not become zero";
    if (!unbound) {
        EXPECT_TRUE((unbound.error().code == CasErrc::UnboundSymbol)) << "unbound x reports UnboundSymbol";
    }

    auto invalid_log = SymbolicExpr::ln(SymbolicExpr::number(-1));
    auto log_result = evaluate_numeric(*invalid_log);
    EXPECT_FALSE((log_result.has_value())) << "ln(-1) is not a real number";
    if (!log_result) {
        EXPECT_TRUE((log_result.error().code == CasErrc::DomainError)) << "ln(-1) reports DomainError";
    }
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);
    auto exhausted = evaluate_numeric(*expression, {{"x", 1.0}}, context);
    EXPECT_FALSE((exhausted.has_value())) << "step budget stops recursive evaluation";
    if (!exhausted) {
        EXPECT_TRUE((exhausted.error().code == CasErrc::ResourceLimit)) << "budget exhaustion reports ResourceLimit";
    }
}

TEST(NumericEvaluation, CertifiedAbsoluteError) {
    auto atom = evaluate_numeric(*SymbolicExpr::number(0.5));
    ASSERT_TRUE(atom);
    EXPECT_DOUBLE_EQ(atom.value().absolute_error, 0.0);

    auto bound = evaluate_numeric(*SymbolicExpr::variable("x"), {{"x", 0.5}});
    ASSERT_TRUE(bound);
    EXPECT_DOUBLE_EQ(bound.value().absolute_error, 0.0);

    auto integer = evaluate_numeric(*SymbolicExpr::number(BigInt(1)));
    ASSERT_TRUE(integer);
    EXPECT_EQ(integer.value().absolute_error, std::numeric_limits<double>::infinity());
    auto fraction = evaluate_numeric(*SymbolicExpr::number(Rational(1, 3)));
    ASSERT_TRUE(fraction);
    EXPECT_EQ(fraction.value().absolute_error, std::numeric_limits<double>::infinity());

    auto expression = detail::make_expression_ptr(detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::node(SymbolicExpr::variable("x")),
            detail::node(SymbolicExpr::variable("y")),
            detail::node(SymbolicExpr::variable("z"))}));
    auto cancelled_terms = evaluate_numeric(
        *expression, {{"x", 1e16}, {"y", 1.0}, {"z", -1e16}});
    ASSERT_TRUE(cancelled_terms);
    EXPECT_DOUBLE_EQ(cancelled_terms.value().value, 0.0);
    EXPECT_EQ(cancelled_terms.value().status, NumericStatus::Finite);
    EXPECT_EQ(cancelled_terms.value().absolute_error, std::numeric_limits<double>::infinity());

    auto overflow = evaluate_numeric(
        *SymbolicExpr::add(SymbolicExpr::variable("x"), SymbolicExpr::variable("y")),
        {{"x", std::numeric_limits<double>::max()},
         {"y", std::numeric_limits<double>::max()}});
    EXPECT_TRUE((!overflow && overflow.error().code == CasErrc::NumericFailure));
    auto exponential_overflow = evaluate_numeric(
        *SymbolicExpr::exp(SymbolicExpr::number(1000.0)));
    EXPECT_TRUE((!exponential_overflow &&
                 exponential_overflow.error().code == CasErrc::NumericFailure));

    auto infinity = evaluate_numeric(*SymbolicExpr::infinity(1));
    ASSERT_TRUE(infinity);
    EXPECT_EQ(infinity.value().status, NumericStatus::PositiveInfinity);
    EXPECT_TRUE(std::isinf(infinity.value().value));
    auto negative_infinity = evaluate_numeric(*SymbolicExpr::infinity(-1));
    ASSERT_TRUE(negative_infinity);
    EXPECT_EQ(negative_infinity.value().status, NumericStatus::NegativeInfinity);
    EXPECT_TRUE(std::isinf(negative_infinity.value().value));

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto stopped = evaluate_numeric(*expression, {{"x", 1.0}, {"y", 2.0}, {"z", 3.0}}, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled));
}

TEST(NumericEvaluation, NumericConstantsAndRoots) {
    auto pi_constant = LMCAS::pi();
    ASSERT_TRUE(pi_constant);
    auto evaluated_pi = evaluate_numeric(*pi_constant.value());
    ASSERT_TRUE(evaluated_pi);
    EXPECT_EQ(evaluated_pi.value().absolute_error, std::numeric_limits<double>::infinity());
    EXPECT_NEAR(evaluated_pi.value().value, std::acos(-1.0), 1e-15);

    auto ascii_pi = evaluate_numeric(*SymbolicExpr::variable("pi"));
    auto unicode_pi = evaluate_numeric(*SymbolicExpr::variable("π"));
    EXPECT_TRUE((!ascii_pi && ascii_pi.error().code == CasErrc::UnboundSymbol));
    EXPECT_TRUE((!unicode_pi && unicode_pi.error().code == CasErrc::UnboundSymbol));
    auto bound_pi = evaluate_numeric(*SymbolicExpr::variable("pi"), {{"pi", 2.0}});
    ASSERT_TRUE(bound_pi);
    EXPECT_DOUBLE_EQ(bound_pi.value().value, 2.0);
    auto constant_with_binding = evaluate_numeric(*pi_constant.value(), {{"pi", 2.0}});
    ASSERT_TRUE(constant_with_binding);
    EXPECT_NEAR(constant_with_binding.value().value, std::acos(-1.0), 1e-15);

    auto root_variable = SymbolicExpr::variable("z");
    auto root_polynomial = SymbolicExpr::add(
        SymbolicExpr::multiply(root_variable, root_variable),
        SymbolicExpr::number(-2));
    auto root = evaluate_numeric(
        *SymbolicExpr::root_of(root_polynomial, "z", 0));
    ASSERT_TRUE((root.has_value())) << "real RootOf evaluates numerically";
    EXPECT_EQ(root.value().absolute_error, std::numeric_limits<double>::infinity());
    if (root) {
        {
            const double actual_value = (root.value().value);
            const double expected_value = (-std::sqrt(2.0));
            const double tolerance = (1e-12);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(NumericEvaluation, RealAxisComplexArgument) {
    auto x = SymbolicExpr::variable("x");
    auto argument_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ComplexArg,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(x)}));
    const struct {
        double input;
        double expected;
    } argument_cases[] = {
        {-4.0, std::acos(-1.0)},
        {4.0, 0.0},
        {-0.0, std::acos(-1.0)},
        {0.0, 0.0}};
    for (const auto &test : argument_cases) {
        auto argument = evaluate_numeric(*argument_expression, {{"x", test.input}});
        EXPECT_TRUE((argument && argument.value().is_finite())) << "bound real-axis argument has a finite principal value";
        if (argument) {
            {
                const double actual_value = (argument.value().value);
                const double expected_value = (test.expected);
                const double tolerance = (1e-14);
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
    auto imaginary_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ImagPart,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(x)}));
    auto imaginary = evaluate_numeric(*imaginary_expression, {{"x", -4.0}});
    EXPECT_TRUE((imaginary && imaginary.value().value == 0.0)) << "negative real values have zero imaginary part";
}

TEST(NumericEvaluation, ExactPowerExponentDomains) {
    auto x = SymbolicExpr::variable("x");
    const BigInt odd_exponent("9007199254740993");
    for (auto exponent : {SymbolicExpr::number(odd_exponent),
                          SymbolicExpr::number(Rational(odd_exponent, BigInt(1)))}) {
        auto odd_power = SymbolicExpr::power(x, exponent);
        auto negative = evaluate_numeric(*odd_power, {{"x", -1.0}});
        EXPECT_TRUE((negative && negative.value().value == -1.0)) << "an exact odd exponent keeps the sign of a negative base";
        auto signed_zero = evaluate_numeric(*odd_power, {{"x", -0.0}});
        EXPECT_TRUE((signed_zero && signed_zero.value().value == 0.0 &&
                     std::signbit(signed_zero.value().value)))
            << "an exact odd exponent preserves negative zero";
    }
    auto even_power = SymbolicExpr::power(
        x, SymbolicExpr::number(odd_exponent + BigInt(1)));
    auto even_value = evaluate_numeric(*even_power, {{"x", -1.0}});
    EXPECT_TRUE((even_value && even_value.value().value == 1.0)) << "an adjacent exact even exponent produces a positive result";

    auto fractional_power = SymbolicExpr::power(
        x, SymbolicExpr::number(Rational(odd_exponent, BigInt(2))));
    auto fractional_value = evaluate_numeric(*fractional_power, {{"x", -1.0}});
    EXPECT_TRUE((!fractional_value &&
                 fractional_value.error().code == CasErrc::DomainError))
        << "a fraction rounded to an integer double remains outside the real power domain";
}

TEST(NumericEvaluation, LambertNumericContract) {
    auto expression = SymbolicExpr::lambertw(SymbolicExpr::variable("x"));
    for (const auto &sample : std::vector<std::pair<double, double>>{
             {1e200, 454.39804503371403}, {1.0, 0.56714329040978384}, {1e-300, 1e-300}, {-1e-300, -1e-300}}) {
        const auto result = evaluate_numeric(*expression, {{"x", sample.first}});
        EXPECT_TRUE((result && result.value().is_finite())) << "Lambert W evaluation succeeds";
        if (result) {
            EXPECT_EQ(result.value().absolute_error, std::numeric_limits<double>::infinity());
            {
                const double actual_value = (result.value().value / sample.second);
                const double expected_value = (1.0);
                const double tolerance = (64.0 * std::numeric_limits<double>::epsilon());
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, tolerance);
            }
        }
    }
    const auto outside = evaluate_numeric(*expression, {{"x", -1.0}});
    EXPECT_TRUE((!outside && outside.error().code == CasErrc::DomainError)) << "finite out-of-domain Lambert input is a domain failure";
    const auto infinite = evaluate_numeric(
        *expression, {{"x", std::numeric_limits<double>::infinity()}});
    EXPECT_TRUE((!infinite && infinite.error().code == CasErrc::NumericFailure)) << "nonfinite numerical input is not a real-domain error";
}
