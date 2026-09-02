#include "../include/symbolic.hpp"
#include "../include/solve_strategies.hpp"
#include "test_common.hpp"
#include "../src/internal/exact_algebraic.hpp"
#include "../include/integration.hpp"
#include "../include/integrator.hpp"
#include "../include/residual_verification.hpp"
#include "../include/lsr_expr.hpp"
#include "../include/quantity.hpp"

int main() {
    TEST_CASE("Result propagation preserves lvalue lifetime");
    auto mutable_result =
        lamina::Result<std::vector<int>>::success({1, 2, 3});
    auto& mutable_value =
        lamina::detail::propagate_result(mutable_result);
    EXPECT_TRUE(&mutable_value == &mutable_result.value(),
                "mutable propagation references the stored value");
    mutable_value[0] = 7;
    EXPECT_TRUE(mutable_result.value()[0] == 7,
                "mutable propagation does not reference a temporary copy");

    const auto const_result =
        lamina::Result<std::vector<int>>::success({4, 5});
    const auto& const_value =
        lamina::detail::propagate_result(const_result);
    EXPECT_TRUE(&const_value == &const_result.value(),
                "const propagation references the stored value");

    auto moved_value = lamina::detail::propagate_result(
        lamina::Result<std::unique_ptr<int>>::success(
            std::make_unique<int>(9)));
    EXPECT_TRUE(moved_value && *moved_value == 9,
                "rvalue propagation moves move-only values");

    TEST_CASE("SolutionSet alternatives are closed and state-specific");

    lamina::SolutionSet empty = lamina::EmptySolutions{};
    EXPECT_TRUE(std::holds_alternative<lamina::EmptySolutions>(empty),
                "empty solution set has one explicit alternative");

    auto x = SymbolicExpr::variable("x");
    lamina::SolutionSet finite = lamina::FiniteSolutions{{
        lamina::FiniteSolution{x, 2, {}}
    }};
    const auto* finite_values =
        std::get_if<lamina::FiniteSolutions>(&finite);
    EXPECT_TRUE(finite_values && finite_values->values.size() == 1,
                "finite solution set stores roots");
    EXPECT_TRUE(finite_values &&
                    finite_values->values[0].multiplicity == 2,
                "finite solution set preserves multiplicity");

    auto k = SymbolicExpr::variable("k");
    lamina::SolutionSet parametric = lamina::ParametricSolutions{{
        lamina::ParametricSolution{
            SymbolicExpr::multiply(k, SymbolicExpr::variable("pi")),
            {"k"},
            {}
        }
    }};
    const auto* parametric_values =
        std::get_if<lamina::ParametricSolutions>(&parametric);
    EXPECT_TRUE(parametric_values &&
                    parametric_values->values.size() == 1 &&
                    parametric_values->values[0].integer_parameters ==
                        std::vector<std::string>{"k"},
                "parametric solution set stores integer parameters");

    TEST_CASE("Integral outcomes cannot contain inactive fields");
    lamina::IntegralOutcome closed = lamina::Verified<lamina::ClosedFormIntegral>{
        lamina::ClosedFormIntegral{
            x, {SymbolicExpr::variable("x_nonzero")}},
        lamina::ByConstructionProof{}};
    EXPECT_TRUE(
        std::holds_alternative<
            lamina::Verified<lamina::ClosedFormIntegral>>(closed),
        "closed integral carries a proof certificate");
    lamina::IntegralOutcome unevaluated = lamina::UnevaluatedIntegral{x};
    EXPECT_TRUE(
        std::holds_alternative<lamina::UnevaluatedIntegral>(unevaluated),
        "unevaluated integral has no verification fields");

    TEST_CASE("Residual verification distinguishes proof from uncertainty");
    lamina::ComputationContext residual_context;
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
        SymbolicExpr::number(1));
    auto factored_square = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        SymbolicExpr::number(2));
    auto equivalent = lamina::check_equivalent(
        expanded_square, factored_square, residual_context);
    EXPECT_TRUE(equivalent &&
                    std::holds_alternative<lamina::ProvedZeroResidual>(
                        equivalent.value()),
                "expanded polynomial identity has an exact certificate");
    auto nonzero = lamina::check_zero_residual(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        residual_context);
    EXPECT_TRUE(nonzero &&
                    std::holds_alternative<lamina::ProvedNonzeroResidual>(
                        nonzero.value()),
                "nonzero rational polynomial is exactly refuted");
    auto unproved = lamina::check_zero_residual(
        SymbolicExpr::sin(x), residual_context);
    EXPECT_TRUE(unproved &&
                    std::holds_alternative<lamina::UnprovedResidual>(
                        unproved.value()),
                "unsupported transcendental residual remains unproved");

    TEST_CASE("Checked solve dispatcher distinguishes mathematical outcomes");

    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    lamina::ComputationContext exact_context;
    auto exact = lamina::solve_equation(
        x2_plus_one, "x", exact_context, lamina::SolveOptions{});
    const auto* exact_finite = exact
        ? std::get_if<lamina::FiniteSolutions>(&exact.value()) : nullptr;
    EXPECT_TRUE(exact_finite && exact_finite->values.size() == 2,
                "exact polynomial preserves algebraic root count");
    if (exact_finite && exact_finite->values.size() == 2) {
        EXPECT_TRUE(
            !exact_finite->values[0].value->to_string().empty(),
            "exact candidates remain symbolic expressions");
    }

    lamina::ComputationContext empty_context;
    auto no_solution = lamina::solve_equation(
        SymbolicExpr::number(1), "x", empty_context, lamina::SolveOptions{});
    EXPECT_TRUE(no_solution &&
                    std::holds_alternative<lamina::EmptySolutions>(
                        no_solution.value()),
                "nonzero constant equation has an empty solution set");

    lamina::ComputationContext universal_context;
    auto universal = lamina::solve_equation(
        SymbolicExpr::number(0), "x", universal_context, lamina::SolveOptions{});
    EXPECT_TRUE(universal &&
                    std::holds_alternative<lamina::UniversalSolutions>(
                        universal.value()),
                "zero equation has the universal solution set");

    auto periodic = lamina::solve_equation(
        SymbolicExpr::sin(x), "x", lamina::SolveOptions{});
    const auto* periodic_values = periodic
        ? std::get_if<lamina::ParametricSolutions>(&periodic.value())
        : nullptr;
    EXPECT_TRUE(periodic_values && periodic_values->values.size() == 1,
                "sin(x)=0 returns a complete integer-parameter family");
    if (periodic_values) {
        EXPECT_TRUE(
            periodic_values->values[0].value->to_string().find("_k") !=
                std::string::npos,
            "periodic family preserves its integer parameter");
    }
    auto sine_equation = SymbolicExpr::add(
        SymbolicExpr::sin(x), SymbolicExpr::number(Rational(-1, 2)));
    auto unsupported_equation = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<std::shared_ptr<const SymbolicNode>>{
                     lamina::detail::node(x)}));
    lamina::ComputationContext unsupported_context;
    auto unsupported = lamina::solve_equation(
        unsupported_equation, "x", unsupported_context, lamina::SolveOptions{});
    EXPECT_TRUE(!unsupported &&
                    unsupported.error().code == lamina::CasErrc::Inconclusive,
                "unsupported symbolic equation returns outer Inconclusive");

    auto default_unsupported = lamina::solve_equation(
        unsupported_equation, "x", lamina::SolveOptions{});
    EXPECT_TRUE(!default_unsupported &&
                    default_unsupported.error().code ==
                        lamina::CasErrc::Inconclusive,
                "default-context solver preserves outer Inconclusive");

    lamina::SolveOptions numeric_options;
    numeric_options.allow_numeric = true;
    numeric_options.has_initial_guess = true;
    numeric_options.initial_guess = 0.5;
    lamina::ComputationContext numeric_context;
    auto numeric = lamina::solve_equation(
        sine_equation, "x", numeric_context, numeric_options);
    const auto* numeric_finite = numeric
        ? std::get_if<lamina::FiniteSolutions>(&numeric.value()) : nullptr;
    EXPECT_TRUE(numeric_finite && !numeric_finite->values.empty(),
                "explicit numeric solving returns verified finite candidates");

    TEST_CASE("Checked solve dispatcher preserves computation errors");

    auto symbolic_parameter_equation =
        SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    lamina::ComputationContext symbolic_parameter_context;
    auto symbolic_parameter = lamina::solve_equation(
        symbolic_parameter_equation, "x",
        symbolic_parameter_context, numeric_options);
    EXPECT_TRUE(symbolic_parameter &&
                    std::holds_alternative<lamina::FiniteSolutions>(
                        symbolic_parameter.value()),
                "symbolic linear coefficients remain exact");

    lamina::CancellationToken cancellation;
    cancellation.cancel();
    lamina::ComputationContext cancelled_context({}, cancellation);
    auto cancelled = lamina::solve_equation(
        x2_plus_one, "x", cancelled_context, lamina::SolveOptions{});
    EXPECT_TRUE(!cancelled && cancelled.error().code == lamina::CasErrc::Cancelled,
                "dispatcher observes cancellation");

    lamina::ResourceLimits limits;
    limits.max_steps = 1;
    lamina::ComputationContext limited_context(limits);
    auto limited = lamina::solve_equation(
        x2_plus_one, "x", limited_context, lamina::SolveOptions{});
    EXPECT_TRUE(!limited && limited.error().code == lamina::CasErrc::ResourceLimit,
                "dispatcher enforces shared step budgets");

    TEST_CASE("Integrator preserves cancellation and recursive budgets");
    lamina::Integrator checked_integrator;
    lamina::CancellationToken integration_cancellation;
    integration_cancellation.cancel();
    lamina::ComputationContext cancelled_integration_context(
        {}, integration_cancellation);
    auto cancelled_integral = checked_integrator.integrate_checked(
        *x2_plus_one, "x", cancelled_integration_context);
    EXPECT_TRUE(!cancelled_integral &&
                    cancelled_integral.error().code ==
                        lamina::CasErrc::Cancelled,
                "integration preserves cancellation");

    lamina::ResourceLimits integration_limits;
    integration_limits.max_steps = 2;
    lamina::ComputationContext limited_integration_context(
        integration_limits);
    auto recursive_integrand = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(x));
    auto limited_integral = checked_integrator.integrate_checked(
        *recursive_integrand, "x", limited_integration_context);
    EXPECT_TRUE(!limited_integral &&
                    limited_integral.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "integration preserves recursive step budgets");



    TEST_CASE("LSR number domains classify exact algebraic expressions");
    auto sqrt_two_expression = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    auto sqrt_two_real = lamina::lsr::domain_contains(
        lamina::lsr::reals(), sqrt_two_expression);
    auto sqrt_two_rational = lamina::lsr::domain_contains(
        lamina::lsr::rationals(), sqrt_two_expression);
    EXPECT_TRUE(sqrt_two_real && sqrt_two_real.value(),
                "sqrt(2) is proven real");
    EXPECT_TRUE(sqrt_two_rational && !sqrt_two_rational.value(),
                "sqrt(2) is proven non-rational");

    TEST_CASE("Quantity powers distinguish dimensioned and dimensionless bases");
    lamina::ComputationContext quantity_context;
    auto dimensionless_power = lamina::quantity_power(
        SymbolicExpr::number(2), SymbolicExpr::number(0.5),
        quantity_context);
    EXPECT_TRUE(dimensionless_power.has_value(),
                "dimensionless quantities accept approximate real exponents");
    auto metre = lamina::attach_unit(
        SymbolicExpr::number(1), "m", quantity_context);
    EXPECT_TRUE(metre.has_value(), "metre quantity construction succeeds");
    if (metre) {
        auto approximate_dimensioned = lamina::quantity_power(
            metre.value(), SymbolicExpr::number(2.0), quantity_context);
        EXPECT_TRUE(!approximate_dimensioned &&
                        approximate_dimensioned.error().code ==
                            lamina::CasErrc::UnitInvalid,
                    "dimensioned quantities require exact rational exponents");
    }

    TEST_CASE("Certified real algebraic values compare without floating point");
    lamina::ComputationContext algebraic_context;
    auto sqrt_two = lamina::detail::make_exact_real_algebraic(
        lamina::Polynomial<Rational>(
            {Rational(-2), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    auto sqrt_three = lamina::detail::make_exact_real_algebraic(
        lamina::Polynomial<Rational>(
            {Rational(-3), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    EXPECT_TRUE(sqrt_two.has_value() && sqrt_three.has_value(),
                "positive square roots have certified isolating intervals");
    if (sqrt_two && sqrt_three) {
        auto ordering = lamina::detail::compare_exact_real_algebraic(
            sqrt_two.value(), sqrt_three.value(), algebraic_context);
        EXPECT_TRUE(ordering.has_value() && ordering.value() < 0,
                    "sqrt(2) is certified less than sqrt(3)");
        auto equality = lamina::detail::equal_exact_real_algebraic(
            sqrt_two.value(), sqrt_two.value(), algebraic_context);
        EXPECT_TRUE(equality.has_value() && equality.value(),
                    "same certified root compares equal");
    }

    TEST_CASE("IntegralNode binds variables and supports calculus");
    auto y = SymbolicExpr::variable("y");
    auto body = SymbolicExpr::multiply(x, y);
    auto integral_expression = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<IntegralNode>(
            lamina::detail::node(body), "x"));
    EXPECT_TRUE(integral_expression->to_string() == "Integral(x*y, x)" ||
                    integral_expression->to_string() == "Integral(y*x, x)",
                "integral has a first-class printed form");
    auto parsed_integral = lamina::lsr::parse_expr(
        integral_expression->to_string());
    EXPECT_TRUE(parsed_integral.has_value() &&
                    lamina::detail::node(parsed_integral.value())->equals(
                        *lamina::detail::node(integral_expression)),
                "integral print/parse round-trip preserves binder structure");
    EXPECT_TRUE(!lamina::Integrator::depends_on(*integral_expression, "x") &&
                    lamina::Integrator::depends_on(*integral_expression, "y"),
                "integral variable is bound in free-variable analysis");
    auto capture_avoiding = integral_expression->substitute("y", x);
    EXPECT_TRUE(lamina::Integrator::depends_on(*capture_avoiding, "x"),
                "substitution does not capture a replacement variable");
    auto primitive = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<IntegralNode>(
            lamina::detail::node(x), "x"));
    auto derivative = primitive->differentiate("x");
    EXPECT_TRUE(derivative && derivative->to_string() == "x",
                "d/dx Integral(x,x) returns the integrand");

    return TEST_REPORT();
}
