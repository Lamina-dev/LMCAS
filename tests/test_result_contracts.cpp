#include "../include/symbolic.hpp"
#include "../include/solve_strategies.hpp"
#include "test_common.hpp"
#include "../src/internal/exact_algebraic.hpp"
#include "../include/integration.hpp"
#include "../include/integrator.hpp"
#include "../include/residual_verification.hpp"
#include "../include/expr.hpp"
#include "../include/quantity.hpp"

using namespace LMCAS;

int main() {
    TEST_CASE("Result carries values without exception propagation");
    auto make_values = []() -> LMCAS::Result<std::vector<int>> {
        return std::vector<int>{1, 2, 3};
    };
    auto mutable_result = make_values();
    auto& mutable_value = mutable_result.value();
    mutable_value[0] = 7;
    EXPECT_TRUE(mutable_result.value()[0] == 7,
                "value access preserves the stored object");

    const auto const_result =
        LMCAS::Result<std::vector<int>>::success({4, 5});
    const auto& const_value = const_result.value();
    EXPECT_TRUE(&const_value == &const_result.value(),
                "const value access preserves the stored object");

    auto make_pointer =
        []() -> LMCAS::Result<std::unique_ptr<int>> {
        return std::make_unique<int>(9);
    };
    auto pointer_result = make_pointer();
    auto moved_value = std::move(pointer_result.value());
    EXPECT_TRUE(moved_value && *moved_value == 9,
                "direct success returns support move-only values");

    TEST_CASE("SolutionSet alternatives are closed and state-specific");

    LMCAS::SolutionSet empty = LMCAS::EmptySolutions{};
    EXPECT_TRUE(std::holds_alternative<LMCAS::EmptySolutions>(empty),
                "empty solution set has one explicit alternative");

    auto x = SymbolicExpr::variable("x");
    LMCAS::SolutionSet finite = LMCAS::FiniteSolutions{{
        LMCAS::FiniteSolution{x, 2, {}}
    }};
    const auto* finite_values =
        std::get_if<LMCAS::FiniteSolutions>(&finite);
    EXPECT_TRUE(finite_values && finite_values->values.size() == 1,
                "finite solution set stores roots");
    EXPECT_TRUE(finite_values &&
                    finite_values->values[0].multiplicity == 2,
                "finite solution set preserves multiplicity");

    auto k = SymbolicExpr::variable("k");
    LMCAS::SolutionSet parametric = LMCAS::ParametricSolutions{{
        LMCAS::ParametricSolution{
            SymbolicExpr::multiply(k, SymbolicExpr::variable("pi")),
            {"k"},
            {}
        }
    }};
    const auto* parametric_values =
        std::get_if<LMCAS::ParametricSolutions>(&parametric);
    EXPECT_TRUE(parametric_values &&
                    parametric_values->values.size() == 1 &&
                    parametric_values->values[0].integer_parameters ==
                        std::vector<std::string>{"k"},
                "parametric solution set stores integer parameters");

    TEST_CASE("Integral outcomes cannot contain inactive fields");
    LMCAS::IntegralOutcome closed = LMCAS::Verified<LMCAS::ClosedFormIntegral>{
        LMCAS::ClosedFormIntegral{
            x, {SymbolicExpr::variable("x_nonzero")}},
        LMCAS::ByConstructionProof{}};
    EXPECT_TRUE(
        std::holds_alternative<
            LMCAS::Verified<LMCAS::ClosedFormIntegral>>(closed),
        "closed integral carries a proof certificate");
    LMCAS::IntegralOutcome unevaluated = LMCAS::UnevaluatedIntegral{x};
    EXPECT_TRUE(
        std::holds_alternative<LMCAS::UnevaluatedIntegral>(unevaluated),
        "unevaluated integral has no verification fields");

    TEST_CASE("Residual verification distinguishes proof from uncertainty");
    LMCAS::ComputationContext residual_context;
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
        SymbolicExpr::number(1));
    auto factored_square = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        SymbolicExpr::number(2));
    auto equivalent = LMCAS::check_equivalent(
        expanded_square, factored_square, residual_context);
    EXPECT_TRUE(equivalent &&
                    std::holds_alternative<LMCAS::ProvedZeroResidual>(
                        equivalent.value()),
                "expanded polynomial identity has an exact certificate");
    auto nonzero = LMCAS::check_zero_residual(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        residual_context);
    EXPECT_TRUE(nonzero &&
                    std::holds_alternative<LMCAS::ProvedNonzeroResidual>(
                        nonzero.value()),
                "nonzero rational polynomial is exactly refuted");
    auto unproved = LMCAS::check_zero_residual(
        SymbolicExpr::sin(x), residual_context);
    EXPECT_TRUE(unproved &&
                    std::holds_alternative<LMCAS::UnprovedResidual>(
                        unproved.value()),
                "unsupported transcendental residual remains unproved");

    TEST_CASE("Checked solve dispatcher distinguishes mathematical outcomes");

    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    LMCAS::ComputationContext exact_context;
    auto exact = LMCAS::solve_equation(
        x2_plus_one, "x", exact_context, LMCAS::SolveOptions{});
    const auto* exact_finite = exact
        ? std::get_if<LMCAS::FiniteSolutions>(&exact.value()) : nullptr;
    EXPECT_TRUE(exact_finite && exact_finite->values.size() == 2,
                "exact polynomial preserves algebraic root count");
    if (exact_finite && exact_finite->values.size() == 2) {
        EXPECT_TRUE(
            !exact_finite->values[0].value->to_string().empty(),
            "exact candidates remain symbolic expressions");
    }

    LMCAS::ComputationContext empty_context;
    auto no_solution = LMCAS::solve_equation(
        SymbolicExpr::number(1), "x", empty_context, LMCAS::SolveOptions{});
    EXPECT_TRUE(no_solution &&
                    std::holds_alternative<LMCAS::EmptySolutions>(
                        no_solution.value()),
                "nonzero constant equation has an empty solution set");

    LMCAS::ComputationContext universal_context;
    auto universal = LMCAS::solve_equation(
        SymbolicExpr::number(0), "x", universal_context, LMCAS::SolveOptions{});
    EXPECT_TRUE(universal &&
                    std::holds_alternative<LMCAS::UniversalSolutions>(
                        universal.value()),
                "zero equation has the universal solution set");

    auto periodic = LMCAS::solve_equation(
        SymbolicExpr::sin(x), "x", LMCAS::SolveOptions{});
    const auto* periodic_values = periodic
        ? std::get_if<LMCAS::ParametricSolutions>(&periodic.value())
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
    auto unsupported_equation = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<std::shared_ptr<const SymbolicNode>>{
                     LMCAS::detail::node(x)}));
    LMCAS::ComputationContext unsupported_context;
    auto unsupported = LMCAS::solve_equation(
        unsupported_equation, "x", unsupported_context, LMCAS::SolveOptions{});
    EXPECT_TRUE(!unsupported &&
                    unsupported.error().code == LMCAS::CasErrc::Inconclusive,
                "unsupported symbolic equation returns outer Inconclusive");

    auto default_unsupported = LMCAS::solve_equation(
        unsupported_equation, "x", LMCAS::SolveOptions{});
    EXPECT_TRUE(!default_unsupported &&
                    default_unsupported.error().code ==
                        LMCAS::CasErrc::Inconclusive,
                "default-context solver preserves outer Inconclusive");

    LMCAS::SolveOptions numeric_options;
    numeric_options.allow_numeric = true;
    numeric_options.has_initial_guess = true;
    numeric_options.initial_guess = 0.5;
    LMCAS::ComputationContext numeric_context;
    auto numeric = LMCAS::solve_equation(
        sine_equation, "x", numeric_context, numeric_options);
    const auto* numeric_finite = numeric
        ? std::get_if<LMCAS::FiniteSolutions>(&numeric.value()) : nullptr;
    EXPECT_TRUE(numeric_finite && !numeric_finite->values.empty(),
                "explicit numeric solving returns verified finite candidates");

    TEST_CASE("Checked solve dispatcher preserves computation errors");

    auto symbolic_parameter_equation =
        SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    LMCAS::ComputationContext symbolic_parameter_context;
    auto symbolic_parameter = LMCAS::solve_equation(
        symbolic_parameter_equation, "x",
        symbolic_parameter_context, numeric_options);
    EXPECT_TRUE(symbolic_parameter &&
                    std::holds_alternative<LMCAS::FiniteSolutions>(
                        symbolic_parameter.value()),
                "symbolic linear coefficients remain exact");

    LMCAS::CancellationToken cancellation;
    cancellation.cancel();
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    auto cancelled = LMCAS::solve_equation(
        x2_plus_one, "x", cancelled_context, LMCAS::SolveOptions{});
    EXPECT_TRUE(!cancelled && cancelled.error().code == LMCAS::CasErrc::Cancelled,
                "dispatcher observes cancellation");

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::solve_equation(
        x2_plus_one, "x", limited_context, LMCAS::SolveOptions{});
    EXPECT_TRUE(!limited && limited.error().code == LMCAS::CasErrc::ResourceLimit,
                "dispatcher enforces shared step budgets");

    TEST_CASE("Integrator preserves cancellation and recursive budgets");
    LMCAS::Integrator checked_integrator;
    LMCAS::CancellationToken integration_cancellation;
    integration_cancellation.cancel();
    LMCAS::ComputationContext cancelled_integration_context(
        {}, integration_cancellation);
    auto cancelled_integral = checked_integrator.integrate_checked(
        *x2_plus_one, "x", cancelled_integration_context);
    EXPECT_TRUE(!cancelled_integral &&
                    cancelled_integral.error().code ==
                        LMCAS::CasErrc::Cancelled,
                "integration preserves cancellation");

    LMCAS::ResourceLimits integration_limits;
    integration_limits.max_steps = 2;
    LMCAS::ComputationContext limited_integration_context(
        integration_limits);
    auto recursive_integrand = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(x));
    auto limited_integral = checked_integrator.integrate_checked(
        *recursive_integrand, "x", limited_integration_context);
    EXPECT_TRUE(!limited_integral &&
                    limited_integral.error().code ==
                        LMCAS::CasErrc::ResourceLimit,
                "integration preserves recursive step budgets");



    TEST_CASE("LMCAS number domains classify exact algebraic expressions");
    auto sqrt_two_expression = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    auto sqrt_two_real = LMCAS::domain_contains(
        LMCAS::reals(), sqrt_two_expression);
    auto sqrt_two_rational = LMCAS::domain_contains(
        LMCAS::rationals(), sqrt_two_expression);
    EXPECT_TRUE(sqrt_two_real && sqrt_two_real.value(),
                "sqrt(2) is proven real");
    EXPECT_TRUE(sqrt_two_rational && !sqrt_two_rational.value(),
                "sqrt(2) is proven non-rational");

    TEST_CASE("Quantity powers distinguish dimensioned and dimensionless bases");
    LMCAS::ComputationContext quantity_context;
    auto dimensionless_power = LMCAS::quantity_power(
        SymbolicExpr::number(2), SymbolicExpr::number(0.5),
        quantity_context);
    EXPECT_TRUE(dimensionless_power.has_value(),
                "dimensionless quantities accept approximate real exponents");
    auto metre = LMCAS::attach_unit(
        SymbolicExpr::number(1), "m", quantity_context);
    EXPECT_TRUE(metre.has_value(), "metre quantity construction succeeds");
    if (metre) {
        auto approximate_dimensioned = LMCAS::quantity_power(
            metre.value(), SymbolicExpr::number(2.0), quantity_context);
        EXPECT_TRUE(!approximate_dimensioned &&
                        approximate_dimensioned.error().code ==
                            LMCAS::CasErrc::UnitInvalid,
                    "dimensioned quantities require exact rational exponents");
    }

    TEST_CASE("Certified real algebraic values compare without floating point");
    LMCAS::ComputationContext algebraic_context;
    auto sqrt_two = LMCAS::detail::make_exact_real_algebraic(
        LMCAS::Polynomial<Rational>(
            {Rational(-2), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    auto sqrt_three = LMCAS::detail::make_exact_real_algebraic(
        LMCAS::Polynomial<Rational>(
            {Rational(-3), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    EXPECT_TRUE(sqrt_two.has_value() && sqrt_three.has_value(),
                "positive square roots have certified isolating intervals");
    if (sqrt_two && sqrt_three) {
        auto ordering = LMCAS::detail::compare_exact_real_algebraic(
            sqrt_two.value(), sqrt_three.value(), algebraic_context);
        EXPECT_TRUE(ordering.has_value() && ordering.value() < 0,
                    "sqrt(2) is certified less than sqrt(3)");
        auto equality = LMCAS::detail::equal_exact_real_algebraic(
            sqrt_two.value(), sqrt_two.value(), algebraic_context);
        EXPECT_TRUE(equality.has_value() && equality.value(),
                    "same certified root compares equal");
    }

    TEST_CASE("IntegralNode binds variables and supports calculus");
    auto y = SymbolicExpr::variable("y");
    auto body = SymbolicExpr::multiply(x, y);
    auto integral_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(body), "x"));
    EXPECT_TRUE(integral_expression->to_string() == "Integral(x*y, x)" ||
                    integral_expression->to_string() == "Integral(y*x, x)",
                "integral has a first-class printed form");
    auto parsed_integral = LMCAS::parse_expr(
        integral_expression->to_string());
    EXPECT_TRUE(parsed_integral.has_value() &&
                    LMCAS::detail::node(parsed_integral.value())->equals(
                        *LMCAS::detail::node(integral_expression)),
                "integral print/parse round-trip preserves binder structure");
    EXPECT_TRUE(!LMCAS::Integrator::depends_on(*integral_expression, "x") &&
                    LMCAS::Integrator::depends_on(*integral_expression, "y"),
                "integral variable is bound in free-variable analysis");
    auto capture_avoiding = integral_expression->substitute("y", x);
    EXPECT_TRUE(LMCAS::Integrator::depends_on(*capture_avoiding, "x"),
                "substitution does not capture a replacement variable");
    auto primitive = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(x), "x"));
    auto derivative = primitive->differentiate("x");
    EXPECT_TRUE(derivative && derivative->to_string() == "x",
                "d/dx Integral(x,x) returns the integrand");

    return TEST_REPORT();
}
