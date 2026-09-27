#include "symbolic.hpp"
#include "solve_strategies.hpp"
#include "test_common.hpp"
#include "internal/exact_algebraic.hpp"
#include "integration.hpp"
#include "integrator.hpp"
#include "residual_verification.hpp"
#include "expr.hpp"
#include "quantity.hpp"
#include "assumption_context.hpp"
#include "internal/pointwise_comparison.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/assumption_simplification.hpp"

using namespace LMCAS;

TEST(ResultContracts, ResultValueAndMoveOwnership) {
    auto make_values = []() -> LMCAS::Result<std::vector<int>> {
        return std::vector<int>{1, 2, 3};
    };
    auto mutable_result = make_values();
    auto &mutable_value = mutable_result.value();
    mutable_value[0] = 7;
    EXPECT_TRUE((mutable_result.value()[0] == 7)) << "value access preserves the stored object";
    auto make_pointer =
        []() -> LMCAS::Result<std::unique_ptr<int>> {
        return std::make_unique<int>(9);
    };
    auto pointer_result = make_pointer();
    auto moved_value = std::move(pointer_result.value());
    EXPECT_TRUE((moved_value && *moved_value == 9)) << "direct success returns support move-only values";
}

TEST(ResultContracts, ResidualCertification) {
    auto x = SymbolicExpr::variable("x");
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
    EXPECT_TRUE((equivalent &&
                 std::holds_alternative<LMCAS::ProvedZeroResidual>(
                     equivalent.value())))
        << "expanded polynomial identity has an exact certificate";
    auto nonzero = LMCAS::check_zero_residual(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        residual_context);
    EXPECT_TRUE((nonzero &&
                 std::holds_alternative<LMCAS::ProvedNonIdentityResidual>(
                     nonzero.value())))
        << "a polynomial with a nonzero coefficient is not an identity";
    auto at_root = LMCAS::check_zero_residual(
        SymbolicExpr::add(x, SymbolicExpr::number(1))->substitute("x", SymbolicExpr::number(-1)), residual_context);
    EXPECT_TRUE((at_root && std::holds_alternative<LMCAS::ProvedZeroResidual>(at_root.value()))) << "a proved nonidentity can still vanish at a particular point";
    auto unproved = LMCAS::check_zero_residual(
        SymbolicExpr::sin(x), residual_context);
    EXPECT_TRUE((unproved &&
                 std::holds_alternative<LMCAS::UnprovedResidual>(
                     unproved.value())))
        << "unsupported transcendental residual remains unproved";
}

TEST(ResultContracts, RationalSquareRootScaling) {
    auto sqrt_two = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    auto half_sqrt_eight = SymbolicExpr::divide(
        SymbolicExpr::sqrt(SymbolicExpr::number(8)), SymbolicExpr::number(2));
    ComputationContext context;
    auto positive = check_equivalent(half_sqrt_eight, sqrt_two, context);
    ASSERT_TRUE(positive);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(positive.value()));
    auto negative = check_equivalent(
        SymbolicExpr::multiply(SymbolicExpr::number(-1), half_sqrt_eight),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), sqrt_two), context);
    ASSERT_TRUE(negative);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(negative.value()));
    auto opposite = check_equivalent(
        half_sqrt_eight,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), sqrt_two), context);
    ASSERT_TRUE(opposite);
    EXPECT_FALSE(std::holds_alternative<ProvedZeroResidual>(opposite.value()));
}

TEST(ResultContracts, RationalRadicandCanonicalization) {
    ComputationContext context;
    auto scaled = check_equivalent(
        SymbolicExpr::sqrt(SymbolicExpr::number(Rational(8, 9))),
        SymbolicExpr::multiply(SymbolicExpr::number(Rational(2, 3)),
                               SymbolicExpr::sqrt(SymbolicExpr::number(2))), context);
    ASSERT_TRUE(scaled);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(scaled.value()));
    auto rational = check_equivalent(
        SymbolicExpr::sqrt(SymbolicExpr::number(Rational(9, 4))),
        SymbolicExpr::number(Rational(3, 2)), context);
    ASSERT_TRUE(rational);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(rational.value()));
    auto power = check_equivalent(
        SymbolicExpr::power(SymbolicExpr::number(8), SymbolicExpr::number(Rational(1, 2))),
        SymbolicExpr::multiply(SymbolicExpr::number(2),
                               SymbolicExpr::sqrt(SymbolicExpr::number(2))), context);
    ASSERT_TRUE(power);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(power.value()));
}

TEST(ResultContracts, SquareRootProofPreservesPrincipalBranch) {
    ComputationContext context;
    auto x = SymbolicExpr::variable("x");
    auto symbolic = check_equivalent(
        SymbolicExpr::sqrt(SymbolicExpr::power(x, SymbolicExpr::number(2))), x, context);
    ASSERT_TRUE(symbolic);
    EXPECT_FALSE(std::holds_alternative<ProvedZeroResidual>(symbolic.value()));
    auto imaginary = check_equivalent(
        SymbolicExpr::sqrt(SymbolicExpr::number(-8)),
        SymbolicExpr::multiply(SymbolicExpr::number(-2),
                               SymbolicExpr::sqrt(SymbolicExpr::number(-2))), context);
    ASSERT_TRUE(imaginary);
    EXPECT_FALSE(std::holds_alternative<ProvedZeroResidual>(imaginary.value()));
}

TEST(ResultContracts, SquareRootProofUsesCallerLimits) {
    auto left = SymbolicExpr::divide(
        SymbolicExpr::sqrt(SymbolicExpr::number(8)), SymbolicExpr::number(2));
    auto right = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    ResourceLimits limits;
    limits.max_integer_bits = 3;
    ComputationContext integer_context(limits);
    auto oversized = check_equivalent(left, right, integer_context);
    ASSERT_FALSE(oversized);
    EXPECT_EQ(oversized.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(integer_context.recursion_depth(), 0u);
    limits = ResourceLimits{};
    limits.max_ast_nodes = 0;
    ComputationContext node_context(limits);
    auto unallocated = check_equivalent(left, right, node_context);
    ASSERT_FALSE(unallocated);
    EXPECT_EQ(unallocated.error().code, CasErrc::ResourceLimit);
    limits = ResourceLimits{};
    limits.max_steps = 128;
    ComputationContext step_context(limits);
    auto unfactored = check_equivalent(
        SymbolicExpr::sqrt(SymbolicExpr::number(1000003)),
        SymbolicExpr::number(1), step_context);
    ASSERT_FALSE(unfactored);
    EXPECT_EQ(unfactored.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(step_context.recursion_depth(), 0u);
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = check_equivalent(left, right, cancelled);
    ASSERT_FALSE(stopped);
    EXPECT_EQ(stopped.error().code, CasErrc::Cancelled);
}

TEST(ResultContracts, FiniteEmptyAndUniversalOutcomes) {
    auto x = SymbolicExpr::variable("x");

    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    LMCAS::ComputationContext exact_context;
    auto exact = LMCAS::solve_equation(
        x2_plus_one, "x", exact_context, LMCAS::SolveOptions{});
    const auto *exact_finite = exact
                                   ? std::get_if<LMCAS::FiniteSolutions>(&exact.value())
                                   : nullptr;
    EXPECT_TRUE((exact_finite && exact_finite->values.size() == 2)) << "exact polynomial preserves algebraic root count";
    LMCAS::ComputationContext empty_context;
    auto no_solution = LMCAS::solve_equation(
        SymbolicExpr::number(1), "x", empty_context, LMCAS::SolveOptions{});
    EXPECT_TRUE((no_solution &&
                 std::holds_alternative<LMCAS::EmptySolutions>(
                     no_solution.value())))
        << "nonzero constant equation has an empty solution set";

    LMCAS::ComputationContext universal_context;
    auto universal = LMCAS::solve_equation(
        SymbolicExpr::number(0), "x", universal_context, LMCAS::SolveOptions{});
    EXPECT_TRUE((universal &&
                 std::holds_alternative<LMCAS::UniversalSolutions>(
                     universal.value())))
        << "zero equation has the universal solution set";
}

TEST(ResultContracts, PeriodicAndUnsupportedOutcomes) {
    auto x = SymbolicExpr::variable("x");
    auto periodic = LMCAS::solve_equation(
        SymbolicExpr::sin(x), "x", LMCAS::SolveOptions{});
    const auto *periodic_values = periodic
                                      ? std::get_if<LMCAS::ParametricSolutions>(&periodic.value())
                                      : nullptr;
    EXPECT_TRUE((periodic_values && periodic_values->values.size() == 1)) << "sin(x)=0 returns a complete integer-parameter family";
    if (periodic_values) {
        EXPECT_TRUE((periodic_values->values[0].value->to_string().find("_k") !=
                     std::string::npos))
            << "periodic family preserves its integer parameter";
    }
    auto unsupported_equation = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<std::shared_ptr<const SymbolicNode>>{
                     LMCAS::detail::node(x)}));
    LMCAS::ComputationContext unsupported_context;
    auto unsupported = LMCAS::solve_equation(
        unsupported_equation, "x", unsupported_context, LMCAS::SolveOptions{});
    EXPECT_TRUE((!unsupported &&
                 unsupported.error().code == LMCAS::CasErrc::Inconclusive))
        << "unsupported symbolic equation returns outer Inconclusive";

    auto default_unsupported = LMCAS::solve_equation(
        unsupported_equation, "x", LMCAS::SolveOptions{});
    EXPECT_TRUE((!default_unsupported &&
                 default_unsupported.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "default-context solver preserves outer Inconclusive";
}

TEST(ResultContracts, SolveComputationErrors) {
    auto x = SymbolicExpr::variable("x");
    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    LMCAS::SolveOptions numeric_options;
    numeric_options.allow_numeric = true;
    numeric_options.has_initial_guess = true;
    numeric_options.initial_guess = 0.5;

    auto symbolic_parameter_equation =
        SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    LMCAS::ComputationContext symbolic_parameter_context;
    auto symbolic_parameter = LMCAS::solve_equation(
        symbolic_parameter_equation, "x",
        symbolic_parameter_context, numeric_options);
    EXPECT_TRUE((symbolic_parameter &&
                 std::holds_alternative<LMCAS::FiniteSolutions>(
                     symbolic_parameter.value())))
        << "symbolic linear coefficients remain exact";

    LMCAS::CancellationToken cancellation;
    cancellation.cancel();
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    auto cancelled = LMCAS::solve_equation(
        x2_plus_one, "x", cancelled_context, LMCAS::SolveOptions{});
    EXPECT_TRUE((!cancelled && cancelled.error().code == LMCAS::CasErrc::Cancelled)) << "dispatcher observes cancellation";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::solve_equation(
        x2_plus_one, "x", limited_context, LMCAS::SolveOptions{});
    EXPECT_TRUE((!limited && limited.error().code == LMCAS::CasErrc::ResourceLimit)) << "dispatcher enforces shared step budgets";
}

TEST(ResultContracts, IntegratorComputationErrors) {
    auto x = SymbolicExpr::variable("x");
    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    LMCAS::Integrator checked_integrator;
    LMCAS::CancellationToken integration_cancellation;
    integration_cancellation.cancel();
    LMCAS::ComputationContext cancelled_integration_context(
        {}, integration_cancellation);
    auto cancelled_integral = checked_integrator.integrate_checked(
        *x2_plus_one, "x", cancelled_integration_context);
    EXPECT_TRUE((!cancelled_integral &&
                 cancelled_integral.error().code ==
                     LMCAS::CasErrc::Cancelled))
        << "integration preserves cancellation";

    LMCAS::ResourceLimits integration_limits;
    integration_limits.max_steps = 2;
    LMCAS::ComputationContext limited_integration_context(
        integration_limits);
    auto recursive_integrand = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(x));
    auto limited_integral = checked_integrator.integrate_checked(
        *recursive_integrand, "x", limited_integration_context);
    EXPECT_TRUE((!limited_integral &&
                 limited_integral.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "integration preserves recursive step budgets";

    auto oversized_exponential = SymbolicExpr::exp(
        SymbolicExpr::power(x, SymbolicExpr::number(1000)));
    LMCAS::ComputationContext conversion_context;
    auto conversion_limited = checked_integrator.integrate_checked(
        *oversized_exponential, "x", conversion_context);
    EXPECT_TRUE((!conversion_limited &&
                 conversion_limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked integration preserves polynomial conversion limits";

    LMCAS::SpecialFunctionStrategy special_functions;
    LMCAS::ComputationContext special_context;
    auto special_limited = special_functions.try_integrate(
        *oversized_exponential, "x", checked_integrator, special_context, 0);
    EXPECT_TRUE((!special_limited &&
                 special_limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "special-function matching preserves conversion limits";
}

TEST(ResultContracts, NumberDomainClassification) {
    auto sqrt_two_expression = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    auto sqrt_two_real = LMCAS::domain_contains(
        LMCAS::reals(), sqrt_two_expression);
    auto sqrt_two_rational = LMCAS::domain_contains(
        LMCAS::rationals(), sqrt_two_expression);
    EXPECT_TRUE((sqrt_two_real && sqrt_two_real.value())) << "sqrt(2) is proven real";
    EXPECT_TRUE((sqrt_two_rational && !sqrt_two_rational.value())) << "sqrt(2) is proven non-rational";
}

TEST(ResultContracts, DimensionedPowerDomains) {
    LMCAS::ComputationContext quantity_context;
    auto dimensionless_power = LMCAS::quantity_power(
        SymbolicExpr::number(2), SymbolicExpr::number(0.5),
        quantity_context);
    EXPECT_TRUE((dimensionless_power.has_value())) << "dimensionless quantities accept approximate real exponents";
    auto metre = LMCAS::attach_unit(
        SymbolicExpr::number(1), "m", quantity_context);
    ASSERT_TRUE((metre.has_value())) << "metre quantity construction succeeds";
    if (metre) {
        auto approximate_dimensioned = LMCAS::quantity_power(
            metre.value(), SymbolicExpr::number(2.0), quantity_context);
        EXPECT_TRUE((!approximate_dimensioned &&
                     approximate_dimensioned.error().code ==
                         LMCAS::CasErrc::UnitInvalid))
            << "dimensioned quantities require exact rational exponents";
    }
}

TEST(ResultContracts, CertifiedAlgebraicOrdering) {
    LMCAS::ComputationContext algebraic_context;
    auto sqrt_two = LMCAS::detail::make_exact_real_algebraic(
        LMCAS::Polynomial<Rational>(
            {Rational(-2), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    auto sqrt_three = LMCAS::detail::make_exact_real_algebraic(
        LMCAS::Polynomial<Rational>(
            {Rational(-3), Rational(0), Rational(1)}),
        1, 1, algebraic_context);
    EXPECT_TRUE((sqrt_two.has_value() && sqrt_three.has_value())) << "positive square roots have certified isolating intervals";
    if (sqrt_two && sqrt_three) {
        auto ordering = LMCAS::detail::compare_exact_real_algebraic(
            sqrt_two.value(), sqrt_three.value(), algebraic_context);
        EXPECT_TRUE((ordering.has_value() && ordering.value() < 0)) << "sqrt(2) is certified less than sqrt(3)";
        auto equality = LMCAS::detail::equal_exact_real_algebraic(
            sqrt_two.value(), sqrt_two.value(), algebraic_context);
        EXPECT_TRUE((equality.has_value() && equality.value())) << "same certified root compares equal";
    }
}

TEST(ResultContracts, IntegralBindingAndCalculus) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto body = SymbolicExpr::multiply(x, y);
    auto integral_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(body), "x"));
    EXPECT_TRUE((integral_expression->to_string() == "Integral(x*y, x)" ||
                 integral_expression->to_string() == "Integral(y*x, x)"))
        << "integral has a first-class printed form";
    auto parsed_integral = LMCAS::parse_expr(
        integral_expression->to_string());
    EXPECT_TRUE((parsed_integral.has_value() &&
                 LMCAS::detail::node(parsed_integral.value())->equals(*LMCAS::detail::node(integral_expression))))
        << "integral print/parse round-trip preserves binder structure";
    EXPECT_TRUE((!LMCAS::Integrator::depends_on(*integral_expression, "x") &&
                 LMCAS::Integrator::depends_on(*integral_expression, "y")))
        << "integral variable is bound in free-variable analysis";
    auto capture_avoiding = integral_expression->substitute("y", x);
    EXPECT_TRUE((LMCAS::Integrator::depends_on(*capture_avoiding, "x"))) << "substitution does not capture a replacement variable";
    auto primitive = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(x), "x"));
    auto derivative = primitive->differentiate("x");
    EXPECT_TRUE((derivative && derivative->to_string() == "x")) << "d/dx Integral(x,x) returns the integrand";
}

TEST(ResultContracts, PointwiseResidualObligations) {
    auto a = SymbolicExpr::variable("a");
    auto zero = SymbolicExpr::number(0);
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "attach real assumptions";
    auto nonidentity = check_equivalent(a, zero, context);
    EXPECT_TRUE((nonidentity && std::holds_alternative<ProvedNonIdentityResidual>(nonidentity.value()))) << "a is formally not the zero polynomial";
    auto unequal = detail::compare_pointwise_values(a, zero, context, Domain::Real);
    EXPECT_TRUE((unequal && unequal.value() == Tribool::Unknown)) << "formal nonidentity cannot prove pointwise inequality for a real parameter";
    auto reciprocal = SymbolicExpr::divide(SymbolicExpr::number(1), a);
    auto identity = check_equivalent(reciprocal, reciprocal, context);
    EXPECT_TRUE((identity && std::holds_alternative<ProvedZeroResidual>(identity.value()))) << "identical reciprocals have a common-domain identity";
    auto equal = detail::compare_pointwise_values(reciprocal, reciprocal, context, Domain::Real);
    EXPECT_TRUE((equal && equal.value() == Tribool::Unknown)) << "even identical expressions require their denominator to be defined";
}

TEST(ResultContracts, PointwiseComparisonDomains) {
    ComputationContext context;
    auto complex = SymbolicExpr::sqrt(SymbolicExpr::number(-1));
    auto real_equal = detail::compare_pointwise_values(complex, complex, context, Domain::Real);
    auto complex_equal = detail::compare_pointwise_values(complex, complex, context, Domain::Complex);
    EXPECT_TRUE((real_equal && real_equal.value() == Tribool::Unknown &&
                 complex_equal && complex_equal.value() == Tribool::True))
        << "the explicit comparison domain controls definedness";
}

TEST(ResultContracts, PointwiseComparisonErrors) {
    auto zero = SymbolicExpr::number(0);
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = detail::compare_pointwise_values(zero, zero, cancelled, Domain::Real);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "structurally equal values cannot bypass cancellation";
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    auto limited = detail::compare_pointwise_values(zero, zero, exhausted, Domain::Real);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "pointwise comparison retains the caller's exhausted proof budget";
}

TEST(ResultContracts, PointwiseTraversalUsesCallerLimits) {
    detail::SymbolicNodePtr node = detail::make_node<VariableNode>("x");
    for (int i = 0; i < 8; ++i) {
        node = detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp, std::vector<detail::SymbolicNodePtr>{node});
    }
    auto expression = detail::make_expression_ptr(node);
    ResourceLimits steps;
    steps.max_steps = 1;
    ComputationContext step_context(steps);
    auto step_result = detail::compare_pointwise_values(
        expression, expression, step_context, Domain::Real);
    ASSERT_FALSE(step_result);
    EXPECT_EQ(step_result.error().code, CasErrc::ResourceLimit);
    ResourceLimits depth;
    depth.max_recursion_depth = 1;
    ComputationContext depth_context(depth);
    auto depth_result = detail::compare_pointwise_values(
        expression, expression, depth_context, Domain::Real);
    ASSERT_FALSE(depth_result);
    EXPECT_EQ(depth_result.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(depth_context.recursion_depth(), 0u);
    auto shallow = detail::query_definedness(
        detail::node(SymbolicExpr::number(1)), detail::no_facts(), Domain::Real, depth_context);
    ASSERT_TRUE(shallow);
    EXPECT_EQ(shallow.value(), Tribool::True);
}

TEST(ResultContracts, ValueFactsPreserveIntegerAndCumulativeLimits) {
    auto expression = parse_expr("(-2)^(2^200)");
    ASSERT_TRUE(expression);
    ResourceLimits limits;
    limits.max_integer_bits = 32;
    ComputationContext integer_context(limits);
    auto value = detail::query_nonzero_value(
        detail::node(expression.value()), detail::no_facts(), Domain::Real, integer_context);
    ASSERT_FALSE(value);
    EXPECT_EQ(value.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(integer_context.recursion_depth(), 0u);
    limits = ResourceLimits{};
    limits.max_steps = 1;
    ComputationContext cumulative(limits);
    const auto literal = detail::node(SymbolicExpr::number(1));
    auto first = detail::query_definedness(literal, detail::no_facts(), Domain::Real, cumulative);
    ASSERT_TRUE(first);
    EXPECT_EQ(first.value(), Tribool::True);
    auto second = detail::query_definedness(literal, detail::no_facts(), Domain::Real, cumulative);
    ASSERT_FALSE(second);
    EXPECT_EQ(second.error().code, CasErrc::ResourceLimit);
}

namespace {
class CancellingFacts final : public FactsQuery {
    static Result<Tribool> query(ComputationContext& context) {
        auto access = context.consume_steps(0, "test.cancel_facts");
        if (!access) { return Result<Tribool>::failure(access.error()); }
        context.cancellation().cancel();
        return Tribool::Unknown;
    }
public:
    Result<Tribool> is_positive(const detail::SymbolicNodePtr&, ComputationContext& context) const override { return query(context); }
    Result<Tribool> is_negative(const detail::SymbolicNodePtr&, ComputationContext& context) const override { return query(context); }
    Result<Tribool> is_nonnegative(const detail::SymbolicNodePtr&, ComputationContext& context) const override { return query(context); }
    Result<Tribool> is_nonzero(const detail::SymbolicNodePtr&, ComputationContext& context) const override { return query(context); }
    Result<Tribool> is_real(const detail::SymbolicNodePtr&, ComputationContext& context) const override { return query(context); }
};
}

TEST(ResultContracts, FactCallbackCancellationRestoresDepth) {
    CancellingFacts facts;
    ComputationContext context;
    auto result = detail::query_positive_value(
        detail::node(SymbolicExpr::variable("x")), facts, context);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Cancelled);
    EXPECT_EQ(context.recursion_depth(), 0u);
    auto subsequent = detail::query_definedness(
        detail::node(SymbolicExpr::number(1)), detail::no_facts(), Domain::Real, context);
    ASSERT_FALSE(subsequent);
    EXPECT_EQ(subsequent.error().code, CasErrc::Cancelled);
}

TEST(ResultContracts, LiteralSimplificationHonorsDepthAndCancellation) {
    auto literal = SymbolicExpr::number(1);
    ResourceLimits limits;
    limits.max_recursion_depth = 0;
    ComputationContext bounded(limits);
    auto exhausted = detail::simplify_expression(literal, bounded);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(bounded.recursion_depth(), 0u);
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = detail::simplify_expression(literal, cancelled);
    ASSERT_FALSE(stopped);
    EXPECT_EQ(stopped.error().code, CasErrc::Cancelled);
    auto defined = detail::query_definedness(detail::node(literal), detail::no_facts(),
                                             Domain::Real, cancelled);
    ASSERT_FALSE(defined);
    EXPECT_EQ(defined.error().code, CasErrc::Cancelled);
    AssumptionContext assumptions;
    ASSERT_TRUE(assumptions.assume_sign_checked("x", Sign::Positive));
    detail::AssumptionFacts facts(assumptions);
    auto known = facts.is_positive(detail::node(SymbolicExpr::variable("x")), cancelled);
    ASSERT_FALSE(known);
    EXPECT_EQ(known.error().code, CasErrc::Cancelled);
}
