#include "../include/symbolic.hpp"
#include "../include/solve_strategies.hpp"
#include "test_common.hpp"

int main() {
    TEST_CASE("SolutionSet distinguishes empty, finite, and inconclusive");

    auto empty = lamina::SolutionSet::empty();
    EXPECT_TRUE(empty.kind() == lamina::SolutionSet::Kind::Empty,
                "empty solution set represents mathematical no-solution");

    auto x = SymbolicExpr::variable("x");
    auto finite = lamina::SolutionSet::finite({
        lamina::FiniteSolution{x, 2, {}}
    });
    EXPECT_TRUE(finite.kind() == lamina::SolutionSet::Kind::Finite,
                "finite solution set preserves explicit kind");
    EXPECT_TRUE(finite.finite_solutions().size() == 1,
                "finite solution set stores roots");
    EXPECT_TRUE(finite.finite_solutions()[0].multiplicity == 2,
                "finite solution set preserves multiplicity");

    auto inconclusive = lamina::SolutionSet::inconclusive("unsupported transcendental equation");
    EXPECT_TRUE(inconclusive.kind() == lamina::SolutionSet::Kind::Inconclusive,
                "unsupported domain is represented as inconclusive");
    EXPECT_TRUE(inconclusive.reason() == "unsupported transcendental equation",
                "inconclusive result carries a reason");

    TEST_CASE("Conditional and verification status are explicit");

    lamina::Conditional<std::shared_ptr<SymbolicExpr>> conditional{
        x,
        {SymbolicExpr::variable("x_nonzero")},
        lamina::VerificationStatus::Verified
    };
    EXPECT_TRUE(conditional.value == x, "conditional stores value");
    EXPECT_TRUE(conditional.conditions.size() == 1, "conditional stores conditions");
    EXPECT_TRUE(conditional.verification == lamina::VerificationStatus::Verified,
                "conditional stores verification status");

    lamina::IntegralResult integral;
    integral.result = conditional;
    integral.singularities.push_back(SymbolicExpr::number(0));
    integral.verification = lamina::VerificationStatus::Inconclusive;
    EXPECT_TRUE(integral.singularities.size() == 1,
                "integral result carries singularities");
    EXPECT_TRUE(integral.verification == lamina::VerificationStatus::Inconclusive,
                "integral result can be explicitly inconclusive");

    TEST_CASE("Checked solve dispatcher distinguishes mathematical outcomes");

    auto x2_plus_one = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    lamina::ComputationContext exact_context;
    auto exact = lamina::solve_equation(
        x2_plus_one, "x", exact_context, lamina::SolveOptions{});
    EXPECT_TRUE(exact && exact.value().kind() == lamina::SolutionSet::Kind::Finite,
                "exact polynomial produces a finite RootOf set");
    EXPECT_TRUE(exact && exact.value().finite_solutions().size() == 2,
                "exact polynomial preserves algebraic root count");
    if (exact && exact.value().finite_solutions().size() == 2) {
        EXPECT_TRUE(exact.value().finite_solutions()[0].value->to_string().find("rootof") !=
                        std::string::npos,
                    "exact candidates remain exact RootOf expressions");
    }

    lamina::ComputationContext empty_context;
    auto no_solution = lamina::solve_equation(
        SymbolicExpr::number(1), "x", empty_context, lamina::SolveOptions{});
    EXPECT_TRUE(no_solution &&
                    no_solution.value().kind() == lamina::SolutionSet::Kind::Empty,
                "nonzero constant equation has an empty solution set");

    lamina::ComputationContext universal_context;
    auto universal = lamina::solve_equation(
        SymbolicExpr::number(0), "x", universal_context, lamina::SolveOptions{});
    EXPECT_TRUE(universal &&
                    universal.value().kind() == lamina::SolutionSet::Kind::Universal,
                "zero equation has the universal solution set");

    auto sine_equation = SymbolicExpr::add(
        SymbolicExpr::sin(x), SymbolicExpr::number(Rational(-1, 2)));
    lamina::ComputationContext unsupported_context;
    auto unsupported = lamina::solve_equation(
        sine_equation, "x", unsupported_context, lamina::SolveOptions{});
    EXPECT_TRUE(unsupported &&
                    unsupported.value().kind() == lamina::SolutionSet::Kind::Inconclusive,
                "unsupported symbolic equation is not reported as empty");

    auto default_unsupported = lamina::solve_equation(
        sine_equation, "x", lamina::SolveOptions{});
    EXPECT_TRUE(default_unsupported &&
                    default_unsupported.value().kind() == lamina::SolutionSet::Kind::Inconclusive,
                "default-context checked dispatcher preserves Inconclusive");

    lamina::SolveOptions numeric_options;
    numeric_options.allow_numeric = true;
    numeric_options.has_initial_guess = true;
    numeric_options.initial_guess = 0.5;
    lamina::ComputationContext numeric_context;
    auto numeric = lamina::solve_equation(
        sine_equation, "x", numeric_context, numeric_options);
    EXPECT_TRUE(numeric && numeric.value().kind() == lamina::SolutionSet::Kind::Finite &&
                    numeric.value().finite_solutions().size() == 1,
                "explicit numeric solving returns a verified finite candidate");

    TEST_CASE("Checked solve dispatcher preserves computation errors");

    auto unbound_equation = SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    lamina::ComputationContext unbound_context;
    auto unbound = lamina::solve_equation(
        unbound_equation, "x", unbound_context, numeric_options);
    EXPECT_TRUE(!unbound && unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "unbound coefficients propagate as UnboundSymbol");

    auto default_unbound = lamina::solve_equation(
        unbound_equation, "x", numeric_options);
    EXPECT_TRUE(!default_unbound &&
                    default_unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "default-context checked dispatcher propagates UnboundSymbol");

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

    TEST_CASE("Checked vector solve dispatcher exposes legacy numeric errors");

    lamina::ComputationContext vector_context;
    auto vector_roots = lamina::solve_dispatch_vector_checked(
        SymbolicExpr::add(x, SymbolicExpr::number(-1)),
        "x",
        vector_context,
        lamina::SolveOptions{});
    EXPECT_TRUE(vector_roots && vector_roots.value().size() == 1,
                "checked vector dispatcher preserves legacy vector results");

    auto abs_arg = SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    auto abs_expr = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(abs_arg)}));
    lamina::SolveOptions vector_numeric_options;
    vector_numeric_options.allow_numeric = true;
    vector_numeric_options.has_initial_guess = true;
    vector_numeric_options.initial_guess = 0.25;
    lamina::ComputationContext vector_unbound_context;
    auto vector_unbound = lamina::solve_dispatch_vector_checked(
        abs_expr, "x", vector_unbound_context, vector_numeric_options);
    EXPECT_TRUE(!vector_unbound &&
                    vector_unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "checked vector dispatcher propagates numeric fallback errors");

    auto default_vector_unbound = lamina::solve_dispatch_vector_checked(
        abs_expr, "x", vector_numeric_options);
    EXPECT_TRUE(!default_vector_unbound &&
                    default_vector_unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "default-context checked vector dispatcher propagates numeric fallback errors");

    auto legacy_suppressed = lamina::solve_dispatch(
        abs_expr, "x", vector_numeric_options);
    EXPECT_TRUE(legacy_suppressed.empty(),
                "legacy vector dispatcher still unwraps errors to an empty vector");

    return TEST_REPORT();
}
