#include "assumption_context.hpp"
#include "poly_utils.hpp"
#include "property_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "solve_strategies.hpp"

#include <iostream>

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(x, SymbolicExpr::number(1));
    if (!expr) {
        std::cerr << "failed to construct expression\n";
        return 1;
    }

    lamina::ComputationContext context;
    auto solved = lamina::solve_dispatch_checked(
        expr, "x", context, lamina::SolveOptions{});
    if (!solved || solved.value().kind() != lamina::SolutionSet::Kind::Finite ||
        solved.value().finite_solutions().size() != 1) {
        std::cerr << "failed to solve expression\n";
        return 2;
    }

    auto polynomial = lamina::symbolic_to_poly<Rational>(expr, "x");
    if (polynomial.coeffs.size() != 2 ||
        polynomial.coeffs[0] != Rational(1) ||
        polynomial.coeffs[1] != Rational(1)) {
        std::cerr << "failed to convert expression to polynomial\n";
        return 3;
    }

    auto interval = lamina::Interval::point(SymbolicExpr::number(0));
    lamina::ComputationContext interval_context;
    auto interval_union = lamina::IntervalUnion::from_intervals_checked(
        {interval}, interval_context);
    if (!interval_union || interval_union.value().intervals().size() != 1) {
        std::cerr << "failed to construct interval union\n";
        return 4;
    }

    lamina::PropertyStore properties;
    lamina::ComputationContext property_context;
    auto declared = properties.declare_continuous_checked(
        "f", interval, property_context);
    if (!declared) {
        std::cerr << "failed to declare property\n";
        return 5;
    }
    lamina::ComputationContext query_context;
    auto continuous = properties.is_continuous_checked(
        "f", interval, query_context);
    if (!continuous || !continuous.value()) {
        std::cerr << "failed to query property\n";
        return 6;
    }

    lamina::AssumptionContext assumptions;
    lamina::QueryInterface queries(assumptions);
    auto one = SymbolicExpr::number(1);
    auto positive = queries.query_positive_checked(*one);
    if (!positive || positive.value() != lamina::Tribool::True) {
        std::cerr << "failed to query positivity\n";
        return 7;
    }

    auto assumed_continuous = assumptions.current_properties().declare_continuous_checked(
        "g", interval);
    if (!assumed_continuous) {
        std::cerr << "failed to declare assumption property\n";
        return 8;
    }
    auto context_continuous = assumptions.is_continuous_checked("g", interval);
    if (!context_continuous || context_continuous.value() != lamina::Tribool::True) {
        std::cerr << "failed to query assumption property\n";
        return 9;
    }

    std::cout << expr->to_string() << '\n';
    if (expr->to_string().empty()) {
        std::cerr << "expression string is empty\n";
        return 10;
    }
    return 0;
}
