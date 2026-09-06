#include "assumption_context.hpp"
#include "poly_utils.hpp"
#include "property_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "solve_strategies.hpp"

#include <iostream>

using namespace LMCAS;

int run_lmmc_linalg_consumer_checks();
int run_expr_consumer_checks();
int run_lmmc_stdlib_consumer_checks();

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(x, SymbolicExpr::number(1));
    if (!expr) {
        std::cerr << "failed to construct expression\n";
        return 1;
    }

    LMCAS::ComputationContext context;
    auto solved = LMCAS::solve_equation(
        expr, "x", context, LMCAS::SolveOptions{});
    const auto* finite = solved
        ? std::get_if<LMCAS::FiniteSolutions>(&solved.value()) : nullptr;
    if (!finite || finite->values.size() != 1) {
        std::cerr << "failed to solve expression\n";
        return 2;
    }

    auto polynomial = LMCAS::symbolic_to_poly<Rational>(expr, "x");
    if (polynomial.coeffs.size() != 2 ||
        polynomial.coeffs[0] != Rational(1) ||
        polynomial.coeffs[1] != Rational(1)) {
        std::cerr << "failed to convert expression to polynomial\n";
        return 3;
    }

    auto interval = LMCAS::Interval::point(SymbolicExpr::number(0));
    LMCAS::ComputationContext interval_context;
    auto interval_union = LMCAS::IntervalUnion::from_intervals_checked(
        {interval}, interval_context);
    if (!interval_union || interval_union.value().intervals().size() != 1) {
        std::cerr << "failed to construct interval union\n";
        return 4;
    }

    LMCAS::PropertyStore properties;
    LMCAS::ComputationContext property_context;
    auto declared = properties.declare_continuous_checked(
        "f", interval, property_context);
    if (!declared) {
        std::cerr << "failed to declare property\n";
        return 5;
    }
    LMCAS::ComputationContext query_context;
    auto continuous = properties.is_continuous_checked(
        "f", interval, query_context);
    if (!continuous || !continuous.value()) {
        std::cerr << "failed to query property\n";
        return 6;
    }

    LMCAS::AssumptionContext assumptions;
    LMCAS::QueryInterface queries(assumptions);
    auto one = SymbolicExpr::number(1);
    auto positive = queries.query_positive_checked(*one);
    if (!positive || positive.value() != LMCAS::Tribool::True) {
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
    if (!context_continuous || context_continuous.value() != LMCAS::Tribool::True) {
        std::cerr << "failed to query assumption property\n";
        return 9;
    }
    if (const int status = run_expr_consumer_checks(); status != 0) {
        return status;
    }
    if (const int status = run_lmmc_stdlib_consumer_checks(); status != 0) {
        return status;
    }

    if (const int linalg_status = run_lmmc_linalg_consumer_checks();
        linalg_status != 0) {
        return linalg_status;
    }
    std::cout << expr->to_string() << '\n';
    return 0;
}
