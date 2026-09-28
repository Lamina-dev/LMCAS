#include "assumption_context.hpp"
#include "expr.hpp"
#include "inequality_solver.hpp"
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

static int check_equation_and_polynomial(const std::shared_ptr<SymbolicExpr>& expr) {
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
    if (!polynomial || polynomial.value().coeffs.size() != 2 ||
        polynomial.value().coeffs[0] != Rational(1) ||
        polynomial.value().coeffs[1] != Rational(1)) {
        std::cerr << "failed to convert expression to polynomial\n";
        return 3;
    }
    return 0;
}

static int check_installed_expression_contracts(const ExprPtr& expr) {
    auto solved = solve_set(expr, "x");
    const auto* finite = solved ? std::get_if<FiniteSolutions>(&solved.value()) : nullptr;
    if (!finite || finite->values.size() != 1 ||
        finite->values[0].value->to_string() != "-1") {
        std::cerr << "installed solve_set lost the root of x+1\n";
        return 12;
    }
    ComputationContext context;
    auto distinct = equivalent(*expr, *SymbolicExpr::number(0), context);
    if (!distinct || distinct.value()) {
        std::cerr << "installed equivalent failed to prove x+1 differs from zero\n";
        return 12;
    }
    auto encoded = serialize_expr(expr);
    if (!encoded) {
        std::cerr << "installed expression encoding failed\n";
        return 12;
    }
    auto decoded = parse_serialized_expr(encoded.value());
    if (!decoded) {
        std::cerr << "installed expression decoding failed\n";
        return 12;
    }
    auto restored = serialize_expr(decoded.value());
    if (!restored || restored.value() != encoded.value()) {
        std::cerr << "installed expression encoding failed its round trip\n";
        return 12;
    }
    auto a = SymbolicExpr::variable("a");
    auto linear = SymbolicExpr::add(SymbolicExpr::variable("x"), a);
    auto positive = InequalitySolver::solve_parametric_inequality_checked(
        linear, InequalityType::GreaterThan, "x", {"a"});
    if (!positive || positive.value().cases.size() != 1 ||
        positive.value().cases[0].solution.intervals().size() != 1 ||
        !positive.value().cases[0].solution.intervals()[0].upper.is_pos_infinity) {
        std::cerr << "installed parametric inequality lost x+a>0\n";
        return 12;
    }
    return 0;
}

static int check_property_store(const LMCAS::Interval& interval) {
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
    return 0;
}

static int check_assumptions(const LMCAS::Interval& interval) {
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
    return 0;
}

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(x, SymbolicExpr::number(1));
    if (!expr) {
        std::cerr << "failed to construct expression\n";
        return 1;
    }

    if (const int status = check_equation_and_polynomial(expr); status != 0) {
        return status;
    }
    if (const int status = check_installed_expression_contracts(expr); status != 0) {
        return status;
    }
    auto interval = LMCAS::Interval::point(SymbolicExpr::number(0));
    LMCAS::ComputationContext interval_context;
    auto interval_union = LMCAS::IntervalUnion::from_intervals_checked(
        {interval}, interval_context);
    if (!interval_union || interval_union.value().intervals().size() != 1) {
        std::cerr << "failed to construct interval union\n";
        return 4;
    }

    if (const int status = check_property_store(interval); status != 0) {
        return status;
    }
    if (const int status = check_assumptions(interval); status != 0) {
        return status;
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
