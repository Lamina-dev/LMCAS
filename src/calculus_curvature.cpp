#include "solve_strategies.hpp"
#include "internal/calculus_utils_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/expression_transform.hpp"

namespace LMCAS {

using namespace calculus_utils_detail;

static ExpressionResult checked_derivative(
    const ExprPtr& expression, const std::string& variable,
    ComputationContext& context, const char* operation)
{
    return detail::checked_transform_expr(
        expression, context, operation, "differentiate",
        [&variable](const SymbolicExpr& value) { return value.differentiate(variable); });
}

ExpressionResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    const std::string operation = "curvature";
    auto input = calculus_utils_validate_expr(f, var, context, operation);
    if (!input) {
        return ExpressionResult::failure(input.error());
    }

    auto first = checked_derivative(f, var, context, operation.c_str());
    if (!first) {
        return ExpressionResult::failure(first.error());
    }
    const auto& f_prime = first.value();
    auto second = checked_derivative(f_prime, var, context, operation.c_str());
    if (!second) {
        return ExpressionResult::failure(second.error());
    }
    const auto& f_double_prime = second.value();

    auto abs_f_pp = calculus_utils_make_abs(f_double_prime);
    if (!abs_f_pp || !LMCAS::detail::node(abs_f_pp)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           "absolute value node construction failed",
                                           operation);
    }

    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);

    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(one_plus_fp_sq, three_half);
    auto result = SymbolicExpr::divide(abs_f_pp, denom);
    auto simplified = result->simplify();
    return ExpressionResult::success(simplified ? simplified : result);
}

ExpressionResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return curvature_checked(f, var, context);
}


ExpressionResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t,
    ComputationContext& context)
{
    const std::string operation = "curvature_parametric";
    auto input = calculus_utils_validate_two_exprs(x_t, y_t, t, context, operation);
    if (!input) {
        return ExpressionResult::failure(input.error());
    }

    auto x_first = checked_derivative(x_t, t, context, operation.c_str());
    if (!x_first) {
        return ExpressionResult::failure(x_first.error());
    }
    const auto& x_prime = x_first.value();
    auto x_second = checked_derivative(x_prime, t, context, operation.c_str());
    if (!x_second) {
        return ExpressionResult::failure(x_second.error());
    }
    const auto& x_double_prime = x_second.value();

    auto y_first = checked_derivative(y_t, t, context, operation.c_str());
    if (!y_first) {
        return ExpressionResult::failure(y_first.error());
    }
    const auto& y_prime = y_first.value();
    auto y_second = checked_derivative(y_prime, t, context, operation.c_str());
    if (!y_second) {
        return ExpressionResult::failure(y_second.error());
    }
    const auto& y_double_prime = y_second.value();

    if (x_prime->is_zero() && y_prime->is_zero()) {
        return ExpressionResult::failure(CasErrc::DomainError,
                                           "parametric curvature is undefined for zero velocity",
                                           operation);
    }

    auto cross_term = SymbolicExpr::add(
        SymbolicExpr::multiply(x_prime, y_double_prime),
        SymbolicExpr::multiply(
            SymbolicExpr::number(-1),
            SymbolicExpr::multiply(y_prime, x_double_prime)));

    auto abs_cross = calculus_utils_make_abs(cross_term);
    if (!abs_cross || !LMCAS::detail::node(abs_cross)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           "absolute value node construction failed",
                                           operation);
    }

    auto x_prime_sq = SymbolicExpr::power(x_prime, SymbolicExpr::number(2));
    auto y_prime_sq = SymbolicExpr::power(y_prime, SymbolicExpr::number(2));
    auto sum_sq = SymbolicExpr::add(x_prime_sq, y_prime_sq);

    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(sum_sq, three_half);
    auto result = SymbolicExpr::divide(abs_cross, denom);
    auto simplified = result->simplify();
    return ExpressionResult::success(simplified ? simplified : result);
}

ExpressionResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t)
{
    ComputationContext context;
    return curvature_parametric_checked(x_t, y_t, t, context);
}



static Result<bool> inflection_curve_defined(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& f_prime,
    const std::shared_ptr<SymbolicExpr>& f_double_prime,
    const std::shared_ptr<SymbolicExpr>& point, const std::string& var,
    const FactsQuery& facts, const std::string& operation, ComputationContext& context) {
    auto value = f->substitute(var, point);
    auto first = f_prime->substitute(var, point);
    auto second = f_double_prime->substitute(var, point);
    for (const auto& expression : {value, first, second}) {
        auto defined = detail::query_definedness(detail::node(expression), facts, Domain::Real, context);
        if (!defined) { return Result<bool>::failure(defined.error()); }
        if (defined.value() == Tribool::False) { return false; }
        auto real_value = detail::query_real_value(detail::node(expression), facts, context);
        if (!real_value) { return Result<bool>::failure(real_value.error()); }
        if (real_value.value() == Tribool::False) { return false; }
        if (defined.value() != Tribool::True || real_value.value() != Tribool::True) {
            return Result<bool>::failure(CasErrc::Inconclusive,
                "real curve domain at an inflection candidate is unproved",
                operation);
        }
    }
    return true;
}

static SymbolicExprVectorResult classify_inflection_candidates(
    const FiniteSolutions& finite, const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& f_prime,
    const std::shared_ptr<SymbolicExpr>& f_double_prime, const std::string& var,
    const FactsQuery& facts, ComputationContext& context, const std::string& operation) {
    std::vector<std::shared_ptr<SymbolicExpr>> points;
    points.reserve(finite.values.size());
    for (const auto& solution : finite.values) {
        auto step = context.consume_steps(1, operation);
        if (!step) { return SymbolicExprVectorResult::failure(step.error()); }
        if (!solution.value || !LMCAS::detail::node(solution.value)) {
            return SymbolicExprVectorResult::failure(CasErrc::InternalInvariant,
                "checked solver returned a null inflection candidate", operation);
        }
        auto real = detail::query_real_value(detail::node(solution.value), facts, context);
        if (!real) { return SymbolicExprVectorResult::failure(real.error()); }
        if (real.value() == Tribool::False) { continue; }
        if (real.value() != Tribool::True || !solution.conditions.empty()) {
            return SymbolicExprVectorResult::failure(CasErrc::Inconclusive,
                "inflection candidate is not proved unconditionally real", operation);
        }
        auto defined = inflection_curve_defined(
            f, f_prime, f_double_prime, solution.value, var, facts, operation, context);
        if (!defined) { return SymbolicExprVectorResult::failure(defined.error()); }
        if (!defined.value()) { continue; }
        points.push_back(solution.value);
    }
    return SymbolicExprVectorResult::success(std::move(points));
}

SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    const std::string operation = "inflection_points";
    auto input = calculus_utils_validate_expr(f, var, context, operation);
    if (!input) {
        return SymbolicExprVectorResult::failure(input.error());
    }

    auto first = checked_derivative(f, var, context, operation.c_str());
    if (!first) {
        return SymbolicExprVectorResult::failure(first.error());
    }
    const auto& f_prime = first.value();

    auto second = checked_derivative(f_prime, var, context, operation.c_str());
    if (!second) {
        return SymbolicExprVectorResult::failure(second.error());
    }
    const auto& f_double_prime = second.value();

    auto equation = f_double_prime->simplify();
    if (!equation || !LMCAS::detail::node(equation)) { equation = f_double_prime; }

    auto solved = solve_equation(equation, var, context, SolveOptions{});
    if (!solved) {
        return SymbolicExprVectorResult::failure(solved.error());
    }

    const auto& solutions = solved.value();
    if (std::holds_alternative<EmptySolutions>(solutions)) {
        return SymbolicExprVectorResult::success({});
    }
    std::optional<detail::AssumptionFacts> assumed_facts;
    if (context.assumptions()) { assumed_facts.emplace(*context.assumptions()); }
    const FactsQuery& facts = assumed_facts
        ? static_cast<const FactsQuery&>(*assumed_facts) : detail::no_facts();
    if (const auto* finite = std::get_if<FiniteSolutions>(&solutions)) {
        return classify_inflection_candidates(
            *finite, f, f_prime, f_double_prime, var, facts, context, operation);
    }


    return SymbolicExprVectorResult::failure(
        CasErrc::Inconclusive,
        "inflection equation is outside the finite exact support domain",
        operation);
}

SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return inflection_points_checked(f, var, context);
}

}
