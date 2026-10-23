#include "internal/inequality_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_transform.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/rewrite_budget.hpp"
#include <algorithm>
#include <limits>
#include <optional>

namespace LMCAS {
namespace {

constexpr const char* kParametricOperation = "solve_parametric_inequality_checked";

using PiecewiseResult = Result<PiecewiseIntervalResult>;

PiecewiseResult checked_piecewise(PiecewiseResult result, detail::RewriteBudget& budget) {
    if (!result) { return result; }
    budget.consume();
    std::size_t nodes = 0;
    for (const auto& branch : result.value().cases) {
        nodes = budget.append_size(nodes, budget.measure(detail::node(branch.condition)));
        for (const auto& interval : branch.solution.intervals()) {
            nodes = budget.append_size(nodes, budget.measure(detail::node(interval.lower.value)));
            nodes = budget.append_size(nodes, budget.measure(detail::node(interval.upper.value)));
        }
    }
    return result;
}

PiecewiseResult inconclusive() {
    return PiecewiseResult::failure(
        CasErrc::Inconclusive,
        "parameter signs or quadratic roots cannot be certified",
        kParametricOperation);
}

IntervalUnion zero_polynomial_solution(InequalityType type) {
    return type == InequalityType::GreaterEqual || type == InequalityType::LessEqual
        ? IntervalUnion::entire_line()
        : IntervalUnion::empty();
}

bool satisfies(int sign, InequalityType type) {
    switch (type) {
        case InequalityType::GreaterThan: return sign > 0;
        case InequalityType::GreaterEqual: return sign >= 0;
        case InequalityType::LessThan: return sign < 0;
        case InequalityType::LessEqual: return sign <= 0;
    }
    return false;
}

std::shared_ptr<SymbolicExpr> sign_condition(
    const std::shared_ptr<SymbolicExpr>& coefficient, int sign) {
    const auto op = sign > 0 ? RelationalNode::Op::GT
        : sign < 0 ? RelationalNode::Op::LT : RelationalNode::Op::EQ;
    return detail::make_expression_ptr(detail::make_node<RelationalNode>(
        detail::node(coefficient), detail::node(SymbolicExpr::number(0)), op));
}

template <typename Node>
Result<void> exact_operands(const Node& node, const std::string& variable,
                            const std::vector<std::string>& parameters,
                            ComputationContext& context, std::size_t depth);

Result<void> exact_polynomial_nodes(
    const detail::SymbolicNodePtr& node, const std::string& variable,
    const std::vector<std::string>& parameters, ComputationContext& context,
    std::size_t depth) {
    if (depth >= context.limits().max_recursion_depth) {
        return Result<void>::failure(CasErrc::ResourceLimit,
            "parameter expression recursion limit exceeded", kParametricOperation);
    }
    auto step = context.consume_steps(1, kParametricOperation);
    if (!step) return step;
    if (const auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (detail::exact_rational_value(*number)) return Result<void>::success();
    } else if (const auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (!symbol->is_constant() &&
            (symbol->name() == variable ||
             std::find(parameters.begin(), parameters.end(), symbol->name()) != parameters.end())) {
            return Result<void>::success();
        }
    } else if (const auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        return exact_operands(*sum, variable, parameters, context, depth);
    } else if (const auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return exact_operands(*product, variable, parameters, context, depth);
    } else if (const auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (exact_small_integer_node(
                power->exponent(), 0, std::numeric_limits<int>::max())) {
            return exact_polynomial_nodes(
                power->base(), variable, parameters, context, depth + 1);
        }
    }
    return Result<void>::failure(CasErrc::Inconclusive,
        "parameter expression is not an exact polynomial", kParametricOperation);
}

template <typename Node>
Result<void> exact_operands(const Node& node, const std::string& variable,
                            const std::vector<std::string>& parameters,
                            ComputationContext& context, std::size_t depth) {
    for (const auto& operand : node.operands()) {
        auto checked = exact_polynomial_nodes(
            operand, variable, parameters, context, depth + 1);
        if (!checked) return checked;
    }
    return Result<void>::success();
}

Result<Interval> linear_interval(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    InequalityType type, detail::RewriteBudget& budget) {
    const auto slope = detail::exact_rational_value(detail::node(coeffs[1]));
    if (!slope) {
        return Result<Interval>::failure(inconclusive().error());
    }
    if (*slope == Rational(0)) {
        return Result<Interval>::failure(inconclusive().error());
    }
    std::optional<detail::AssumptionFacts> facts;
    if (budget.context().assumptions()) { facts.emplace(*budget.context().assumptions()); }
    auto root = detail::make_expression_ptr(detail::normalize_node(
        detail::node(SymbolicExpr::divide(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), coeffs[0]), coeffs[1])),
        budget, facts ? *facts : detail::no_facts()));

    const bool greater = satisfies(1, type);
    const bool strict = type == InequalityType::GreaterThan ||
                        type == InequalityType::LessThan;
    const bool above = greater != slope->get_numerator().is_negative();
    const auto endpoint = strict ? Endpoint::open(root)
                                 : Endpoint::closed(root);
    return above ? Interval{endpoint, Endpoint::pos_inf()}
                 : Interval{Endpoint::neg_inf(), endpoint};
}

PiecewiseResult solved_inequality(
    const std::shared_ptr<SymbolicExpr>& expr, InequalityType type,
    const std::string& variable, ComputationContext& context) {
    auto solved = InequalitySolver::solve_inequality_checked(expr, type, variable, context);
    if (!solved) return PiecewiseResult::failure(solved.error());
    PiecewiseIntervalResult result;
    result.cases.push_back({nullptr, std::move(solved.value())});
    return result;
}

PiecewiseResult constant_solution(const std::shared_ptr<SymbolicExpr>& coefficient,
                                  InequalityType type) {
    PiecewiseIntervalResult result;
    for (int sign : {1, 0, -1}) {
        result.cases.push_back({
            sign_condition(coefficient, sign),
            satisfies(sign, type) ? IntervalUnion::entire_line()
                                  : IntervalUnion::empty()});
    }
    return result;
}

PiecewiseResult quadratic_solution(const std::shared_ptr<SymbolicExpr>& coefficient,
                                   InequalityType type, ComputationContext& context) {
    PiecewiseIntervalResult result;
    for (int sign : {1, 0, -1}) {
        auto solution = sign == 0 ? zero_polynomial_solution(type)
            : detail::inequality_support::build_parametric_intervals(
                {SymbolicExpr::number(0)}, {2}, sign, type, context);
        result.cases.push_back({sign_condition(coefficient, sign), std::move(solution)});
    }
    return result;
}
PiecewiseResult nonlinear_solution(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    int degree, InequalityType type, ComputationContext& context) {
    if (degree == 0) return constant_solution(coeffs[0], type);
    if (!coeffs[0]->is_zero() || !coeffs[1]->is_zero()) return inconclusive();
    return quadratic_solution(coeffs[2], type, context);
}

Result<void> check_arguments(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& variable,
    const std::vector<std::string>& parameters) {
    if (parameters.empty()) return Result<void>::success();
    if (!expr || !detail::node(expr) || variable.empty()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "inequality expression and variable must not be empty",
            kParametricOperation);
    }
    if (std::find(parameters.begin(), parameters.end(), variable) != parameters.end()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "inequality variable cannot also be a parameter",
            kParametricOperation);
    }
    return Result<void>::success();
}

Result<std::vector<std::shared_ptr<SymbolicExpr>>> simplify_coefficients(
    const Polynomial<SymbolicPolyCoeff>& poly, bool& symbolic,
    detail::RewriteBudget& budget) {
    std::vector<std::shared_ptr<SymbolicExpr>> coeffs;
    coeffs.reserve(poly.coeffs.size());
    std::optional<detail::AssumptionFacts> facts;
    if (budget.context().assumptions()) { facts.emplace(*budget.context().assumptions()); }
    std::size_t nodes = 0;
    for (const auto& coefficient : poly.coeffs) {
        auto value = detail::make_expression_ptr(detail::normalize_node(
            detail::node(coefficient.val), budget, facts ? *facts : detail::no_facts()));
        nodes = budget.append_size(nodes, budget.measure(detail::node(value)));
        if (!detail::exact_rational_value(detail::node(value))) symbolic = true;
        coeffs.push_back(std::move(value));
    }
    return coeffs;
}
} // namespace

Result<PiecewiseIntervalResult> InequalitySolver::solve_parametric_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters,
    ComputationContext& context) {
    auto arguments = check_arguments(expr, variable, parameters);
    if (!arguments) {
        return PiecewiseResult::failure(arguments.error());
    }

    try {
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                     context.limits().max_ast_nodes, kParametricOperation);
        budget.measure(detail::node(expr));
        if (parameters.empty()) {
            return checked_piecewise(solved_inequality(expr, type, variable, context), budget);
        }
        auto converted = detail::symbolic_to_poly_checked(*expr, variable, context);
        if (!converted) {
            return PiecewiseResult::failure(converted.error());
        }
        const auto& poly = converted.value();
        if (poly.degree() > 2) {
            return inconclusive();
        }

        auto exact = exact_polynomial_nodes(
            detail::node(expr), variable, parameters, context, 0);
        if (!exact) {
            return PiecewiseResult::failure(exact.error());
        }

        if (poly.is_zero()) {
            PiecewiseIntervalResult result;
            result.cases.push_back({nullptr, zero_polynomial_solution(type)});
            return checked_piecewise(std::move(result), budget);
        }

        bool symbolic = false;
        auto simplified = simplify_coefficients(poly, symbolic, budget);
        if (!simplified) {
            return PiecewiseResult::failure(simplified.error());
        }
        const auto& coeffs = simplified.value();

        if (!symbolic) {
            return checked_piecewise(solved_inequality(expr, type, variable, context), budget);
        }

        if (poly.degree() != 1) {
            return checked_piecewise(nonlinear_solution(coeffs, poly.degree(), type, context), budget);
        }
        auto interval = linear_interval(coeffs, type, budget);
        if (!interval) {
            return PiecewiseResult::failure(interval.error());
        }
        PiecewiseIntervalResult result;
        result.cases.push_back({nullptr, IntervalUnion::from_checked_normalized(
            {std::move(interval.value())})});
        return checked_piecewise(std::move(result), budget);
    } catch (const CasError& error) {
        return PiecewiseResult::failure(error);
    } catch (const std::length_error& error) {
        return PiecewiseResult::failure(CasErrc::ResourceLimit, error.what(), kParametricOperation);
    } catch (const std::bad_alloc&) {
        return PiecewiseResult::failure(
            CasErrc::ResourceLimit, "parametric inequality allocation failed",
            kParametricOperation);
    } catch (const std::exception& error) {
        return PiecewiseResult::failure(
            CasErrc::InternalInvariant, error.what(), kParametricOperation);
    }
}

Result<PiecewiseIntervalResult> InequalitySolver::solve_parametric_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters) {
    ComputationContext context;
    return solve_parametric_inequality_checked(expr, type, variable, parameters, context);
}

} // namespace LMCAS
