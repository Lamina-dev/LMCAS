#include "internal/inequality_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "expr.hpp"
#include <algorithm>
#include <limits>

namespace LMCAS {
namespace {

constexpr const char* kParametricOperation = "solve_parametric_inequality_checked";

using PiecewiseResult = Result<PiecewiseIntervalResult>;

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
        for (const auto& term : sum->operands()) {
            auto checked = exact_polynomial_nodes(term, variable, parameters, context, depth + 1);
            if (!checked) return checked;
        }
        return Result<void>::success();
    } else if (const auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& factor : product->operands()) {
            auto checked = exact_polynomial_nodes(factor, variable, parameters, context, depth + 1);
            if (!checked) return checked;
        }
        return Result<void>::success();
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

} // namespace

Result<PiecewiseIntervalResult> InequalitySolver::solve_parametric_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters,
    ComputationContext& context) {
    if (!parameters.empty()) {
        if (!expr || !detail::node(expr) || variable.empty()) {
            return PiecewiseResult::failure(
                CasErrc::InvalidArgument, "inequality expression and variable must not be empty",
                kParametricOperation);
        }
        if (std::find(parameters.begin(), parameters.end(), variable) != parameters.end()) {
            return PiecewiseResult::failure(
                CasErrc::InvalidArgument, "inequality variable cannot also be a parameter",
                kParametricOperation);
        }
    }

    try {
        if (parameters.empty()) {
            auto solved = solve_inequality_checked(expr, type, variable, context);
            if (!solved) return PiecewiseResult::failure(solved.error());
            PiecewiseIntervalResult result;
            result.cases.push_back({nullptr, std::move(solved.value())});
            return result;
        }
        auto converted = detail::symbolic_to_poly_checked(*expr, variable, context);
        if (!converted) return PiecewiseResult::failure(converted.error());
        const auto& poly = converted.value();
        if (poly.degree() > 2) return inconclusive();

        auto exact = exact_polynomial_nodes(
            detail::node(expr), variable, parameters, context, 0);
        if (!exact) return PiecewiseResult::failure(exact.error());

        if (poly.is_zero()) {
            PiecewiseIntervalResult result;
            result.cases.push_back({nullptr, zero_polynomial_solution(type)});
            return result;
        }

        std::vector<std::shared_ptr<SymbolicExpr>> coeffs;
        coeffs.reserve(poly.coeffs.size());
        bool symbolic = false;
        for (const auto& coefficient : poly.coeffs) {
            auto simplified = simplify(coefficient.val, context);
            if (!simplified) return PiecewiseResult::failure(simplified.error());
            auto value = std::move(simplified.value());
            if (!detail::exact_rational_value(detail::node(value))) symbolic = true;
            coeffs.push_back(std::move(value));
        }

        if (!symbolic) {
            auto solved = solve_inequality_checked(expr, type, variable, context);
            if (!solved) return PiecewiseResult::failure(solved.error());
            PiecewiseIntervalResult result;
            result.cases.push_back({nullptr, std::move(solved.value())});
            return result;
        }

        auto nodes = context.reserve_nodes(16, kParametricOperation);
        if (!nodes) return PiecewiseResult::failure(nodes.error());
        PiecewiseIntervalResult result;

        if (poly.degree() == 0) {
            for (int sign : {1, 0, -1}) {
                result.cases.push_back({
                    sign_condition(coeffs[0], sign),
                    satisfies(sign, type) ? IntervalUnion::entire_line()
                                          : IntervalUnion::empty()});
            }
            return result;
        }

        if (poly.degree() == 1) {
            const auto slope = detail::exact_rational_value(detail::node(coeffs[1]));
            if (!slope || *slope == Rational(0)) return inconclusive();
            auto root = simplify(SymbolicExpr::divide(
                SymbolicExpr::multiply(SymbolicExpr::number(-1), coeffs[0]),
                coeffs[1]), context);
            if (!root) return PiecewiseResult::failure(root.error());

            const bool greater = type == InequalityType::GreaterThan ||
                                 type == InequalityType::GreaterEqual;
            const bool strict = type == InequalityType::GreaterThan ||
                                type == InequalityType::LessThan;
            const bool above = greater != slope->get_numerator().is_negative();
            const auto endpoint = strict ? Endpoint::open(root.value())
                                         : Endpoint::closed(root.value());
            result.cases.push_back({nullptr, IntervalUnion::from_checked_normalized({
                above ? Interval{endpoint, Endpoint::pos_inf()}
                      : Interval{Endpoint::neg_inf(), endpoint}})});
            return result;
        }

        if (!coeffs[0]->is_zero() || !coeffs[1]->is_zero()) return inconclusive();
        for (int sign : {1, 0, -1}) {
            auto solution = sign == 0 ? zero_polynomial_solution(type)
                : detail::inequality_support::build_parametric_intervals(
                    {SymbolicExpr::number(0)}, {2}, sign, type);
            result.cases.push_back({sign_condition(coeffs[2], sign), std::move(solution)});
        }
        return result;
    } catch (const CasError& error) {
        return PiecewiseResult::failure(error);
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
