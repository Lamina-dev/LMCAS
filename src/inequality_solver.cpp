#include "inequality_solver.hpp"
#include "internal/symbolic_ast.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "newton_raphson.hpp"
#include "root_of_utils.hpp"
#include "internal/exact_algebraic.hpp"
#include "internal/inequality_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <functional>
#include <optional>

namespace LMCAS {
using detail::inequality_support::exact_numeric_sign;
using detail::inequality_support::find_roots_with_multiplicity;
using detail::inequality_support::root_less_than;
using detail::inequality_support::roots_equal;
using detail::inequality_support::solve_exact_affine_inequality;
using detail::inequality_support::solve_exact_polynomial_inequality;


namespace {

constexpr const char* kCheckedInequalityOperation = "solve_inequality_checked";


}



Result<IntervalUnion> InequalitySolver::solve_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    ComputationContext& context) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality expression cannot be null",
            kCheckedInequalityOperation);
    }
    if (variable.empty()) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality variable cannot be empty",
            kCheckedInequalityOperation);
    }

    try {
        auto recognized = recognize_rational_polynomial(*expr, variable, context);
        if (!recognized) {
            return Result<IntervalUnion>::failure(recognized.error());
        }
        if (!recognized.value()) {
            return Result<IntervalUnion>::failure(
                CasErrc::Inconclusive,
                "expression is not an exact rational polynomial in the requested variable",
                kCheckedInequalityOperation);
        }
        if (recognized.value()->degree() == 2 &&
            recognized.value()->coeffs.size() >= 3) {
            return solve_exact_quadratic_inequality(
                *recognized.value(), type, context);
        }
        if (recognized.value()->degree() <= 1) {
            return solve_exact_affine_inequality(
                *recognized.value(), type, context);
        }
        return solve_exact_polynomial_inequality(
            *recognized.value(), type, context);
    } catch (const std::bad_alloc&) {
        return Result<IntervalUnion>::failure(
            CasErrc::ResourceLimit, "inequality allocation failed",
            kCheckedInequalityOperation);
    } catch (const std::exception& error) {
        return Result<IntervalUnion>::failure(
            CasErrc::InternalInvariant, error.what(), kCheckedInequalityOperation);
    }
}

Result<IntervalUnion> InequalitySolver::solve_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {
    ComputationContext context;
    return solve_inequality_checked(expr, type, variable, context);
}

Result<IntervalUnion> InequalitySolver::solve_inequalities_checked(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable,
    ComputationContext& context) {
    if (variable.empty()) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality variable cannot be empty",
            "solve_inequalities_checked");
    }

    auto initial_step = context.consume_steps(1, "solve_inequalities_checked");
    if (!initial_step) {
        return Result<IntervalUnion>::failure(initial_step.error());
    }
    if (inequalities.empty()) {
        return Result<IntervalUnion>::success(IntervalUnion::entire_line());
    }

    std::optional<IntervalUnion> aggregate;
    for (const auto& inequality : inequalities) {
        auto solved = solve_inequality_checked(
            inequality.first, inequality.second, variable, context);
        if (!solved) {
            return Result<IntervalUnion>::failure(solved.error());
        }
        if (!aggregate) {
            aggregate = std::move(solved.value());
            continue;
        }
        auto intersection = aggregate->intersect_checked(solved.value(), context);
        if (!intersection) {
            return Result<IntervalUnion>::failure(intersection.error());
        }
        aggregate = std::move(intersection.value());
        if (aggregate->is_empty()) {
            break;
        }
    }
    return Result<IntervalUnion>::success(std::move(*aggregate));
}

Result<IntervalUnion> InequalitySolver::solve_inequalities_checked(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {
    ComputationContext context;
    return solve_inequalities_checked(inequalities, variable, context);
}



static IntervalUnion solve_constant_inequality(
    const Polynomial<SymbolicPolyCoeff>& poly, InequalityType type) {
    auto lc = poly.lead_coeff().val;
    if (!lc) {
        return IntervalUnion::empty();
    }
    int sign = exact_numeric_sign(lc);
    if (sign == 0 && !lc->is_zero()) {
        auto simplified = lc->simplify();
        if (!simplified || !simplified->is_zero()) {
            return IntervalUnion::empty();
        }
    }

    bool satisfies = false;
    switch (type) {
        case InequalityType::GreaterThan: satisfies = (sign > 0); break;
        case InequalityType::GreaterEqual: satisfies = (sign >= 0); break;
        case InequalityType::LessThan: satisfies = (sign < 0); break;
        case InequalityType::LessEqual: satisfies = (sign <= 0); break;
        default:
            return IntervalUnion::empty();
    }
    return satisfies ? IntervalUnion::entire_line() : IntervalUnion::empty();
}

IntervalUnion InequalitySolver::solve_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {

    if (!expr) {
        return IntervalUnion::empty();
    }

    auto converted = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (!converted) throw std::invalid_argument(converted.error().message);
    const auto& poly = converted.value();
    if (poly.is_zero()) {
        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    if (poly.degree() <= 0) {
        return solve_constant_inequality(poly, type);
    }

    {
        auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
        if (!poly_rat) throw std::invalid_argument(poly_rat.error().message);

    }

    auto roots_with_mult = find_roots_with_multiplicity(expr, variable);

    std::sort(roots_with_mult.begin(), roots_with_mult.end(),
        [](const auto& a, const auto& b) {
            return root_less_than(a.first, b.first);
        });

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> unique_roots;
    for (const auto& [root, mult] : roots_with_mult) {
        if (!unique_roots.empty() && roots_equal(unique_roots.back().first, root)) {

            unique_roots.back().second = std::max(unique_roots.back().second, mult);
        } else {
            unique_roots.push_back({root, mult});
        }
    }

    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    std::vector<int> multiplicities;
    for (const auto& [root, mult] : unique_roots) {
        roots.push_back(root);
        multiplicities.push_back(mult);
    }

    auto chart = build_sign_chart(expr, variable, roots, multiplicities);

    return select_intervals(chart, type, roots, multiplicities);
}


IntervalUnion InequalitySolver::solve_inequalities(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {

    if (inequalities.empty()) {
        return IntervalUnion::entire_line();
    }

    IntervalUnion result = solve_inequality(inequalities[0].first, inequalities[0].second, variable);

    for (size_t i = 1; i < inequalities.size(); ++i) {
        auto solution = solve_inequality(inequalities[i].first, inequalities[i].second, variable);
        result = result.intersect(solution);
        if (result.is_empty()) {
            break;
        }
    }

    return result;
}




}
