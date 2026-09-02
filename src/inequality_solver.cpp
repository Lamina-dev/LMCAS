#include "../include/inequality_solver.hpp"
#include "symbolic_ast.hpp"
#include "../include/numeric_evaluation.hpp"
#include "../include/poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_strategies.hpp"
#include "../include/newton_raphson.hpp"
#include "../include/root_of_utils.hpp"
#include "internal/exact_algebraic.hpp"
#include "internal/inequality_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <functional>
#include <optional>

namespace lamina {
using detail::inequality_support::depends_on_any_param;
using detail::inequality_support::determine_leading_sign;
using detail::inequality_support::exact_numeric_sign;
using detail::inequality_support::find_roots_with_multiplicity;
using detail::inequality_support::root_less_than;
using detail::inequality_support::roots_equal;
using detail::inequality_support::solve_exact_affine_inequality;
using detail::inequality_support::solve_exact_polynomial_inequality;
using detail::inequality_support::solve_symbolic_poly;
using detail::inequality_support::try_checked_numeric_constant;


namespace {

constexpr const char* kCheckedInequalityOperation = "solve_inequality_checked";


} // namespace



Result<IntervalUnion> InequalitySolver::solve_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    ComputationContext& context) {
    if (!expr || !lamina::detail::node(expr)) {
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
        if (!recognized) return Result<IntervalUnion>::failure(recognized.error());
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
    if (!initial_step) return Result<IntervalUnion>::failure(initial_step.error());
    if (inequalities.empty()) {
        return Result<IntervalUnion>::success(IntervalUnion::entire_line());
    }

    std::optional<IntervalUnion> aggregate;
    for (const auto& inequality : inequalities) {
        auto solved = solve_inequality_checked(
            inequality.first, inequality.second, variable, context);
        if (!solved) return Result<IntervalUnion>::failure(solved.error());
        if (!aggregate) {
            aggregate = std::move(solved.value());
            continue;
        }
        auto intersection = aggregate->intersect_checked(solved.value(), context);
        if (!intersection) {
            CasError error = intersection.error();
            error.operation = "solve_inequalities_checked";
            return Result<IntervalUnion>::failure(std::move(error));
        }
        aggregate = std::move(intersection.value());
        if (aggregate->is_empty()) break;
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



IntervalUnion InequalitySolver::solve_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {

    if (!expr) return IntervalUnion::empty();

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (poly.is_zero()) {

        if (expression_depends_on_variable(lamina::detail::node(expr), variable)) {
            return IntervalUnion::empty();
        }

        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    if (poly.degree() <= 0) {

        auto lc = poly.lead_coeff().val;
        if (!lc) return IntervalUnion::empty();
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

    {
        auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
        if (poly_rat.is_zero() && expression_depends_on_variable(lamina::detail::node(expr), variable)) {

            return IntervalUnion::empty();
        }

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

IntervalUnion InequalitySolver::solve_rational_inequality(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    InequalityType type,
    const std::string& variable) {

    if (!numerator || !denominator) return IntervalUnion::empty();

    auto num_poly = symbolic_to_poly<SymbolicPolyCoeff>(numerator, variable);
    auto den_poly = symbolic_to_poly<SymbolicPolyCoeff>(denominator, variable);

    if (den_poly.is_zero()) return IntervalUnion::empty();

    if (num_poly.is_zero()) {
        if (type == InequalityType::GreaterThan || type == InequalityType::LessThan) {
            return IntervalUnion::empty();
        }

        auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);
        if (den_roots_with_mult.empty()) {
            return IntervalUnion::entire_line();
        }

        std::vector<std::shared_ptr<SymbolicExpr>> den_roots;
        for (const auto& [root, mult] : den_roots_with_mult) {
            den_roots.push_back(root);
        }
        std::sort(den_roots.begin(), den_roots.end(), root_less_than);

        std::vector<Interval> intervals;

        {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(den_roots[0]);
            intervals.push_back(iv);
        }
        for (size_t i = 0; i + 1 < den_roots.size(); ++i) {
            Interval iv;
            iv.lower = Endpoint::open(den_roots[i]);
            iv.upper = Endpoint::open(den_roots[i + 1]);
            intervals.push_back(iv);
        }

        {
            Interval iv;
            iv.lower = Endpoint::open(den_roots.back());
            iv.upper = Endpoint::pos_inf();
            intervals.push_back(iv);
        }
        return IntervalUnion(intervals);
    }

    {
        auto num_poly_rat = symbolic_to_poly<Rational>(numerator, variable);
        if (num_poly_rat.is_zero() && expression_depends_on_variable(lamina::detail::node(numerator), variable)) {
            return IntervalUnion::empty();
        }

    }

    {
        auto den_poly_rat = symbolic_to_poly<Rational>(denominator, variable);
        if (den_poly_rat.is_zero() && expression_depends_on_variable(lamina::detail::node(denominator), variable)) {
            return IntervalUnion::empty();
        }

    }

    auto num_roots_with_mult = find_roots_with_multiplicity(numerator, variable);

    auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);

    struct CriticalPoint {
        std::shared_ptr<SymbolicExpr> value;
        int num_multiplicity;
        int den_multiplicity;
    };

    std::vector<CriticalPoint> critical_points;

    for (const auto& [root, mult] : num_roots_with_mult) {
        critical_points.push_back({root, mult, 0});
    }
    for (const auto& [root, mult] : den_roots_with_mult) {

        bool found = false;
        for (auto& cp : critical_points) {
            if (roots_equal(cp.value, root)) {
                cp.den_multiplicity = mult;
                found = true;
                break;
            }
        }
        if (!found) {
            critical_points.push_back({root, 0, mult});
        }
    }

    std::sort(critical_points.begin(), critical_points.end(),
        [](const CriticalPoint& a, const CriticalPoint& b) {
            return root_less_than(a.value, b.value);
        });

    int num_leading_sign = determine_leading_sign(num_poly);
    int den_leading_sign = determine_leading_sign(den_poly);
    int combined_leading_sign = num_leading_sign * den_leading_sign;

    size_t n_cp = critical_points.size();
    std::vector<int> interval_signs(n_cp + 1);

    interval_signs[n_cp] = combined_leading_sign;

    for (int i = (int)n_cp - 1; i >= 0; --i) {
        int total_mult = critical_points[i].num_multiplicity + critical_points[i].den_multiplicity;
        interval_signs[i] = interval_signs[i + 1];
        if (total_mult % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    std::vector<Interval> result_intervals;

    if (n_cp == 0) {

        if (interval_signs[0] == target_sign) {
            result_intervals.push_back(Interval::entire_line());
        }
    } else {

        if (interval_signs[0] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(critical_points[0].value);
            result_intervals.push_back(iv);
        }

        for (size_t i = 0; i + 1 < n_cp; ++i) {
            if (interval_signs[i + 1] == target_sign) {
                Interval iv;
                iv.lower = Endpoint::open(critical_points[i].value);
                iv.upper = Endpoint::open(critical_points[i + 1].value);
                result_intervals.push_back(iv);
            }
        }

        if (interval_signs[n_cp] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::open(critical_points[n_cp - 1].value);
            iv.upper = Endpoint::pos_inf();
            result_intervals.push_back(iv);
        }
    }

    if (!is_strict) {
        for (const auto& cp : critical_points) {

            if (cp.num_multiplicity > 0 && cp.den_multiplicity == 0) {

                bool merged = false;
                for (auto& iv : result_intervals) {

                    if (!iv.upper.is_pos_infinity && iv.upper.value) {
                        if (roots_equal(iv.upper.value, cp.value)) {
                            iv.upper.is_open = false;
                            merged = true;
                        }
                    }

                    if (!iv.lower.is_neg_infinity && iv.lower.value) {
                        if (roots_equal(iv.lower.value, cp.value)) {
                            iv.lower.is_open = false;
                            merged = true;
                        }
                    }
                }

                if (!merged) {
                    result_intervals.push_back(Interval::point(cp.value));
                }
            }
        }
    }

    return IntervalUnion(result_intervals);
}

IntervalUnion InequalitySolver::solve_inequalities(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {

    if (inequalities.empty()) return IntervalUnion::entire_line();

    IntervalUnion result = solve_inequality(inequalities[0].first, inequalities[0].second, variable);

    for (size_t i = 1; i < inequalities.size(); ++i) {
        auto solution = solve_inequality(inequalities[i].first, inequalities[i].second, variable);
        result = result.intersect(solution);
        if (result.is_empty()) break;
    }

    return result;
}



PiecewiseIntervalResult InequalitySolver::solve_parametric_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters) {

    PiecewiseIntervalResult result;

    if (!expr) return result;

    if (parameters.empty()) {
        auto solution = solve_inequality(expr, type, variable);
        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;
        single_case.solution = solution;
        result.cases.push_back(single_case);
        return result;
    }

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);

    if (poly.is_zero()) {

        if (expression_depends_on_variable(lamina::detail::node(expr), variable)) {

            return result;
        }

        PiecewiseIntervalResult::Case zero_case;
        zero_case.condition = nullptr;
        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            zero_case.solution = IntervalUnion::entire_line();
        } else {
            zero_case.solution = IntervalUnion::empty();
        }
        result.cases.push_back(zero_case);
        return result;
    }

    int deg = poly.degree();
    auto leading_coeff = poly.coeffs[deg].val;
    if (!leading_coeff) leading_coeff = SymbolicExpr::number(0);
    leading_coeff = leading_coeff->simplify();

    bool lc_depends_on_params = depends_on_any_param(leading_coeff, parameters);

    if (!lc_depends_on_params) {

        int leading_sign = exact_numeric_sign(leading_coeff);
        if (leading_sign == 0 && !leading_coeff->is_zero()) {
            if (auto val = try_checked_numeric_constant(*leading_coeff)) {
                leading_sign = (*val > 0) ? 1 : ((*val < 0) ? -1 : 0);
            }
        }

        auto symbolic_roots = solve_symbolic_poly(poly, variable);

        // Sort roots in ascending order
        /// 当根含参数时,root_less_than 通过差值符号判断排序.
        /// 对于二次公式根 r1=(-b+sqrt(d))/(2a), r2=(-b-sqrt(d))/(2a),
        /// 当 a>0 时 r1>r2,需要交换为 [r2, r1].
        if (symbolic_roots.size() == 2) {
            /// 尝试判断 root[0] > root[1],若是则交换
            bool swapped = false;
            auto d = SymbolicExpr::add(symbolic_roots[0],
                SymbolicExpr::multiply(symbolic_roots[1], SymbolicExpr::number(-1)))->simplify();
            if (!d) {
                std::sort(symbolic_roots.begin(), symbolic_roots.end(), root_less_than);
                swapped = true;
            } else {
            /// 如果差值可以求值为正数,说明 root[0] > root[1],需要交换
            if (auto dv = try_checked_numeric_constant(*d)) {
                if (*dv > 0) { std::swap(symbolic_roots[0], symbolic_roots[1]); swapped = true; }
            } else {
                /// 尝试结构化符号判断
                /// diff 为 (disc)^0.5 形式(PowerNode with exp=0.5)或含 sqrt 的乘积
                auto check_positive = [](const std::shared_ptr<SymbolicExpr>& e) -> bool {
                    if (!e || !lamina::detail::node(e)) return false;
                    /// PowerNode with exponent in (0,1) -> non-negative
                    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(e))) {
                        auto exp_e = lamina::detail::make_expression_ptr(pw->exponent());
                        if (auto ev = try_checked_numeric_constant(*exp_e)) {
                            if (*ev > 0 && *ev < 1.0) return true;
                        }
                    }
                    /// FunctionNode::Sqrt -> non-negative
                    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(e))) {
                        if (fn->type() == FunctionNode::FuncType::Sqrt) return true;
                    }
                    /// MultiplyNode: all factors positive
                    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(e))) {
                        int sign = 1;
                        for (const auto& op : mul->operands()) {
                            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                                if (std::holds_alternative<BigInt>(n->value())) {
                                    if (std::get<BigInt>(n->value()).IsNegative()) sign *= -1;
                                    else if (std::get<BigInt>(n->value()).is_zero()) return false;
                                } else if (std::holds_alternative<Rational>(n->value())) {
                                    if (std::get<Rational>(n->value()).get_numerator().IsNegative()) sign *= -1;
                                    else if (std::get<Rational>(n->value()).get_numerator().is_zero()) return false;
                                } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
                                    if (std::get<lmmc_real_t>(n->value()) < 0) sign *= -1;
                                    else if (std::get<lmmc_real_t>(n->value()) == 0) return false;
                                }
                            } else if (auto pw2 = std::dynamic_pointer_cast<const PowerNode>(op)) {
                                auto exp_e2 = lamina::detail::make_expression_ptr(pw2->exponent());
                                auto ev2 = try_checked_numeric_constant(*exp_e2);
                                if (ev2 && *ev2 > 0 && *ev2 < 1.0) { /* positive, sign unchanged */ }
                                else return false; // can't determine
                            } else if (auto fn2 = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                                if (fn2->type() == FunctionNode::FuncType::Sqrt) { /* positive */ }
                                else return false;
                            } else {
                                return false; // unknown factor
                            }
                        }
                        return sign > 0;
                    }
                    return false;
                };
                auto check_negative = [](const std::shared_ptr<SymbolicExpr>& e) -> bool {
                    if (!e || !lamina::detail::node(e)) return false;
                    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(e))) {
                        int sign = 1;
                        for (const auto& op : mul->operands()) {
                            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                                if (std::holds_alternative<BigInt>(n->value())) {
                                    if (std::get<BigInt>(n->value()).IsNegative()) sign *= -1;
                                    else if (std::get<BigInt>(n->value()).is_zero()) return false;
                                } else if (std::holds_alternative<Rational>(n->value())) {
                                    if (std::get<Rational>(n->value()).get_numerator().IsNegative()) sign *= -1;
                                    else if (std::get<Rational>(n->value()).get_numerator().is_zero()) return false;
                                } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
                                    if (std::get<lmmc_real_t>(n->value()) < 0) sign *= -1;
                                    else if (std::get<lmmc_real_t>(n->value()) == 0) return false;
                                }
                            } else if (auto pw2 = std::dynamic_pointer_cast<const PowerNode>(op)) {
                                auto exp_e2 = lamina::detail::make_expression_ptr(pw2->exponent());
                                auto ev2 = try_checked_numeric_constant(*exp_e2);
                                if (ev2 && *ev2 > 0 && *ev2 < 1.0) { /* positive */ }
                                else return false;
                            } else if (auto fn2 = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                                if (fn2->type() == FunctionNode::FuncType::Sqrt) { /* positive */ }
                                else return false;
                            } else {
                                return false;
                            }
                        }
                        return sign < 0;
                    }
                    return false;
                };
                if (check_positive(d)) {
                    std::swap(symbolic_roots[0], symbolic_roots[1]);
                    swapped = true;
                } else if (check_negative(d)) {
                    /// already in correct order
                    swapped = false;
                }
            }
            }
            if (!swapped) {
                /// Fallback: 对于 a>0 的二次多项式,solve_quadratic_internal 返回
                /// [大根, 小根],需要交换.对于 a<0 则已经是 [小根, 大根].
                /// 这里利用 leading_sign 直接判断.
                if (leading_sign > 0 && deg == 2) {
                    std::swap(symbolic_roots[0], symbolic_roots[1]);
                }
            }
        } else {
            std::sort(symbolic_roots.begin(), symbolic_roots.end(), root_less_than);
        }

        std::vector<int> multiplicities(symbolic_roots.size(), 1);

        auto solution = build_parametric_solution(symbolic_roots, multiplicities, leading_sign, type);

        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;
        single_case.solution = solution;
        result.cases.push_back(single_case);
    } else {

        {
            PiecewiseIntervalResult::Case pos_case;
            pos_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
                    RelationalNode::Op::GT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);

            // Sort roots in ascending order and apply same permutation to multiplicities
            std::vector<size_t> indices(symbolic_roots.size());
            for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
            std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                return root_less_than(symbolic_roots[a], symbolic_roots[b]);
            });
            std::vector<std::shared_ptr<SymbolicExpr>> sorted_roots;
            std::vector<int> sorted_mults;
            for (size_t idx : indices) {
                sorted_roots.push_back(symbolic_roots[idx]);
                sorted_mults.push_back(multiplicities[idx]);
            }

            pos_case.solution = build_parametric_solution(sorted_roots, sorted_mults, 1, type);
            result.cases.push_back(pos_case);
        }

        {
            PiecewiseIntervalResult::Case neg_case;
            neg_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
                    RelationalNode::Op::LT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);

            // Sort roots in ascending order and apply same permutation to multiplicities
            std::vector<size_t> indices(symbolic_roots.size());
            for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
            std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                return root_less_than(symbolic_roots[a], symbolic_roots[b]);
            });
            std::vector<std::shared_ptr<SymbolicExpr>> sorted_roots;
            std::vector<int> sorted_mults;
            for (size_t idx : indices) {
                sorted_roots.push_back(symbolic_roots[idx]);
                sorted_mults.push_back(multiplicities[idx]);
            }

            neg_case.solution = build_parametric_solution(sorted_roots, sorted_mults, -1, type);
            result.cases.push_back(neg_case);
        }

        {
            PiecewiseIntervalResult::Case degen_case;
            degen_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
                    RelationalNode::Op::EQ));

            /// 降阶多项式仍按参数分段时,expanded_into_subcases 标记各子分支
            /// 已直接加入 result.cases.
            bool expanded_into_subcases = false;

            if (deg >= 1) {
                std::vector<SymbolicPolyCoeff> reduced_coeffs;
                for (int i = 0; i < deg; ++i) {
                    reduced_coeffs.push_back(poly.coeffs[i]);
                }
                Polynomial<SymbolicPolyCoeff> reduced_poly(reduced_coeffs, variable);

                if (reduced_poly.is_zero()) {

                    if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
                        degen_case.solution = IntervalUnion::entire_line();
                    } else {
                        degen_case.solution = IntervalUnion::empty();
                    }
                } else {

                    auto new_lc = reduced_poly.lead_coeff().val;
                    if (new_lc) new_lc = new_lc->simplify();

                    if (new_lc && depends_on_any_param(new_lc, parameters)) {

                        auto reduced_expr = SymbolicExpr::number(0);
                        auto var_expr = SymbolicExpr::variable(variable);
                        for (int i = reduced_poly.degree(); i >= 0; --i) {
                            auto coeff_val = reduced_poly.coeffs[i].val;
                            if (!coeff_val) continue;
                            if (i == 0) {
                                reduced_expr = SymbolicExpr::add(reduced_expr, coeff_val);
                            } else if (i == 1) {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val, var_expr));
                            } else {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val,
                                        SymbolicExpr::power(var_expr, SymbolicExpr::number(i))));
                            }
                        }
                        reduced_expr = reduced_expr->simplify();
                        auto sub_result = solve_parametric_inequality(reduced_expr, type, variable, parameters);

                        /// 将每个子分支的参数条件与父级 leading coeff == 0 条件合取,
                        /// 再作为独立分支加入结果,完整保留降阶后的参数分段.
                        if (sub_result.cases.empty()) {
                            degen_case.solution = IntervalUnion::empty();
                        } else {
                            for (const auto& sub_case : sub_result.cases) {
                                PiecewiseIntervalResult::Case merged;
                                if (sub_case.condition) {
                                    merged.condition = lamina::detail::make_expression_ptr(
                                        lamina::detail::make_node<LogicalNode>(
                                            lamina::detail::node(degen_case.condition),
                                            lamina::detail::node(sub_case.condition),
                                            LogicalNode::Op::And));
                                } else {
                                    merged.condition = degen_case.condition;
                                }
                                merged.solution = sub_case.solution;
                                result.cases.push_back(merged);
                            }
                            expanded_into_subcases = true;
                        }
                    } else {

                        int reduced_leading_sign = 0;
                        if (new_lc) {
                            reduced_leading_sign = exact_numeric_sign(new_lc);
                            if (reduced_leading_sign == 0 && !new_lc->is_zero()) {
                                if (auto val = try_checked_numeric_constant(*new_lc)) {
                                    reduced_leading_sign = (*val > 0) ? 1 : ((*val < 0) ? -1 : 0);
                                }
                            }
                        }

                        auto symbolic_roots = solve_symbolic_poly(reduced_poly, variable);
                        std::vector<int> multiplicities(symbolic_roots.size(), 1);
                        degen_case.solution = build_parametric_solution(
                            symbolic_roots, multiplicities, reduced_leading_sign, type);
                    }
                }
            } else {
                degen_case.solution = IntervalUnion::empty();
            }

            if (!expanded_into_subcases) {
                result.cases.push_back(degen_case);
            }
        }
    }

    return result;
}

}
