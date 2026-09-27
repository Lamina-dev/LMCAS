#include "internal/inequality_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>

namespace LMCAS {
using namespace detail::inequality_support;
namespace {

struct CriticalPoint {
    std::shared_ptr<SymbolicExpr> value;
    int num_multiplicity;
    int den_multiplicity;
};

IntervalUnion zero_numerator_domain(const std::shared_ptr<SymbolicExpr>& denominator,
                                    InequalityType type, const std::string& variable) {
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

std::vector<CriticalPoint> rational_critical_points(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    const std::string& variable) {
    auto num_roots_with_mult = find_roots_with_multiplicity(numerator, variable);

    auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);
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
    return critical_points;
}

std::vector<Interval> rational_sign_intervals(
    const std::vector<CriticalPoint>& critical_points,
    int combined_leading_sign, int target_sign) {
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
    std::vector<Interval> result_intervals;
    for (size_t i = 0; i <= n_cp; ++i) {
        if (interval_signs[i] != target_sign) {
            continue;
        }
        Interval interval;
        interval.lower = i == 0 ? Endpoint::neg_inf() : Endpoint::open(critical_points[i - 1].value);
        interval.upper = i == n_cp ? Endpoint::pos_inf() : Endpoint::open(critical_points[i].value);
        result_intervals.push_back(interval);
    }
    return result_intervals;
}

bool has_unrecognized_variable_polynomial(
    const std::shared_ptr<SymbolicExpr>& expression, const std::string& variable) {
    auto polynomial = symbolic_to_poly<Rational>(expression, variable);
    return !polynomial;
}

void include_rational_zeros(
    std::vector<Interval>& intervals, const std::vector<CriticalPoint>& points) {
    for (const auto& point : points) {
        if (point.num_multiplicity > 0 && point.den_multiplicity == 0) {
            include_numeric_root(intervals, point.value);
        }
    }
}

}

IntervalUnion InequalitySolver::solve_rational_inequality(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    InequalityType type,
    const std::string& variable) {

    if (!numerator || !denominator) {
        return IntervalUnion::empty();
    }

    auto num_conversion = symbolic_to_poly<SymbolicPolyCoeff>(numerator, variable);
    auto den_conversion = symbolic_to_poly<SymbolicPolyCoeff>(denominator, variable);
    if (!num_conversion) throw std::invalid_argument(num_conversion.error().message);
    if (!den_conversion) throw std::invalid_argument(den_conversion.error().message);
    const auto& num_poly = num_conversion.value();
    const auto& den_poly = den_conversion.value();

    if (den_poly.is_zero()) {
        return IntervalUnion::empty();
    }

    if (num_poly.is_zero()) {
        return zero_numerator_domain(denominator, type, variable);
    }
    if (has_unrecognized_variable_polynomial(numerator, variable)) {
        throw std::invalid_argument("Numerator is not a supported rational polynomial");
    }
    if (has_unrecognized_variable_polynomial(denominator, variable)) {
        throw std::invalid_argument("Denominator is not a supported rational polynomial");
    }
    auto critical_points = rational_critical_points(numerator, denominator, variable);
    int num_leading_sign = determine_leading_sign(num_poly);
    int den_leading_sign = determine_leading_sign(den_poly);
    int combined_leading_sign = num_leading_sign * den_leading_sign;
    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;
    auto result_intervals = rational_sign_intervals(critical_points, combined_leading_sign, target_sign);
    if (!is_strict) {
        include_rational_zeros(result_intervals, critical_points);
    }
    return IntervalUnion(result_intervals);
}

}
