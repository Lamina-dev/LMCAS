#include "internal/inequality_solver_support.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>

namespace LMCAS {
using namespace detail::inequality_support;

namespace detail::inequality_support {

void include_numeric_root(std::vector<Interval>& intervals,
                          const std::shared_ptr<SymbolicExpr>& root) {
    bool merged = false;
    for (auto& interval : intervals) {
        if (!interval.upper.is_pos_infinity && interval.upper.value &&
            roots_equal(interval.upper.value, root)) {
            interval.upper.is_open = false;
            merged = true;
        }
        if (!interval.lower.is_neg_infinity && interval.lower.value &&
            roots_equal(interval.lower.value, root)) {
            interval.lower.is_open = false;
            merged = true;
        }
    }
    if (!merged) {
        intervals.push_back(Interval::point(root));
    }
}

static bool close_symbolic_endpoint(Endpoint& endpoint,
                                   const std::shared_ptr<SymbolicExpr>& root) {
    if (detail::node(endpoint.value)->equals(*detail::node(root))) {
        endpoint.is_open = false;
        return true;
    }
    auto diff = SymbolicExpr::add(endpoint.value,
        SymbolicExpr::multiply(root, SymbolicExpr::number(-1)));
    if (!diff->simplify()->is_zero()) {
        return false;
    }
    endpoint.is_open = false;
    return true;
}

static void include_symbolic_root(std::vector<Interval>& intervals,
                                 const std::shared_ptr<SymbolicExpr>& root) {
    bool merged = false;
    for (auto& interval : intervals) {
        if (!interval.upper.is_pos_infinity && interval.upper.value) {
            if (close_symbolic_endpoint(interval.upper, root)) {
                merged = true;
            }
        }
        if (!interval.lower.is_neg_infinity && interval.lower.value) {
            if (close_symbolic_endpoint(interval.lower, root)) {
                merged = true;
            }
        }
    }
    if (!merged) {
        intervals.push_back(Interval::point(root));
    }
}

static std::vector<Interval> select_open_intervals(
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& signs, int target_sign) {
    std::vector<Interval> intervals;
    for (size_t index = 0; index < signs.size(); ++index) {
        if (signs[index] != target_sign) {
            continue;
        }
        Interval interval;
        interval.lower = index == 0 ? Endpoint::neg_inf() : Endpoint::open(roots[index - 1]);
        interval.upper = index == roots.size() ? Endpoint::pos_inf() : Endpoint::open(roots[index]);
        intervals.push_back(interval);
    }
    return intervals;
}

IntervalUnion build_parametric_intervals(
    const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
    const std::vector<int>& multiplicities,
    int leading_sign,
    InequalityType type) {

    if (symbolic_roots.empty()) {

        bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
        int target_sign = want_positive ? 1 : -1;
        if (leading_sign == target_sign) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    size_t n = symbolic_roots.size();
    std::vector<int> interval_signs(n + 1);

    interval_signs[n] = leading_sign;

    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    auto result_intervals = select_open_intervals(symbolic_roots, interval_signs, target_sign);

    if (!is_strict) {
        for (const auto& root : symbolic_roots) include_symbolic_root(result_intervals, root);
    }
    return IntervalUnion(result_intervals);
}

}

std::vector<SignChartEntry> InequalitySolver::build_sign_chart(
    const std::shared_ptr<SymbolicExpr>& poly,
    const std::string& variable,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities) {

    std::vector<SignChartEntry> chart;
    auto converted = symbolic_to_poly<SymbolicPolyCoeff>(poly, variable);
    if (!converted) throw std::invalid_argument(converted.error().message);
    const auto& p = converted.value();

    if (roots.empty()) {

        int sign = determine_leading_sign(p);
        chart.push_back({Interval::entire_line(), sign});
        return chart;
    }

    int leading_sign = determine_leading_sign(p);

    size_t n = roots.size();
    std::vector<int> interval_signs(n + 1);

    interval_signs[n] = leading_sign;

    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    {
        Interval iv;
        iv.lower = Endpoint::neg_inf();
        iv.upper = Endpoint::open(roots[0]);
        chart.push_back({iv, interval_signs[0]});
    }

    for (size_t i = 0; i + 1 < n; ++i) {
        Interval iv;
        iv.lower = Endpoint::open(roots[i]);
        iv.upper = Endpoint::open(roots[i + 1]);
        chart.push_back({iv, interval_signs[i + 1]});
    }

    {
        Interval iv;
        iv.lower = Endpoint::open(roots[n - 1]);
        iv.upper = Endpoint::pos_inf();
        chart.push_back({iv, interval_signs[n]});
    }

    return chart;
}

IntervalUnion InequalitySolver::select_intervals(
    const std::vector<SignChartEntry>& chart,
    InequalityType type,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>&) {

    std::vector<Interval> result_intervals;

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    for (const auto& entry : chart) {
        if (entry.sign == target_sign) {
            result_intervals.push_back(entry.interval);
        }
    }

    if (!is_strict) {
        for (const auto& root : roots) include_numeric_root(result_intervals, root);
    }
    return IntervalUnion(result_intervals);
}

IntervalUnion InequalitySolver::build_parametric_solution(
    const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
    const std::vector<int>& multiplicities,
    int leading_sign,
    InequalityType type) {
    return build_parametric_intervals(symbolic_roots, multiplicities, leading_sign, type);
}

}
