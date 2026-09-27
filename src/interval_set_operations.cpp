#include "internal/interval_endpoint.hpp"
#include <utility>

namespace LMCAS {
using detail::CheckedInterval;
using detail::comparable_endpoint;
using detail::compare_comparable;

namespace {

Result<std::vector<CheckedInterval>> checked_interval_views(
    const std::vector<Interval>& intervals,
    ComputationContext& context,
    const std::string& operation) {
    auto normalized = normalize_intervals_checked(intervals, context);
    if (!normalized) {
        return Result<std::vector<CheckedInterval>>::failure(
            normalized.error().code, normalized.error().message, operation);
    }

    std::vector<CheckedInterval> checked;
    auto normalized_intervals = std::move(normalized.value());
    checked.reserve(normalized_intervals.size());
    for (auto& interval : normalized_intervals) {
        auto lower = comparable_endpoint(interval.lower, context, operation);
        if (!lower) {
            return Result<std::vector<CheckedInterval>>::failure(lower.error());
        }
        auto upper = comparable_endpoint(interval.upper, context, operation);
        if (!upper) {
            return Result<std::vector<CheckedInterval>>::failure(upper.error());
        }
        checked.push_back(CheckedInterval{
            std::move(interval), std::move(lower.value()), std::move(upper.value())});
    }
    return Result<std::vector<CheckedInterval>>::success(std::move(checked));
}

Endpoint complement_lower_from_upper(const Endpoint& upper) {
    Endpoint lower;
    lower.value = upper.value;
    lower.is_open = !upper.is_open;
    lower.is_neg_infinity = false;
    lower.is_pos_infinity = false;
    return lower;
}

Endpoint complement_upper_from_lower(const Endpoint& lower) {
    Endpoint upper;
    upper.value = lower.value;
    upper.is_open = !lower.is_open;
    upper.is_neg_infinity = false;
    upper.is_pos_infinity = false;
    return upper;
}

Interval intersection_candidate(const CheckedInterval& a, const CheckedInterval& b,
                                int lower_order, int upper_order) {
    Interval candidate{
        lower_order >= 0 ? a.interval.lower : b.interval.lower,
        upper_order <= 0 ? a.interval.upper : b.interval.upper
    };
    if (lower_order == 0) {
        candidate.lower.is_open = a.interval.lower.is_open || b.interval.lower.is_open;
    }
    if (upper_order == 0) {
        candidate.upper.is_open = a.interval.upper.is_open || b.interval.upper.is_open;
    }
    return candidate;
}

}

Result<IntervalUnion> IntervalUnion::intersect_checked(
    const IntervalUnion& other,
    ComputationContext& context) const {
    constexpr const char* operation = "interval_union_intersect";
    auto left = checked_interval_views(intervals_, context, operation);
    if (!left) {
        return Result<IntervalUnion>::failure(left.error());
    }
    auto right = checked_interval_views(other.intervals_, context, operation);
    if (!right) {
        return Result<IntervalUnion>::failure(right.error());
    }

    const auto& a_intervals = left.value();
    const auto& b_intervals = right.value();
    std::vector<Interval> result;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a_intervals.size() && j < b_intervals.size()) {
        const auto& a = a_intervals[i];
        const auto& b = b_intervals[j];
        auto lower_order = compare_comparable(a.lower, b.lower, context);
        if (!lower_order) {
            return Result<IntervalUnion>::failure(lower_order.error());
        }
        auto upper_order = compare_comparable(a.upper, b.upper, context);
        if (!upper_order) {
            return Result<IntervalUnion>::failure(upper_order.error());
        }
        auto candidate = intersection_candidate(
            a, b, lower_order.value(), upper_order.value());
        auto empty = interval_is_empty_checked(candidate, context);
        if (!empty) {
            return Result<IntervalUnion>::failure(
                empty.error().code, empty.error().message, operation);
        }
        if (!empty.value()) {
            result.push_back(std::move(candidate));
        }

        if (upper_order.value() < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    auto normalized = normalize_intervals_checked(std::move(result), context);
    if (!normalized) {
        return Result<IntervalUnion>::failure(
            normalized.error().code, normalized.error().message, operation);
    }
    return Result<IntervalUnion>::success(
        IntervalUnion::from_checked_normalized(std::move(normalized.value())));
}

Result<IntervalUnion> IntervalUnion::intersect_checked(
    const IntervalUnion& other) const {
    ComputationContext context;
    return intersect_checked(other, context);
}

Result<IntervalUnion> IntervalUnion::unite_checked(
    const IntervalUnion& other,
    ComputationContext& context) const {
    constexpr const char* operation = "interval_union_unite";
    std::vector<Interval> all;
    all.reserve(intervals_.size() + other.intervals_.size());
    all.insert(all.end(), intervals_.begin(), intervals_.end());
    all.insert(all.end(), other.intervals_.begin(), other.intervals_.end());

    auto step = context.consume_steps(all.size(), operation);
    if (!step) {
        return Result<IntervalUnion>::failure(step.error());
    }
    auto normalized = normalize_intervals_checked(std::move(all), context);
    if (!normalized) {
        return Result<IntervalUnion>::failure(normalized.error());
    }
    return Result<IntervalUnion>::success(
        IntervalUnion::from_checked_normalized(std::move(normalized.value())));
}

Result<IntervalUnion> IntervalUnion::unite_checked(
    const IntervalUnion& other) const {
    ComputationContext context;
    return unite_checked(other, context);
}

Result<IntervalUnion> IntervalUnion::complement_checked(
    ComputationContext& context) const {
    constexpr const char* operation = "interval_union_complement";
    auto checked = checked_interval_views(intervals_, context, operation);
    if (!checked) {
        return Result<IntervalUnion>::failure(checked.error());
    }
    const auto& intervals = checked.value();

    if (intervals.empty()) {
        return Result<IntervalUnion>::success(
            IntervalUnion::from_checked_normalized({Interval::entire_line()}));
    }
    if (intervals.size() == 1 &&
        intervals[0].interval.lower.is_neg_infinity &&
        intervals[0].interval.upper.is_pos_infinity) {
        return Result<IntervalUnion>::success(IntervalUnion::from_checked_normalized({}));
    }

    std::vector<Interval> result;
    const Interval& first = intervals.front().interval;
    if (!first.lower.is_neg_infinity) {
        result.push_back(Interval{
            Endpoint::neg_inf(),
            complement_upper_from_lower(first.lower)
        });
    }

    for (std::size_t i = 0; i + 1 < intervals.size(); ++i) {
        const Interval& current = intervals[i].interval;
        const Interval& next = intervals[i + 1].interval;
        result.push_back(Interval{
            complement_lower_from_upper(current.upper),
            complement_upper_from_lower(next.lower)
        });
    }

    const Interval& last = intervals.back().interval;
    if (!last.upper.is_pos_infinity) {
        result.push_back(Interval{
            complement_lower_from_upper(last.upper),
            Endpoint::pos_inf()
        });
    }

    auto normalized = normalize_intervals_checked(std::move(result), context);
    if (!normalized) {
        return Result<IntervalUnion>::failure(normalized.error());
    }
    return Result<IntervalUnion>::success(
        IntervalUnion::from_checked_normalized(std::move(normalized.value())));
}

Result<IntervalUnion> IntervalUnion::complement_checked() const {
    ComputationContext context;
    return complement_checked(context);
}

}
