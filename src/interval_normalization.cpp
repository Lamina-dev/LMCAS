#include "internal/interval_endpoint.hpp"
#include <utility>

namespace LMCAS {
using detail::CheckedInterval;
using detail::comparable_endpoint;
using detail::compare_comparable;

namespace {

constexpr const char* kCheckedIntervalOperation = "normalize_intervals";

Result<std::vector<CheckedInterval>> collect_nonempty_intervals(
    std::vector<Interval> intervals, ComputationContext& context) {
    std::vector<CheckedInterval> checked;
    checked.reserve(intervals.size());
    for (auto& interval : intervals) {
        auto lower = comparable_endpoint(
            interval.lower, context, kCheckedIntervalOperation);
        if (!lower) return Result<std::vector<CheckedInterval>>::failure(lower.error());
        auto upper = comparable_endpoint(
            interval.upper, context, kCheckedIntervalOperation);
        if (!upper) return Result<std::vector<CheckedInterval>>::failure(upper.error());
        auto comparison =
            compare_comparable(lower.value(), upper.value(), context);
        if (!comparison) {
            return Result<std::vector<CheckedInterval>>::failure(comparison.error());
        }
        const bool empty = comparison.value() > 0 ||
            (comparison.value() == 0 &&
             (interval.lower.is_open || interval.upper.is_open));
        if (!empty) {
            checked.push_back(CheckedInterval{
                std::move(interval), std::move(lower.value()), std::move(upper.value())});
        }
    }
    return Result<std::vector<CheckedInterval>>::success(std::move(checked));
}

Result<void> sort_checked_intervals(
    std::vector<CheckedInterval>& checked, ComputationContext& context) {
    for (std::size_t i = 1; i < checked.size(); ++i) {
        std::size_t position = i;
        while (position > 0) {
            auto comparison = compare_comparable(
                checked[position - 1].lower, checked[position].lower, context);
            if (!comparison) {
                return Result<void>::failure(comparison.error());
            }
            const bool out_of_order = comparison.value() > 0 ||
                (comparison.value() == 0 &&
                 checked[position - 1].interval.lower.is_open &&
                 !checked[position].interval.lower.is_open);
            if (!out_of_order) break;
            std::swap(checked[position - 1], checked[position]);
            --position;
        }
    }
    return Result<void>::success();
}

Result<std::vector<Interval>> merge_checked_intervals(
    std::vector<CheckedInterval> checked, ComputationContext& context) {
    std::vector<CheckedInterval> merged;
    merged.reserve(checked.size());
    for (auto& next : checked) {
        if (merged.empty()) {
            merged.push_back(std::move(next));
            continue;
        }

        CheckedInterval& current = merged.back();
        auto boundary =
            compare_comparable(current.upper, next.lower, context);
        if (!boundary) {
            return Result<std::vector<Interval>>::failure(boundary.error());
        }
        const bool overlaps = boundary.value() > 0 ||
            (boundary.value() == 0 &&
             (!current.interval.upper.is_open || !next.interval.lower.is_open));
        if (!overlaps) {
            merged.push_back(std::move(next));
            continue;
        }

        auto upper_comparison =
            compare_comparable(current.upper, next.upper, context);
        if (!upper_comparison) {
            return Result<std::vector<Interval>>::failure(upper_comparison.error());
        }
        if (upper_comparison.value() < 0 ||
            (upper_comparison.value() == 0 && current.interval.upper.is_open &&
             !next.interval.upper.is_open)) {
            current.interval.upper = std::move(next.interval.upper);
            current.upper = std::move(next.upper);
        }
    }

    std::vector<Interval> result;
    result.reserve(merged.size());
    for (auto& interval : merged) {
        result.push_back(std::move(interval.interval));
    }
    return Result<std::vector<Interval>>::success(std::move(result));
}

}

Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals,
    ComputationContext& context) {
    auto checked = collect_nonempty_intervals(std::move(intervals), context);
    if (!checked) return Result<std::vector<Interval>>::failure(checked.error());
    auto sorted = sort_checked_intervals(checked.value(), context);
    if (!sorted) return Result<std::vector<Interval>>::failure(sorted.error());
    return merge_checked_intervals(std::move(checked.value()), context);
}

Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals) {
    ComputationContext context;
    return normalize_intervals_checked(std::move(intervals), context);
}

}
