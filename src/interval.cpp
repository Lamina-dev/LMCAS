#include "internal/interval_endpoint.hpp"
#include <utility>

namespace LMCAS {
using detail::ComparableEndpoint;
using detail::comparable_endpoint;
using detail::compare_comparable;
using detail::exact_double_rational;

Endpoint Endpoint::neg_inf() {
    return Endpoint{nullptr, true, true, false};
}

Endpoint Endpoint::pos_inf() {
    return Endpoint{nullptr, true, false, true};
}

Endpoint Endpoint::closed(std::shared_ptr<SymbolicExpr> val) {
    return Endpoint{std::move(val), false, false, false};
}

Endpoint Endpoint::open(std::shared_ptr<SymbolicExpr> val) {
    return Endpoint{std::move(val), true, false, false};
}

bool Interval::contains(double value) const {
    auto result = interval_contains_checked(*this, value);
    return result ? result.value() : false;
}

bool Interval::is_empty() const {
    auto result = interval_is_empty_checked(*this);
    return result ? result.value() : false;
}

bool Interval::is_entire_line() const {
    return lower.is_neg_infinity && upper.is_pos_infinity;
}

Interval Interval::empty() {

    return Interval{Endpoint::open(SymbolicExpr::number(1)), Endpoint::open(SymbolicExpr::number(0))};
}

Interval Interval::entire_line() {
    return Interval{Endpoint::neg_inf(), Endpoint::pos_inf()};
}

Interval Interval::point(std::shared_ptr<SymbolicExpr> val) {
    return Interval{Endpoint::closed(val), Endpoint::closed(val)};
}

IntervalUnion::IntervalUnion() : intervals_() {}

IntervalUnion::IntervalUnion(std::vector<Interval> intervals) : intervals_(std::move(intervals)) {
    normalize();
}

IntervalUnion IntervalUnion::from_single(const Interval& iv) {
    return IntervalUnion(std::vector<Interval>{iv});
}

Result<IntervalUnion> IntervalUnion::from_intervals_checked(
    std::vector<Interval> intervals,
    ComputationContext& context) {
    auto normalized = normalize_intervals_checked(std::move(intervals), context);
    if (!normalized) {
        return Result<IntervalUnion>::failure(normalized.error());
    }
    return Result<IntervalUnion>::success(
        from_checked_normalized(std::move(normalized.value())));
}

Result<IntervalUnion> IntervalUnion::from_intervals_checked(
    std::vector<Interval> intervals) {
    ComputationContext context;
    return from_intervals_checked(std::move(intervals), context);
}

IntervalUnion IntervalUnion::empty() {
    return IntervalUnion();
}

IntervalUnion IntervalUnion::entire_line() {
    return IntervalUnion::from_single(Interval::entire_line());
}

bool IntervalUnion::contains(double value) const {
    for (const auto& iv : intervals_) {
        if (iv.contains(value)) {
            return true;
        }
    }
    return false;
}

bool IntervalUnion::is_empty() const {
    return intervals_.empty();
}

bool IntervalUnion::is_entire_line() const {
    return intervals_.size() == 1 && intervals_[0].is_entire_line();
}

const std::vector<Interval>& IntervalUnion::intervals() const {
    return intervals_;
}

IntervalUnion IntervalUnion::from_checked_normalized(std::vector<Interval> intervals) {
    IntervalUnion result;
    result.intervals_ = std::move(intervals);
    return result;
}

void IntervalUnion::normalize() {
    auto normalized = normalize_intervals_checked(intervals_);
    if (normalized) {
        intervals_ = std::move(normalized.value());
    }
}

Result<bool> interval_contains_checked(
    const Interval& interval,
    double value,
    ComputationContext& context) {
    constexpr const char* operation = "interval_contains";
    auto point_step = context.consume_steps(1, operation);
    if (!point_step) {
        return Result<bool>::failure(point_step.error());
    }
    auto point = exact_double_rational(value, context, operation);
    if (!point) {
        return Result<bool>::failure(point.error());
    }
    const ComparableEndpoint point_key{
        0, std::move(point.value()), Rational(0), Rational(0)};

    auto lower = comparable_endpoint(interval.lower, context, operation);
    if (!lower) {
        return Result<bool>::failure(lower.error());
    }
    auto upper = comparable_endpoint(interval.upper, context, operation);
    if (!upper) {
        return Result<bool>::failure(upper.error());
    }

    auto lower_cmp = compare_comparable(point_key, lower.value(), context);
    if (!lower_cmp) {
        return Result<bool>::failure(lower_cmp.error());
    }
    if (lower_cmp.value() < 0 ||
        (lower_cmp.value() == 0 && interval.lower.is_open)) {
        return Result<bool>::success(false);
    }
    auto upper_cmp = compare_comparable(point_key, upper.value(), context);
    if (!upper_cmp) {
        return Result<bool>::failure(upper_cmp.error());
    }
    if (upper_cmp.value() > 0 ||
        (upper_cmp.value() == 0 && interval.upper.is_open)) {
        return Result<bool>::success(false);
    }
    return Result<bool>::success(true);
}

Result<bool> interval_contains_checked(
    const Interval& interval,
    double value) {
    ComputationContext context;
    return interval_contains_checked(interval, value, context);
}

Result<bool> interval_is_empty_checked(
    const Interval& interval,
    ComputationContext& context) {
    constexpr const char* operation = "interval_is_empty";
    auto lower = comparable_endpoint(interval.lower, context, operation);
    if (!lower) {
        return Result<bool>::failure(lower.error());
    }
    auto upper = comparable_endpoint(interval.upper, context, operation);
    if (!upper) {
        return Result<bool>::failure(upper.error());
    }
    auto comparison = compare_comparable(lower.value(), upper.value(), context);
    if (!comparison) {
        return Result<bool>::failure(comparison.error());
    }
    return Result<bool>::success(
        comparison.value() > 0 ||
        (comparison.value() == 0 &&
         (interval.lower.is_open || interval.upper.is_open)));
}

Result<bool> interval_is_empty_checked(const Interval& interval) {
    ComputationContext context;
    return interval_is_empty_checked(interval, context);
}

IntervalUnion IntervalUnion::intersect(const IntervalUnion& other) const {
    auto result = intersect_checked(other);
    return result ? result.value() : IntervalUnion::empty();
}

IntervalUnion IntervalUnion::unite(const IntervalUnion& other) const {
    auto result = unite_checked(other);
    return result ? result.value() : IntervalUnion::empty();
}

IntervalUnion IntervalUnion::complement() const {
    auto result = complement_checked();
    return result ? result.value() : IntervalUnion::empty();
}

}
