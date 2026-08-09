#include "../include/interval.hpp"
#include "../include/symbolic_ast.hpp"
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <cctype>

namespace lamina {

namespace {

constexpr const char* kCheckedIntervalOperation = "normalize_intervals";

struct ComparableEndpoint {
    int infinity = 0;
    Rational rational{};
    Rational radical_coefficient{};
    Rational radicand{};
};

struct CheckedInterval {
    Interval interval;
    ComparableEndpoint lower;
    ComparableEndpoint upper;
};

Result<Rational> exact_double_rational(double value,
                                       ComputationContext& context,
                                       const std::string& operation) {
    if (!std::isfinite(value)) {
        return Result<Rational>::failure(
            CasErrc::NumericFailure,
            "finite interval endpoint evaluated to NaN or infinity",
            operation);
    }
    if (value == 0.0) return Result<Rational>::success(Rational(0));

    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    const auto mantissa = static_cast<std::int64_t>(std::ldexp(fraction, 53));
    const int binary_exponent = exponent - 53;
    const std::size_t required_bits = binary_exponent < 0
        ? static_cast<std::size_t>(-binary_exponent) + 1
        : static_cast<std::size_t>(binary_exponent) + 54;
    if (required_bits > context.limits().max_integer_bits) {
        return Result<Rational>::failure(
            CasErrc::ResourceLimit,
            "IEEE endpoint conversion exceeds the integer bit budget",
            operation);
    }

    BigInt numerator(static_cast<long long>(mantissa));
    BigInt denominator(1);
    if (binary_exponent >= 0) {
        numerator <<= static_cast<mp_size_t>(binary_exponent);
    } else {
        denominator <<= static_cast<mp_size_t>(-binary_exponent);
    }
    return Result<Rational>::success(Rational(numerator, denominator));
}

std::optional<Rational> exact_number_value(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;
    if (std::holds_alternative<BigInt>(number->value())) {
        return Rational(std::get<BigInt>(number->value()));
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value());
    }
    return std::nullopt;
}

std::optional<ComparableEndpoint> parse_quadratic_surd(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return std::nullopt;
    if (auto number = exact_number_value(node)) {
        return ComparableEndpoint{0, *number, Rational(0), Rational(0)};
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() != FunctionNode::FuncType::Sqrt ||
            function->arguments().size() != 1) {
            return std::nullopt;
        }
        auto radicand = exact_number_value(function->arguments()[0]);
        if (!radicand || *radicand <= Rational(0)) return std::nullopt;
        return ComparableEndpoint{0, Rational(0), Rational(1), *radicand};
    }

    auto combine_add = [](const ComparableEndpoint& left,
                          const ComparableEndpoint& right)
        -> std::optional<ComparableEndpoint> {
        if (left.infinity != 0 || right.infinity != 0) return std::nullopt;
        if (left.radical_coefficient != Rational(0) &&
            right.radical_coefficient != Rational(0) &&
            left.radicand != right.radicand) {
            return std::nullopt;
        }
        const Rational radicand = left.radical_coefficient != Rational(0)
            ? left.radicand : right.radicand;
        return ComparableEndpoint{
            0,
            left.rational + right.rational,
            left.radical_coefficient + right.radical_coefficient,
            radicand};
    };

    auto combine_multiply = [](const ComparableEndpoint& left,
                               const ComparableEndpoint& right)
        -> std::optional<ComparableEndpoint> {
        if (left.infinity != 0 || right.infinity != 0) return std::nullopt;
        const bool left_radical = left.radical_coefficient != Rational(0);
        const bool right_radical = right.radical_coefficient != Rational(0);
        if (left_radical && right_radical && left.radicand != right.radicand) {
            return std::nullopt;
        }
        const Rational radicand = left_radical ? left.radicand : right.radicand;
        const Rational rational = left.rational * right.rational +
            left.radical_coefficient * right.radical_coefficient * radicand;
        const Rational coefficient =
            left.rational * right.radical_coefficient +
            left.radical_coefficient * right.rational;
        return ComparableEndpoint{0, rational, coefficient, radicand};
    };

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        ComparableEndpoint result{0, Rational(0), Rational(0), Rational(0)};
        for (const auto& operand : add->operands()) {
            auto parsed = parse_quadratic_surd(operand);
            if (!parsed) return std::nullopt;
            auto combined = combine_add(result, *parsed);
            if (!combined) return std::nullopt;
            result = std::move(*combined);
        }
        return result;
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        ComparableEndpoint result{0, Rational(1), Rational(0), Rational(0)};
        for (const auto& operand : multiply->operands()) {
            auto parsed = parse_quadratic_surd(operand);
            if (!parsed) return std::nullopt;
            auto combined = combine_multiply(result, *parsed);
            if (!combined) return std::nullopt;
            result = std::move(*combined);
        }
        return result;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = exact_number_value(power->exponent());
        if (!exponent) return std::nullopt;
        if (*exponent == Rational(1, 2)) {
            auto radicand = exact_number_value(power->base());
            if (!radicand || *radicand <= Rational(0)) return std::nullopt;
            return ComparableEndpoint{0, Rational(0), Rational(1), *radicand};
        }
        auto base = parse_quadratic_surd(power->base());
        if (!base || *exponent != Rational(-1)) return std::nullopt;
        const Rational norm = base->rational * base->rational -
            base->radical_coefficient * base->radical_coefficient * base->radicand;
        if (norm == Rational(0)) return std::nullopt;
        return ComparableEndpoint{
            0,
            base->rational / norm,
            (Rational(0) - base->radical_coefficient) / norm,
            base->radicand};
    }

    return std::nullopt;
}

int rational_sign(const Rational& value) {
    if (value == Rational(0)) return 0;
    return value.get_numerator().IsNegative() ? -1 : 1;
}

int quadratic_surd_sign(const Rational& rational,
                        const Rational& coefficient,
                        const Rational& radicand) {
    const int rational_part_sign = rational_sign(rational);
    const int radical_part_sign = rational_sign(coefficient);
    if (radical_part_sign == 0) return rational_part_sign;
    if (rational_part_sign == 0 || rational_part_sign == radical_part_sign) {
        return radical_part_sign;
    }

    const Rational rational_square = rational * rational;
    const Rational radical_square = coefficient * coefficient * radicand;
    if (rational_square == radical_square) return 0;
    if (rational_part_sign > 0) {
        return rational_square > radical_square ? 1 : -1;
    }
    return radical_square > rational_square ? 1 : -1;
}

Result<int> compare_comparable(const ComparableEndpoint& left,
                               const ComparableEndpoint& right,
                               const std::string& operation) {
    if (left.infinity < right.infinity) return Result<int>::success(-1);
    if (left.infinity > right.infinity) return Result<int>::success(1);
    if (left.infinity != 0) return Result<int>::success(0);

    const bool left_radical = left.radical_coefficient != Rational(0);
    const bool right_radical = right.radical_coefficient != Rational(0);
    if (!left_radical && !right_radical) {
        if (left.rational < right.rational) return Result<int>::success(-1);
        if (left.rational > right.rational) return Result<int>::success(1);
        return Result<int>::success(0);
    }
    if (left_radical && right_radical && left.radicand != right.radicand) {
        return Result<int>::failure(
            CasErrc::Inconclusive,
            "exact comparison across distinct quadratic extensions is not proven",
            operation);
    }
    const Rational radicand = left_radical ? left.radicand : right.radicand;
    return Result<int>::success(quadratic_surd_sign(
        left.rational - right.rational,
        left.radical_coefficient - right.radical_coefficient,
        radicand));
}

Result<ComparableEndpoint> comparable_endpoint(
    const Endpoint& endpoint,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<ComparableEndpoint>::failure(step.error());

    if (endpoint.is_neg_infinity && endpoint.is_pos_infinity) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InvalidArgument,
            "an endpoint cannot be both negative and positive infinity",
            operation);
    }
    if (endpoint.is_neg_infinity || endpoint.is_pos_infinity) {
        if (endpoint.value) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::InvalidArgument,
                "infinite endpoints cannot also contain a finite value",
                operation);
        }
        if (!endpoint.is_open) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::InvalidArgument,
                "infinite interval endpoints must be open",
                operation);
        }
        return Result<ComparableEndpoint>::success(
            ComparableEndpoint{endpoint.is_neg_infinity ? -1 : 1,
                               Rational(0), Rational(0), Rational(0)});
    }
    if (!endpoint.value || !lamina::detail::node(endpoint.value)) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InvalidArgument,
            "finite interval endpoint must contain an expression",
            operation);
    }

    auto simplified = endpoint.value->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InternalInvariant,
            "interval endpoint simplification produced a null expression",
            operation);
    }

    if (auto number = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(number->value())) {
            return Result<ComparableEndpoint>::success(ComparableEndpoint{
                0, Rational(std::get<BigInt>(number->value())),
                Rational(0), Rational(0)});
        }
        if (std::holds_alternative<Rational>(number->value())) {
            return Result<ComparableEndpoint>::success(ComparableEndpoint{
                0, std::get<Rational>(number->value()), Rational(0), Rational(0)});
        }
        auto rational = exact_double_rational(
            std::get<lmmc_real_t>(number->value()), context, operation);
        if (!rational) return Result<ComparableEndpoint>::failure(rational.error());
        return Result<ComparableEndpoint>::success(
            ComparableEndpoint{0, std::move(rational.value()),
                               Rational(0), Rational(0)});
    }

    if (auto surd = parse_quadratic_surd(lamina::detail::node(simplified))) {
        return Result<ComparableEndpoint>::success(std::move(*surd));
    }

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(simplified))) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::UnboundSymbol,
            "interval endpoint contains unbound symbol '" + variable->name() + "'",
            operation);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(simplified))) {
        if (function->arguments().size() == 1) {
            auto argument = exact_number_value(function->arguments()[0]);
            if (argument &&
                ((function->type() == FunctionNode::FuncType::Ln &&
                  *argument <= Rational(0)) ||
                 (function->type() == FunctionNode::FuncType::Sqrt &&
                  *argument < Rational(0)))) {
                return Result<ComparableEndpoint>::failure(
                    CasErrc::DomainError,
                    "interval endpoint is outside the real function domain",
                    operation);
            }
        }
    }
    return Result<ComparableEndpoint>::failure(
        CasErrc::Inconclusive,
        "exact ordering of a symbolic interval endpoint is not proven",
        operation);
}

Result<std::vector<CheckedInterval>> checked_interval_views(
    const std::vector<Interval>& intervals,
    ComputationContext& context,
    const std::string& operation) {
    auto normalized = normalize_intervals_checked(intervals, context);
    if (!normalized) {
        return Result<std::vector<CheckedInterval>>::failure(normalized.error());
    }

    std::vector<CheckedInterval> checked;
    auto normalized_intervals = std::move(normalized.value());
    checked.reserve(normalized_intervals.size());
    for (auto& interval : normalized_intervals) {
        auto lower = comparable_endpoint(interval.lower, context, operation);
        if (!lower) return Result<std::vector<CheckedInterval>>::failure(lower.error());
        auto upper = comparable_endpoint(interval.upper, context, operation);
        if (!upper) return Result<std::vector<CheckedInterval>>::failure(upper.error());
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

} // namespace

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
        if (iv.contains(value)) return true;
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

Result<IntervalUnion> IntervalUnion::intersect_checked(
    const IntervalUnion& other,
    ComputationContext& context) const {
    constexpr const char* operation = "interval_union_intersect";
    auto left = checked_interval_views(intervals_, context, operation);
    if (!left) return Result<IntervalUnion>::failure(left.error());
    auto right = checked_interval_views(other.intervals_, context, operation);
    if (!right) return Result<IntervalUnion>::failure(right.error());

    const auto& a_intervals = left.value();
    const auto& b_intervals = right.value();
    std::vector<Interval> result;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a_intervals.size() && j < b_intervals.size()) {
        const auto& a = a_intervals[i];
        const auto& b = b_intervals[j];
        auto lower_order = compare_comparable(a.lower, b.lower, operation);
        if (!lower_order) return Result<IntervalUnion>::failure(lower_order.error());
        auto upper_order = compare_comparable(a.upper, b.upper, operation);
        if (!upper_order) return Result<IntervalUnion>::failure(upper_order.error());
        const bool use_a_lower = lower_order.value() >= 0;
        const bool use_a_upper = upper_order.value() <= 0;

        Interval candidate{
            use_a_lower ? a.interval.lower : b.interval.lower,
            use_a_upper ? a.interval.upper : b.interval.upper
        };
        auto empty = interval_is_empty_checked(candidate, context);
        if (!empty) return Result<IntervalUnion>::failure(empty.error());
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
    if (!normalized) return Result<IntervalUnion>::failure(normalized.error());
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
    if (!step) return Result<IntervalUnion>::failure(step.error());
    auto normalized = normalize_intervals_checked(std::move(all), context);
    if (!normalized) return Result<IntervalUnion>::failure(normalized.error());
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
    if (!checked) return Result<IntervalUnion>::failure(checked.error());
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
    if (!normalized) return Result<IntervalUnion>::failure(normalized.error());
    return Result<IntervalUnion>::success(
        IntervalUnion::from_checked_normalized(std::move(normalized.value())));
}

Result<IntervalUnion> IntervalUnion::complement_checked() const {
    ComputationContext context;
    return complement_checked(context);
}

Result<bool> interval_contains_checked(
    const Interval& interval,
    double value,
    ComputationContext& context) {
    constexpr const char* operation = "interval_contains";
    auto point_step = context.consume_steps(1, operation);
    if (!point_step) return Result<bool>::failure(point_step.error());
    auto point = exact_double_rational(value, context, operation);
    if (!point) return Result<bool>::failure(point.error());
    const ComparableEndpoint point_key{
        0, std::move(point.value()), Rational(0), Rational(0)};

    auto lower = comparable_endpoint(interval.lower, context, operation);
    if (!lower) return Result<bool>::failure(lower.error());
    auto upper = comparable_endpoint(interval.upper, context, operation);
    if (!upper) return Result<bool>::failure(upper.error());

    auto lower_cmp = compare_comparable(point_key, lower.value(), operation);
    if (!lower_cmp) return Result<bool>::failure(lower_cmp.error());
    if (lower_cmp.value() < 0 ||
        (lower_cmp.value() == 0 && interval.lower.is_open)) {
        return Result<bool>::success(false);
    }
    auto upper_cmp = compare_comparable(point_key, upper.value(), operation);
    if (!upper_cmp) return Result<bool>::failure(upper_cmp.error());
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
    if (!lower) return Result<bool>::failure(lower.error());
    auto upper = comparable_endpoint(interval.upper, context, operation);
    if (!upper) return Result<bool>::failure(upper.error());
    auto comparison = compare_comparable(
        lower.value(), upper.value(), operation);
    if (!comparison) return Result<bool>::failure(comparison.error());
    return Result<bool>::success(
        comparison.value() > 0 ||
        (comparison.value() == 0 &&
         (interval.lower.is_open || interval.upper.is_open)));
}

Result<bool> interval_is_empty_checked(const Interval& interval) {
    ComputationContext context;
    return interval_is_empty_checked(interval, context);
}

Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals,
    ComputationContext& context) {
    struct CheckedInterval {
        Interval interval;
        ComparableEndpoint lower;
        ComparableEndpoint upper;
    };

    std::vector<CheckedInterval> checked;
    checked.reserve(intervals.size());
    for (auto& interval : intervals) {
        auto lower = comparable_endpoint(
            interval.lower, context, kCheckedIntervalOperation);
        if (!lower) return Result<std::vector<Interval>>::failure(lower.error());
        auto upper = comparable_endpoint(
            interval.upper, context, kCheckedIntervalOperation);
        if (!upper) return Result<std::vector<Interval>>::failure(upper.error());
        auto comparison = compare_comparable(
            lower.value(), upper.value(), kCheckedIntervalOperation);
        if (!comparison) {
            return Result<std::vector<Interval>>::failure(comparison.error());
        }
        const bool empty = comparison.value() > 0 ||
            (comparison.value() == 0 &&
             (interval.lower.is_open || interval.upper.is_open));
        if (!empty) {
            checked.push_back(CheckedInterval{
                std::move(interval), std::move(lower.value()), std::move(upper.value())});
        }
    }

    for (std::size_t i = 1; i < checked.size(); ++i) {
        std::size_t position = i;
        while (position > 0) {
            auto comparison = compare_comparable(
                checked[position - 1].lower,
                checked[position].lower,
                kCheckedIntervalOperation);
            if (!comparison) {
                return Result<std::vector<Interval>>::failure(comparison.error());
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

    std::vector<CheckedInterval> merged;
    merged.reserve(checked.size());
    for (auto& next : checked) {
        if (merged.empty()) {
            merged.push_back(std::move(next));
            continue;
        }

        CheckedInterval& current = merged.back();
        auto boundary = compare_comparable(
            current.upper, next.lower, kCheckedIntervalOperation);
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

        auto upper_comparison = compare_comparable(
            current.upper, next.upper, kCheckedIntervalOperation);
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

Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals) {
    ComputationContext context;
    return normalize_intervals_checked(std::move(intervals), context);
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

std::string IntervalUnion::to_string() const {

    if (intervals_.empty()) {
        return "\xe2\x88\x85";
    }

    if (intervals_.size() == 1 && intervals_[0].is_entire_line()) {
        return "(-\xe2\x88\x9e, +\xe2\x88\x9e)";
    }

    std::string result;
    for (size_t i = 0; i < intervals_.size(); ++i) {
        if (i > 0) {
            result += " \xe2\x88\xaa ";
        }

        const Interval& iv = intervals_[i];

        if (iv.lower.is_open) {
            result += "(";
        } else {
            result += "[";
        }

        if (iv.lower.is_neg_infinity) {
            result += "-\xe2\x88\x9e";
        } else if (iv.lower.value) {
            result += iv.lower.value->to_string();
        } else {
            result += "0";
        }

        result += ", ";

        if (iv.upper.is_pos_infinity) {
            result += "+\xe2\x88\x9e";
        } else if (iv.upper.value) {
            result += iv.upper.value->to_string();
        } else {
            result += "0";
        }

        if (iv.upper.is_open) {
            result += ")";
        } else {
            result += "]";
        }
    }

    return result;
}

static bool starts_with_at(const std::string& str, size_t pos, const std::string& prefix) {
    if (pos + prefix.size() > str.size()) return false;
    return str.compare(pos, prefix.size(), prefix) == 0;
}

static size_t skip_ws(const std::string& str, size_t pos) {
    while (pos < str.size() && (str[pos] == ' ' || str[pos] == '\t')) {
        ++pos;
    }
    return pos;
}

static std::shared_ptr<SymbolicExpr> parse_numeric_value(const std::string& str, size_t& pos) {
    size_t start = pos;
    bool has_digit = false;
    bool has_dot = false;

    if (pos < str.size() && (str[pos] == '-' || str[pos] == '+')) {
        ++pos;
    }

    while (pos < str.size() && std::isdigit(static_cast<unsigned char>(str[pos]))) {
        has_digit = true;
        ++pos;
    }

    if (pos < str.size() && str[pos] == '.') {
        has_dot = true;
        ++pos;
        while (pos < str.size() && std::isdigit(static_cast<unsigned char>(str[pos]))) {
            has_digit = true;
            ++pos;
        }
    }

    if (!has_digit) {
        pos = start;
        return nullptr;
    }

    std::string num_str = str.substr(start, pos - start);
    if (has_dot) {
        double val = std::stod(num_str);
        return SymbolicExpr::number(val);
    } else {

        long long val = std::stoll(num_str);
        return SymbolicExpr::number(val);
    }
}

static bool parse_endpoint_value(const std::string& str, size_t& pos, Endpoint& ep) {
    pos = skip_ws(str, pos);

    if (starts_with_at(str, pos, "-\xe2\x88\x9e")) {
        ep.is_neg_infinity = true;
        ep.is_pos_infinity = false;
        ep.is_open = true;
        ep.value = nullptr;
        pos += 4;
        return true;
    }

    if (starts_with_at(str, pos, "+\xe2\x88\x9e")) {
        ep.is_pos_infinity = true;
        ep.is_neg_infinity = false;
        ep.is_open = true;
        ep.value = nullptr;
        pos += 4;
        return true;
    }

    auto val = parse_numeric_value(str, pos);
    if (!val) return false;

    ep.value = val;
    ep.is_neg_infinity = false;
    ep.is_pos_infinity = false;

    return true;
}

static bool parse_single_interval(const std::string& str, size_t& pos, Interval& iv) {
    pos = skip_ws(str, pos);
    if (pos >= str.size()) return false;

    char lower_bracket = str[pos];
    if (lower_bracket != '(' && lower_bracket != '[') return false;
    bool lower_open = (lower_bracket == '(');
    ++pos;

    Endpoint lower;
    if (!parse_endpoint_value(str, pos, lower)) return false;

    if (!lower.is_neg_infinity && !lower.is_pos_infinity) {
        lower.is_open = lower_open;
    }

    pos = skip_ws(str, pos);
    if (pos >= str.size() || str[pos] != ',') return false;
    ++pos;

    Endpoint upper;
    if (!parse_endpoint_value(str, pos, upper)) return false;

    pos = skip_ws(str, pos);
    if (pos >= str.size()) return false;
    char upper_bracket = str[pos];
    if (upper_bracket != ')' && upper_bracket != ']') return false;
    bool upper_open = (upper_bracket == ')');
    ++pos;

    if (!upper.is_neg_infinity && !upper.is_pos_infinity) {
        upper.is_open = upper_open;
    }

    iv.lower = lower;
    iv.upper = upper;
    return true;
}

std::optional<IntervalUnion> IntervalUnion::parse(const std::string& str) {
    if (str.empty()) return std::nullopt;

    if (str == "\xe2\x88\x85") {
        return IntervalUnion::empty();
    }

    std::vector<Interval> intervals;
    size_t pos = 0;

    Interval iv;
    if (!parse_single_interval(str, pos, iv)) return std::nullopt;
    intervals.push_back(iv);

    while (pos < str.size()) {
        pos = skip_ws(str, pos);
        if (pos >= str.size()) break;

        if (!starts_with_at(str, pos, "\xe2\x88\xaa")) {
            return std::nullopt;
        }
        pos += 3;

        pos = skip_ws(str, pos);

        Interval next_iv;
        if (!parse_single_interval(str, pos, next_iv)) return std::nullopt;
        intervals.push_back(next_iv);
    }

    return IntervalUnion(std::move(intervals));
}

std::shared_ptr<SymbolicExpr> IntervalUnion::to_expr(const std::string& var) const {

    if (intervals_.empty()) {
        return nullptr;
    }

    auto var_node = lamina::detail::node(SymbolicExpr::variable(var));

    auto interval_to_expr = [&](const Interval& iv) -> std::shared_ptr<const SymbolicNode> {
        std::shared_ptr<const SymbolicNode> lower_cond = nullptr;
        std::shared_ptr<const SymbolicNode> upper_cond = nullptr;

        if (!iv.lower.is_neg_infinity) {
            auto bound = iv.lower.value ? lamina::detail::node(iv.lower.value) : lamina::detail::make_node<NumberNode>(0.0);
            RelationalNode::Op op = iv.lower.is_open ? RelationalNode::Op::GT : RelationalNode::Op::GEQ;
            lower_cond = lamina::detail::make_node<RelationalNode>(var_node, bound, op);
        }

        if (!iv.upper.is_pos_infinity) {
            auto bound = iv.upper.value ? lamina::detail::node(iv.upper.value) : lamina::detail::make_node<NumberNode>(0.0);
            RelationalNode::Op op = iv.upper.is_open ? RelationalNode::Op::LT : RelationalNode::Op::LEQ;
            upper_cond = lamina::detail::make_node<RelationalNode>(var_node, bound, op);
        }

        if (lower_cond && upper_cond) {
            return lamina::detail::make_node<LogicalNode>(lower_cond, upper_cond, LogicalNode::Op::And);
        } else if (lower_cond) {
            return lower_cond;
        } else if (upper_cond) {
            return upper_cond;
        } else {

            return nullptr;
        }
    };

    if (intervals_.size() == 1 && intervals_[0].is_entire_line()) {

        return nullptr;
    }

    if (intervals_.size() == 1) {
        auto node = interval_to_expr(intervals_[0]);
        if (!node) return nullptr;
        return lamina::detail::make_expression_ptr(node);
    }

    auto result = interval_to_expr(intervals_[0]);
    for (size_t i = 1; i < intervals_.size(); ++i) {
        auto next = interval_to_expr(intervals_[i]);
        if (result && next) {
            result = lamina::detail::make_node<LogicalNode>(result, next, LogicalNode::Op::Or);
        } else if (next) {
            result = next;
        }

    }

    if (!result) return nullptr;
    return lamina::detail::make_expression_ptr(result);
}

}
