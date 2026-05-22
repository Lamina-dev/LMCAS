#include "../include/interval.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <cctype>

namespace lamina {

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

static double endpoint_numeric_value(const Endpoint& ep) {
    if (ep.is_neg_infinity) return -std::numeric_limits<double>::infinity();
    if (ep.is_pos_infinity) return std::numeric_limits<double>::infinity();
    if (ep.value) return ep.value->to_numeric();
    return 0.0;
}

bool Interval::contains(double value) const {

    if (lower.is_neg_infinity) {

    } else {
        double lo = endpoint_numeric_value(lower);
        if (lower.is_open) {
            if (value <= lo) return false;
        } else {
            if (value < lo) return false;
        }
    }

    if (upper.is_pos_infinity) {

    } else {
        double hi = endpoint_numeric_value(upper);
        if (upper.is_open) {
            if (value >= hi) return false;
        } else {
            if (value > hi) return false;
        }
    }

    return true;
}

bool Interval::is_empty() const {

    if (lower.is_neg_infinity || upper.is_pos_infinity) return false;
    if (lower.is_pos_infinity) return true;
    if (upper.is_neg_infinity) return true;

    double lo = endpoint_numeric_value(lower);
    double hi = endpoint_numeric_value(upper);

    if (lo > hi) return true;
    if (lo == hi && (lower.is_open || upper.is_open)) return true;

    return false;
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

static int compare_endpoints_lower(const Endpoint& a, const Endpoint& b) {

    if (a.is_neg_infinity && b.is_neg_infinity) return 0;
    if (a.is_neg_infinity) return -1;
    if (b.is_neg_infinity) return 1;

    if (a.is_pos_infinity && b.is_pos_infinity) return 0;
    if (a.is_pos_infinity) return 1;
    if (b.is_pos_infinity) return -1;

    double va = endpoint_numeric_value(a);
    double vb = endpoint_numeric_value(b);

    if (va < vb) return -1;
    if (va > vb) return 1;

    if (a.is_open && !b.is_open) return 1;
    if (!a.is_open && b.is_open) return -1;
    return 0;
}

static int compare_endpoints_upper(const Endpoint& a, const Endpoint& b) {
    if (a.is_pos_infinity && b.is_pos_infinity) return 0;
    if (a.is_pos_infinity) return 1;
    if (b.is_pos_infinity) return -1;

    if (a.is_neg_infinity && b.is_neg_infinity) return 0;
    if (a.is_neg_infinity) return -1;
    if (b.is_neg_infinity) return 1;

    double va = endpoint_numeric_value(a);
    double vb = endpoint_numeric_value(b);

    if (va < vb) return -1;
    if (va > vb) return 1;

    if (!a.is_open && b.is_open) return 1;
    if (a.is_open && !b.is_open) return -1;
    return 0;
}

static bool can_merge(const Interval& a, const Interval& b) {

    if (a.upper.is_pos_infinity) return true;
    if (b.lower.is_neg_infinity) return true;

    if (a.upper.is_neg_infinity || b.lower.is_pos_infinity) return false;

    double au = endpoint_numeric_value(a.upper);
    double bl = endpoint_numeric_value(b.lower);

    if (au > bl) return true;
    if (au == bl) {

        return !a.upper.is_open || !b.lower.is_open;
    }
    return false;
}

static Endpoint max_upper(const Endpoint& a, const Endpoint& b) {
    int cmp = compare_endpoints_upper(a, b);
    return (cmp >= 0) ? a : b;
}

void IntervalUnion::normalize() {

    intervals_.erase(
        std::remove_if(intervals_.begin(), intervals_.end(),
                       [](const Interval& iv) { return iv.is_empty(); }),
        intervals_.end());

    if (intervals_.empty()) return;

    std::sort(intervals_.begin(), intervals_.end(),
              [](const Interval& a, const Interval& b) {
                  return compare_endpoints_lower(a.lower, b.lower) < 0;
              });

    std::vector<Interval> merged;
    merged.push_back(intervals_[0]);

    for (size_t i = 1; i < intervals_.size(); ++i) {
        Interval& current = merged.back();
        const Interval& next = intervals_[i];

        if (can_merge(current, next)) {

            current.upper = max_upper(current.upper, next.upper);
        } else {
            merged.push_back(next);
        }
    }

    intervals_ = std::move(merged);
}

IntervalUnion IntervalUnion::intersect(const IntervalUnion& other) const {
    if (intervals_.empty() || other.intervals_.empty()) {
        return IntervalUnion::empty();
    }

    std::vector<Interval> result;
    size_t i = 0, j = 0;

    while (i < intervals_.size() && j < other.intervals_.size()) {
        const Interval& a = intervals_[i];
        const Interval& b = other.intervals_[j];

        Endpoint lo = (compare_endpoints_lower(a.lower, b.lower) >= 0) ? a.lower : b.lower;

        Endpoint hi = (compare_endpoints_upper(a.upper, b.upper) <= 0) ? a.upper : b.upper;

        Interval candidate{lo, hi};
        if (!candidate.is_empty()) {
            result.push_back(candidate);
        }

        if (compare_endpoints_upper(a.upper, b.upper) < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    IntervalUnion res;
    res.intervals_ = std::move(result);
    return res;
}

IntervalUnion IntervalUnion::unite(const IntervalUnion& other) const {
    std::vector<Interval> all;
    all.reserve(intervals_.size() + other.intervals_.size());
    all.insert(all.end(), intervals_.begin(), intervals_.end());
    all.insert(all.end(), other.intervals_.begin(), other.intervals_.end());
    return IntervalUnion(std::move(all));
}

IntervalUnion IntervalUnion::complement() const {

    if (intervals_.empty()) {
        return IntervalUnion::entire_line();
    }

    if (is_entire_line()) {
        return IntervalUnion::empty();
    }

    std::vector<Interval> result;

    const Interval& first = intervals_[0];
    if (!first.lower.is_neg_infinity) {
        Endpoint gap_upper;
        gap_upper.value = first.lower.value;
        gap_upper.is_open = !first.lower.is_open;
        gap_upper.is_neg_infinity = false;
        gap_upper.is_pos_infinity = false;

        Interval gap{Endpoint::neg_inf(), gap_upper};
        if (!gap.is_empty()) {
            result.push_back(gap);
        }
    }

    for (size_t i = 0; i + 1 < intervals_.size(); ++i) {
        const Interval& curr = intervals_[i];
        const Interval& next = intervals_[i + 1];

        Endpoint gap_lower;
        gap_lower.value = curr.upper.value;
        gap_lower.is_open = !curr.upper.is_open;
        gap_lower.is_neg_infinity = false;
        gap_lower.is_pos_infinity = false;

        Endpoint gap_upper;
        gap_upper.value = next.lower.value;
        gap_upper.is_open = !next.lower.is_open;
        gap_upper.is_neg_infinity = false;
        gap_upper.is_pos_infinity = false;

        Interval gap{gap_lower, gap_upper};
        if (!gap.is_empty()) {
            result.push_back(gap);
        }
    }

    const Interval& last = intervals_.back();
    if (!last.upper.is_pos_infinity) {
        Endpoint gap_lower;
        gap_lower.value = last.upper.value;
        gap_lower.is_open = !last.upper.is_open;
        gap_lower.is_neg_infinity = false;
        gap_lower.is_pos_infinity = false;

        Interval gap{gap_lower, Endpoint::pos_inf()};
        if (!gap.is_empty()) {
            result.push_back(gap);
        }
    }

    IntervalUnion res;
    res.intervals_ = std::move(result);
    return res;
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
    bool has_sign = false;
    bool has_digit = false;
    bool has_dot = false;

    if (pos < str.size() && (str[pos] == '-' || str[pos] == '+')) {
        has_sign = true;
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

    auto var_node = SymbolicExpr::variable(var)->root;

    auto interval_to_expr = [&](const Interval& iv) -> std::shared_ptr<SymbolicNode> {
        std::shared_ptr<SymbolicNode> lower_cond = nullptr;
        std::shared_ptr<SymbolicNode> upper_cond = nullptr;

        if (!iv.lower.is_neg_infinity) {
            auto bound = iv.lower.value ? iv.lower.value->root : std::make_shared<NumberNode>(0.0);
            RelationalNode::Op op = iv.lower.is_open ? RelationalNode::Op::GT : RelationalNode::Op::GEQ;
            lower_cond = std::make_shared<RelationalNode>(var_node, bound, op);
        }

        if (!iv.upper.is_pos_infinity) {
            auto bound = iv.upper.value ? iv.upper.value->root : std::make_shared<NumberNode>(0.0);
            RelationalNode::Op op = iv.upper.is_open ? RelationalNode::Op::LT : RelationalNode::Op::LEQ;
            upper_cond = std::make_shared<RelationalNode>(var_node, bound, op);
        }

        if (lower_cond && upper_cond) {
            return std::make_shared<LogicalNode>(lower_cond, upper_cond, LogicalNode::Op::And);
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
        return std::make_shared<SymbolicExpr>(node);
    }

    auto result = interval_to_expr(intervals_[0]);
    for (size_t i = 1; i < intervals_.size(); ++i) {
        auto next = interval_to_expr(intervals_[i]);
        if (result && next) {
            result = std::make_shared<LogicalNode>(result, next, LogicalNode::Op::Or);
        } else if (next) {
            result = next;
        }

    }

    if (!result) return nullptr;
    return std::make_shared<SymbolicExpr>(result);
}

}
