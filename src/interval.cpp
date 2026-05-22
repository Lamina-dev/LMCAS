#include "../include/interval.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <cctype>

namespace lamina {

// ============================================================================
// Endpoint static factory methods
// ============================================================================

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

// ============================================================================
// Helper: extract numeric value from an Endpoint
// ============================================================================

static double endpoint_numeric_value(const Endpoint& ep) {
    if (ep.is_neg_infinity) return -std::numeric_limits<double>::infinity();
    if (ep.is_pos_infinity) return std::numeric_limits<double>::infinity();
    if (ep.value) return ep.value->to_numeric();
    return 0.0;
}

// ============================================================================
// Interval methods
// ============================================================================

bool Interval::contains(double value) const {
    // Check lower bound
    if (lower.is_neg_infinity) {
        // -∞ always passes lower check
    } else {
        double lo = endpoint_numeric_value(lower);
        if (lower.is_open) {
            if (value <= lo) return false;
        } else {
            if (value < lo) return false;
        }
    }

    // Check upper bound
    if (upper.is_pos_infinity) {
        // +∞ always passes upper check
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
    // Infinity endpoints can never form an empty interval (unless inverted)
    if (lower.is_neg_infinity || upper.is_pos_infinity) return false;
    if (lower.is_pos_infinity) return true;  // (+∞, ...) is empty
    if (upper.is_neg_infinity) return true;  // (..., -∞) is empty

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
    // An interval where lower > upper: (1, 0) is empty
    return Interval{Endpoint::open(SymbolicExpr::number(1)), Endpoint::open(SymbolicExpr::number(0))};
}

Interval Interval::entire_line() {
    return Interval{Endpoint::neg_inf(), Endpoint::pos_inf()};
}

Interval Interval::point(std::shared_ptr<SymbolicExpr> val) {
    return Interval{Endpoint::closed(val), Endpoint::closed(val)};
}

// ============================================================================
// IntervalUnion constructors and factories
// ============================================================================

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

// ============================================================================
// IntervalUnion query methods
// ============================================================================

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

// ============================================================================
// IntervalUnion::normalize()
// ============================================================================

// Helper: compare two endpoints numerically for sorting
// Returns negative if a < b, 0 if equal, positive if a > b
static int compare_endpoints_lower(const Endpoint& a, const Endpoint& b) {
    // -∞ is less than everything
    if (a.is_neg_infinity && b.is_neg_infinity) return 0;
    if (a.is_neg_infinity) return -1;
    if (b.is_neg_infinity) return 1;

    // +∞ is greater than everything
    if (a.is_pos_infinity && b.is_pos_infinity) return 0;
    if (a.is_pos_infinity) return 1;
    if (b.is_pos_infinity) return -1;

    double va = endpoint_numeric_value(a);
    double vb = endpoint_numeric_value(b);

    if (va < vb) return -1;
    if (va > vb) return 1;

    // Same numeric value: closed < open for lower bounds (closed comes first)
    if (a.is_open && !b.is_open) return 1;
    if (!a.is_open && b.is_open) return -1;
    return 0;
}

// Helper: for upper bounds, closed > open (closed extends further)
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

    // Same numeric value: closed > open for upper bounds (closed extends further)
    if (!a.is_open && b.is_open) return 1;
    if (a.is_open && !b.is_open) return -1;
    return 0;
}

// Helper: check if two intervals overlap or are adjacent (can be merged)
static bool can_merge(const Interval& a, const Interval& b) {
    // a's upper must reach b's lower for them to overlap/be adjacent
    // They overlap if a.upper >= b.lower (considering open/closed)
    // They are adjacent if a.upper == b.lower and at least one is closed at that point

    if (a.upper.is_pos_infinity) return true;
    if (b.lower.is_neg_infinity) return true;

    if (a.upper.is_neg_infinity || b.lower.is_pos_infinity) return false;

    double au = endpoint_numeric_value(a.upper);
    double bl = endpoint_numeric_value(b.lower);

    if (au > bl) return true;
    if (au == bl) {
        // Adjacent or overlapping at a single point
        // Can merge if at least one endpoint is closed at this value
        return !a.upper.is_open || !b.lower.is_open;
    }
    return false;
}

// Helper: merge the upper endpoint (take the larger one)
static Endpoint max_upper(const Endpoint& a, const Endpoint& b) {
    int cmp = compare_endpoints_upper(a, b);
    return (cmp >= 0) ? a : b;
}

void IntervalUnion::normalize() {
    // 1. Remove empty intervals
    intervals_.erase(
        std::remove_if(intervals_.begin(), intervals_.end(),
                       [](const Interval& iv) { return iv.is_empty(); }),
        intervals_.end());

    if (intervals_.empty()) return;

    // 2. Sort by lower bound
    std::sort(intervals_.begin(), intervals_.end(),
              [](const Interval& a, const Interval& b) {
                  return compare_endpoints_lower(a.lower, b.lower) < 0;
              });

    // 3. Merge overlapping/adjacent intervals
    std::vector<Interval> merged;
    merged.push_back(intervals_[0]);

    for (size_t i = 1; i < intervals_.size(); ++i) {
        Interval& current = merged.back();
        const Interval& next = intervals_[i];

        if (can_merge(current, next)) {
            // Merge: extend current's upper bound
            current.upper = max_upper(current.upper, next.upper);
        } else {
            merged.push_back(next);
        }
    }

    intervals_ = std::move(merged);
}

// ============================================================================
// IntervalUnion::intersect() — Two-pointer sweep algorithm
// Requirement 5.4: intersection contains only points present in both operands
// ============================================================================

IntervalUnion IntervalUnion::intersect(const IntervalUnion& other) const {
    if (intervals_.empty() || other.intervals_.empty()) {
        return IntervalUnion::empty();
    }

    std::vector<Interval> result;
    size_t i = 0, j = 0;

    while (i < intervals_.size() && j < other.intervals_.size()) {
        const Interval& a = intervals_[i];
        const Interval& b = other.intervals_[j];

        // Compute overlap lower bound: max of the two lower bounds
        Endpoint lo = (compare_endpoints_lower(a.lower, b.lower) >= 0) ? a.lower : b.lower;
        // Compute overlap upper bound: min of the two upper bounds
        Endpoint hi = (compare_endpoints_upper(a.upper, b.upper) <= 0) ? a.upper : b.upper;

        // Check if the resulting interval is non-empty
        Interval candidate{lo, hi};
        if (!candidate.is_empty()) {
            result.push_back(candidate);
        }

        // Advance the pointer with the smaller upper bound
        if (compare_endpoints_upper(a.upper, b.upper) < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    // Result is already sorted and disjoint since inputs are sorted/disjoint
    IntervalUnion res;
    res.intervals_ = std::move(result);
    return res;
}

// ============================================================================
// IntervalUnion::unite() — Concatenate all intervals and normalize
// Requirement 5.5: union contains all points present in either operand
// ============================================================================

IntervalUnion IntervalUnion::unite(const IntervalUnion& other) const {
    std::vector<Interval> all;
    all.reserve(intervals_.size() + other.intervals_.size());
    all.insert(all.end(), intervals_.begin(), intervals_.end());
    all.insert(all.end(), other.intervals_.begin(), other.intervals_.end());
    return IntervalUnion(std::move(all));  // constructor calls normalize()
}

// ============================================================================
// IntervalUnion::complement() — Gap-based complement algorithm
// Requirement 5.6: complement contains all reals not in the original set
// Requirement 5.8: A ∩ complement(A) = ∅, A ∪ complement(A) = ℝ
// ============================================================================

IntervalUnion IntervalUnion::complement() const {
    // Empty set's complement is the entire real line
    if (intervals_.empty()) {
        return IntervalUnion::entire_line();
    }

    // Entire line's complement is empty
    if (is_entire_line()) {
        return IntervalUnion::empty();
    }

    std::vector<Interval> result;

    // Gap before the first interval: (-∞, first.lower) with flipped openness
    const Interval& first = intervals_[0];
    if (!first.lower.is_neg_infinity) {
        Endpoint gap_upper;
        gap_upper.value = first.lower.value;
        gap_upper.is_open = !first.lower.is_open;  // flip open/closed
        gap_upper.is_neg_infinity = false;
        gap_upper.is_pos_infinity = false;

        Interval gap{Endpoint::neg_inf(), gap_upper};
        if (!gap.is_empty()) {
            result.push_back(gap);
        }
    }

    // Gaps between consecutive intervals
    for (size_t i = 0; i + 1 < intervals_.size(); ++i) {
        const Interval& curr = intervals_[i];
        const Interval& next = intervals_[i + 1];

        // Lower of gap: flip curr.upper openness
        Endpoint gap_lower;
        gap_lower.value = curr.upper.value;
        gap_lower.is_open = !curr.upper.is_open;  // flip
        gap_lower.is_neg_infinity = false;
        gap_lower.is_pos_infinity = false;

        // Upper of gap: flip next.lower openness
        Endpoint gap_upper;
        gap_upper.value = next.lower.value;
        gap_upper.is_open = !next.lower.is_open;  // flip
        gap_upper.is_neg_infinity = false;
        gap_upper.is_pos_infinity = false;

        Interval gap{gap_lower, gap_upper};
        if (!gap.is_empty()) {
            result.push_back(gap);
        }
    }

    // Gap after the last interval: (last.upper, +∞) with flipped openness
    const Interval& last = intervals_.back();
    if (!last.upper.is_pos_infinity) {
        Endpoint gap_lower;
        gap_lower.value = last.upper.value;
        gap_lower.is_open = !last.upper.is_open;  // flip
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
    // Empty set
    if (intervals_.empty()) {
        return "\xe2\x88\x85";  // ∅ (U+2205)
    }

    // Entire real line shortcut
    if (intervals_.size() == 1 && intervals_[0].is_entire_line()) {
        return "(-\xe2\x88\x9e, +\xe2\x88\x9e)";  // (-∞, +∞)
    }

    std::string result;
    for (size_t i = 0; i < intervals_.size(); ++i) {
        if (i > 0) {
            result += " \xe2\x88\xaa ";  // " ∪ " (U+222A with spaces)
        }

        const Interval& iv = intervals_[i];

        // Lower bracket
        if (iv.lower.is_open) {
            result += "(";
        } else {
            result += "[";
        }

        // Lower value
        if (iv.lower.is_neg_infinity) {
            result += "-\xe2\x88\x9e";  // -∞
        } else if (iv.lower.value) {
            result += iv.lower.value->to_string();
        } else {
            result += "0";
        }

        result += ", ";

        // Upper value
        if (iv.upper.is_pos_infinity) {
            result += "+\xe2\x88\x9e";  // +∞
        } else if (iv.upper.value) {
            result += iv.upper.value->to_string();
        } else {
            result += "0";
        }

        // Upper bracket
        if (iv.upper.is_open) {
            result += ")";
        } else {
            result += "]";
        }
    }

    return result;
}

// ============================================================================
// IntervalUnion::parse() helpers
// ============================================================================

// Helper: check if string starts with a given prefix at position pos
static bool starts_with_at(const std::string& str, size_t pos, const std::string& prefix) {
    if (pos + prefix.size() > str.size()) return false;
    return str.compare(pos, prefix.size(), prefix) == 0;
}

// Helper: skip whitespace
static size_t skip_ws(const std::string& str, size_t pos) {
    while (pos < str.size() && (str[pos] == ' ' || str[pos] == '\t')) {
        ++pos;
    }
    return pos;
}

// Helper: parse a numeric value (integer or decimal, possibly negative)
// Returns the parsed SymbolicExpr and advances pos past the number.
// Returns nullptr if no valid number found.
static std::shared_ptr<SymbolicExpr> parse_numeric_value(const std::string& str, size_t& pos) {
    size_t start = pos;
    bool has_sign = false;
    bool has_digit = false;
    bool has_dot = false;

    // Optional sign
    if (pos < str.size() && (str[pos] == '-' || str[pos] == '+')) {
        has_sign = true;
        ++pos;
    }

    // Digits before decimal point
    while (pos < str.size() && std::isdigit(static_cast<unsigned char>(str[pos]))) {
        has_digit = true;
        ++pos;
    }

    // Optional decimal point and digits after
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
        // Try as integer
        long long val = std::stoll(num_str);
        return SymbolicExpr::number(val);
    }
}

// Helper: parse an endpoint value (either -∞, +∞, or a numeric value)
// Returns true on success, sets the endpoint fields.
static bool parse_endpoint_value(const std::string& str, size_t& pos, Endpoint& ep) {
    pos = skip_ws(str, pos);

    // Check for -∞ : "-" followed by ∞ (UTF-8: 0xE2 0x88 0x9E)
    if (starts_with_at(str, pos, "-\xe2\x88\x9e")) {
        ep.is_neg_infinity = true;
        ep.is_pos_infinity = false;
        ep.is_open = true;  // infinity is always open
        ep.value = nullptr;
        pos += 4;  // "-" (1 byte) + "∞" (3 bytes)
        return true;
    }

    // Check for +∞ : "+" followed by ∞
    if (starts_with_at(str, pos, "+\xe2\x88\x9e")) {
        ep.is_pos_infinity = true;
        ep.is_neg_infinity = false;
        ep.is_open = true;  // infinity is always open
        ep.value = nullptr;
        pos += 4;  // "+" (1 byte) + "∞" (3 bytes)
        return true;
    }

    // Try to parse a numeric value
    auto val = parse_numeric_value(str, pos);
    if (!val) return false;

    ep.value = val;
    ep.is_neg_infinity = false;
    ep.is_pos_infinity = false;
    // is_open will be set by the bracket parsing
    return true;
}

// Helper: parse a single interval like "[3, +∞)" or "(-∞, -2)"
// Returns true on success.
static bool parse_single_interval(const std::string& str, size_t& pos, Interval& iv) {
    pos = skip_ws(str, pos);
    if (pos >= str.size()) return false;

    // Parse lower bracket
    char lower_bracket = str[pos];
    if (lower_bracket != '(' && lower_bracket != '[') return false;
    bool lower_open = (lower_bracket == '(');
    ++pos;

    // Parse lower value
    Endpoint lower;
    if (!parse_endpoint_value(str, pos, lower)) return false;

    // If it's not infinity, set the open/closed from bracket
    if (!lower.is_neg_infinity && !lower.is_pos_infinity) {
        lower.is_open = lower_open;
    }
    // -∞ is always open regardless of bracket

    // Expect comma
    pos = skip_ws(str, pos);
    if (pos >= str.size() || str[pos] != ',') return false;
    ++pos;

    // Parse upper value
    Endpoint upper;
    if (!parse_endpoint_value(str, pos, upper)) return false;

    // Parse upper bracket
    pos = skip_ws(str, pos);
    if (pos >= str.size()) return false;
    char upper_bracket = str[pos];
    if (upper_bracket != ')' && upper_bracket != ']') return false;
    bool upper_open = (upper_bracket == ')');
    ++pos;

    // If it's not infinity, set the open/closed from bracket
    if (!upper.is_neg_infinity && !upper.is_pos_infinity) {
        upper.is_open = upper_open;
    }
    // +∞ is always open regardless of bracket

    iv.lower = lower;
    iv.upper = upper;
    return true;
}

std::optional<IntervalUnion> IntervalUnion::parse(const std::string& str) {
    if (str.empty()) return std::nullopt;

    // Check for empty set: ∅ (UTF-8: 0xE2 0x88 0x85)
    if (str == "\xe2\x88\x85") {
        return IntervalUnion::empty();
    }

    std::vector<Interval> intervals;
    size_t pos = 0;

    // Parse first interval
    Interval iv;
    if (!parse_single_interval(str, pos, iv)) return std::nullopt;
    intervals.push_back(iv);

    // Parse additional intervals separated by " ∪ " (UTF-8: space + 0xE2 0x88 0xAA + space)
    while (pos < str.size()) {
        pos = skip_ws(str, pos);
        if (pos >= str.size()) break;

        // Check for union separator: " ∪ " — we already skipped leading whitespace
        // The ∪ character is UTF-8: 0xE2 0x88 0xAA
        if (!starts_with_at(str, pos, "\xe2\x88\xaa")) {
            return std::nullopt;  // unexpected character
        }
        pos += 3;  // skip ∪ (3 bytes)

        // Skip trailing space after ∪
        pos = skip_ws(str, pos);

        // Parse next interval
        Interval next_iv;
        if (!parse_single_interval(str, pos, next_iv)) return std::nullopt;
        intervals.push_back(next_iv);
    }

    return IntervalUnion(std::move(intervals));
}

std::shared_ptr<SymbolicExpr> IntervalUnion::to_expr(const std::string& var) const {
    // Empty set: no solution, return nullptr
    if (intervals_.empty()) {
        return nullptr;
    }

    auto var_node = SymbolicExpr::variable(var)->root;

    // Helper: convert a single interval to a SymbolicExpr using relational nodes and And
    auto interval_to_expr = [&](const Interval& iv) -> std::shared_ptr<SymbolicNode> {
        std::shared_ptr<SymbolicNode> lower_cond = nullptr;
        std::shared_ptr<SymbolicNode> upper_cond = nullptr;

        // Lower bound condition: var > a or var >= a
        if (!iv.lower.is_neg_infinity) {
            auto bound = iv.lower.value ? iv.lower.value->root : std::make_shared<NumberNode>(0.0);
            RelationalNode::Op op = iv.lower.is_open ? RelationalNode::Op::GT : RelationalNode::Op::GEQ;
            lower_cond = std::make_shared<RelationalNode>(var_node, bound, op);
        }

        // Upper bound condition: var < b or var <= b
        if (!iv.upper.is_pos_infinity) {
            auto bound = iv.upper.value ? iv.upper.value->root : std::make_shared<NumberNode>(0.0);
            RelationalNode::Op op = iv.upper.is_open ? RelationalNode::Op::LT : RelationalNode::Op::LEQ;
            upper_cond = std::make_shared<RelationalNode>(var_node, bound, op);
        }

        // Combine lower and upper with And
        if (lower_cond && upper_cond) {
            return std::make_shared<LogicalNode>(lower_cond, upper_cond, LogicalNode::Op::And);
        } else if (lower_cond) {
            return lower_cond;
        } else if (upper_cond) {
            return upper_cond;
        } else {
            // Entire line: no constraints (should not happen for individual intervals in a union
            // unless the union is the entire line itself)
            return nullptr;
        }
    };

    // Entire line: (-∞, +∞) has no constraints
    if (intervals_.size() == 1 && intervals_[0].is_entire_line()) {
        // Return a tautology: var = var (or we could return nullptr to indicate "all reals")
        // Per the design, we return nullptr for "no constraint" / entire line
        return nullptr;
    }

    // Single interval
    if (intervals_.size() == 1) {
        auto node = interval_to_expr(intervals_[0]);
        if (!node) return nullptr;
        return std::make_shared<SymbolicExpr>(node);
    }

    // Multiple intervals: connect with Or
    auto result = interval_to_expr(intervals_[0]);
    for (size_t i = 1; i < intervals_.size(); ++i) {
        auto next = interval_to_expr(intervals_[i]);
        if (result && next) {
            result = std::make_shared<LogicalNode>(result, next, LogicalNode::Op::Or);
        } else if (next) {
            result = next;
        }
        // If next is nullptr (entire line interval), the whole thing is trivially true
        // but this shouldn't happen in a normalized IntervalUnion with multiple intervals
    }

    if (!result) return nullptr;
    return std::make_shared<SymbolicExpr>(result);
}

} // namespace lamina
