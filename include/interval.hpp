#pragma once
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>
#include <optional>

namespace lamina {

enum class InequalityType {
    GreaterThan,
    GreaterEqual,
    LessThan,
    LessEqual
};

struct LAMINA_API Endpoint {
    std::shared_ptr<SymbolicExpr> value;
    bool is_open;
    bool is_neg_infinity;
    bool is_pos_infinity;

    static Endpoint neg_inf();
    static Endpoint pos_inf();
    static Endpoint closed(std::shared_ptr<SymbolicExpr> val);
    static Endpoint open(std::shared_ptr<SymbolicExpr> val);
};

struct LAMINA_API Interval {
    Endpoint lower;
    Endpoint upper;

    bool contains(double value) const;
    bool is_empty() const;
    bool is_entire_line() const;

    static Interval empty();
    static Interval entire_line();
    static Interval point(std::shared_ptr<SymbolicExpr> val);
};

class LAMINA_API IntervalUnion {
public:
    IntervalUnion();
    explicit IntervalUnion(std::vector<Interval> intervals);
    static IntervalUnion from_single(const Interval& iv);
    static IntervalUnion empty();
    static IntervalUnion entire_line();

    IntervalUnion intersect(const IntervalUnion& other) const;
    IntervalUnion unite(const IntervalUnion& other) const;
    IntervalUnion complement() const;

    bool contains(double value) const;
    bool is_empty() const;
    bool is_entire_line() const;
    const std::vector<Interval>& intervals() const;

    std::string to_string() const;
    static std::optional<IntervalUnion> parse(const std::string& str);

    std::shared_ptr<SymbolicExpr> to_expr(const std::string& var) const;

private:
    std::vector<Interval> intervals_;
    void normalize();
};

}
