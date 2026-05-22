#pragma once
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>
#include <optional>

namespace lamina {

// 不等式类型枚举
enum class InequalityType {
    GreaterThan,   // >
    GreaterEqual,  // >=
    LessThan,      // <
    LessEqual      // <=
};

// 区间端点
struct LAMINA_API Endpoint {
    std::shared_ptr<SymbolicExpr> value;  // 数值或符号表达式
    bool is_open;                          // true = 开, false = 闭
    bool is_neg_infinity;                  // -∞
    bool is_pos_infinity;                  // +∞

    static Endpoint neg_inf();
    static Endpoint pos_inf();
    static Endpoint closed(std::shared_ptr<SymbolicExpr> val);
    static Endpoint open(std::shared_ptr<SymbolicExpr> val);
};

// 单个区间 [a, b], (a, b), [a, b), (a, b]
struct LAMINA_API Interval {
    Endpoint lower;
    Endpoint upper;

    bool contains(double value) const;
    bool is_empty() const;
    bool is_entire_line() const;

    static Interval empty();
    static Interval entire_line();
    static Interval point(std::shared_ptr<SymbolicExpr> val);  // [a, a]
};

// 区间并集: 不相交区间的有序集合
class LAMINA_API IntervalUnion {
public:
    IntervalUnion();
    explicit IntervalUnion(std::vector<Interval> intervals);
    static IntervalUnion from_single(const Interval& iv);
    static IntervalUnion empty();
    static IntervalUnion entire_line();

    // 集合运算
    IntervalUnion intersect(const IntervalUnion& other) const;
    IntervalUnion unite(const IntervalUnion& other) const;
    IntervalUnion complement() const;

    // 查询
    bool contains(double value) const;
    bool is_empty() const;
    bool is_entire_line() const;
    const std::vector<Interval>& intervals() const;

    // 序列化
    std::string to_string() const;
    static std::optional<IntervalUnion> parse(const std::string& str);

    // 转为符号表达式 (And/Or 逻辑连接)
    std::shared_ptr<SymbolicExpr> to_expr(const std::string& var) const;

private:
    std::vector<Interval> intervals_;  // 按 lower bound 排序, 不相交
    void normalize();  // 合并重叠/相邻区间, 排序
};

} // namespace lamina
