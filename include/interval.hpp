/**
 * @file interval.hpp
 * @brief 区间与区间并集 IntervalUnion，用于不等式解集表示。
 */
#pragma once
#include "computation_context.hpp"
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>
#include <optional>

namespace lamina {

class InequalitySolver;

/** @brief 不等式类型 */
enum class InequalityType {
    GreaterThan,  ///< 严格大于
    GreaterEqual, ///< 大于等于
    LessThan,     ///< 严格小于
    LessEqual     ///< 小于等于
};

/** @brief 区间端点，支持有限值与正负无穷 */
struct LAMINA_API Endpoint {
    std::shared_ptr<SymbolicExpr> value{}; ///< 端点值（无穷时为空）
    bool is_open = false;          ///< 是否为开端点
    bool is_neg_infinity = false;  ///< 是否为负无穷
    bool is_pos_infinity = false;  ///< 是否为正无穷

    /**
     * @brief 构造负无穷端点
     * @return 负无穷端点
     */
    static Endpoint neg_inf();

    /**
     * @brief 构造正无穷端点
     * @return 正无穷端点
     */
    static Endpoint pos_inf();

    /**
     * @brief 构造闭端点
     * @param val 端点值
     * @return 闭端点
     */
    static Endpoint closed(std::shared_ptr<SymbolicExpr> val);

    /**
     * @brief 构造开端点
     * @param val 端点值
     * @return 开端点
     */
    static Endpoint open(std::shared_ptr<SymbolicExpr> val);
};

/** @brief 单个区间 [a, b] 或 (a, b) 等 */
struct LAMINA_API Interval {
    Endpoint lower; ///< 下端点
    Endpoint upper; ///< 上端点

    /**
     * @brief 判断数值是否在区间内
     * @param value 待判断的数值
     * @return 包含返回 true
     */
    bool contains(double value) const;

    /**
     * @brief 判断区间是否为空
     * @return 空区间返回 true
     */
    bool is_empty() const;

    /**
     * @brief 判断区间是否为整条实数轴
     * @return 是整条实数轴返回 true
     */
    bool is_entire_line() const;

    /**
     * @brief 构造空区间
     * @return 空区间
     */
    static Interval empty();

    /**
     * @brief 构造整条实数轴 (-∞, +∞)
     * @return 全实数轴区间
     */
    static Interval entire_line();

    /**
     * @brief 构造单点区间 [val, val]
     * @param val 点值
     * @return 单点区间
     */
    static Interval point(std::shared_ptr<SymbolicExpr> val);
};

/** @brief 区间并集，表示不等式的解集 */
class LAMINA_API IntervalUnion {
public:
    IntervalUnion();

    /**
     * @brief 从区间列表构造并集
     * @param intervals 区间列表
     */
    explicit IntervalUnion(std::vector<Interval> intervals);

    /**
     * @brief 从单个区间构造并集
     * @param iv 区间
     * @return 包含单个区间的并集
     */
    static IntervalUnion from_single(const Interval& iv);

    /** @brief Validate and construct an interval union with an explicit context. */
    static Result<IntervalUnion> from_intervals_checked(
        std::vector<Interval> intervals,
        ComputationContext& context);

    /** @brief Validate and construct an interval union with a default context. */
    static Result<IntervalUnion> from_intervals_checked(
        std::vector<Interval> intervals);

    /**
     * @brief 构造空集
     * @return 空的区间并集
     */
    static IntervalUnion empty();

    /**
     * @brief 构造全实数轴
     * @return 表示 (-∞, +∞) 的并集
     */
    static IntervalUnion entire_line();

    /**
     * @brief 与另一个并集求交
     * @param other 另一个区间并集
     * @return 交集结果
     */
    IntervalUnion intersect(const IntervalUnion& other) const;

    /**
     * @brief Checked intersection with explicit endpoint-validation errors.
     * @param other 另一个区间并集
     * @param context 计算上下文和资源预算
     * @return 成功时返回交集；端点不可比较或资源耗尽时返回错误
     */
    Result<IntervalUnion> intersect_checked(
        const IntervalUnion& other,
        ComputationContext& context) const;

    Result<IntervalUnion> intersect_checked(const IntervalUnion& other) const;

    /**
     * @brief 与另一个并集求并
     * @param other 另一个区间并集
     * @return 并集结果
     */
    IntervalUnion unite(const IntervalUnion& other) const;

    /**
     * @brief Checked union with exact endpoint ordering.
     * @param other 另一个区间并集
     * @param context 计算上下文和资源预算
     * @return 成功时返回并集；端点不可比较或资源耗尽时返回错误
     */
    Result<IntervalUnion> unite_checked(
        const IntervalUnion& other,
        ComputationContext& context) const;

    Result<IntervalUnion> unite_checked(const IntervalUnion& other) const;

    /**
     * @brief 求补集
     * @return 补集结果
     */
    IntervalUnion complement() const;

    /**
     * @brief Checked complement with exact endpoint ordering.
     * @param context 计算上下文和资源预算
     * @return 成功时返回补集；端点不可比较或资源耗尽时返回错误
     */
    Result<IntervalUnion> complement_checked(ComputationContext& context) const;

    Result<IntervalUnion> complement_checked() const;

    /**
     * @brief 判断数值是否在并集内
     * @param value 待判断的数值
     * @return 包含返回 true
     */
    bool contains(double value) const;

    /**
     * @brief 判断并集是否为空
     * @return 空集返回 true
     */
    bool is_empty() const;

    /**
     * @brief 判断并集是否为全实数轴
     * @return 是全实数轴返回 true
     */
    bool is_entire_line() const;

    /**
     * @brief 获取内部区间列表
     * @return 区间列表的常引用
     */
    const std::vector<Interval>& intervals() const;

    /**
     * @brief 转换为字符串表示
     * @return 区间并集的字符串形式
     */
    std::string to_string() const;

    /**
     * @brief 从字符串解析区间并集
     * @param str 字符串表示
     * @return 解析成功返回并集，失败返回 nullopt
     */
    static std::optional<IntervalUnion> parse(const std::string& str);

    /**
     * @brief 转换为符号表达式
     * @param var 变量名
     * @return 表示该区间并集的符号表达式
     */
    std::shared_ptr<SymbolicExpr> to_expr(const std::string& var) const;

private:
    friend class InequalitySolver;

    std::vector<Interval> intervals_;
    void normalize();
    static IntervalUnion from_checked_normalized(std::vector<Interval> intervals);
};

/** @brief Checked numeric membership test for an interval. */
LAMINA_API Result<bool> interval_contains_checked(
    const Interval& interval,
    double value,
    ComputationContext& context);

LAMINA_API Result<bool> interval_contains_checked(
    const Interval& interval,
    double value);

/** @brief Checked emptiness decision for an interval. */
LAMINA_API Result<bool> interval_is_empty_checked(
    const Interval& interval,
    ComputationContext& context);

LAMINA_API Result<bool> interval_is_empty_checked(const Interval& interval);

/**
 * @brief Validate, sort, remove empty intervals, and merge proven overlaps.
 *
 * Exact endpoints remain exact. Approximate endpoints are ordered by their
 * exact IEEE binary values. Any incomparable or invalid endpoint returns an
 * error before a normalized list is produced.
 */
LAMINA_API Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals,
    ComputationContext& context);

LAMINA_API Result<std::vector<Interval>> normalize_intervals_checked(
    std::vector<Interval> intervals);

}
