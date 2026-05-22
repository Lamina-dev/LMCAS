/**
 * @file integration.hpp
 * @brief 符号积分引擎：策略模式，支持查表、幂律、换元、部分分式、分部积分。
 */
#pragma once
#include "symbolic.hpp"
#include "matcher.hpp"
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <unordered_map>

namespace lamina {

/** @brief 积分表条目，描述一条积分规则（模式 → 结果） */
struct LAMINA_API IntegrationEntry {

    std::string name;                    ///< 规则名称

    SymbolicExpr pattern;                ///< 被积表达式的匹配模式

    SymbolicExpr result;                 ///< 积分结果模板

    std::unordered_set<std::string> wildcards; ///< 模式中的通配符集合

    std::function<bool(const MatchMap&, const std::string& var)> condition; ///< 附加匹配条件

    int priority = 100;                  ///< 优先级，数值越小越优先

    IntegrationEntry() = default;
    IntegrationEntry(std::string name, SymbolicExpr pat, SymbolicExpr res,
                     std::unordered_set<std::string> wc,
                     std::function<bool(const MatchMap&, const std::string& var)> cond = nullptr,
                     int prio = 100)
        : name(std::move(name)), pattern(std::move(pat)), result(std::move(res)),
          wildcards(std::move(wc)), condition(std::move(cond)), priority(prio) {}
};

/** @brief 积分查找表，按类别管理积分规则 */
class LAMINA_API IntegrationTable {
public:
    /** @brief 积分规则类别 */
    enum class Category {
        Polynomial,    ///< 多项式
        Exponential,   ///< 指数函数
        Logarithmic,   ///< 对数函数
        Trigonometric, ///< 三角函数
        InverseTrig,   ///< 反三角函数
        Hyperbolic,    ///< 双曲函数
        Algebraic,     ///< 代数函数
        Special,       ///< 特殊函数
        UserDefined    ///< 用户自定义
    };

    IntegrationTable();

    /**
     * @brief 向指定类别添加积分规则
     * @param cat 规则类别
     * @param entry 积分条目
     */
    void add_entry(Category cat, const IntegrationEntry& entry);

    /**
     * @brief 清空指定类别的所有规则
     * @param cat 规则类别
     */
    void clear_category(Category cat);

    /**
     * @brief 获取指定类别的所有规则
     * @param cat 规则类别
     * @return 该类别下的积分条目列表
     */
    const std::vector<IntegrationEntry>& get_entries(Category cat) const;

    /**
     * @brief 获取所有规则，按优先级排序
     * @return 排序后的积分条目指针列表
     */
    std::vector<const IntegrationEntry*> get_all_sorted() const;

    /** @brief 加载默认积分规则表 */
    void load_defaults();

private:
    std::unordered_map<int, std::vector<IntegrationEntry>> entries_;
    static const std::vector<IntegrationEntry> empty_entries_;
};

class Integrator;

/** @brief 积分策略基类，定义策略接口 */
class LAMINA_API IntegrationStrategy {
public:
    virtual ~IntegrationStrategy() = default;

    /**
     * @brief 尝试对表达式进行积分
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param ctx 积分器上下文
     * @param depth 当前递归深度
     * @return 积分结果，失败返回 nullptr
     */
    virtual std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        int depth = 0) = 0;

    /**
     * @brief 获取策略名称
     * @return 策略名称字符串
     */
    virtual std::string name() const = 0;
};

/** @brief 查表积分策略 */
class LAMINA_API TableLookupStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "TableLookup"; }
};

/** @brief 幂律积分策略 */
class LAMINA_API PowerRuleStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PowerRule"; }
};

/** @brief 换元积分策略 */
class LAMINA_API SubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "Substitution"; }
};

/** @brief 部分分式积分策略 */
class LAMINA_API PartialFractionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PartialFraction"; }
};

/** @brief 分部积分策略 */
class LAMINA_API IBPStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "IntegrationByParts"; }
};

/** @brief 符号积分器，协调各策略完成积分运算 */
class LAMINA_API Integrator {
public:
    Integrator();

    /**
     * @brief 计算不定积分
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @return 积分结果表达式
     */
    SymbolicExpr integrate(const SymbolicExpr& expr, const std::string& var_name);

    /**
     * @brief 计算定积分
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @param lower 积分下限
     * @param upper 积分上限
     * @return 定积分结果表达式
     */
    SymbolicExpr integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                               const SymbolicExpr& lower, const SymbolicExpr& upper);

    /**
     * @brief 添加积分策略
     * @param strategy 策略对象
     * @param position 插入位置，-1 表示追加到末尾
     */
    void add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position = -1);

    /** @brief 获取积分表（可修改） */
    IntegrationTable& table() { return table_; }
    /** @brief 获取积分表（只读） */
    const IntegrationTable& table() const { return table_; }

    /**
     * @brief 递归积分入口
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param depth 当前递归深度
     * @return 积分结果，失败返回 nullptr
     */
    std::shared_ptr<SymbolicExpr> integrate_recursive(
        const SymbolicExpr& expr, const std::string& var, int depth = 0);

    /**
     * @brief 判断表达式是否依赖指定变量
     * @param expr 表达式
     * @param var 变量名
     * @return 依赖返回 true
     */
    static bool depends_on(const SymbolicExpr& expr, const std::string& var);

    /** @brief 获取最大递归深度 */
    int max_depth() const { return max_depth_; }
    /**
     * @brief 设置最大递归深度
     * @param d 深度值
     */
    void set_max_depth(int d) { max_depth_ = d; }

private:
    IntegrationTable table_;
    std::vector<std::unique_ptr<IntegrationStrategy>> strategies_;

    struct CycleState {
        std::vector<SymbolicExpr> history;
    };
    CycleState cycle_state_;

    int max_depth_ = 8;

    std::shared_ptr<SymbolicExpr> apply_linearity(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> dispatch_strategies(
        const SymbolicExpr& expr, const std::string& var, int depth);

    static std::shared_ptr<SymbolicExpr> make_integral_node(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr& expr, const std::string& var);
    void resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx);
};

}
