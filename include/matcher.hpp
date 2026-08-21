/**
 * @file matcher.hpp
 * @brief 模式匹配引擎 Matcher 与重写规则系统 RewriteEngine。
 */
#pragma once
#include "computation_context.hpp"
#include "symbolic.hpp"
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <functional>

namespace lamina {

// Forward declaration for optional assumption context integration
class AssumptionContext;

/** @brief 匹配结果映射，通配符名称 → 绑定的表达式 */
using MatchMap = std::unordered_map<std::string, SymbolicExpr>;

/** @brief 符号表达式模式匹配器 */
class LAMINA_API Matcher {
public:

    /**
     * @brief 将模式与目标表达式进行匹配
     * @param pattern 模式表达式
     * @param target 目标表达式
     * @param wildcards 模式中的通配符名称集合
     * @param results 匹配成功时输出绑定结果
     * @return 匹配成功返回 true
     */
    static bool match(const SymbolicExpr& pattern, const SymbolicExpr& target,
                      const std::unordered_set<std::string>& wildcards,
                      MatchMap& results);

    /**
     * @brief 将绑定结果代入模板表达式
     * @param template_expr 模板表达式
     * @param bindings 通配符绑定映射
     * @param use_rest 是否使用剩余项替换
     * @return 替换后的表达式
     */
    static SymbolicExpr replace(const SymbolicExpr& template_expr, const MatchMap& bindings, bool use_rest = true);
};

/** @brief 重写规则，定义模式到替换表达式的变换 */
struct LAMINA_API Rule {
    SymbolicExpr pattern;                          ///< 匹配模式
    SymbolicExpr replacement;                      ///< 替换模板
    std::unordered_set<std::string> wildcards;     ///< 通配符集合

    std::function<bool(const MatchMap&)> condition; ///< 附加匹配条件（可选）

    /// Assumption-aware condition: receives bindings and the assumption context.
    /// If set, this is preferred over `condition` when a context is available.
    std::function<bool(const MatchMap&, const AssumptionContext*)> assumption_condition;

    Rule(SymbolicExpr p, SymbolicExpr r, std::unordered_set<std::string> w,
         std::function<bool(const MatchMap&)> c = nullptr)
        : pattern(p), replacement(r), wildcards(w), condition(c) {}

    Rule(SymbolicExpr p, SymbolicExpr r, std::unordered_set<std::string> w,
         std::function<bool(const MatchMap&, const AssumptionContext*)> ac)
        : pattern(p), replacement(r), wildcards(w), assumption_condition(ac) {}
};

/** @brief 基于规则的表达式重写引擎 */
class LAMINA_API RewriteEngine {
    std::vector<Rule> rules;

public:
    /**
     * @brief 添加重写规则
     * @param rule 规则对象
     */
    void add_rule(const Rule& rule);

    /**
     * @brief 反复应用规则直到不动点或达到最大迭代次数
     * @param expr 输入表达式
     * @param max_iterations 最大迭代次数
     * @return 重写后的表达式
     */
    SymbolicExpr apply(const SymbolicExpr& expr, int max_iterations = 100) const;
    Result<SymbolicExpr> apply_checked(
        const SymbolicExpr& expr,
        ComputationContext& context,
        int max_iterations = 100) const;
    /**
     * @brief 对表达式应用一步重写
     * @param expr 输入表达式
     * @return 重写后的表达式，若无规则匹配则返回原表达式
     */
    SymbolicExpr apply_step(const SymbolicExpr& expr) const;
    Result<SymbolicExpr> apply_step_checked(
        const SymbolicExpr& expr,
        ComputationContext& context) const;
    /**
     * @brief 获取当前所有规则
     * @return 规则列表的常引用
     */
    const std::vector<Rule>& get_rules() const { return rules; }
};

/**
 * @brief 创建通配符表达式
 * @param name 通配符名称
 * @return 表示通配符的符号表达式
 */
LAMINA_API SymbolicExpr wildcard(const std::string& name);

}
