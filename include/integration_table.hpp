/** @file integration_table.hpp */ #pragma once
#include "symbolic.hpp"
#include "matcher.hpp"
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace lamina {

/** @brief 积分表条目，描述一条积分规则（模式 → 结果） */
struct LAMINA_API IntegrationEntry {

    std::string name;                    ///< 规则名称

    SymbolicExpr pattern;                ///< 被积表达式的匹配模式

    SymbolicExpr result;                 ///< 积分结果模板

    std::unordered_set<std::string> wildcards; ///< 模式中的通配符集合

    std::function<bool(const MatchMap&, const std::string& var)> condition; ///< 附加匹配条件

    int priority = 100;                  ///< 优先级，数值越小越优先

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

} // namespace lamina
