#pragma once
#include "interval.hpp"
#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace lamina {

// 符号图表中的一个条目
struct SignChartEntry {
    Interval interval;
    int sign;  // +1, -1, 或 0
};

// 含参数不等式的分段解: 根据参数条件给出不同的区间解集
struct LAMINA_API PiecewiseIntervalResult {
    struct Case {
        std::shared_ptr<SymbolicExpr> condition;  // 参数条件 (如 a > 0)
        IntervalUnion solution;                    // 该条件下的解集
    };
    std::vector<Case> cases;

    // 是否为空 (无任何分段)
    bool is_empty() const { return cases.empty(); }

    // 是否为单一解 (非分段)
    bool is_single() const { return cases.size() == 1; }

    // 获取单一解 (当 is_single() 为 true 时)
    IntervalUnion single_solution() const {
        if (cases.empty()) return IntervalUnion::empty();
        return cases[0].solution;
    }
};

// 不等式求解器
class LAMINA_API InequalitySolver {
public:
    // 求解单个多项式不等式: p(x) ⊳ 0
    static IntervalUnion solve_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable);

    // 求解不等式组: {p₁(x) ⊳₁ 0, p₂(x) ⊳₂ 0, ...}
    static IntervalUnion solve_inequalities(
        const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                     InequalityType>>& inequalities,
        const std::string& variable);

    // 求解有理不等式: p(x)/q(x) ⊳ 0
    static IntervalUnion solve_rational_inequality(
        const std::shared_ptr<SymbolicExpr>& numerator,
        const std::shared_ptr<SymbolicExpr>& denominator,
        InequalityType type,
        const std::string& variable);

    // 求解含参数的多项式不等式: p(x; params) ⊳ 0
    // 当根依赖参数时，以符号表达式表示区间边界
    // 当首项系数符号依赖参数时，返回分段解
    // 当参数值导致次数下降时，包含退化情况
    static PiecewiseIntervalResult solve_parametric_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable,
        const std::vector<std::string>& parameters);

private:
    // 构造符号图表
    static std::vector<SignChartEntry> build_sign_chart(
        const std::shared_ptr<SymbolicExpr>& poly,
        const std::string& variable,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    // 根据符号图表和不等式类型选取满足条件的区间
    static IntervalUnion select_intervals(
        const std::vector<SignChartEntry>& chart,
        InequalityType type,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    // 构造含参数的符号区间解集 (根为符号表达式)
    // leading_sign: +1 或 -1 表示首项系数的符号
    static IntervalUnion build_parametric_solution(
        const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
        const std::vector<int>& multiplicities,
        int leading_sign,
        InequalityType type);
};

} // namespace lamina
