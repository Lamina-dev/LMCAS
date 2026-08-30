/**
 * @file inequality_solver.hpp
 * @brief 不等式求解器 InequalitySolver,支持多项式,有理式,参数化不等式.
 */
#pragma once
#include "computation_context.hpp"
#include "interval.hpp"
#include "polynomial.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace lamina {

/** @brief 符号图条目,记录区间及其上的符号 */
struct SignChartEntry {
    Interval interval; ///< 区间
    int sign;          ///< 该区间上表达式的符号(+1, 0, -1)
};

/** @brief 参数化不等式的分段解集 */
struct LAMINA_API PiecewiseIntervalResult {
    /** @brief 单个分段:条件 + 对应解集 */
    struct Case {
        std::shared_ptr<SymbolicExpr> condition; ///< 参数满足的条件
        IntervalUnion solution;                  ///< 该条件下的解集
    };
    std::vector<Case> cases; ///< 所有分段

    /**
     * @brief 判断结果是否为空
     * @return 无分段时返回 true
     */
    bool is_empty() const { return cases.empty(); }

    /**
     * @brief 判断是否只有一个分段
     * @return 恰好一个分段时返回 true
     */
    bool is_single() const { return cases.size() == 1; }

    /**
     * @brief 获取唯一分段的解集
     *
     * is_single() 为真时返回该分段;多分段解依赖参数取值,
     * 因而映射为空集以完整保留分段语义.
     * @return 单一分段的解集;多分段结果映射为空集
     */
    IntervalUnion single_solution() const {
        if (!is_single()) return IntervalUnion::empty();
        return cases[0].solution;
    }
};

/** @brief 不等式求解器,提供多项式,有理式及参数化不等式的求解 */
class LAMINA_API InequalitySolver {
public:

    /**
     * @brief Solve an inequality in the checked exact support domain.
     *
     * The current complete support domain is an exact rational polynomial of
     * degree at most two in @p variable. Approximate coefficients, parameters,
     * and higher-degree expressions return CasErrc::Inconclusive.
     */
    static Result<IntervalUnion> solve_inequality_checked(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable,
        ComputationContext& context);

    static Result<IntervalUnion> solve_inequality_checked(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable);

    /**
     * @brief Solve a conjunction of inequalities through the checked solver.
     *
     * Every component must lie in the checked scalar support domain and every
     * required endpoint comparison must be exact. An empty conjunction denotes
     * the entire real line.
     */
    static Result<IntervalUnion> solve_inequalities_checked(
        const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                     InequalityType>>& inequalities,
        const std::string& variable,
        ComputationContext& context);

    static Result<IntervalUnion> solve_inequalities_checked(
        const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                     InequalityType>>& inequalities,
        const std::string& variable);

    /**
     * @brief 求解单个不等式
     * @param expr 不等式左端表达式(右端为 0)
     * @param type 不等式类型
     * @param variable 求解变量名
     * @return 解集
     */
    static IntervalUnion solve_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable);

    /**
     * @brief 求解不等式组(取交集)
     * @param inequalities 不等式列表,每项为 (表达式, 不等式类型)
     * @param variable 求解变量名
     * @return 解集
     */
    static IntervalUnion solve_inequalities(
        const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                     InequalityType>>& inequalities,
        const std::string& variable);

    /**
     * @brief 求解有理不等式 numerator/denominator >/< 0
     * @param numerator 分子表达式
     * @param denominator 分母表达式
     * @param type 不等式类型
     * @param variable 求解变量名
     * @return 解集
     */
    static IntervalUnion solve_rational_inequality(
        const std::shared_ptr<SymbolicExpr>& numerator,
        const std::shared_ptr<SymbolicExpr>& denominator,
        InequalityType type,
        const std::string& variable);

    /**
     * @brief 求解含参数的不等式,返回分段结果
     * @param expr 不等式左端表达式
     * @param type 不等式类型
     * @param variable 求解变量名
     * @param parameters 参数名列表
     * @return 分段解集
     */
    static PiecewiseIntervalResult solve_parametric_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable,
        const std::vector<std::string>& parameters);

private:

    static Result<IntervalUnion> solve_exact_quadratic_inequality(
        const Polynomial<Rational>& polynomial,
        InequalityType type,
        ComputationContext& context);

    static std::vector<SignChartEntry> build_sign_chart(
        const std::shared_ptr<SymbolicExpr>& poly,
        const std::string& variable,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    static IntervalUnion select_intervals(
        const std::vector<SignChartEntry>& chart,
        InequalityType type,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    static IntervalUnion build_parametric_solution(
        const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
        const std::vector<int>& multiplicities,
        int leading_sign,
        InequalityType type);
};

}
