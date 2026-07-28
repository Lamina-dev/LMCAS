/**
 * @file solve_strategies.hpp
 * @brief 方程求解策略调度器 solve_dispatch 及相关选项结构体。
 */
#pragma once

#include "symbolic.hpp"
#include "computation_context.hpp"
#include "conditional_result.hpp"
#include <vector>
#include <memory>
#include <string>

namespace lamina {

/** @brief 方程求解选项，控制数值求解行为和结果形式。 */
struct SolveOptions {
    bool allow_numeric = false;           ///< 是否允许数值求解
    int max_newton_iterations = 100;      ///< Newton 迭代最大次数
    lmmc_real_t tolerance = 1e-12;        ///< 收敛容差
    int max_roots = -1;                   ///< 最大返回根数，-1 表示不限制
    bool return_rootof = true;            ///< 无闭式解时是否返回 RootOf 表达式
    lmmc_real_t initial_guess = 0.0;      ///< 数值迭代初始猜测值
    bool has_initial_guess = false;       ///< 是否指定了初始猜测值

    lmmc_real_t search_lo = -10.0;        ///< 搜索区间下界（默认 -10）
    lmmc_real_t search_hi = 10.0;         ///< 搜索区间上界（默认 10）
    bool has_search_interval = false;     ///< 是否指定了搜索区间
};

/** @brief 数值根结果，包含根值、残差和迭代次数。 */
struct NumericRoot {
    lmmc_real_t value;      ///< 根的数值
    lmmc_real_t residual;   ///< 残差（将根代入方程后的绝对值）
    int iterations;         ///< 收敛所用迭代次数
};

/** @brief 求解策略枚举，标识所使用的求解方法。 */
enum class SolveStrategy {
    ClosedForm,       ///< 闭式公式求解
    Preprocessing,    ///< 预处理（化简、换元等）
    Transcendental,   ///< 超越方程求解
    Numerical,        ///< 数值迭代求解
    RootOf            ///< 以 RootOf 表达式表示
};

/**
 * @brief 方程求解策略调度器，根据表达式类型选择合适的求解方法。
 * @param expr 待求解的表达式（视为等于零）
 * @param var 求解变量名
 * @param opts 求解选项
 * @return 所有根的符号表达式列表
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_dispatch(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts);

using SolveVectorResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

/**
 * @brief Compatibility checked dispatcher for callers still expecting vectors.
 *
 * This preserves the legacy vector result shape while returning CasError for
 * invalid input and checked numeric fallback failures. Prefer
 * solve_dispatch_checked for new 2.0 code because it can represent Empty,
 * Universal, Conditional, and Inconclusive solution sets.
 */
LAMINA_API SolveVectorResult solve_dispatch_vector_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts = {});

/**
 * @brief Checked legacy-vector dispatcher with a default computation context.
 *
 * Use this as the migration path from solve_dispatch() when the caller needs
 * explicit CasError propagation but does not need custom budgets or
 * cancellation.
 */
LAMINA_API SolveVectorResult solve_dispatch_vector_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

using SolveResult = Result<SolutionSet>;

/**
 * @brief Solve an equation using the verified 2.0 support domain.
 *
 * Exact rational polynomials produce exact RootOf solutions. Approximate
 * candidates are produced only when `allow_numeric` is set and are verified by
 * checked numeric evaluation. Unsupported symbolic domains return an
 * inconclusive SolutionSet rather than an empty set.
 */
LAMINA_API SolveResult solve_dispatch_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts = {});

/**
 * @brief Checked solution-set dispatcher with a default computation context.
 */
LAMINA_API SolveResult solve_dispatch_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

}
