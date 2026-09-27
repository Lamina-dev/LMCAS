/**
 * @file solve_strategies.hpp
 * @brief 方程求解策略调度器 solve_dispatch 及相关选项结构体.
 */
#pragma once

#include "symbolic.hpp"
#include "computation_context.hpp"
#include "conditional_result.hpp"
#include <vector>
#include <memory>
#include <string>

namespace LMCAS {

/** @brief 方程求解选项,控制数值求解行为和结果形式. */
struct SolveOptions {
    bool allow_numeric = false;           ///< 是否允许数值求解
    int max_newton_iterations = 100;      ///< Newton 迭代最大次数
    lmmc_real_t tolerance = 1e-12;        ///< 收敛容差
    int max_roots = -1;                   /**< 数值候选数量预算；精确解集和整数解族保持完整。 */
    bool return_rootof = true;            /**< 为 false 时，任何不可表示的根使整个求解返回 Inconclusive。 */
    lmmc_real_t initial_guess = 0.0;      ///< 数值迭代初始猜测值
    bool has_initial_guess = false;       ///< 是否指定了初始猜测值

    lmmc_real_t search_lo = -10.0;        ///< 搜索区间下界(默认 -10)
    lmmc_real_t search_hi = 10.0;         ///< 搜索区间上界(默认 10)
    bool has_search_interval = false;     ///< 是否指定了搜索区间
};

/** @brief 数值根结果,包含根值,残差和迭代次数. */
struct NumericRoot {
    lmmc_real_t value;      ///< 根的数值
    lmmc_real_t residual;   ///< 残差(将根代入方程后的绝对值)
    int iterations;         ///< 收敛所用迭代次数
};

/** @brief 求解策略枚举,标识所使用的求解方法. */
enum class SolveStrategy {
    ClosedForm,       ///< 闭式公式求解
    Preprocessing,    ///< 预处理(化简,换元等)
    Transcendental,   ///< 超越方程求解
    Numerical,        ///< 数值迭代求解
    RootOf            ///< 以 RootOf 表达式表示
};




/**
 * @brief 在已验证支持域内求方程的完整解集。
 * 精确有理多项式保留全部复根及重数，以精确根式或 RootOf 表示，保持系数精确。
 * 实超越反演保留原操作数定义域、整数参数和未决条件；系数退化未知时，
 * 返回完整的 ConditionalSolutions，并以条件限定通用根。
 * 数值搜索可辅助求解；候选残差全部通过仍不足以证明完备性。
 * 完备性未证或完整定义域无法表示时返回 CasErrc::Inconclusive。
 * 成功表示在上下文假设及所有返回条件下，得到操作定义域内的完整解集。
 * UniversalSolutions 仅在该域内普遍成立，原式极点始终排除。
 * 候选残差恒等式仅证明公共定义域上的正确性；可定义性、根的互异性及覆盖性须另证。
 * 未决逐点条件须保留为条件或 Inconclusive，形式非恒等式不足以消除这些条件。
 * 资源耗尽与取消仍作为外层错误返回。
 */
LMCAS_API SolveResult solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts = {});

/**
 * @brief Solves an equation using a default computation context.
 */
LMCAS_API SolveResult solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});


using FiniteSolveResult =
    Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

/**
 * @brief 投影完整、无条件的有限解集，并保留重数。
 * 空集投影为空向量；无限解集或条件解集返回 Inconclusive。
 */
LMCAS_API FiniteSolveResult solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    ComputationContext& context,
    const SolveOptions& options = {});

LMCAS_API FiniteSolveResult solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const SolveOptions& options = {});
}
