/**
 * @file solve_mixed_transcendental.hpp
 * @brief 混合超越方程求解：符号预处理 + 自适应根隔离 + Newton-Raphson 精化。
 *
 * 本模块处理同时包含超越函数（sin、cos、tan、exp、ln）和代数项的方程，
 * 这类方程无法通过换元化为纯多项式或纯超越形式。求解流程：
 * 1. 通过 factor_transcendental 尝试符号分解；
 * 2. 对不可约因子执行自适应根隔离（符号变化检测 + 单调性确认）；
 * 3. 使用 Newton-Raphson（带区间约束）精化每个隔离根，失败时回退到二分法。
 */
#pragma once

#include "symbolic.hpp"
#include "solve_strategies.hpp"
#include "transcendental_factor.hpp"

#include <memory>
#include <vector>
#include <optional>
#include <string>

namespace lamina {

/// 搜索区间
struct SearchInterval {
    lmmc_real_t lo;  ///< 下界
    lmmc_real_t hi;  ///< 上界
};

/// 隔离子区间，包含至多一个根
struct IsolatedInterval {
    lmmc_real_t lo;   ///< 子区间下界
    lmmc_real_t hi;   ///< 子区间上界
    bool confirmed;   ///< 若导数单调性确认恰含一个根则为 true
};

/**
 * @brief 混合超越方程求解主入口：符号预处理 + 数值根隔离 + 精化。
 *
 * 对输入表达式先调用 factor_transcendental 尝试分解为多个因子：
 * - 可约为 ≤4 次多项式的因子委托给多项式求解器；
 * - 可通过 Lambert W 或直接反演求解的因子委托给 solve_transcendental；
 * - 不可约因子进入数值路径（根隔离 + 精化）。
 * 最终对所有因子的根集取并集，去重排序后返回。
 *
 * @param[in] expr 待求解表达式（视为 = 0）
 * @param[in] var  求解变量名
 * @param[in] opts 求解选项（容差、最大迭代次数、搜索区间等）
 * @return 数值根列表，以 NumberNode 表达式形式返回，按升序排列；
 *         无根时返回空向量
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_mixed_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts);

/**
 * @brief 判断表达式是否包含依赖指定变量的超越函数。
 *
 * 遍历 AST，检测 sin/cos/tan/exp/ln 节点的参数是否依赖 var。
 * 若超越函数的参数不含 var（即为常数），则不计入。
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 包含依赖 var 的超越函数时返回 true，否则返回 false
 */
LAMINA_API bool contains_transcendental_of_var(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

/**
 * @brief 判断换元结果是否为关于换元变量的多项式。
 *
 * 检查 detect_trans_substitutions 的输出：若映射为空或换元后表达式
 * 仍非多项式形式，则返回 false。
 *
 * @param[in] sub_result 换元分析结果
 * @return 换元后为多项式时返回 true，否则返回 false
 */
bool is_polynomial_after_substitution(
    const TransSubstitutionResult& sub_result);

/**
 * @brief 确定搜索区间：用户指定 > 周期扩展 > 默认 [-10, 10]。
 *
 * 优先级：
 * 1. 若 opts.has_search_interval 为 true，使用 [search_lo, search_hi]；
 * 2. 若表达式含周期函数（sin/cos/tan）且参数为 k*x+c 形式，
 *    扩展默认区间覆盖至少两个完整周期，上限 [-100, 100]；
 * 3. 否则使用默认区间 [-10, 10]。
 *
 * @param[in] expr 表达式
 * @param[in] var  变量名
 * @param[in] opts 求解选项
 * @return 有效搜索区间；区间无效（lo >= hi 或宽度 <= tolerance）时返回 nullopt
 */
LAMINA_API std::optional<SearchInterval> determine_search_interval(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts);

/**
 * @brief 自适应根隔离：将搜索区间划分为含至多一个根的子区间。
 *
 * 通过递归二分检测符号变化，定位根所在子区间。当导数可用时，
 * 利用导数符号确认子区间内的单调性（即恰含一个根）。
 *
 * 特殊处理：
 * - 子区间宽度 < 1e-6 时停止细分；
 * - 采样点产生 NaN/Inf 时跳过并缩小区间重新采样；
 * - 超过 max_roots 限制时截断返回。
 *
 * @param[in] expr       表达式
 * @param[in] derivative 表达式的导数（可为 nullptr，跳过单调性确认）
 * @param[in] var        变量名
 * @param[in] interval   搜索区间
 * @param[in] opts       求解选项（max_roots、tolerance）
 * @return 隔离子区间列表
 */
LAMINA_API std::vector<IsolatedInterval> isolate_roots(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const SearchInterval& interval,
    const SolveOptions& opts);

/**
 * @brief 去重并排序数值根。
 *
 * 对候选根进行两两比较，差值小于 10 * tolerance 的视为重复，
 * 保留残差较小者。去重后按升序排列，并应用 max_roots 限制。
 *
 * @param[in,out] roots     候选根列表（会被排序修改）
 * @param[in]     tolerance 去重容差（实际阈值为 10 * tolerance）
 * @param[in]     max_roots 最大返回根数（-1 表示不限制）
 * @return 去重、排序后的根值列表
 */
LAMINA_API std::vector<lmmc_real_t> deduplicate_roots(
    std::vector<NumericRoot>& roots,
    lmmc_real_t tolerance,
    int max_roots);

/**
 * @brief 对隔离区间执行根精化，返回高精度数值根。
 *
 * 使用带区间约束的 Newton-Raphson 迭代精化根，导数不可用时回退到纯二分法。
 * 收敛条件为 |f(x)| < tolerance，达到最大迭代次数后返回残差最小的结果
 * （仅当残差 < tolerance），否则返回 nullopt。
 *
 * @param[in] expr       待求解表达式
 * @param[in] derivative 表达式的导数（可为 nullptr，此时使用纯二分法）
 * @param[in] var        求解变量名
 * @param[in] interval   隔离子区间
 * @param[in] opts       求解选项（tolerance, max_newton_iterations）
 * @return 精化后的数值根；未收敛或残差超限时返回 nullopt
 */
LAMINA_API std::optional<NumericRoot> refine_root(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const IsolatedInterval& interval,
    const SolveOptions& opts);

} // namespace lamina
