/**
 * @file series_engine.hpp
 * @brief 级数引擎模块：收敛性判定、幂级数运算、傅里叶级数、洛朗级数、渐近展开、符号求和与乘积。
 */
#pragma once

#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>

namespace lamina {

using SeriesExprResult = Result<std::shared_ptr<SymbolicExpr>>;
using PowerSeriesResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;


/**
 * @brief 收敛性结果枚举。
 */
enum class ConvergenceResult {
    Convergent,   ///< 收敛
    Divergent,    ///< 发散
    Inconclusive  ///< 无法判定
};

/**
 * @brief 收敛性判定信息。
 */
struct ConvergenceInfo {
    ConvergenceResult result;  ///< 判定结果
    std::string test_used;     ///< 使用的判定方法（"ratio", "root", "comparison"）
};

using ConvergenceInfoResult = Result<ConvergenceInfo>;

/**
 * @brief 计算幂级数的收敛半径。
 *
 * 先尝试比值判别法 R = lim|aₙ/aₙ₊₁|，若不确定则使用根值判别法 R = 1/lim sup |aₙ|^(1/n)。
 *
 * @param[in] coefficients 幂级数系数列表 {a₀, a₁, a₂, ...}
 * @param[in] var 级数变量名
 * @return 收敛半径表达式（可能为数值、∞ 或符号表达式）
 */
LAMINA_API SeriesExprResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var,
    ComputationContext& context);

LAMINA_API SeriesExprResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var);

LAMINA_API std::shared_ptr<SymbolicExpr> convergence_radius(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients, const std::string& var);

/**
 * @brief 判定无穷级数的收敛性。
 *
 * 依次尝试比值判别法、根值判别法和比较判别法。
 *
 * @param[in] general_term 通项表达式 aₙ
 * @param[in] index_var 求和指标变量名
 * @return 收敛性判定信息
 */
LAMINA_API ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var,
    ComputationContext& context);

LAMINA_API ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var);

LAMINA_API ConvergenceInfo convergence_test(
    const std::shared_ptr<SymbolicExpr>& general_term, const std::string& index_var);


/**
 * @brief 幂级数加法：逐项相加。
 *
 * @param[in] a 第一个幂级数系数列表
 * @param[in] b 第二个幂级数系数列表
 * @return 和级数的系数列表
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> power_series_add(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b);

/**
 * @brief 幂级数乘法：柯西乘积（卷积）。
 *
 * @param[in] a 第一个幂级数系数列表
 * @param[in] b 第二个幂级数系数列表
 * @param[in] order 截断阶数
 * @return 乘积级数的系数列表（截断到 order 项）
 */
LAMINA_API PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order,
    ComputationContext& context);

LAMINA_API PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order);

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> power_series_multiply(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order);

/**
 * @brief 幂级数复合：计算 f(g(x)) 的系数，要求 g(0) = 0。
 *
 * @param[in] f 外层幂级数系数列表
 * @param[in] g 内层幂级数系数列表（g[0] 必须为 0）
 * @param[in] order 截断阶数
 * @return 复合级数的系数列表
 */
LAMINA_API PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order,
    ComputationContext& context);

LAMINA_API PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order);

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> power_series_compose(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order);


/**
 * @brief 计算函数的傅里叶级数展开。
 *
 * 计算截断傅里叶级数：a₀/2 + ∑(aₖcos(2πkx/T) + bₖsin(2πkx/T))，k=1..n_terms。
 * 系数公式：aₖ = (2/T)∫f(x)cos(2πkx/T)dx，bₖ = (2/T)∫f(x)sin(2πkx/T)dx。
 * 当函数为偶函数时，所有 bₖ = 0（余弦级数）；
 * 当函数为奇函数时，所有 aₖ = 0（正弦级数）。
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 自变量名
 * @param[in] period 周期 T
 * @param[in] n_terms 展开项数
 * @return 截断傅里叶级数表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> fourier_series(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms);


/**
 * @brief 奇点类型枚举。
 */
enum class SingularityType {
    Removable,  ///< 可去奇点
    Pole,       ///< 极点
    Essential   ///< 本性奇点
};

/**
 * @brief 洛朗级数展开结果。
 */
struct LaurentResult {
    std::shared_ptr<SymbolicExpr> series;  ///< 洛朗级数表达式
    SingularityType singularity;           ///< 奇点类型
    int pole_order;                        ///< 极点阶数（非极点时为 0）
    std::shared_ptr<SymbolicExpr> residue; ///< 留数（(z-z₀)⁻¹ 的系数）
};

using LaurentSeriesResult = Result<LaurentResult>;

/**
 * @brief 计算函数在指定点的洛朗级数展开。
 *
 * 展开包含从 (z-z₀)^(-order_neg) 到 (z-z₀)^(order_pos) 的项。
 * 通过计算 (z-z₀)^k · f(z) 在 z₀ 处的极限和导数来确定负幂次系数。
 * 根据洛朗系数判定奇点类型：
 *   - 所有负幂次系数为零 → 可去奇点
 *   - 有限个非零负幂次系数 → m 阶极点
 *   - 无穷多非零负幂次系数 → 本性奇点
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] center 展开中心 z₀
 * @param[in] order_neg 负幂次最大阶数
 * @param[in] order_pos 正幂次最大阶数
 * @return 洛朗级数表达式
 */
LAMINA_API SeriesExprResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context);

LAMINA_API SeriesExprResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos);

LAMINA_API std::shared_ptr<SymbolicExpr> laurent_series(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos);

/**
 * @brief 计算函数在指定点的洛朗级数展开（含详细信息）。
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] center 展开中心 z₀
 * @param[in] order_neg 负幂次最大阶数
 * @param[in] order_pos 正幂次最大阶数
 * @return 洛朗级数展开结果（含奇点分类和留数）
 */
LAMINA_API LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context);

LAMINA_API LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos);

LAMINA_API LaurentResult laurent_series_full(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos);


/**
 * @brief 计算函数在 x→∞ 时的渐近展开。
 *
 * 将函数展开为递减幂次的和：∑ cₖ · x^αₖ，其中 α₁ > α₂ > ... > αₙ。
 * 对有理函数使用多项式长除法；对含 exp、ln、幂函数的复合表达式使用已知渐近恒等式。
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] order 展开项数
 * @return 渐近展开表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> asymptotic_expand(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var, int order);


/**
 * @brief 计算符号求和 ∑_{k=lower}^{upper} body(k) 的闭合形式。
 *
 * 对多项式通项使用 Faulhaber 公式；对有理函数通项尝试部分分式分解和望远镜求和；
 * 对几何/指数通项使用等比级数公式。无法求得闭合形式时返回未求值的 SummationNode。
 *
 * @param[in] body 通项表达式
 * @param[in] index 求和指标变量名
 * @param[in] lower 下界
 * @param[in] upper 上界
 * @return 求和结果表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> symbolic_sum(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper);


/**
 * @brief 计算符号乘积 ∏_{k=lower}^{upper} body(k) 的闭合形式。
 *
 * 对望远镜乘积和阶乘/Pochhammer 符号使用已知公式。
 * 无法求得闭合形式时返回未求值的 ProductNode_Op。
 *
 * @param[in] body 通项表达式
 * @param[in] index 乘积指标变量名
 * @param[in] lower 下界
 * @param[in] upper 上界
 * @return 乘积结果表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> symbolic_product(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper);


/**
 * @brief 计算数列的上极限 lim sup aₙ。
 * @param[in] a_n 通项表达式
 * @param[in] n 指标变量名
 * @return 上极限表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> lim_sup(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n);

/**
 * @brief 计算数列的下极限 lim inf aₙ。
 * @param[in] a_n 通项表达式
 * @param[in] n 指标变量名
 * @return 下极限表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> lim_inf(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n);

} // namespace lamina
