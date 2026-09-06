/**
 * @file series_engine.hpp
 * @brief 级数引擎模块:收敛性判定,幂级数运算,傅里叶级数,洛朗级数,渐近展开,符号求和与乘积.
 */
#pragma once

#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>

namespace LMCAS {

using PowerSeriesResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;


/**
 * @brief 收敛性结果枚举.
 */
enum class ConvergenceResult {
    Convergent,   ///< 收敛
    Divergent,    ///< 发散
    Inconclusive  ///< 无法判定
};

/**
 * @brief 收敛性判定信息.
 */
struct ConvergenceInfo {
    ConvergenceResult result;  ///< 判定结果
    std::string test_used;     ///< 使用的判定方法("ratio", "root", "comparison")
};

using ConvergenceInfoResult = Result<ConvergenceInfo>;

/**
 * @brief 从无限级数通项计算幂级数收敛半径.
 *
 * 当前对常数通项,b^n 和 n^p 给出精确结果;无法证明尾项行为时返回
 * ::CasErrc::Inconclusive,不会从有限样本外推极限.
 *
 * @param[in] general_coefficient 系数通项 a_n
 * @param[in] index_var 系数指标变量名
 * @return 已证明的收敛半径表达式
 */
LMCAS_API ExpressionResult convergence_radius_checked(
    const std::shared_ptr<SymbolicExpr>& general_coefficient,
    const std::string& index_var,
    ComputationContext& context);

LMCAS_API ExpressionResult convergence_radius_checked(
    const std::shared_ptr<SymbolicExpr>& general_coefficient,
    const std::string& index_var);

/**
 * @brief 返回有限系数多项式的收敛半径.
 *
 * coefficients 完整表示有限多项式,未列出的高阶系数均为零,因此非空
 * 有效输入的收敛半径恒为 infinity.若要表示无限级数,应使用通项重载;
 * 有限前缀不足以证明无限尾项的收敛半径.
 */
LMCAS_API ExpressionResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var,
    ComputationContext& context);

LMCAS_API ExpressionResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var);


/**
 * @brief 判定无穷级数的收敛性.
 *
 * 依次尝试比值判别法,根值判别法和比较判别法.
 *
 * @param[in] general_term 通项表达式 a_n
 * @param[in] index_var 求和指标变量名
 * @return 收敛性判定信息
 */
LMCAS_API ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var,
    ComputationContext& context);

LMCAS_API ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var);



/**
 * @brief 幂级数加法:逐项相加.
 *
 * @param[in] a 第一个幂级数系数列表
 * @param[in] b 第二个幂级数系数列表
 * @return 和级数的系数列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> power_series_add(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b);

/**
 * @brief 幂级数乘法:柯西乘积(卷积).
 *
 * @param[in] a 第一个幂级数系数列表
 * @param[in] b 第二个幂级数系数列表
 * @param[in] order 截断阶数
 * @return 乘积级数的系数列表(截断到 order 项)
 */
LMCAS_API PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order,
    ComputationContext& context);

LMCAS_API PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order);


/**
 * @brief 幂级数复合:计算 f(g(x)) 的系数,要求 g(0) = 0.
 *
 * @param[in] f 外层幂级数系数列表
 * @param[in] g 内层幂级数系数列表(g[0] 必须为 0)
 * @param[in] order 截断阶数
 * @return 复合级数的系数列表
 */
LMCAS_API PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order,
    ComputationContext& context);

LMCAS_API PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order);



/**
 * @brief Computes a bounded truncated Fourier expansion.
 *
 * All coefficient integrations share @p context. Resource exhaustion,
 * cancellation, invalid inputs, and nested integration failures are returned
 * without being collapsed to a null expression.
 */
LMCAS_API ExpressionResult fourier_series_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms,
    ComputationContext& context);

LMCAS_API ExpressionResult fourier_series_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms);


/**
 * @brief 奇点类型枚举.
 */
enum class SingularityType {
    Removable,  ///< 可去奇点
    Pole,       ///< 极点
    Essential   ///< 本性奇点
};

/**
 * @brief 洛朗级数展开结果.
 */
struct LaurentResult {
    std::shared_ptr<SymbolicExpr> series;  ///< 洛朗级数表达式
    SingularityType singularity;           ///< 奇点类型
    int pole_order;                        ///< 极点阶数(非极点时为 0)
    std::shared_ptr<SymbolicExpr> residue; ///< 留数((z-z_0)-¹ 的系数)
};

using LaurentSeriesResult = Result<LaurentResult>;

/**
 * @brief 计算函数在指定点的洛朗级数展开.
 *
 * 展开包含从 (z-z_0)^(-order_neg) 到 (z-z_0)^(order_pos) 的项.
 * 通过计算 (z-z_0)^k * f(z) 在 z_0 处的极限和导数来确定负幂次系数.
 * 根据洛朗系数判定奇点类型:
 *   - 所有负幂次系数为零 -> 可去奇点
 *   - 有限个非零负幂次系数 -> m 阶极点
 *   - 无穷多非零负幂次系数 -> 本性奇点
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] center 展开中心 z_0
 * @param[in] order_neg 负幂次最大阶数
 * @param[in] order_pos 正幂次最大阶数
 * @return 洛朗级数表达式
 */
LMCAS_API ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context);

LMCAS_API ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos);


/**
 * @brief 计算函数在指定点的洛朗级数展开(含详细信息).
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] center 展开中心 z_0
 * @param[in] order_neg 负幂次最大阶数
 * @param[in] order_pos 正幂次最大阶数
 * @return 洛朗级数展开结果(含奇点分类和留数)
 */
LMCAS_API LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context);

LMCAS_API LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos);



/**
 * @brief 计算函数在 x->infinity 时的渐近展开.
 *
 * 将函数展开为递减幂次的和:sum cₖ * x^alphaₖ,其中 alpha_1 > alpha_2 > ... > alpha_n.
 * 对有理函数使用多项式长除法;对含 exp,ln,幂函数的复合表达式使用已知渐近恒等式.
 *
 * @param[in] f 待展开的函数表达式
 * @param[in] var 变量名
 * @param[in] order 展开项数
 * @return 渐近展开表达式
 */
LMCAS_API std::shared_ptr<SymbolicExpr> asymptotic_expand(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var, int order);


/**
 * @brief 计算符号求和 sum_{k=lower}^{upper} body(k) 的闭合形式.
 *
 * 对多项式通项使用 Faulhaber 公式;对有理函数通项尝试部分分式分解和望远镜求和;
 * 对几何/指数通项使用等比级数公式;其余通项映射为未求值 SummationNode.
 *
 * @param[in] body 通项表达式
 * @param[in] index 求和指标变量名
 * @param[in] lower 下界
 * @param[in] upper 上界
 * @return 求和结果表达式
 */
LMCAS_API std::shared_ptr<SymbolicExpr> symbolic_sum(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper);


/**
 * @brief 计算符号乘积 prod_{k=lower}^{upper} body(k) 的闭合形式.
 *
 * 对望远镜乘积和阶乘/Pochhammer 符号使用已知公式;
 * 其余通项映射为未求值 ProductNode.
 *
 * @param[in] body 通项表达式
 * @param[in] index 乘积指标变量名
 * @param[in] lower 下界
 * @param[in] upper 上界
 * @return 乘积结果表达式
 */
LMCAS_API std::shared_ptr<SymbolicExpr> symbolic_product(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper);


/**
 * @brief 计算数列的上极限.
 */
LMCAS_API ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context);

LMCAS_API ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n);

/**
 * @brief 计算数列的下极限.
 */
LMCAS_API ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context);

LMCAS_API ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n);

} // namespace LMCAS
