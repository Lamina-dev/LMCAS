/**
 * @file transcendental_factor.hpp
 * @brief 混合超越方程不可约因式分解器公共 API.
 *
 * 将 Berlekamp/Zassenhaus 风格的因式分解算法从纯多项式域推广到含超越函数的
 * 表达式空间.核心流程:换元 -> 多项式构造 -> 模分解 -> Hensel 提升 -> 因子组合 -> 逆换元.
 */
#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "result.hpp"
#include "symbolic.hpp"

#include <memory>
#include <string>
#include <vector>
#include <cstdint>

namespace LMCAS {

/// 换元映射条目:超越子表达式 -> 代数不定元名
struct TransSubstitution {
    std::shared_ptr<SymbolicExpr> trans_expr;  ///< 原始超越子表达式 (e.g., sin(x))
    std::string indeterminate;                  ///< 代数不定元名 (e.g., "u0")
};

/// 换元结果
struct TransSubstitutionResult {
    std::vector<TransSubstitution> mappings;    ///< 所有换元映射
    std::shared_ptr<SymbolicExpr> poly_expr;   ///< 换元后的多项式表达式
    std::vector<std::shared_ptr<SymbolicExpr>> constraints; ///< 代数约束 (e.g., u0^2+u1^2-1=0)
};


/**
 * @brief 分解含 sin、cos、exp、ln 等超越函数的表达式。
 *
 * 不支持的转换保留原表达式；资源限制、取消和内部错误通过 Result 传播。
 *
 * @param[in] expr 待分解的符号表达式
 * @param[in] var  目标变量名
 * @param[in,out] context 共享计算预算和取消状态
 * @return 不可约因子的列表(乘积等于原表达式,可能含常数因子)
 */
LMCAS_API Result<std::vector<std::shared_ptr<SymbolicExpr>>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context);

LMCAS_API Result<std::vector<std::shared_ptr<SymbolicExpr>>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

/**
 * @brief 检测表达式中的超越函数换元模式.
 *
 * 遍历表达式 AST,识别所有依赖目标变量的超越子表达式,为每个分配代数不定元,
 * 并记录不定元之间的代数约束(如三角恒等式 u_sin^2 + u_cos^2 = 1).
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 换元结果,包含映射表,换元后多项式表达式及约束列表
 */
LMCAS_API TransSubstitutionResult detect_trans_substitutions(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);


/// 多项式构造结果
struct TfPolyBuildResult {
    bool success;                              ///< 转换是否成功
    Polynomial<Rational> poly;                 ///< 主变量多项式(单变量或主变量策略)
    std::string main_variable;                 ///< 选定的主变量名
    std::vector<std::string> param_variables;  ///< 参数变量列表(非主变量的不定元)

    TfPolyBuildResult() : success(false), poly("x") {}
};

/**
 * @brief 将换元后的表达式构造为有理系数多项式。
 *
 * 选择次数最高的不定元为主变量，经 symbolic_to_poly 转为 Polynomial<Rational>。
 * 非主变量系数须可表示为有理数，否则转换失败，保持系数完整并以错误区别于零。
 * success=false 仅表示策略不匹配；其他转换错误由外层 Result 原样传播。
 *
 * @param[in] poly_expr       换元后的符号表达式
 * @param[in] indeterminates  不定元名称列表(如 {"u0", "u1"})
 * @param[in] original_var    原始目标变量名(如 "x")
 * @return 多项式构造结果
 *
 * @internal
 */
LMCAS_API Result<TfPolyBuildResult> tf_build_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& indeterminates,
    const std::string& original_var);

/**
 * @brief 对逆换元后的因子列表执行化简与常数乘子提取.
 *
 * 对每个因子调用 simplify() 规范化,提取数值前导系数,
 * 将所有常数乘子合并为单一数值因子.
 *
 * @param[in,out] factors 因子列表
 * @return 化简后的因子列表(可能含首位常数因子)
 *
 * @internal
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors(
    std::vector<std::shared_ptr<SymbolicExpr>>& factors);

} // namespace LMCAS
