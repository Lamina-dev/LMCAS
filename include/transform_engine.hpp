/**
 * @file transform_engine.hpp
 * @brief 积分变换引擎:Laplace 变换,逆 Laplace 变换,Fourier 变换,Z 变换.
 *
 * 提供基于内置规则匹配,线性性质和位移定理的符号积分变换计算;
 * 当前规则集之外的表达式由 checked API 返回 Inconclusive.
 */
#pragma once

#include "computation_context.hpp"
#include "conditional_result.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>


namespace LMCAS {


/**
 * @brief checked 变换仅返回已求值并认证的表达式。
 * 条件与收敛域约束须在给定上下文假设下同时成立，并非全局已证事实。
 * 支持范围未知或收敛域无法表示时返回 Inconclusive；资源与取消错误原样传播。
 * 成功结果仅在保留的定义域内有效，规则集可不完备。
 */



/**
 * @brief 计算函数的 Laplace 变换 L{f(t)} = F(s).
 *
 * 算法:
 * 1. 线性性:L{af + bg} = aL{f} + bL{g}
 * 2. 内置规则匹配:识别已支持的变换对(t^n, e^at, sin, cos, sinh, cosh)
 * 3. 位移定理:L{e^(at)*f(t)} = F(s-a)
 * 4. 当前规则集之外的表达式返回 Inconclusive
 *
 * 多项式次数与阶乘系数保持为精确 BigInt。
 * 阶乘计算消耗上下文的步数和整数位数预算；耗尽时返回 checked 阶乘错误，保持精确语义。
 *
 * @param[in] f 时域函数表达式
 * @param[in] t 时域变量名
 * @param[in] s 频域变量名
 * @return Laplace 变换结果 F(s)
 */
LMCAS_API TransformEngineResult laplace_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算 Laplace 变换,并返回条件/ROC 容器.
 */
LMCAS_API TransformEngineResult laplace_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s);

/**
 * @brief 计算逆 Laplace 变换 L-¹{F(s)} = f(t).
 *
 * 算法:
 * 1. 线性性:L-¹{aF + bG} = aL-¹{F} + bL-¹{G}
 * 2. 部分分式分解 F(s)
 * 3. 对每个部分分式项查表求逆变换
 * 4. 重极点处理:L-¹{1/(s-a)^n} = t^(n-1)*e^(at)/(n-1)!
 * 5. 逆位移定理:L-¹{F(s-a)} = e^(at)*f(t)
 * 6. 当前规则集之外的表达式返回 Inconclusive
 *
 * @param[in] F 频域函数表达式
 * @param[in] s 频域变量名
 * @param[in] t 时域变量名
 * @return 逆 Laplace 变换结果 f(t)
 */
LMCAS_API TransformEngineResult inverse_laplace_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& s,
    const std::string& t,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算逆 Laplace 变换,并返回条件/ROC 容器.
 */
LMCAS_API TransformEngineResult inverse_laplace_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& s,
    const std::string& t);


/**
 * @brief 计算函数的 Fourier 变换 F{f(t)} = F(omega).
 *
 * 指数衰减率须由上下文假设证明为实数；无约束符号参数须显式声明为实数。
 * exp(-a*abs(t)) 的变换为 2*a/(a*a+omega*omega)，保留所有与 t 无关的因子；
 * a 的正性未知时附带 a > 0 条件。
 * 已知非正、复数或实性未证的衰减率返回 Inconclusive。
 * 线性组合与缩放保留各已求值项的条件。
 * @param[in] f 时域函数表达式
 * @param[in] t 时域变量名
 * @param[in] omega 频域变量名
 * @return Fourier 变换结果
 */
LMCAS_API TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算 Fourier 变换,并返回条件/ROC 容器.
 */
LMCAS_API TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega);

/**
 * @brief 计算逆 Fourier 变换 F-¹{F(omega)} = f(t).
 * @param[in] F 频域函数表达式
 * @param[in] omega 频域变量名
 * @param[in] t 时域变量名
 * @return 逆 Fourier 变换结果
 */
LMCAS_API TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算逆 Fourier 变换,并返回条件/ROC 容器.
 */
LMCAS_API TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t);

/**
 * @brief 计算两个函数的卷积 (f * g)(t) = integralf(tau)g(t-tau)dtau.
 * @param[in] f 第一个函数
 * @param[in] g 第二个函数
 * @param[in] var 卷积变量名
 * @return 卷积结果
 */
LMCAS_API TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算卷积,并显式报告无效输入.
 */
LMCAS_API TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var);


/**
 * @brief 计算序列的 Z 变换 Z{f[n]} = F(z).
 * @param[in] f_n 序列通项表达式
 * @param[in] n 序列指标变量名
 * @param[in] z Z 域变量名
 * @return Z 变换结果
 */
LMCAS_API TransformEngineResult z_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算 Z 变换,并返回条件/ROC 容器.
 */
LMCAS_API TransformEngineResult z_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z);

} // namespace LMCAS
