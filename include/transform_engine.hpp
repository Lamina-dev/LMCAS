/**
 * @file transform_engine.hpp
 * @brief 积分变换引擎：Laplace 变换、逆 Laplace 变换、Fourier 变换、Z 变换。
 *
 * 提供基于变换表查找、线性性质和位移定理的符号积分变换计算。
 * 无法求得闭合形式时返回未求值的 TransformNode。
 */
#pragma once

#include "symbolic_ast.hpp"
#include <memory>
#include <string>

class SymbolicExpr;

#include <vector>
#include <functional>

namespace lamina {

// ============================================================
// 变换表 (Transform Table)
// ============================================================

/**
 * @brief 变换表条目，存储一对时域/频域对应关系。
 */
struct TransformTableEntry {
    std::string name;                          ///< 条目名称
    std::shared_ptr<SymbolicExpr> time_domain; ///< 时域表达式
    std::shared_ptr<SymbolicExpr> freq_domain; ///< 频域表达式
};

/**
 * @brief Laplace 变换表，管理已知变换对。
 *
 * 内置常见函数的 Laplace 变换对：
 * - 多项式 t^n ↔ n!/s^(n+1)
 * - 指数 e^(at) ↔ 1/(s-a)
 * - 正弦 sin(at) ↔ a/(s²+a²)
 * - 余弦 cos(at) ↔ s/(s²+a²)
 * - 双曲正弦 sinh(at) ↔ a/(s²-a²)
 * - 双曲余弦 cosh(at) ↔ s/(s²-a²)
 * - 阶跃函数 u(t) ↔ 1/s
 */
class LAMINA_API TransformTable {
public:
    TransformTable();

    /**
     * @brief 获取所有变换对条目。
     * @return 变换对列表的常引用
     */
    const std::vector<TransformTableEntry>& entries() const { return entries_; }

    /**
     * @brief 添加自定义变换对。
     * @param[in] entry 变换表条目
     */
    void add_entry(TransformTableEntry entry);

private:
    std::vector<TransformTableEntry> entries_;
    void init_laplace_pairs();
};

// ============================================================
// Laplace 变换 (Requirements 37, 86)
// ============================================================

/**
 * @brief 计算函数的 Laplace 变换 ℒ{f(t)} = F(s)。
 *
 * 算法：
 * 1. 线性性：ℒ{af + bg} = aℒ{f} + bℒ{g}
 * 2. 变换表查找：匹配已知变换对（t^n, e^at, sin, cos, sinh, cosh）
 * 3. 位移定理：ℒ{e^(at)·f(t)} = F(s-a)
 * 4. 无法求得闭合形式时返回未求值的 LaplaceNode
 *
 * @param[in] f 时域函数表达式
 * @param[in] t 时域变量名
 * @param[in] s 频域变量名
 * @return Laplace 变换结果 F(s)
 */
LAMINA_API std::shared_ptr<SymbolicExpr> laplace_transform(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s);

/**
 * @brief 计算逆 Laplace 变换 ℒ⁻¹{F(s)} = f(t)。
 *
 * 算法：
 * 1. 线性性：ℒ⁻¹{aF + bG} = aℒ⁻¹{F} + bℒ⁻¹{G}
 * 2. 部分分式分解 F(s)
 * 3. 对每个部分分式项查表求逆变换
 * 4. 重极点处理：ℒ⁻¹{1/(s-a)^n} = t^(n-1)·e^(at)/(n-1)!
 * 5. 逆位移定理：ℒ⁻¹{F(s-a)} = e^(at)·f(t)
 * 6. 无法求得闭合形式时返回未求值的逆 Laplace 节点
 *
 * @param[in] F 频域函数表达式
 * @param[in] s 频域变量名
 * @param[in] t 时域变量名
 * @return 逆 Laplace 变换结果 f(t)
 */
LAMINA_API std::shared_ptr<SymbolicExpr> inverse_laplace(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& s,
    const std::string& t);

// ============================================================
// Fourier 变换 (Requirements 63, 64)
// ============================================================

/**
 * @brief 计算函数的 Fourier 变换 F{f(t)} = F(ω)。
 * @param[in] f 时域函数表达式
 * @param[in] t 时域变量名
 * @param[in] omega 频域变量名
 * @return Fourier 变换结果
 */
LAMINA_API std::shared_ptr<SymbolicExpr> fourier_transform(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega);

/**
 * @brief 计算逆 Fourier 变换 F⁻¹{F(ω)} = f(t)。
 * @param[in] F 频域函数表达式
 * @param[in] omega 频域变量名
 * @param[in] t 时域变量名
 * @return 逆 Fourier 变换结果
 */
LAMINA_API std::shared_ptr<SymbolicExpr> inverse_fourier_transform(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t);

/**
 * @brief 计算两个函数的卷积 (f * g)(t) = ∫f(τ)g(t-τ)dτ。
 * @param[in] f 第一个函数
 * @param[in] g 第二个函数
 * @param[in] var 卷积变量名
 * @return 卷积结果
 */
LAMINA_API std::shared_ptr<SymbolicExpr> convolve(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var);

// ============================================================
// Z 变换 (Requirement 65)
// ============================================================

/**
 * @brief 计算序列的 Z 变换 Z{f[n]} = F(z)。
 * @param[in] f_n 序列通项表达式
 * @param[in] n 序列指标变量名
 * @param[in] z Z 域变量名
 * @return Z 变换结果
 */
LAMINA_API std::shared_ptr<SymbolicExpr> z_transform(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z);

} // namespace lamina
