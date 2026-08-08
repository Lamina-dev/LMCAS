/**
 * @file calculus_utils.hpp
 * @brief 微积分工具函数：连续性判定、渐近线分析、对数微分、微分、全微分、反函数导数、反函数求解等。
 */
#pragma once

#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <utility>
#include <memory>

namespace lamina {


/**
 * @brief 连续性类型枚举。
 */
enum class ContinuityType {
    Continuous,  ///< 连续：左极限 = 右极限 = 函数值
    Removable,   ///< 可去间断点：左极限 = 右极限 ≠ 函数值（或函数无定义）
    Jump,        ///< 跳跃间断点：左极限 ≠ 右极限，但两侧极限均存在
    Essential    ///< 本性间断点：至少一侧极限不存在或为无穷
};

/**
 * @brief 判断函数在指定点的连续性类型。
 *
 * 通过计算左极限、右极限和函数值来确定连续性分类。
 *
 * @param[in] f     待分析的函数表达式
 * @param[in] var   自变量名
 * @param[in] point 待检测的点
 * @return 连续性类型
 */
LAMINA_API ContinuityType continuity_at(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point);


/**
 * @brief 渐近线分析结果。
 */
struct AsymptoteResult {
    std::vector<std::shared_ptr<SymbolicExpr>> vertical;     ///< 垂直渐近线 x = a
    std::vector<std::shared_ptr<SymbolicExpr>> horizontal;   ///< 水平渐近线 y = L
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                          std::shared_ptr<SymbolicExpr>>> oblique; ///< 斜渐近线 (斜率, 截距)
};

using AsymptoteAnalysisResult = Result<AsymptoteResult>;
using SymbolicExprResult = Result<std::shared_ptr<SymbolicExpr>>;
using SymbolicExprVectorResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

/**
 * @brief 计算函数的渐近线。
 *
 * 垂直渐近线：分母为零且极限为 ±∞ 的点。
 * 水平渐近线：x→±∞ 时的极限值。
 * 斜渐近线：当无水平渐近线时，计算 slope = lim(f/x) 和 intercept = lim(f - slope·x)。
 *
 * @param[in] f   待分析的函数表达式
 * @param[in] var 自变量名
 * @return 渐近线分析结果
 */
LAMINA_API AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算函数的渐近线，并显式报告无效输入。
 */
LAMINA_API AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

LAMINA_API AsymptoteResult asymptotes(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);


/**
 * @brief 对数微分法：计算 f 关于 var 的导数。
 *
 * 利用公式 d/dx[f] = f · d/dx[ln(f)]，对 ln(f) 先应用对数化简规则
 * （乘积变求和、幂次提前）再求导，适用于 x^x 或多项乘积等复杂表达式。
 *
 * @param[in] f   待求导的表达式
 * @param[in] var 求导变量名
 * @return 对数微分结果表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> log_differentiate(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

/**
 * @brief 计算单变量微分 dy = f'(x)dx 中的系数 f'(x)。
 *
 * 返回 f 关于 var 的导数，即微分 df 中 dx 的系数。
 *
 * @param[in] f   函数表达式
 * @param[in] var 微分变量名
 * @return f'(var) 表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

/**
 * @brief 计算全微分 df = ∑(∂f/∂xᵢ)·dxᵢ。
 *
 * 对多元函数 f 关于每个变量求偏导，返回 (偏导数, 变量名) 对的列表，
 * 表示全微分中各项。
 *
 * @param[in] f    多元函数表达式
 * @param[in] vars 变量名列表
 * @return (∂f/∂xᵢ, "xᵢ") 对的向量
 */
LAMINA_API std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> total_differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars);

/**
 * @brief 计算反函数导数 (f⁻¹)'(point) = 1 / f'(f⁻¹(point))。
 *
 * Checked API 只在反函数候选唯一且 f'(x₀) 可构造、非零时成功。
 */
LAMINA_API SymbolicExprResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算反函数导数，并显式报告失败语义。
 */
LAMINA_API SymbolicExprResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point);

/**
 * @brief 计算反函数导数 (f⁻¹)'(point) = 1 / f'(f⁻¹(point))。
 *
 * Legacy wrapper 先求解 f(x) = point 得到第一个 x₀，再计算
 * 1/f'(x₀)。若 f'(x₀) = 0 则返回 nullptr。
 *
 * @param[in] f     函数表达式
 * @param[in] var   变量名
 * @param[in] point 求导点（y 值）
 * @return 反函数导数表达式，不存在时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> inverse_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point);

/**
 * @brief 求解反函数 f⁻¹(y)，即求解方程 f(x) = y 关于 x 的所有解。
 *
 * Checked API 使用 solution-set dispatcher；只在有限精确解或可证明空集
 * 时成功，超出支持域返回 `CasErrc::Inconclusive`。
 */
LAMINA_API SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y, ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解反函数，并显式报告失败语义。
 */
LAMINA_API SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y);

/**
 * @brief 求解反函数 f⁻¹(y)，即求解方程 f(x) = y 关于 x 的所有解。
 *
 * Legacy wrapper 利用现有求解器将 f(var) - y = 0 求解，返回所有分支。
 *
 * @param[in] f   函数表达式
 * @param[in] var 变量名
 * @param[in] y   目标值表达式
 * @return 所有解的列表
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> inverse_function(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y);


/**
 * @brief 计算显式曲线 y = f(x) 的曲率 κ = |f''| / (1 + f'²)^(3/2)
 * @param[in] f 曲线函数表达式 f(x)
 * @param[in] var 自变量名称
 * @return 曲率的符号表达式
 */
LAMINA_API SymbolicExprResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算显式曲线曲率，并显式报告无效输入。
 */
LAMINA_API SymbolicExprResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

LAMINA_API std::shared_ptr<SymbolicExpr> curvature(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

/**
 * @brief 计算参数曲线 (x(t), y(t)) 的曲率 κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
 * @param[in] x_t x 分量关于参数 t 的表达式
 * @param[in] y_t y 分量关于参数 t 的表达式
 * @param[in] t 参数变量名称
 * @return 曲率的符号表达式
 */
LAMINA_API SymbolicExprResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算参数曲线曲率，并显式报告无效输入和零速度。
 */
LAMINA_API SymbolicExprResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t);

LAMINA_API std::shared_ptr<SymbolicExpr> curvature_parametric(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t);


/**
 * @brief 计算函数 y = f(x) 的拐点（二阶导数为零的点）。
 * 
 * @param[in] f 曲线函数表达式 f(x)
 * @param[in] var 自变量名称
 * @return 拐点的 x 坐标列表
 */
LAMINA_API SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算拐点，并显式报告无效输入。
 */
LAMINA_API SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> inflection_points(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var);


/**
 * @brief 计算曲线 y = f(x) 绕 x 轴旋转所得旋转面的面积
 *
 * 公式: S = 2π ∫ₐᵇ |f(x)| · √(1 + f'(x)²) dx
 * Checked API 只接受可精确求值的符号定积分；超出支持域返回
 * `CasErrc::Inconclusive`。
 */
LAMINA_API SymbolicExprResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算绕 x 轴旋转面积，并显式报告失败语义。
 */
LAMINA_API SymbolicExprResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

/**
 * @brief 计算曲线 y = f(x) 绕 x 轴旋转所得旋转面的面积。
 *
 * Legacy wrapper 先尝试符号积分，若失败则回退到数值积分。
 *
 * @param[in] f 曲线函数表达式 f(x)
 * @param[in] var 自变量名称
 * @param[in] a 积分下限
 * @param[in] b 积分上限
 * @return 旋转面面积的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> surface_area_revolution_x(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

/**
 * @brief 计算曲线 y = f(x) 绕 y 轴旋转所得旋转面的面积
 *
 * 公式: S = 2π ∫ₐᵇ |x| · √(1 + f'(x)²) dx
 * Checked API 只接受可精确求值的符号定积分；超出支持域返回
 * `CasErrc::Inconclusive`。
 */
LAMINA_API SymbolicExprResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算绕 y 轴旋转面积，并显式报告失败语义。
 */
LAMINA_API SymbolicExprResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

/**
 * @brief 计算曲线 y = f(x) 绕 y 轴旋转所得旋转面的面积。
 *
 * Legacy wrapper 先尝试符号积分，若失败则回退到数值积分。
 *
 * @param[in] f 曲线函数表达式 f(x)
 * @param[in] var 自变量名称
 * @param[in] a 积分下限
 * @param[in] b 积分上限
 * @return 旋转面面积的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> surface_area_revolution_y(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

} // namespace lamina
