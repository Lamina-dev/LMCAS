/**
 * @file symbolic_geometry.hpp
 * @brief 解析几何应用：旋转体体积、弧长公式。
 */
#pragma once
#include "symbolic.hpp"
#include <string>
#include <memory>

namespace lamina {

/**
 * @brief 计算函数 f(x) 绕 x 轴旋转所得旋转体的体积
 * @param fx 被旋转的函数表达式 f(x)
 * @param a 积分下限
 * @param b 积分上限
 * @return 旋转体体积的符号表达式
 */
std::shared_ptr<SymbolicExpr> volume_of_revolution_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

/**
 * @brief 计算函数 f(x) 在 [a, b] 上的弧长
 * @param fx 曲线函数表达式 f(x)
 * @param a 区间左端点
 * @param b 区间右端点
 * @return 弧长的符号表达式
 */
std::shared_ptr<SymbolicExpr> arc_length_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

/**
 * @brief 计算函数 f(y) 绕 y 轴旋转所得旋转体的体积
 * @param fy 被旋转的函数表达式 f(y)
 * @param a 积分下限
 * @param b 积分上限
 * @return 旋转体体积的符号表达式
 */
std::shared_ptr<SymbolicExpr> volume_of_revolution_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

/**
 * @brief 计算函数 f(y) 在 [a, b] 上的弧长
 * @param fy 曲线函数表达式 f(y)
 * @param a 区间左端点
 * @param b 区间右端点
 * @return 弧长的符号表达式
 */
std::shared_ptr<SymbolicExpr> arc_length_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

}
