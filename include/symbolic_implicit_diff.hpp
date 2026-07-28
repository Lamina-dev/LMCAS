/**
 * @file symbolic_implicit_diff.hpp
 * @brief 隐函数微分 dy/dx = -F_x / F_y。
 */
#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

/**
 * @brief 对隐函数 F(x, y) = 0 求 dy/dx
 *
 * 利用隐函数定理：dy/dx = -F_x / F_y
 *
 * @param F 隐函数表达式 F(x, y)
 * @param x 自变量名
 * @param y 因变量名
 * @return dy/dx 的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> implicit_diff(
    std::shared_ptr<SymbolicExpr> F,
    const std::string& x,
    const std::string& y
);

}
