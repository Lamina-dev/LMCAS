#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {
// 一阶微分方程分离变量法求解：dy/dx = f(x, y)
// 返回 y(x) 的隐式表达式（如 F(y) = G(x) + C）
std::shared_ptr<SymbolicExpr> solve_separable_ode(
    std::shared_ptr<SymbolicExpr> rhs, // f(x, y)
    const std::string& x,
    const std::string& y
);

// 一阶线性微分方程：y' + P(x)y = Q(x)
// 返回通解表达式
std::shared_ptr<SymbolicExpr> solve_linear1_ode(
    std::shared_ptr<SymbolicExpr> Px, // P(x)
    std::shared_ptr<SymbolicExpr> Qx, // Q(x)
    const std::string& x,
    const std::string& y
);

// 二阶常系数齐次/非齐次微分方程 ay''+by'+cy=f(x)
// 返回通解表达式
std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx, // f(x)
    const std::string& x,
    const std::string& y
);

}
