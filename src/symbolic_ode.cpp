#include "../include/symbolic_ode.hpp"
#include "../include/symbolic.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <cmath>
#include <memory>
#include <string>

namespace lamina {
// 分离变量法：dy/dx = f(x, y)
std::shared_ptr<SymbolicExpr> solve_separable_ode(
    std::shared_ptr<SymbolicExpr> rhs,
    const std::string& x,
    const std::string& y
) {
    // 目标：∫ 1/f_y dy = ∫ 1 dx
    // 1. 试图将 rhs 拆分为 X(x)*Y(y) 形式
    // 2. 若能分离，构造两边积分表达式
    // 这里只做形式化（不做自动分离），用户需传入已分离的表达式
    // 例：rhs = X(x)*Y(y) 传入 1/Y(y) dy = X(x) dx
    auto inv_y = SymbolicExpr::divide(SymbolicExpr::number(1), rhs);
    auto int_y = inv_y->integrate(y);
    auto int_x = SymbolicExpr::number(1)->integrate(x);
    // 返回 int_y = int_x + C
    // 用 AddNode 表示 int_y - int_x = C
    return SymbolicExpr::add(int_y, SymbolicExpr::multiply(SymbolicExpr::number(-1), int_x));
}

// 一阶线性微分方程 y' + P(x)y = Q(x)
std::shared_ptr<SymbolicExpr> solve_linear1_ode(
    std::shared_ptr<SymbolicExpr> Px,
    std::shared_ptr<SymbolicExpr> Qx,
    const std::string& x,
    const std::string& y
) {
    // 积分因子 μ(x) = exp(∫P(x)dx)
    auto intP = Px->integrate(x);
    auto mu = SymbolicExpr::exp(intP);
    // y*μ(x) = ∫Q(x)μ(x)dx + C
    auto Qmu = SymbolicExpr::multiply(Qx, mu);
    auto intQmu = Qmu->integrate(x);
    // 返回 y = (intQmu + C)/mu
    // 形式上 y = (∫Q(x)μ(x)dx + C)/μ(x)
    // 用 AddNode 表示 intQmu + C
    auto C = SymbolicExpr::variable("C");
    auto num = SymbolicExpr::add(intQmu, C);
    return SymbolicExpr::divide(num, mu);
}

// 二阶常系数齐次/非齐次 ay''+by'+cy=f(x)
std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y
) {
    // 齐次解：ay''+by'+cy=0
    // 特征方程：ar^2+br+c=0
    double D = b*b - 4*a*c;
    auto C1 = SymbolicExpr::variable("C1");
    auto C2 = SymbolicExpr::variable("C2");
    std::shared_ptr<SymbolicExpr> yh;
    int eq;
    lmmc_double_nearly_equal_tol(D, 0.0, 1e-12, 1e-12, &eq);
    if (!eq && D > 0) {
        double r1 = (-b + std::sqrt(D)) / (2*a);
        double r2 = (-b - std::sqrt(D)) / (2*a);
        yh = SymbolicExpr::add(
            SymbolicExpr::multiply(C1, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r1), SymbolicExpr::variable(x)))),
            SymbolicExpr::multiply(C2, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r2), SymbolicExpr::variable(x))))
        );
    } else if (eq) {
        double r = -b / (2*a);
        yh = SymbolicExpr::add(
            SymbolicExpr::multiply(C1, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r), SymbolicExpr::variable(x)))),
            SymbolicExpr::multiply(C2, SymbolicExpr::multiply(SymbolicExpr::variable(x), SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r), SymbolicExpr::variable(x)))))
        );
    } else {
        double real = -b / (2*a);
        double imag = std::sqrt(-D) / (2*a);
        // yh = e^{real x} (C1 cos(imag x) + C2 sin(imag x))
        auto exp_part = SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(real), SymbolicExpr::variable(x)));
        auto cos_part = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        auto sin_part = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        yh = SymbolicExpr::multiply(exp_part, SymbolicExpr::add(SymbolicExpr::multiply(C1, cos_part), SymbolicExpr::multiply(C2, sin_part)));
    }
    // 非齐次：常数变易法
    if (!fx->is_zero()) {
        // 只做形式化，实际需用常数变易法求 particular solution
        // 这里只返回 yh + Yp 形式，Yp = 0
        return SymbolicExpr::add(yh, SymbolicExpr::number(0));
    } else {
        return yh;
    }
}

} // namespace lamina
