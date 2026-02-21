#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {
// 复数符号类
struct ComplexSymbolic {
    std::shared_ptr<SymbolicExpr> real;
    std::shared_ptr<SymbolicExpr> imag;
};

// 构造复数 z = x + iy
inline ComplexSymbolic make_complex(
    std::shared_ptr<SymbolicExpr> real,
    std::shared_ptr<SymbolicExpr> imag
) {
    return ComplexSymbolic{real, imag};
}

// 四则运算
ComplexSymbolic complex_add(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_sub(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_mul(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_div(const ComplexSymbolic& a, const ComplexSymbolic& b);

// 共轭
ComplexSymbolic complex_conj(const ComplexSymbolic& z);

// 模长
std::shared_ptr<SymbolicExpr> complex_abs(const ComplexSymbolic& z);

// 辐角 Arg
std::shared_ptr<SymbolicExpr> complex_arg(const ComplexSymbolic& z);

// 指数形式 re^{iθ}
ComplexSymbolic complex_exp_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

// 三角形式 r(cosθ + i sinθ)
ComplexSymbolic complex_trig_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

// 复数方程求解 z^n = c
std::vector<ComplexSymbolic> solve_complex_nth_root(
    std::shared_ptr<SymbolicExpr> c,
    int n
);

// 二次复数方程 az^2 + bz + c = 0
std::vector<ComplexSymbolic> solve_complex_quadratic(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c
);

// 复数轨迹
// |z - a| = r 圆
// |z - a| = |z - b| 中垂线
std::shared_ptr<SymbolicExpr> complex_locus_circle(
    const ComplexSymbolic& a,
    std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var = "z"
);
std::shared_ptr<SymbolicExpr> complex_locus_perpendicular_bisector(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    const std::string& z_var = "z"
);

}
