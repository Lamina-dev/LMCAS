#include "../include/symbolic_ode.hpp"
#include "../include/symbolic.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <cmath>
#include <memory>
#include <string>

namespace lamina {

std::shared_ptr<SymbolicExpr> solve_separable_ode(
    std::shared_ptr<SymbolicExpr> rhs,
    const std::string& x,
    const std::string& y
) {

    auto inv_y = SymbolicExpr::divide(SymbolicExpr::number(1), rhs);
    auto int_y = inv_y->integrate(y);
    auto int_x = SymbolicExpr::number(1)->integrate(x);

    return SymbolicExpr::add(int_y, SymbolicExpr::multiply(SymbolicExpr::number(-1), int_x));
}

std::shared_ptr<SymbolicExpr> solve_linear1_ode(
    std::shared_ptr<SymbolicExpr> Px,
    std::shared_ptr<SymbolicExpr> Qx,
    const std::string& x,
    const std::string& y
) {

    auto intP = Px->integrate(x);
    auto mu = SymbolicExpr::exp(intP);

    auto Qmu = SymbolicExpr::multiply(Qx, mu);
    auto intQmu = Qmu->integrate(x);

    auto C = SymbolicExpr::variable("C");
    auto num = SymbolicExpr::add(intQmu, C);
    return SymbolicExpr::divide(num, mu);
}

std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y
) {

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

        auto exp_part = SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(real), SymbolicExpr::variable(x)));
        auto cos_part = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        auto sin_part = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        yh = SymbolicExpr::multiply(exp_part, SymbolicExpr::add(SymbolicExpr::multiply(C1, cos_part), SymbolicExpr::multiply(C2, sin_part)));
    }

    if (!fx->is_zero()) {

        return SymbolicExpr::add(yh, SymbolicExpr::number(0));
    } else {
        return yh;
    }
}

}
