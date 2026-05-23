#include "../include/symbolic_ode.hpp"
#include "../include/symbolic.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <cmath>
#include <memory>
#include <stdexcept>
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

    // Guard a == 0: equation degenerates to first-order (or constant) form.
    int a_is_zero = 0;
    lmmc_double_nearly_equal_tol(a, 0.0, 1e-12, 1e-12, &a_is_zero);
    if (a_is_zero) {
        int b_is_zero = 0;
        lmmc_double_nearly_equal_tol(b, 0.0, 1e-12, 1e-12, &b_is_zero);
        if (b_is_zero) {
            // a == 0 and b == 0: not a proper second-order ODE.
            throw std::invalid_argument(
                "solve_linear2_ode: leading and first-derivative coefficients are both zero");
        }
        // a == 0, b != 0: degenerates to b*y' + c*y = f(x), i.e. y' + (c/b)*y = f/b.
        auto Px = SymbolicExpr::number(c / b);
        auto Qx = fx->is_zero()
                      ? SymbolicExpr::number(0)
                      : SymbolicExpr::divide(fx, SymbolicExpr::number(b));
        return solve_linear1_ode(Px, Qx, x, y);
    }

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
        // Particular-solution computation for non-homogeneous case is not yet
        // implemented. Fail loudly instead of returning the homogeneous solution
        // disguised as the full general solution.
        throw std::logic_error(
            "solve_linear2_ode: non-homogeneous case (f(x) != 0) is not implemented");
    }
    return yh;
}

}
