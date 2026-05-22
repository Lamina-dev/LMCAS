#include "../include/symbolic_geometry.hpp"
#include "../include/symbolic.hpp"
#include <cmath>

namespace lamina {

std::shared_ptr<SymbolicExpr> volume_of_revolution_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    auto pi = SymbolicExpr::number(M_PI);
    auto integrand = SymbolicExpr::multiply(pi, SymbolicExpr::power(fx, SymbolicExpr::number(2)));
    auto integral = integrand->integrate("x");
    auto upper = integral->substitute("x", b);
    auto lower = integral->substitute("x", a);
    return SymbolicExpr::add(upper, SymbolicExpr::multiply(SymbolicExpr::number(-1), lower));
}

std::shared_ptr<SymbolicExpr> arc_length_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    auto dfdx = fx->differentiate("x");
    auto integrand = SymbolicExpr::sqrt(SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::power(dfdx, SymbolicExpr::number(2))));
    auto integral = integrand->integrate("x");
    auto upper = integral->substitute("x", b);
    auto lower = integral->substitute("x", a);
    return SymbolicExpr::add(upper, SymbolicExpr::multiply(SymbolicExpr::number(-1), lower));
}

std::shared_ptr<SymbolicExpr> volume_of_revolution_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    auto pi = SymbolicExpr::number(M_PI);
    auto integrand = SymbolicExpr::multiply(pi, SymbolicExpr::power(fy, SymbolicExpr::number(2)));
    auto integral = integrand->integrate("y");
    auto upper = integral->substitute("y", b);
    auto lower = integral->substitute("y", a);
    return SymbolicExpr::add(upper, SymbolicExpr::multiply(SymbolicExpr::number(-1), lower));
}

std::shared_ptr<SymbolicExpr> arc_length_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    auto dfdy = fy->differentiate("y");
    auto integrand = SymbolicExpr::sqrt(SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::power(dfdy, SymbolicExpr::number(2))));
    auto integral = integrand->integrate("y");
    auto upper = integral->substitute("y", b);
    auto lower = integral->substitute("y", a);
    return SymbolicExpr::add(upper, SymbolicExpr::multiply(SymbolicExpr::number(-1), lower));
}

}
