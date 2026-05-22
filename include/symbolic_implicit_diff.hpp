#include "symbolic.hpp"
#include "visitors/differentiation_visitor.hpp"
#include <memory>
#include <string>

namespace lamina {

std::shared_ptr<SymbolicExpr> implicit_diff(
    std::shared_ptr<SymbolicExpr> F,
    const std::string& x,
    const std::string& y
) {

    DifferentiationVisitor d_dx(x);
    F->root->accept(d_dx);
    auto dFdx = d_dx.result;

    DifferentiationVisitor d_dy(x, y);
    F->root->accept(d_dy);
    auto dFdy = d_dy.result;

    auto minus = SymbolicExpr::number(-1);
    auto num = SymbolicExpr::multiply(minus, std::make_shared<SymbolicExpr>(dFdx));
    auto denom = std::make_shared<SymbolicExpr>(dFdy);
    return SymbolicExpr::divide(num, denom);
}
}
