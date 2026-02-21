#include "symbolic.hpp"
#include "visitors/differentiation_visitor.hpp"
#include <memory>
#include <string>

namespace lamina {
// 隐函数微分接口
// 返回 dF/dx, 其中 F(x, y(x)) = 0，y 隐含 x
std::shared_ptr<SymbolicExpr> implicit_diff(
    std::shared_ptr<SymbolicExpr> F,
    const std::string& x,
    const std::string& y
) {
    // dF/dx + dF/dy * y' = 0  =>  y' = - (dF/dx) / (dF/dy)
    DifferentiationVisitor d_dx(x);
    F->root->accept(d_dx);
    auto dFdx = d_dx.result;

    DifferentiationVisitor d_dy(x, y); // 隐函数模式
    F->root->accept(d_dy);
    auto dFdy = d_dy.result;

    auto minus = SymbolicExpr::number(-1);
    auto num = SymbolicExpr::multiply(minus, std::make_shared<SymbolicExpr>(dFdx));
    auto denom = std::make_shared<SymbolicExpr>(dFdy);
    return SymbolicExpr::divide(num, denom);
}
} // namespace lamina
