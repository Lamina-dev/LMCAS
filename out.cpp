#include <iostream>
#include <memory>
#include "../include/symbolic.hpp"

int main() {
    using Expr = std::shared_ptr<SymbolicExpr>;
    Expr expr = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)), SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::variable("x"))), SymbolicExpr::sin(SymbolicExpr::variable("x")));
    std::cout << "expr = " << expr->to_string() << "\n";
    auto simplified = expr->simplify();
    std::cout << "simplified = " << simplified->to_string() << "\n";
    // example: derivative with respect to x
    auto deriv = expr->differentiate("x");
    std::cout << "derivative = " << (deriv ? deriv->to_string() : "<null>") << "\n";
    return 0;
}