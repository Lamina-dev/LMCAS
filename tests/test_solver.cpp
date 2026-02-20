#include "../include/symbolic.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

void print_sols(const std::string& eq, const std::vector<std::shared_ptr<SymbolicExpr>>& sols) {
    std::cout << eq << " -> { ";
    for(size_t i=0; i<sols.size(); ++i) {
        
        auto sim = sols[i]->simplify();
        std::cout << sim->to_string();
        if(i < sols.size()-1) std::cout << ", ";
    }
    std::cout << " }\n";
}

int main() {
    auto x = SymbolicExpr::variable("x");
    auto n2 = SymbolicExpr::number(2);
    auto n6 = SymbolicExpr::number(6);
    auto n3 = SymbolicExpr::number(3);
    auto nm9 = SymbolicExpr::number(-9);
    auto nm5 = SymbolicExpr::number(-5);
    auto nm4 = SymbolicExpr::number(-4);
    auto n1 = SymbolicExpr::number(1);
    
    
    auto expr1 = SymbolicExpr::add(SymbolicExpr::multiply(n2, x), n6);
    auto sol1 = SymbolicExpr::solve(expr1, "x");
    print_sols("2x + 6 = 0", sol1);
    
    
    auto expr2 = SymbolicExpr::add(SymbolicExpr::multiply(n3, x), nm9);
    auto sol2 = SymbolicExpr::solve(expr2, "x");
    print_sols("3x - 9 = 0", sol2);

    
    auto expr3 = SymbolicExpr::add(x, nm5);
    auto sol3 = SymbolicExpr::solve(expr3, "x");
    print_sols("x - 5 = 0", sol3);
    
    
    auto x2 = SymbolicExpr::power(x, n2);
    auto expr4 = SymbolicExpr::add(x2, nm4);
    auto sol4 = SymbolicExpr::solve(expr4, "x");
    print_sols("x^2 - 4 = 0", sol4);
    
    
    auto term2 = SymbolicExpr::multiply(n2, x);
    auto part1 = SymbolicExpr::add(x2, term2);
    auto expr5 = SymbolicExpr::add(part1, n1);
    auto sol5 = SymbolicExpr::solve(expr5, "x");
    print_sols("x^2 + 2x + 1 = 0", sol5);

    return 0;
}
