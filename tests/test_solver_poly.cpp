#include "../include/solver.hpp"
#include <iostream>
#include <cassert>

using namespace lamina;

SymbolicExpr var(const std::string& name) {
    return SymbolicExpr(SymbolicFactory::create_variable(name));
}

SymbolicExpr num(int n) {
    return SymbolicExpr(SymbolicFactory::create_number(BigInt(n)));
}

SymbolicExpr operator+(const SymbolicExpr& a, const SymbolicExpr& b) {

    std::vector<std::shared_ptr<SymbolicNode>> ops = {a.root, b.root};
    return SymbolicExpr(SymbolicFactory::create_add(ops));
}

SymbolicExpr operator*(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops = {a.root, b.root};
    return SymbolicExpr(SymbolicFactory::create_multiply(ops));
}

SymbolicExpr operator-(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops = {SymbolicFactory::create_number(BigInt(-1)), b.root};
    auto neg = SymbolicFactory::create_multiply(ops);
    std::vector<std::shared_ptr<SymbolicNode>> aops = {a.root, neg};
    return SymbolicExpr(SymbolicFactory::create_add(aops));
}

SymbolicExpr pow(const SymbolicExpr& a, int n) {
    return SymbolicExpr(SymbolicFactory::create_power(a.root, SymbolicFactory::create_number(BigInt(n))));
}

void test_groebner_basis_simple() {
    std::cout << "Testing Groebner Basis simple..." << std::endl;

    auto x = var("x");
    auto y = var("y");

    auto f1 = x + y;
    auto f2 = x - y;

    std::vector<SymbolicExpr> F = {f1, f2};
    std::vector<std::string> vars = {"x", "y"};

    auto G = Solver::groebner_basis(F, vars);

    std::cout << "Basis size: " << G.size() << std::endl;
    for (const auto& g : G) {
        std::cout << "  " << g.to_string() << std::endl;
    }

}

void test_groebner_basis_circle() {
     std::cout << "Testing Groebner Basis circle line..." << std::endl;

     auto x = var("x");
     auto y = var("y");

     auto f1 = (pow(x, 2) + pow(y, 2)) - num(1);
     auto f2 = x - y;

     std::vector<SymbolicExpr> F = {f1, f2};
     std::vector<std::string> vars = {"x", "y"};

     auto G = Solver::groebner_basis(F, vars);

     std::cout << "Basis size: " << G.size() << std::endl;
     for (const auto& g : G) {
         std::cout << "  " << g.to_string() << std::endl;
     }
}

int main() {
    try {
        test_groebner_basis_simple();
        test_groebner_basis_circle();
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
