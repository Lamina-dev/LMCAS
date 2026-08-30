#include "../include/solver.hpp"
#include "test_common.hpp"
#include <iostream>

using namespace lamina;

SymbolicExpr var(const std::string& name) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

SymbolicExpr num(int n) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(n)));
}

static SymbolicExpr operator+(const SymbolicExpr& a, const SymbolicExpr& b) {

    std::vector<std::shared_ptr<const SymbolicNode>> ops = {lamina::detail::node(a), lamina::detail::node(b)};
    return lamina::detail::expression_from_node(SymbolicFactory::create_add(ops));
}

static SymbolicExpr operator*(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops = {lamina::detail::node(a), lamina::detail::node(b)};
    return lamina::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
}

static SymbolicExpr operator-(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops = {SymbolicFactory::create_number(BigInt(-1)), lamina::detail::node(b)};
    auto neg = SymbolicFactory::create_multiply(ops);
    std::vector<std::shared_ptr<const SymbolicNode>> aops = {lamina::detail::node(a), neg};
    return lamina::detail::expression_from_node(SymbolicFactory::create_add(aops));
}

SymbolicExpr pow(const SymbolicExpr& a, int n) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_power(lamina::detail::node(a), SymbolicFactory::create_number(BigInt(n))));
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
    EXPECT_TRUE(!G.empty(), "simple Groebner basis is non-empty");
    for (const auto& g : G) {
        std::cout << "  " << g.to_string() << std::endl;
        EXPECT_TRUE(!g.to_string().empty(), "simple Groebner basis element is printable");
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
     EXPECT_TRUE(!G.empty(), "circle-line Groebner basis is non-empty");
     for (const auto& g : G) {
         std::cout << "  " << g.to_string() << std::endl;
         EXPECT_TRUE(!g.to_string().empty(), "circle-line Groebner basis element is printable");
     }
}

int main() {
    try {
        test_groebner_basis_simple();
        test_groebner_basis_circle();
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
