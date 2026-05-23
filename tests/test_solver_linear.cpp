#include "../include/solver.hpp"
#include <iostream>
#include <cassert>

using namespace lamina;

SymbolicExpr create_var(const std::string& name) {
    return SymbolicExpr(SymbolicFactory::create_variable(name));
}

SymbolicExpr create_num(int n) {
    return SymbolicExpr(SymbolicFactory::create_number(BigInt(n)));
}

SymbolicExpr operator+(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(a.root);
    ops.push_back(b.root);
    return SymbolicExpr(SymbolicFactory::create_add(ops));
}

SymbolicExpr operator-(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
    ops.push_back(b.root);
    auto neg = SymbolicFactory::create_multiply(ops);

    std::vector<std::shared_ptr<SymbolicNode>> aops;
    aops.push_back(a.root);
    aops.push_back(neg);
    return SymbolicExpr(SymbolicFactory::create_add(aops));
}

SymbolicExpr operator*(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(a.root);
    ops.push_back(b.root);
    return SymbolicExpr(SymbolicFactory::create_multiply(ops));
}

void test_linear_solver_2x2() {
    std::cout << "Testing linear solver 2x2..." << std::endl;

    auto x = create_var("x");
    auto y = create_var("y");
    auto three = create_num(3);
    auto one = create_num(1);

    auto eq1 = (x + y) - three;
    auto eq2 = (x - y) - one;

    std::vector<SymbolicExpr> equations = {eq1, eq2};
    std::vector<std::string> variables = {"x", "y"};

    auto result = Solver::solve_linear_system(equations, variables);

    if (result.count("x")) {

        std::cout << "x = " << result["x"].to_string() << std::endl;
    }
    if (result.count("y")) {
        std::cout << "y = " << result["y"].to_string() << std::endl;
    }

    auto x_val = result["x"];
    auto y_val = result["y"];

    assert(x_val.is_number() && x_val.get_number().index() == 1 && std::get<BigInt>(x_val.get_number()).to_int() == 2);
    assert(y_val.is_number() && y_val.get_number().index() == 1 && std::get<BigInt>(y_val.get_number()).to_int() == 1);

    std::cout << "2x2 test passed." << std::endl;
}

void test_linear_solver_3x3() {
     std::cout << "Testing linear solver 3x3..." << std::endl;

    auto x = create_var("x");
    auto y = create_var("y");
    auto z = create_var("z");
    auto num6 = create_num(6);
    auto num1 = create_num(1);
    auto num2 = create_num(2);

    auto eq1 = (x + y + z) - num6;
    auto eq2 = (create_num(2)*x + y - z) - num1;
    auto eq3 = (x - y + z) - num2;

    std::vector<SymbolicExpr> equations = {eq1, eq2, eq3};
    std::vector<std::string> variables = {"x", "y", "z"};

    auto result = Solver::solve_linear_system(equations, variables);

    assert(result.count("x"));
    assert(result.count("y"));
    assert(result.count("z"));

    std::cout << "x = " << result["x"].to_string() << std::endl;
    std::cout << "y = " << result["y"].to_string() << std::endl;
    std::cout << "z = " << result["z"].to_string() << std::endl;

}

int main() {
    try {
        test_linear_solver_2x2();
        test_linear_solver_3x3();
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
