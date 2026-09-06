#include "../include/solver.hpp"
#include "test_common.hpp"
#include <iostream>

using namespace LMCAS;

SymbolicExpr create_var(const std::string& name) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

SymbolicExpr create_num(int n) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_number(BigInt(n)));
}

static SymbolicExpr operator+(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(LMCAS::detail::node(a));
    ops.push_back(LMCAS::detail::node(b));
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(ops));
}

static SymbolicExpr operator-(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
    ops.push_back(LMCAS::detail::node(b));
    auto neg = SymbolicFactory::create_multiply(ops);

    std::vector<std::shared_ptr<const SymbolicNode>> aops;
    aops.push_back(LMCAS::detail::node(a));
    aops.push_back(neg);
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(aops));
}

static SymbolicExpr operator*(const SymbolicExpr& a, const SymbolicExpr& b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(LMCAS::detail::node(a));
    ops.push_back(LMCAS::detail::node(b));
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
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

        std::cout << "x = " << result.at("x").to_string() << std::endl;
    }
    if (result.count("y")) {
        std::cout << "y = " << result.at("y").to_string() << std::endl;
    }

    auto x_val = result.at("x");
    auto y_val = result.at("y");

    EXPECT_TRUE(x_val.is_number() && x_val.get_number().index() == 1 &&
                    std::get<BigInt>(x_val.get_number()).to_int() == 2,
                "2x2 linear solver returns x = 2");
    EXPECT_TRUE(y_val.is_number() && y_val.get_number().index() == 1 &&
                    std::get<BigInt>(y_val.get_number()).to_int() == 1,
                "2x2 linear solver returns y = 1");

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

    EXPECT_TRUE(result.count("x") == 1, "3x3 linear solver returns x");
    EXPECT_TRUE(result.count("y") == 1, "3x3 linear solver returns y");
    EXPECT_TRUE(result.count("z") == 1, "3x3 linear solver returns z");

    std::cout << "x = " << result.at("x").to_string() << std::endl;
    std::cout << "y = " << result.at("y").to_string() << std::endl;
    std::cout << "z = " << result.at("z").to_string() << std::endl;

}

int main() {
    try {
        test_linear_solver_2x2();
        test_linear_solver_3x3();
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
