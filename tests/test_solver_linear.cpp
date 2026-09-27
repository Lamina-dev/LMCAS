#include "solver.hpp"
#include "test_common.hpp"

using namespace LMCAS;

SymbolicExpr create_var(const std::string &name) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

SymbolicExpr create_num(int n) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_number(BigInt(n)));
}

static SymbolicExpr operator+(const SymbolicExpr &a, const SymbolicExpr &b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(LMCAS::detail::node(a));
    ops.push_back(LMCAS::detail::node(b));
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(ops));
}

static SymbolicExpr operator-(const SymbolicExpr &a, const SymbolicExpr &b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
    ops.push_back(LMCAS::detail::node(b));
    auto neg = SymbolicFactory::create_multiply(ops);

    std::vector<std::shared_ptr<const SymbolicNode>> aops;
    aops.push_back(LMCAS::detail::node(a));
    aops.push_back(neg);
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(aops));
}

static SymbolicExpr operator*(const SymbolicExpr &a, const SymbolicExpr &b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(LMCAS::detail::node(a));
    ops.push_back(LMCAS::detail::node(b));
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
}

TEST(LmcasSolverLinear, LinearSolver2x2) {
    auto x = create_var("x");
    auto y = create_var("y");
    auto three = create_num(3);
    auto one = create_num(1);

    auto eq1 = (x + y) - three;
    auto eq2 = (x - y) - one;

    std::vector<SymbolicExpr> equations = {eq1, eq2};
    std::vector<std::string> variables = {"x", "y"};

    auto result = Solver::solve_linear_system(equations, variables);

    auto x_val = result.at("x");
    auto y_val = result.at("y");

    EXPECT_TRUE((x_val.is_number() && std::holds_alternative<BigInt>(x_val.get_number()) &&
                 std::get<BigInt>(x_val.get_number()) == BigInt(2)))
        << ("2x2 linear solver returns exact x = 2");
    EXPECT_TRUE((y_val.is_number() && std::holds_alternative<BigInt>(y_val.get_number()) &&
                 std::get<BigInt>(y_val.get_number()) == BigInt(1)))
        << ("2x2 linear solver returns exact y = 1");
}

TEST(LmcasSolverLinear, LinearSolver3x3) {
    auto x = create_var("x");
    auto y = create_var("y");
    auto z = create_var("z");
    auto num6 = create_num(6);
    auto num1 = create_num(1);
    auto num2 = create_num(2);

    auto eq1 = (x + y + z) - num6;
    auto eq2 = (create_num(2) * x + y - z) - num1;
    auto eq3 = (x - y + z) - num2;

    std::vector<SymbolicExpr> equations = {eq1, eq2, eq3};
    std::vector<std::string> variables = {"x", "y", "z"};

    auto result = Solver::solve_linear_system(equations, variables);

    EXPECT_TRUE((result.count("x") == 1)) << ("3x3 linear solver returns x");
    EXPECT_TRUE((result.count("y") == 1)) << ("3x3 linear solver returns y");
    EXPECT_TRUE((result.count("z") == 1)) << ("3x3 linear solver returns z");

}
