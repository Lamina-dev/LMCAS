#include "solver.hpp"
#include "test_common.hpp"

using namespace LMCAS;

SymbolicExpr var(const std::string &name) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

SymbolicExpr num(int n) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_number(BigInt(n)));
}

static SymbolicExpr operator+(const SymbolicExpr &a, const SymbolicExpr &b) {

    std::vector<std::shared_ptr<const SymbolicNode>> ops = {LMCAS::detail::node(a), LMCAS::detail::node(b)};
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(ops));
}

static SymbolicExpr operator-(const SymbolicExpr &a, const SymbolicExpr &b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops = {SymbolicFactory::create_number(BigInt(-1)), LMCAS::detail::node(b)};
    auto neg = SymbolicFactory::create_multiply(ops);
    std::vector<std::shared_ptr<const SymbolicNode>> aops = {LMCAS::detail::node(a), neg};
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_add(aops));
}

SymbolicExpr pow(const SymbolicExpr &a, int n) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_power(LMCAS::detail::node(a), SymbolicFactory::create_number(BigInt(n))));
}

TEST(LmcasSolverPoly, GroebnerBasisSimple) {
    auto x = var("x");
    auto y = var("y");

    auto f1 = x + y;
    auto f2 = x - y;

    std::vector<SymbolicExpr> F = {f1, f2};
    std::vector<std::string> vars = {"x", "y"};

    auto G = Solver::groebner_basis(F, vars);

    ASSERT_FALSE(G.empty());
    for (const auto &generator : F) {
        EXPECT_TRUE(Solver::ideal_membership(generator, G, vars));
    }
}

TEST(LmcasSolverPoly, GroebnerBasisCircle) {
    auto x = var("x");
    auto y = var("y");

    auto f1 = (pow(x, 2) + pow(y, 2)) - num(1);
    auto f2 = x - y;

    std::vector<SymbolicExpr> F = {f1, f2};
    std::vector<std::string> vars = {"x", "y"};

    auto G = Solver::groebner_basis(F, vars);

    ASSERT_FALSE(G.empty());
    for (const auto &generator : F) {
        EXPECT_TRUE(Solver::ideal_membership(generator, G, vars));
    }
}
