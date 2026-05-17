#include "test_common.hpp"
#include "solver.hpp"

int main() {
    TEST_CASE("Solve Higher Degree Polynomial (RootOf)");
    {
        auto x = SymbolicExpr::variable("x");
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto eq = SymbolicExpr::add(x3, SymbolicExpr::number(-2));

        auto sols = SymbolicExpr::solve(eq, "x");
        EXPECT_TRUE(sols.size() == 3, "rootof count should be 3");
        for (const auto& s : sols) {
            EXPECT_CONTAINS(s->to_string(), {"rootof"}, "rootof token");
        }
    }

    TEST_CASE("Solve Linear+Exp (LambertW)");
    {
        auto x = SymbolicExpr::variable("x");
        auto eq = SymbolicExpr::add(x, SymbolicExpr::exp(x));

        auto sols = SymbolicExpr::solve(eq, "x");
        EXPECT_TRUE(sols.size() == 1, "lambertw solution count");
        if (!sols.empty()) {
            EXPECT_CONTAINS(sols[0]->to_string(), {"lambertw"}, "lambertw token");
        }
    }

    TEST_CASE("Solve Rational System (Denominator Filter)");
    {
        auto x = SymbolicExpr::variable("x");
        auto denom = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto frac = SymbolicExpr::divide(x, denom);
        auto eq = SymbolicExpr::add(frac, SymbolicExpr::number(-2));

        std::vector<SymbolicExpr> eqs = {*eq};
        auto sols = lamina::Solver::solve_polynomial_system(eqs, {"x"});
        EXPECT_TRUE(sols.size() == 1, "rational system solutions size");
        if (!sols.empty()) {
            auto x_val = std::make_shared<SymbolicExpr>(sols[0]["x"]);
            EXPECT_EQ_EXPR(x_val, SymbolicExpr::number(2), "rational system x=2");
        }
    }

    return TEST_REPORT();
}
