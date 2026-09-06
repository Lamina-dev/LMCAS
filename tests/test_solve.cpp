#include "test_common.hpp"
#include "inequality_solver.hpp"
#include "symbolic_matrix.hpp"

using namespace LMCAS;

int main() {
    TEST_CASE("Solve System (Direct & Linear)");

    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");

        auto eq1 = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::add(y, SymbolicExpr::number(-5))
        );

        auto eq2 = SymbolicExpr::add(
            x,
            SymbolicExpr::add(SymbolicExpr::multiply(y, SymbolicExpr::number(-1)), SymbolicExpr::number(-1))
        );

        auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        EXPECT_TRUE(solutions.size() == 1, "Solutions size should be 1");
        if (solutions.size() > 0) {
            EXPECT_EQ_EXPR(solutions[0].at("x"), SymbolicExpr::number(2), "x should be 2");
            EXPECT_EQ_EXPR(solutions[0].at("y"), SymbolicExpr::number(1), "y should be 1");
        }
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto a = SymbolicExpr::variable("a");

        auto eq1 = SymbolicExpr::add(x, SymbolicExpr::multiply(a, y));
        auto eq2 = SymbolicExpr::add(x, SymbolicExpr::add(SymbolicExpr::multiply(a, SymbolicExpr::multiply(y, SymbolicExpr::number(-1))), SymbolicExpr::number(-2)));

        auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});

        EXPECT_TRUE(
            solutions.empty(),
            "unproved symbolic pivot does not fabricate a conditional solution");
    }

    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(2)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto checked_eigenvectors = LMCAS::matrix_eigenvectors_checked(m);
        EXPECT_TRUE(checked_eigenvectors.has_value(),
                    "checked eigenvector solve succeeds");
        if (checked_eigenvectors) {
            const auto& eigs = checked_eigenvectors.value();
            EXPECT_TRUE(eigs.size() == 2, "Should return 2 eigenvectors");
            for (const auto& vector : eigs) {
                EXPECT_TRUE(vector.size() == 2,
                            "each eigenvector has two components");
            }
        }
    }

    {

        auto x = SymbolicExpr::variable("x");
        auto left = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(-6));

        auto solution = LMCAS::InequalitySolver::solve_inequality_checked(
            left, LMCAS::InequalityType::GreaterThan, "x");
        EXPECT_TRUE(
            solution && !solution.value().intervals().empty(),
            "checked inequality solver returns a nonempty interval");
    }

    return TEST_REPORT();
}
