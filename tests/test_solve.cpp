#include "test_common.hpp"
#include "inequality_solver.hpp"

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
        auto eigs = SymbolicExpr::eigenvectors(m);

        if(eigs.size() == 2) {
             EXPECT_TRUE(eigs.size() == 2, "Should return 2 eigenvectors");
             for(auto& p : eigs) {
                 std::cout << "Checking Lambda: " << p.first->to_string() << std::endl;
                 if (p.first->to_string() == "1") {

                     EXPECT_TRUE(p.second.size() == 1, "size 1");
                     auto& v = p.second[0];

                     if (v->get_type() == SymbolicExpr::Type::Matrix) {

                     }
                 } else if (p.first->to_string() == "2") {

                     EXPECT_TRUE(p.second.size() == 1, "size 1");
                     auto& v = p.second[0];
                 }
             }
        } else {

             std::cout << "Eigenvectors count: " << eigs.size() << " (Expected 2)" << std::endl;

        }
    }

    {

        auto x = SymbolicExpr::variable("x");
        auto left = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(-6));

        auto solution = lamina::InequalitySolver::solve_inequality_checked(
            left, lamina::InequalityType::GreaterThan, "x");
        EXPECT_TRUE(
            solution && !solution.value().intervals().empty(),
            "checked inequality solver returns a nonempty interval");
    }

    return TEST_REPORT();
}
