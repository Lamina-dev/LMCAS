#include "test_common.hpp"
#include "inequality_solver.hpp"
#include "symbolic_matrix.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

namespace {

TEST(Solve, DirectLinearSystem) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto eq1 = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::add(y, SymbolicExpr::number(-5)));

    auto eq2 = SymbolicExpr::add(
        x,
        SymbolicExpr::add(SymbolicExpr::multiply(y, SymbolicExpr::number(-1)), SymbolicExpr::number(-1)));

    auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
    ASSERT_EQ(solutions.size(), 1u) << "Solutions size should be 1";
    EXPECT_TRUE(test_same_expression(solutions[0].at("x"), SymbolicExpr::number(2)))
        << "x should be 2";
    EXPECT_TRUE(test_same_expression(solutions[0].at("y"), SymbolicExpr::number(1)))
        << "y should be 1";
}

TEST(Solve, SymbolicPivotSystem) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto a = SymbolicExpr::variable("a");

    auto eq1 = SymbolicExpr::add(x, SymbolicExpr::multiply(a, y));
    auto eq2 = SymbolicExpr::add(x, SymbolicExpr::add(SymbolicExpr::multiply(a, SymbolicExpr::multiply(y, SymbolicExpr::number(-1))), SymbolicExpr::number(-2)));

    auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});

    EXPECT_TRUE(solutions.empty()) << "unproved symbolic pivot does not fabricate a conditional solution";
}

TEST(Solve, DiagonalEigenvectorSystem) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
        {SymbolicExpr::number(1), SymbolicExpr::number(0)},
        {SymbolicExpr::number(0), SymbolicExpr::number(2)}};
    auto m = SymbolicExpr::matrix(m_data);
    auto checked_eigenvectors = LMCAS::matrix_eigenvectors_checked(m);
    if (!checked_eigenvectors) {
        std::cerr << "matrix_eigenvectors: code=" << static_cast<int>(checked_eigenvectors.error().code)
                  << " operation=" << checked_eigenvectors.error().operation
                  << " message=" << checked_eigenvectors.error().message << '\n';
    }
    EXPECT_TRUE(checked_eigenvectors.has_value()) << "checked eigenvector solve succeeds";
    if (checked_eigenvectors) {
        const auto &eigs = checked_eigenvectors.value();
        EXPECT_TRUE(eigs.size() == 2) << "Should return 2 eigenvectors";
        if (eigs.size() == 2 && eigs[0].size() == 2 && eigs[1].size() == 2) {
            auto basis = SymbolicExpr::matrix({{eigs[0][0], eigs[1][0]}, {eigs[0][1], eigs[1][1]}});
            auto rank = matrix_rank_checked(basis);
            EXPECT_TRUE(rank && rank.value() == 2) << "returned eigenvectors are independent";
            for (const auto &vector : eigs) {
                ComputationContext context;
                auto residual = check_zero_residual(
                    SymbolicExpr::multiply(vector[0], vector[1]), context);
                if (!residual) {
                    std::cerr << "eigenvector residual: code=" << static_cast<int>(residual.error().code)
                              << " operation=" << residual.error().operation
                              << " message=" << residual.error().message << '\n';
                }
                EXPECT_TRUE(residual && std::holds_alternative<ProvedZeroResidual>(residual.value())) << "vector lies in an actual eigenspace";
            }
        } else {
            ADD_FAILURE() << "eigenbasis dimensions match the matrix";
        }
    }
}

TEST(Solve, LinearInequality) {
    auto x = SymbolicExpr::variable("x");
    auto left = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(-6));

    auto solution = LMCAS::InequalitySolver::solve_inequality_checked(
        left, LMCAS::InequalityType::GreaterThan, "x");
    EXPECT_TRUE(solution && !solution.value().intervals().empty()) << "checked inequality solver returns a nonempty interval";
}

} // namespace
