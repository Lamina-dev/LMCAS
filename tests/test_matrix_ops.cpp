#include "test_common.hpp"
#include "symbolic_matrix.hpp"
#include "assumption_context.hpp"

int main() {
    TEST_CASE("Matrix Transpose & Inverse & RREF");

    {
         std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m1_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}
        };
        auto m1 = SymbolicExpr::matrix(m1_data);
        auto t1 = SymbolicExpr::transpose(m1);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m1_t_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4)}
        };
        auto m1_t_expected = SymbolicExpr::matrix(m1_t_data);

        EXPECT_EQ_EXPR(t1, m1_t_expected, "Transpose 2x2");
    }

    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m2_data = {
            {SymbolicExpr::number(4), SymbolicExpr::number(7)},
            {SymbolicExpr::number(2), SymbolicExpr::number(6)}
        };
        auto m2 = SymbolicExpr::matrix(m2_data);
        auto inv2 = lamina::matrix_inverse_checked(m2).value();

        auto prod = SymbolicExpr::multiply(m2, inv2)->simplify();

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto id = SymbolicExpr::matrix(id_data);
        EXPECT_EQ_EXPR(prod, id, "Inverse * Original == Identity");
    }

    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto id = SymbolicExpr::matrix(id_data);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m3_data = {
            {SymbolicExpr::number(2), SymbolicExpr::number(4)},
            {SymbolicExpr::number(1), SymbolicExpr::number(5)}
        };
        auto m3 = SymbolicExpr::matrix(m3_data);
        auto rref3 = lamina::matrix_rref_checked(m3).value();
        EXPECT_EQ_EXPR(rref3, id, "RREF of 2x2 full rank matrix");
    }

    TEST_CASE("Exact matrix kernel proves rank, inverse, nullspace, and solve");
    auto rational_matrix = SymbolicExpr::matrix({
        {SymbolicExpr::number(Rational(1, 2)), SymbolicExpr::number(1)},
        {SymbolicExpr::number(1), SymbolicExpr::number(3)}
    });
    auto determinant = lamina::matrix_determinant_checked(rational_matrix);
    EXPECT_TRUE(determinant && determinant.value()->to_string() == "1/2",
                "Bareiss rational determinant is exact");
    auto inverse = lamina::matrix_inverse_checked(rational_matrix);
    EXPECT_TRUE(inverse.has_value(), "exact rational inverse succeeds");
    if (inverse) {
        auto reconstruction = lamina::matrix_multiply_checked(
            rational_matrix, inverse.value());
        auto identity = SymbolicExpr::matrix({
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        });
        auto simplified_reconstruction = reconstruction
            ? reconstruction.value()->simplify() : nullptr;
        EXPECT_EQ_EXPR(simplified_reconstruction, identity,
                       "A times exact inverse reconstructs identity");
    }
    auto algebraic_matrix = SymbolicExpr::matrix({
        {SymbolicExpr::sqrt(SymbolicExpr::number(2)), SymbolicExpr::number(0)},
        {SymbolicExpr::number(0), SymbolicExpr::sqrt(SymbolicExpr::number(3))}
    });
    auto algebraic_rank = lamina::matrix_rank_checked(algebraic_matrix);
    EXPECT_TRUE(algebraic_rank && algebraic_rank.value() == 2,
                "exact algebraic pivots prove full rank");


    auto rectangular = SymbolicExpr::matrix({
        {SymbolicExpr::number(1), SymbolicExpr::number(2)},
        {SymbolicExpr::number(2), SymbolicExpr::number(4)}
    });
    auto rank = lamina::matrix_rank_checked(rectangular);
    EXPECT_TRUE(rank && rank.value() == 1,
                "rectangular rank-deficient matrix has exact rank one");
    auto nullspace = lamina::matrix_nullspace_checked(rectangular);
    EXPECT_TRUE(nullspace && nullspace.value().size() == 1,
                "rank-deficient matrix has one nullspace basis vector");
    if (nullspace && !nullspace.value().empty()) {
        EXPECT_TRUE(nullspace.value()[0][0]->to_string() == "-2" &&
                        nullspace.value()[0][1]->to_string() == "1",
                    "nullspace basis is exact");
    }

    auto rhs = SymbolicExpr::matrix({
        {SymbolicExpr::number(3)},
        {SymbolicExpr::number(6)}
    });
    auto parametric = lamina::matrix_solve_linear_checked(rectangular, rhs);
    EXPECT_TRUE(parametric &&
                    std::holds_alternative<
                        lamina::MatrixParametricLinearSolution>(
                        parametric.value()),
                "rank-deficient consistent system is parametric");

    auto inconsistent_rhs = SymbolicExpr::matrix({
        {SymbolicExpr::number(3)},
        {SymbolicExpr::number(7)}
    });
    auto inconsistent = lamina::matrix_solve_linear_checked(
        rectangular, inconsistent_rhs);
    EXPECT_TRUE(inconsistent &&
                    std::holds_alternative<
                        lamina::MatrixInconsistentLinearSolution>(
                        inconsistent.value()),
                "proved inconsistent system has typed outcome");

    TEST_CASE("Exact pivots require proof or assumptions");
    auto symbolic = SymbolicExpr::matrix({
        {SymbolicExpr::variable("a")}
    });
    lamina::ComputationContext unknown_context;
    auto unknown_rank = lamina::matrix_rank_checked(
        symbolic, unknown_context);
    EXPECT_TRUE(!unknown_rank &&
                    unknown_rank.error().code == lamina::CasErrc::Inconclusive,
                "unproved symbolic pivot is Inconclusive");

    auto assumptions = std::make_shared<lamina::AssumptionContext>();
    auto assumption = assumptions->assume_sign(
        "a", lamina::Sign::NonZero);
    EXPECT_TRUE(assumption.has_value(), "nonzero pivot assumption is accepted");
    lamina::ComputationContext assumed_context;
    auto installed = assumed_context.set_assumptions(assumptions);
    EXPECT_TRUE(installed.has_value(), "pivot assumptions are installed");
    auto assumed_rank = lamina::matrix_rank_checked(
        symbolic, assumed_context);
    EXPECT_TRUE(assumed_rank && assumed_rank.value() == 1,
                "nonzero assumption proves symbolic pivot");

    auto approximate = SymbolicExpr::matrix({
        {SymbolicExpr::number(1.0)}
    });
    auto approximate_rank = lamina::matrix_rank_checked(approximate);
    EXPECT_TRUE(!approximate_rank &&
                    approximate_rank.error().code ==
                        lamina::CasErrc::Inconclusive,
                "approximate pivot is never exact proof");

    lamina::CancellationToken cancellation;
    cancellation.cancel();
    lamina::ComputationContext cancelled_context({}, cancellation);
    auto cancelled_rank = lamina::matrix_rank_checked(
        rational_matrix, cancelled_context);
    EXPECT_TRUE(!cancelled_rank &&
                    cancelled_rank.error().code == lamina::CasErrc::Cancelled,
                "exact elimination preserves cancellation");
    return TEST_REPORT();
}
