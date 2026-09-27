#include "test_common.hpp"
#include "symbolic_matrix.hpp"
#include "residual_verification.hpp"

#include <iostream>
#include <memory>
#include <variant>
#include <vector>

using namespace LMCAS;

TEST(SymbolicMatrixAdvanced, MatrixEigenvalues) {
    auto matrix = SymbolicExpr::matrix({
        {SymbolicExpr::number(2), SymbolicExpr::number(1)},
        {SymbolicExpr::number(0), SymbolicExpr::number(3)}});
    auto result = matrix_eigenvalues_checked(matrix);

    ASSERT_TRUE(result);
    const auto &eigenvalues = result.value();
    ASSERT_EQ(eigenvalues.size(), 2U);
    bool found_two = false;
    bool found_three = false;
    for (const auto &eigenvalue : eigenvalues) {
        found_two = found_two || test_proved_equivalent(
                                     eigenvalue, SymbolicExpr::number(2));
        found_three = found_three || test_proved_equivalent(
                                         eigenvalue, SymbolicExpr::number(3));
    }
    EXPECT_TRUE(found_two);
    EXPECT_TRUE(found_three);
}

static void expect_eigenvector_reconstruction(const MatrixNode &image,
                                              const std::vector<ExprPtr> &vector,
                                              const std::vector<int> &eigenvalues,
                                              std::size_t dimension) {
    bool matched = false;
    for (int lambda : eigenvalues) {
        bool zero = true;
        ComputationContext context;
        for (std::size_t row = 0; row < dimension; ++row) {
            auto residual = SymbolicExpr::add(
                detail::make_expression_ptr(image.get(row, 0)),
                SymbolicExpr::multiply(SymbolicExpr::number(-lambda), vector[row]));
            auto proof = check_zero_residual(residual, context);
            if (!proof) {
                std::cerr << "eigenvector residual: code=" << static_cast<int>(proof.error().code)
                          << " operation=" << proof.error().operation
                          << " message=" << proof.error().message << '\n';
            }
            zero = zero && proof && std::holds_alternative<ProvedZeroResidual>(proof.value());
        }
        matched = matched || zero;
    }
    EXPECT_TRUE((matched)) << "ordinary vector satisfies Av=lambda*v for an actual eigenvalue";
}

static bool check_eigenvector_image(const ExprPtr &matrix,
                                    const std::vector<std::vector<ExprPtr>> &column,
                                    const std::vector<ExprPtr> &vector,
                                    const std::vector<int> &eigenvalues,
                                    std::size_t dimension) {
    auto v = SymbolicExpr::matrix(column);
    auto av = matrix_multiply_checked(matrix, v);
    EXPECT_TRUE((av.has_value())) << "eigenvector image is defined";
    if (!av) {
        return false;
    }
    auto image = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(av.value()->simplify()));
    EXPECT_TRUE((image && image->rows() == dimension && image->cols() == 1)) << "eigenvector image has the required dimensions";
    if (!image || image->rows() != dimension || image->cols() != 1) {
        return false;
    }
    expect_eigenvector_reconstruction(*image, vector, eigenvalues, dimension);
    return true;
}

static void require_ordinary_eigenbasis(const ExprPtr &matrix, std::size_t dimension,
                                        std::size_t expected_rank,
                                        const std::vector<int> &eigenvalues) {
    auto result = matrix_eigenvectors_checked(matrix);
    if (!result) {
        std::cerr << "ordinary eigenspaces: code=" << static_cast<int>(result.error().code)
                  << " operation=" << result.error().operation
                  << " message=" << result.error().message << '\n';
    }
    ASSERT_TRUE((result.has_value())) << "ordinary eigenspaces are complete";
    if (!result) {
        return;
    }
    const auto &vectors = result.value();
    EXPECT_TRUE((vectors.size() == expected_rank)) << "each eigenspace is returned once";
    if (vectors.size() != expected_rank) {
        return;
    }
    std::vector<std::vector<ExprPtr>> columns(dimension, std::vector<ExprPtr>(vectors.size()));
    for (std::size_t index = 0; index < vectors.size(); ++index) {
        EXPECT_TRUE((vectors[index].size() == dimension)) << "eigenvector dimension matches input";
        if (vectors[index].size() != dimension) {
            return;
        }
        std::vector<std::vector<ExprPtr>> column(dimension, std::vector<ExprPtr>(1));
        for (std::size_t row = 0; row < dimension; ++row) {
            column[row][0] = columns[row][index] = vectors[index][row];
        }
        if (!check_eigenvector_image(matrix, column, vectors[index], eigenvalues, dimension)) {
            return;
        }
    }
    auto rank = matrix_rank_checked(SymbolicExpr::matrix(columns));
    EXPECT_TRUE((rank && rank.value() == expected_rank)) << "ordinary eigenvectors are independent";
}

TEST(SymbolicMatrixAdvanced, MatrixEigenvectors) {
    auto n = [](int value) { return SymbolicExpr::number(value); };
    require_ordinary_eigenbasis(SymbolicExpr::matrix({{n(2), n(1)}, {n(0), n(3)}}), 2, 2, {2, 3});
    require_ordinary_eigenbasis(SymbolicExpr::matrix({{n(0), n(0)}, {n(0), n(0)}}), 2, 2, {0});
    require_ordinary_eigenbasis(SymbolicExpr::matrix({{n(1), n(1), n(0)}, {n(0), n(1), n(0)}, {n(0), n(0), n(2)}}), 3, 2, {1, 2});
    require_ordinary_eigenbasis(SymbolicExpr::matrix({{n(1), n(0), n(1)}, {n(0), n(1), n(0)}, {n(0), n(0), n(2)}}), 3, 3, {1, 2});
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);
    auto exhausted = matrix_eigenvectors_checked(SymbolicExpr::matrix({{n(1)}}), context);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "eigenspace errors are not skipped";
}

TEST(SymbolicMatrixAdvanced, MatrixCheckedContracts) {
    auto x = SymbolicExpr::variable("x");
    auto not_matrix = LMCAS::matrix_determinant_checked(x);
    EXPECT_TRUE((!not_matrix &&
                 not_matrix.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked determinant rejects non-matrix input";

    auto null_mul = LMCAS::matrix_multiply_checked(nullptr, nullptr);
    EXPECT_TRUE((!null_mul &&
                 null_mul.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked matrix multiply rejects null inputs";

    auto A = SymbolicExpr::matrix({{SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
                                   {SymbolicExpr::number(4), SymbolicExpr::number(5), SymbolicExpr::number(6)}});
    auto B = SymbolicExpr::matrix({{SymbolicExpr::number(1), SymbolicExpr::number(2)},
                                   {SymbolicExpr::number(3), SymbolicExpr::number(4)}});
    auto bad_mul = LMCAS::matrix_multiply_checked(A, B);
    EXPECT_TRUE((!bad_mul &&
                 bad_mul.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked matrix multiply rejects incompatible dimensions";

    auto non_square_det = LMCAS::matrix_determinant_checked(A);
    EXPECT_TRUE((!non_square_det &&
                 non_square_det.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked determinant rejects non-square matrix";

    auto singular = SymbolicExpr::matrix({{SymbolicExpr::number(1), SymbolicExpr::number(2)},
                                          {SymbolicExpr::number(2), SymbolicExpr::number(4)}});
    auto singular_inverse = LMCAS::matrix_inverse_checked(singular);
    EXPECT_TRUE((!singular_inverse &&
                 singular_inverse.error().code == LMCAS::CasErrc::DomainError))
        << "checked inverse rejects provably singular matrices";

    auto symbolic_det = SymbolicExpr::matrix({{SymbolicExpr::variable("x"), SymbolicExpr::number(0)},
                                              {SymbolicExpr::number(0), SymbolicExpr::number(1)}});
    auto conditional_inverse = LMCAS::matrix_inverse_checked(symbolic_det);
    EXPECT_TRUE((!conditional_inverse &&
                 conditional_inverse.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked inverse reports inconclusive when determinant nonzero cannot be verified";

    auto unsupported_rot = LMCAS::matrix_rotation_checked(0.0, 3);
    EXPECT_TRUE((!unsupported_rot &&
                 unsupported_rot.error().code ==
                     LMCAS::CasErrc::DimensionMismatch))
        << "checked rotation reports a dimension mismatch";

    auto bad_eigenvalues = LMCAS::matrix_eigenvalues_checked(A);
    EXPECT_TRUE((!bad_eigenvalues &&
                 bad_eigenvalues.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked eigenvalues reject non-square matrix";

    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext context({}, token);
    auto cancelled = LMCAS::matrix_scaling_checked(1.0, 1.0, 2, context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked scaling observes cancelled context";
}
