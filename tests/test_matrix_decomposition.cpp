#include "matrix_decomposition.hpp"
#include "matrix_quadratic_form.hpp"
#include "symbolic_matrix.hpp"
#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"
#include <string>
#include <variant>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> bigint_num(const BigInt &n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(
            std::variant<BigInt, Rational, lmmc_real_t>{std::in_place_type<BigInt>, n}));
}

static std::shared_ptr<SymbolicExpr> exact_dense_matrix(
    size_t rows,
    size_t cols,
    const std::vector<std::shared_ptr<SymbolicExpr>> &entries) {
    MatrixNode::DenseStorage storage;
    storage.reserve(entries.size());
    for (const auto &entry : entries) {
        storage.push_back(LMCAS::detail::node(entry));
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MatrixNode>(rows, cols, std::move(storage)));
}

static std::shared_ptr<SymbolicExpr> mat2(int a, int b, int c, int d) {
    return SymbolicExpr::matrix({{num(a), num(b)}, {num(c), num(d)}});
}

static std::shared_ptr<SymbolicExpr> mat3(int a, int b, int c, int d, int e, int f, int g, int h, int i) {
    return SymbolicExpr::matrix({{num(a), num(b), num(c)}, {num(d), num(e), num(f)}, {num(g), num(h), num(i)}});
}

static bool require_jordan_result(const JordanDecompositionResult &result,
                                  const std::string &message) {
    if (!result)
        std::cerr << message << ": code=" << static_cast<int>(result.error().code)
                  << " operation=" << result.error().operation
                  << " message=" << result.error().message << '\n';
    EXPECT_TRUE((result.has_value())) << (message);
    return result.has_value();
}

static void expect_exact_matrix_entries(const MatrixNode &actual, const MatrixNode &expected,
                                        const std::string &message) {
    ComputationContext context;
    for (std::size_t row = 0; row < actual.rows(); ++row) {
        for (std::size_t col = 0; col < actual.cols(); ++col) {
            auto proof = check_equivalent(detail::make_expression_ptr(actual.get(row, col)),
                                          detail::make_expression_ptr(expected.get(row, col)), context);
            if (!proof) {
                std::cerr << message << ": code=" << static_cast<int>(proof.error().code)
                          << " operation=" << proof.error().operation
                          << " message=" << proof.error().message << '\n';
            }
            EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << (message + " exact entry (" + std::to_string(row) + "," +
                                                                                                  std::to_string(col) + ")");
        }
    }
}

static void expect_exact_matrix(const ExprPtr &actual, const ExprPtr &expected,
                                const std::string &message) {
    auto a = actual ? std::dynamic_pointer_cast<const MatrixNode>(
                          detail::node(actual->simplify()))
                    : nullptr;
    auto b = expected ? std::dynamic_pointer_cast<const MatrixNode>(
                            detail::node(expected->simplify()))
                      : nullptr;
    if (!a || !b) {
        ADD_FAILURE() << (message + " dimensions agree");
        return;
    }
    ASSERT_EQ(a->rows(), b->rows()) << message + " row count";
    ASSERT_EQ(a->cols(), b->cols()) << message + " column count";
    expect_exact_matrix_entries(*a, *b, message);
}

TEST(LmcasMatrixDecomposition, Trace) {
    auto A = mat3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    auto tr = matrix_trace(A);
    EXPECT_TRUE((tr != nullptr)) << ("trace not null");
    EXPECT_TRUE(test_expression_text(tr, num(15))) << ("trace([[1..9]]) = 1+5+9 = 15");
}

TEST(LmcasMatrixDecomposition, Rank) {
    auto I = mat2(1, 0, 0, 1);
    auto identity_rank = matrix_rank_checked(I);
    EXPECT_TRUE((identity_rank && identity_rank.value() == 2)) << ("rank(I2) = 2");
    auto singular = mat2(1, 2, 2, 4);
    auto singular_rank = matrix_rank_checked(singular);
    EXPECT_TRUE((singular_rank && singular_rank.value() == 1)) << ("rank([[1,2],[2,4]]) = 1");
    auto zero = mat2(0, 0, 0, 0);
    auto zero_rank = matrix_rank_checked(zero);
    EXPECT_TRUE((zero_rank && zero_rank.value() == 0)) << ("rank(0) = 0");
}

TEST(LmcasMatrixDecomposition, LuRoundtrip) {
    auto A = mat2(4, 3, 6, 3);
    auto decomposition = lu_decomposition_checked(A);
    ASSERT_TRUE(decomposition.has_value()) << ("checked LU succeeds");
    auto prod = SymbolicExpr::multiply(
                    decomposition.value().L, decomposition.value().U)
                    ->simplify();
    EXPECT_TRUE(test_expression_text(prod, A->simplify())) << ("L*U == A");
}

TEST(LmcasMatrixDecomposition, GramSchmidtOrthogonality) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> vs = {
        {num(1), num(1), num(0)},
        {num(1), num(0), num(1)}};
    auto basis = gram_schmidt(vs, false);
    ASSERT_EQ(basis.size(), 2u) << ("gram_schmidt returns 2 vectors");
    // dot(b0, b1) == 0
    auto dot = SymbolicExpr::number(0);
    for (size_t k = 0; k < 3; ++k)
        dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(basis[0][k], basis[1][k]));
    {
        const auto actual_expr = dot->simplify();
        const auto expected_expr = num(0);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("gram_schmidt vectors orthogonal");
    }
}

TEST(LmcasMatrixDecomposition, GramSchmidtDependence) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> vs = {
        {num(1), num(0)},
        {num(2), num(0)} // dependent
    };
    auto basis = gram_schmidt(vs, false);
    EXPECT_TRUE((basis.size() == 1)) << ("gram_schmidt drops linearly dependent vector");
}

TEST(LmcasMatrixDecomposition, Kronecker) {
    auto A = mat2(1, 2, 3, 4);
    auto B = mat2(0, 1, 1, 0);
    auto K = kronecker(A, B);
    auto kn = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(K));
    EXPECT_TRUE((kn && kn->rows() == 4 && kn->cols() == 4)) << ("kron(2x2,2x2) is 4x4");
    // top-left block = 1*B => [[0,1],[1,0]]; element (0,1) = 1
    {
        const auto actual_expr = LMCAS::detail::make_expression_ptr(kn->get(0, 1));
        const auto expected_expr = num(1);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("kron element (0,1)=1");
    }
    {
        const auto actual_expr = LMCAS::detail::make_expression_ptr(kn->get(0, 0));
        const auto expected_expr = num(0);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("kron element (0,0)=0");
    }
}

TEST(LmcasMatrixDecomposition, FrobeniusNorm) {
    auto A = mat2(3, 0, 0, 4);
    auto fn = matrix_norm(A, "frobenius");
    {
        const auto actual_expr = fn->simplify();
        const auto expected_expr = num(5);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("frobenius([[3,0],[0,4]]) = 5");
    }
}

TEST(LmcasMatrixDecomposition, InducedNorms) {
    auto A = mat2(1, 2, 3, 4);
    auto n1 = matrix_norm(A, "1");
    {
        const auto actual_expr = n1->simplify();
        const auto expected_expr = num(6);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("1-norm = max col sum = 2+4 = 6");
    }
    auto ninf = matrix_norm(A, "inf");
    {
        const auto actual_expr = ninf->simplify();
        const auto expected_expr = num(7);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("inf-norm = max row sum = 3+4 = 7");
    }
}

TEST(LmcasMatrixDecomposition, ExactNormOrdering) {
    const BigInt two_to_53("9007199254740992");
    const BigInt next_integer = two_to_53 + BigInt(1);
    auto A = exact_dense_matrix(2, 2, {bigint_num(two_to_53), bigint_num(next_integer), num(0), num(0)});
    auto n1 = matrix_norm(A, "1");
    EXPECT_TRUE(test_expression_text(n1, bigint_num(next_integer))) << ("1-norm preserves exact large ordering");

    auto B = exact_dense_matrix(2, 2, {bigint_num(two_to_53), num(0), bigint_num(next_integer), num(0)});
    auto ninf = matrix_norm(B, "inf");
    EXPECT_TRUE(test_expression_text(ninf, bigint_num(next_integer))) << ("inf-norm preserves exact large ordering");
}

TEST(LmcasMatrixDecomposition, MatrixExponential) {
    auto Z = mat2(0, 0, 0, 0);
    auto E = matrix_exp(Z);
    auto en = E
        ? std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(E))
        : nullptr;
    ASSERT_TRUE(en && en->rows() == 2 && en->cols() == 2)
        << "matrix_exp of zero must return an explicit 2x2 matrix";
    expect_exact_matrix(
        E, mat2(1, 0, 0, 1), "exp(0) is the complete identity matrix");
}

TEST(LmcasMatrixDecomposition, QuadraticFormClassification) {
    // x^2 + y^2 -> positive definite
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto q = SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y));
    auto A = quadratic_form_matrix(q, {"x", "y"});
    EXPECT_TRUE((A != nullptr)) << ("quadratic_form_matrix not null");
    EXPECT_TRUE((classify_quadratic_form(A) == "positive_definite")) << ("x^2+y^2 positive definite");

    // x^2 - y^2 -> indefinite
    auto q2 = SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                                SymbolicExpr::multiply(num(-1), SymbolicExpr::multiply(y, y)));
    auto A2 = quadratic_form_matrix(q2, {"x", "y"});
    EXPECT_TRUE((classify_quadratic_form(A2) == "indefinite")) << ("x^2-y^2 indefinite");
}

TEST(LmcasMatrixDecomposition, QuadraticFormMatrixValidatesHomogeneousQuadratics) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto a = SymbolicExpr::variable("a");
    auto xy = SymbolicExpr::multiply(x, y);
    auto quadratic = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(num(2), xy),
            SymbolicExpr::multiply(
                num(3), SymbolicExpr::power(y, num(2)))));

    auto matrix_expression =
        quadratic_form_matrix(quadratic, {"x", "y"});
    ASSERT_NE(matrix_expression, nullptr);
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(matrix_expression));
    ASSERT_NE(matrix, nullptr);
    ASSERT_EQ(matrix->rows(), 2u);
    ASSERT_EQ(matrix->cols(), 2u);
    const std::shared_ptr<SymbolicExpr> expected[2][2] = {
        {num(1), num(1)}, {num(1), num(3)}};
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 2; ++col) {
            EXPECT_TRUE(test_proved_equivalent(
                detail::make_expression_ptr(matrix->get(row, col)),
                expected[row][col]));
        }
    }

    auto parameterized = quadratic_form_matrix(
        SymbolicExpr::multiply(
            a, SymbolicExpr::power(x, num(2))),
        {"x"});
    ASSERT_NE(parameterized, nullptr);
    auto parameter_matrix = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(parameterized));
    ASSERT_NE(parameter_matrix, nullptr);
    ASSERT_EQ(parameter_matrix->rows(), 1u);
    ASSERT_EQ(parameter_matrix->cols(), 1u);
    EXPECT_TRUE(test_proved_equivalent(
        detail::make_expression_ptr(parameter_matrix->get(0, 0)), a));

    EXPECT_EQ(quadratic_form_matrix(nullptr, {"x"}), nullptr);
    EXPECT_EQ(quadratic_form_matrix(quadratic, {}), nullptr);
    EXPECT_EQ(quadratic_form_matrix(quadratic, {"x", "x"}), nullptr);
    EXPECT_EQ(
        quadratic_form_matrix(
            SymbolicExpr::power(x, num(3)), {"x"}),
        nullptr);
    EXPECT_EQ(
        quadratic_form_matrix(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(2)), x),
            {"x"}),
        nullptr);
}

TEST(LmcasMatrixDecomposition, DiagonalJordan) {
    auto A = mat2(2, 0, 0, 3);
    auto decomposition = jordan_form_checked(A);
    if (require_jordan_result(decomposition,
                              "checked Jordan form exists for a diagonal matrix")) {
        auto inverse =
            matrix_inverse_checked(decomposition.value().P).value();
        auto recon = SymbolicExpr::multiply(
                         decomposition.value().P,
                         SymbolicExpr::multiply(
                             decomposition.value().J, inverse))
                         ->simplify();
        expect_exact_matrix(recon, A, "P J P^-1 == A");
    }
}

TEST(LmcasMatrixDecomposition, CheckedLu) {
    auto A = mat2(4, 3, 6, 3);
    auto lu = lu_decomposition_checked(A);
    ASSERT_TRUE(lu.has_value()) << "checked LU succeeds";
    auto permuted = SymbolicExpr::multiply(lu.value().P, A)->simplify();
    auto prod =
        SymbolicExpr::multiply(lu.value().L, lu.value().U)->simplify();
    EXPECT_TRUE(test_expression_text(prod, permuted))
        << "checked PLU reconstructs P*A";
}

TEST(LmcasMatrixDecomposition, CheckedQrCholesky) {
    auto qr = qr_decomposition_checked(mat2(1, 0, 0, 1));
    ASSERT_TRUE(qr.has_value()) << "checked QR succeeds on identity";
    EXPECT_TRUE(qr.value().Q != nullptr && qr.value().R != nullptr)
        << "checked QR returns Q and R";

    auto chol = cholesky_decomposition_checked(mat2(4, 0, 0, 9));
    ASSERT_TRUE(chol.has_value())
        << "checked Cholesky succeeds on diagonal SPD matrix";
    EXPECT_TRUE(chol.value().L != nullptr)
        << "checked Cholesky returns L";
}

TEST(LmcasMatrixDecomposition, CheckedJordan) {
    auto jordan = jordan_form_checked(mat2(2, 0, 0, 3));
    if (require_jordan_result(
            jordan, "checked Jordan succeeds on diagonal matrix")) {
        ASSERT_TRUE(
            jordan.value().J != nullptr && jordan.value().P != nullptr)
            << "checked Jordan returns J and P";
        auto inverse = matrix_inverse_checked(jordan.value().P);
        ASSERT_TRUE(inverse.has_value())
            << "checked Jordan returns invertible P";
        ASSERT_TRUE(inverse.value() != nullptr);
        auto reconstructed =
            SymbolicExpr::multiply(
                jordan.value().P,
                SymbolicExpr::multiply(jordan.value().J, inverse.value()))
                ->simplify();
        expect_exact_matrix(
            reconstructed, mat2(2, 0, 0, 3),
            "checked Jordan reconstructs exact diagonal input");
    }
}

TEST(LmcasMatrixDecomposition, CheckedJordanBlock) {
    auto jordan_block = mat2(2, 1, 0, 2);
    auto block_form = jordan_form_checked(jordan_block);
    if (require_jordan_result(block_form,
                              "checked Jordan accepts a non-diagonal Jordan block")) {
        auto Pinv = matrix_inverse_checked(block_form.value().P).value();
        auto reconstructed = Pinv
                                 ? SymbolicExpr::multiply(
                                       block_form.value().P,
                                       SymbolicExpr::multiply(block_form.value().J, Pinv))
                                       ->simplify()
                                 : nullptr;
        expect_exact_matrix(reconstructed, jordan_block,
                            "Jordan block satisfies A=P*J*P^-1");
    }
}

TEST(LmcasMatrixDecomposition, CheckedSvd) {
    auto svd = svd_decomposition_checked(mat2(2, 0, 0, 3));
    ASSERT_TRUE((svd.has_value())) << ("checked SVD succeeds on exact nonnegative diagonal matrix");
    if (svd) {
        auto reconstructed = SymbolicExpr::multiply(
                                 svd.value().U,
                                 SymbolicExpr::multiply(svd.value().S, SymbolicExpr::transpose(svd.value().V)))
                                 ->simplify();
        EXPECT_TRUE(test_expression_text(reconstructed, mat2(2, 0, 0, 3)->simplify())) << ("checked SVD reconstructs exact diagonal input");
    }
}

TEST(LmcasMatrixDecomposition, LuDomain) {
    auto non_square = SymbolicExpr::matrix({{num(1), num(2)}});
    auto bad_lu = lu_decomposition_checked(non_square);
    EXPECT_TRUE((!bad_lu && bad_lu.error().code == CasErrc::InvalidArgument)) << ("checked LU rejects non-square matrix");

    auto needs_pivot_lu = lu_decomposition_checked(mat2(0, 1, 1, 0));
    ASSERT_TRUE((needs_pivot_lu.has_value())) << ("checked PLU succeeds when row pivoting is required");
    if (needs_pivot_lu) {
        auto pivoted = SymbolicExpr::multiply(
                           needs_pivot_lu.value().P, mat2(0, 1, 1, 0))
                           ->simplify();
        auto reconstructed = SymbolicExpr::multiply(
                                 needs_pivot_lu.value().L,
                                 needs_pivot_lu.value().U)
                                 ->simplify();
        EXPECT_TRUE(test_expression_text(reconstructed, pivoted)) << ("pivoting PLU satisfies P*A=L*U");
    }
}

TEST(LmcasMatrixDecomposition, DecompositionContextErrors) {
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled_context({}, token);
    auto cancelled = lu_decomposition_checked(mat2(1, 0, 0, 1), cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << ("checked LU observes cancellation");

    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext limited_context(limits);
    auto limited = jordan_form_checked(mat2(2, 0, 0, 3), limited_context);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << ("checked Jordan observes exhausted step budget");
}

TEST(LmcasMatrixDecomposition, SymbolicLu) {
    auto x = SymbolicExpr::variable("x");
    auto symbolic_lu_input = SymbolicExpr::matrix({{x, num(1)}, {num(1), num(1)}});
    auto symbolic_lu = lu_decomposition_checked(symbolic_lu_input);
    EXPECT_TRUE((!symbolic_lu && symbolic_lu.error().code == CasErrc::Inconclusive)) << ("checked PLU requires exact rational entries or proved pivots");
}

TEST(LmcasMatrixDecomposition, QrDomain) {
    auto x = SymbolicExpr::variable("x");
    auto bad_qr = qr_decomposition_checked(num(1));
    EXPECT_TRUE((!bad_qr && bad_qr.error().code == CasErrc::InvalidArgument)) << ("checked QR rejects non-matrix input");

    auto dependent_qr = qr_decomposition_checked(mat2(1, 0, 2, 0));
    EXPECT_TRUE((!dependent_qr && dependent_qr.error().code == CasErrc::Inconclusive)) << ("checked QR reports Inconclusive for rank-deficient columns");

    auto wide_qr = qr_decomposition_checked(
        SymbolicExpr::matrix({{num(1), num(0), num(0)}, {num(0), num(1), num(0)}}));
    EXPECT_TRUE((!wide_qr && wide_qr.error().code == CasErrc::Inconclusive)) << ("checked QR reports Inconclusive outside tall/full-column-rank support");

    auto symbolic_qr_input = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
    auto symbolic_qr = qr_decomposition_checked(symbolic_qr_input);
    EXPECT_TRUE((!symbolic_qr && symbolic_qr.error().code == CasErrc::Inconclusive)) << ("checked QR requires proven exact rational full-column-rank support");
}

TEST(LmcasMatrixDecomposition, CholeskyDomain) {
    auto x = SymbolicExpr::variable("x");
    auto symbolic_spd = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
    auto null_chol = cholesky_decomposition_checked(nullptr);
    EXPECT_TRUE((!null_chol && null_chol.error().code == CasErrc::InvalidArgument)) << ("checked Cholesky rejects null input");

    auto non_spd_chol = cholesky_decomposition_checked(mat2(-1, 0, 0, 1));
    EXPECT_TRUE((!non_spd_chol && non_spd_chol.error().code == CasErrc::DomainError)) << ("checked Cholesky rejects proven non-SPD matrices");

    auto semidefinite_chol = cholesky_decomposition_checked(mat2(1, 0, 0, 0));
    EXPECT_TRUE((!semidefinite_chol && semidefinite_chol.error().code == CasErrc::DomainError)) << ("checked Cholesky rejects positive semidefinite matrices");
    auto symbolic_chol = cholesky_decomposition_checked(symbolic_spd);
    EXPECT_TRUE((!symbolic_chol && symbolic_chol.error().code == CasErrc::Inconclusive)) << ("checked Cholesky requires a proven exact rational SPD matrix");
}

TEST(LmcasMatrixDecomposition, JordanDomain) {
    auto non_square = SymbolicExpr::matrix({{num(1), num(2)}});
    auto x = SymbolicExpr::variable("x");
    auto symbolic_spd = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
    auto bad_jordan = jordan_form_checked(non_square);
    EXPECT_TRUE((!bad_jordan && bad_jordan.error().code == CasErrc::InvalidArgument)) << ("checked Jordan rejects non-square matrix");

    auto non_diagonal_jordan = jordan_form_checked(mat2(1, 1, 0, 1));
    require_jordan_result(non_diagonal_jordan,
                          "checked Jordan supports an exact rational Jordan block");

    auto symbolic_jordan = jordan_form_checked(symbolic_spd);
    EXPECT_TRUE((!symbolic_jordan && symbolic_jordan.error().code == CasErrc::Inconclusive)) << ("checked Jordan requires exact rational entries or proved chains");
}

TEST(LmcasMatrixDecomposition, SvdDomain) {
    auto x = SymbolicExpr::variable("x");
    auto symbolic_spd = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
    auto non_diagonal_svd = svd_decomposition_checked(mat2(1, 1, 0, 1));
    EXPECT_TRUE((!non_diagonal_svd &&
                 non_diagonal_svd.error().code == CasErrc::Inconclusive))
        << ("checked SVD remains explicit when a complete singular basis is unproved");

    auto negative_diagonal_svd = svd_decomposition_checked(mat2(-1, 0, 0, 1));
    EXPECT_TRUE((negative_diagonal_svd.has_value())) << ("checked SVD absorbs diagonal signs into singular vectors");

    auto symbolic_svd = svd_decomposition_checked(symbolic_spd);
    EXPECT_TRUE((!symbolic_svd && symbolic_svd.error().code == CasErrc::Inconclusive)) << ("checked SVD requires exact rational entries or proved eigenspaces");
}

static void require_jordan_certificate(const ExprPtr &input, const JordanDecomposition &value,
                                       std::size_t dimension) {
    auto rank = matrix_rank_checked(value.P);
    ASSERT_TRUE(rank.has_value());
    EXPECT_EQ(rank.value(), dimension) << "Jordan basis has full rank";
    auto ap = matrix_multiply_checked(input, value.P);
    auto pj = matrix_multiply_checked(value.P, value.J);
    ASSERT_TRUE(ap.has_value() && pj.has_value())
        << "Jordan certificate products are defined";
    expect_exact_matrix(
        ap.value(), pj.value(), "Jordan certificate AP=PJ");
}

TEST(LmcasMatrixDecomposition, RepeatedEigenspaces) {
    for (const auto &item : std::vector<std::pair<ExprPtr, std::size_t>>{
             {mat3(1, 0, 1, 0, 1, 0, 0, 0, 2), 3}, {mat3(1, 1, 0, 0, 2, 0, 0, 1, 1), 3}, {mat2(0, 0, 0, 0), 2}, {mat2(1, 0, 0, 1), 2}, {mat2(2, 1, 0, 3), 2}}) {
        const auto &input = item.first;
        auto result = jordan_form_checked(input);
        if (require_jordan_result(result, "complete repeated eigenspaces diagonalize"))
            require_jordan_certificate(input, result.value(), item.second);
    }
    auto defective = mat3(1, 1, 0, 0, 1, 0, 0, 0, 2);
    auto result = jordan_form_checked(defective);
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << ("larger defective matrix does not return a singular basis");
    auto block3 = jordan_form_checked(mat3(1, 1, 0, 0, 1, 1, 0, 0, 1));
    EXPECT_TRUE((!block3 && block3.error().code == CasErrc::Inconclusive)) << ("unsupported length-three chain is explicit");
    auto x = SymbolicExpr::variable("x");
    auto symbolic = jordan_form_checked(SymbolicExpr::matrix({{num(0), num(0)}, {num(0), x}}));
    EXPECT_TRUE((!symbolic && symbolic.error().code == CasErrc::Inconclusive)) << ("unproved eigenvalue distinction is not a complete basis");
    for (const auto &scale : {num(1), num(2), SymbolicExpr::divide(num(-3), num(2))}) {
        auto input = SymbolicExpr::matrix({{num(2), scale}, {num(0), num(2)}});
        auto normalized = jordan_form_checked(input);
        if (require_jordan_result(normalized, "rational two-dimensional chain is supported")) {
            require_jordan_certificate(input, normalized.value(), 2);
            expect_exact_matrix(normalized.value().J, mat2(2, 1, 0, 2), "Jordan superdiagonal is one");
        }
    }
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = jordan_form_checked(defective, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << ("Jordan cancellation remains Cancelled");
}

TEST(LmcasMatrixDecomposition, PairedMatrixFunctions) {
    auto input = mat3(1, 0, 1, 0, 1, 0, 0, 0, 2);
    auto e = SymbolicExpr::exp(num(1));
    auto e2 = SymbolicExpr::exp(num(2));
    auto expected_exp = SymbolicExpr::matrix({{e, num(0), SymbolicExpr::add(e2, SymbolicExpr::multiply(num(-1), e))},
                                              {num(0), e, num(0)},
                                              {num(0), num(0), e2}});
    expect_exact_matrix(matrix_exp(input), expected_exp, "exp uses paired repeated eigenspaces");
    auto log2 = SymbolicExpr::ln(num(2));
    auto expected_log = SymbolicExpr::matrix({{num(0), num(0), log2}, {num(0), num(0), num(0)}, {num(0), num(0), log2}});
    expect_exact_matrix(matrix_log(input), expected_log, "log uses paired repeated eigenspaces");
    auto defect = mat3(1, 1, 0, 0, 1, 0, 0, 0, 2);
    {
        const auto actual_expr = matrix_exp(defect);
        const auto expected_expr = SymbolicExpr::exp(defect);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("unsupported defect exponential remains unevaluated");
    }
    EXPECT_TRUE((matrix_log(defect) == nullptr)) << ("unsupported defect logarithm fails");
    EXPECT_TRUE((matrix_log(mat2(-1, 0, 0, 2)) == nullptr)) << ("real logarithm needs positive spectrum");
}

TEST(LmcasMatrixDecomposition, SymbolicNormCandidates) {
    auto radical = SymbolicExpr::sqrt(num(2));
    auto x = SymbolicExpr::variable("x");
    for (const auto &type : {std::string("1"), std::string("inf")}) {
        for (const auto &item : std::vector<std::pair<ExprPtr, double>>{
                 {SymbolicExpr::matrix({{num(1), num(0)}, {num(0), radical}}), std::sqrt(2.0)},
                 {SymbolicExpr::matrix({{radical, num(0)}, {num(0), num(1)}}), std::sqrt(2.0)},
                 {SymbolicExpr::matrix({{num(2), num(0)}, {num(0), radical}}), 2.0},
                 {SymbolicExpr::matrix({{radical, num(1)}, {num(1), num(0)}}), 1.0 + std::sqrt(2.0)}}) {
            auto norm = matrix_norm(item.first, type);
            ASSERT_TRUE(norm != nullptr) << "norm candidate expression exists";
            auto value = evaluate_numeric(*norm);
            ASSERT_TRUE(value.has_value()) << "radical norm evaluates";
            EXPECT_TRUE(std::isfinite(value.value().value))
                << "complete axis sums determine the norm";
            EXPECT_NEAR(value.value().value, item.second, 1e-12)
                << "complete axis sums determine the norm";
        }
        auto symbolic = matrix_norm(SymbolicExpr::matrix({{num(1), num(0)}, {num(0), x}}), type);
        for (const auto point : {0, -2}) {
            auto value = evaluate_numeric(*symbolic->substitute("x", num(point)));
            ASSERT_TRUE((value.has_value())) << ("symbolic norm specializes");
            const auto actual_value = value.value().value;
            const auto expected_value = point == 0 ? 1.0 : 2.0;
            EXPECT_TRUE(std::isfinite(actual_value))
                << "all symbolic norm candidates survive";
            EXPECT_NEAR(actual_value, expected_value, 1e-12)
                << "all symbolic norm candidates survive";
        }
    }
}

TEST(LmcasMatrixDecomposition, ExactQuadraticSigns) {
    for (const auto &item : std::vector<std::pair<ExprPtr, std::string>>{
             {mat2(0, 0, 0, 0), "positive_semidefinite"}, {SymbolicExpr::matrix({{num(0)}}), "positive_semidefinite"}, {mat2(1, 0, 0, 0), "positive_semidefinite"}, {mat2(-1, 0, 0, 0), "negative_semidefinite"}, {mat2(1, 0, 0, 1), "positive_definite"}, {mat2(-1, 0, 0, -1), "negative_definite"}, {mat2(1, 0, 0, -1), "indefinite"}, {SymbolicExpr::matrix({{SymbolicExpr::variable("x"), num(0)}, {num(0), num(1)}}), "unknown"}}) {
        EXPECT_TRUE((classify_quadratic_form(item.first) == item.second)) << ("proved quadratic sign classification");
    }
    auto tiny = SymbolicExpr::number(Rational(BigInt(1), BigInt("1000000000000")));
    EXPECT_TRUE((classify_quadratic_form(SymbolicExpr::matrix({{tiny, num(0)}, {num(0), tiny}})) ==
                 "positive_definite"))
        << ("exact tiny positive spectrum is not zero");
    auto radical = SymbolicExpr::sqrt(num(2));
    EXPECT_TRUE((classify_quadratic_form(SymbolicExpr::matrix({{radical, num(0)}, {num(0), num(1)}})) ==
                 "positive_definite"))
        << ("radical signs are proved without approximate ordering");
    for (double value : {0.0, 1e-12, -1e-12}) {
        auto approximate = SymbolicExpr::number(value);
        auto expected = value > 0 ? "positive_definite" : (value < 0 ? "negative_definite" : "positive_semidefinite");
        EXPECT_TRUE((classify_quadratic_form(SymbolicExpr::matrix({{approximate}})) == expected)) << ("approximate spectrum uses its stored sign, not a zero tolerance");
    }
}
