#include "matrix_decomposition.hpp"
#include "symbolic_matrix.hpp"
#include "test_common.hpp"
#include <string>
#include <variant>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> bigint_num(const BigInt& n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(
            std::variant<BigInt, Rational, lmmc_real_t>{std::in_place_type<BigInt>, n}));
}

static std::shared_ptr<SymbolicExpr> exact_dense_matrix(
    size_t rows,
    size_t cols,
    const std::vector<std::shared_ptr<SymbolicExpr>>& entries) {
    MatrixNode::DenseStorage storage;
    storage.reserve(entries.size());
    for (const auto& entry : entries) {
        storage.push_back(LMCAS::detail::node(entry));
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MatrixNode>(rows, cols, std::move(storage)));
}

static std::shared_ptr<SymbolicExpr> mat2(int a, int b, int c, int d) {
    return SymbolicExpr::matrix({{num(a), num(b)}, {num(c), num(d)}});
}

static std::shared_ptr<SymbolicExpr> mat3(int a,int b,int c,int d,int e,int f,int g,int h,int i){
    return SymbolicExpr::matrix({{num(a),num(b),num(c)},{num(d),num(e),num(f)},{num(g),num(h),num(i)}});
}

int main() {
    // ---- trace ----
    {
        auto A = mat3(1,2,3, 4,5,6, 7,8,9);
        auto tr = matrix_trace(A);
        EXPECT_TRUE(tr != nullptr, "trace not null");
        EXPECT_EQ_EXPR(tr, num(15), "trace([[1..9]]) = 1+5+9 = 15");
    }

    // ---- rank ----
    {
        auto I = mat2(1,0, 0,1);
        auto identity_rank = matrix_rank_checked(I);
        EXPECT_TRUE(identity_rank && identity_rank.value() == 2,
                    "rank(I2) = 2");
        auto singular = mat2(1,2, 2,4);
        auto singular_rank = matrix_rank_checked(singular);
        EXPECT_TRUE(singular_rank && singular_rank.value() == 1,
                    "rank([[1,2],[2,4]]) = 1");
        auto zero = mat2(0,0, 0,0);
        auto zero_rank = matrix_rank_checked(zero);
        EXPECT_TRUE(zero_rank && zero_rank.value() == 0, "rank(0) = 0");
    }

    // ---- LU round-trip P*A = L*U (here no pivoting): A = L*U ----
    {
        auto A = mat2(4,3, 6,3);
        auto decomposition = lu_decomposition_checked(A);
        EXPECT_TRUE(decomposition.has_value(), "checked LU succeeds");
        auto prod = SymbolicExpr::multiply(
            decomposition.value().L, decomposition.value().U)->simplify();
        EXPECT_EQ_EXPR(prod, A->simplify(), "L*U == A");
    }

    // ---- Gram-Schmidt orthogonality ----
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> vs = {
            {num(1), num(1), num(0)},
            {num(1), num(0), num(1)}
        };
        auto basis = gram_schmidt(vs, false);
        EXPECT_TRUE(basis.size() == 2, "gram_schmidt returns 2 vectors");
        // dot(b0, b1) == 0
        auto dot = SymbolicExpr::number(0);
        for (size_t k = 0; k < 3; ++k)
            dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(basis[0][k], basis[1][k]));
        EXPECT_EQ_EXPR(dot->simplify(), num(0), "gram_schmidt vectors orthogonal");
    }

    // ---- Gram-Schmidt skips dependent vector ----
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> vs = {
            {num(1), num(0)},
            {num(2), num(0)}  // dependent
        };
        auto basis = gram_schmidt(vs, false);
        EXPECT_TRUE(basis.size() == 1, "gram_schmidt drops linearly dependent vector");
    }

    // ---- Kronecker product dimensions and elements ----
    {
        auto A = mat2(1,2, 3,4);
        auto B = mat2(0,1, 1,0);
        auto K = kronecker(A, B);
        auto kn = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(K));
        EXPECT_TRUE(kn && kn->rows() == 4 && kn->cols() == 4, "kron(2x2,2x2) is 4x4");
        // top-left block = 1*B => [[0,1],[1,0]]; element (0,1) = 1
        EXPECT_EQ_EXPR(LMCAS::detail::make_expression_ptr(kn->get(0,1)), num(1), "kron element (0,1)=1");
        EXPECT_EQ_EXPR(LMCAS::detail::make_expression_ptr(kn->get(0,0)), num(0), "kron element (0,0)=0");
    }

    // ---- Frobenius norm ----
    {
        auto A = mat2(3,0, 0,4);
        auto fn = matrix_norm(A, "frobenius");
        EXPECT_EQ_EXPR(fn->simplify(), num(5), "frobenius([[3,0],[0,4]]) = 5");
    }

    // ---- 1-norm (max col sum) and inf-norm (max row sum) ----
    {
        auto A = mat2(1,2, 3,4);
        auto n1 = matrix_norm(A, "1");
        EXPECT_EQ_EXPR(n1->simplify(), num(6), "1-norm = max col sum = 2+4 = 6");
        auto ninf = matrix_norm(A, "inf");
        EXPECT_EQ_EXPR(ninf->simplify(), num(7), "inf-norm = max row sum = 3+4 = 7");
    }

    // ---- exact large norm comparisons ----
    {
        const BigInt two_to_53("9007199254740992");
        const BigInt next_integer = two_to_53 + BigInt(1);
        auto A = exact_dense_matrix(2, 2, {
            bigint_num(two_to_53), bigint_num(next_integer),
            num(0), num(0)
        });
        auto n1 = matrix_norm(A, "1");
        std::string n1_text = n1 ? n1->to_string() : "<null>";
        EXPECT_TRUE(n1 && n1->to_string().find(next_integer.to_string()) != std::string::npos,
                    "1-norm keeps exact large column-sum ordering, got " + n1_text);

        auto B = exact_dense_matrix(2, 2, {
            bigint_num(two_to_53), num(0),
            bigint_num(next_integer), num(0)
        });
        auto ninf = matrix_norm(B, "inf");
        std::string ninf_text = ninf ? ninf->to_string() : "<null>";
        EXPECT_TRUE(ninf && ninf->to_string().find(next_integer.to_string()) != std::string::npos,
                    "inf-norm keeps exact large row-sum ordering, got " + ninf_text);
    }

    // ---- matrix_exp of zero is identity ----
    {
        auto Z = mat2(0,0, 0,0);
        auto E = matrix_exp(Z);
        auto en = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(E));
        if (en) {
            EXPECT_EQ_EXPR(LMCAS::detail::make_expression_ptr(en->get(0,0))->simplify(), num(1), "exp(0)[0,0]=1");
            EXPECT_EQ_EXPR(LMCAS::detail::make_expression_ptr(en->get(0,1))->simplify(), num(0), "exp(0)[0,1]=0");
        } else {
            EXPECT_TRUE(false, "matrix_exp of zero must return an explicit matrix");
        }
    }

    // ---- quadratic form classification ----
    {
        // x^2 + y^2 -> positive definite
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto q = SymbolicExpr::add(SymbolicExpr::multiply(x, x), SymbolicExpr::multiply(y, y));
        auto A = quadratic_form_matrix(q, {"x", "y"});
        EXPECT_TRUE(A != nullptr, "quadratic_form_matrix not null");
        EXPECT_TRUE(classify_quadratic_form(A) == "positive_definite", "x^2+y^2 positive definite");

        // x^2 - y^2 -> indefinite
        auto q2 = SymbolicExpr::add(SymbolicExpr::multiply(x, x),
            SymbolicExpr::multiply(num(-1), SymbolicExpr::multiply(y, y)));
        auto A2 = quadratic_form_matrix(q2, {"x", "y"});
        EXPECT_TRUE(classify_quadratic_form(A2) == "indefinite", "x^2-y^2 indefinite");
    }

    // ---- Jordan form of diagonalizable matrix: A = P J P^-1 ----
    {
        auto A = mat2(2,0, 0,3);
        auto decomposition = jordan_form_checked(A);
        EXPECT_TRUE(decomposition.has_value(),
                    "checked Jordan form exists for a diagonal matrix");
        if (decomposition) {
            auto inverse =
                matrix_inverse_checked(decomposition.value().P).value();
            auto recon = SymbolicExpr::multiply(
                decomposition.value().P,
                SymbolicExpr::multiply(
                    decomposition.value().J, inverse))->simplify();
            EXPECT_EQ_EXPR(recon, A->simplify(), "P J P^-1 == A");
        }
    }

    // ---- checked decomposition APIs: success and explicit errors ----
    {
        auto A = mat2(4,3, 6,3);
        auto lu = lu_decomposition_checked(A);
        EXPECT_TRUE(lu.has_value(), "checked LU succeeds");
        if (lu) {
            auto permuted = SymbolicExpr::multiply(lu.value().P, A)->simplify();
            auto prod = SymbolicExpr::multiply(
                lu.value().L, lu.value().U)->simplify();
            EXPECT_EQ_EXPR(prod, permuted, "checked PLU reconstructs P*A");
        }

        auto qr = qr_decomposition_checked(mat2(1,0, 0,1));
        EXPECT_TRUE(qr.has_value(), "checked QR succeeds on identity");
        if (qr) {
            EXPECT_TRUE(qr.value().Q != nullptr && qr.value().R != nullptr,
                        "checked QR returns Q and R");
        }

        auto chol = cholesky_decomposition_checked(mat2(4,0, 0,9));
        EXPECT_TRUE(chol.has_value(), "checked Cholesky succeeds on diagonal SPD matrix");
        if (chol) {
            EXPECT_TRUE(chol.value().L != nullptr, "checked Cholesky returns L");
        }

        auto jordan = jordan_form_checked(mat2(2,0, 0,3));
        EXPECT_TRUE(jordan.has_value(), "checked Jordan succeeds on diagonal matrix");
        if (jordan) {
            EXPECT_TRUE(jordan.value().J != nullptr && jordan.value().P != nullptr,
                        "checked Jordan returns J and P");
            auto Pinv = matrix_inverse_checked(jordan.value().P).value();
            EXPECT_TRUE(Pinv != nullptr, "checked Jordan returns invertible P");
            if (Pinv) {
                auto reconstructed = SymbolicExpr::multiply(
                    jordan.value().P,
                    SymbolicExpr::multiply(jordan.value().J, Pinv))->simplify();
                EXPECT_EQ_EXPR(reconstructed, mat2(2,0, 0,3)->simplify(),
                               "checked Jordan reconstructs exact diagonal input");
            }
        }

        auto jordan_block = mat2(2,1, 0,2);
        auto block_form = jordan_form_checked(jordan_block);
        EXPECT_TRUE(block_form.has_value(),
                    "checked Jordan accepts a non-diagonal Jordan block");
        if (block_form) {
            auto Pinv = matrix_inverse_checked(block_form.value().P).value();
            auto reconstructed = Pinv
                ? SymbolicExpr::multiply(
                      block_form.value().P,
                      SymbolicExpr::multiply(block_form.value().J, Pinv))->simplify()
                : nullptr;
            EXPECT_EQ_EXPR(reconstructed, jordan_block->simplify(),
                           "Jordan block satisfies A=P*J*P^-1");
        }

        auto svd = svd_decomposition_checked(mat2(2,0, 0,3));
        EXPECT_TRUE(svd.has_value(), "checked SVD succeeds on exact nonnegative diagonal matrix");
        if (svd) {
            auto reconstructed = SymbolicExpr::multiply(
                svd.value().U,
                SymbolicExpr::multiply(svd.value().S, SymbolicExpr::transpose(svd.value().V)))->simplify();
            EXPECT_EQ_EXPR(reconstructed, mat2(2,0, 0,3)->simplify(),
                           "checked SVD reconstructs exact diagonal input");
        }
    }

    {
        auto non_square = SymbolicExpr::matrix({{num(1), num(2)}});
        auto bad_lu = lu_decomposition_checked(non_square);
        EXPECT_TRUE(!bad_lu && bad_lu.error().code == CasErrc::InvalidArgument,
                    "checked LU rejects non-square matrix");

        auto needs_pivot_lu = lu_decomposition_checked(mat2(0,1, 1,0));
        EXPECT_TRUE(needs_pivot_lu.has_value(),
                    "checked PLU succeeds when row pivoting is required");
        if (needs_pivot_lu) {
            auto pivoted = SymbolicExpr::multiply(
                needs_pivot_lu.value().P, mat2(0,1, 1,0))->simplify();
            auto reconstructed = SymbolicExpr::multiply(
                needs_pivot_lu.value().L,
                needs_pivot_lu.value().U)->simplify();
            EXPECT_EQ_EXPR(reconstructed, pivoted,
                           "pivoting PLU satisfies P*A=L*U");
        }

        auto x = SymbolicExpr::variable("x");
        auto symbolic_lu_input = SymbolicExpr::matrix({{x, num(1)}, {num(1), num(1)}});
        auto symbolic_lu = lu_decomposition_checked(symbolic_lu_input);
        EXPECT_TRUE(!symbolic_lu && symbolic_lu.error().code == CasErrc::Inconclusive,
                    "checked PLU requires exact rational entries or proved pivots");

        auto bad_qr = qr_decomposition_checked(num(1));
        EXPECT_TRUE(!bad_qr && bad_qr.error().code == CasErrc::InvalidArgument,
                    "checked QR rejects non-matrix input");

        auto dependent_qr = qr_decomposition_checked(mat2(1,0, 2,0));
        EXPECT_TRUE(!dependent_qr && dependent_qr.error().code == CasErrc::Inconclusive,
                    "checked QR reports Inconclusive for rank-deficient columns");

        auto wide_qr = qr_decomposition_checked(
            SymbolicExpr::matrix({{num(1), num(0), num(0)}, {num(0), num(1), num(0)}}));
        EXPECT_TRUE(!wide_qr && wide_qr.error().code == CasErrc::Inconclusive,
                    "checked QR reports Inconclusive outside tall/full-column-rank support");

        auto symbolic_qr_input = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
        auto symbolic_qr = qr_decomposition_checked(symbolic_qr_input);
        EXPECT_TRUE(!symbolic_qr && symbolic_qr.error().code == CasErrc::Inconclusive,
                    "checked QR requires proven exact rational full-column-rank support");

        auto null_chol = cholesky_decomposition_checked(nullptr);
        EXPECT_TRUE(!null_chol && null_chol.error().code == CasErrc::InvalidArgument,
                    "checked Cholesky rejects null input");

        auto non_spd_chol = cholesky_decomposition_checked(mat2(-1,0, 0,1));
        EXPECT_TRUE(!non_spd_chol && non_spd_chol.error().code == CasErrc::DomainError,
                    "checked Cholesky rejects proven non-SPD matrices");

        auto semidefinite_chol = cholesky_decomposition_checked(mat2(1,0, 0,0));
        EXPECT_TRUE(!semidefinite_chol && semidefinite_chol.error().code == CasErrc::DomainError,
                    "checked Cholesky rejects positive semidefinite matrices");

        auto symbolic_spd = SymbolicExpr::matrix({{x, num(0)}, {num(0), num(1)}});
        auto symbolic_chol = cholesky_decomposition_checked(symbolic_spd);
        EXPECT_TRUE(!symbolic_chol && symbolic_chol.error().code == CasErrc::Inconclusive,
                    "checked Cholesky requires a proven exact rational SPD matrix");

        auto bad_jordan = jordan_form_checked(non_square);
        EXPECT_TRUE(!bad_jordan && bad_jordan.error().code == CasErrc::InvalidArgument,
                    "checked Jordan rejects non-square matrix");

        auto non_diagonal_jordan = jordan_form_checked(mat2(1,1, 0,1));
        EXPECT_TRUE(non_diagonal_jordan.has_value(),
                    "checked Jordan supports an exact rational Jordan block");

        auto symbolic_jordan = jordan_form_checked(symbolic_spd);
        EXPECT_TRUE(!symbolic_jordan && symbolic_jordan.error().code == CasErrc::Inconclusive,
                    "checked Jordan requires exact rational entries or proved chains");

        auto non_diagonal_svd = svd_decomposition_checked(mat2(1,1, 0,1));
        EXPECT_TRUE(!non_diagonal_svd &&
                        non_diagonal_svd.error().code == CasErrc::Inconclusive,
                    "checked SVD remains explicit when a complete singular basis is unproved");

        auto negative_diagonal_svd = svd_decomposition_checked(mat2(-1,0, 0,1));
        EXPECT_TRUE(negative_diagonal_svd.has_value(),
                    "checked SVD absorbs diagonal signs into singular vectors");

        auto symbolic_svd = svd_decomposition_checked(symbolic_spd);
        EXPECT_TRUE(!symbolic_svd && symbolic_svd.error().code == CasErrc::Inconclusive,
                    "checked SVD requires exact rational entries or proved eigenspaces");
    }

    {
        CancellationToken token;
        token.cancel();
        ComputationContext cancelled_context({}, token);
        auto cancelled = lu_decomposition_checked(mat2(1,0, 0,1), cancelled_context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked LU observes cancellation");

        ResourceLimits limits;
        limits.max_steps = 1;
        ComputationContext limited_context(limits);
        auto limited = jordan_form_checked(mat2(2,0, 0,3), limited_context);
        EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                    "checked Jordan observes exhausted step budget");
    }

    return TEST_REPORT();
}
