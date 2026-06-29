/**
 * @file matrix_decomposition.cpp
 * @brief 矩阵高级分解算法实现。
 */
#include "matrix_decomposition.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <cmath>
#include <stdexcept>

namespace lamina {

bool lu_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L,
    std::shared_ptr<SymbolicExpr>& U) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return false;
    size_t n = mat->rows;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> L_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    
    for (size_t i = 0; i < n; i++) {
        for (size_t k = i; k < n; k++) {
            auto sum = SymbolicExpr::number(0);
            for (size_t j = 0; j < i; j++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[i][j], U_grid[j][k]));
            auto a_ik = std::make_shared<SymbolicExpr>(mat->get(i, k));
            U_grid[i][k] = SymbolicExpr::add(a_ik, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum));
        }

        for (size_t k = i; k < n; k++) {
            if (i == k) {
                L_grid[i][i] = SymbolicExpr::number(1);
            } else {
                auto sum = SymbolicExpr::number(0);
                for (size_t j = 0; j < i; j++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[k][j], U_grid[j][i]));
                auto a_ki = std::make_shared<SymbolicExpr>(mat->get(k, i));
                L_grid[k][i] = SymbolicExpr::divide(SymbolicExpr::add(a_ki, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum)), U_grid[i][i]);
            }
        }
    }
    
    L = SymbolicExpr::matrix(L_grid);
    U = SymbolicExpr::matrix(U_grid);
    return true;
}

bool qr_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& Q,
    std::shared_ptr<SymbolicExpr>& R) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat) return false;
    size_t m = mat->rows;
    size_t n = mat->cols;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> Q_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> R_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_cols(n, std::vector<std::shared_ptr<SymbolicExpr>>(m));
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < m; i++) A_cols[j][i] = std::make_shared<SymbolicExpr>(mat->get(i, j));
    }
    
    for (size_t j = 0; j < n; j++) {
        std::vector<std::shared_ptr<SymbolicExpr>> u_j = A_cols[j];
        for (size_t i = 0; i < j; i++) {
            auto dot = SymbolicExpr::number(0);
            for (size_t k = 0; k < m; k++) dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(Q_grid[k][i], A_cols[j][k]));
            R_grid[i][j] = dot;
            
            for (size_t k = 0; k < m; k++) {
                u_j[k] = SymbolicExpr::add(u_j[k], SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(R_grid[i][j], Q_grid[k][i])));
            }
        }
        
        auto norm_sq = SymbolicExpr::number(0);
        for (size_t k = 0; k < m; k++) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(u_j[k], u_j[k]));
        R_grid[j][j] = SymbolicExpr::sqrt(norm_sq);
        
        for (size_t k = 0; k < m; k++) Q_grid[k][j] = SymbolicExpr::divide(u_j[k], R_grid[j][j]);
    }
    
    Q = SymbolicExpr::matrix(Q_grid);
    R = SymbolicExpr::matrix(R_grid);
    return true;
}

bool cholesky_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return false;
    size_t n = mat->rows;
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> L_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            auto sum = SymbolicExpr::number(0);
            for (size_t k = 0; k < j; k++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[i][k], L_grid[j][k]));
            auto a_ij = std::make_shared<SymbolicExpr>(mat->get(i, j));
            auto diff = SymbolicExpr::add(a_ij, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum));
            if (i == j) L_grid[i][j] = SymbolicExpr::sqrt(diff);
            else L_grid[i][j] = SymbolicExpr::divide(diff, L_grid[j][j]);
        }
    }
    L = SymbolicExpr::matrix(L_grid);
    return true;
}

static std::vector<std::shared_ptr<SymbolicExpr>> normalize_vec(std::vector<std::shared_ptr<SymbolicExpr>> v) {
    auto norm_sq = SymbolicExpr::number(0);
    for (auto& x : v) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
    auto norm = SymbolicExpr::sqrt(norm_sq);
    std::vector<std::shared_ptr<SymbolicExpr>> res;
    for (auto& x : v) res.push_back(SymbolicExpr::divide(x, norm));
    return res;
}

bool svd_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& U,
    std::shared_ptr<SymbolicExpr>& S,
    std::shared_ptr<SymbolicExpr>& V) {
    
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat) return false;
    
    auto AT = SymbolicExpr::transpose(A);
    auto ATA = SymbolicExpr::multiply(AT, A);
    
    auto eigen_V = SymbolicExpr::eigenvectors(ATA);
    if (eigen_V.empty()) return false; // Fallback fails if matrix > 4x4 or no symbolic roots
    
    size_t n = mat->cols;
    size_t m = mat->rows;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> V_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> sigmas;
    
    for (auto& pair : eigen_V) {
        auto sigma = SymbolicExpr::sqrt(pair.first);
        for (auto& vec_mat : pair.second) {
            auto v_mat_node = std::dynamic_pointer_cast<MatrixNode>(vec_mat->root);
            std::vector<std::shared_ptr<SymbolicExpr>> v_col;
            for (size_t i = 0; i < n; i++) v_col.push_back(std::make_shared<SymbolicExpr>(v_mat_node->get(i, 0)));
            V_cols.push_back(normalize_vec(v_col));
            sigmas.push_back(sigma);
        }
    }
    
    if (V_cols.size() != n) return false; // Incomplete basis
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_cols;
    for (size_t i = 0; i < n; i++) {
        if (!sigmas[i]->root->is_zero()) {
            /// u_i = A v_i / sigma_i
            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> v_col_mat(n, std::vector<std::shared_ptr<SymbolicExpr>>(1));
            for (size_t j = 0; j < n; j++) v_col_mat[j][0] = V_cols[i][j];
            auto Avi = SymbolicExpr::multiply(A, SymbolicExpr::matrix(v_col_mat));
            auto Avi_node = std::dynamic_pointer_cast<MatrixNode>(Avi->root);
            std::vector<std::shared_ptr<SymbolicExpr>> u_col;
            for (size_t j = 0; j < m; j++) {
                u_col.push_back(SymbolicExpr::divide(std::make_shared<SymbolicExpr>(Avi_node->get(j, 0)), sigmas[i]));
            }
            U_cols.push_back(u_col);
        }
    }
    
    /// Complete U basis with standard basis via Gram-Schmidt if needed (simplified for stub replacement)
    for (size_t i = 0; i < m && U_cols.size() < m; i++) {
        std::vector<std::shared_ptr<SymbolicExpr>> ei(m, SymbolicExpr::number(0));
        ei[i] = SymbolicExpr::number(1);
        auto ui = ei;
        for (auto& uj : U_cols) {
            auto dot = SymbolicExpr::number(0);
            for (size_t k = 0; k < m; k++) dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(uj[k], ei[k]));
            for (size_t k = 0; k < m; k++) ui[k] = SymbolicExpr::add(ui[k], SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(dot, uj[k])));
        }
        auto norm_sq = SymbolicExpr::number(0);
        for (auto& x : ui) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
        if (!norm_sq->root->is_zero()) {
            U_cols.push_back(normalize_vec(ui));
        }
    }
    
    if (U_cols.size() != m) return false;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(m));
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) U_grid[i][j] = U_cols[j][i];
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> V_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) V_grid[i][j] = V_cols[j][i];
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> S_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < std::min(m, n); i++) S_grid[i][i] = sigmas[i];
    
    U = SymbolicExpr::matrix(U_grid);
    S = SymbolicExpr::matrix(S_grid);
    V = SymbolicExpr::matrix(V_grid);
    return true;
}

std::shared_ptr<SymbolicExpr> matrix_exp(
    const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return SymbolicExpr::exp(A);
    size_t n = mat->rows;
    
    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return SymbolicExpr::exp(A); // Fallback to AST node
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (auto& pair : eigen_V) {
        for (auto& vec : pair.second) {
            auto v_mat_node = std::dynamic_pointer_cast<MatrixNode>(vec->root);
            std::vector<std::shared_ptr<SymbolicExpr>> v_col;
            for (size_t i = 0; i < n; i++) v_col.push_back(std::make_shared<SymbolicExpr>(v_mat_node->get(i, 0)));
            P_cols.push_back(v_col);
            evals.push_back(pair.first);
        }
    }
    
    if (P_cols.size() != n) return SymbolicExpr::exp(A); // Not diagonalizable
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) P_grid[i][j] = P_cols[j][i];
    }
    auto P = SymbolicExpr::matrix(P_grid);
    auto P_inv = SymbolicExpr::inverse(P);
    if (!P_inv) return SymbolicExpr::exp(A);
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> expD_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; i++) expD_grid[i][i] = SymbolicExpr::exp(evals[i]);
    auto expD = SymbolicExpr::matrix(expD_grid);
    
    return SymbolicExpr::multiply(P, SymbolicExpr::multiply(expD, P_inv));
}

// ============================================================
/// 迹与秩 (Requirements 52, 58)
// ============================================================

std::shared_ptr<SymbolicExpr> matrix_trace(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return nullptr;
    auto sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < mat->rows; ++i) {
        sum = SymbolicExpr::add(sum, std::make_shared<SymbolicExpr>(mat->get(i, i)));
    }
    return sum->simplify();
}

int matrix_rank(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat) return -1;
    size_t m = mat->rows, n = mat->cols;
    if (m == 0 || n == 0) return 0;

    /// 构造可变副本（数值化为 double，符号项无法判零时视为非零主元）
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> g(m,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            g[i][j] = std::make_shared<SymbolicExpr>(mat->get(i, j))->simplify();

    auto is_zero_expr = [](const std::shared_ptr<SymbolicExpr>& e) {
        return e && e->root && e->root->is_zero();
    };

    int rank = 0;
    size_t row = 0;
    for (size_t col = 0; col < n && row < m; ++col) {
        /// 寻找主元
        size_t pivot = row;
        while (pivot < m && is_zero_expr(g[pivot][col])) ++pivot;
        if (pivot == m) continue;
        std::swap(g[row], g[pivot]);

        auto pval = g[row][col];
        /// 消元
        for (size_t i = 0; i < m; ++i) {
            if (i == row || is_zero_expr(g[i][col])) continue;
            auto factor = SymbolicExpr::divide(g[i][col], pval);
            for (size_t j = col; j < n; ++j) {
                auto sub = SymbolicExpr::multiply(factor, g[row][j]);
                g[i][j] = SymbolicExpr::add(g[i][j],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), sub))->simplify();
            }
        }
        ++row;
        ++rank;
    }
    return rank;
}

// ============================================================
/// Gram-Schmidt 正交化 (Requirement 57)
// ============================================================

std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gram_schmidt(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& vectors,
    bool normalize) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> basis;
    for (const auto& v : vectors) {
        std::vector<std::shared_ptr<SymbolicExpr>> u = v;
        for (const auto& b : basis) {
            /// 投影系数 <v,b>/<b,b>
            auto vb = SymbolicExpr::number(0);
            auto bb = SymbolicExpr::number(0);
            for (size_t k = 0; k < u.size(); ++k) {
                vb = SymbolicExpr::add(vb, SymbolicExpr::multiply(v[k], b[k]));
                bb = SymbolicExpr::add(bb, SymbolicExpr::multiply(b[k], b[k]));
            }
            if (bb->root && bb->root->is_zero()) continue;
            auto coeff = SymbolicExpr::divide(vb, bb);
            for (size_t k = 0; k < u.size(); ++k) {
                u[k] = SymbolicExpr::add(u[k],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1),
                        SymbolicExpr::multiply(coeff, b[k])))->simplify();
            }
        }
        /// 检查 u 是否为零向量（线性相关）
        auto norm_sq = SymbolicExpr::number(0);
        for (auto& x : u) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
        norm_sq = norm_sq->simplify();
        if (norm_sq->root && norm_sq->root->is_zero()) continue;
        basis.push_back(u);
    }
    if (normalize) {
        for (auto& u : basis) {
            auto norm_sq = SymbolicExpr::number(0);
            for (auto& x : u) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
            auto norm = SymbolicExpr::sqrt(norm_sq);
            for (auto& x : u) x = SymbolicExpr::divide(x, norm)->simplify();
        }
    }
    return basis;
}

// ============================================================
/// 矩阵对数 (Requirement 69)
// ============================================================

std::shared_ptr<SymbolicExpr> matrix_log(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return nullptr;
    size_t n = mat->rows;

    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (auto& pr : eigen_V) {
        for (auto& vec : pr.second) {
            auto vnode = std::dynamic_pointer_cast<MatrixNode>(vec->root);
            if (!vnode) return nullptr;
            std::vector<std::shared_ptr<SymbolicExpr>> col;
            for (size_t i = 0; i < n; ++i) col.push_back(std::make_shared<SymbolicExpr>(vnode->get(i, 0)));
            P_cols.push_back(col);
            evals.push_back(pr.first);
        }
    }
    if (P_cols.size() != n) return nullptr;

    /// 检查特征值为正（实数对数存在性）
    for (auto& ev : evals) {
        auto se = ev->simplify();
        if (se->is_number()) {
            double d = se->to_numeric();
            if (d <= 0) return nullptr;
        }
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) P_grid[i][j] = P_cols[j][i];
    auto P = SymbolicExpr::matrix(P_grid);
    auto P_inv = SymbolicExpr::inverse(P);
    if (!P_inv) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> logD(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) logD[i][i] = SymbolicExpr::ln(evals[i]);
    auto logDm = SymbolicExpr::matrix(logD);

    return SymbolicExpr::multiply(P, SymbolicExpr::multiply(logDm, P_inv));
}

// ============================================================
/// Kronecker 积 (Requirement 70)
// ============================================================

std::shared_ptr<SymbolicExpr> kronecker(const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B) {
    auto a = std::dynamic_pointer_cast<MatrixNode>(A->root);
    auto b = std::dynamic_pointer_cast<MatrixNode>(B->root);
    if (!a || !b) return nullptr;
    size_t m = a->rows, n = a->cols, p = b->rows, q = b->cols;
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(m * p,
        std::vector<std::shared_ptr<SymbolicExpr>>(n * q));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j) {
            auto aij = std::make_shared<SymbolicExpr>(a->get(i, j));
            for (size_t k = 0; k < p; ++k)
                for (size_t l = 0; l < q; ++l) {
                    auto bkl = std::make_shared<SymbolicExpr>(b->get(k, l));
                    grid[i * p + k][j * q + l] = SymbolicExpr::multiply(aij, bkl)->simplify();
                }
        }
    return SymbolicExpr::matrix(grid);
}

// ============================================================
/// 矩阵范数 (Requirement 71)
// ============================================================

std::shared_ptr<SymbolicExpr> matrix_norm(const std::shared_ptr<SymbolicExpr>& A,
    const std::string& type) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat) return nullptr;
    size_t m = mat->rows, n = mat->cols;

    auto abs_of = [](const std::shared_ptr<SymbolicExpr>& e) {
        auto node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<SymbolicNode>>{e->root});
        return std::make_shared<SymbolicExpr>(node);
    };

    if (type == "frobenius") {
        auto sum = SymbolicExpr::number(0);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j) {
                auto e = std::make_shared<SymbolicExpr>(mat->get(i, j));
                sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(e, e));
            }
        return SymbolicExpr::sqrt(sum)->simplify();
    }
    if (type == "1") {
        /// 最大列和
        std::shared_ptr<SymbolicExpr> best = nullptr;
        double best_val = -1;
        for (size_t j = 0; j < n; ++j) {
            auto colsum = SymbolicExpr::number(0);
            for (size_t i = 0; i < m; ++i)
                colsum = SymbolicExpr::add(colsum, abs_of(std::make_shared<SymbolicExpr>(mat->get(i, j))));
            colsum = colsum->simplify();
            double v = 0;
            if (colsum->is_number()) v = colsum->to_numeric();
            if (!best || v > best_val) { best = colsum; best_val = v; }
        }
        return best;
    }
    if (type == "inf") {
        /// 最大行和
        std::shared_ptr<SymbolicExpr> best = nullptr;
        double best_val = -1;
        for (size_t i = 0; i < m; ++i) {
            auto rowsum = SymbolicExpr::number(0);
            for (size_t j = 0; j < n; ++j)
                rowsum = SymbolicExpr::add(rowsum, abs_of(std::make_shared<SymbolicExpr>(mat->get(i, j))));
            rowsum = rowsum->simplify();
            double v = 0;
            if (rowsum->is_number()) v = rowsum->to_numeric();
            if (!best || v > best_val) { best = rowsum; best_val = v; }
        }
        return best;
    }
    return nullptr;
}

// ============================================================
/// 二次型 (Requirement 59)
// ============================================================

std::shared_ptr<SymbolicExpr> quadratic_form_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& vars) {
    if (!expr) return nullptr;
    size_t n = vars.size();
    auto e = expr->expand();
    if (!e) e = expr;

    /// A[i][j] = (1/2) ∂²(expr)/∂xᵢ∂xⱼ （对称矩阵，对角为 ∂²/2 即系数本身的一半*2）
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) {
        auto di = e->differentiate(vars[i]);
        for (size_t j = 0; j < n; ++j) {
            auto dij = di->differentiate(vars[j]);
            /// 对二次型，二阶偏导为常数；A_ij = (1/2)·∂²/∂xi∂xj
            grid[i][j] = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), dij)->simplify();
        }
    }
    return SymbolicExpr::matrix(grid);
}

std::string classify_quadratic_form(const std::shared_ptr<SymbolicExpr>& A) {
    /// 使用特征值列表（来自特征多项式求根），避免脆弱的特征向量求解。
    auto evals_expr = SymbolicExpr::eigenvalues(A);
    if (!evals_expr) return "unknown";
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(evals_expr->root);
    if (!mat_node || mat_node->cols == 0) return "unknown";

    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (size_t i = 0; i < mat_node->cols; ++i) {
        evals.push_back(std::make_shared<SymbolicExpr>(mat_node->get(0, i)));
    }
    if (evals.empty()) return "unknown";

    /// 递归数值求值，处理未化简的根式（如 (1/2)*(4^0.5)）。
    std::function<bool(const std::shared_ptr<SymbolicNode>&, double&)> num_eval =
        [&](const std::shared_ptr<SymbolicNode>& n, double& out) -> bool {
        if (!n) return false;
        if (auto nn = std::dynamic_pointer_cast<NumberNode>(n)) {
            SymbolicExpr tmp(n);
            out = tmp.to_numeric();
            return true;
        }
        if (auto a = std::dynamic_pointer_cast<AddNode>(n)) {
            double s = 0, t;
            for (auto& op : a->operands) { if (!num_eval(op, t)) return false; s += t; }
            out = s; return true;
        }
        if (auto m = std::dynamic_pointer_cast<MultiplyNode>(n)) {
            double p = 1, t;
            for (auto& op : m->operands) { if (!num_eval(op, t)) return false; p *= t; }
            out = p; return true;
        }
        if (auto pw = std::dynamic_pointer_cast<PowerNode>(n)) {
            double b, e;
            if (!num_eval(pw->base, b) || !num_eval(pw->exponent, e)) return false;
            out = std::pow(b, e); return true;
        }
        if (auto f = std::dynamic_pointer_cast<FunctionNode>(n)) {
            if (f->type == FunctionNode::FuncType::Sqrt && f->arguments.size() == 1) {
                double a; if (!num_eval(f->arguments[0], a)) return false;
                out = std::sqrt(a); return true;
            }
        }
        return false;
    };

    bool all_pos = true, all_neg = true, any_zero = false;
    bool has_pos = false, has_neg = false;
    const double eps = 1e-9;
    for (auto& ev : evals) {
        auto se = ev->simplify();
        double d;
        if (!num_eval(se->root, d)) return "unknown";
        if (d > eps) { all_neg = false; has_pos = true; }
        else if (d < -eps) { all_pos = false; has_neg = true; }
        else { any_zero = true; all_pos = false; all_neg = false; }
    }
    if (all_pos) return "positive_definite";
    if (all_neg) return "negative_definite";
    if (has_pos && has_neg) return "indefinite";
    if (has_pos && any_zero) return "positive_semidefinite";
    if (has_neg && any_zero) return "negative_semidefinite";
    return "indefinite";
}

// ============================================================
/// Jordan 标准型 (Requirement 56)
// ============================================================

bool jordan_form(const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& J, std::shared_ptr<SymbolicExpr>& P) {
    auto mat = std::dynamic_pointer_cast<MatrixNode>(A->root);
    if (!mat || mat->rows != mat->cols) return false;
    size_t n = mat->rows;

    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return false;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> diag;
    for (auto& pr : eigen_V) {
        for (auto& vec : pr.second) {
            auto vnode = std::dynamic_pointer_cast<MatrixNode>(vec->root);
            if (!vnode) return false;
            std::vector<std::shared_ptr<SymbolicExpr>> col;
            for (size_t i = 0; i < n; ++i) col.push_back(std::make_shared<SymbolicExpr>(vnode->get(i, 0)));
            P_cols.push_back(col);
            diag.push_back(pr.first);
        }
    }
    /// 当特征向量数等于 n 时，矩阵可对角化，Jordan 型即对角阵。
    /// 缺陷情形（重根但特征向量不足）当前不支持广义特征向量链，返回 false。
    if (P_cols.size() != n) return false;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) P_grid[i][j] = P_cols[j][i];
    P = SymbolicExpr::matrix(P_grid);

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> J_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) J_grid[i][i] = diag[i];
    J = SymbolicExpr::matrix(J_grid);
    return true;
}

} // namespace lamina
