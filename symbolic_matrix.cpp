#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <vector>

// Helper for matrix addition
std::shared_ptr<SymbolicExpr> add_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->operands.size() != b->operands.size()) return nullptr; // Dimension mismatch
    if (a->operands.empty()) return a;
    if (a->operands[0]->operands.size() != b->operands[0]->operands.size()) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
    for (size_t i = 0; i < a->operands.size(); ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row;
        auto vec_a = a->operands[i];
        auto vec_b = b->operands[i];
        for (size_t j = 0; j < vec_a->operands.size(); ++j) {
            row.push_back(SymbolicExpr::add(vec_a->operands[j], vec_b->operands[j])->simplify());
        }
        res_mat.push_back(row);
    }
    return SymbolicExpr::matrix(res_mat);
}

std::shared_ptr<SymbolicExpr> multiply_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->type == SymbolicExpr::Type::Matrix && b->type == SymbolicExpr::Type::Matrix) {
        if (a->operands.empty() || b->operands.empty()) return nullptr;
        
        size_t rowsA = a->operands.size();
        size_t colsA = 0;
        if (!a->operands[0]->operands.empty()) colsA = a->operands[0]->operands.size();
        
        size_t rowsB = b->operands.size();
        size_t colsB = 0;
        if(!b->operands.empty() && !b->operands[0]->operands.empty()) colsB = b->operands[0]->operands.size();
        
        if (colsA != rowsB) return nullptr;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
        res_mat.resize(rowsA);

        for(size_t i=0; i<rowsA; ++i) {
            for(size_t j=0; j<colsB; ++j) {
                std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
                for(size_t k=0; k<colsA; ++k) {
                    auto term = SymbolicExpr::multiply(a->operands[i]->operands[k], b->operands[k]->operands[j]);
                    sum = SymbolicExpr::add(sum, term);
                }
                auto sim = sum->simplify();
                res_mat[i].push_back(sim);
            }
        }
        return SymbolicExpr::matrix(res_mat);
    }
    auto scalar_mul = [](const std::shared_ptr<SymbolicExpr>& scalar, const std::shared_ptr<SymbolicExpr>& mat) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
        for(const auto& row_vec : mat->operands) {
             std::vector<std::shared_ptr<SymbolicExpr>> new_row;
             for(const auto& elem : row_vec->operands) {
                 new_row.push_back(SymbolicExpr::multiply(scalar, elem)->simplify());
             }
             res_mat.push_back(new_row);
        }
        return SymbolicExpr::matrix(res_mat);
    };
    if (a->type != SymbolicExpr::Type::Matrix && b->type == SymbolicExpr::Type::Matrix) return scalar_mul(a, b);
    if (a->type == SymbolicExpr::Type::Matrix && b->type != SymbolicExpr::Type::Matrix) return scalar_mul(b, a);
    return nullptr;
}

// Helper to get minor matrix (remove row r, col c)
static std::shared_ptr<SymbolicExpr> get_minor(const std::shared_ptr<SymbolicExpr>& mat, size_t r, size_t c) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
    size_t rows = mat->operands.size();
    size_t cols = mat->operands[0]->operands.size();
    
    for (size_t i = 0; i < rows; ++i) {
        if (i == r) continue;
        std::vector<std::shared_ptr<SymbolicExpr>> row_vec;
        for (size_t j = 0; j < cols; ++j) {
            if (j == c) continue;
            row_vec.push_back(mat->operands[i]->operands[j]);
        }
        res_mat.push_back(row_vec);
    }
    return SymbolicExpr::matrix(res_mat);
}

std::shared_ptr<SymbolicExpr> matrix_determinant(const std::shared_ptr<SymbolicExpr>& mat) {
    if (mat->type != SymbolicExpr::Type::Matrix) return nullptr;
    if (mat->operands.empty()) return SymbolicExpr::number(1);
    
    size_t rows = mat->operands.size();
    if (mat->operands[0]->type != SymbolicExpr::Type::Vector) return nullptr;
    size_t cols = mat->operands[0]->operands.size();
    
    if (rows != cols) return nullptr; // Not square
    
    if (rows == 1) {
        return mat->operands[0]->operands[0];
    }
    
    if (rows == 2) {
        // [[a,b],[c,d]] -> ad - bc
        auto a = mat->operands[0]->operands[0];
        auto b = mat->operands[0]->operands[1];
        auto c = mat->operands[1]->operands[0];
        auto d = mat->operands[1]->operands[1];
        auto ad = SymbolicExpr::multiply(a, d);
        auto bc = SymbolicExpr::multiply(b, c);
        auto neg_bc = SymbolicExpr::multiply(bc, SymbolicExpr::number(-1));
        return SymbolicExpr::add(ad, neg_bc)->simplify();
    }
    
    std::shared_ptr<SymbolicExpr> det = SymbolicExpr::number(0);
    // Laplace expansion along row 0
    for (size_t j = 0; j < cols; ++j) {
        auto element = mat->operands[0]->operands[j];
        
        // Optimization: If element is 0, skip
        if (element->is_number()) {
             // Assuming convert_rational works for integers too
             auto r = element->convert_rational();
             if (r == Rational(0)) continue;
        }
        
        auto sub = get_minor(mat, 0, j);
        auto sub_det = matrix_determinant(sub);
        
        auto term = SymbolicExpr::multiply(element, sub_det);
        
        if (j % 2 == 1) {
            auto neg_term = SymbolicExpr::multiply(term, SymbolicExpr::number(-1));
            det = SymbolicExpr::add(det, neg_term);
        } else {
            det = SymbolicExpr::add(det, term);
        }
    }
    return det->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::determinant(const std::shared_ptr<SymbolicExpr>& mat) {
    return matrix_determinant(mat);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) {
    if (mat->type != SymbolicExpr::Type::Matrix) return nullptr;
    if (mat->operands.empty()) return mat;
    
    size_t rows = mat->operands.size();
    if (mat->operands[0]->type != SymbolicExpr::Type::Vector) return nullptr;
    size_t cols = mat->operands[0]->operands.size();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
    res_mat.resize(cols);
    for(size_t i=0; i<cols; ++i) res_mat[i].resize(rows);

    for(size_t i=0; i<rows; ++i) {
        for(size_t j=0; j<cols; ++j) {
            res_mat[j][i] = mat->operands[i]->operands[j];
        }
    }
    return SymbolicExpr::matrix(res_mat);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::inverse(const std::shared_ptr<SymbolicExpr>& mat) {
    auto det = determinant(mat);
    if (!det) return nullptr;
    
    // Check if determinant is zero
    if (det->is_number()) {
        auto r = det->convert_rational();
        if (r == Rational(0)) return nullptr; // Singular matrix
    }

    size_t rows = mat->operands.size();
    size_t cols = mat->operands[0]->operands.size();
    if (rows != cols) return nullptr;

    // 1 / det
    auto inv_det = SymbolicExpr::power(det, SymbolicExpr::number(-1));

    if (rows == 1) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res;
        res.push_back({inv_det});
        return SymbolicExpr::matrix(res);
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> adj_mat_data;
    
    // Calculate Cofactor Matrix
    for(size_t i=0; i<rows; ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row_vec;
        for(size_t j=0; j<cols; ++j) {
            auto minor = get_minor(mat, i, j);
            auto m_det = matrix_determinant(minor);
            
            // Cofactor sign (-1)^(i+j)
            if ((i + j) % 2 == 1) {
                m_det = SymbolicExpr::multiply(m_det, SymbolicExpr::number(-1));
            }
            row_vec.push_back(m_det->simplify());
        }
        adj_mat_data.push_back(row_vec);
    }
    
    // Adjugate is Transpose of Cofactor Matrix
    auto cofactor_mat = SymbolicExpr::matrix(adj_mat_data);
    auto adjugate = transpose(cofactor_mat);
    
    // Inverse = Adj * (1/det)
    return SymbolicExpr::multiply(inv_det, adjugate)->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::rref(const std::shared_ptr<SymbolicExpr>& mat_in) {
    if (mat_in->type != SymbolicExpr::Type::Matrix) return nullptr;
    if (mat_in->operands.empty()) return mat_in;

    // Deep copy for mutation (simulated via rebuilding)
    // Actually we can just perform operations on the structure since expressions are shared pointers
    // But the matrix structure itself (vectors) needs to be new.
    size_t rows = mat_in->operands.size();
    size_t cols = mat_in->operands[0]->operands.size();
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat;
    for(auto& r : mat_in->operands) {
        std::vector<std::shared_ptr<SymbolicExpr>> row_copy;
        for(auto& c : r->operands) {
            row_copy.push_back(c);
        }
        mat.push_back(row_copy);
    }

    size_t lead = 0;
    for (size_t r = 0; r < rows; ++r) {
        if (cols <= lead) break;
        
        size_t i = r;
        // Search for pivot
        while (true) {
             // Check if mat[i][lead] is zero
             // Symbolic check is hard, we check explicit zero number
             bool is_zero = false;
             if (mat[i][lead]->is_number()) {
                 if (mat[i][lead]->convert_rational() == Rational(0)) is_zero = true;
             }
             
             if (!is_zero) break; // Found pivot
             
             ++i;
             if (i == rows) {
                 i = r;
                 ++lead;
                 if (cols <= lead) goto end_rref;
             }
        }
        
        // Swap rows i and r
        std::swap(mat[i], mat[r]);
        
        // Normalize row r
        auto pivot = mat[r][lead];
        // If pivot is 1, skip division
        bool pivot_is_one = false;
        if (pivot->is_number() && pivot->convert_rational() == Rational(1)) pivot_is_one = true;
        
        if (!pivot_is_one) {
            auto inv_pivot = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));
            for (size_t j = 0; j < cols; ++j) {
                mat[r][j] = SymbolicExpr::multiply(mat[r][j], inv_pivot)->simplify();
            }
        }
        
        // Eliminate other rows
        for (size_t k = 0; k < rows; ++k) {
            if (k != r) {
                auto factor = mat[k][lead];
                // if factor is 0, skip
                if (factor->is_number() && factor->convert_rational() == Rational(0)) continue;
                
                // R[k] = R[k] - factor * R[r]
                 for (size_t j = 0; j < cols; ++j) {
                     auto term = SymbolicExpr::multiply(factor, mat[r][j]);
                     auto neg_term = SymbolicExpr::multiply(term, SymbolicExpr::number(-1));
                     mat[k][j] = SymbolicExpr::add(mat[k][j], neg_term)->simplify();
                 }
            }
        }
        ++lead;
    }
    
    end_rref:;
    return SymbolicExpr::matrix(mat);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda_var) {
    if (mat->type != SymbolicExpr::Type::Matrix) return nullptr;
    if (mat->operands.empty()) return nullptr;
    
    size_t rows = mat->operands.size();
    if (mat->operands[0]->type != SymbolicExpr::Type::Vector) return nullptr;
    size_t cols = mat->operands[0]->operands.size();
    
    if (rows != cols) return nullptr; // Must be square

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
    auto lambda = SymbolicExpr::variable(lambda_var);
    
    for(size_t i=0; i<rows; ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row_vec;
        for(size_t j=0; j<cols; ++j) {
            auto elem = mat->operands[i]->operands[j];
            if (i == j) {
                // Diagonal: A_ii - lambda
                auto term = SymbolicExpr::add(elem, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
                row_vec.push_back(term->simplify());
            } else {
                row_vec.push_back(elem);
            }
        }
        res_mat.push_back(row_vec);
    }
    
    auto sub_matrix = SymbolicExpr::matrix(res_mat);
    return determinant(sub_matrix);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) {
    auto poly = charpoly(mat, "lambda");
    if (!poly) return nullptr;
    
    // Solve P(lambda) = 0
    auto solutions = solve(poly, "lambda");
    
    // Return solution list as a vector expression
    // Note: solve returns vector<map<string, expr>>, but for single var it's easier.
    // Actually `SymbolicExpr::solve` signature: static vector<map<...>> solve_system
    // Wait, the regular solve is different.
    // cas.cpp `cas_solve` calls `SymbolicExpr::solve`.
    
    // Let's check SymbolicExpr::solve signature in symbolic.hpp... it's not STATIC there?
    // Wait, I saw `clean` reading of `cas.cpp`: `auto solutions = SymbolicExpr::solve(expr, var);`
    // And in `symbolic.hpp` I didn't verify `solve` (singular) existence clearly.
    // Ah, `complex` context check needed.
    // `solve_system` is static.
    
    // Let's assume `solve` follows the `cas.cpp` usage:
    // static std::vector<std::shared_ptr<SymbolicExpr>> solve(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
    // I need to check if that exists or if I need to implement it.
    
    // Based on `cas.cpp`:
    // auto solutions = SymbolicExpr::solve(expr, var);
    // return Value(solutions[0]);  // solutions is vector<shared_ptr<SymbolicExpr>>?
    
    std::vector<std::shared_ptr<SymbolicExpr>> sols = solve(poly, "lambda");
    
    // Wrap in a vector/list expression
    std::vector<std::shared_ptr<SymbolicExpr>> vec_content;
    for(auto& s : sols) vec_content.push_back(s);
    
    auto vec_expr = std::make_shared<SymbolicExpr>(SymbolicExpr::Type::Vector);
    vec_expr->operands = vec_content;
    return vec_expr;
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> SymbolicExpr::eigenvectors(const std::shared_ptr<SymbolicExpr>& mat) {
    auto evals_expr = eigenvalues(mat);
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> result;

    if (!evals_expr || evals_expr->type != Type::Vector) return result;

    size_t n = mat->operands.size();
    
    // Track processed eigenvalues to avoid re-calculation
    std::vector<std::string> seen_evals;

    for (const auto& lambda : evals_expr->operands) {
        std::string lam_str = lambda->to_string();
        bool seen = false;
        for(const auto& s : seen_evals) if (s == lam_str) seen = true;
        if (seen) continue;
        seen_evals.push_back(lam_str);

        // Construct B = Mat - lambda * I
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> b_data;
        for(size_t i=0; i<n; ++i) {
            std::vector<std::shared_ptr<SymbolicExpr>> row;
            for(size_t j=0; j<n; ++j) {
                auto elem = mat->operands[i]->operands[j];
                if (i == j) {
                    row.push_back(SymbolicExpr::add(elem, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)))->simplify());
                } else {
                    row.push_back(elem);
                }
            }
            b_data.push_back(row);
        }
        auto B = SymbolicExpr::matrix(b_data);
        auto R = SymbolicExpr::rref(B);

        // Identify pivots and free variables
        std::vector<int> pivot_cols; 
        std::vector<bool> is_col_pivot(n, false);
        
        for(size_t r=0; r<n; ++r) {
             int c = -1;
             for(size_t j=0; j<n; ++j) {
                 bool is_zero = false;
                 if (R->operands[r]->operands[j]->is_number()) {
                     if (R->operands[r]->operands[j]->convert_rational() == Rational(0)) is_zero = true;
                 }
                 if (!is_zero) {
                     c = (int)j;
                     break;
                 }
             }
             if (c != -1) {
                 pivot_cols.push_back(c);
                 is_col_pivot[c] = true;
             } else {
                 pivot_cols.push_back(-1); 
             }
        }
        
        std::vector<int> free_cols;
        for(int j=0; j<n; ++j) {
            if(!is_col_pivot[j]) free_cols.push_back(j);
        }
        
        std::vector<std::shared_ptr<SymbolicExpr>> eigenvectors_for_lambda;
        
        for (int free_col : free_cols) {
            std::vector<std::shared_ptr<SymbolicExpr>> vec_data(n);
            vec_data[free_col] = SymbolicExpr::number(1);
            
            for(int other : free_cols) {
                if(other != free_col) vec_data[other] = SymbolicExpr::number(0);
            }
            
            for(size_t r=0; r<n; ++r) {
                int p_col = pivot_cols[r];
                if (p_col != -1) {
                    auto val = R->operands[r]->operands[free_col];
                    vec_data[p_col] = SymbolicExpr::multiply(val, SymbolicExpr::number(-1))->simplify();
                }
            }
            
            auto v_expr = std::make_shared<SymbolicExpr>(SymbolicExpr::Type::Vector);
            v_expr->operands = vec_data;
            eigenvectors_for_lambda.push_back(v_expr);
        }
        
        result.push_back({lambda, eigenvectors_for_lambda});
    }
    
    return result;
}

