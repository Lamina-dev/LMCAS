#include "../include/symbolic.hpp"
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/poly_utils.hpp"

using namespace lamina;






std::shared_ptr<SymbolicExpr> get_coeff(const Polynomial<SymbolicPolyCoeff>& p, int deg) {
    if (deg < 0 || deg > p.degree()) return SymbolicExpr::number(0);
    return p.coeffs[deg].val;
}

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
    if (!eq) return {};
    auto simplified_eq = eq->simplify();
    
    
    if (auto rel = std::dynamic_pointer_cast<RelationalNode>(simplified_eq->root)) {
        auto L = std::make_shared<SymbolicExpr>(rel->left);
        auto R = std::make_shared<SymbolicExpr>(rel->right);
        auto diff = SymbolicExpr::add(L, SymbolicExpr::multiply(R, SymbolicExpr::number(-1)));
        
        
        auto poly = symbolic_to_poly<SymbolicPolyCoeff>(diff, var_name);
        
        if (!poly.is_zero()) {
            int max_deg = poly.degree();
            if (max_deg == 1) {
                auto a = get_coeff(poly, 1);
                auto b = get_coeff(poly, 0);
                
                auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
                auto a_inv = SymbolicExpr::power(a, SymbolicExpr::number(-1));
                
                bool flip = false;
                
                try {
                    auto a_simp = a->simplify();
                    if (auto num_a = std::dynamic_pointer_cast<NumberNode>(a_simp->root)) {
                        if (std::holds_alternative<double>(num_a->value)) {
                             if (std::get<double>(num_a->value) < 0) flip = true;
                        } else if (std::holds_alternative<BigInt>(num_a->value)) {
                             if (std::get<BigInt>(num_a->value).IsNegative()) flip = true;
                        } else if (std::holds_alternative<Rational>(num_a->value)) {
                             if (std::get<Rational>(num_a->value).get_numerator().IsNegative()) flip = true;
                        }
                    }
                    
                } catch(...) {}
                
                RelationalNode::Op new_op = rel->op;
                if (flip) {
                    switch(rel->op) {
                        case RelationalNode::Op::LT: new_op = RelationalNode::Op::GT; break;
                        case RelationalNode::Op::GT: new_op = RelationalNode::Op::LT; break;
                        case RelationalNode::Op::LEQ: new_op = RelationalNode::Op::GEQ; break;
                        case RelationalNode::Op::GEQ: new_op = RelationalNode::Op::LEQ; break;
                        default: break;
                    }
                }
                
                auto rhs = SymbolicExpr::multiply(neg_b, a_inv)->simplify();
                auto var_node = std::make_shared<VariableNode>(var_name);
                auto res_node = std::make_shared<RelationalNode>(var_node, rhs->root, new_op);
                
                return { std::make_shared<SymbolicExpr>(res_node) };
            }
        }
        return {}; 
    }
    
    std::shared_ptr<SymbolicExpr> diff_expr = simplified_eq;
    if (auto rel = std::dynamic_pointer_cast<RelationalNode>(simplified_eq->root)) {
        if (rel->op == RelationalNode::Op::EQ) {
             auto L = std::make_shared<SymbolicExpr>(rel->left);
             auto R = std::make_shared<SymbolicExpr>(rel->right);
             diff_expr = SymbolicExpr::add(L, SymbolicExpr::multiply(R, SymbolicExpr::number(-1)));
        }
    }
    
    
    
    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(diff_expr, var_name);
    
    if (poly.is_zero()) return {}; 

    int max_deg = poly.degree();

    if (max_deg == 1) {
        
        auto a = get_coeff(poly, 1); 
        auto b = get_coeff(poly, 0); 
        
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto result = SymbolicExpr::divide(neg_b, a);
        return {result->simplify()};
    } else if (max_deg == 2) {
        
        auto a = get_coeff(poly, 2);
        auto b = get_coeff(poly, 1);
        auto c = get_coeff(poly, 0);

        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto ac4 = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(b2, SymbolicExpr::multiply(ac4, SymbolicExpr::number(-1)));

        
        auto sqrt_delta = SymbolicExpr::power(delta, SymbolicExpr::number(0.5));
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);
        
        auto numerator1 = SymbolicExpr::add(neg_b, sqrt_delta);
        auto numerator2 = SymbolicExpr::add(neg_b, SymbolicExpr::multiply(sqrt_delta, SymbolicExpr::number(-1)));
        
        auto x1 = SymbolicExpr::divide(numerator1, two_a);
        auto x2 = SymbolicExpr::divide(numerator2, two_a);
        
        return {x1->simplify(), x2->simplify()};
    }

    return {};
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(const std::vector<std::shared_ptr<SymbolicExpr>>& equations, const std::vector<std::string>& vars) {
    size_t n = vars.size();
    size_t m = equations.size();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(m, std::vector<std::shared_ptr<SymbolicExpr>>(n + 1));
    
    for(size_t i=0; i<m; ++i) {
        auto eq = equations[i]->expand();
        
        auto C = eq;
        for(const auto& v : vars) C = C->substitute(v, SymbolicExpr::number(0));
        C = C->simplify();
        A[i][n] = SymbolicExpr::multiply(C, SymbolicExpr::number(-1)); 
        
        for(size_t j=0; j<n; ++j) {
            auto p = symbolic_to_poly<SymbolicPolyCoeff>(eq, vars[j]);
            if (p.degree() >= 1) {
                 A[i][j] = get_coeff(p, 1);
            } else {
                 A[i][j] = SymbolicExpr::number(0);
            }
        }
    }

    std::vector<size_t> pivot_col_for_row;
    int sign;
    gaussian_eliminate(A, m, n, pivot_col_for_row, sign);
    
    std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
    for(const auto& v : vars) solution[v] = SymbolicExpr::number(0);

    for(size_t r=0; r<m; ++r) {
        size_t pivot_col = pivot_col_for_row[r];
        if (pivot_col != (size_t)-1) {
            solution[vars[pivot_col]] = A[r][n];
        }
    }
    return { solution };
}






std::shared_ptr<SymbolicExpr> SymbolicExpr::determinant(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat || mat->get_type() != Type::Matrix) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    if (!mat_node || mat_node->rows != mat_node->cols) return SymbolicExpr::number(0);

    size_t n = mat_node->rows;
    if (n == 1) return std::make_shared<SymbolicExpr>(mat_node->get(0,0));
    if (n == 2) {
        auto a = std::make_shared<SymbolicExpr>(mat_node->get(0,0));
        auto b = std::make_shared<SymbolicExpr>(mat_node->get(0,1));
        auto c = std::make_shared<SymbolicExpr>(mat_node->get(1,0));
        auto d = std::make_shared<SymbolicExpr>(mat_node->get(1,1));
        return SymbolicExpr::add(SymbolicExpr::multiply(a,d), SymbolicExpr::multiply(SymbolicExpr::multiply(b,c), SymbolicExpr::number(-1)))->simplify();
    }

    if (n > 3) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
        for(size_t i=0; i<n; ++i) {
            for(size_t j=0; j<n; ++j) {
                A[i][j] = std::make_shared<SymbolicExpr>(mat_node->get(i,j));
            }
        }
        
        int sign = 1;
        auto det = SymbolicExpr::number(1);
        
        for (size_t col = 0; col < n; ++col) {
            size_t pivot_row = col;
            while (pivot_row < n && A[pivot_row][col]->is_zero()) {
                pivot_row++;
            }
            if (pivot_row == n) return SymbolicExpr::number(0);
            
            if (pivot_row != col) {
                std::swap(A[col], A[pivot_row]);
                sign = -sign;
            }
            
            auto pivot = A[col][col];
            det = SymbolicExpr::multiply(det, pivot)->simplify();
            auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));
            
            for (size_t r = col + 1; r < n; ++r) {
                auto factor = A[r][col];
                if (!factor->is_zero()) {
                    auto mult = SymbolicExpr::multiply(factor, pivot_inv)->simplify();
                    auto neg_mult = SymbolicExpr::multiply(mult, SymbolicExpr::number(-1));
                    for (size_t c = col + 1; c < n; ++c) {
                        auto term = SymbolicExpr::multiply(neg_mult, A[col][c]);
                        A[r][c] = SymbolicExpr::add(A[r][c], term)->simplify();
                    }
                    A[r][col] = SymbolicExpr::number(0);
                }
            }
        }
        
        if (sign == -1) {
            det = SymbolicExpr::multiply(det, SymbolicExpr::number(-1))->simplify();
        }
        return det;
    }

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for(size_t c=0; c<n; ++c) {
        auto elem = std::make_shared<SymbolicExpr>(mat_node->get(0,c));
        if (elem->is_zero()) continue;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> minor_data;
        for(size_t r=1; r<n; ++r) {
            std::vector<std::shared_ptr<SymbolicExpr>> row;
            for(size_t k=0; k<n; ++k) {
                if (k == c) continue;
                row.push_back(std::make_shared<SymbolicExpr>(mat_node->get(r,k)));
            }
            minor_data.push_back(row);
        }
        auto minor_mat = SymbolicExpr::matrix(minor_data);
        auto minor_det = SymbolicExpr::determinant(minor_mat);
        
        auto term = SymbolicExpr::multiply(elem, minor_det);
        if (c % 2 == 1) term = SymbolicExpr::multiply(term, SymbolicExpr::number(-1));
        terms.push_back(term);
    }
    
    if (terms.empty()) return SymbolicExpr::number(0);
    auto result = terms[0];
    for(size_t k=1; k<terms.size(); ++k) result = SymbolicExpr::add(result, terms[k]);
    return result->simplify();
}



std::shared_ptr<SymbolicExpr> SymbolicExpr::charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda_name) {
    if (!mat || mat->get_type() != Type::Matrix) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t n = mat_node->rows;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    auto lambda = SymbolicExpr::variable(lambda_name);
    
    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
            auto val = std::make_shared<SymbolicExpr>(mat_node->get(i,j));
            if (i == j) {
                data[i][j] = SymbolicExpr::add(val, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
            } else {
                data[i][j] = val;
            }
        }
    }
    
    auto poly_mat = SymbolicExpr::matrix(data);
    return SymbolicExpr::determinant(poly_mat);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) {
    auto cp = charpoly(mat, "lambda");
    auto solutions = solve(cp, "lambda");
    
    std::vector<std::shared_ptr<SymbolicExpr>> distinct_solutions;

    std::set<std::string> seen;
    for(auto& s : solutions) {
        auto str = s->to_string();
        if (seen.find(str) == seen.end()) {
            seen.insert(str);
            distinct_solutions.push_back(s);
        }
    }
    
    std::vector<std::shared_ptr<SymbolicNode>> vec_nodes;
    for(auto& s : distinct_solutions) vec_nodes.push_back(s->root);
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_data;
    mat_data.push_back(distinct_solutions);
    return SymbolicExpr::matrix(mat_data);
}


std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> SymbolicExpr::eigenvectors(const std::shared_ptr<SymbolicExpr>& mat) {
    auto evals_expr = eigenvalues(mat);

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> result;
    
    if (!evals_expr || (evals_expr->get_type() != SymbolicExpr::Type::Matrix && evals_expr->get_type() != SymbolicExpr::Type::Vector)) {
        return {}; 
    }
    
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(evals_expr->root);
    size_t num_evals = mat_node->cols;
    
    auto A_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t n = A_node->rows;
    
    for(size_t i=0; i<num_evals; ++i) {
        auto lambda_node = mat_node->get(0, i);
        auto lambda = std::make_shared<SymbolicExpr>(lambda_node); 
        
        std::vector<std::shared_ptr<SymbolicExpr>> equations;
        std::vector<std::string> vars;
        for(size_t k=0; k<n; ++k) vars.push_back("v" + std::to_string(k));
        
        for(size_t i=0; i<n; ++i) {
            std::vector<std::shared_ptr<SymbolicExpr>> terms;
            for(size_t j=0; j<n; ++j) {
                auto a_ij = std::make_shared<SymbolicExpr>(A_node->get(i,j));
                std::shared_ptr<SymbolicExpr> coeff;
                if (i == j) {
                    coeff = SymbolicExpr::add(a_ij, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
                } else {
                    coeff = a_ij;
                }
                
                auto var = SymbolicExpr::variable(vars[j]);
                terms.push_back(SymbolicExpr::multiply(coeff, var));
            }
            auto row_eq = terms[0];
            for(size_t k=1; k<terms.size(); ++k) row_eq = SymbolicExpr::add(row_eq, terms[k]);
            equations.push_back(row_eq->simplify());
        }
        
        auto sols = solve_system(equations, vars);
        
        std::vector<std::shared_ptr<SymbolicExpr>> eigenvec;
        for(const auto& v : vars) eigenvec.push_back(sols[0].at(v));
        
        bool is_non_zero = false;
        for(auto& x : eigenvec) if(!x->is_zero()) is_non_zero = true;
        
        if (is_non_zero) {
             std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data;
             for(auto& val : eigenvec) col_vec_data.push_back({val});
             result.push_back({lambda, {SymbolicExpr::matrix(col_vec_data)}});
        } else {
             if (n > 1) {
                  bool found_vec = false;
                  for(int free_var_idx = n-1; free_var_idx >= 0 && !found_vec; --free_var_idx) {
                       std::vector<std::shared_ptr<SymbolicExpr>> sub_eqs;
                       std::vector<std::string> sub_vars;
                       
                       for(size_t k=0; k<n; ++k) {
                           auto eq_sub = equations[k]->substitute(vars[free_var_idx], SymbolicExpr::number(1))->simplify();
                           if (!eq_sub->is_zero()) {
                               sub_eqs.push_back(eq_sub);
                           }
                       }
                       
                       for(size_t k=0; k<n; ++k) {
                           if (k != (size_t)free_var_idx) sub_vars.push_back(vars[k]);
                       }
                       
                       auto sub_sols = solve_system(sub_eqs, sub_vars);
                       if (!sub_sols.empty()) {
                            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data(n);
                            for(size_t k=0; k<n; ++k) {
                                if (k == (size_t)free_var_idx) col_vec_data[k] = {SymbolicExpr::number(1)};
                                else col_vec_data[k] = {sub_sols[0].at(vars[k])};
                            }
                            
                            auto vec_expr = SymbolicExpr::matrix(col_vec_data); 
                            result.push_back({lambda, {vec_expr}});
                            found_vec = true;
                       }
                  }
             }
        }
    }
    
    return result;
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat || (mat->get_type() != Type::Matrix && mat->get_type() != Type::Vector)) return mat;
    
    auto m_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t r = m_node->rows;
    size_t c = m_node->cols;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> new_data(c, std::vector<std::shared_ptr<SymbolicExpr>>(r));
    
    for(size_t i=0; i<r; ++i) {
        for(size_t j=0; j<c; ++j) {
            new_data[j][i] = std::make_shared<SymbolicExpr>(m_node->get(i,j));
        }
    }
    
    return SymbolicExpr::matrix(new_data);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::rref(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat || (mat->get_type() != Type::Matrix && mat->get_type() != Type::Vector)) return mat;

    auto m_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t rows = m_node->rows;
    size_t cols = m_node->cols;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> data(rows, std::vector<std::shared_ptr<SymbolicExpr>>(cols));
    for(size_t i=0; i<rows; ++i) {
        for(size_t j=0; j<cols; ++j) {
            data[i][j] = std::make_shared<SymbolicExpr>(m_node->get(i,j));
        }
    }
    
    size_t lead = 0;
    for (size_t r = 0; r < rows; ++r) {
        if (cols <= lead) break;
        
        size_t i = r;
        while (true) {
             if (i >= rows) {
                 i = r;
                 lead++;
                 if (cols == lead) goto end_loops;
                 continue;
             }
             auto val = data[i][lead]->simplify();
             if (!val->get_number_value_is_zero()) {
                 break; 
             }
             i++;
        }
        
        std::swap(data[i], data[r]);
        
        auto pivot = data[r][lead];
        auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));
        
        for (size_t j = 0; j < cols; ++j) {
            data[r][j] = SymbolicExpr::multiply(data[r][j], pivot_inv)->simplify();
        }
        
        for (size_t k = 0; k < rows; ++k) {
            if (k != r) {
                auto factor = data[k][lead];
                if (factor->simplify()->get_number_value_is_zero()) continue;
                
                auto neg_factor = SymbolicExpr::multiply(factor, SymbolicExpr::number(-1));
                for (size_t j = 0; j < cols; ++j) {
                    auto term = SymbolicExpr::multiply(neg_factor, data[r][j]);
                    data[k][j] = SymbolicExpr::add(data[k][j], term)->simplify();
                }
            }
        }
        lead++;
    }
    end_loops:;

    return SymbolicExpr::matrix(data);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::inverse(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat || mat->get_type() != Type::Matrix) return nullptr;
    
    auto m_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t n = m_node->rows;
    if (n != m_node->cols) return nullptr; 
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> aug_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(2*n));
    
    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
            aug_data[i][j] = std::make_shared<SymbolicExpr>(m_node->get(i,j));
        }
        for(size_t j=0; j<n; ++j) {
            aug_data[i][n+j] = (i==j ? SymbolicExpr::number(1) : SymbolicExpr::number(0));
        }
    }
    
    auto aug_mat = SymbolicExpr::matrix(aug_data);
    auto rref_mat = rref(aug_mat);
    
    if (!rref_mat) return nullptr;
    auto rref_node = std::dynamic_pointer_cast<MatrixNode>(rref_mat->root);
    
    for(size_t i=0; i<n; ++i) {
        auto diag = std::make_shared<SymbolicExpr>(rref_node->get(i,i))->simplify();
        if (!diag->root->is_one()) return nullptr;
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> inv_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
             inv_data[i][j] = std::make_shared<SymbolicExpr>(rref_node->get(i, n+j));
        }
    }
    
    return SymbolicExpr::matrix(inv_data);
}
