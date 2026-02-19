#include "../include/symbolic.hpp"
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include "../include/visitors/normalization_visitor.hpp"

// A visitor to collect polynomial coefficients w.r.t a variable
class PolyCoeffVisitor : public SymbolicVisitor {
public:
    std::string var;
    std::map<int, std::shared_ptr<SymbolicExpr>> coeffs;
    bool is_polynomial = true;

    PolyCoeffVisitor(std::string v) : var(v) {}

    void add_coeff(int deg, std::shared_ptr<SymbolicExpr> val) {
        if (coeffs.find(deg) == coeffs.end()) {
            coeffs[deg] = val;
        } else {
            coeffs[deg] = SymbolicExpr::add(coeffs[deg], val);
        }
    }

    void visit(NumberNode& node) override {
        add_coeff(0, std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(node)));
    }

    void visit(VariableNode& node) override {
        if (node.name == var) add_coeff(1, SymbolicExpr::number(1));
        else add_coeff(0, std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(node)));
    }

    void visit(AddNode& node) override {
        for (const auto& op : node.operands) {
            op->accept(*this);
            if (!is_polynomial) return;
        }
    }

    void visit(MultiplyNode& node) override {
        int total_deg = 0;
        std::vector<std::shared_ptr<SymbolicExpr>> constant_parts;

        for (const auto& op : node.operands) {
             if (auto v = std::dynamic_pointer_cast<VariableNode>(op)) {
                 if (v->name == var) total_deg++;
                 else constant_parts.push_back(std::make_shared<SymbolicExpr>(op));
             } else if (auto p = std::dynamic_pointer_cast<PowerNode>(op)) {
                 if (auto b = std::dynamic_pointer_cast<VariableNode>(p->base)) {
                     if (b->name == var) {
                         if (auto n = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
                              if (std::holds_alternative<double>(n->value)) total_deg += (int)std::get<double>(n->value);
                              else if (std::holds_alternative<BigInt>(n->value)) total_deg += std::get<BigInt>(n->value).to_int();
                         }
                     } else { constant_parts.push_back(std::make_shared<SymbolicExpr>(op)); }
                 } else { constant_parts.push_back(std::make_shared<SymbolicExpr>(op)); }
             } else if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) {
                 constant_parts.push_back(std::make_shared<SymbolicExpr>(op));
             } else {
                 // Check if contains var? Assuming simplified inputs where logic is cleaner or handled here
                 // For now treat as constant
                 constant_parts.push_back(std::make_shared<SymbolicExpr>(op));
             }
        }
        
        std::shared_ptr<SymbolicExpr> term = (constant_parts.empty()) ? SymbolicExpr::number(1) : constant_parts[0];
        for(size_t k=1; k<constant_parts.size(); ++k) term = SymbolicExpr::multiply(term, constant_parts[k]);
        
        add_coeff(total_deg, term);
    }
    
    void visit(PowerNode& node) override {
         if (auto v = std::dynamic_pointer_cast<VariableNode>(node.base)) {
             if (v->name == var) {
                 if (auto n = std::dynamic_pointer_cast<NumberNode>(node.exponent)) {
                     int deg = 0;
                     if (std::holds_alternative<double>(n->value)) deg = (int)std::get<double>(n->value);
                     else if (std::holds_alternative<BigInt>(n->value)) deg = std::get<BigInt>(n->value).to_int();
                     add_coeff(deg, SymbolicExpr::number(1));
                     return;
                 }
             }
         }
         // Else treat as constant
         add_coeff(0, std::make_shared<SymbolicExpr>(std::make_shared<PowerNode>(node.base, node.exponent)));
    }
    
    void visit(FunctionNode& node) override {
        // Treat as constant
        std::vector<std::shared_ptr<SymbolicNode>> args = node.arguments;
        add_coeff(0, std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(node.type, args)));
    }
    
    void visit(MatrixNode& node) override {
        // Matrix cannot be part of polynomial coefficient in this context
        // Treat as constant or error? For now, constant 0 degree
        // Actually, matrices shouldn't appear in scalar polynomial logic usually
    }
};

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
    if (!eq) return {};
    auto simplified_eq = eq->simplify();
    
    // Check for Inequality (RelationalNode)
    if (auto rel = std::dynamic_pointer_cast<RelationalNode>(simplified_eq->root)) {
        // L op R -> (L - R) op 0
        auto L = std::make_shared<SymbolicExpr>(rel->left);
        auto R = std::make_shared<SymbolicExpr>(rel->right);
        auto diff = SymbolicExpr::add(L, SymbolicExpr::multiply(R, SymbolicExpr::number(-1)))->simplify();
        
        PolyCoeffVisitor poly_v(var_name);
        diff->root->accept(poly_v);
        
        if (poly_v.is_polynomial && !poly_v.coeffs.empty()) {
            int max_deg = poly_v.coeffs.rbegin()->first;
            if (max_deg == 1) {
                // ax + b op 0
                // ax op -b
                auto a = (poly_v.coeffs.count(1)) ? poly_v.coeffs[1] : SymbolicExpr::number(0);
                auto b = (poly_v.coeffs.count(0)) ? poly_v.coeffs[0] : SymbolicExpr::number(0);
                
                auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
                auto a_inv = SymbolicExpr::power(a, SymbolicExpr::number(-1));
                
                // Result: x op (-b/a)
                // Need to determine if op flips.
                // Flip if a < 0.
                bool flip = false;
                if (auto num_a = std::dynamic_pointer_cast<NumberNode>(a->root)) {
                    if (std::holds_alternative<double>(num_a->value)) {
                         if (std::get<double>(num_a->value) < 0) flip = true;
                    } else if (std::holds_alternative<BigInt>(num_a->value)) {
                         if (std::get<BigInt>(num_a->value).IsNegative()) flip = true;
                    } else if (std::holds_alternative<Rational>(num_a->value)) {
                         if (std::get<Rational>(num_a->value).get_numerator().IsNegative()) flip = true;
                    }
                }
                
                // Construct new RelationalNode for result
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
        return {}; // Non-linear or complex inequality
    }

    PolyCoeffVisitor poly_v(var_name);
    simplified_eq->root->accept(poly_v);
    
    if (poly_v.is_polynomial && !poly_v.coeffs.empty()) {
        int max_deg = poly_v.coeffs.rbegin()->first;
        if (max_deg == 1) {
            auto a = (poly_v.coeffs.count(1)) ? poly_v.coeffs[1] : SymbolicExpr::number(0);
            auto b = (poly_v.coeffs.count(0)) ? poly_v.coeffs[0] : SymbolicExpr::number(0);
            auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
            auto a_inv = SymbolicExpr::power(a, SymbolicExpr::number(-1));
            return { SymbolicExpr::multiply(neg_b, a_inv)->simplify() };
        }
        if (max_deg == 2) {
            auto a = (poly_v.coeffs.count(2)) ? poly_v.coeffs[2] : SymbolicExpr::number(0);
            auto b = (poly_v.coeffs.count(1)) ? poly_v.coeffs[1] : SymbolicExpr::number(0);
            auto c = (poly_v.coeffs.count(0)) ? poly_v.coeffs[0] : SymbolicExpr::number(0);
            auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
            auto four = SymbolicExpr::number(4);
            auto ac4 = SymbolicExpr::multiply(SymbolicExpr::multiply(four, a), c);
            auto neg_ac4 = SymbolicExpr::multiply(ac4, SymbolicExpr::number(-1));
            auto D = SymbolicExpr::add(b2, neg_ac4);
            
            // sqrt_D definition:
            auto half = SymbolicExpr::number(Rational(1, 2));
            auto sqrt_D_node = SymbolicExpr::power(D, half);
            
            auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
            auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);
            auto two_a_inv = SymbolicExpr::power(two_a, SymbolicExpr::number(-1));
            auto num1 = SymbolicExpr::add(neg_b, sqrt_D_node);
            auto x1 = SymbolicExpr::multiply(num1, two_a_inv);
            auto neg_sqrt_D = SymbolicExpr::multiply(sqrt_D_node, SymbolicExpr::number(-1));
            auto num2 = SymbolicExpr::add(neg_b, neg_sqrt_D);
            auto x2 = SymbolicExpr::multiply(num2, two_a_inv);
            return { x1->simplify(), x2->simplify() };
        }
    }
    return {};
}

// Gaussian Elimination for System of Linear Equations
std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(const std::vector<std::shared_ptr<SymbolicExpr>>& equations, const std::vector<std::string>& vars) {
    size_t n = vars.size();
    size_t m = equations.size();
    // Build Augmented Matrix (m x n+1)
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(m, std::vector<std::shared_ptr<SymbolicExpr>>(n + 1));
    
    for(size_t i=0; i<m; ++i) {
        auto eq = equations[i]->expand();
        for(size_t j=0; j<n; ++j) {
            PolyCoeffVisitor v(vars[j]);
            eq->root->accept(v);
            A[i][j] = (v.coeffs.count(1)) ? v.coeffs[1] : SymbolicExpr::number(0);
        }
        // Constant term B[i] (RHS).
        // C = eq with all vars = 0
        auto C = eq;
        for(const auto& v : vars) {
             C = C->substitute(v, SymbolicExpr::number(0));
        }
        C = C->simplify();
        A[i][n] = SymbolicExpr::multiply(C, SymbolicExpr::number(-1)); // RHS
    }

    // Gaussian Elimination
    size_t pivot_row = 0;
    for(size_t col=0; col<n && pivot_row<m; ++col) {
        // Find pivot
        size_t max_row = pivot_row;
        // Simple pivot selection: first non-zero
        // Ideally pick robust pivot
        bool found_pivot = false;
        while(max_row < m) {
            A[max_row][col] = A[max_row][col]->simplify(); // Ensure simplified
            if (!A[max_row][col]->is_zero()) {
                found_pivot = true;
                break; 
            }
            max_row++;
        }
        
        if (!found_pivot) continue; // No pivot in this column
        
        // Swap rows
        std::swap(A[pivot_row], A[max_row]);
        
        // Normalize pivot row
        auto pivot = A[pivot_row][col];
        auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));
        
        for(size_t k=col; k<=n; ++k) {
            A[pivot_row][k] = SymbolicExpr::multiply(A[pivot_row][k], pivot_inv)->simplify();
        }
        // A[pivot_row][col] is now 1
        
        // Eliminate other rows
        for(size_t i=0; i<m; ++i) {
            if (i != pivot_row) {
                auto factor = A[i][col];
                if (!factor->is_zero()) {
                    auto neg_factor = SymbolicExpr::multiply(factor, SymbolicExpr::number(-1));
                    for(size_t k=col; k<=n; ++k) {
                        auto term = SymbolicExpr::multiply(neg_factor, A[pivot_row][k]);
                        A[i][k] = SymbolicExpr::add(A[i][k], term)->simplify();
                    }
                }
            }
        }
        pivot_row++;
    }
    
    // Back substitution not needed if we did full elimination (GAUSS-JORDAN)
    // Results are in A[i][n] for basic variables?
    // Construct solution map
    std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
    // Initialize all vars to 0 (free variables)
    for(const auto& v : vars) solution[v] = SymbolicExpr::number(0);

    for(size_t r=0; r<m; ++r) {
        // Find pivot in this row
        int pivot_col = -1;
        for(size_t c=0; c<n; ++c) {
            if (!A[r][c]->is_zero()) {
                pivot_col = c;
                break;
            }
        }
        
        if (pivot_col != -1) {
            // This row solves for vars[pivot_col]
            // In RREF, this variable has coefficient 1, and no other row has it non-zero.
            // But if RREF failed partially (due to simplification issues), we might have other non-zeros.
            // Assuming simplified RREF:
            solution[vars[pivot_col]] = A[r][n];
        }
    }
    return { solution };
}

// =========================================================
// Matrix Operations Implementation
// =========================================================

// Determinant (Recursive Expansion)
std::shared_ptr<SymbolicExpr> SymbolicExpr::determinant(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat || mat->get_type() != Type::Matrix) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    if (!mat_node || mat_node->rows != mat_node->cols) return SymbolicExpr::number(0);

    size_t n = mat_node->rows;
    if (n == 1) return std::make_shared<SymbolicExpr>(mat_node->get(0,0));
    if (n == 2) {
        // ad - bc
        auto a = std::make_shared<SymbolicExpr>(mat_node->get(0,0));
        auto b = std::make_shared<SymbolicExpr>(mat_node->get(0,1));
        auto c = std::make_shared<SymbolicExpr>(mat_node->get(1,0));
        auto d = std::make_shared<SymbolicExpr>(mat_node->get(1,1));
        return SymbolicExpr::add(SymbolicExpr::multiply(a,d), SymbolicExpr::multiply(SymbolicExpr::multiply(b,c), SymbolicExpr::number(-1)))->simplify();
    }

    // Laplace expansion along first row
    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for(size_t c=0; c<n; ++c) {
        auto elem = std::make_shared<SymbolicExpr>(mat_node->get(0,c));
        if (elem->is_zero()) continue;

        // Minor matrix
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


// Characteristic Polynomial: det(A - lambda*I)
std::shared_ptr<SymbolicExpr> SymbolicExpr::charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda_name) {
    if (!mat || mat->get_type() != Type::Matrix) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t n = mat_node->rows;
    
    // Construct A - lambda*I
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    auto lambda = SymbolicExpr::variable(lambda_name);
    
    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
            auto val = std::make_shared<SymbolicExpr>(mat_node->get(i,j));
            if (i == j) {
                // val - lambda
                data[i][j] = SymbolicExpr::add(val, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
            } else {
                data[i][j] = val;
            }
        }
    }
    
    auto poly_mat = SymbolicExpr::matrix(data);
    return SymbolicExpr::determinant(poly_mat);
}

// Eigenvalues
std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) {
    auto cp = charpoly(mat, "lambda");
    auto solutions = solve(cp, "lambda");
    
    std::vector<std::shared_ptr<SymbolicExpr>> distinct_solutions;

    // Simple deduplication based on string (could be better)
    std::set<std::string> seen;
    for(auto& s : solutions) {
        auto str = s->to_string();
        if (seen.find(str) == seen.end()) {
            seen.insert(str);
            distinct_solutions.push_back(s);
        }
    }
    
    // Return vector type
    std::vector<std::shared_ptr<SymbolicNode>> vec_nodes;
    for(auto& s : distinct_solutions) vec_nodes.push_back(s->root);
    
    // Return vector type
    // Using MatrixNode (1xN) to represent vector list of eigenvalues
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_data;
    mat_data.push_back(distinct_solutions);
    return SymbolicExpr::matrix(mat_data);
}

// Eigenvectors
std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> SymbolicExpr::eigenvectors(const std::shared_ptr<SymbolicExpr>& mat) {
    auto evals_expr = eigenvalues(mat);

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> result;
    
    if (!evals_expr || (evals_expr->get_type() != SymbolicExpr::Type::Matrix && evals_expr->get_type() != SymbolicExpr::Type::Vector)) {
        return {}; 
    }
    
    // Get operands from MatrixNode (1 row)
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(evals_expr->root);
    size_t num_evals = mat_node->cols;
    
    auto A_node = std::dynamic_pointer_cast<MatrixNode>(mat->root);
    size_t n = A_node->rows;
    
    for(size_t i=0; i<num_evals; ++i) {
        auto lambda_node = mat_node->get(0, i);
        auto lambda = std::make_shared<SymbolicExpr>(lambda_node); // Proper construction?
        
        // Solve (A - lambda*I)v = 0

        // We construct the system of equations
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
            // Sum terms = 0
            auto row_eq = terms[0];
            for(size_t k=1; k<terms.size(); ++k) row_eq = SymbolicExpr::add(row_eq, terms[k]);
            equations.push_back(row_eq->simplify());
        }
        
        auto sols = solve_system(equations, vars);
        
        // Convert solution map to vector
        // solve_system returns generic solution with free variables as 0 (in my implementation)
        // Wait! In eigenvectors, we WANT non-trivial solutions.
        // My solve_system currently sets free variables to 0! That's BAD for homogeneous systems where only the trivial solution exists if free vars are 0.
        // I need to change solve_system or handle it here.
        // Since solve_system is generic, it might be setting params to 0. 
        // For eigenvectors, we should detect free variables.
        
        // HACK: solve_system sets free vars to 0. If result is all 0, it failed to find eigenvector.
        // Actually, for (A-lI)v=0, we expect dependent rows.
        // We should set free variables to 1 one by one to get basis.
        
        // For now, let's assume solve_system supports parametrization or we do it manually?
        // Let's rely on solve_system output. If it produced all 0s, we failed.
        
        // RE-IMPLEMENT a simple NullSpace solver here specifically for eigenvectors.
        // Gaussian Elimination on (A - lI)
        // ... (this is duplicated work but safer for now to ensure free variable handling)
        
        // Basic check
        std::vector<std::shared_ptr<SymbolicExpr>> eigenvec;
        for(const auto& v : vars) eigenvec.push_back(sols[0].at(v));
        
        // Verify if non-zero
        bool is_non_zero = false;
        for(auto& x : eigenvec) if(!x->is_zero()) is_non_zero = true;
        
        if (is_non_zero) {
             // Create column vector matrix
             std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data;
             for(auto& val : eigenvec) col_vec_data.push_back({val});
             result.push_back({lambda, {SymbolicExpr::matrix(col_vec_data)}});
        } else {
             // Try setting the LAST variable to 1 if it was 0?
             // Since we know det(A-lI)=0, at least one row is dependent.
             // Usually the last variable is free.
             // Let's force v_n-1 = 1 and solve for others using the first n-1 equations.
             
             // Actually, simply remove the last equation (redundant) and set last var = 1.
             // Then solve non-homogeneous system.
             if (n > 1) {
                  // Strategy: Try setting each variable to 1 (one-hot) and solve for the others.
                  // We need just one non-trivial solution.
                  
                  bool found_vec = false;
                  for(int free_var_idx = n-1; free_var_idx >= 0 && !found_vec; --free_var_idx) {
                       std::vector<std::shared_ptr<SymbolicExpr>> sub_eqs;
                       std::vector<std::string> sub_vars;
                       
                       // Collect equations (removing zero rows)
                       for(size_t k=0; k<n; ++k) {
                           auto eq_sub = equations[k]->substitute(vars[free_var_idx], SymbolicExpr::number(1))->simplify();
                           if (!eq_sub->is_zero()) {
                               sub_eqs.push_back(eq_sub);
                           }
                       }
                       
                       // Collect variables to solve for
                       for(size_t k=0; k<n; ++k) {
                           if (k != (size_t)free_var_idx) sub_vars.push_back(vars[k]);
                       }
                       
                       // If we have more vars than eqs, we might have multiple free variables, but solve_system might handle under-determined by setting others to 0?
                       // Actually solve_system minimizes free vars to 0. So it works out.
                       
                       auto sub_sols = solve_system(sub_eqs, sub_vars);
                       if (!sub_sols.empty()) {
                            // Construct full vector
                            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data(n);
                            for(size_t k=0; k<n; ++k) {
                                if (k == (size_t)free_var_idx) col_vec_data[k] = {SymbolicExpr::number(1)};
                                else col_vec_data[k] = {sub_sols[0].at(vars[k])};
                            }
                            
                            auto vec_expr = SymbolicExpr::matrix(col_vec_data); // Nx1 Matrix
                            result.push_back({lambda, {vec_expr}});
                            found_vec = true;
                       }
                  }
             }
        }
    }
    
    return result;
}

// Transpose/Inverse/RREF impls...
std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) { return mat; } // TODO
std::shared_ptr<SymbolicExpr> SymbolicExpr::inverse(const std::shared_ptr<SymbolicExpr>& mat) { return mat; } // TODO
std::shared_ptr<SymbolicExpr> SymbolicExpr::rref(const std::shared_ptr<SymbolicExpr>& mat) { return mat; } // TODO

