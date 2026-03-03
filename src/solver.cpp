// Force rebuild for SymbolicFactory::create_add update
#include "solver.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <optional>

namespace lamina {


static bool contains(const SymbolicExpr& expr, const std::string& var) {
    if (!expr.root) return false;
    
    struct ContainsVisitor : public SymbolicVisitor {
        bool found = false;
        std::string target;
        ContainsVisitor(const std::string& t) : target(t) {}
        
        void visit(NumberNode& n) override {}
        void visit(VariableNode& n) override { if (n.name == target) found = true; }
        void visit(AddNode& n) override { for (auto& op : n.operands) if (!found) op->accept(*this); }
        void visit(MultiplyNode& n) override { for (auto& op : n.operands) if (!found) op->accept(*this); }
        void visit(PowerNode& n) override { 
            if (!found) n.base->accept(*this);
            if (!found) n.exponent->accept(*this); 
        }
        void visit(FunctionNode& n) override { 
            for (auto& arg : n.arguments) if (!found) arg->accept(*this);
        }
        void visit(MatrixNode& n) override {} 
        void visit(RelationalNode& n) override {} 
    };
    
    ContainsVisitor v(var);
    expr.root->accept(v);
    return v.found;
}


static std::shared_ptr<SymbolicExpr> to_ptr(const SymbolicExpr& expr) {
    return std::make_shared<SymbolicExpr>(expr);
}







static std::pair<SymbolicExpr, SymbolicExpr> isolate_linear_coeff(const SymbolicExpr& expr, const std::string& var) {
    
    
    
    

    auto expr_ptr = to_ptr(expr); 
    auto A_ptr = expr_ptr->differentiate(var);
    
    
    std::vector<std::shared_ptr<SymbolicNode>> mops;
    mops.push_back(A_ptr->root);
    mops.push_back(SymbolicFactory::create_variable(var));
    auto term_Ax = SymbolicFactory::create_multiply(mops);
    
    
    std::vector<std::shared_ptr<SymbolicNode>> nops;
    nops.push_back(SymbolicFactory::create_number(BigInt(-1)));
    nops.push_back(term_Ax);
    auto neg_term = SymbolicFactory::create_multiply(nops);
    
    
    std::vector<std::shared_ptr<SymbolicNode>> aops;
    aops.push_back(expr.root);
    aops.push_back(neg_term);
    auto B_node = SymbolicFactory::create_add(aops);
    
    
    
    SymbolicExpr B_expr(B_node);
    auto B_simp = to_ptr(B_expr)->simplify(); 
    
    SymbolicExpr A_expr(A_ptr->root);
    auto A_simp = to_ptr(A_expr)->simplify();

    return {SymbolicExpr(A_simp->root), SymbolicExpr(B_simp->root)};
}


std::map<std::string, SymbolicExpr> Solver::solve_linear_system(
    const std::vector<SymbolicExpr>& equations_in, 
    const std::vector<std::string>& variables) 
{
    
    
    
    
    
    
    size_t num_vars = variables.size();
    size_t num_eqs = equations_in.size();
    
    
    std::vector<std::vector<SymbolicExpr>> matrix(num_eqs, std::vector<SymbolicExpr>(num_vars + 1));
    
    for (size_t i = 0; i < num_eqs; ++i) {
        SymbolicExpr current_eq = equations_in[i]; 
        
        
        
        SymbolicExpr constant_part = current_eq;
        
        for (size_t j = 0; j < num_vars; ++j) {
            auto [coeff, remainder] = isolate_linear_coeff(constant_part, variables[j]);
            matrix[i][j] = coeff;
            
            
            
            
            constant_part = remainder;
        }
        
        
        
        
        
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
        ops.push_back(constant_part.root);
        auto neg_const = SymbolicFactory::create_multiply(ops);
        matrix[i][num_vars] = SymbolicExpr(neg_const);
    }
    
    
    size_t pivot_row = 0;
    std::vector<size_t> pivot_col_for_row(num_eqs, -1);
    
    for (size_t col = 0; col < num_vars && pivot_row < num_eqs; ++col) {
        
        size_t sel = pivot_row;
        
        
        while (sel < num_eqs && matrix[sel][col].root->is_zero()) {
            sel++;
        }
        
        if (sel == num_eqs) continue; 
        
        
        std::swap(matrix[sel], matrix[pivot_row]);
        pivot_col_for_row[pivot_row] = col;
        
        
        SymbolicExpr pivot = matrix[pivot_row][col];
        
        auto inv_pivot_node = std::make_shared<PowerNode>(pivot.root, SymbolicFactory::create_number(BigInt(-1)));
        
        for (size_t j = col; j <= num_vars; ++j) {
            
            std::vector<std::shared_ptr<SymbolicNode>> ops;
            ops.push_back(matrix[pivot_row][j].root);
            ops.push_back(inv_pivot_node); 
            auto mult = SymbolicFactory::create_multiply(ops);
            
            matrix[pivot_row][j] = SymbolicExpr(to_ptr(SymbolicExpr(mult))->simplify()->root);
        }
        
        
        for (size_t i = 0; i < num_eqs; ++i) {
            if (i != pivot_row) {
                SymbolicExpr fac = matrix[i][col];
                if (!fac.root->is_zero()) {
                    for (size_t j = col; j <= num_vars; ++j) {
                        
                        std::vector<std::shared_ptr<SymbolicNode>> mops;
                        mops.push_back(fac.root);
                        mops.push_back(matrix[pivot_row][j].root);
                        auto term = SymbolicFactory::create_multiply(mops);
                        
                        
                        std::vector<std::shared_ptr<SymbolicNode>> nops;
                        nops.push_back(SymbolicFactory::create_number(BigInt(-1)));
                        nops.push_back(term);
                        auto neg_term = SymbolicFactory::create_multiply(nops);
                        
                        std::vector<std::shared_ptr<SymbolicNode>> aops;
                        aops.push_back(matrix[i][j].root);
                        aops.push_back(neg_term);
                        auto res = SymbolicFactory::create_add(aops);
                        matrix[i][j] = SymbolicExpr(to_ptr(SymbolicExpr(res))->simplify()->root);
                    }
                }
            }
        }
        
        pivot_row++;
    }
    
    
    
    std::map<std::string, SymbolicExpr> solution;
    for (size_t i = 0; i < num_eqs; ++i) {
        size_t pcol = pivot_col_for_row[i];
        if (pcol != (size_t)-1) {
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            SymbolicExpr val = matrix[i][num_vars];
            for (size_t j = pcol + 1; j < num_vars; ++j) {
                SymbolicExpr c = matrix[i][j];
                if (!c.root->is_zero()) {
                    
                    std::vector<std::shared_ptr<SymbolicNode>> mops;
                    mops.push_back(c.root);
                    mops.push_back(SymbolicFactory::create_variable(variables[j]));
                    auto term = SymbolicFactory::create_multiply(mops);
                    
                    std::vector<std::shared_ptr<SymbolicNode>> nops;
                    nops.push_back(SymbolicFactory::create_number(BigInt(-1)));
                    nops.push_back(term);
                    
                    std::vector<std::shared_ptr<SymbolicNode>> aops;
                    aops.push_back(val.root);
                    aops.push_back(SymbolicFactory::create_multiply(nops));
                    val = SymbolicExpr(SymbolicFactory::create_add(aops));
                }
            }
            solution[variables[pcol]] = val;
        }
    }
    
    return solution;
}





namespace {

    using Monomial = std::vector<int>;

    
    struct MonomialLess {
        bool operator()(const Monomial& a, const Monomial& b) const {
            
            
            
            
            
            
            
            
            
            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i] != b[i]) {
                    return a[i] > b[i]; 
                }
            }
            return false;
        }
    };

    struct Term {
        Rational coeff;
        Monomial mono;
    };
    
    
    
    
    using PolyTerms = std::map<Monomial, Rational, MonomialLess>;

    struct Poly {
        PolyTerms terms;
        size_t num_vars;

        Poly() : num_vars(0) {}
        Poly(size_t n) : num_vars(n) {}

        bool is_zero() const { return terms.empty(); }

        Monomial LM() const {
            if (terms.empty()) return std::vector<int>(num_vars, 0); 
            return terms.begin()->first;
        }

        Rational LC() const {
            if (terms.empty()) return Rational(0);
            return terms.begin()->second;
        }
        
        Term LT() const {
             if (terms.empty()) return {Rational(0), std::vector<int>(num_vars, 0)};
             return {terms.begin()->second, terms.begin()->first};
        }
        
        void add_term(const Monomial& m, const Rational& c) {
            if (c == Rational(0)) return;
            auto it = terms.find(m);
            if (it != terms.end()) {
                
                
                
                Rational new_c = it->second + c; 
                if (new_c == Rational(0)) {
                    terms.erase(it);
                } else {
                    it->second = new_c;
                }
            } else {
                terms[m] = c;
            }
        }
    };

    Monomial mul_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] + b[i];
        return res;
    }
    
    Monomial lcm_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
        return res;
    }
    
    
    bool divides_mono(const Monomial& a, const Monomial& b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return false;
        }
        return true;
    }
    
    Monomial div_mono(const Monomial& num, const Monomial& den) {
        Monomial res(num.size());
        for (size_t i = 0; i < num.size(); ++i) res[i] = num[i] - den[i];
        return res;
    }
    
    Poly add_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {
            res.add_term(m, c);
        }
        return res;
    }

    Poly sub_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {
            
            Rational c_neg = c * Rational(-1);
            res.add_term(m, c_neg);
        }
        return res;
    }
    
    Poly mul_poly_term(const Poly& p, const Term& t) {
        Poly res(p.num_vars);
        if (t.coeff == Rational(0)) return res;
        for (auto const& [m, c] : p.terms) {
            res.add_term(mul_mono(m, t.mono), c * t.coeff);
        }
        return res;
    }
    
    Poly mul_poly(const Poly& a, const Poly& b) {
        Poly res(a.num_vars);
        for (auto const& [ma, ca] : a.terms) {
            for (auto const& [mb, cb] : b.terms) {
                res.add_term(mul_mono(ma, mb), ca * cb);
            }
        }
        return res;
    }
    
    
    class PolyBuilder : public SymbolicVisitor {
        std::vector<std::string> vars;
        Poly result;
    public:
        PolyBuilder(const std::vector<std::string>& v) : vars(v), result(v.size()) {}
        
        Poly get_result() const { return result; }
        
        void visit(NumberNode& node) override {
            result = Poly(vars.size());
            
            if (std::holds_alternative<Rational>(node.value)) {
                result.add_term(std::vector<int>(vars.size(), 0), std::get<Rational>(node.value));
            } else if (std::holds_alternative<BigInt>(node.value)) {
                result.add_term(std::vector<int>(vars.size(), 0), Rational(std::get<BigInt>(node.value)));
            } else if (std::holds_alternative<double>(node.value)) {
                
                result.add_term(std::vector<int>(vars.size(), 0), Rational((long long)std::get<double>(node.value)));
            } 
        }
        
        void visit(VariableNode& node) override {
            result = Poly(vars.size());
            auto it = std::find(vars.begin(), vars.end(), node.name);
            if (it != vars.end()) {
                Monomial m(vars.size(), 0);
                m[std::distance(vars.begin(), it)] = 1;
                result.add_term(m, Rational(1));
            } else {
                
                
                
                
                
            }
        }
        
        void visit(AddNode& node) override {
            Poly sum(vars.size());
            for (auto& op : node.operands) {
                PolyBuilder b(vars);
                op->accept(b);
                sum = add_poly(sum, b.get_result());
            }
            result = sum;
        }
        
        void visit(MultiplyNode& node) override {
            Poly prod(vars.size());
            prod.add_term(std::vector<int>(vars.size(), 0), Rational(1)); 
            
            for (auto& op : node.operands) {
                PolyBuilder b(vars);
                op->accept(b);
                prod = mul_poly(prod, b.get_result());
            }
            result = prod;
        }
        
        void visit(PowerNode& node) override {
            PolyBuilder b_base(vars);
            node.base->accept(b_base);
            Poly base = b_base.get_result();
            
            
            
            
            long long exp = 0;
            if (node.exponent->is_number()) {
                 
                 auto num_node = std::dynamic_pointer_cast<NumberNode>(node.exponent);
                 if (num_node) {
                     const auto& val = num_node->value;
                     if (std::holds_alternative<BigInt>(val)) exp = std::get<BigInt>(val).to_int();
                     else if (std::holds_alternative<Rational>(val)) exp = (long long)std::get<Rational>(val).to_double(); 
                     else exp = (long long)std::get<double>(val);
                 }
            }
            
            if (exp == 0) {
                result = Poly(vars.size());
                result.add_term(std::vector<int>(vars.size(), 0), Rational(1));
            } else {
                Poly res(vars.size());
                res.add_term(std::vector<int>(vars.size(), 0), Rational(1));
                for (int i = 0; i < exp; ++i) {
                    res = mul_poly(res, base);
                }
                result = res;
            }
        }
        
        void visit(FunctionNode& node) override { result = Poly(vars.size()); } 
        void visit(MatrixNode& node) override { result = Poly(vars.size()); }
        void visit(RelationalNode& node) override { result = Poly(vars.size()); }
    };

    Poly to_poly(const SymbolicExpr& expr, const std::vector<std::string>& vars) {
        PolyBuilder b(vars);
        expr.root->accept(b);
        return b.get_result();
    }
    
    SymbolicExpr from_poly(const Poly& p, const std::vector<std::string>& vars) {
        if (p.terms.empty()) return SymbolicExpr(SymbolicFactory::create_number(BigInt(0)));
        
        std::vector<std::shared_ptr<SymbolicNode>> add_ops;
        
        for (auto const& [m, c] : p.terms) {
            
            std::vector<std::shared_ptr<SymbolicNode>> mul_ops;
            
            
            if (c.get_denominator() == BigInt(1)) {
                mul_ops.push_back(SymbolicFactory::create_number(c.get_numerator()));
            } else {
                mul_ops.push_back(SymbolicFactory::create_number(c));
            }
            
            for (size_t i = 0; i < m.size(); ++i) {
                if (m[i] > 0) {
                    auto var = SymbolicFactory::create_variable(vars[i]);
                    if (m[i] == 1) {
                         mul_ops.push_back(var);
                    } else {
                         auto pow = SymbolicFactory::create_power(var, SymbolicFactory::create_number(BigInt(m[i])));
                         mul_ops.push_back(pow);
                    }
                }
            }
            
            if (mul_ops.size() == 0) { 
                
            } else if (mul_ops.size() == 1) {
                add_ops.push_back(mul_ops[0]);
            } else {
                add_ops.push_back(SymbolicFactory::create_multiply(mul_ops));
            }
        }
        
        if (add_ops.empty()) return SymbolicExpr(SymbolicFactory::create_number(BigInt(0)));
        if (add_ops.size() == 1) return SymbolicExpr(add_ops[0]);
        return SymbolicExpr(SymbolicFactory::create_add(add_ops));
    }
    
    
    Poly reduce(Poly p, const std::vector<Poly>& G) {
        Poly r(p.num_vars);
        
        while (!p.is_zero()) {
            bool reduced = false;
            Monomial p_lm = p.LM();
            
            for (const auto& g : G) {
                if (divides_mono(g.LM(), p_lm)) {
                    
                    Term factor;
                    factor.coeff = p.LC() / g.LC();
                    factor.mono = div_mono(p_lm, g.LM());
                    
                    
                    Poly term_g = mul_poly_term(g, factor);
                    p = sub_poly(p, term_g);
                    reduced = true;
                    break;
                }
            }
            
            if (!reduced) {
                Term lt = p.LT();
                r.add_term(lt.mono, lt.coeff);
                
                
                
                p.terms.erase(p.terms.begin());
            }
        }
        return r;
    }
    
    Poly s_poly(const Poly& f, const Poly& g) {
        Monomial L = lcm_mono(f.LM(), g.LM());
        
        
        
        Term t1;
        t1.mono = div_mono(L, f.LM());
        t1.coeff = Rational(1) / f.LC();
        
        Term t2;
        t2.mono = div_mono(L, g.LM());
        t2.coeff = Rational(1) / g.LC();
        
        return sub_poly(mul_poly_term(f, t1), mul_poly_term(g, t2));
    }

} 


std::vector<SymbolicExpr> Solver::groebner_basis(
    const std::vector<SymbolicExpr>& polynomials,
    const std::vector<std::string>& variables) 
{
    std::vector<Poly> G;
    for (const auto& expr : polynomials) {
        Poly p = to_poly(expr, variables);
        if (!p.is_zero()) G.push_back(p);
    }
    
    
    
    
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < G.size(); ++i) {
        for (size_t j = i + 1; j < G.size(); ++j) {
            pairs.push_back({i, j});
        }
    }
    
    while (!pairs.empty()) {
        auto [i, j] = pairs.back();
        pairs.pop_back();
        
        
        
        
        Poly S = s_poly(G[i], G[j]);
        Poly r = reduce(S, G);
        
        if (!r.is_zero()) {
            size_t new_idx = G.size();
            G.push_back(r);
            for (size_t k = 0; k < new_idx; ++k) {
                pairs.push_back({k, new_idx});
            }
        }
    }
    
    
    std::vector<SymbolicExpr> result;
    for (const auto& p : G) {
        result.push_back(from_poly(p, variables));
    }
    return result;
}

std::vector<std::map<std::string, SymbolicExpr>> Solver::solve_polynomial_system(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables) 
{
    
    
    auto G_basis = groebner_basis(equations, variables);
    
    return {};
}

} 
