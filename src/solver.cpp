// Force rebuild for SymbolicFactory::create_add update
#include "solver.hpp"
#include "poly_utils.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <optional>
#include <limits>

namespace lamina {

void gaussian_eliminate(std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& A, size_t m, size_t n, std::vector<size_t>& pivot_col_for_row, int& sign) {
    size_t pivot_row = 0;
    sign = 1;
    pivot_col_for_row.assign(m, (size_t)-1);

    for (size_t col = 0; col < n && pivot_row < m; ++col) {
        size_t max_row = pivot_row;
        bool found_pivot = false;
        while (max_row < m) {
            A[max_row][col] = A[max_row][col]->simplify();
            if (!A[max_row][col]->is_zero()) {
                found_pivot = true;
                break;
            }
            max_row++;
        }
        
        if (!found_pivot) continue;
        
        if (pivot_row != max_row) {
            std::swap(A[pivot_row], A[max_row]);
            sign = -sign;
        }
        
        pivot_col_for_row[pivot_row] = col;
        auto pivot = A[pivot_row][col];
        auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));
        
        for (size_t k = col; k < A[pivot_row].size(); ++k) {
            A[pivot_row][k] = SymbolicExpr::multiply(A[pivot_row][k], pivot_inv)->simplify();
        }
        
        for (size_t i = 0; i < m; ++i) {
            if (i != pivot_row) {
                auto factor = A[i][col];
                if (!factor->is_zero()) {
                    auto neg_factor = SymbolicExpr::multiply(factor, SymbolicExpr::number(-1));
                    for (size_t k = col; k < A[i].size(); ++k) {
                        auto term = SymbolicExpr::multiply(neg_factor, A[pivot_row][k]);
                        A[i][k] = SymbolicExpr::add(A[i][k], term)->simplify();
                    }
                }
            }
        }
        pivot_row++;
    }
}

static std::shared_ptr<SymbolicExpr> to_ptr(const SymbolicExpr& expr) {
    return std::make_shared<SymbolicExpr>(expr);
}

static bool get_integer_value(const std::shared_ptr<SymbolicNode>& node, long long& value) {
    auto num = std::dynamic_pointer_cast<NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value)) {
        value = std::get<BigInt>(num->value).to_int();
        return true;
    }
    if (std::holds_alternative<Rational>(num->value)) {
        const auto& r = std::get<Rational>(num->value);
        if (!r.is_integer()) return false;
        value = r.to_BigInt().to_int();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value)) {
        lmmc_real_t d = std::get<lmmc_real_t>(num->value);
        int eq;
        lmmc_double_nearly_equal_tol(d, std::round(d), 1e-12, 1e-12, &eq);
        if (!eq) return false;
        value = static_cast<long long>(std::llround(d));
        return true;
    }
    return false;
}

static std::shared_ptr<NumberNode> add_number_nodes(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
    if (std::holds_alternative<lmmc_real_t>(a->value) || std::holds_alternative<lmmc_real_t>(b->value)) {
        auto to_real = [](const auto& v) {
            if (std::holds_alternative<lmmc_real_t>(v)) return std::get<lmmc_real_t>(v);
            if (std::holds_alternative<Rational>(v)) return (lmmc_real_t)std::get<Rational>(v).to_double();
            return (lmmc_real_t)std::get<BigInt>(v).to_double();
        };
        lmmc_real_t r1 = to_real(a->value);
        lmmc_real_t r2 = to_real(b->value);
        return std::make_shared<NumberNode>(r1 + r2);
    }

    if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
        Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) : Rational(std::get<BigInt>(a->value));
        Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) : Rational(std::get<BigInt>(b->value));
        return std::make_shared<NumberNode>(r1 + r2);
    }

    return std::make_shared<NumberNode>(std::get<BigInt>(a->value) + std::get<BigInt>(b->value));
}

static std::shared_ptr<NumberNode> multiply_number_nodes(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
    if (std::holds_alternative<lmmc_real_t>(a->value) || std::holds_alternative<lmmc_real_t>(b->value)) {
        auto to_real = [](const auto& v) {
            if (std::holds_alternative<lmmc_real_t>(v)) return std::get<lmmc_real_t>(v);
            if (std::holds_alternative<Rational>(v)) return (lmmc_real_t)std::get<Rational>(v).to_double();
            return (lmmc_real_t)std::get<BigInt>(v).to_double();
        };
        lmmc_real_t r1 = to_real(a->value);
        lmmc_real_t r2 = to_real(b->value);
        return std::make_shared<NumberNode>(r1 * r2);
    }

    if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
        Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) : Rational(std::get<BigInt>(a->value));
        Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) : Rational(std::get<BigInt>(b->value));
        return std::make_shared<NumberNode>(r1 * r2);
    }

    return std::make_shared<NumberNode>(std::get<BigInt>(a->value) * std::get<BigInt>(b->value));
}

struct NodeLess {
    bool operator()(const std::shared_ptr<SymbolicNode>& a, const std::shared_ptr<SymbolicNode>& b) const {
        if (!a || !b) return a < b;
        return a->compare(*b) < 0;
    }
};

static std::shared_ptr<SymbolicExpr> multiply_no_expand(
    const std::shared_ptr<SymbolicNode>& term,
    const std::vector<std::shared_ptr<SymbolicNode>>& den_factors
) {
    std::vector<std::shared_ptr<SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(term)) {
        factors.insert(factors.end(), mul->operands.begin(), mul->operands.end());
    } else if (term) {
        factors.push_back(term);
    }
    factors.insert(factors.end(), den_factors.begin(), den_factors.end());

    auto const_acc = std::make_shared<NumberNode>(BigInt(1));
    std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeLess> bases;

    for (const auto& op : factors) {
        if (!op) continue;
        if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
            const_acc = multiply_number_nodes(const_acc, num);
            continue;
        }

        std::shared_ptr<SymbolicNode> base = op;
        std::shared_ptr<NumberNode> exp = std::make_shared<NumberNode>(BigInt(1));
        if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
            base = pow->base;
            if (auto e_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                exp = e_num;
            }
        }

        auto it = bases.find(base);
        if (it == bases.end()) {
            bases[base] = exp;
        } else {
            bases[base] = add_number_nodes(it->second, exp);
        }
    }

    std::vector<std::shared_ptr<SymbolicNode>> final_ops;
    if (!const_acc->is_one()) final_ops.push_back(const_acc);

    for (const auto& [base, exp] : bases) {
        if (exp->is_zero()) continue;
        if (exp->is_one()) final_ops.push_back(base);
        else final_ops.push_back(std::make_shared<PowerNode>(base, exp));
    }

    if (final_ops.empty()) return SymbolicExpr::number(1);
    if (final_ops.size() == 1) return std::make_shared<SymbolicExpr>(final_ops[0]);
    return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(final_ops));
}

static bool is_polynomial_node(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (std::dynamic_pointer_cast<NumberNode>(node)) return true;
    if (std::dynamic_pointer_cast<VariableNode>(node)) return true;

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (const auto& op : add->operands) {
            if (!is_polynomial_node(op)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (const auto& op : mul->operands) {
            if (!is_polynomial_node(op)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        long long exp = 0;
        if (!get_integer_value(pow->exponent, exp)) return false;
        if (exp < 0) return false;
        return is_polynomial_node(pow->base);
    }

    return false;
}

static std::shared_ptr<SymbolicExpr> multiply_factors(const std::vector<std::shared_ptr<SymbolicNode>>& factors) {
    if (factors.empty()) return SymbolicExpr::number(1);
    auto res = std::make_shared<SymbolicExpr>(factors[0]);
    for (size_t i = 1; i < factors.size(); ++i) {
        res = SymbolicExpr::multiply(res, std::make_shared<SymbolicExpr>(factors[i]));
    }
    return res->simplify();
}

static bool collect_denominator_factors(
    const std::shared_ptr<SymbolicNode>& node,
    std::vector<std::shared_ptr<SymbolicNode>>& den_factors,
    std::vector<std::shared_ptr<SymbolicExpr>>& den_constraints
) {
    if (!node) return false;

    if (std::dynamic_pointer_cast<NumberNode>(node) || std::dynamic_pointer_cast<VariableNode>(node)) {
        return true;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (const auto& op : add->operands) {
            if (!collect_denominator_factors(op, den_factors, den_constraints)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (const auto& op : mul->operands) {
            if (!collect_denominator_factors(op, den_factors, den_constraints)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        long long exp = 0;
        if (!get_integer_value(pow->exponent, exp)) return false;

        if (exp < 0) {
            if (!is_polynomial_node(pow->base)) return false;
            long long k = -exp;
            if (k == 1) {
                den_factors.push_back(pow->base);
            } else {
                den_factors.push_back(SymbolicFactory::create_power(pow->base, SymbolicFactory::create_number(BigInt(k))));
            }
            den_constraints.push_back(std::make_shared<SymbolicExpr>(pow->base));
            return true;
        }

        if (!collect_denominator_factors(pow->base, den_factors, den_constraints)) return false;
        return true;
    }

    // FunctionNode: traverse into arguments (functions don't introduce denominators directly)
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        for (const auto& arg : func->arguments) {
            if (!collect_denominator_factors(arg, den_factors, den_constraints)) return false;
        }
        return true;
    }

    return false;
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
    
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix(num_eqs, std::vector<std::shared_ptr<SymbolicExpr>>(num_vars + 1));
    
    for (size_t i = 0; i < num_eqs; ++i) {
        SymbolicExpr current_eq = equations_in[i]; 
        
        SymbolicExpr constant_part = current_eq;
        
        for (size_t j = 0; j < num_vars; ++j) {
            auto [coeff, remainder] = isolate_linear_coeff(constant_part, variables[j]);
            matrix[i][j] = std::make_shared<SymbolicExpr>(coeff);
            
            constant_part = remainder;
        }
        
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
        ops.push_back(constant_part.root);
        auto neg_const = SymbolicFactory::create_multiply(ops);
        matrix[i][num_vars] = std::make_shared<SymbolicExpr>(neg_const);
    }
    
    std::vector<size_t> pivot_col_for_row;
    int sign;
    gaussian_eliminate(matrix, num_eqs, num_vars, pivot_col_for_row, sign);
    
    std::map<std::string, SymbolicExpr> solution;
    for (size_t i = 0; i < num_eqs; ++i) {
        size_t pcol = pivot_col_for_row[i];
        if (pcol != (size_t)-1) {
            SymbolicExpr val = *matrix[i][num_vars];
            for (size_t j = pcol + 1; j < num_vars; ++j) {
                SymbolicExpr c = *matrix[i][j];
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
            } else if (std::holds_alternative<lmmc_real_t>(node.value)) {
                result.add_term(std::vector<int>(vars.size(), 0), Rational((long long)std::get<lmmc_real_t>(node.value)));
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
                      else exp = (long long)std::get<lmmc_real_t>(val);
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
    std::vector<SymbolicExpr> cleared_equations;
    std::vector<std::shared_ptr<SymbolicExpr>> denom_constraints;
    cleared_equations.reserve(equations.size());

    for (const auto& eq : equations) {
        if (!eq.root) return {};

        std::vector<std::shared_ptr<SymbolicNode>> den_factors;
        std::vector<std::shared_ptr<SymbolicExpr>> den_local_constraints;
        if (!collect_denominator_factors(eq.root, den_factors, den_local_constraints)) return {};

        auto denom_expr = multiply_factors(den_factors);
        auto cleared = to_ptr(eq);
        if (!den_factors.empty()) {
            if (auto add = std::dynamic_pointer_cast<AddNode>(cleared->root)) {
                std::vector<std::shared_ptr<SymbolicNode>> new_ops;
                new_ops.reserve(add->operands.size());
                for (const auto& op : add->operands) {
                    auto prod = multiply_no_expand(op, den_factors);
                    new_ops.push_back(prod->root);
                }
                cleared = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(new_ops));
            } else {
                cleared = multiply_no_expand(cleared->root, den_factors);
            }
        } else {
            cleared = cleared->simplify();
        }

        if (!cleared || !cleared->root || !is_polynomial_node(cleared->root)) return {};
        cleared_equations.push_back(*cleared);

        for (const auto& c : den_local_constraints) {
            denom_constraints.push_back(c);
        }
    }

    if (cleared_equations.size() == 1 && variables.size() == 1) {
        auto roots = SymbolicExpr::solve(std::make_shared<SymbolicExpr>(cleared_equations[0]), variables[0]);
        std::vector<std::map<std::string, SymbolicExpr>> single_solutions;
        for (const auto& r : roots) {
            single_solutions.push_back({{variables[0], *r}});
        }
        if (denom_constraints.empty()) return single_solutions;

        std::vector<std::map<std::string, SymbolicExpr>> filtered;
        filtered.reserve(single_solutions.size());
        for (const auto& sol : single_solutions) {
            bool ok = true;
            for (const auto& den : denom_constraints) {
                auto sub = den;
                for (const auto& [name, val] : sol) {
                    sub = sub->substitute(name, std::make_shared<SymbolicExpr>(val));
                    if (!sub) break;
                }
                if (!sub) continue;
                sub = sub->simplify();
                if (sub && sub->is_zero()) {
                    ok = false;
                    break;
                }
            }
            if (ok) filtered.push_back(sol);
        }
        return filtered;
    }

    auto G_basis = groebner_basis(cleared_equations, variables);

    std::vector<std::shared_ptr<SymbolicExpr>> basis;
    basis.reserve(G_basis.size());
    for (const auto& g : G_basis) {
        auto g_ptr = std::make_shared<SymbolicExpr>(g);
        auto simp = g_ptr->simplify();
        if (simp && !simp->is_zero()) {
            basis.push_back(simp);
        }
    }

    if (variables.empty()) {
        for (const auto& p : basis) {
            if (p->is_number() && !p->is_zero()) return {};
        }
        return {{} };
    }

    auto substitute_all = [&](const std::shared_ptr<SymbolicExpr>& expr,
                              const std::map<std::string, SymbolicExpr>& subs) {
        auto res = expr;
        for (const auto& [name, val] : subs) {
            res = res->substitute(name, std::make_shared<SymbolicExpr>(val));
            if (!res) return std::shared_ptr<SymbolicExpr>(nullptr);
        }
        return res->simplify();
    };

    auto solve_rec = [&](auto&& self, int var_pos,
                         const std::map<std::string, SymbolicExpr>& partial)
        -> std::vector<std::map<std::string, SymbolicExpr>> {
        std::vector<std::shared_ptr<SymbolicExpr>> reduced;
        reduced.reserve(basis.size());

        for (const auto& p : basis) {
            auto r = substitute_all(p, partial);
            if (!r) continue;
            if (r->is_zero()) continue;

            bool depends = false;
            for (int i = 0; i <= var_pos && i < (int)variables.size(); ++i) {
                if (contains(*r, variables[i])) {
                    depends = true;
                    break;
                }
            }

            if (!depends) {
                if (r->is_number() && !r->is_zero()) return std::vector<std::map<std::string, SymbolicExpr>>{};
                continue;
            }

            reduced.push_back(r);
        }

        if (var_pos < 0) {
            return {partial};
        }

        const auto& curr_var = variables[var_pos];
        bool curr_var_appears = false;
        std::shared_ptr<SymbolicExpr> target = nullptr;
        int best_deg = std::numeric_limits<int>::max();

        for (const auto& r : reduced) {
            if (!contains(*r, curr_var)) continue;
            curr_var_appears = true;

            bool has_other = false;
            for (int i = 0; i < var_pos; ++i) {
                if (contains(*r, variables[i])) {
                    has_other = true;
                    break;
                }
            }
            if (has_other) continue;

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(r, curr_var);
            int deg = poly.degree();
            if (deg >= 1 && deg < best_deg) {
                best_deg = deg;
                target = r;
            }
        }

        if (!curr_var_appears) {
            auto next_partial = partial;
            next_partial[curr_var] = *SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        if (!target) return {};

        auto roots = SymbolicExpr::solve(target, curr_var);
        if (roots.empty()) return {};

        std::vector<std::map<std::string, SymbolicExpr>> results;
        for (const auto& r : roots) {
            auto next_partial = partial;
            next_partial[curr_var] = *r;
            auto sub_res = self(self, var_pos - 1, next_partial);
            results.insert(results.end(), sub_res.begin(), sub_res.end());
        }
        return results;
    };

    std::map<std::string, SymbolicExpr> empty;
    auto candidates = solve_rec(solve_rec, static_cast<int>(variables.size()) - 1, empty);

    if (denom_constraints.empty()) return candidates;

    std::vector<std::map<std::string, SymbolicExpr>> filtered;
    filtered.reserve(candidates.size());

    for (const auto& sol : candidates) {
        bool ok = true;
        for (const auto& den : denom_constraints) {
            auto sub = den;
            for (const auto& [name, val] : sol) {
                sub = sub->substitute(name, std::make_shared<SymbolicExpr>(val));
                if (!sub) break;
            }
            if (!sub) continue;
            sub = sub->simplify();
            if (sub && sub->is_zero()) {
                ok = false;
                break;
            }
        }
        if (ok) filtered.push_back(sol);
    }

    return filtered;
}

} 
