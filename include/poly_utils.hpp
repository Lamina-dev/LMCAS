#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "symbolic_ast.hpp"

namespace lamina {






template <typename T>
Polynomial<T> symbolic_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);


template <typename T>
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const Polynomial<T>& poly);






struct SymbolicPolyCoeff {
    std::shared_ptr<SymbolicExpr> val;

    SymbolicPolyCoeff() : val(SymbolicExpr::number(0)) {}
    
    
    explicit SymbolicPolyCoeff(int v) : val(SymbolicExpr::number(v)) {}
    
    
    SymbolicPolyCoeff(std::shared_ptr<SymbolicExpr> v) : val(std::move(v)) {}

    bool operator==(const SymbolicPolyCoeff& other) const {
        
        
        if (!val || !other.val) return false;
        
        if (val == other.val) return true;
        
        
        auto diff = SymbolicExpr::add(val, SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)));
        return diff->simplify()->is_zero();
    }
    
    bool operator!=(const SymbolicPolyCoeff& other) const {
        return !(*this == other);
    }
    
    SymbolicPolyCoeff operator+(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(val, other.val));
    }
    
    SymbolicPolyCoeff operator-(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(
            SymbolicExpr::add(val, SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)))
        );
    }
    
    SymbolicPolyCoeff operator*(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, other.val));
    }

    SymbolicPolyCoeff operator/(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::divide(val, other.val));
    }
    
    SymbolicPolyCoeff operator-() const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, SymbolicExpr::number(-1)));
    }
    
    
    std::string ToString() const {
        return val ? val->to_string() : "0";
    }
    
    
    friend SymbolicPolyCoeff abs(const SymbolicPolyCoeff& s) {
        
        
        
        return s;
    }
};










template<typename T>
T extract_coeff_value(const std::shared_ptr<SymbolicExpr>& c);


template<>
inline SymbolicPolyCoeff extract_coeff_value<SymbolicPolyCoeff>(const std::shared_ptr<SymbolicExpr>& c) {
    return SymbolicPolyCoeff(c);
}

template<>
inline BigInt extract_coeff_value<BigInt>(const std::shared_ptr<SymbolicExpr>& c) {
    auto simp = c->simplify(); 
    if (auto n = std::dynamic_pointer_cast<NumberNode>(simp->root)) {
        if (std::holds_alternative<BigInt>(n->value)) return std::get<BigInt>(n->value);
        if (std::holds_alternative<Rational>(n->value)) {
            
            Rational r = std::get<Rational>(n->value);
            if (r.is_integer()) return r.to_BigInt();
        }
        if (std::holds_alternative<double>(n->value)) return BigInt((long long)std::get<double>(n->value));
    }
    
    if (simp->is_zero()) return BigInt(0);
    
    if (simp->is_one()) return BigInt(1);
    
    
    
    
    
    return BigInt(0); 
}

template<>
inline Rational extract_coeff_value<Rational>(const std::shared_ptr<SymbolicExpr>& c) {
    auto simp = c->simplify();
    if (auto n = std::dynamic_pointer_cast<NumberNode>(simp->root)) {
        if (std::holds_alternative<Rational>(n->value)) return std::get<Rational>(n->value);
        if (std::holds_alternative<BigInt>(n->value)) return Rational(std::get<BigInt>(n->value));
        if (std::holds_alternative<double>(n->value)) return Rational::from_double(std::get<double>(n->value));
    }
    if (simp->is_zero()) return Rational(0);
    if (simp->is_one()) return Rational(1);
    return Rational(0);
}

inline bool depends_on_var(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    if (!node) return false;
    
    struct DepVisitor : public SymbolicVisitor {
        std::string v;
        bool found = false;
        void visit(NumberNode&) override {}
        void visit(VariableNode& n) override { if (n.name == v) found = true; }
        void visit(AddNode& n) override { for(auto& op : n.operands) { if(found) return; op->accept(*this); } }
        void visit(MultiplyNode& n) override { for(auto& op : n.operands) { if(found) return; op->accept(*this); } }
        void visit(PowerNode& n) override { n.base->accept(*this); if(!found) n.exponent->accept(*this); }
        void visit(FunctionNode& n) override { for(auto& arg : n.arguments) { if(found) return; arg->accept(*this); } }
        void visit(MatrixNode& n) override {
            if (std::holds_alternative<MatrixNode::DenseStorage>(n.storage)) {
                for (auto& item : std::get<MatrixNode::DenseStorage>(n.storage)) {
                    if (item) { item->accept(*this); if (found) return; }
                }
            } else {
                for (auto& [idx, item] : std::get<MatrixNode::SparseStorage>(n.storage)) {
                    item->accept(*this); if (found) return;
                }
            }
        }
    } visitor;
    visitor.v = var;
    node->accept(visitor);
    return visitor.found;
}

inline bool contains(const SymbolicExpr& expr, const std::string& var) {
    return depends_on_var(expr.root, var);
}

template <typename T>
Polynomial<T> symbolic_to_poly_recursive(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    if (!node) return Polynomial<T>(var);

    // If it doesn't contain the variable, it's a coefficient (T)
    if (!depends_on_var(node, var)) {
        return Polynomial<T>({extract_coeff_value<T>(std::make_shared<SymbolicExpr>(node))}, var);
    }

    // Case: Variable
    if (auto v = std::dynamic_pointer_cast<VariableNode>(node)) {
        if (v->name == var) {
            return Polynomial<T>({T(0), T(1)}, var);
        }
    }

    // Case: Add
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        Polynomial<T> res(var);
        for (auto& op : add->operands) {
            res = res + symbolic_to_poly_recursive<T>(op, var);
        }
        return res;
    }

    // Case: Multiply
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        Polynomial<T> res({T(1)}, var);
        for (auto& op : mul->operands) {
            res = res * symbolic_to_poly_recursive<T>(op, var);
        }
        return res;
    }

    // Case: Power
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            int e_val = 0;
            if (std::holds_alternative<BigInt>(exp_num->value)) e_val = std::get<BigInt>(exp_num->value).to_int();
            else if (std::holds_alternative<double>(exp_num->value)) e_val = (int)std::get<double>(exp_num->value);
            else if (std::holds_alternative<Rational>(exp_num->value)) e_val = (int)std::get<Rational>(exp_num->value).to_double();
            
            if (e_val == 0) return Polynomial<T>({T(1)}, var);
            if (e_val > 0) {
                auto base_poly = symbolic_to_poly_recursive<T>(pow->base, var);
                Polynomial<T> res({T(1)}, var);
                for (int i = 0; i < e_val; ++i) res = res * base_poly;
                return res;
            }
        }
    }

    // Fallback or non-polynomial term
    return Polynomial<T>(var);
}

template <typename T>
Polynomial<T> symbolic_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !expr->root) return Polynomial<T>(var);
    return symbolic_to_poly_recursive<T>(expr->root, var);
}


template <typename T>
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const Polynomial<T>& poly) {
    if (poly.is_zero()) return SymbolicExpr::number(0);
    
    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        if (poly.coeffs[i] == T(0)) continue;
        
        auto coeff_node = SymbolicExpr::number(0);
        
        
        if constexpr (std::is_same_v<T, BigInt>) {
            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        } else if constexpr (std::is_same_v<T, Rational>) {
            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        } else {
            
            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        }
        
        if (i == 0) {
            terms.push_back(coeff_node);
        } else {
            
            auto var_node = std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(poly.variable_name));
            std::shared_ptr<SymbolicExpr> var_part;
            if (i == 1) {
                var_part = var_node;
            } else {
                var_part = SymbolicExpr::power(var_node, SymbolicExpr::number((int)i));
            }
            
            
            if (poly.coeffs[i] == T(1)) {
                terms.push_back(var_part);
            } else if (poly.coeffs[i] == T(-1)) {
                 terms.push_back(SymbolicExpr::multiply(SymbolicExpr::number(-1), var_part));
            } else {
                terms.push_back(SymbolicExpr::multiply(coeff_node, var_part));
            }
        }
    }
    
    
    
    std::reverse(terms.begin(), terms.end());
    
    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];
    
    
    
    auto res = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) {
        res = SymbolicExpr::add(res, terms[i]);
    }
    return res;
}

void gaussian_eliminate(std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& A, size_t m, size_t n, std::vector<size_t>& pivot_col_for_row, int& sign);

} 
