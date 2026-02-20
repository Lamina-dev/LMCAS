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
        return SymbolicPolyCoeff(SymbolicExpr::add(val, other.val)->simplify());
    }
    
    SymbolicPolyCoeff operator-(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(
            SymbolicExpr::add(val, SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)))->simplify()
        );
    }
    
    SymbolicPolyCoeff operator*(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, other.val)->simplify());
    }

    SymbolicPolyCoeff operator/(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::divide(val, other.val)->simplify());
    }
    
    SymbolicPolyCoeff operator-() const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, SymbolicExpr::number(-1))->simplify());
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

template <typename T>
Polynomial<T> symbolic_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    
    auto expanded = expr->expand(); 
    
    std::map<int, T> coeff_map;
    
    
    auto process_term = [&](std::shared_ptr<SymbolicNode> node) {
        
        
        
        
        int degree = 0;
        
        
        
        T coeff_val = T(1);
        
        std::vector<std::shared_ptr<SymbolicNode>> factors;
        if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
            for(auto& op : mul->operands) factors.push_back(op);
        } else {
            factors.push_back(node);
        }
        
        for (const auto& f : factors) {
            
            
            bool is_var_term = false;
            
            if (auto v = std::dynamic_pointer_cast<VariableNode>(f)) {
                if (v->name == var) {
                    degree += 1;
                    is_var_term = true;
                }
            } else if (auto p = std::dynamic_pointer_cast<PowerNode>(f)) {
                
                if (auto v = std::dynamic_pointer_cast<VariableNode>(p->base)) {
                    if (v->name == var) {
                        
                        if (auto exp = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
                            
                             if (std::holds_alternative<BigInt>(exp->value)) 
                                 degree += std::get<BigInt>(exp->value).to_int(); 
                             else if (std::holds_alternative<double>(exp->value)) 
                                 degree += (int)std::get<double>(exp->value);
                             else if (std::holds_alternative<Rational>(exp->value))
                                 degree += (int)std::get<Rational>(exp->value).to_double();
                             is_var_term = true;
                        }
                    }
                }
            }
            
            if (!is_var_term) {
                
                
                
                
                
                
                auto node_expr = std::make_shared<SymbolicExpr>(f);
                T val = extract_coeff_value<T>(node_expr);
                coeff_val = coeff_val * val;
            }
        }
        
        
        if (coeff_map.find(degree) == coeff_map.end()) {
            coeff_map[degree] = coeff_val;
        } else {
            coeff_map[degree] = coeff_map[degree] + coeff_val;
        }
    };
    
    if (auto add = std::dynamic_pointer_cast<AddNode>(expanded->root)) {
        for (const auto& op : add->operands) process_term(op);
    } else {
        process_term(expanded->root);
    }
    
    
    if (coeff_map.empty()) return Polynomial<T>(var);
    
    int max_deg = coeff_map.rbegin()->first;
    if (max_deg < 0) max_deg = 0; 
    
    std::vector<T> coeffs(max_deg + 1, T(0));
    for (const auto& [deg, val] : coeff_map) {
        if (deg >= 0) coeffs[deg] = val;
    }
    
    return Polynomial<T>(coeffs, var);
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

} 
