#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <numeric> 
#include "bigint.hpp"
#include "rational.hpp"

namespace lamina {


template<typename T>
T gcd_coeff_impl(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        return std::abs(std::gcd(a, b));
    } 
    else if constexpr (std::is_same_v<T, BigInt>) {
        return BigInt::gcd(a, b);
    }
    else if constexpr (std::is_same_v<T, Rational>) {
        return Rational(1); 
    }
    else {
        return T(1); 
    }
}

template <typename CoeffType>
class Polynomial {
public:
    
    
    std::vector<CoeffType> coeffs;
    std::string variable_name;

    Polynomial(const std::string& var = "x") : variable_name(var) {}

    Polynomial(const std::vector<CoeffType>& c, const std::string& var = "x") 
        : coeffs(c), variable_name(var) {
        trim();
    }
    
    Polynomial(const CoeffType& c, const std::string& var = "x") 
        : coeffs{c}, variable_name(var) {
        trim();
    }

    
    void trim() {
        while (coeffs.size() > 0 && coeffs.back() == CoeffType(0)) {
            coeffs.pop_back();
        }
    }

    bool is_zero() const {
        return coeffs.empty();
    }

    
    int degree() const {
        if (coeffs.empty()) return -1;
        return static_cast<int>(coeffs.size()) - 1;
    }

    CoeffType lead_coeff() const {
        if (coeffs.empty()) return CoeffType(0);
        return coeffs.back();
    }

    Polynomial operator+(const Polynomial& other) const {
        if (variable_name != other.variable_name) {
            
            
            
            
            if (!is_zero() && !other.is_zero()) {
                
            }
        }
        
        Polynomial res(variable_name);
        size_t n = std::max(coeffs.size(), other.coeffs.size());
        res.coeffs.resize(n);
        
        for (size_t i = 0; i < n; ++i) {
            CoeffType a = (i < coeffs.size()) ? coeffs[i] : CoeffType(0);
            CoeffType b = (i < other.coeffs.size()) ? other.coeffs[i] : CoeffType(0);
            res.coeffs[i] = a + b;
        }
        res.trim();
        return res;
    }

    Polynomial operator-(const Polynomial& other) const {
         Polynomial res(variable_name);
         size_t n = std::max(coeffs.size(), other.coeffs.size());
         res.coeffs.resize(n);
         
         for (size_t i = 0; i < n; ++i) {
             CoeffType a = (i < coeffs.size()) ? coeffs[i] : CoeffType(0);
             CoeffType b = (i < other.coeffs.size()) ? other.coeffs[i] : CoeffType(0);
             res.coeffs[i] = a - b;
         }
         res.trim();
         return res;
    }

    Polynomial operator*(const Polynomial& other) const {
        if (is_zero() || other.is_zero()) return Polynomial(variable_name);
        
        Polynomial res(variable_name);
        res.coeffs.resize(coeffs.size() + other.coeffs.size() - 1, CoeffType(0));
        
        for (size_t i = 0; i < coeffs.size(); ++i) {
            for (size_t j = 0; j < other.coeffs.size(); ++j) {
                res.coeffs[i + j] = res.coeffs[i + j] + coeffs[i] * other.coeffs[j];
            }
        }
        res.trim();
        return res;
    }
    
    bool operator==(const Polynomial& other) const {
         if (degree() != other.degree()) return false;
         if (variable_name != other.variable_name && !is_zero() && !other.is_zero()) return false;
         for (size_t i = 0; i < coeffs.size(); ++i) {
             if (coeffs[i] != other.coeffs[i]) return false;
         }
         return true;
    }

    
    
    std::pair<Polynomial, Polynomial> div_mod(const Polynomial& other) const {
        if (other.is_zero()) throw std::runtime_error("Division by zero polynomial");
        
        Polynomial quotient(variable_name);
        Polynomial remainder = *this;
        
        int deg_rem = remainder.degree();
        int deg_div = other.degree();
        CoeffType lc_div = other.lead_coeff();

        if (deg_rem < deg_div) {
             return {quotient, remainder};
        }

        
        quotient.coeffs.resize(deg_rem - deg_div + 1, CoeffType(0));

        while (deg_rem >= deg_div && !remainder.is_zero()) {
             
             CoeffType factor = remainder.lead_coeff() / lc_div;
             int deg_diff = deg_rem - deg_div;
             
             
             quotient.coeffs[deg_diff] = factor;
             
             
             
             Polynomial term(variable_name);
             term.coeffs.resize(deg_diff + 1, CoeffType(0));
             term.coeffs[deg_diff] = factor;
             
             Polynomial subtrahend = other * term;
             remainder = remainder - subtrahend; 
             deg_rem = remainder.degree();
        }
        
        quotient.trim();
        return {quotient, remainder};
    }
    

    
    static CoeffType gcd_coeff(const CoeffType& a, const CoeffType& b) {
        if constexpr (std::is_same_v<CoeffType, BigInt>) {
            return BigInt::gcd(a, b);
        } else if constexpr (std::is_same_v<CoeffType, Rational>) {
            return CoeffType(1);
        } else {
             
             if constexpr (std::is_integral_v<CoeffType>) {
                 return std::gcd(a, b);
             }
             
             return CoeffType(1);
        }
    }

    
    
    
    CoeffType content() const {
        if (is_zero()) return CoeffType(0);
        CoeffType g = coeffs[0];
        
        
        if constexpr (std::is_same_v<CoeffType, Rational>) {
             return CoeffType(1); 
        } else {
            
            
            
            
            for (size_t i = 1; i < coeffs.size(); ++i) {
                g = gcd_coeff(g, coeffs[i]);
                if (g == CoeffType(1)) break; 
            }
            return g; 
        }
    }

    
    Polynomial primitive_part() const {
        if (is_zero()) return *this;
        CoeffType c = content();
        if (c == CoeffType(0) || c == CoeffType(1)) return *this; 
        
        Polynomial res(variable_name);
        res.coeffs.reserve(coeffs.size());
        for (const auto& val : coeffs) {
            
            res.coeffs.push_back(val / c); 
        }
        
        
        if (res.lead_coeff() < CoeffType(0)) {
             
             for (auto& val : res.coeffs) {
                 
                 
                 
                 val = CoeffType(0) - val; 
             }
        }
        return res;
    }

    
    
    std::pair<Polynomial, Polynomial> pseudo_div_mod(const Polynomial& other) const {
        if (other.is_zero()) throw std::runtime_error("Division by zero polynomial");
        
        int deg_rem = degree();
        int deg_div = other.degree();
        
        if (deg_rem < deg_div) {
             return {Polynomial(variable_name), *this};
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        Polynomial remainder = *this;
        Polynomial quotient(variable_name); 
        CoeffType lc_div = other.lead_coeff();
        
        
        
        
        
        
        
        
        int delta = deg_rem - deg_div;
        quotient.coeffs.resize(delta + 1, CoeffType(0));
        
        while (deg_rem >= deg_div && !remainder.is_zero()) {
             int current_diff = deg_rem - deg_div;
             CoeffType lc_rem = remainder.lead_coeff();
             
             
             
             
             
             for (auto& c : quotient.coeffs) c = c * lc_div;
             
             quotient.coeffs[current_diff] = quotient.coeffs[current_diff] + lc_rem;

             
             
             
             
             
             
             
             
             Polynomial term(variable_name);
             term.coeffs.resize(current_diff + 1, CoeffType(0)); 
             term.coeffs[current_diff] = lc_rem; 
             
             
             
             
             
             for(auto& c : remainder.coeffs) c = c * lc_div;
             
             
             
             for (size_t i = 0; i < other.coeffs.size(); ++i) {
                 size_t target_idx = i + current_diff;
                 
                 if (target_idx < remainder.coeffs.size()) {
                    remainder.coeffs[target_idx] = remainder.coeffs[target_idx] - other.coeffs[i] * lc_rem;
                 }
                 
             }
             
             remainder.trim();
             deg_rem = remainder.degree();
        }
        
        return {quotient, remainder};
    }

    
    
    CoeffType eval(const CoeffType& val) const {
        if (coeffs.empty()) return CoeffType(0);
        CoeffType res = coeffs.back();
        for (int i = (int)coeffs.size() - 2; i >= 0; --i) {
            res = res * val + coeffs[i];
        }
        return res;
    }
    
    
    Polynomial differentiate() const {
        if (degree() < 1) return Polynomial(variable_name);
        Polynomial res(variable_name);
        res.coeffs.resize(coeffs.size() - 1);
        for (size_t i = 1; i < coeffs.size(); ++i) {
            CoeffType c = coeffs[i]; 
            
            
            
            res.coeffs[i-1] = c * CoeffType(i); 
        }
        res.trim();
        return res;
    }

    
    
    Polynomial make_monic() const {
        if (is_zero()) return *this;
        CoeffType lc = lead_coeff();
        if (lc == CoeffType(1)) return *this;
        
        Polynomial res(variable_name);
        res.coeffs.reserve(coeffs.size());
        for (const auto& c : coeffs) {
            res.coeffs.push_back(c / lc);
        }
        return res;
    }
    
    Polynomial pseudo_div_mod_rem(const Polynomial& other) const {
        
        return pseudo_div_mod(other).second;
    }
    
    
    
    
    static Polynomial gcd(Polynomial a, Polynomial b) {
        if (a.is_zero()) return b;
        if (b.is_zero()) return a;

        if constexpr (std::is_same_v<CoeffType, Rational>) {
            
            while (!b.is_zero()) {
                auto [q, r] = a.div_mod(b);
                a = b;
                b = r;
            }
            
            return a.make_monic();
        } else {
            
            
            CoeffType cA = a.content();
            CoeffType cB = b.content();
            CoeffType c  = gcd_coeff_impl(cA, cB);
            
            a = a.primitive_part();
            b = b.primitive_part();
            
            while (!b.is_zero()) {
                
                Polynomial r = a.pseudo_div_mod_rem(b);
                
                if (r.is_zero()) {
                    a = b;
                    b = r; 
                } else {
                    
                    a = b;
                    b = r.primitive_part(); 
                }
            }
            
            
            if (c == CoeffType(1)) return a;
            
            
            for (auto& val : a.coeffs) val = val * c;
            return a;
        }
    }

    
    Polynomial square_free_part() const {
        
        if (degree() <= 0) return *this;
        
        Polynomial deriv = differentiate();
        Polynomial g = gcd(*this, deriv);
        
        
        auto [q, r] = div_mod(g);
        
        
        
        if constexpr (std::is_same_v<CoeffType, Rational>) {
            return q.make_monic();
        }
        return q;
    }


    std::string to_string() const {
        if (is_zero()) return "0";
        std::string s = "";
        for (int i = degree(); i >= 0; --i) {
            CoeffType c = coeffs[i];
            if (c == CoeffType(0)) continue;
            
            bool positive = (c > CoeffType(0));
            if (!s.empty()) {
                s += (positive ? " + " : " - ");
                
                c = (c < CoeffType(0)) ? (c * CoeffType(-1)) : c;
            } else {
                if (!positive) {
                    s += "-";
                    c = c * CoeffType(-1);
                }
            }
            
            bool print_coeff = (c != CoeffType(1)) || (i == 0);
            if (print_coeff) {
                
                
                
                
                
                
                 
                 
                 
                 
            }

            
            
            
            
            
        }
        
        return "Poly(" + std::to_string(degree()) + ")";
    }
};


template<typename T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p) {
    if (p.is_zero()) return os << "0";
    for (int i = p.degree(); i >= 0; --i) {
        T c = p.coeffs[i];
        if (c == T(0)) continue;
        
        if (i < p.degree()) {
            if (c > T(0)) os << " + ";
            else os << " - ";
        } else {
            if (c < T(0)) os << "-";
        }
        
        T abs_c = (c < T(0)) ? (c * T(-1)) : c;
        if (abs_c != T(1) || i == 0) os << abs_c;
        
        if (i > 0) os << p.variable_name;
        if (i > 1) os << "^" << i;
    }
    return os;
}

} 
