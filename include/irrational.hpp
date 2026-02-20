#pragma once
#define _USE_MATH_DEFINES
#include "symbolic.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif


class Irrational {
public:
    

    std::shared_ptr<SymbolicExpr> to_symbolic() const {
        switch (type) {
            case Type::SQRT: {
                
                auto sqrtExpr = SymbolicExpr::sqrt(SymbolicExpr::number(static_cast<int>(radicand)));
                if (std::abs(coefficient) < 1e-15) {
                    return SymbolicExpr::number(0);
                } else if (std::abs(coefficient - 1.0) < 1e-15) {
                    return sqrtExpr;
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational(coefficient)), sqrtExpr);
                }
            }
            case Type::PI:
                
                if (std::abs(coefficient) < 1e-15) {
                    return SymbolicExpr::variable("π");
                } else if (std::abs(coefficient - 1.0) < 1e-15) {
                    return SymbolicExpr::variable("π");
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational(coefficient)), SymbolicExpr::variable("π"));
                }
            case Type::E:
                if (std::abs(coefficient) < 1e-15) {
                    return SymbolicExpr::number(0);
                } else if (std::abs(coefficient - 1.0) < 1e-15) {
                    return SymbolicExpr::variable("e");
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational(coefficient)), SymbolicExpr::variable("e"));
                }
            case Type::LOG:
                
                if (std::abs(coefficient) < 1e-15) {
                    return SymbolicExpr::number(0);
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational(coefficient)), SymbolicExpr::variable("log(" + std::to_string(radicand) + ")"));
                }
            case Type::COMPLEX:
                
                return SymbolicExpr::number(::Rational(constant_term));
            default:
                return SymbolicExpr::number(0);
        }
    }
    enum class Type {
        SQRT,  
        PI,    
        E,     
        LOG,   
        COMPLEX
    };

private:
    Type type;

    
    double coefficient;
    long long radicand;

    
    std::map<std::string, double> coefficients;
    double constant_term;

    
    static std::pair<long long, long long> simplify_sqrt(long long n) {
        long long perfect_square = 1;
        long long remainder = n;

        for (long long i = 2; i * i <= n; ++i) {
            while (remainder % (i * i) == 0) {
                perfect_square *= i;
                remainder /= (i * i);
            }
        }
        return {perfect_square, remainder};
    }

public:
    
    Irrational() : type(Type::COMPLEX), coefficient(0), radicand(1), constant_term(0) {}

    
    static Irrational sqrt(long long n, double coeff = 1.0) {
        Irrational result;
        result.type = Type::SQRT;

        auto [perfect, remainder] = simplify_sqrt(n);
        result.coefficient = coeff * perfect;
        result.radicand = remainder;
        result.constant_term = 0;

        return result;
    }

    
    static Irrational pi(double coeff = 1.0) {
        Irrational result;
        result.type = Type::PI;
        result.coefficient = coeff;
        result.radicand = 1;
        result.constant_term = 0;
        return result;
    }

    
    static Irrational e(double coeff = 1.0) {
        Irrational result;
        result.type = Type::E;
        result.coefficient = coeff;
        result.radicand = 1;
        result.constant_term = 0;
        return result;
    }

    
    static Irrational constant(double value) {
        Irrational result;
        result.type = Type::COMPLEX;
        result.coefficient = 0;
        result.radicand = 1;
        result.constant_term = value;
        return result;
    }

    
    void to_complex() {
        if (type == Type::COMPLEX) return;

        coefficients.clear();
        constant_term = 0;

        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {
                    constant_term = coefficient;
                } else {
                    coefficients["sqrt" + std::to_string(radicand)] = coefficient;
                }
                break;
            case Type::PI:
                coefficients["pi"] = coefficient;
                break;
            case Type::E:
                coefficients["e"] = coefficient;
                break;
            default:
                break;
        }
        type = Type::COMPLEX;
    }

    
    Irrational operator+(const Irrational& other) const {
        Irrational result = *this;
        Irrational other_copy = other;

        result.to_complex();
        other_copy.to_complex();

        result.constant_term += other_copy.constant_term;

        for (const auto& [key, coeff]: other_copy.coefficients) {
            result.coefficients[key] += coeff;
        }

        return result;
    }

    
    Irrational operator-(const Irrational& other) const {
        Irrational result = *this;
        Irrational other_copy = other;

        result.to_complex();
        other_copy.to_complex();

        result.constant_term -= other_copy.constant_term;

        for (const auto& [key, coeff]: other_copy.coefficients) {
            result.coefficients[key] -= coeff;
        }

        return result;
    }

    
    Irrational operator*(double scalar) const {
        Irrational result = *this;

        if (type == Type::COMPLEX) {
            result.constant_term *= scalar;
            for (auto& [key, coeff]: result.coefficients) {
                coeff *= scalar;
            }
        } else {
            result.coefficient *= scalar;
        }

        return result;
    }

    
    Irrational operator*(const Irrational& other) const {
        
        if (type == Type::COMPLEX && coefficients.empty()) {
            return other * constant_term;
        }
        if (other.type == Type::COMPLEX && other.coefficients.empty()) {
            return *this * other.constant_term;
        }

        
        if (type == Type::SQRT && other.type == Type::SQRT) {
            return Irrational::sqrt(radicand * other.radicand,
                                    coefficient * other.coefficient);
        }

        
        return Irrational::constant(to_double() * other.to_double());
    }

    
    Irrational operator/(const Irrational& other) const {
        
        if (other.type == Type::COMPLEX && other.coefficients.empty() && other.constant_term != 0) {
            return *this * (1.0 / other.constant_term);
        }

        
        double other_val = other.to_double();
        if (std::abs(other_val) < 1e-15) {
            throw std::runtime_error("Irrational: division by zero");
        }
        return Irrational::constant(to_double() / other_val);
    }

    
    Irrational operator-() const {
        return *this * (-1.0);
    }

    
    bool operator==(const Irrational& other) const {
        return std::abs(to_double() - other.to_double()) < 1e-12;
    }

    bool operator<(const Irrational& other) const {
        return to_double() < other.to_double();
    }

    bool operator<=(const Irrational& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const Irrational& other) const {
        return other < *this;
    }

    bool operator>=(const Irrational& other) const {
        return *this > other || *this == other;
    }

    
    double to_double() const {
        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {
                    return coefficient;
                }
                return coefficient * std::sqrt(radicand);
            case Type::PI:
                return coefficient * M_PI;
            case Type::E:
                return coefficient * M_E;
            case Type::LOG:
                return coefficient * std::log(radicand);
            case Type::COMPLEX: {
                double result = constant_term;
                for (const auto& [key, coeff]: coefficients) {
                    if (key == "pi") {
                        result += coeff * M_PI;
                    } else if (key == "e") {
                        result += coeff * M_E;
                    } else if (key.substr(0, 4) == "sqrt") {
                        long long n = std::stoll(key.substr(4));
                        result += coeff * std::sqrt(n);
                    }
                }
                return result;
            }
            default:
                return 0.0;
        }
    }

    
    std::string to_string() const {
        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {
                    
                    if (coefficient == static_cast<int>(coefficient)) {
                        return std::to_string(static_cast<int>(coefficient));
                    }
                    return std::to_string(coefficient);
                }
                if (coefficient == 1.0) {
                    return "√" + std::to_string(radicand);
                }
                if (coefficient == -1.0) {
                    return "-√" + std::to_string(radicand);
                }
                
                if (std::abs(coefficient - std::round(coefficient)) < 1e-15) {
                    return std::to_string(static_cast<int>(std::round(coefficient))) + "√" + std::to_string(radicand);
                } else {
                    
                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "√" + std::to_string(radicand);
                }

            case Type::PI:
                if (coefficient == 1.0) {
                    return "π";
                }
                if (coefficient == -1.0) {
                    return "-π";
                }
                
                if (std::abs(coefficient - std::round(coefficient)) < 1e-15) {
                    return std::to_string(static_cast<int>(std::round(coefficient))) + "π";
                } else {
                    
                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "π";
                }

            case Type::E:
                if (coefficient == 1.0) {
                    return "e";
                }
                if (coefficient == -1.0) {
                    return "-e";
                }
                
                if (std::abs(coefficient - std::round(coefficient)) < 1e-15) {
                    return std::to_string(static_cast<int>(std::round(coefficient))) + "e";
                } else {
                    
                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "e";
                }

            case Type::LOG:
                if (coefficient == 1.0) {
                    return "log(" + std::to_string(radicand) + ")";
                }
                if (coefficient == -1.0) {
                    return "-log(" + std::to_string(radicand) + ")";
                }
                
                if (std::abs(coefficient - std::round(coefficient)) < 1e-15) {
                    return std::to_string(static_cast<int>(std::round(coefficient))) + "log(" + std::to_string(radicand) + ")";
                } else {
                    
                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "log(" + std::to_string(radicand) + ")";
                }

            case Type::COMPLEX: {
                std::string result;
                bool first = true;

                
                if (std::abs(constant_term) > 1e-15) {
                    if (std::abs(constant_term - std::round(constant_term)) < 1e-15) {
                        
                        result += std::to_string(static_cast<int>(std::round(constant_term)));
                    } else {
                        
                        std::ostringstream oss;
                        oss << std::fixed << std::setprecision(6) << constant_term;
                        std::string temp = oss.str();
                        
                        temp.erase(temp.find_last_not_of('0') + 1);
                        if (temp.back() == '.') temp.pop_back();
                        result += temp;
                    }
                    first = false;
                }

                
                for (const auto& [key, coeff]: coefficients) {
                    if (std::abs(coeff) < 1e-15) continue;

                    if (!first && coeff > 0) result += " + ";
                    else if (!first && coeff < 0)
                        result += " - ";

                    double abs_coeff = std::abs(coeff);
                    std::string term;

                    if (key == "pi") {
                        if (abs_coeff == 1.0) {
                            term = "π";
                        } else if (abs_coeff == static_cast<int>(abs_coeff)) {
                            term = std::to_string(static_cast<int>(abs_coeff)) + "π";
                        } else {
                            term = std::to_string(abs_coeff) + "π";
                        }
                    } else if (key == "e") {
                        if (abs_coeff == 1.0) {
                            term = "e";
                        } else if (abs_coeff == static_cast<int>(abs_coeff)) {
                            term = std::to_string(static_cast<int>(abs_coeff)) + "e";
                        } else {
                            term = std::to_string(abs_coeff) + "e";
                        }
                    } else if (key.substr(0, 4) == "sqrt") {
                        long long n = std::stoll(key.substr(4));
                        if (abs_coeff == 1.0) {
                            term = "√" + std::to_string(n);
                        } else if (abs_coeff == static_cast<int>(abs_coeff)) {
                            term = std::to_string(static_cast<int>(abs_coeff)) + "√" + std::to_string(n);
                        } else {
                            term = std::to_string(abs_coeff) + "√" + std::to_string(n);
                        }
                    }

                    if (first && coeff < 0) result += "-";
                    result += term;
                    first = false;
                }

                return result.empty() ? "0" : result;
            }
            default:
                return "0";
        }
    }

    
    bool is_zero() const {
        return std::abs(to_double()) < 1e-15;
    }

    
    bool is_rational() const {
        if (type == Type::COMPLEX) {
            return coefficients.empty();
        }
        return false;
    }

    
    void simplify() {
        if (type == Type::COMPLEX) {
            auto it = coefficients.begin();
            while (it != coefficients.end()) {
                if (std::abs(it->second) < 1e-15) {
                    it = coefficients.erase(it);
                } else {
                    ++it;
                }
            }

            
            if (coefficients.empty() && std::abs(constant_term) < 1e-15) {
                constant_term = 0.0;
            }
        }
    }

    
    bool is_positive() const {
        return to_double() > 1e-15;
    }

    
    bool is_negative() const {
        return to_double() < -1e-15;
    }

    
    Irrational abs() const {
        if (is_negative()) {
            return -*this;
        }
        return *this;
    }

    
    Irrational pow(int exponent) const {
        if (exponent == 0) {
            return Irrational::constant(1.0);
        }
        if (exponent == 1) {
            return *this;
        }
        if (exponent == 2 && type == Type::SQRT) {
            
            return Irrational::constant(coefficient * coefficient * radicand);
        }

        
        return Irrational::constant(std::pow(to_double(), exponent));
    }

    
    friend std::ostream& operator<<(std::ostream& os, const Irrational& ir) {
        os << ir.to_string();
        return os;
    }

    
    Type get_type() const { return type; }
};
