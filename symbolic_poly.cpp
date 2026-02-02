#include "symbolic.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include <map>
#include <cstdio>

//#define DEBUG_FACTOR 1
#if 0
#define DBG(fmt, ...) do { fprintf(stderr, "FACTOR_DBG: " fmt "\n", ##__VA_ARGS__); fflush(stderr); } while(0)
#else
#define DBG(fmt, ...)
#endif

// 辅助函数：计算两个整数的最大公约数
static long long gcd(long long a, long long b) {
    return std::gcd(a, b);
}

// Helper to check if a Rational is a perfect square and return its root
static bool rational_sqrt_exact(const Rational& r, Rational& out_root) {
    if (r.get_numerator().IsNegative() || r.get_denominator().IsNegative()) return false; 
    
    BigInt num = r.get_numerator();
    BigInt den = r.get_denominator();
    
    if (num.is_perfect_square() && den.is_perfect_square()) {
        out_root = Rational(num.sqrt(), den.sqrt());
        return true;
    }
    return false;
}

// 辅助函数：提取公因式
static std::shared_ptr<SymbolicExpr> factor_common(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || expr->type != SymbolicExpr::Type::Add) {
        return expr; 
    }

    const auto& terms = expr->operands;
    if (terms.empty()) return expr;

    // 1. 提取数值公因子
    long long common_coeff = 0;
    bool first = true;
    
    // 我们主要处理整数系数的情况
    for (const auto& term : terms) {
        if (!term) continue;
        
        long long coeff = 1; // Default coeff is 1 (implicit)
        bool has_coeff = false; // Whether we successfully extracted a numeric coefficient
        
        // Check for (Number * ...) pattern
        if (term->type == SymbolicExpr::Type::Multiply && !term->operands.empty() && term->operands[0] && term->operands[0]->is_number()) {
            try {
                Rational r = term->operands[0]->convert_rational();
                if (r.is_integer()) {
                    // Use string conversion to avoid potential issues with direct BigInt->int conversion
                    coeff = std::abs(std::stoll(r.get_numerator().ToString()));
                } else {
                    double d = r.to_double();
                    coeff = std::abs((long long)d);
                }
                has_coeff = true;
            } catch (...) {
                // Ignore conversion errors, treat as coeff 1
            }
        } 
        // Check for Number pattern
        else if (term->is_number()) {
             try {
                Rational r = term->convert_rational();
                if (r.is_integer()) {
                    coeff = std::abs(std::stoll(r.get_numerator().ToString()));
                } else {
                     double d = r.to_double();
                     coeff = std::abs((long long)d);
                }
                has_coeff = true;
             } catch (...) {}
        }
        
        if (!has_coeff) coeff = 1;

        if (first) {
            common_coeff = coeff;
            first = false;
        } else {
            common_coeff = gcd(common_coeff, coeff);
        }
    }
    
    if (common_coeff == 0) common_coeff = 1;

    // 如果公因子为1，直接返回
    if (common_coeff == 1) {
        return expr; 
    }

    auto common_factor = SymbolicExpr::number((int)common_coeff);

    // 构建剩余部分
    std::vector<std::shared_ptr<SymbolicExpr>> new_terms;
    for (const auto& term : terms) {
        std::shared_ptr<SymbolicExpr> new_term;
        bool divided = false;

        // Try to divide by common_coeff
        if (term->is_number()) {
             Rational r = term->convert_rational();
             new_term = SymbolicExpr::number(r / Rational(common_coeff));
             divided = true;
        } else if (term->type == SymbolicExpr::Type::Multiply && term->operands[0]->is_number()) {
            Rational r = term->operands[0]->convert_rational();
            Rational new_c = r / Rational(common_coeff);
            
            std::vector<std::shared_ptr<SymbolicExpr>> mult_ops = term->operands;
            mult_ops[0] = SymbolicExpr::number(new_c);
            
            // Normalize: if 1, remove it
            if (new_c == Rational(1) && mult_ops.size() > 1) {
                mult_ops.erase(mult_ops.begin());
            }
            
            // Rebuild chain
            if (mult_ops.size() == 1) {
                 new_term = mult_ops[0];
            }
            else {
                auto acc = mult_ops[0];
                for(size_t i=1; i<mult_ops.size(); ++i) {
                    acc = SymbolicExpr::multiply(acc, mult_ops[i]);
                }
                new_term = acc;
            }
            divided = true;
        }

        if (!divided) {
             Rational recip(1, common_coeff);
             new_term = SymbolicExpr::multiply(SymbolicExpr::number(recip), term);
        }
        new_terms.push_back(new_term);
    }

    // Reconstruct Add
    if (new_terms.empty()) return SymbolicExpr::number(0);
    auto rest = new_terms[0];
    for(size_t i=1; i<new_terms.size(); ++i) {
        rest = SymbolicExpr::add(rest, new_terms[i]);
    }
    
    return SymbolicExpr::multiply(common_factor, rest);
}

// Helper to flatten nested Add expressions
static void collect_add_terms(const std::shared_ptr<SymbolicExpr>& expr, std::vector<std::shared_ptr<SymbolicExpr>>& out_terms) {
    if (expr && expr->type == SymbolicExpr::Type::Add) {
        for (const auto& op : expr->operands) {
            collect_add_terms(op, out_terms);
        }
    } else if (expr) {
        out_terms.push_back(expr);
    }
}

// 占位符：二次多项式因式分解
static std::shared_ptr<SymbolicExpr> factor_quadratic(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || expr->type != SymbolicExpr::Type::Add) return expr;

    // Flatten terms from nested Adds
    std::vector<std::shared_ptr<SymbolicExpr>> flat_terms;
    collect_add_terms(expr, flat_terms);
    DBG("factor_quadratic: flat_terms size %zu", flat_terms.size());

    // 1. Identify coefficients a, b, c for ax^2 + bx + c
    Rational a(0), b(0), c(0);
    std::string var_name;
    bool is_quadratic = false;

    // Temporary storage to track if terms are accounted for.
    // However, if we fail to form a quadratic, we just return expr.
    // If we succeed, we replace it entirely.
    // If there are left over terms (e.g. x^2 + 5x + 6 + y), we haven't handled that mixing yet.
    // Ideally we should isolate the quadratic part.
    // For now, assume the WHOLE expression is the quadratic in one variable + constant.

    for (const auto& term : flat_terms) {
        bool handled = false;
        if (!term) continue;
        DBG("Processing term type %d", (int)term->type);
        
        // Check for x^2
        if (term->type == SymbolicExpr::Type::Power && term->operands.size() == 2 && 
            term->operands[0]->type == SymbolicExpr::Type::Variable && 
            term->operands[1]->is_number()) {
            try {
                if (term->operands[1]->convert_rational() == Rational(2)) {
                    // It is x^2
                    if (!is_quadratic || var_name == term->operands[0]->identifier) {
                        a = a + Rational(1);
                        var_name = term->operands[0]->identifier;
                        is_quadratic = true;
                        handled = true;
                        DBG("Found x^2, var=%s", var_name.c_str());
                    }
                }
            } catch (...){}
        }
        
        // Check for ax^2 or bx
        if (!handled && term->type == SymbolicExpr::Type::Multiply && term->operands.size() == 2) {
            auto op1 = term->operands[0];
            auto op2 = term->operands[1];
            // Normalize so op1 is number if possible
            if (!op1->is_number() && op2->is_number()) std::swap(op1, op2);
            
            if (op1->is_number()) {
                try {
                     Rational coeff = op1->convert_rational();
                     // Check op2 is x^2
                     if (op2->type == SymbolicExpr::Type::Power && op2->operands.size() == 2 &&
                         op2->operands[0]->type == SymbolicExpr::Type::Variable &&
                         op2->operands[1]->is_number() && op2->operands[1]->convert_rational() == Rational(2)) {
                             if (!is_quadratic || var_name == op2->operands[0]->identifier) {
                                 a = a + coeff;
                                 var_name = op2->operands[0]->identifier;
                                 is_quadratic = true;
                                 handled = true;
                                 DBG("Found ax^2, a=%s", coeff.to_string().c_str());
                             }
                     }
                     // Check op2 is x
                     else if (op2->type == SymbolicExpr::Type::Variable) {
                         if (var_name.empty() || var_name == op2->identifier) {
                             b = b + coeff;
                             if (var_name.empty()) var_name = op2->identifier;
                             handled = true;
                             DBG("Found bx, b=%s", coeff.to_string().c_str());
                         }
                     }
                } catch(...) {}
            }
        }
        
        // Check for x
        if (!handled && term->type == SymbolicExpr::Type::Variable) {
             if (var_name.empty() || var_name == term->identifier) {
                 b = b + Rational(1);
                 if (var_name.empty()) var_name = term->identifier;
                 handled = true;
                 DBG("Found x, b+=1");
             }
        }
        
        // Check for c (number)
        if (!handled && term->is_number()) {
            try {
                c = c + term->convert_rational();
                handled = true;
                DBG("Found c, c=%s", c.to_string().c_str());
            } catch (...) {}
        }
        
        // If we found a term that doesn't fit the single-variable quadratic model (e.g. variable y when x is selected),
        // then this simple quadratic factorizer should probably abort or be smarter.
        // For this task, we assume the inputs are clean quadratics or mixed simple sums.
        // If we have 'x^2 + 5x + 6 + y', treating 'y' as part of constant 'c' only works if 'y' is considered constant relative to 'x'.
        if (!handled) {
            DBG("Term unhandled: type %d", (int)term->type);
            // For safety, if we encounter something we don't understand, abort.
            return expr;
        }
    }

    if (!is_quadratic || a.is_zero()) {
        DBG("Not quadratic or a=0");
        return expr; 
    }

    // 2. Calculate discriminant: b^2 - 4ac
    Rational discriminant = b.power(BigInt(2)) - Rational(4) * a * c;
    DBG("Discriminant: %s", discriminant.to_string().c_str());

    // 3. Check if discriminant is a perfect square
    Rational disc_root;
    if (rational_sqrt_exact(discriminant, disc_root)) {
        DBG("Perfect square discriminant, root: %s", disc_root.to_string().c_str());
        
        Rational two_a = Rational(2) * a;
        if (two_a.is_zero()) {
            DBG("Division by zero in quadratic factor: 2a is zero. a=%s", a.to_string().c_str());
            return expr;
        }

        // Roots are rational: (-b ± sqrt(D)) / 2a
        Rational root1 = (-b + disc_root) / two_a;
        Rational root2 = (-b - disc_root) / two_a;

        // Factor is a * (x - root1) * (x - root2)
        // Note: root = r => factor is (x - r)
        // If root is -2, factor is (x - (-2)) = (x + 2)
        
        auto x = SymbolicExpr::variable(var_name);
        
        // Helper to create (x - r)
        auto make_factor = [&](const Rational& r) -> std::shared_ptr<SymbolicExpr> {
             // x - r = x + (-r)
             if (r.is_zero()) return x;
             Rational neg_r = -r;
             if (neg_r.get_denominator() == BigInt(1) && neg_r.get_numerator().IsNegative()) {
                  // x + (-3)
                  return SymbolicExpr::add(x, SymbolicExpr::number(neg_r));
             } else {
                  return SymbolicExpr::add(x, SymbolicExpr::number(neg_r));
             }
        };

        auto factor1 = make_factor(root1);
        auto factor2 = make_factor(root2);
        
        std::shared_ptr<SymbolicExpr> result;
        if (a == Rational(1)) {
            result = SymbolicExpr::multiply(factor1, factor2);
        } else {
            result = SymbolicExpr::multiply(SymbolicExpr::number(a), SymbolicExpr::multiply(factor1, factor2));
        }

        return result; 
    }

    return expr; 
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::factor() const {
    // 1. Recursive factor on operands (Depth First)
    std::vector<std::shared_ptr<SymbolicExpr>> new_ops;
    bool changed = false;
    
    // Base case: no operands (Leaf nodes like Number, Variable)
    if (operands.empty()) {
        auto copy = std::make_shared<SymbolicExpr>(type);
        copy->number_value = number_value;
        copy->identifier = identifier;
        return copy;
    }
    
    for (const auto& op : operands) {
        if (!op) continue;
        auto new_op = op->factor();
        new_ops.push_back(new_op);
        if (new_op != op) changed = true;
    }

    // Preserve structure with factored operands
    std::shared_ptr<SymbolicExpr> current = std::make_shared<SymbolicExpr>(type);
    current->operands = new_ops;
    current->number_value = number_value;
    current->identifier = identifier;
    
    if (current->type == Type::Add) {
        // Try to extract common factors
        auto fact = factor_common(current);
        if (fact != current) return fact;
        
        // Try quadratic factorization
        auto quad = factor_quadratic(current);
        if (quad != current) return quad;
    }

    return current;
}
