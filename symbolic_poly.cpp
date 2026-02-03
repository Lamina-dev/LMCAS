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

// ==========================================
// Extension: Expansion, GCD, Solver
// ==========================================

// Helper for distributing sums: (a+b)*(c+d)
static std::shared_ptr<SymbolicExpr> expand_distribute(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    std::function<void(const std::shared_ptr<SymbolicExpr>&, std::vector<std::shared_ptr<SymbolicExpr>>&)> collect_terms;
    collect_terms = [&](const std::shared_ptr<SymbolicExpr>& e, std::vector<std::shared_ptr<SymbolicExpr>>& terms) {
        if (e->type == SymbolicExpr::Type::Add) {
            for(const auto& op : e->operands) collect_terms(op, terms);
        } else {
            terms.push_back(e);
        }
    };
    
    std::vector<std::shared_ptr<SymbolicExpr>> terms_a;
    collect_terms(a, terms_a);
    
    std::vector<std::shared_ptr<SymbolicExpr>> terms_b;
    collect_terms(b, terms_b);
    
    std::shared_ptr<SymbolicExpr> res = SymbolicExpr::number(0);
    
    // (a1...an) * (b1...bm) -> sum(ai * bj)
    for(const auto& ta : terms_a) {
        for(const auto& tb : terms_b) {
            auto prod = SymbolicExpr::multiply(ta, tb); 
            res = SymbolicExpr::add(res, prod); 
        }
    }
    // simplify should flatten adds
    auto sim = res->simplify();
    return sim;
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::expand() const {
    // 1. Expand operands first
    std::vector<std::shared_ptr<SymbolicExpr>> exp_ops;
    for(const auto& op : operands) {
        if(op) exp_ops.push_back(op->expand());
    }

    if (type == Type::Add) {
        auto res = SymbolicExpr::number(0);
        for(const auto& op : exp_ops) {
            res = SymbolicExpr::add(res, op);
        }
        return res->simplify();
    }
    
    if (type == Type::Multiply) {
        if (exp_ops.empty()) return SymbolicExpr::number(1);
        auto current = exp_ops[0];
        for(size_t i=1; i<exp_ops.size(); ++i) {
            auto next_op = exp_ops[i];
            current = expand_distribute(current, next_op);
        }
        return current;
    }
    
    if (type == Type::Power) {
        auto base = exp_ops[0];
        auto exp = exp_ops[1]; 
        
        if (exp->is_number()) {
            Rational r = exp->convert_rational();
            if (r.is_integer() && r > Rational(1) && r < Rational(20)) {
                // (a+b)^n -> expand
                long long n = (long long)r.to_double();
                auto current = base;
                for(int i=1; i<n; ++i) {
                    current = expand_distribute(current, base);
                }
                return current;
            }
        }
        return SymbolicExpr::power(base, exp);
    }
    
    // Default: restructure
    std::shared_ptr<SymbolicExpr> current = std::make_shared<SymbolicExpr>(type);
    current->identifier = identifier;
    current->number_value = number_value;
    current->operands = exp_ops;
    
    return current; 
}

// ==========================================
// Polynomial GCD
// ==========================================

using PolyMap = std::map<int, std::shared_ptr<SymbolicExpr>>;

static PolyMap expr_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    PolyMap poly;
    
    std::function<void(const std::shared_ptr<SymbolicExpr>&)> add_term;
    add_term = [&](const std::shared_ptr<SymbolicExpr>& term) {
        if (term->type == SymbolicExpr::Type::Add) {
            for(const auto& op : term->operands) add_term(op);
            return;
        }

        int deg = 0;
        std::shared_ptr<SymbolicExpr> coeff = SymbolicExpr::number(1);
        
        // Decompose monomial term
        std::function<void(const std::shared_ptr<SymbolicExpr>&)> analyze = 
            [&](const std::shared_ptr<SymbolicExpr>& t) {
            if (t->type == SymbolicExpr::Type::Multiply) {
                for(auto& op : t->operands) analyze(op);
            } else if (t->type == SymbolicExpr::Type::Power && 
                       t->operands.size() == 2 &&
                       t->operands[0]->type == SymbolicExpr::Type::Variable && 
                       t->operands[0]->identifier == var && 
                       t->operands[1]->is_number()) {
                 Rational r = t->operands[1]->convert_rational();
                 deg += (int)r.to_double();
            } else if (t->type == SymbolicExpr::Type::Variable && t->identifier == var) {
                deg += 1;
            } else {
                coeff = SymbolicExpr::multiply(coeff, t);
            }
        };
        analyze(term);
        
        if (poly.find(deg) == poly.end()) poly[deg] = coeff;
        else poly[deg] = SymbolicExpr::add(poly[deg], coeff);
    };
    
    add_term(expr);
    
    for(auto& [d, c] : poly) c = c->simplify();
    auto it = poly.begin();
    while (it != poly.end()) {
        if (it->second->is_number() && it->second->get_number_value_is_zero()) {
            it = poly.erase(it);
        } else {
            ++it;
        }
    }
    return poly;
}

static std::shared_ptr<SymbolicExpr> poly_to_expr(const PolyMap& poly, const std::string& var) {
    if (poly.empty()) return SymbolicExpr::number(0);
    std::shared_ptr<SymbolicExpr> res = SymbolicExpr::number(0);
    auto x = SymbolicExpr::variable(var);
    
    for(auto& [deg, coeff] : poly) {
        std::shared_ptr<SymbolicExpr> term;
        if (deg == 0) term = coeff;
        else if (deg == 1) term = SymbolicExpr::multiply(coeff, x);
        else term = SymbolicExpr::multiply(coeff, SymbolicExpr::power(x, SymbolicExpr::number(deg)));
        res = SymbolicExpr::add(res, term);
    }
    return res->simplify();
}

static int poly_degree(const PolyMap& p) {
    if (p.empty()) return 0;
    return p.rbegin()->first;
}

static PolyMap poly_div(PolyMap r, PolyMap d, PolyMap& q_out) {
    q_out.clear();
    if (d.empty()) return r;
    
    int deg_d = poly_degree(d);
    auto lead_d = d.rbegin()->second; 
    
    // Limit iterations to prevent hanging on non-divisible symbolic coeffs
    int iter = 0;
    while (!r.empty() && poly_degree(r) >= deg_d && iter++ < 100) {
         int deg_r = poly_degree(r);
         auto lead_r = r.rbegin()->second;
         
         int diff = deg_r - deg_d;
         auto factor = SymbolicExpr::multiply(lead_r, SymbolicExpr::power(lead_d, SymbolicExpr::number(-1)))->simplify(); 
         
         if (q_out.count(diff)) q_out[diff] = SymbolicExpr::add(q_out[diff], factor);
         else q_out[diff] = factor;
         
         for(auto& [deg, coeff] : d) {
             int target_deg = deg + diff;
             auto term = SymbolicExpr::multiply(coeff, factor)->simplify();
             if (r.count(target_deg)) {
                 r[target_deg] = SymbolicExpr::add(r[target_deg], SymbolicExpr::multiply(term, SymbolicExpr::number(-1)))->simplify();
                 if (r[target_deg]->is_number() && r[target_deg]->get_number_value_is_zero()) r.erase(target_deg);
             } else {
                 r[target_deg] = SymbolicExpr::multiply(term, SymbolicExpr::number(-1))->simplify();
             }
         }
         
         // Force remove leading term to ensure progress if symbolic cancellation failed slightly 
         // logic ensures it SHOULD cancel if arithmetic is exact.
         if (r.count(deg_r)) r.erase(deg_r);
    }
    return r;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    std::string var = "x";
    
    // Find variables
    // For now try 'x'
    
    auto pa = expr_to_poly(a->expand(), var);
    auto pb = expr_to_poly(b->expand(), var);
    
    while (!pb.empty()) {
        PolyMap q;
        auto r = poly_div(pa, pb, q);
        pa = pb;
        pb = r;
    }
    
    // Normalize
    if (!pa.empty()) {
        auto lc = pa.rbegin()->second;
        for(auto& [d, c] : pa) {
            c = SymbolicExpr::multiply(c, SymbolicExpr::power(lc, SymbolicExpr::number(-1)))->simplify();
        }
    }
    return poly_to_expr(pa, var);
}

// ==========================================
// Solver
// ==========================================


std::shared_ptr<SymbolicExpr> SymbolicExpr::substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const {
    if (type == SymbolicExpr::Type::Variable && identifier == var_name) {
        return value;
    }
    
    // Recursive copy with substitution
    std::vector<std::shared_ptr<SymbolicExpr>> new_ops;
    bool changed = false;
    for(const auto& op : operands) {
        if(!op) { new_ops.push_back(nullptr); continue; }
        auto new_op = op->substitute(var_name, value);
        new_ops.push_back(new_op);
        if (new_op != op) changed = true; // Pointer comparison
    }
    // Deep check if pointer check isn't sufficient for shared_ptr
    
    // Minimal optimization: if no operands changed, return copy of self (or self if const)
    // Here we return a new simplified expression usually
    
    if (type == SymbolicExpr::Type::Number) return std::make_shared<SymbolicExpr>(*this);
    if (type == SymbolicExpr::Type::Variable) return std::make_shared<SymbolicExpr>(*this);

    auto res = std::make_shared<SymbolicExpr>(type);
    res->number_value = number_value;
    res->identifier = identifier;
    res->operands = new_ops;
    
    return res->simplify(); 
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_single_poly(const std::shared_ptr<SymbolicExpr>& eq, const std::string& var) {
    auto poly = expr_to_poly(eq->expand(), var);
    int deg = poly_degree(poly);
    
    if (deg == 1) {
        // ax + b = 0 => x = -b/a
        auto a = poly[1];
        auto b = poly.count(0) ? poly[0] : SymbolicExpr::number(0);
        return {SymbolicExpr::multiply(SymbolicExpr::multiply(b, SymbolicExpr::number(-1)), SymbolicExpr::power(a, SymbolicExpr::number(-1)))->simplify()};
    }
    
    if (deg == 2) {
        // ax^2 + bx + c = 0
        auto a = poly.count(2) ? poly[2] : SymbolicExpr::number(0);
        auto b = poly.count(1) ? poly[1] : SymbolicExpr::number(0);
        auto c = poly.count(0) ? poly[0] : SymbolicExpr::number(0);
        
        auto delta = SymbolicExpr::add(SymbolicExpr::power(b, SymbolicExpr::number(2)), 
                     SymbolicExpr::multiply(SymbolicExpr::number(-4), SymbolicExpr::multiply(a, c)));
        
        auto den = SymbolicExpr::multiply(SymbolicExpr::number(2), a);
        auto inv_den = SymbolicExpr::power(den, SymbolicExpr::number(-1));

        // Root 1: (-b + sqrt(delta)) / 2a
        auto num1 = SymbolicExpr::add(SymbolicExpr::multiply(b, SymbolicExpr::number(-1)), SymbolicExpr::sqrt(delta));
        auto root1 = SymbolicExpr::multiply(num1, inv_den)->simplify();
        
        // Root 2: (-b - sqrt(delta)) / 2a
        auto num2 = SymbolicExpr::add(SymbolicExpr::multiply(b, SymbolicExpr::number(-1)), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::sqrt(delta)));
        auto root2 = SymbolicExpr::multiply(num2, inv_den)->simplify();

        return {root1, root2};
    }
    
    return {}; 
}

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var) {
    return solve_single_poly(eq, var);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var) {
    auto pa = expr_to_poly(a->expand(), var);
    auto pb = expr_to_poly(b->expand(), var);
    
    std::shared_ptr<SymbolicExpr> res = SymbolicExpr::number(1);
    
    while(true) {
        int deg_a = poly_degree(pa);
        int deg_b = poly_degree(pb);
        
        if (deg_b == -1 || (pb.empty())) { // Zero polynomial
             return SymbolicExpr::number(0);
        }

        if (deg_b == 0) {
            auto const_b = pb.at(0); // Coeff at 0
            res = SymbolicExpr::multiply(res, SymbolicExpr::power(const_b, SymbolicExpr::number(deg_a)))->simplify();
            break;
        }
        
        if (deg_a < deg_b) {
            std::swap(pa, pb);
            std::swap(deg_a, deg_b);
            if ((deg_a % 2 == 1) && (deg_b % 2 == 1)) {
                 res = SymbolicExpr::multiply(res, SymbolicExpr::number(-1));
            }
            continue;
        }
        
        PolyMap q; // quotient
        PolyMap r = poly_div(pa, pb, q);
        int deg_r = poly_degree(r);
        
        auto lc_b = pb.rbegin()->second;
        
        // Res(A,B) = (-1)^(mn) * b_n^(m-deg(R)) * Res(B,R)
        // Note: using basic Euclidean remainder
        
        // Sign
        if ((deg_a % 2 == 1) && (deg_b % 2 == 1)) {
             res = SymbolicExpr::multiply(res, SymbolicExpr::number(-1));
        }
        
        // Leading coeff power
        auto factor = SymbolicExpr::power(lc_b, SymbolicExpr::number(deg_a - deg_r));
        res = SymbolicExpr::multiply(res, factor);
        
        pa = pb;
        pb = r;
    }
    
    return res->simplify(); 
}

// 辅助：检查表达式是否依赖于某些变量
static bool depends_on(const std::shared_ptr<SymbolicExpr>& expr, const std::vector<std::string>& vars) {
    if (expr->type == SymbolicExpr::Type::Variable) {
        for (const auto& v : vars) {
            if (expr->identifier == v) return true;
        }
        return false;
    }
    for (const auto& op : expr->operands) {
        if (depends_on(op, vars)) return true;
    }
    return false;
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations, 
        const std::vector<std::string>& vars) {

    size_t m = equations.size();
    size_t n = vars.size();

    // 尝试构建线性方程组矩阵 [A|b]
    // 方程形式： sum(a_ij * x_j) + C = 0
    // 移项后： sum(a_ij * x_j) = -C
    // 增广矩阵最后一列存 -C
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> augmented_matrix(m);
    bool is_linear = true;

    for (size_t i = 0; i < m; ++i) {
        // 先展开以分离项
        auto eq = equations[i]->expand()->simplify(); // 确保 fully expanded
        augmented_matrix[i].resize(n + 1, SymbolicExpr::number(0));
        
        // 收集所有加法项
        std::vector<std::shared_ptr<SymbolicExpr>> terms;
        auto collect_terms = [&](auto&& self, const std::shared_ptr<SymbolicExpr>& e) -> void {
            if (e->type == Type::Add) {
                for (const auto& op : e->operands) self(self, op);
            } else {
                terms.push_back(e);
            }
        };
        collect_terms(collect_terms, eq);
        
        std::shared_ptr<SymbolicExpr> constant_part = SymbolicExpr::number(0);

        for (const auto& term : terms) {
            // 分析每一项包含的变量
            std::vector<int> found_indices;
            for (size_t j = 0; j < n; ++j) {
                if (depends_on(term, {vars[j]})) {
                    found_indices.push_back(j);
                }
            }
            
            if (found_indices.empty()) {
                // 常数项
                constant_part = SymbolicExpr::add(constant_part, term);
            } else if (found_indices.size() == 1) {
                // 只含有一个变量，检查线性性
                int var_idx = found_indices[0];
                std::string var_name = vars[var_idx];
                
                // 为了验证线性性，求导数。如果是线性项 c * x，导数为 c，且 c 不依赖于任何 vars
                auto coeff = term->differentiate(var_name)->simplify();
                
                if (depends_on(coeff, vars)) {
                    // 系数依赖于任何待解变量（包括 x 自身，例如 x^2 -> 2x），则非线性
                    is_linear = false;
                    break;
                }
                
                // 还要检查是否真的是线性项。
                // 比如 sin(x) 对 x 求导是 cos(x)，此时 depends_on(coeff, vars) 会为 true。
                // 但如果 floor(x) 这种？ Cas 不支持。
                // 还有一个特例： x * y (两个变量)，这里 found_indices.size() 会是 2，已经被上面排除了。
                // 所以这里应该是安全的。
                
                augmented_matrix[i][var_idx] = SymbolicExpr::add(augmented_matrix[i][var_idx], coeff);
            } else {
                // 包含多个变量（如 x*y），非线性
                is_linear = false;
                break;
            }
        }
        
        if (!is_linear) break;
        
        // b = -constant
        augmented_matrix[i][n] = SymbolicExpr::multiply(constant_part, SymbolicExpr::number(-1))->simplify();
    }
    
    // 如果是线性，使用 RREF 求解
    if (is_linear) {
        // 化简矩阵中的每个元素
        for(auto& row : augmented_matrix) {
            for(auto& cell : row) cell = cell->simplify();
        }

        auto mat = SymbolicExpr::matrix(augmented_matrix);
        // 使用 RREF
        auto solved_mat_expr = SymbolicExpr::rref(mat);
        
        if (!solved_mat_expr || solved_mat_expr->type != Type::Matrix) return {};
        
        auto& solved_rows = solved_mat_expr->operands;
        if (solved_rows.empty()) return {};

        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        
        // 解析 RREF 结果 (Back substitution if not diagonal, but RREF should be mostly diagonal)
        // RREF 形式：
        // [1 0 2 | 5] -> x + 2z = 5 -> x = 5 - 2z
        // [0 1 3 | 6] -> y + 3z = 6 -> y = 6 - 3z
        // [0 0 0 | 0]
        
        // 识别 pivot 变量和自由变量
        // 每一行寻找第一个非零元素（pivot）
        // 如果 pivot 是 1，且该列其他为 0 (RREF定义)
        
        // 由于是符号计算，判断“非零”比较困难 (depends on simplify)。
        // 假设 simplify 足够强，能把 0 化简为 Number(0)
        
        std::vector<int> pivot_col_for_row(m, -1);
        std::vector<bool> is_free_var(n, true);
        
        for (size_t i = 0; i < solved_rows.size(); ++i) {
            auto& row_vec = solved_rows[i]->operands;
            for (size_t j = 0; j < n; ++j) {
                // 检查是否非零
                bool is_zero = false;
                if (row_vec[j]->is_number()) {
                    auto val = row_vec[j]->convert_rational();
                    if (val == 0) is_zero = true;
                }
                
                if (!is_zero) {
                    pivot_col_for_row[i] = j;
                    is_free_var[j] = false;
                    
                    // 检查无解情况: [0 0 ... 0 | 1]
                    // 但这里 j < n，所以这是系数部分的非零
                    // 如果整行系数部分都是 0，但最后部分非0，则无解
                    break;
                }
            }
            
            // 检查无解
            if (pivot_col_for_row[i] == -1) {
                auto rhs = row_vec[n];
                bool rhs_zero = false;
                 if (rhs->is_number()) {
                    auto val = rhs->convert_rational();
                    if (val == 0) rhs_zero = true;
                }
                
                if (!rhs_zero) {
                    // 0 = non_zero -> 无解
                    return {}; 
                }
            }
        }
        
        // 构建解
        // 对于自由变量，我们暂时无法用 "t" 表示返回一般解 (TODO)
        // 目前如果存在自由变量，我们可能返回空，或者仅仅返回确定解的部分？
        // 标准行为通常是返回参数化解。这里暂时只支持唯一解情况，或者若包含自由变量则返回部分绑定
        
        // 从下往上回代（RREF 其实不需要回代，直接移项即可）
        for (int i = m - 1; i >= 0; --i) {
            int p = pivot_col_for_row[i];
            if (p == -1) continue;
            
            // x_p + sum(c_k * x_k) = rhs
            // x_p = rhs - sum(c_k * x_k)
            auto val = solved_rows[i]->operands[n];
            
            for (size_t j = p + 1; j < n; ++j) {
                auto coeff = solved_rows[i]->operands[j];
                bool is_zero = false;
                if (coeff->is_number() && coeff->convert_rational() == 0) is_zero = true;
                
                if (!is_zero) {
                    // 这是一个依赖自由变量的解
                    // 如果 user 没有提供自由变量的值，结果将依赖于该变量
                    // 我们可以保留变量名在表达式中
                    auto term = SymbolicExpr::multiply(coeff, SymbolicExpr::variable(vars[j]));
                    val = SymbolicExpr::add(val, SymbolicExpr::multiply(term, SymbolicExpr::number(-1)));
                }
            }
            
            // Pivot 系数归一化 (RREF 应该已经是 1，除非符号计算导致未能除尽)
            auto pivot_val = solved_rows[i]->operands[p];
            bool is_one = false;
            if (pivot_val->is_number() && pivot_val->convert_rational() == 1) is_one = true;
            
            if (!is_one) {
                val = SymbolicExpr::multiply(val, SymbolicExpr::power(pivot_val, SymbolicExpr::number(-1)));
            }
            
            solution[vars[p]] = val->simplify();
        }
        
        // 填充自由变量： x_k = x_k
        // 这一步对于 CAS 很重要，告诉用户哪些是自由的
        for (size_t j = 0; j < n; ++j) {
            if (is_free_var[j]) {
                // solution[vars[j]] = SymbolicExpr::variable(vars[j]); 
                // 或者不放入 map，表示它是自由的
            }
        }
        
        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> res;
        res.push_back(solution);
        return res;
    }

     // Substitution Solver Strategy (Fallback)
     std::vector<std::shared_ptr<SymbolicExpr>> current_eqs = equations;
     std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
     
     for(size_t step = 0; step < vars.size(); ++step) {
         bool progress = false;
         
         for(int i = 0; i < (int)current_eqs.size(); ++i) {
            auto& eq = current_eqs[i];
            
            for(const auto& v : vars) {
                if (solution.count(v)) continue;
                
                auto res_list = solve(eq, v);
                if (!res_list.empty()) {
                    auto sol_val = res_list[0];
                    solution[v] = sol_val;
                    progress = true;
                    
                    current_eqs.erase(current_eqs.begin() + i);
                    
                    for(auto& other : current_eqs) {
                         other = other->substitute(v, sol_val);
                    }
                    
                    for(auto& [sv, sval] : solution) {
                        if (sv != v) {
                            solution[sv] = sval->substitute(v, sol_val);
                        }
                    }
                    
                    break;
                }
            }
            if (progress) break;
         }
         if (!progress) break;
     }
     
     std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> res;
     if (!solution.empty()) res.push_back(solution);
     return res;
}
