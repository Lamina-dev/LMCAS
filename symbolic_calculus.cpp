#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <vector>
#include <string>

// --- Helper Declarations ---

// Check if expression tends to 0 at var -> val (Forward declared for limit)
static bool tends_to_zero(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& val) {
    auto sub = expr->substitute(var, val)->simplify();
    return sub->is_number() && sub->convert_rational() == Rational(0);
}

// Helper to check dependency
static bool has_variable(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr) return false;
    if (expr->type == SymbolicExpr::Type::Variable && expr->identifier == var) return true;
    for (const auto& op : expr->operands) {
        if (has_variable(op, var)) return true;
    }
    return false;
}

// --- Additional Symbolic API Implementations ---

std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    // num / den -> num * den^-1
    if (!den) return num; // Error handling usually
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    auto res = std::make_shared<SymbolicExpr>(Type::Integral);
    res->operands.push_back(expr);
    res->identifier = var;
    return res;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point) {
    auto res = std::make_shared<SymbolicExpr>(Type::Limit);
    res->operands.push_back(expr);
    res->operands.push_back(point);
    res->identifier = var;
    return res;
}

// --- Calculus Core ---

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point) const {
    // 1. Check for L'Hopital Rule Candidate Indeterminate Form (0/0)
    if (type == Type::Multiply) {
         auto num = operands[0];
         std::shared_ptr<SymbolicExpr> den = nullptr;
         
         if (operands[1]->type == Type::Power && operands[1]->operands[1]->is_number()) {
             auto exp_val = operands[1]->operands[1]->convert_rational();
             if (exp_val == Rational(-1)) {
                 den = operands[1]->operands[0];
             }
         }
         
         if (den) {
             if (tends_to_zero(num, var, point) && tends_to_zero(den, var, point)) {
                 auto d_num = num->differentiate(var);
                 auto d_den = den->differentiate(var);
                 auto new_limit_val = SymbolicExpr::divide(d_num, d_den)->simplify();
                 return new_limit_val->limit(var, point); 
             }
         }
    }

    // 2. Direct Substitution
    return this->substitute(var, point)->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integrate(const std::string& var_name) const {
    // Simplest: reconstruct 'this' into a shared_ptr if needed.
    // For Number/Variable, we generate new nodes.
    
    switch (type) {
        case Type::Number:
             // int c dx = cx
             // Use std::visit to handle the variant number_value
             {
                 auto num_node = std::visit([](const auto& val) { return SymbolicExpr::number(val); }, number_value);
                 return SymbolicExpr::multiply(num_node, SymbolicExpr::variable(var_name))->simplify();
             }
             
        case Type::Variable:
             if (identifier == var_name) {
                 // int x dx = x^2/2
                 auto x = SymbolicExpr::variable(identifier);
                 auto two = SymbolicExpr::number(2);
                 return SymbolicExpr::divide(SymbolicExpr::power(x, two), two)->simplify();
             }
             // int y dx = yx
             return SymbolicExpr::multiply(SymbolicExpr::variable(identifier), SymbolicExpr::variable(var_name))->simplify();
             
        case Type::Add:
             // int (u+v) = int u + int v
             return SymbolicExpr::add(operands[0]->integrate(var_name), operands[1]->integrate(var_name))->simplify();
             
        case Type::Power: {
             auto base = operands[0];
             auto exp = operands[1];
             if (base->type == Type::Variable && base->identifier == var_name) {
                 if (exp->is_number()) {
                     // int x^n = x^(n+1)/(n+1)
                     // Check if n = -1
                     if (exp->convert_rational() == Rational(-1)) {
                         return SymbolicExpr::ln(base);
                     }
                     auto n_plus_1 = SymbolicExpr::add(exp, SymbolicExpr::number(1));
                     return SymbolicExpr::divide(SymbolicExpr::power(base, n_plus_1), n_plus_1)->simplify();
                 }
             }
             // Fallback
             break;
        }
    }
    
    // Fallback: Reconstruct current node type and operands for the integral wrapper
    auto self_copy = std::make_shared<SymbolicExpr>(type);
    self_copy->number_value = number_value;
    self_copy->identifier = identifier;
    self_copy->operands = operands; // Shared ptrs copy is cheap
    
    return SymbolicExpr::make_integral(self_copy, var_name); 
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate(const std::string& var_name) const {
    switch (type) {
        case Type::Number:
        case Type::Infinity:
             return SymbolicExpr::number(0);
             
        case Type::Variable:
             if (identifier == var_name) return SymbolicExpr::number(1);
             return SymbolicExpr::number(0);
             
        case Type::Add:
             return SymbolicExpr::add(operands[0]->differentiate(var_name), operands[1]->differentiate(var_name))->simplify();
             
        case Type::Multiply: {
             auto u = operands[0];
             auto v = operands[1];
             auto du = u->differentiate(var_name);
             auto dv = v->differentiate(var_name);
             // product rule: u'v + uv'
             return SymbolicExpr::add(
                 SymbolicExpr::multiply(du, v),
                 SymbolicExpr::multiply(u, dv)
             )->simplify();
        }
        
        case Type::Power: {
             auto u = operands[0];
             auto v = operands[1];
             auto du = u->differentiate(var_name);
             auto dv = v->differentiate(var_name);
             // if v is constant number
             if (v->is_number() || !has_variable(v, var_name)) {
                 // n * u^(n-1) * u'
                 auto n_minus_1 = SymbolicExpr::add(v, SymbolicExpr::number(-1));
                 return SymbolicExpr::multiply(
                     SymbolicExpr::multiply(v, SymbolicExpr::power(u, n_minus_1)),
                     du
                 )->simplify();
             }
             // u^v * (v'ln(u) + v*u'/u)
             return SymbolicExpr::multiply(
                 SymbolicExpr::power(u, v),
                 SymbolicExpr::add(
                     SymbolicExpr::multiply(dv, SymbolicExpr::ln(u)),
                     SymbolicExpr::multiply(SymbolicExpr::multiply(v, du), SymbolicExpr::power(u, SymbolicExpr::number(-1)))
                 )
             )->simplify();
        }
        
        case Type::Sqrt: {
             // sqrt(u)' = u' / 2sqrt(u)
             auto u = operands[0];
             auto du = u->differentiate(var_name);
             return SymbolicExpr::multiply(
                 du,
                 SymbolicExpr::power(
                    SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::sqrt(u)),
                    SymbolicExpr::number(-1)
                 )
             )->simplify();
        }
        
        case Type::Sin: {
             // cos(u)u'
             auto u = operands[0];
             return SymbolicExpr::multiply(SymbolicExpr::cos(u), u->differentiate(var_name))->simplify();
        }
        
        case Type::Cos: {
             // -sin(u)u'
             auto u = operands[0];
             return SymbolicExpr::multiply(
                 SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::sin(u)),
                 u->differentiate(var_name)
             )->simplify();
        }
        
        case Type::Tan: {
             // (1+tan^2(u))u' or sec^2(u)u'
             // Let's use 1/cos^2(u) * u'
             auto u = operands[0];
             auto sec2 = SymbolicExpr::power(SymbolicExpr::cos(u), SymbolicExpr::number(-2));
             return SymbolicExpr::multiply(sec2, u->differentiate(var_name))->simplify();
        }
        
        case Type::ArcSin: {
             // u' / sqrt(1 - u^2)
             auto u = operands[0];
             auto du = u->differentiate(var_name);
             auto denom = SymbolicExpr::sqrt(
                 SymbolicExpr::add(SymbolicExpr::number(1), 
                     SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(u, SymbolicExpr::number(2))))
             );
             return SymbolicExpr::multiply(du, SymbolicExpr::power(denom, SymbolicExpr::number(-1)))->simplify();
        }

        case Type::ArcCos: {
             // -u' / sqrt(1 - u^2)
             auto u = operands[0];
             auto du = u->differentiate(var_name);
             auto denom = SymbolicExpr::sqrt(
                 SymbolicExpr::add(SymbolicExpr::number(1), 
                     SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(u, SymbolicExpr::number(2))))
             );
             auto neg_du = SymbolicExpr::multiply(SymbolicExpr::number(-1), du);
             return SymbolicExpr::multiply(neg_du, SymbolicExpr::power(denom, SymbolicExpr::number(-1)))->simplify();
        }

        case Type::ArcTan: {
             // u' / (1 + u^2)
             auto u = operands[0];
             auto du = u->differentiate(var_name);
             auto denom = SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::power(u, SymbolicExpr::number(2)));
             return SymbolicExpr::multiply(du, SymbolicExpr::power(denom, SymbolicExpr::number(-1)))->simplify();
        }

        case Type::Sinh: {
             // cosh(u) * u'
             // Note: SymbolicExpr::cosh helper might not be declared in symbolic.hpp yet, use make generic
             // Assuming make generic logic or wait for helper.
             // For now, let's comment out or implement generic constructor
             // return SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(Type::Cosh, u), u->differentiate(var_name))->simplify();
             return SymbolicExpr::number(0); // TODO: implement Hyperbolic helpers in symbolic.hpp
        }

        case Type::Cosh: {
             // sinh(u) * u'
             // return SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(Type::Sinh, u), u->differentiate(var_name))->simplify();
             return SymbolicExpr::number(0); // TODO
        }

        case Type::Tanh: {
             // sech^2(u) * u' = (1 - tanh^2(u)) * u'
             // return ...
             return SymbolicExpr::number(0); // TODO
        }

        case Type::Ln: {
             // u'/u
             auto u = operands[0];
             return SymbolicExpr::multiply(u->differentiate(var_name), SymbolicExpr::power(u, SymbolicExpr::number(-1)))->simplify();
        }

        case Type::Log: {
             // log_b(u) = ln(u) / ln(b)
             // Differentiate the equivalent expression using chain rule/product rule logic already implemented
             auto u = operands[0];
             auto b = operands[1];
             auto ln_u = SymbolicExpr::ln(u);
             auto ln_b = SymbolicExpr::ln(b);
             auto expr = SymbolicExpr::multiply(ln_u, SymbolicExpr::power(ln_b, SymbolicExpr::number(-1)));
             return expr->differentiate(var_name);
        }

        default:
             return SymbolicExpr::number(0);
    }
}

// Solver implementation (basic)
// Solves LHS = RHS for var
static std::vector<std::shared_ptr<SymbolicExpr>> solve_internal(
    std::shared_ptr<SymbolicExpr> lhs, 
    std::shared_ptr<SymbolicExpr> rhs, 
    const std::string& var) 
{
    // Try to isolate var
    // Case: x = val
    if (lhs->type == SymbolicExpr::Type::Variable && lhs->identifier == var) {
        if (!has_variable(rhs, var)) return { rhs }; // Found solution
        return {}; // Implicit equation x = f(x), too hard for now
    }
    
    // Case: val = x
    if (rhs->type == SymbolicExpr::Type::Variable && rhs->identifier == var) {
        if (!has_variable(lhs, var)) return { lhs };
        return {};
    }
    
    // Case: u + v = rhs (assuming u has var, v doesn't)
    if (lhs->type == SymbolicExpr::Type::Add) {
        auto u = lhs->operands[0];
        auto v = lhs->operands[1];
        bool u_has = has_variable(u, var);
        bool v_has = has_variable(v, var);
        
        if (u_has && !v_has) {
            // u = rhs - v
            return solve_internal(u, SymbolicExpr::add(rhs, SymbolicExpr::multiply(SymbolicExpr::number(-1), v))->simplify(), var);
        }
        if (!u_has && v_has) {
            // v = rhs - u
            return solve_internal(v, SymbolicExpr::add(rhs, SymbolicExpr::multiply(SymbolicExpr::number(-1), u))->simplify(), var);
        }
    }
    
    // Case: u * v = rhs
    if (lhs->type == SymbolicExpr::Type::Multiply) {
        auto u = lhs->operands[0];
        auto v = lhs->operands[1];
        bool u_has = has_variable(u, var);
        bool v_has = has_variable(v, var);
        
        if (u_has && !v_has) {
            // u = rhs / v
            if (v->is_number()) {
                 // Optimization: numerical division
                 return solve_internal(u, SymbolicExpr::multiply(rhs, SymbolicExpr::number(v->convert_rational().reciprocal()))->simplify(), var);
            }
            return solve_internal(u, SymbolicExpr::multiply(rhs, SymbolicExpr::power(v, SymbolicExpr::number(-1)))->simplify(), var);
        }
        if (!u_has && v_has) {
             if (u->is_number()) {
                 return solve_internal(v, SymbolicExpr::multiply(rhs, SymbolicExpr::number(u->convert_rational().reciprocal()))->simplify(), var);
             }
             return solve_internal(v, SymbolicExpr::multiply(rhs, SymbolicExpr::power(u, SymbolicExpr::number(-1)))->simplify(), var);
        }
    }
    
    // Case u^n = rhs (n const) -> u = rhs^(1/n) (basic only, ignores negative roots etc)
    if (lhs->type == SymbolicExpr::Type::Power) {
        auto u = lhs->operands[0];
        auto n = lhs->operands[1];
        if (has_variable(u, var) && !has_variable(n, var) && n->is_number()) {
             return solve_internal(u, SymbolicExpr::power(rhs, SymbolicExpr::power(n, SymbolicExpr::number(-1)))->simplify(), var);
        }
    }
    
    return {};
}

/*
std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
    // Moved to symbolic_poly.cpp
    return {};
}
*/
