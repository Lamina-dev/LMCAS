#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <vector>
#include <string>

// Helper to check dependency
static bool has_variable(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr) return false;
    if (expr->type == SymbolicExpr::Type::Variable && expr->identifier == var) return true;
    for (const auto& op : expr->operands) {
        if (has_variable(op, var)) return true;
    }
    return false;
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

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
    // Assuming eq is an expression equal to 0
    return solve_internal(eq, SymbolicExpr::number(0), var_name);
}
