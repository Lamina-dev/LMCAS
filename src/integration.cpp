#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <variant>

// Note: SymbolicExpr, SymbolicNode, FunctionNode etc are in global namespace
// Integrator is in lamina namespace

namespace lamina {

using IntegralHistory = std::vector<SymbolicExpr>;

// Forward declarations of helper functions
static std::shared_ptr<SymbolicExpr> integrate_strategy(const SymbolicExpr& expr, const std::string& var, IntegralHistory& history, int depth);
static std::shared_ptr<SymbolicExpr> integrate_rec(const SymbolicExpr& expr, const std::string& var, IntegralHistory& history, int depth = 0);
static std::shared_ptr<SymbolicExpr> integrate_power_rule(const SymbolicExpr& expr, const std::string& var);
static std::shared_ptr<SymbolicExpr> integrate_basic_funcs(const SymbolicExpr& expr, const std::string& var);
static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_substitution(const SymbolicExpr& expr, const std::string& var);
static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_partial_fraction(const SymbolicExpr& expr, const std::string& var);
static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_heuristic_ibp(const SymbolicExpr& expr, const std::string& var, IntegralHistory& history, int depth);

// ----------------------------------------------------------------------
// Helper Functions (Local to this file)
// ----------------------------------------------------------------------

// make_expr_ptr: Helper to create a shared_ptr for SymbolicExpr
static std::shared_ptr<SymbolicExpr> make_expr_ptr(const SymbolicExpr& e) {
    return std::make_shared<SymbolicExpr>(e);
}

// valid_dependency: Check if expr depends on var (by checking if derivative is non-zero)
static bool valid_dependency(const SymbolicExpr& expr, const std::string& var) {
    auto diff = expr.differentiate(var);
    if (!diff) return false;
    auto simp_diff = diff->simplify();
    return !simp_diff->is_zero();
}

// sym_sub: Symbolic subtraction (a - b)
static std::shared_ptr<SymbolicExpr> sym_sub(const SymbolicExpr& a, const SymbolicExpr& b) {
    auto neg_b = SymbolicExpr::multiply(SymbolicExpr::number(-1), std::make_shared<SymbolicExpr>(b));
    return SymbolicExpr::add(std::make_shared<SymbolicExpr>(a), neg_b);
}

// sym_rational: Create a rational number SymbolicExpr
static std::shared_ptr<SymbolicExpr> sym_rational(long long num, long long den) {
    return SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
}

// make_arctan: Create an ArcTan function node manually since SymbolicExpr doesn't wrap it
static std::shared_ptr<SymbolicExpr> make_arctan(const std::shared_ptr<SymbolicExpr>& op) {
     return std::make_shared<SymbolicExpr>(
         std::make_shared<FunctionNode>(
             FunctionNode::FuncType::ArcTan, 
             std::vector<std::shared_ptr<SymbolicNode>>{op->root}
         )
     );
}

// has_integral_node_check: Recursively check if an expression contains an unresolved integral
static bool has_integral_node_check(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Calculus_Integral) return true;
        for (auto& arg : fn->arguments) {
            if (has_integral_node_check(arg)) return true;
        }
    } else if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands) if (has_integral_node_check(op)) return true;
    } else if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
         for (auto& op : mul->operands) if (has_integral_node_check(op)) return true;
    } else if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (has_integral_node_check(pow->base)) return true;
        if (has_integral_node_check(pow->exponent)) return true;
    }
    return false;
}

// ----------------------------------------------------------------------
// Main Integration Logic
// ----------------------------------------------------------------------

Integrator::Integrator() {}

SymbolicExpr Integrator::integrate(const SymbolicExpr& expr, const std::string& var_name) {
    // 1. Linearity: Integral(a*f + b*g) = a*Integral(f) + b*Integral(g)
    auto simp_expr = expr.simplify();
    
    
    // Check for constant multiply: Integral(c * f(x)) = c * Integral(f(x))
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> constants;
        std::vector<std::shared_ptr<SymbolicNode>> dependents;
        
        for (auto& op : mul->operands) {
            SymbolicExpr term(op);
            if (!valid_dependency(term, var_name)) {
                constants.push_back(op);
            } else {
                dependents.push_back(op);
            }
        }
        
        // Only split if we found constants and at least one dependent term
        if (!constants.empty() && dependents.size() < mul->operands.size()) {
            SymbolicExpr const_part = (constants.size() == 1) ? 
                SymbolicExpr(constants[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(constants));
                
            SymbolicExpr dep_part = (dependents.empty()) ? 
                *SymbolicExpr::number(1) : 
                ((dependents.size() == 1) ? SymbolicExpr(dependents[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(dependents)));
                
            auto int_part = integrate(dep_part, var_name);
            
            // Recombine: const * integral(dep)
            return *(SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(const_part), std::make_shared<SymbolicExpr>(int_part)));
        }
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> results;
        for (auto& op : add->operands) {
            SymbolicExpr term(op);
            auto int_term = integrate(term, var_name);
            results.push_back(int_term.root); // Extract root
        }
        return SymbolicExpr(std::make_shared<AddNode>(results));
    }

    // Delegate to recursive handler
    IntegralHistory history;
    return *integrate_rec(*simp_expr, var_name, history);
}

// Check for singularity in [lower, upper]
// Return value: detected singularity point (or nullptr if none found)
std::shared_ptr<SymbolicExpr> check_singularity(const SymbolicExpr& expr, const std::string& var_name, 
                                                const SymbolicExpr& lower, const SymbolicExpr& upper) {
    // Simple check: Look for denominator roots
    // Or argument checks for specific functions (ln(x) -> x=0)
    // Current heuristic: Denominator 0
    
    // Decompose into Num / Den
    std::shared_ptr<SymbolicExpr> den_expr = nullptr;

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        std::vector<std::shared_ptr<SymbolicNode>> den_nodes;
        for(const auto& op : mul->operands) {
             if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                 if (auto sn = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                      if ((std::holds_alternative<double>(sn->value) && std::get<double>(sn->value) < 0) ||
                          (std::holds_alternative<BigInt>(sn->value) && std::get<BigInt>(sn->value).to_double() < 0) ||
                          (std::holds_alternative<Rational>(sn->value) && std::get<Rational>(sn->value).to_double() < 0)) {
                          den_nodes.push_back(pow->base);
                      }
                 }
             }
        }
        if (!den_nodes.empty()) {
            den_expr = std::make_shared<SymbolicExpr>(
                den_nodes.size() == 1 ? den_nodes[0] : std::make_shared<MultiplyNode>(den_nodes)
            );
        }
    } else if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
         // Single power term check
         if (auto sn = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
              if ((std::holds_alternative<double>(sn->value) && std::get<double>(sn->value) < 0) ||
                  (std::holds_alternative<BigInt>(sn->value) && std::get<BigInt>(sn->value).to_double() < 0) ||
                  (std::holds_alternative<Rational>(sn->value) && std::get<Rational>(sn->value).to_double() < 0)) {
                  den_expr = std::make_shared<SymbolicExpr>(pow->base);
              }
         }
    }

    if (den_expr) {
        // Solve den(x) = 0
        // Currently we only check if den(0) == 0 or den(lower/upper) == 0?
        // Let's assume linear factors x-c for now.
        if (auto v = std::dynamic_pointer_cast<VariableNode>(den_expr->root)) {
            if (v->name == var_name) return SymbolicExpr::number(0); // x=0 is singularity
        }
        
        // Check if [lower, upper] constains 0
        // We need numerical check for [lower, upper]
        double l = lower.to_double();
        double u = upper.to_double();
        
        if (l <= 0 && u >= 0 && l != u) {
             // 0 is inside
             // Return 0
             // But only if den actually has root 0. 
             // We assumed den IS 'x' above.
        }
    }
    
    // Hardcoded check for 1/x type
    // If expr is x^-1
    // We can rely on user to be sensible OR implement solve(den=0)
    // Let's implement a very basic Solver call? 
    // Or just check if 0 is in the interval for 1/x type integrals.
    
    return nullptr; 
}


SymbolicExpr Integrator::integrate_def(const SymbolicExpr& expr, const std::string& var_name, 
                              const SymbolicExpr& lower, const SymbolicExpr& upper) {
    // 0. Check for singularities (Initial attempt for 1/x at 0)
    // Heuristic: Check if integrand is x^-1 (or c*x^-1) and interval crosses 0
    // This is very limited but serves the "1/x from -1 to 1" example.
    
    SymbolicExpr simp_expr_val = *expr.simplify();
    bool is_inv_x = false;
    
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(simp_expr_val.root)) {
        if (auto v = std::dynamic_pointer_cast<VariableNode>(pow->base)) {
             if (v->name == var_name) {
                 if (auto en = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                     if ((std::holds_alternative<double>(en->value) && std::abs(std::get<double>(en->value) + 1.0) < 1e-9) ||
                         (std::holds_alternative<BigInt>(en->value) && std::get<BigInt>(en->value).to_int() == -1) || 
                         (std::holds_alternative<Rational>(en->value) && std::get<Rational>(en->value).to_double() == -1.0)) {
                         is_inv_x = true;
                     }
                 }
             }
        }
    }
    
    // Check bounds
    // We need numerical values for bounds to check if 0 is inside
    // If symbols, we can't reliably know without constraints solver.
    double l_val = lower.to_double();
    double u_val = upper.to_double();
    
    // Check if bounds are purely numeric (non-zero implies success of to_double usually, but 0 is valid too)
    // If to_double() returns 0 for a symbol 'a', this logic is flawed for symbols.
    // So we check if roots are NumberNodes.
    
    bool numeric_bounds = (lower.root && std::dynamic_pointer_cast<NumberNode>(lower.root)) && 
                          (upper.root && std::dynamic_pointer_cast<NumberNode>(upper.root));

    if (is_inv_x && numeric_bounds) {
         if (l_val < -1e-9 && u_val > 1e-9) {
             // Split around 0
             // Limit_{t->0-} Int(l, t) + Limit_{t->0+} Int(t, u)
             
             auto t = std::make_shared<SymbolicExpr>(*SymbolicExpr::variable("t"));
             auto zero = std::make_shared<SymbolicExpr>(*SymbolicExpr::number(0));
             
             // Int(l, t) = ln(|t|) - ln(|l|)
             // Limit t->0- of ln(|t|) -> -Inf
             
             // Int(t, u) = ln(|u|) - ln(|t|)
             // Limit t->0+ of -ln(|t|) -> -(-Inf) = +Inf
             
             // -Inf + Inf -> Indeterminate / Divergent (Cauchy Principal Value might exist = 0)
             // Standard Riemann integral diverges.
             
             // Our limit function returns Infinity node.
             
             auto int_left = integrate_def(expr, var_name, lower, *t);
             auto lim_left = int_left.limit("t", zero, "-");
             
             auto int_right = integrate_def(expr, var_name, *t, upper);
             auto lim_right = int_right.limit("t", zero, "+");
             
             if (lim_left && lim_right) {
                  return *SymbolicExpr::add(lim_left, lim_right);
             }
         }
    }

    // 1. Calculate Indefinite Integral
    SymbolicExpr indefinite = integrate(expr, var_name);
    
    // 2. Check if integration failed (returns a symbolic Integral node)
    // Note: The recursive integrator returns an integral node if it fails.
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(indefinite.root)) {
        if (func->type == FunctionNode::FuncType::Calculus_Integral) {
             // Return symbolic definite integral node: Int(f, x, a, b)
             std::vector<std::shared_ptr<SymbolicNode>> args;
             args.push_back(expr.root);
             args.push_back(SymbolicExpr::variable(var_name)->root);
             args.push_back(lower.root);
             args.push_back(upper.root); 
             return SymbolicExpr(std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
        }
    }

    // 3. Apply Fundamental Theorem of Calculus: F(b) - F(a)
    // We attempt simple substitution. 
    // Ideally we should use limits (lim_{x->b} F(x)) to handle singularities,
    // but simple substitution covers standard cases.
    
    auto F_b = indefinite.substitute(var_name, make_expr_ptr(upper));
    auto F_a = indefinite.substitute(var_name, make_expr_ptr(lower));
    
    auto result = sym_sub(*F_b, *F_a);
    
    return *result->simplify();
}


// Recursive integration wrapper with cycle detection
static std::shared_ptr<SymbolicExpr> integrate_rec(const SymbolicExpr& expr, const std::string& var_name, IntegralHistory& history, int depth) {
    // Hard recursion limit check
    if (depth > 8) {
        std::vector<std::shared_ptr<SymbolicNode>> args;
        args.push_back(expr.root);
        args.push_back(SymbolicExpr::variable(var_name)->root);
        return std::make_shared<SymbolicExpr>(
            std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args)
        );
    }

    // Cycle check
    for (size_t i = 0; i < history.size(); ++i) {
        // Try to divide current expr by history expr
        auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), make_expr_ptr(history[i]))->simplify();
        
        // If ratio is constant w.r.t var_name, we found a cycle
        if (!valid_dependency(*ratio, var_name)) {
            // Return k * INT_CYCLE_i
            return SymbolicExpr::multiply(ratio, SymbolicExpr::variable("INT_CYCLE_" + std::to_string(i)));
        }
    }

    history.push_back(expr);
    size_t my_idx = history.size() - 1;
    
    // Pass depth + 1 to strategy? No, strategy decides which sub-calls increment depth.
    // Actually, strategy just passes 'depth' along, but if strategy calls integrate_rec recursively, it should increment.
    // EXCEPT: integrate_strategy itself is not recursive, it calls helpers.
    // The helpers (like IBP) call integrate_rec with depth + 1.
    auto res = integrate_strategy(expr, var_name, history, depth);
    
    history.pop_back();

    if (res) {
        std::string cycle_var = "INT_CYCLE_" + std::to_string(my_idx);
        // Check if result contains the cycle variable
        if (valid_dependency(*res, cycle_var)) {
             // Solve for cycle_var: I = A + B*I  => I = A / (1-B)
             // Here res is (A + B*I).
             // Wait, usually res is linear in I? Yes.
             
             auto B = res->differentiate(cycle_var)->simplify();
             auto A = res->substitute(cycle_var, SymbolicExpr::number(0))->simplify();
             
             auto one_minus_B = sym_sub(*SymbolicExpr::number(1), *B)->simplify();
             if (!one_minus_B->is_zero()) {
                 res = SymbolicExpr::divide(A, one_minus_B);
             }
        }
    }
    return res;
}

// Actual integration strategy
static std::shared_ptr<SymbolicExpr> integrate_strategy(const SymbolicExpr& expr, const std::string& var_name, IntegralHistory& history, int depth) {
    // 0. Check if constant (Integral(c dx) = c*x)
    if (!valid_dependency(expr, var_name)) {
        return SymbolicExpr::multiply(make_expr_ptr(expr), SymbolicExpr::variable(var_name));
    }
    
    // 1. Basic Power Rule
    auto pr_res = integrate_power_rule(expr, var_name);
    if (pr_res) return pr_res;

    // 2. Basic Transcendental Functions (sin, cos, exp, etc.)
    auto bf_res = integrate_basic_funcs(expr, var_name);
    if (bf_res) return bf_res;

    // 3. Substitution
    auto sub_res = try_substitution(expr, var_name);
    if (sub_res.first) return sub_res.second;

    // 4. Partial Fractions
    auto pf_res = try_partial_fraction(expr, var_name);
    if (pf_res.first) return pf_res.second;

    // 5. Integration by Parts (Heuristic)
    auto ibp_res = try_heuristic_ibp(expr, var_name, history, depth);
    if (ibp_res.first) return ibp_res.second;

    // Fallback: Return unintegrated Integral form
    // Create Integral(expr, var) node
    std::vector<std::shared_ptr<SymbolicNode>> args;
    args.push_back(expr.root);
    args.push_back(SymbolicExpr::variable(var_name)->root);
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args)
    );
}

// ----------------------------------------------------------------------
// Specific Integration Strategies
// ----------------------------------------------------------------------

static std::shared_ptr<SymbolicExpr> integrate_power_rule(const SymbolicExpr& expr, const std::string& var) {
    // Case 1: x -> x^2/2
    if (auto v_node = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        if (v_node->name == var) {
             return SymbolicExpr::multiply(
                SymbolicExpr::power(make_expr_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2)
             );
        }
    }

    // Case 2: x^n
    if (auto p_node = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        SymbolicExpr base(p_node->base);
        SymbolicExpr exp(p_node->exponent);

        if (auto b_var = std::dynamic_pointer_cast<VariableNode>(base.root)) {
            if (b_var->name == var && !valid_dependency(exp, var)) {
                // x^n -> x^(n+1)/(n+1)
                auto n_plus_1 = SymbolicExpr::add(make_expr_ptr(exp), SymbolicExpr::number(1))->simplify();
                
                if (n_plus_1->is_zero()) { // n = -1 -> Integral(1/x) = ln(x)
                    return SymbolicExpr::ln(make_expr_ptr(base)); 
                }
                
                return SymbolicExpr::divide(
                    SymbolicExpr::power(make_expr_ptr(base), n_plus_1),
                    n_plus_1
                );
            }
        }
    }
    return nullptr; 
}

static std::shared_ptr<SymbolicExpr> integrate_basic_funcs(const SymbolicExpr& expr, const std::string& var) {
    // Special case for 1/x -> ln(x) handled in power rule usually, but check just in case
    // Handle forms like x^-1 if not caught by power rule? (Should be caught)

    // Handle standard functions
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        // We only handle simple f(x) here. f(2x) etc should be handled by substitution.
        if (func->arguments.empty()) return nullptr;
        SymbolicExpr arg(func->arguments[0]);
        
        // Ensure argument is exactly the variable
        if (arg.to_string() != var) return nullptr; 
        
        if (func->type == FunctionNode::FuncType::Sin) {
            return SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(arg)));
        }
        if (func->type == FunctionNode::FuncType::Cos) {
            return SymbolicExpr::sin(make_expr_ptr(arg));
        }
        if (func->type == FunctionNode::FuncType::Tan) {
             return SymbolicExpr::multiply(SymbolicExpr::number(-1), 
                SymbolicExpr::ln(SymbolicExpr::cos(make_expr_ptr(arg))));
        }
        if (func->type == FunctionNode::FuncType::Exp) {
            return make_expr_ptr(expr); // exp(x) -> exp(x)
        }
    }
    return nullptr;
}

static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_substitution(const SymbolicExpr& expr, const std::string& var) {
    // Attempt simple u-substitutions
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        ops.push_back(expr.root);
    }
    
    // Try each factor as the potential "function of u" where u is inside it
    for (size_t i = 0; i < ops.size(); ++i) {
        SymbolicExpr candidate_term(ops[i]);
        SymbolicExpr u;
        bool possible = false;

        if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
            u = SymbolicExpr(pow->base);
            possible = true;
        } 
        else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
            if (!func->arguments.empty()) {
                u = SymbolicExpr(func->arguments[0]);
                possible = true;
            }
        }

        if (possible && valid_dependency(u, var)) {
            // Check if du exists in the expression
            auto d_ptr = u.differentiate(var);
            if (!d_ptr) continue;
            auto du = d_ptr->simplify();
            if (du->is_zero()) continue;

            // Compute "remainder" = expr / (f(u) * du)
            // Ideally expr is f(u) * du or f(u) * du * constant
            
            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(make_expr_ptr(f_u), make_expr_ptr(*du));
            
            // Check if expr is divisible by term_times_du up to a constant
            auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), term_times_du)->simplify();
            
            // If the ratio does not depend on var, substitution is valid!
            if (!valid_dependency(*ratio, var)) {
                // Integrate f(u) du -> F(u)
                std::shared_ptr<SymbolicExpr> prim = nullptr;
                
                if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
                    SymbolicExpr n(pow->exponent);
                    // u^n du -> u^(n+1)/(n+1)
                    auto np1 = SymbolicExpr::add(make_expr_ptr(n), SymbolicExpr::number(1))->simplify();
                    if (np1->is_zero()) {
                         prim = SymbolicExpr::ln(make_expr_ptr(u));
                    } else {
                         prim = SymbolicExpr::divide(
                            SymbolicExpr::power(make_expr_ptr(u), np1),
                            np1
                         );
                    }
                }
                else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
                    if (func->type == FunctionNode::FuncType::Cos) {
                         prim = SymbolicExpr::sin(make_expr_ptr(u));
                    } else if (func->type == FunctionNode::FuncType::Sin) {
                         prim = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
                    } else if (func->type == FunctionNode::FuncType::Exp) {
                         prim = SymbolicExpr::exp(make_expr_ptr(u));
                    }
                }
                
                if (prim) {
                    // Result = constant_ratio * F(u)
                    return {true, SymbolicExpr::multiply(ratio, prim)};
                }
            }
        }
    }
    return {false, nullptr};
}


static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_partial_fraction(const SymbolicExpr& expr, const std::string& var) {
// Handling 1/(ax^2+bx+c) or like (ax^2+bx+c)^-1
    std::shared_ptr<SymbolicExpr> den = nullptr;
    
    if (auto p = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
         double exp_val = 0;
         bool is_inv = false;
         try {
             if (auto num_node = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
                  if(std::holds_alternative<double>(num_node->value)) exp_val = std::get<double>(num_node->value);
                  else if(std::holds_alternative<Rational>(num_node->value)) exp_val = std::get<Rational>(num_node->value).to_double();
                  else if(std::holds_alternative<BigInt>(num_node->value)) exp_val = std::get<BigInt>(num_node->value).to_double();
                  if (std::abs(exp_val + 1.0) < 1e-9) is_inv = true;
             }
         } catch(...) {}
         
         if (is_inv) {
             den = make_expr_ptr(SymbolicExpr(p->base));
         }
    }
    
    if (!den) return {false, nullptr};

    try {
        // Use lamina::Polynomial
        lamina::Polynomial<lamina::SymbolicPolyCoeff> Q = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(make_expr_ptr(*den), var);
        
        // Only handling quadratics for now
        if (Q.degree() == 2) {
            SymbolicExpr c_expr = *(Q.coeffs[0].val); // x^0
            SymbolicExpr b_expr = *(Q.coeffs[1].val); // x^1
            SymbolicExpr a_expr = *(Q.coeffs[2].val); // x^2
            
            double a = a_expr.is_number() ? a_expr.to_double() : 1.0; 
            double b = b_expr.is_number() ? b_expr.to_double() : 0.0;
            double c = c_expr.is_number() ? c_expr.to_double() : 0.0;
            
            if (std::abs(a) < 1e-9) return {false, nullptr};
            
            double delta = b*b - 4*a*c;
            
            if (delta > 1e-9) {
                // Positive discriminant: two real roots: 1 / (a(x-r1)(x-r2))
                double sqrt_delta = std::sqrt(delta);
                
                // Result formula: (1/sqrt(delta)) * ln |(2ax + b - sqrt(delta))/(2ax + b + sqrt(delta))|
                // Simplified version
                
                auto scalar = SymbolicExpr::number(1.0 / sqrt_delta);
                auto two_a = SymbolicExpr::number(2.0*a);
                auto b_num = SymbolicExpr::number(b);
                
                auto two_a_x = SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var));
                auto two_a_x_plus_b = SymbolicExpr::add(make_expr_ptr(*two_a_x), make_expr_ptr(*b_num));
                
                auto term1_arg = sym_sub(*two_a_x_plus_b, *SymbolicExpr::number(sqrt_delta));
                auto term2_arg = SymbolicExpr::add(make_expr_ptr(*two_a_x_plus_b), SymbolicExpr::number(sqrt_delta));
                
                auto term1 = SymbolicExpr::ln(make_expr_ptr(*term1_arg));
                auto term2 = SymbolicExpr::ln(make_expr_ptr(*term2_arg));
                
                return {true, SymbolicExpr::multiply(scalar, sym_sub(*term1, *term2))};
                
            } else if (delta < -1e-9) {
                // Negative discriminant: ArcTan
                double neg_delta = -delta;
                double sqrt_neg_delta = std::sqrt(neg_delta);
                
                // Result: (2/sqrt(-D)) * arctan((2ax+b)/sqrt(-D))
                
                auto scalar = SymbolicExpr::number(2.0 / sqrt_neg_delta);
                auto two_a = SymbolicExpr::number(2.0*a);
                auto b_num = SymbolicExpr::number(b);
                
                auto num = SymbolicExpr::add(
                    SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var)),
                    make_expr_ptr(*b_num)
                );
                auto inner = SymbolicExpr::divide(make_expr_ptr(*num), SymbolicExpr::number(sqrt_neg_delta));
                
                return {true, SymbolicExpr::multiply(scalar, make_arctan(inner))};
            }
        }
    } catch(...) {
        // Fallback or conversion error
    }

    return {false, nullptr}; 
}


// Heuristic IBP
static std::pair<bool, std::shared_ptr<SymbolicExpr>> try_heuristic_ibp(const SymbolicExpr& expr, const std::string& var, IntegralHistory& history, int depth) {
     // IBP: Integral(u dv) = uv - Integral(v du)
     std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        // Treat as ln(x) * 1
        ops.push_back(expr.root);
    }
    
    // Strategy: Identify u using LIATE
    // L-I-A-T-E rule: Log, Inverse trig, Algebraic, Trig, Exponential.
    // We pick "u" as the one that appears earlier in LIATE.
    
    int best_u_idx = -1;
    int best_score = 100; // lower is better
    
    auto get_score = [&](const std::shared_ptr<SymbolicNode>& node) -> int {
        SymbolicExpr e(node);
        if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
            if (fn->type == FunctionNode::FuncType::Ln || fn->type == FunctionNode::FuncType::Log) return 1;
            if (fn->type == FunctionNode::FuncType::ArcSin || fn->type == FunctionNode::FuncType::ArcTan) return 2;
            if (fn->type == FunctionNode::FuncType::Sin || fn->type == FunctionNode::FuncType::Cos) return 4;
            if (fn->type == FunctionNode::FuncType::Exp) return 5;
        }
        
        if (!valid_dependency(e, var)) return 10;
        
        // Variables and Powers generally algebraic
        if (std::dynamic_pointer_cast<VariableNode>(node)) return 3;
        if (std::dynamic_pointer_cast<PowerNode>(node)) return 3; 

        return 10;
    };
    
    // Scan for best u
    if (ops.size() == 1) {
        int s = get_score(ops[0]);
        if (s <= 2) { 
             best_u_idx = 0; // Only try IBP if it's Log or Inverse Trig alone
        }
    } else {
        for(size_t i=0; i<ops.size(); ++i) {
            // u must depend on var? Not strictly, but typically yes.
            if (!valid_dependency(SymbolicExpr(ops[i]), var)) continue;
            int s = get_score(ops[i]);
            if (s < best_score) {
                best_score = s;
                best_u_idx = (int)i;
            }
        }
    }
    
    if (best_u_idx == -1) return {false, nullptr};
    
    // Assign u
    SymbolicExpr u(ops[best_u_idx]);
    
    // Assign dv (everything else)
    std::vector<std::shared_ptr<SymbolicNode>> dv_ops;
    for(size_t i=0; i<ops.size(); ++i) {
        if ((int)i != best_u_idx) dv_ops.push_back(ops[i]);
    }
    
    std::shared_ptr<SymbolicExpr> dv;
    if (dv_ops.empty()) dv = SymbolicExpr::number(1);
    else if (dv_ops.size() == 1) dv = make_expr_ptr(SymbolicExpr(dv_ops[0]));
    else dv = std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(dv_ops));
    
    // V = Integral(dv)
    // Avoid infinite recursion: if dv is complex, we might loop. 
    // Ideally we want dv to be "easier" to integrate.
    auto v = integrate_rec(*dv, var, history, depth + 1);
    
    // If v is still an integral, IBP failed to simplify start
    if (has_integral_node_check(v->root)) return {false, nullptr};

    // du = differentiate(u)
    auto du_ptr = u.differentiate(var);
    if (!du_ptr) return {false, nullptr};
    auto du = make_expr_ptr(*du_ptr->simplify());

    // Result = uv - Integral(v du)
    auto uv = SymbolicExpr::multiply(make_expr_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    auto int_vdu = integrate_rec(*vdu, var, history, depth + 1);
    
    // If residuals still have integrals, we permit it (it might be simpler)
    return {true, sym_sub(*uv, *int_vdu)};
}

} // namespace lamina
