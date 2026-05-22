// solve_transcendental.cpp - Transcendental equation solver (inverse function method)
// Implements inverse-function pattern matching for sin, cos, tan, exp, ln, Lambert W, and a^g(x)

#include "solve_transcendental.hpp"
#include "solve_strategies.hpp"
#include "poly_utils.hpp"
#include "lmmc/config.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace lamina {

// ============================================================================
// Helper utilities
// ============================================================================

// Create an arcsin expression: arcsin(c)
static std::shared_ptr<SymbolicExpr> make_arcsin(const std::shared_ptr<SymbolicExpr>& c) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::ArcSin,
            std::vector<std::shared_ptr<SymbolicNode>>{c->root}));
}

// Create an arccos expression: arccos(c)
static std::shared_ptr<SymbolicExpr> make_arccos(const std::shared_ptr<SymbolicExpr>& c) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<SymbolicNode>>{c->root}));
}

// Create an arctan expression: arctan(c)
static std::shared_ptr<SymbolicExpr> make_arctan(const std::shared_ptr<SymbolicExpr>& c) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::ArcTan,
            std::vector<std::shared_ptr<SymbolicNode>>{c->root}));
}

// Create a numeric constant for pi
static std::shared_ptr<SymbolicExpr> make_pi() {
    return SymbolicExpr::number(LMMC_CONST_PI);
}

// Try to evaluate an expression as a numeric constant.
// Returns true if the expression is a pure number (no variables), false otherwise.
static bool try_evaluate_numeric(const std::shared_ptr<SymbolicExpr>& expr, lmmc_real_t& out) {
    if (!expr || !expr->root) return false;
    // Check if expression contains any variables
    struct VarChecker : public SymbolicVisitor {
        bool has_var = false;
        void visit(NumberNode&) override {}
        void visit(VariableNode&) override { has_var = true; }
        void visit(AddNode& n) override { for (auto& op : n.operands) { if (has_var) return; op->accept(*this); } }
        void visit(MultiplyNode& n) override { for (auto& op : n.operands) { if (has_var) return; op->accept(*this); } }
        void visit(PowerNode& n) override { n.base->accept(*this); if (!has_var) n.exponent->accept(*this); }
        void visit(FunctionNode& n) override { for (auto& arg : n.arguments) { if (has_var) return; arg->accept(*this); } }
        void visit(MatrixNode&) override {}
    } checker;
    expr->root->accept(checker);
    if (checker.has_var) return false;
    
    out = expr->to_numeric();
    return true;
}

// ============================================================================
// Pattern decomposition: expr = f(g(x)) - c  where c does not depend on var
// Given an expression (which equals 0), try to decompose it as:
//   transcendental_func(inner_expr) + constant_part = 0
// i.e., transcendental_func(inner_expr) = -constant_part
// ============================================================================

// Structure to hold a decomposed pattern: func(g) = c
struct InversePattern {
    FunctionNode::FuncType func_type;
    std::shared_ptr<SymbolicExpr> inner;  // g(x) - the argument of the function
    std::shared_ptr<SymbolicExpr> rhs;    // c - the constant value
};

// Try to decompose expr into: func(inner) - c = 0
// The expression is an AddNode with a function call and constant terms.
// Returns nullopt if no pattern is found.
static std::optional<InversePattern> decompose_trig_exp_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    
    if (!expr || !expr->root) return std::nullopt;
    
    // Case 1: expr is directly a function call (func(g) = 0, so c = 0)
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr->root)) {
        if (func->arguments.size() == 1 && depends_on_var(func->arguments[0], var)) {
            auto ft = func->type;
            if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                ft == FunctionNode::FuncType::Ln) {
                return InversePattern{
                    ft,
                    std::make_shared<SymbolicExpr>(func->arguments[0]),
                    SymbolicExpr::number(0)
                };
            }
        }
    }
    
    // Case 2: expr is an AddNode: look for a function term + constant terms
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr->root)) {
        // Separate terms into: function terms (containing var) and constant terms
        std::shared_ptr<FunctionNode> func_term = nullptr;
        std::shared_ptr<SymbolicNode> func_coeff = nullptr; // coefficient of the function term
        std::vector<std::shared_ptr<SymbolicNode>> const_terms;
        
        for (auto& op : add->operands) {
            // Check if this operand is a transcendental function of var
            if (auto f = std::dynamic_pointer_cast<FunctionNode>(op)) {
                if (f->arguments.size() == 1 && depends_on_var(f->arguments[0], var)) {
                    auto ft = f->type;
                    if ((ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                         ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                         ft == FunctionNode::FuncType::Ln) && !func_term) {
                        func_term = f;
                        func_coeff = std::make_shared<NumberNode>(BigInt(1));
                        continue;
                    }
                }
            }
            // Check if this is coeff * func(g(x)) - e.g. MultiplyNode with a number and a function
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                std::shared_ptr<FunctionNode> mf = nullptr;
                std::vector<std::shared_ptr<SymbolicNode>> coeff_parts;
                for (auto& mop : mul->operands) {
                    if (auto f = std::dynamic_pointer_cast<FunctionNode>(mop)) {
                        if (f->arguments.size() == 1 && depends_on_var(f->arguments[0], var)) {
                            auto ft = f->type;
                            if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                                ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                                ft == FunctionNode::FuncType::Ln) {
                                if (!mf) { mf = f; continue; }
                            }
                        }
                    }
                    coeff_parts.push_back(mop);
                }
                if (mf && !func_term && !depends_on_var(SymbolicFactory::create_multiply(coeff_parts), var)) {
                    func_term = mf;
                    func_coeff = SymbolicFactory::create_multiply(coeff_parts);
                    continue;
                }
            }
            
            // Otherwise it's a constant term (or a term we can't handle)
            if (!depends_on_var(op, var)) {
                const_terms.push_back(op);
            } else {
                // Contains var but isn't a simple transcendental function - can't decompose
                return std::nullopt;
            }
        }
        
        if (func_term) {
            // We have: coeff * func(inner) + sum(const_terms) = 0
            // => func(inner) = -sum(const_terms) / coeff
            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = std::make_shared<SymbolicExpr>(const_terms[0]);
            } else {
                const_sum = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(const_terms));
            }
            // rhs = -const_sum / coeff
            auto neg_const = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1));
            auto coeff_expr = std::make_shared<SymbolicExpr>(func_coeff);
            auto rhs = SymbolicExpr::divide(neg_const, coeff_expr)->simplify();
            
            return InversePattern{
                func_term->type,
                std::make_shared<SymbolicExpr>(func_term->arguments[0]),
                rhs
            };
        }
    }
    
    // Case 3: expr is coeff * func(g(x)) (MultiplyNode at top level, c = 0)
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr->root)) {
        for (auto& op : mul->operands) {
            if (auto f = std::dynamic_pointer_cast<FunctionNode>(op)) {
                if (f->arguments.size() == 1 && depends_on_var(f->arguments[0], var)) {
                    auto ft = f->type;
                    if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                        ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                        ft == FunctionNode::FuncType::Ln) {
                        return InversePattern{
                            ft,
                            std::make_shared<SymbolicExpr>(f->arguments[0]),
                            SymbolicExpr::number(0)
                        };
                    }
                }
            }
        }
    }
    
    return std::nullopt;
}

// ============================================================================
// Lambert W pattern: g(x) * exp(g(x)) = c
// Detect: g*exp(g) - c = 0 in the expression
// ============================================================================

struct LambertWPattern {
    std::shared_ptr<SymbolicExpr> inner;  // g(x)
    std::shared_ptr<SymbolicExpr> rhs;    // c
};

// Check if two AST nodes are structurally equal
static bool nodes_equal(const std::shared_ptr<SymbolicNode>& a, const std::shared_ptr<SymbolicNode>& b) {
    if (!a || !b) return a == b;
    return a->equals(*b);
}

static std::optional<LambertWPattern> decompose_lambert_w_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    
    if (!expr || !expr->root) return std::nullopt;
    
    // Look for pattern: g(x) * exp(g(x)) + const_terms = 0
    // In the AST this could be:
    //   AddNode[ MultiplyNode[g, exp(g)], const_terms... ]
    //   or just MultiplyNode[g, exp(g)] (when c = 0)
    
    auto check_mul_is_g_exp_g = [&](const std::shared_ptr<MultiplyNode>& mul)
        -> std::optional<std::shared_ptr<SymbolicExpr>> {
        // Look for two factors: one is exp(something), the other equals that something
        for (size_t i = 0; i < mul->operands.size(); ++i) {
            if (auto f = std::dynamic_pointer_cast<FunctionNode>(mul->operands[i])) {
                if (f->type == FunctionNode::FuncType::Exp && f->arguments.size() == 1) {
                    // Found exp(something). Check if remaining factors equal that something.
                    auto exp_arg = f->arguments[0];
                    std::vector<std::shared_ptr<SymbolicNode>> remaining;
                    for (size_t j = 0; j < mul->operands.size(); ++j) {
                        if (j != i) remaining.push_back(mul->operands[j]);
                    }
                    std::shared_ptr<SymbolicNode> g_node;
                    if (remaining.size() == 1) {
                        g_node = remaining[0];
                    } else {
                        g_node = std::make_shared<MultiplyNode>(remaining);
                    }
                    
                    if (nodes_equal(g_node, exp_arg) && depends_on_var(g_node, var)) {
                        return std::make_shared<SymbolicExpr>(g_node);
                    }
                }
            }
        }
        return std::nullopt;
    };
    
    // Case 1: expr is directly MultiplyNode (g*exp(g) = 0)
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr->root)) {
        auto g = check_mul_is_g_exp_g(mul);
        if (g) {
            return LambertWPattern{*g, SymbolicExpr::number(0)};
        }
    }
    
    // Case 2: expr is AddNode with a g*exp(g) term and constant terms
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr->root)) {
        std::optional<std::shared_ptr<SymbolicExpr>> g_found;
        std::vector<std::shared_ptr<SymbolicNode>> const_terms;
        bool has_other_var_terms = false;
        
        for (auto& op : add->operands) {
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                auto g = check_mul_is_g_exp_g(mul);
                if (g && !g_found) {
                    g_found = g;
                    continue;
                }
            }
            if (!depends_on_var(op, var)) {
                const_terms.push_back(op);
            } else {
                has_other_var_terms = true;
            }
        }
        
        if (g_found && !has_other_var_terms) {
            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = std::make_shared<SymbolicExpr>(const_terms[0]);
            } else {
                const_sum = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(const_terms));
            }
            auto rhs = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1))->simplify();
            return LambertWPattern{*g_found, rhs};
        }
    }
    
    return std::nullopt;
}

// ============================================================================
// Exponential base pattern: a^g(x) = c  (where a is a constant, a>0, a!=1)
// Detect: a^g(x) - c = 0 in the expression
// This is represented as PowerNode(a, g(x)) where a is constant
// ============================================================================

struct ExpBasePattern {
    std::shared_ptr<SymbolicExpr> base;   // a (constant)
    std::shared_ptr<SymbolicExpr> inner;  // g(x)
    std::shared_ptr<SymbolicExpr> rhs;    // c
};

static std::optional<ExpBasePattern> decompose_exp_base_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    
    if (!expr || !expr->root) return std::nullopt;
    
    auto check_power_is_a_to_g = [&](const std::shared_ptr<PowerNode>& pow)
        -> std::optional<std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>> {
        // a^g(x): base is constant, exponent depends on var
        if (!depends_on_var(pow->base, var) && depends_on_var(pow->exponent, var)) {
            auto base_expr = std::make_shared<SymbolicExpr>(pow->base);
            auto exp_expr = std::make_shared<SymbolicExpr>(pow->exponent);
            // Verify base is a positive constant != 1
            lmmc_real_t base_val;
            if (try_evaluate_numeric(base_expr, base_val)) {
                if (base_val > 0 && std::abs(base_val - 1.0) > LMMC_REAL_EPSILON) {
                    return std::make_pair(base_expr, exp_expr);
                }
            }
        }
        return std::nullopt;
    };
    
    // Case 1: expr is directly PowerNode (a^g(x) = 0 - but a^g > 0 always, so no solution)
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr->root)) {
        auto result = check_power_is_a_to_g(pow);
        if (result) {
            // a^g(x) = 0 has no solution (a^g > 0 for a > 0)
            return ExpBasePattern{result->first, result->second, SymbolicExpr::number(0)};
        }
    }
    
    // Case 2: expr is AddNode with a^g(x) + const_terms = 0
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr->root)) {
        std::optional<std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>> found;
        std::vector<std::shared_ptr<SymbolicNode>> const_terms;
        bool has_other_var_terms = false;
        
        for (auto& op : add->operands) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                auto result = check_power_is_a_to_g(pow);
                if (result && !found) {
                    found = result;
                    continue;
                }
            }
            if (!depends_on_var(op, var)) {
                const_terms.push_back(op);
            } else {
                has_other_var_terms = true;
            }
        }
        
        if (found && !has_other_var_terms) {
            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = std::make_shared<SymbolicExpr>(const_terms[0]);
            } else {
                const_sum = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(const_terms));
            }
            auto rhs = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1))->simplify();
            return ExpBasePattern{found->first, found->second, rhs};
        }
    }
    
    return std::nullopt;
}


// ============================================================================
// Inversion logic: given a pattern, produce the inverted g(x) = value expressions
// ============================================================================

// For sin(g) = c: g = arcsin(c) + 2k*pi, g = pi - arcsin(c) + 2k*pi
static std::vector<std::shared_ptr<SymbolicExpr>> invert_sin(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {
    
    // Domain check: |c| <= 1
    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (std::abs(c_val) > 1.0 + LMMC_REAL_EPSILON) {
            return {};  // Domain violation
        }
    }
    
    auto pi = make_pi();
    auto asin_c = make_arcsin(c);
    
    // Determine how many periods to enumerate
    int num_periods = 1; // default: k=0 only (principal period)
    if (max_roots > 0) {
        // Each period gives 2 solutions, so num_periods = ceil(max_roots / 2)
        num_periods = (max_roots + 1) / 2;
    }
    
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {
        // g = arcsin(c) + 2k*pi
        auto two_k_pi = SymbolicExpr::multiply(SymbolicExpr::number(2 * k), pi);
        auto sol1 = SymbolicExpr::add(asin_c, two_k_pi)->simplify();
        results.push_back(sol1);
        
        // g = pi - arcsin(c) + 2k*pi
        auto pi_minus_asin = SymbolicExpr::add(pi, SymbolicExpr::multiply(asin_c, SymbolicExpr::number(-1)));
        auto sol2 = SymbolicExpr::add(pi_minus_asin, two_k_pi)->simplify();
        results.push_back(sol2);
    }
    
    return results;
}

// For cos(g) = c: g = +/- arccos(c) + 2k*pi
static std::vector<std::shared_ptr<SymbolicExpr>> invert_cos(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {
    
    // Domain check: |c| <= 1
    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (std::abs(c_val) > 1.0 + LMMC_REAL_EPSILON) {
            return {};  // Domain violation
        }
    }
    
    auto pi = make_pi();
    auto acos_c = make_arccos(c);
    
    int num_periods = 1;
    if (max_roots > 0) {
        num_periods = (max_roots + 1) / 2;
    }
    
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {
        auto two_k_pi = SymbolicExpr::multiply(SymbolicExpr::number(2 * k), pi);
        
        // g = arccos(c) + 2k*pi
        auto sol1 = SymbolicExpr::add(acos_c, two_k_pi)->simplify();
        results.push_back(sol1);
        
        // g = -arccos(c) + 2k*pi
        auto neg_acos = SymbolicExpr::multiply(acos_c, SymbolicExpr::number(-1));
        auto sol2 = SymbolicExpr::add(neg_acos, two_k_pi)->simplify();
        results.push_back(sol2);
    }
    
    return results;
}

// For tan(g) = c: g = arctan(c) + k*pi
static std::vector<std::shared_ptr<SymbolicExpr>> invert_tan(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {
    
    // No domain restriction for tan inverse
    auto pi = make_pi();
    auto atan_c = make_arctan(c);
    
    int num_periods = 1;
    if (max_roots > 0) {
        num_periods = max_roots;
    }
    
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {
        // g = arctan(c) + k*pi
        auto k_pi = SymbolicExpr::multiply(SymbolicExpr::number(k), pi);
        auto sol = SymbolicExpr::add(atan_c, k_pi)->simplify();
        results.push_back(sol);
    }
    
    return results;
}

// For exp(g) = c: g = ln(c), requires c > 0
static std::vector<std::shared_ptr<SymbolicExpr>> invert_exp(
    const std::shared_ptr<SymbolicExpr>& c) {
    
    // Domain check: c > 0
    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (c_val <= 0.0) {
            return {};  // Domain violation
        }
    }
    
    // g = ln(c)
    auto sol = SymbolicExpr::ln(c)->simplify();
    return {sol};
}

// For ln(g) = c: g = exp(c), with domain constraint g(x) > 0
static std::vector<std::shared_ptr<SymbolicExpr>> invert_ln(
    const std::shared_ptr<SymbolicExpr>& c) {
    
    // g = exp(c) - always valid (exp(c) > 0 for all c, satisfying domain)
    auto sol = SymbolicExpr::exp(c)->simplify();
    return {sol};
}

// For g*exp(g) = c: g = W(c) (Lambert W function)
static std::vector<std::shared_ptr<SymbolicExpr>> invert_lambert_w(
    const std::shared_ptr<SymbolicExpr>& c) {
    
    // g = W(c)
    auto sol = SymbolicExpr::lambertw(c);
    if (!sol) return {};
    return {sol->simplify()};
}

// For a^g = c: g = ln(c)/ln(a), requires a>0, a!=1, c>0
static std::vector<std::shared_ptr<SymbolicExpr>> invert_exp_base(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::shared_ptr<SymbolicExpr>& c) {
    
    // Domain check: c > 0
    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (c_val <= 0.0) {
            return {};  // Domain violation
        }
    }
    
    // g = ln(c) / ln(a)
    auto ln_c = SymbolicExpr::ln(c);
    auto ln_a = SymbolicExpr::ln(base);
    auto sol = SymbolicExpr::divide(ln_c, ln_a)->simplify();
    return {sol};
}

// ============================================================================
// Substitution inversion helpers: given u = h(x) and a u-root, invert h to
// recover x values, discarding solutions outside h⁻¹'s domain.
// ============================================================================

// Invert u = exp(x): x = ln(u_root), requires u_root > 0
static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_exp(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (u_val <= 0.0) return {};  // Domain: exp(x) > 0, so u must be > 0
    }
    auto x_val = SymbolicExpr::ln(u_root)->simplify();
    return {x_val};
}

// Invert u = sin(x): x = arcsin(u_root) (principal value), requires |u_root| <= 1
static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_sin(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (std::abs(u_val) > 1.0 + LMMC_REAL_EPSILON) return {};  // Domain violation
    }
    // Principal value: x = arcsin(u_root)
    auto x_val = make_arcsin(u_root)->simplify();
    return {x_val};
}

// Invert u = cos(x): x = arccos(u_root) (principal value), requires |u_root| <= 1
static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_cos(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (std::abs(u_val) > 1.0 + LMMC_REAL_EPSILON) return {};  // Domain violation
    }
    auto x_val = make_arccos(u_root)->simplify();
    return {x_val};
}

// Invert u = tan(x): x = arctan(u_root) (principal value), no domain restriction
static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_tan(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {
    auto x_val = make_arctan(u_root)->simplify();
    return {x_val};
}

// Invert u = x^n: x = u_root^(1/n), with domain filtering
// For even n: u_root must be >= 0, returns +/- root
// For odd n: returns the real n-th root
static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_power(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var,
    int n) {
    if (n <= 0) return {};
    
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    auto exponent = SymbolicExpr::divide(SymbolicExpr::number(1), SymbolicExpr::number(n));
    
    if (n % 2 == 0) {
        // Even power: u_root must be >= 0
        lmmc_real_t u_val;
        if (try_evaluate_numeric(u_root, u_val)) {
            if (u_val < 0.0) return {};  // No real solution for even root of negative
        }
        // x = +/- u_root^(1/n)
        auto pos_root = SymbolicExpr::power(u_root, exponent)->simplify();
        auto neg_root = SymbolicExpr::multiply(pos_root, SymbolicExpr::number(-1))->simplify();
        results.push_back(pos_root);
        results.push_back(neg_root);
    } else {
        // Odd power: always has a real root
        auto root_val = SymbolicExpr::power(u_root, exponent)->simplify();
        results.push_back(root_val);
    }
    
    return results;
}

// Determine what type of function h is and invert it for a given u_root.
// h is the substitution expression (e.g., exp(x), sin(x), cos(x), tan(x), or x^n).
// Returns x values such that h(x) = u_root, filtering by domain.
static std::vector<std::shared_ptr<SymbolicExpr>> invert_substitution_h(
    const std::shared_ptr<SymbolicExpr>& h_expr,
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {
    
    if (!h_expr || !h_expr->root) return {};
    
    // Case: h = func(var) where func is exp/sin/cos/tan
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(h_expr->root)) {
        if (func->arguments.size() == 1) {
            // Check that the argument is just the variable
            if (auto v = std::dynamic_pointer_cast<VariableNode>(func->arguments[0])) {
                if (v->name == var) {
                    switch (func->type) {
                        case FunctionNode::FuncType::Exp:
                            return invert_h_exp(u_root, var);
                        case FunctionNode::FuncType::Sin:
                            return invert_h_sin(u_root, var);
                        case FunctionNode::FuncType::Cos:
                            return invert_h_cos(u_root, var);
                        case FunctionNode::FuncType::Tan:
                            return invert_h_tan(u_root, var);
                        default:
                            break;
                    }
                }
            }
            // If the argument is more complex (e.g., exp(2x)), we need to solve
            // func(g(x)) = u_root, which means g(x) = func_inverse(u_root)
            // This is handled by the recursive solver below
        }
    }
    
    // Case: h = var (the variable itself, for x^n substitution where the candidate was x)
    if (auto v = std::dynamic_pointer_cast<VariableNode>(h_expr->root)) {
        if (v->name == var) {
            // u = x, so x = u_root directly
            return {u_root};
        }
    }
    
    // Case: h = x^n (PowerNode with base=var, exponent=constant integer)
    // Note: detect_substitution may return the base (x) as the candidate for x^n patterns,
    // but it could also be a more complex expression. Check for PowerNode form.
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(h_expr->root)) {
        if (auto v = std::dynamic_pointer_cast<VariableNode>(pow->base)) {
            if (v->name == var && !depends_on_var(pow->exponent, var)) {
                if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    int n = 0;
                    if (std::holds_alternative<BigInt>(exp_num->value))
                        n = std::get<BigInt>(exp_num->value).to_int();
                    else if (std::holds_alternative<double>(exp_num->value))
                        n = (int)std::get<double>(exp_num->value);
                    else if (std::holds_alternative<Rational>(exp_num->value))
                        n = (int)std::get<Rational>(exp_num->value).to_double();
                    if (n >= 2) {
                        return invert_h_power(u_root, var, n);
                    }
                }
            }
        }
    }
    
    // Fallback: solve h(x) = u_root as an equation
    // Build h(x) - u_root = 0 and try to solve
    auto eq = SymbolicExpr::add(h_expr,
        SymbolicExpr::multiply(u_root, SymbolicExpr::number(-1)))->simplify();
    
    // Try polynomial solve
    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(eq, var);
    if (!poly.is_zero() && poly.degree() >= 1) {
        return SymbolicExpr::solve(eq, var);
    }
    
    return {};
}

// ============================================================================
// Internal recursive implementation with depth tracking
// ============================================================================

static constexpr int MAX_TRANSCENDENTAL_DEPTH = 5;

// Forward declaration of the internal recursive helper
static std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental_impl(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    int depth);

// Helper: recursively solve an inner equation g(x) = value for x.
// Builds the equation g(x) - value = 0 and attempts polynomial solve first,
// then falls back to recursive transcendental solve.
static std::vector<std::shared_ptr<SymbolicExpr>> solve_inner_equation(
    const std::shared_ptr<SymbolicExpr>& inner,
    const std::shared_ptr<SymbolicExpr>& value,
    const std::string& var,
    int depth) {
    
    // If inner is just the variable, the solution is value directly
    if (auto v = std::dynamic_pointer_cast<VariableNode>(inner->root)) {
        if (v->name == var) {
            return {value};
        }
    }
    
    // Build equation: inner - value = 0
    auto inner_eq = SymbolicExpr::add(inner,
        SymbolicExpr::multiply(value, SymbolicExpr::number(-1)))->simplify();
    
    // Try to solve as a polynomial first
    auto inner_poly = symbolic_to_poly<SymbolicPolyCoeff>(inner_eq, var);
    if (!inner_poly.is_zero() && inner_poly.degree() >= 1) {
        auto inner_solutions = SymbolicExpr::solve(inner_eq, var);
        if (!inner_solutions.empty()) return inner_solutions;
    }
    
    // Fall back to recursive transcendental solve
    return solve_transcendental_impl(inner_eq, var, depth + 1);
}

// The main recursive implementation
static std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental_impl(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    int depth) {
    
    if (!expr || !expr->root) return {};
    
    // Depth limit: abort this branch if recursion exceeds MAX_TRANSCENDENTAL_DEPTH
    if (depth > MAX_TRANSCENDENTAL_DEPTH) return {};
    
    // Use default options for max_roots
    SolveOptions opts;
    int max_roots = opts.max_roots; // -1 means default (principal period k=0)
    
    // Try Lambert W pattern first (more specific: g*exp(g) = c)
    auto lambert_pat = decompose_lambert_w_pattern(expr, var);
    if (lambert_pat) {
        auto inverted = invert_lambert_w(lambert_pat->rhs);
        if (inverted.empty()) return {};
        
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(lambert_pat->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }
    
    // Try a^g(x) = c pattern
    auto exp_base_pat = decompose_exp_base_pattern(expr, var);
    if (exp_base_pat) {
        auto inverted = invert_exp_base(exp_base_pat->base, exp_base_pat->rhs);
        if (inverted.empty()) return {};
        
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(exp_base_pat->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }
    
    // Try standard trig/exp/ln patterns
    auto pattern = decompose_trig_exp_pattern(expr, var);
    if (pattern) {
        std::vector<std::shared_ptr<SymbolicExpr>> inverted;
        
        switch (pattern->func_type) {
            case FunctionNode::FuncType::Sin:
                inverted = invert_sin(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Cos:
                inverted = invert_cos(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Tan:
                inverted = invert_tan(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Exp:
                inverted = invert_exp(pattern->rhs);
                break;
            case FunctionNode::FuncType::Ln:
                inverted = invert_ln(pattern->rhs);
                break;
            default:
                return {};
        }
        
        if (inverted.empty()) return {};
        
        // For each inverted value, recursively solve g(x) = value
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(pattern->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }
    
    // Try substitution method: detect u = h(x) that makes the equation polynomial in u
    auto subst = detect_substitution(expr, var);
    if (subst) {
        // Solve the polynomial in u via the polynomial pipeline
        auto u_solutions = SymbolicExpr::solve(subst->poly_in_u, subst->u_var);
        if (u_solutions.empty()) return {};
        
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& u_root : u_solutions) {
            // Invert h to recover x from each u-root, with domain filtering
            auto x_values = invert_substitution_h(subst->u_expr, u_root, var);
            results.insert(results.end(), x_values.begin(), x_values.end());
        }
        return results;
    }
    
    return {};
}

// ============================================================================
// Main entry point: solve_transcendental (public API)
// ============================================================================

std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    
    return solve_transcendental_impl(expr, var, 0);
}

// ============================================================================
// Substitution detection
// ============================================================================

// Visitor to collect transcendental subexpressions that depend on var.
// Collects: exp(g(x)), sin(g(x)), cos(g(x)), tan(g(x)), and x^n (power with
// base depending on var and constant positive integer exponent >= 2).
static void collect_transcendental_subexprs(
    const std::shared_ptr<SymbolicNode>& node,
    const std::string& var,
    std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {
    
    if (!node) return;
    
    // Check FunctionNode: exp, sin, cos, tan
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (func->arguments.size() == 1 && depends_on_var(func->arguments[0], var)) {
            auto ft = func->type;
            if (ft == FunctionNode::FuncType::Exp ||
                ft == FunctionNode::FuncType::Sin ||
                ft == FunctionNode::FuncType::Cos ||
                ft == FunctionNode::FuncType::Tan) {
                candidates.push_back(std::make_shared<SymbolicExpr>(node));
            }
        }
        // Recurse into arguments
        for (auto& arg : func->arguments) {
            collect_transcendental_subexprs(arg, var, candidates);
        }
        return;
    }
    
    // Check PowerNode: x^n where base depends on var and exponent is a constant integer >= 2
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (depends_on_var(pow->base, var) && !depends_on_var(pow->exponent, var)) {
            if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                int e_val = 0;
                if (std::holds_alternative<BigInt>(exp_num->value))
                    e_val = std::get<BigInt>(exp_num->value).to_int();
                else if (std::holds_alternative<double>(exp_num->value))
                    e_val = (int)std::get<double>(exp_num->value);
                else if (std::holds_alternative<Rational>(exp_num->value))
                    e_val = (int)std::get<Rational>(exp_num->value).to_double();
                
                // Only consider as a candidate if the base itself is transcendental or variable
                // and exponent >= 2 (x^n pattern)
                if (e_val >= 2) {
                    // The candidate is the base (e.g., for x^2 + x + 1, candidate is x)
                    auto base_expr = std::make_shared<SymbolicExpr>(pow->base);
                    candidates.push_back(base_expr);
                }
            }
        }
        // Recurse into base and exponent
        collect_transcendental_subexprs(pow->base, var, candidates);
        collect_transcendental_subexprs(pow->exponent, var, candidates);
        return;
    }
    
    // Recurse into AddNode
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }
    
    // Recurse into MultiplyNode
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }
}

// Remove duplicate candidates (by string representation)
static void deduplicate_candidates(std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {
    std::vector<std::shared_ptr<SymbolicExpr>> unique;
    std::vector<std::string> seen;
    for (auto& c : candidates) {
        std::string s = c->to_string();
        bool found = false;
        for (auto& existing : seen) {
            if (existing == s) { found = true; break; }
        }
        if (!found) {
            seen.push_back(s);
            unique.push_back(c);
        }
    }
    candidates = std::move(unique);
}

// Try to extract the multiplier k from exp(k*var) or exp(var) expressions.
// Returns 0 if the node is not of the form exp(k*var).
static int extract_exp_multiplier(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    auto func = std::dynamic_pointer_cast<FunctionNode>(node);
    if (!func || func->type != FunctionNode::FuncType::Exp || func->arguments.size() != 1)
        return 0;
    
    auto arg = func->arguments[0];
    
    // Case: exp(var) → k = 1
    if (auto v = std::dynamic_pointer_cast<VariableNode>(arg)) {
        if (v->name == var) return 1;
    }
    
    // Case: exp(k * var) where k is a positive integer
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(arg)) {
        // Look for a number and the variable in the multiply operands
        std::shared_ptr<NumberNode> num_part = nullptr;
        bool has_var = false;
        bool has_other = false;
        
        for (auto& op : mul->operands) {
            if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) {
                if (!num_part) num_part = n;
                else has_other = true;
            } else if (auto v = std::dynamic_pointer_cast<VariableNode>(op)) {
                if (v->name == var) has_var = true;
                else has_other = true;
            } else {
                has_other = true;
            }
        }
        
        if (has_var && num_part && !has_other) {
            int k = 0;
            if (std::holds_alternative<BigInt>(num_part->value))
                k = std::get<BigInt>(num_part->value).to_int();
            else if (std::holds_alternative<double>(num_part->value))
                k = (int)std::get<double>(num_part->value);
            else if (std::holds_alternative<Rational>(num_part->value))
                k = (int)std::get<Rational>(num_part->value).to_double();
            if (k > 0) return k;
        }
    }
    
    return 0;
}

// Rewrite expression by replacing exp(k*var) with u^k (where u = exp(var)).
// This handles the special case: exp(2x) → u^2, exp(3x) → u^3, etc.
static std::shared_ptr<SymbolicNode> rewrite_exp_as_u_power(
    const std::shared_ptr<SymbolicNode>& node,
    const std::string& var,
    const std::string& u_var) {
    
    if (!node) return node;
    
    // Check if this node is exp(k*var)
    int k = extract_exp_multiplier(node, var);
    if (k > 0) {
        auto u_node = std::make_shared<VariableNode>(u_var);
        if (k == 1) return u_node;
        return SymbolicFactory::create_power(u_node, std::make_shared<NumberNode>(BigInt(k)));
    }
    
    // Recurse into AddNode
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(add->operands.size());
        for (auto& op : add->operands) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_add(std::move(new_ops));
    }
    
    // Recurse into MultiplyNode
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands.size());
        for (auto& op : mul->operands) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_multiply(std::move(new_ops));
    }
    
    // Recurse into PowerNode
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        auto new_base = rewrite_exp_as_u_power(pow->base, var, u_var);
        auto new_exp = rewrite_exp_as_u_power(pow->exponent, var, u_var);
        return SymbolicFactory::create_power(new_base, new_exp);
    }
    
    // Recurse into FunctionNode (non-exp functions)
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        new_args.reserve(func->arguments.size());
        for (auto& arg : func->arguments) {
            new_args.push_back(rewrite_exp_as_u_power(arg, var, u_var));
        }
        return std::make_shared<FunctionNode>(func->type, std::move(new_args));
    }
    
    // Leaf nodes (Number, Variable) - return as-is
    return node;
}

// Check if the expression contains any exp(k*var) terms (for the special case handling)
static bool has_exp_k_var_terms(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    if (!node) return false;
    
    if (extract_exp_multiplier(node, var) > 0) return true;
    
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands)
            if (has_exp_k_var_terms(op, var)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands)
            if (has_exp_k_var_terms(op, var)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (has_exp_k_var_terms(pow->base, var)) return true;
        if (has_exp_k_var_terms(pow->exponent, var)) return true;
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        for (auto& arg : func->arguments)
            if (has_exp_k_var_terms(arg, var)) return true;
    }
    
    return false;
}

std::optional<SubstitutionResult> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    
    if (!expr || !expr->root) return std::nullopt;
    
    // Don't attempt substitution if expression doesn't depend on var
    if (!depends_on_var(expr->root, var)) return std::nullopt;
    
    const std::string u_var = "_u_subst";
    
    // Step 1: Collect transcendental subexpressions
    std::vector<std::shared_ptr<SymbolicExpr>> candidates;
    collect_transcendental_subexprs(expr->root, var, candidates);
    deduplicate_candidates(candidates);
    
    // Step 2: For each candidate h(x), try substituting u for h(x)
    for (auto& h : candidates) {
        // Substitute h(x) → u in the expression
        // We use the string-based substitute: replace the variable pattern
        // Since h(x) might be complex (e.g., exp(x)), we need to substitute
        // by replacing the subexpression structurally.
        
        // Use SymbolicExpr::substitute which replaces a variable name with a value.
        // For transcendental functions, we need a different approach:
        // Replace h(x) with u by substituting in the expression.
        
        // Strategy: Convert h to string, create a fresh variable u,
        // and try to express the equation as polynomial in u.
        // We do this by substituting var_name in h with a known pattern.
        
        // Actually, the simplest approach: if h is a function like exp(x), sin(x), etc.,
        // we can try to substitute the entire subexpression.
        // Since SymbolicExpr::substitute only replaces variables by name,
        // we need to use a structural replacement approach.
        
        // For simple cases where h(x) = func(var) (e.g., exp(x), sin(x)):
        // We can substitute var → ln(u) for exp case, var → arcsin(u) for sin case, etc.
        // But that's complex. Instead, let's use a direct approach:
        // Replace the function node with u_var in the AST.
        
        // Simpler approach: if h = exp(x), sin(x), cos(x), tan(x), or just x (for x^n),
        // try to see if the expression can be written as polynomial in h.
        
        // For h = variable x itself (x^n pattern), use symbolic_to_poly directly
        if (auto v = std::dynamic_pointer_cast<VariableNode>(h->root)) {
            if (v->name == var) {
                // The expression might already be a polynomial in x
                auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, var);
                if (!poly.is_zero() && poly.degree() >= 2) {
                    // This is already a polynomial in x - not really a "substitution"
                    // Skip this candidate as it's trivial
                    continue;
                }
            }
        }
        
        // For transcendental h(x) = func(x):
        // Try substituting by replacing h(x) with a fresh variable u in the expression.
        // We do this by creating a modified expression where h(x) is replaced by u.
        auto u_expr = SymbolicExpr::variable(u_var);
        
        // Use string-based matching: get h's string representation and try to
        // match/replace in the expression. This is fragile but works for simple cases.
        // Better approach: use the substitute method with the function's argument.
        
        // For h = func(var) where func is exp/sin/cos/tan:
        // If the expression only contains h(var) and powers of h(var), we can
        // substitute by replacing the function call with u.
        
        // Direct structural substitution: walk the AST and replace nodes equal to h
        // with the u variable.
        struct SubstNodeVisitor {
            std::shared_ptr<SymbolicNode> target;
            std::string u_name;
            
            std::shared_ptr<SymbolicNode> replace(const std::shared_ptr<SymbolicNode>& node) {
                if (!node) return node;
                
                // Check if this node equals the target
                if (node->equals(*target)) {
                    return std::make_shared<VariableNode>(u_name);
                }
                
                // Recurse
                if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
                    std::vector<std::shared_ptr<SymbolicNode>> new_ops;
                    new_ops.reserve(add->operands.size());
                    for (auto& op : add->operands) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_add(std::move(new_ops));
                }
                if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
                    std::vector<std::shared_ptr<SymbolicNode>> new_ops;
                    new_ops.reserve(mul->operands.size());
                    for (auto& op : mul->operands) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_multiply(std::move(new_ops));
                }
                if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
                    auto new_base = replace(pow->base);
                    auto new_exp = replace(pow->exponent);
                    return SymbolicFactory::create_power(new_base, new_exp);
                }
                if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
                    std::vector<std::shared_ptr<SymbolicNode>> new_args;
                    new_args.reserve(func->arguments.size());
                    for (auto& arg : func->arguments) {
                        new_args.push_back(replace(arg));
                    }
                    return std::make_shared<FunctionNode>(func->type, std::move(new_args));
                }
                
                return node;
            }
        };
        
        SubstNodeVisitor visitor;
        visitor.target = h->root;
        visitor.u_name = u_var;
        
        auto substituted_node = visitor.replace(expr->root);
        auto substituted_expr = std::make_shared<SymbolicExpr>(substituted_node);
        substituted_expr = substituted_expr->simplify();
        
        // Check if the substituted expression no longer depends on var
        if (!depends_on_var(substituted_expr->root, var)) {
            // Try to convert to polynomial in u
            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(substituted_expr, u_var);
            if (!poly.is_zero() && poly.degree() >= 2) {
                return SubstitutionResult{h, substituted_expr, u_var};
            }
        }
    }
    
    // Step 3: Special case - exp(k*x) pattern
    // Detect expressions like exp(2x) + 3*exp(x) + 2 = 0
    // Unify all exp(k*x) as u^k where u = exp(x)
    if (has_exp_k_var_terms(expr->root, var)) {
        auto rewritten_node = rewrite_exp_as_u_power(expr->root, var, u_var);
        auto rewritten_expr = std::make_shared<SymbolicExpr>(rewritten_node);
        rewritten_expr = rewritten_expr->simplify();
        
        // Check if the rewritten expression no longer depends on var
        if (!depends_on_var(rewritten_expr->root, var)) {
            // Try to convert to polynomial in u
            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(rewritten_expr, u_var);
            if (!poly.is_zero() && poly.degree() >= 2) {
                // u = exp(x)
                auto exp_x = SymbolicExpr::exp(SymbolicExpr::variable(var));
                return SubstitutionResult{exp_x, rewritten_expr, u_var};
            }
        }
    }
    
    return std::nullopt;
}

} // namespace lamina
