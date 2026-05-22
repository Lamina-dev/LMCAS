#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <variant>

namespace lamina {

// ============================================================================
// Utility helpers (file-local)
// ============================================================================

static std::shared_ptr<SymbolicExpr> make_expr_ptr(const SymbolicExpr& e) {
    return std::make_shared<SymbolicExpr>(e);
}

static bool valid_dependency(const SymbolicExpr& expr, const std::string& var) {
    auto diff = expr.differentiate(var);
    if (!diff) return false;
    auto simp_diff = diff->simplify();
    return !simp_diff->is_zero();
}

static std::shared_ptr<SymbolicExpr> sym_sub(const SymbolicExpr& a, const SymbolicExpr& b) {
    auto neg_b = SymbolicExpr::multiply(SymbolicExpr::number(-1), std::make_shared<SymbolicExpr>(b));
    return SymbolicExpr::add(std::make_shared<SymbolicExpr>(a), neg_b);
}

static std::shared_ptr<SymbolicExpr> sym_rational(long long num, long long den) {
    return SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
}

static std::shared_ptr<SymbolicExpr> make_arctan(const std::shared_ptr<SymbolicExpr>& op) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            FunctionNode::FuncType::ArcTan,
            std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
}

static bool has_integral_node_check(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Calculus_Integral) return true;
        for (auto& arg : fn->arguments)
            if (has_integral_node_check(arg)) return true;
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

// ============================================================================
// IntegrationTable
// ============================================================================

const std::vector<IntegrationEntry> IntegrationTable::empty_entries_;

IntegrationTable::IntegrationTable() {
    load_defaults();
}

void IntegrationTable::add_entry(Category cat, const IntegrationEntry& entry) {
    entries_[static_cast<int>(cat)].push_back(entry);
    // Keep sorted by priority
    auto& vec = entries_[static_cast<int>(cat)];
    std::sort(vec.begin(), vec.end(), [](const IntegrationEntry& a, const IntegrationEntry& b) {
        return a.priority < b.priority;
    });
}

void IntegrationTable::clear_category(Category cat) {
    entries_[static_cast<int>(cat)].clear();
}

const std::vector<IntegrationEntry>& IntegrationTable::get_entries(Category cat) const {
    auto it = entries_.find(static_cast<int>(cat));
    if (it == entries_.end()) return empty_entries_;
    return it->second;
}

std::vector<const IntegrationEntry*> IntegrationTable::get_all_sorted() const {
    std::vector<const IntegrationEntry*> all;
    for (const auto& [cat, vec] : entries_) {
        for (const auto& entry : vec) {
            all.push_back(&entry);
        }
    }
    std::sort(all.begin(), all.end(), [](const IntegrationEntry* a, const IntegrationEntry* b) {
        return a->priority < b->priority;
    });
    return all;
}

void IntegrationTable::load_defaults() {
    // Helper lambda: check that a wildcard binding does NOT depend on the integration variable
    auto independent_of_var = [](const std::string& wc_name) {
        return [wc_name](const MatchMap& m, const std::string& var) -> bool {
            auto it = m.find(wc_name);
            if (it == m.end()) return true;
            return !depends_on_var(it->second.root, var);
        };
    };

    // --- Exponential rules ---
    {
        // exp(a*x) -> (1/a)*exp(a*x), where a is independent of x
        auto a = wildcard("_a");
        auto u = wildcard("_u"); // represents x
        auto pat = *SymbolicExpr::exp(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1)),
            SymbolicExpr::exp(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u))));
        add_entry(Category::Exponential, IntegrationEntry(
            "exp(a*x)", pat, res, {"_a", "_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it_a = m.find("_a");
                auto it_u = m.find("_u");
                if (it_a == m.end() || it_u == m.end()) return false;
                // _a must be independent of var, _u must BE the var
                if (depends_on_var(it_a->second.root, var)) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it_u->second.root);
                return v && v->name == var;
            }, 50));
    }
    {
        // exp(x) -> exp(x)
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::exp(make_expr_ptr(u));
        auto res = *SymbolicExpr::exp(make_expr_ptr(u));
        add_entry(Category::Exponential, IntegrationEntry(
            "exp(x)", pat, res, {"_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it = m.find("_u");
                if (it == m.end()) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
                return v && v->name == var;
            }, 60));
    }

    // --- Trigonometric rules ---
    {
        // sin(x) -> -cos(x)
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::sin(make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)", pat, res, {"_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it = m.find("_u");
                if (it == m.end()) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
                return v && v->name == var;
            }, 60));
    }
    {
        // cos(x) -> sin(x)
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::cos(make_expr_ptr(u));
        auto res = *SymbolicExpr::sin(make_expr_ptr(u));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)", pat, res, {"_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it = m.find("_u");
                if (it == m.end()) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
                return v && v->name == var;
            }, 60));
    }
    {
        // tan(x) -> -ln(cos(x))
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::tan(make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::ln(SymbolicExpr::cos(make_expr_ptr(u))));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "tan(x)", pat, res, {"_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it = m.find("_u");
                if (it == m.end()) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
                return v && v->name == var;
            }, 60));
    }

    // --- Logarithmic rules ---
    {
        // ln(x) -> x*ln(x) - x
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::ln(make_expr_ptr(u));
        // x*ln(x) - x
        auto x_ln_x = SymbolicExpr::multiply(make_expr_ptr(u), SymbolicExpr::ln(make_expr_ptr(u)));
        auto res = *sym_sub(*x_ln_x, u);
        add_entry(Category::Logarithmic, IntegrationEntry(
            "ln(x)", pat, res, {"_u"},
            [](const MatchMap& m, const std::string& var) {
                auto it = m.find("_u");
                if (it == m.end()) return false;
                auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
                return v && v->name == var;
            }, 60));
    }
}

// ============================================================================
// TableLookupStrategy
// ============================================================================

std::shared_ptr<SymbolicExpr> TableLookupStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx) {
    
    auto all_entries = ctx.table().get_all_sorted();
    
    for (const auto* entry : all_entries) {
        MatchMap bindings;
        if (Matcher::match(entry->pattern, expr, entry->wildcards, bindings)) {
            // Check condition
            if (entry->condition && !entry->condition(bindings, var)) {
                continue;
            }
            // Apply replacement
            SymbolicExpr result = Matcher::replace(entry->result, bindings, false);
            auto simplified = result.simplify();
            if (simplified && !has_integral_node_check(simplified->root)) {
                return simplified;
            }
        }
    }
    return nullptr;
}

// ============================================================================
// PowerRuleStrategy
// ============================================================================

std::shared_ptr<SymbolicExpr> PowerRuleStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx) {
    
    // Case 1: x -> x^2/2
    if (auto v_node = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        if (v_node->name == var) {
            return SymbolicExpr::multiply(
                SymbolicExpr::power(make_expr_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2));
        }
    }

    // Case 2: x^n where n is constant w.r.t. var
    if (auto p_node = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        SymbolicExpr base(p_node->base);
        SymbolicExpr exp_expr(p_node->exponent);

        if (auto b_var = std::dynamic_pointer_cast<VariableNode>(base.root)) {
            if (b_var->name == var && !valid_dependency(exp_expr, var)) {
                auto n_plus_1 = SymbolicExpr::add(make_expr_ptr(exp_expr), SymbolicExpr::number(1))->simplify();
                if (n_plus_1->is_zero()) {
                    // n = -1 -> ln(x)
                    return SymbolicExpr::ln(make_expr_ptr(base));
                }
                return SymbolicExpr::divide(
                    SymbolicExpr::power(make_expr_ptr(base), n_plus_1),
                    n_plus_1);
            }
        }
    }
    return nullptr;
}

// ============================================================================
// SubstitutionStrategy
// ============================================================================

std::shared_ptr<SymbolicExpr> SubstitutionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx) {
    
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        ops.push_back(expr.root);
    }
    
    for (size_t i = 0; i < ops.size(); ++i) {
        SymbolicExpr candidate_term(ops[i]);
        SymbolicExpr u;
        bool possible = false;

        if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
            u = SymbolicExpr(pow->base);
            possible = true;
        } else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
            if (!func->arguments.empty()) {
                u = SymbolicExpr(func->arguments[0]);
                possible = true;
            }
        }

        if (possible && valid_dependency(u, var)) {
            auto d_ptr = u.differentiate(var);
            if (!d_ptr) continue;
            auto du = d_ptr->simplify();
            if (du->is_zero()) continue;

            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(make_expr_ptr(f_u), make_expr_ptr(*du));
            auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), term_times_du)->simplify();
            
            if (!valid_dependency(*ratio, var)) {
                std::shared_ptr<SymbolicExpr> prim = nullptr;
                
                if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
                    SymbolicExpr n(pow->exponent);
                    auto np1 = SymbolicExpr::add(make_expr_ptr(n), SymbolicExpr::number(1))->simplify();
                    if (np1->is_zero()) {
                        prim = SymbolicExpr::ln(make_expr_ptr(u));
                    } else {
                        prim = SymbolicExpr::divide(
                            SymbolicExpr::power(make_expr_ptr(u), np1), np1);
                    }
                } else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
                    if (func->type == FunctionNode::FuncType::Cos) {
                        prim = SymbolicExpr::sin(make_expr_ptr(u));
                    } else if (func->type == FunctionNode::FuncType::Sin) {
                        prim = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
                    } else if (func->type == FunctionNode::FuncType::Exp) {
                        prim = SymbolicExpr::exp(make_expr_ptr(u));
                    }
                }
                
                if (prim) {
                    return SymbolicExpr::multiply(ratio, prim);
                }
            }
        }
    }
    return nullptr;
}

// ============================================================================
// PartialFractionStrategy
// ============================================================================

std::shared_ptr<SymbolicExpr> PartialFractionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx) {
    
    std::shared_ptr<SymbolicExpr> den = nullptr;
    
    if (auto p = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        double exp_val = 0;
        bool is_inv = false;
        if (auto num_node = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
            if (std::holds_alternative<lmmc_real_t>(num_node->value))
                exp_val = std::get<lmmc_real_t>(num_node->value);
            else if (std::holds_alternative<Rational>(num_node->value))
                exp_val = std::get<Rational>(num_node->value).to_double();
            else if (std::holds_alternative<BigInt>(num_node->value))
                exp_val = std::get<BigInt>(num_node->value).to_double();
            int eq;
            lmmc_double_nearly_equal_tol(exp_val, -1.0, 1e-9, 1e-9, &eq);
            if (eq) is_inv = true;
        }
        if (is_inv) {
            den = make_expr_ptr(SymbolicExpr(p->base));
        }
    }
    
    if (!den) return nullptr;

    try {
        Polynomial<SymbolicPolyCoeff> Q = symbolic_to_poly<SymbolicPolyCoeff>(den, var);
        
        if (Q.degree() == 2) {
            SymbolicExpr c_expr = *(Q.coeffs[0].val);
            SymbolicExpr b_expr = *(Q.coeffs[1].val);
            SymbolicExpr a_expr = *(Q.coeffs[2].val);
            
            // All three coefficients must be numeric, otherwise we can't
            // compute the discriminant correctly. Symbolic coefficients
            // require a more general approach (factoring/CAS) which is
            // beyond this strategy's scope.
            if (!a_expr.is_number() || !b_expr.is_number() || !c_expr.is_number()) {
                return nullptr;
            }
            
            double a = a_expr.to_numeric();
            double b = b_expr.to_numeric();
            double c = c_expr.to_numeric();
            
            int eq_a;
            lmmc_double_nearly_equal_tol(a, 0.0, 1e-9, 1e-9, &eq_a);
            if (eq_a) return nullptr;
            
            double delta = b * b - 4 * a * c;
            int eq_delta;
            lmmc_double_nearly_equal_tol(delta, 0.0, 1e-9, 1e-9, &eq_delta);
            
            if (!eq_delta && delta > 0) {
                double sqrt_delta = std::sqrt(delta);
                auto scalar = SymbolicExpr::number(1.0 / sqrt_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto two_a_x = SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var));
                auto two_a_x_plus_b = SymbolicExpr::add(make_expr_ptr(*two_a_x), make_expr_ptr(*b_num));
                auto term1_arg = sym_sub(*two_a_x_plus_b, *SymbolicExpr::number(sqrt_delta));
                auto term2_arg = SymbolicExpr::add(make_expr_ptr(*two_a_x_plus_b), SymbolicExpr::number(sqrt_delta));
                auto term1 = SymbolicExpr::ln(make_expr_ptr(*term1_arg));
                auto term2 = SymbolicExpr::ln(make_expr_ptr(*term2_arg));
                return SymbolicExpr::multiply(scalar, sym_sub(*term1, *term2));
            } else if (!eq_delta && delta < 0) {
                double sqrt_neg_delta = std::sqrt(-delta);
                auto scalar = SymbolicExpr::number(2.0 / sqrt_neg_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto num = SymbolicExpr::add(
                    SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var)),
                    make_expr_ptr(*b_num));
                auto inner = SymbolicExpr::divide(make_expr_ptr(*num), SymbolicExpr::number(sqrt_neg_delta));
                return SymbolicExpr::multiply(scalar, make_arctan(inner));
            }
        }
    } catch (...) {}
    
    return nullptr;
}

// ============================================================================
// IBPStrategy
// ============================================================================

std::shared_ptr<SymbolicExpr> IBPStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx) {
    
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        ops.push_back(expr.root);
    }
    
    int best_u_idx = -1;
    int best_score = 100;
    
    auto get_score = [&](const std::shared_ptr<SymbolicNode>& node) -> int {
        SymbolicExpr e(node);
        if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
            if (fn->type == FunctionNode::FuncType::Ln || fn->type == FunctionNode::FuncType::Log) return 1;
            if (fn->type == FunctionNode::FuncType::ArcSin || fn->type == FunctionNode::FuncType::ArcTan) return 2;
            if (fn->type == FunctionNode::FuncType::Sin || fn->type == FunctionNode::FuncType::Cos) return 4;
            if (fn->type == FunctionNode::FuncType::Exp) return 5;
        }
        if (!valid_dependency(e, var)) return 10;
        if (std::dynamic_pointer_cast<VariableNode>(node)) return 3;
        if (std::dynamic_pointer_cast<PowerNode>(node)) return 3;
        return 10;
    };
    
    if (ops.size() == 1) {
        int s = get_score(ops[0]);
        if (s <= 2) best_u_idx = 0;
    } else {
        for (size_t i = 0; i < ops.size(); ++i) {
            if (!valid_dependency(SymbolicExpr(ops[i]), var)) continue;
            int s = get_score(ops[i]);
            if (s < best_score) {
                best_score = s;
                best_u_idx = (int)i;
            }
        }
    }
    
    if (best_u_idx == -1) return nullptr;
    
    SymbolicExpr u(ops[best_u_idx]);
    
    std::vector<std::shared_ptr<SymbolicNode>> dv_ops;
    for (size_t i = 0; i < ops.size(); ++i) {
        if ((int)i != best_u_idx) dv_ops.push_back(ops[i]);
    }
    
    std::shared_ptr<SymbolicExpr> dv;
    if (dv_ops.empty()) dv = SymbolicExpr::number(1);
    else if (dv_ops.size() == 1) dv = make_expr_ptr(SymbolicExpr(dv_ops[0]));
    else dv = std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(dv_ops));
    
    auto v = ctx.integrate_recursive(*dv, var, ctx.max_depth() - 2);
    if (has_integral_node_check(v->root)) return nullptr;

    auto du_ptr = u.differentiate(var);
    if (!du_ptr) return nullptr;
    auto du = make_expr_ptr(*du_ptr->simplify());

    auto uv = SymbolicExpr::multiply(make_expr_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    auto int_vdu = ctx.integrate_recursive(*vdu, var, ctx.max_depth() - 1);
    
    return sym_sub(*uv, *int_vdu);
}

// ============================================================================
// Integrator
// ============================================================================

Integrator::Integrator() {
    // Register default strategies in order of priority
    strategies_.push_back(std::make_unique<TableLookupStrategy>());
    strategies_.push_back(std::make_unique<PowerRuleStrategy>());
    strategies_.push_back(std::make_unique<SubstitutionStrategy>());
    strategies_.push_back(std::make_unique<PartialFractionStrategy>());
    strategies_.push_back(std::make_unique<IBPStrategy>());
}

void Integrator::add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position) {
    if (!strategy) {
        throw std::invalid_argument("Integrator::add_strategy: strategy must not be null");
    }
    if (position < 0 || position >= (int)strategies_.size()) {
        strategies_.push_back(std::move(strategy));
    } else {
        strategies_.insert(strategies_.begin() + position, std::move(strategy));
    }
}

bool Integrator::depends_on(const SymbolicExpr& expr, const std::string& var) {
    return depends_on_var(expr.root, var);
}

std::shared_ptr<SymbolicExpr> Integrator::make_integral_node(
    const SymbolicExpr& expr, const std::string& var) {
    std::vector<std::shared_ptr<SymbolicNode>> args;
    args.push_back(expr.root);
    args.push_back(SymbolicExpr::variable(var)->root);
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
}

std::shared_ptr<SymbolicExpr> Integrator::check_cycle(
    const SymbolicExpr& expr, const std::string& var) {
    for (size_t i = 0; i < cycle_state_.history.size(); ++i) {
        auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), make_expr_ptr(cycle_state_.history[i]))->simplify();
        if (!valid_dependency(*ratio, var)) {
            return SymbolicExpr::multiply(ratio, SymbolicExpr::variable("INT_CYCLE_" + std::to_string(i)));
        }
    }
    return nullptr;
}

void Integrator::resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx) {
    std::string cycle_var = "INT_CYCLE_" + std::to_string(cycle_idx);
    if (valid_dependency(*result, cycle_var)) {
        auto B = result->differentiate(cycle_var)->simplify();
        auto A = result->substitute(cycle_var, SymbolicExpr::number(0))->simplify();
        auto one_minus_B = sym_sub(*SymbolicExpr::number(1), *B)->simplify();
        if (!one_minus_B->is_zero()) {
            result = SymbolicExpr::divide(A, one_minus_B);
        }
    }
}

std::shared_ptr<SymbolicExpr> Integrator::apply_linearity(
    const SymbolicExpr& expr, const std::string& var) {
    
    auto simp_expr = expr.simplify();
    
    // Factor out constants from multiplication
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> constants;
        std::vector<std::shared_ptr<SymbolicNode>> dependents;
        
        for (auto& op : mul->operands) {
            SymbolicExpr term(op);
            if (!valid_dependency(term, var)) {
                constants.push_back(op);
            } else {
                dependents.push_back(op);
            }
        }
        
        if (!constants.empty() && dependents.size() < mul->operands.size()) {
            SymbolicExpr const_part = (constants.size() == 1) ?
                SymbolicExpr(constants[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(constants));
            SymbolicExpr dep_part = (dependents.empty()) ?
                *SymbolicExpr::number(1) :
                ((dependents.size() == 1) ? SymbolicExpr(dependents[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(dependents)));
            
            auto int_part = integrate(dep_part, var);
            return SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(const_part), std::make_shared<SymbolicExpr>(int_part));
        }
    }
    
    // Split sums
    if (auto add = std::dynamic_pointer_cast<AddNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> results;
        for (auto& op : add->operands) {
            SymbolicExpr term(op);
            auto int_term = integrate(term, var);
            results.push_back(int_term.root);
        }
        return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(results));
    }
    
    return nullptr; // linearity not applicable
}

std::shared_ptr<SymbolicExpr> Integrator::dispatch_strategies(
    const SymbolicExpr& expr, const std::string& var, int depth) {
    
    for (auto& strategy : strategies_) {
        auto result = strategy->try_integrate(expr, var, *this);
        if (result) {
            err_stream << "[Integration] Strategy '" << strategy->name() << "' succeeded\n";
            return result;
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> Integrator::integrate_recursive(
    const SymbolicExpr& expr, const std::string& var, int depth) {
    
    if (depth > max_depth_) {
        return make_integral_node(expr, var);
    }
    
    // Cycle detection
    auto cycle_result = check_cycle(expr, var);
    if (cycle_result) return cycle_result;
    
    cycle_state_.history.push_back(expr);
    size_t my_idx = cycle_state_.history.size() - 1;
    
    // Constant check
    if (expr.is_number() || !valid_dependency(expr, var)) {
        cycle_state_.history.pop_back();
        return SymbolicExpr::multiply(make_expr_ptr(expr), SymbolicExpr::variable(var));
    }
    
    auto result = dispatch_strategies(expr, var, depth);
    
    cycle_state_.history.pop_back();
    
    if (!result) {
        return make_integral_node(expr, var);
    }
    
    // Resolve any cycle variables
    resolve_cycle(result, my_idx);
    
    return result;
}

SymbolicExpr Integrator::integrate(const SymbolicExpr& expr, const std::string& var_name) {
    // Reset cycle state for this top-level call
    cycle_state_.history.clear();
    
    // Try linearity first (sums and constant factors)
    auto linear_result = apply_linearity(expr, var_name);
    if (linear_result) return *linear_result;
    
    // Delegate to recursive handler
    return *integrate_recursive(expr, var_name, 0);
}

SymbolicExpr Integrator::integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                                        const SymbolicExpr& lower, const SymbolicExpr& upper) {
    // Check for singularities (basic 1/x case)
    SymbolicExpr simp_expr_val = *expr.simplify();
    bool is_inv_x = false;
    
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(simp_expr_val.root)) {
        if (auto v = std::dynamic_pointer_cast<VariableNode>(pow->base)) {
            if (v->name == var_name) {
                if (auto en = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    int eq_minus_one = 0;
                    if (std::holds_alternative<lmmc_real_t>(en->value)) {
                        lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(en->value), -1.0, 1e-9, 1e-9, &eq_minus_one);
                    }
                    if ((std::holds_alternative<lmmc_real_t>(en->value) && eq_minus_one != 0) ||
                        (std::holds_alternative<BigInt>(en->value) && std::get<BigInt>(en->value).to_int() == -1) ||
                        (std::holds_alternative<Rational>(en->value) && std::get<Rational>(en->value).to_double() == -1.0)) {
                        is_inv_x = true;
                    }
                }
            }
        }
    }

    // Check numeric bounds first; only call to_numeric() if bounds are numbers.
    bool numeric_bounds = (lower.root && std::dynamic_pointer_cast<NumberNode>(lower.root)) &&
                          (upper.root && std::dynamic_pointer_cast<NumberNode>(upper.root));

    if (is_inv_x && numeric_bounds) {
        double l_val = lower.to_numeric();
        double u_val = upper.to_numeric();
        int eq_l, eq_u;
        lmmc_double_nearly_equal_tol(l_val, 0.0, 1e-9, 1e-9, &eq_l);
        lmmc_double_nearly_equal_tol(u_val, 0.0, 1e-9, 1e-9, &eq_u);
        if (!eq_l && l_val < 0 && !eq_u && u_val > 0) {
            auto t = std::make_shared<SymbolicExpr>(*SymbolicExpr::variable("t"));
            auto zero = std::make_shared<SymbolicExpr>(*SymbolicExpr::number(0));
            auto int_left = integrate_def(expr, var_name, lower, *t);
            auto lim_left = int_left.limit("t", zero, "-");
            auto int_right = integrate_def(expr, var_name, *t, upper);
            auto lim_right = int_right.limit("t", zero, "+");
            if (lim_left && lim_right) {
                return *SymbolicExpr::add(lim_left, lim_right);
            }
        }
    }

    // Calculate indefinite integral
    SymbolicExpr indefinite = integrate(expr, var_name);
    
    // Check if integration failed
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(indefinite.root)) {
        if (func->type == FunctionNode::FuncType::Calculus_Integral) {
            std::vector<std::shared_ptr<SymbolicNode>> args;
            args.push_back(expr.root);
            args.push_back(SymbolicExpr::variable(var_name)->root);
            args.push_back(lower.root);
            args.push_back(upper.root);
            return SymbolicExpr(std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
        }
    }

    // FTC: F(b) - F(a)
    auto F_b = indefinite.substitute(var_name, make_expr_ptr(upper));
    auto F_a = indefinite.substitute(var_name, make_expr_ptr(lower));
    auto result = sym_sub(*F_b, *F_a);
    return *result->simplify();
}

} // namespace lamina
