// solve_polynomial.cpp - Cubic (Cardano) and Quartic (Ferrari) closed-form solvers

#include "solve_polynomial.hpp"
#include "root_of_utils.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace lamina {

// Helper: check if a SymbolicExpr is purely numeric (no variables at all)
static bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !expr->root) return true;
    
    struct VarDetector : public SymbolicVisitor {
        bool has_var = false;
        void visit(NumberNode&) override {}
        void visit(VariableNode&) override { has_var = true; }
        void visit(AddNode& n) override { for (auto& op : n.operands) { if (has_var) return; op->accept(*this); } }
        void visit(MultiplyNode& n) override { for (auto& op : n.operands) { if (has_var) return; op->accept(*this); } }
        void visit(PowerNode& n) override { n.base->accept(*this); if (!has_var) n.exponent->accept(*this); }
        void visit(FunctionNode& n) override { for (auto& arg : n.arguments) { if (has_var) return; arg->accept(*this); } }
        void visit(MatrixNode&) override {}
    } detector;
    
    expr->root->accept(detector);
    return !detector.has_var;
}

// Helper: create cbrt(x) = power(x, 1/3)
static std::shared_ptr<SymbolicExpr> cbrt_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3)));
}

// Helper: create acos(x) using FunctionNode
static std::shared_ptr<SymbolicExpr> acos_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<SymbolicNode>>{x->root}));
}

// Helper: negate an expression
static std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
}

// Helper: subtract b from a
static std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    return SymbolicExpr::add(a, negate(b));
}

// Helper constants
static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> num_r(int p, int q) { return SymbolicExpr::number(Rational(p, q)); }

// Forward declaration for lower-degree solving
static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

// Solve linear equation ax + b = 0
static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b) {
    // x = -b/a
    auto neg_b = negate(b);
    return { SymbolicExpr::divide(neg_b, a)->simplify() };
}

// Solve quadratic equation ax^2 + bx + c = 0
static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {
    // discriminant = b^2 - 4ac
    auto b2 = SymbolicExpr::power(b, num(2));
    auto four_ac = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, c));
    auto delta = sub(b2, four_ac);
    
    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto neg_b = negate(b);
    auto two_a = SymbolicExpr::multiply(num(2), a);
    
    auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
    auto x2 = SymbolicExpr::divide(sub(neg_b, sqrt_delta), two_a)->simplify();
    
    return { x1, x2 };
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string& var) {
    
    // Check if a == 0: delegate to lower-degree solver
    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {
        // Not a cubic, delegate to quadratic or linear
        auto b_simp = b->simplify();
        if (b_simp->get_number_value_is_zero()) {
            // Linear: cx + d = 0
            return solve_linear_internal(c, d);
        }
        // Quadratic: bx^2 + cx + d = 0
        return solve_quadratic_internal(b, c, d);
    }
    
    // Step 1: Compute depressed cubic coefficients
    // p = (3ac - b^2) / (3a^2)
    // q = (2b^3 - 9abc + 27a^2*d) / (27a^3)
    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    
    // p = (3ac - b^2) / (3a^2)
    auto three_ac = SymbolicExpr::multiply(num(3), SymbolicExpr::multiply(a, c));
    auto p_num = sub(three_ac, b2);
    auto p_den = SymbolicExpr::multiply(num(3), a2);
    auto p = SymbolicExpr::divide(p_num, p_den)->simplify();
    
    // q = (2b^3 - 9abc + 27a^2*d) / (27a^3)
    auto two_b3 = SymbolicExpr::multiply(num(2), b3);
    auto nine_abc = SymbolicExpr::multiply(num(9), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto twentyseven_a2d = SymbolicExpr::multiply(num(27), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(two_b3, nine_abc), twentyseven_a2d);
    auto q_den = SymbolicExpr::multiply(num(27), a3);
    auto q = SymbolicExpr::divide(q_num, q_den)->simplify();
    
    // Shift: x = t - b/(3a)
    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(3), a))->simplify();
    // shift = b/(3a), so x = t - shift means we subtract shift from each t
    
    // Check if coefficients are all numeric (no symbolic parameters)
    bool all_numeric = is_purely_numeric(p) && is_purely_numeric(q);
    
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    
    if (all_numeric) {
        // Numeric path: evaluate discriminant to choose branch
        double p_val = p->to_numeric();
        double q_val = q->to_numeric();
        
        // D = (q/2)^2 + (p/3)^3
        double D = (q_val / 2.0) * (q_val / 2.0) + (p_val / 3.0) * (p_val / 3.0) * (p_val / 3.0);
        
        const double eps = 1e-12;
        
        if (std::abs(p_val) < eps && std::abs(q_val) < eps) {
            // Triple root: t = 0, so x = -b/(3a)
            double shift_val = shift->to_numeric();
            auto root = SymbolicExpr::number(-shift_val);
            roots = { root, root, root };
        } else if (std::abs(D) < eps) {
            // D = 0: repeated roots
            // t1 = -2*cbrt(q/2), t2 = t3 = cbrt(q/2)
            double q_half_val = q_val / 2.0;
            double cbrt_q_half_val = std::cbrt(q_half_val);
            double shift_val = shift->to_numeric();
            
            // t1 = -2 * cbrt(q/2)
            double t1_val = -2.0 * cbrt_q_half_val;
            // t2 = t3 = cbrt(q/2)
            double t2_val = cbrt_q_half_val;
            
            auto x1 = SymbolicExpr::number(t1_val - shift_val);
            auto x2 = SymbolicExpr::number(t2_val - shift_val);
            
            roots = { x1, x2, x2 };
        } else if (D > eps) {
            // D > 0: one real root + two complex conjugate roots
            // Compute numerically: u = cbrt(-q/2 + sqrt(D)), v = cbrt(-q/2 - sqrt(D))
            double neg_q_half_val = -q_val / 2.0;
            double sqrt_D_val = std::sqrt(D);
            
            double u_arg_val = neg_q_half_val + sqrt_D_val;
            double v_arg_val = neg_q_half_val - sqrt_D_val;
            
            double u_val = std::cbrt(u_arg_val);
            double v_val = std::cbrt(v_arg_val);
            
            // Compute shift numerically
            double shift_val = shift->to_numeric();
            
            // Real root: t1 = u + v, x1 = t1 - shift
            double t1_val = u_val + v_val;
            double x1_val = t1_val - shift_val;
            
            // Complex roots: t2,t3 = -(u+v)/2 ± i*sqrt(3)*(u-v)/2
            double real_part = -(u_val + v_val) / 2.0 - shift_val;
            double imag_part_val = std::sqrt(3.0) * (u_val - v_val) / 2.0;
            
            auto x1 = SymbolicExpr::number(x1_val);
            
            // Express complex roots symbolically for proper representation
            auto i_unit = SymbolicExpr::sqrt(num(-1));
            auto x2 = SymbolicExpr::add(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part_val), i_unit))->simplify();
            auto x3 = sub(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part_val), i_unit))->simplify();
            
            roots = { x1, x2, x3 };
        } else {
            // D < 0: Casus irreducibilis - three distinct real roots
            // Use trigonometric form:
            // r = sqrt(-p^3/27)
            // θ = arccos(-q / (2r))
            // t_k = 2*cbrt(r) * cos((θ + 2kπ)/3), k=0,1,2
            
            double neg_p3_over_27 = -(p_val * p_val * p_val) / 27.0;
            double r_val = std::sqrt(neg_p3_over_27);
            double cos_arg = -q_val / (2.0 * r_val);
            // Clamp to [-1, 1] for numerical safety
            if (cos_arg > 1.0) cos_arg = 1.0;
            if (cos_arg < -1.0) cos_arg = -1.0;
            double theta_val = std::acos(cos_arg);
            double two_cbrt_r = 2.0 * std::cbrt(r_val);
            
            // Compute shift numerically
            double shift_val = shift->to_numeric();
            
            for (int k = 0; k < 3; ++k) {
                double angle = (theta_val + 2.0 * k * M_PI) / 3.0;
                double t_k_val = two_cbrt_r * std::cos(angle);
                double x_k_val = t_k_val - shift_val;
                roots.push_back(SymbolicExpr::number(x_k_val));
            }
        }
    } else {
        // Symbolic path: coefficients contain parameters
        // Return general Cardano form with cbrt/sqrt without evaluating discriminant sign
        
        // D = (q/2)^2 + (p/3)^3
        auto q_half = SymbolicExpr::divide(q, num(2));
        auto p_third = SymbolicExpr::divide(p, num(3));
        auto D_expr = SymbolicExpr::add(
            SymbolicExpr::power(q_half, num(2)),
            SymbolicExpr::power(p_third, num(3)))->simplify();
        
        // u = cbrt(-q/2 + sqrt(D))
        auto neg_q_half = negate(q_half)->simplify();
        auto sqrt_D = SymbolicExpr::sqrt(D_expr);
        
        auto u = cbrt_expr(SymbolicExpr::add(neg_q_half, sqrt_D))->simplify();
        auto v = cbrt_expr(sub(neg_q_half, sqrt_D))->simplify();
        
        // Real root: t1 = u + v, x1 = t1 - b/(3a)
        auto t1 = SymbolicExpr::add(u, v);
        auto x1 = sub(t1, shift)->simplify();
        
        // Complex roots using cube roots of unity:
        // t2 = ω*u + ω²*v, t3 = ω²*u + ω*v
        // where ω = (-1 + i*sqrt(3))/2, ω² = (-1 - i*sqrt(3))/2
        // Equivalently:
        // x2 = -(u+v)/2 + i*sqrt(3)*(u-v)/2 - shift
        // x3 = -(u+v)/2 - i*sqrt(3)*(u-v)/2 - shift
        auto neg_half_sum = SymbolicExpr::divide(negate(SymbolicExpr::add(u, v)), num(2));
        auto diff_uv = sub(u, v);
        auto sqrt3_half = SymbolicExpr::divide(SymbolicExpr::sqrt(num(3)), num(2));
        auto imag_coeff = SymbolicExpr::multiply(sqrt3_half, diff_uv);
        
        auto i_unit = SymbolicExpr::sqrt(num(-1));
        auto x2 = sub(SymbolicExpr::add(neg_half_sum, SymbolicExpr::multiply(i_unit, imag_coeff)), shift)->simplify();
        auto x3 = sub(sub(neg_half_sum, SymbolicExpr::multiply(i_unit, imag_coeff)), shift)->simplify();
        
        roots = { x1, x2, x3 };
    }
    
    return roots;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string& var) {
    
    // ax^4 + bx^2 + c = 0, substitute u = x^2
    // Check if all coefficients are numeric
    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) && is_purely_numeric(c);
    
    if (all_numeric) {
        // Numeric path
        double av = a->to_numeric();
        double bv = b->to_numeric();
        double cv = c->to_numeric();
        
        // Solve au^2 + bu + c = 0
        double disc = bv * bv - 4.0 * av * cv;
        double sqrt_disc = std::sqrt(std::abs(disc));
        
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        const double eps = 1e-12;
        
        if (disc >= -eps) {
            double u1 = (-bv + sqrt_disc) / (2.0 * av);
            double u2 = (-bv - sqrt_disc) / (2.0 * av);
            
            double u_vals[2] = {u1, u2};
            for (int i = 0; i < 2; ++i) {
                if (u_vals[i] >= -eps) {
                    double x_pos = std::sqrt(std::max(0.0, u_vals[i]));
                    results.push_back(SymbolicExpr::number(x_pos));
                    results.push_back(SymbolicExpr::number(-x_pos));
                } else {
                    // u < 0: complex roots x = ±i*sqrt(|u|)
                    auto u_expr = SymbolicExpr::number(u_vals[i]);
                    results.push_back(SymbolicExpr::sqrt(u_expr)->simplify());
                    results.push_back(negate(SymbolicExpr::sqrt(u_expr))->simplify());
                }
            }
        } else {
            // Complex discriminant - both u values are complex
            // Return symbolic form
            auto u_roots = solve_quadratic_internal(a, b, c);
            for (const auto& u : u_roots) {
                results.push_back(SymbolicExpr::sqrt(u)->simplify());
                results.push_back(negate(SymbolicExpr::sqrt(u))->simplify());
            }
        }
        
        return results;
    }
    
    // Symbolic path: Solve au^2 + bu + c = 0 for u
    auto u_roots = solve_quadratic_internal(a, b, c);
    
    // For each u root, x = ±sqrt(u)
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& u : u_roots) {
        auto pos_root = SymbolicExpr::sqrt(u)->simplify();
        auto neg_root = negate(SymbolicExpr::sqrt(u))->simplify();
        results.push_back(pos_root);
        results.push_back(neg_root);
    }
    
    return results;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var) {
    
    // Check if a == 0: delegate to lower-degree solver
    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {
        // Not a quartic, delegate to cubic
        return solve_cubic(b, c, d, e, var);
    }
    
    // Step 0: Check for biquadratic (b=0, d=0)
    auto b_simp = b->simplify();
    auto d_simp = d->simplify();
    if (b_simp->get_number_value_is_zero() && d_simp->get_number_value_is_zero()) {
        return solve_biquadratic(a, c, e, var);
    }
    
    // Full Ferrari's method for general quartic ax^4 + bx^3 + cx^2 + dx + e = 0
    
    // Check if all input coefficients are purely numeric
    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) &&
                       is_purely_numeric(c) && is_purely_numeric(d) && is_purely_numeric(e);
    
    if (all_numeric) {
        // ===== Numeric path =====
        double av = a->to_numeric();
        double bv = b->to_numeric();
        double cv = c->to_numeric();
        double dv = d->to_numeric();
        double ev = e->to_numeric();
        
        // Step 1: Compute depressed quartic y^4 + py^2 + qy + r = 0
        // x = y - bv/(4av)
        double shift_val = bv / (4.0 * av);
        double p_val = (8.0*av*cv - 3.0*bv*bv) / (8.0*av*av);
        double q_val = (bv*bv*bv - 4.0*av*bv*cv + 8.0*av*av*dv) / (8.0*av*av*av);
        double r_val = (-3.0*bv*bv*bv*bv + 256.0*av*av*av*ev - 64.0*av*av*bv*dv + 16.0*av*bv*bv*cv) / (256.0*av*av*av*av);
        
        const double eps = 1e-12;
        
        // Step 2: If q == 0, degenerate to biquadratic y^4 + py^2 + r = 0
        if (std::abs(q_val) < eps) {
            // Solve u^2 + p*u + r = 0 where u = y^2
            double disc = p_val * p_val - 4.0 * r_val;
            double sqrt_disc = std::sqrt(std::abs(disc));
            double u1, u2;
            if (disc >= 0) {
                u1 = (-p_val + sqrt_disc) / 2.0;
                u2 = (-p_val - sqrt_disc) / 2.0;
            } else {
                // Complex u values - use symbolic form
                u1 = (-p_val + sqrt_disc) / 2.0;
                u2 = (-p_val - sqrt_disc) / 2.0;
            }
            
            std::vector<std::shared_ptr<SymbolicExpr>> results;
            double u_vals[2] = {u1, u2};
            for (int i = 0; i < 2; ++i) {
                if (u_vals[i] >= -eps) {
                    double y_pos = std::sqrt(std::max(0.0, u_vals[i]));
                    results.push_back(SymbolicExpr::number(y_pos - shift_val));
                    results.push_back(SymbolicExpr::number(-y_pos - shift_val));
                } else {
                    // u < 0: complex roots y = ±i*sqrt(|u|)
                    // Return symbolic form: sqrt(u) where u is negative
                    auto u_expr = SymbolicExpr::number(u_vals[i]);
                    auto shift_expr = SymbolicExpr::number(shift_val);
                    results.push_back(sub(SymbolicExpr::sqrt(u_expr), shift_expr)->simplify());
                    results.push_back(sub(negate(SymbolicExpr::sqrt(u_expr)), shift_expr)->simplify());
                }
            }
            return results;
        }
        
        // Step 3: Solve resolvent cubic 8m^3 + 8pm^2 + (2p^2 - 8r)m - q^2 = 0
        double rc_a = 8.0;
        double rc_b = 8.0 * p_val;
        double rc_c = 2.0 * p_val * p_val - 8.0 * r_val;
        double rc_d = -(q_val * q_val);
        
        auto cubic_roots = solve_cubic(
            SymbolicExpr::number(rc_a),
            SymbolicExpr::number(rc_b),
            SymbolicExpr::number(rc_c),
            SymbolicExpr::number(rc_d),
            var);
        
        // Pick a positive real root m (required so that sqrt(2m) is real).
        // The resolvent cubic always has at least one positive real root when
        // the original quartic has real coefficients and q != 0.
        // Strategy: collect all finite real roots, prefer the largest positive one.
        double m_val = 0.0;
        bool found_m = false;
        
        // First pass: find the largest positive root
        for (const auto& root : cubic_roots) {
            double val = root->to_numeric();
            if (std::isfinite(val) && val > eps) {
                if (!found_m || val > m_val) {
                    m_val = val;
                    found_m = true;
                }
            }
        }
        
        // Second pass: if no positive root found, try any non-zero finite root
        if (!found_m) {
            for (const auto& root : cubic_roots) {
                double val = root->to_numeric();
                if (std::isfinite(val) && std::abs(val) > eps) {
                    if (!found_m || std::abs(val) > std::abs(m_val)) {
                        m_val = val;
                        found_m = true;
                    }
                }
            }
        }
        
        if (!found_m && !cubic_roots.empty()) {
            m_val = cubic_roots[0]->to_numeric();
            found_m = true;
        }
        
        if (!found_m) return {};
        
        // Step 4: Factor into two quadratics
        // s = sqrt(2m) - if m is negative, we need to handle this carefully
        if (m_val < 0) {
            // When m < 0, sqrt(2m) is imaginary. This shouldn't happen for quartics
            // with all real roots, but can happen when the quartic has complex roots.
            // Use |m| and adjust signs accordingly.
            double abs_2m = std::abs(2.0 * m_val);
            double s_abs = std::sqrt(abs_2m);
            // In this case the factoring produces complex quadratics.
            // Fall back to using the absolute value and accepting complex results.
            // Actually, re-derive: if m < 0, let's just use m as-is with complex arithmetic.
            // For simplicity, return symbolic expressions.
            auto shift_expr = SymbolicExpr::number(shift_val);
            auto m_expr = SymbolicExpr::number(m_val);
            auto p_expr = SymbolicExpr::number(p_val);
            auto q_expr = SymbolicExpr::number(q_val);
            auto two_m_expr = SymbolicExpr::number(2.0 * m_val);
            auto s_expr = SymbolicExpr::sqrt(two_m_expr);
            auto p_half_expr = SymbolicExpr::number(p_val / 2.0);
            auto m_plus_p_half = SymbolicExpr::number(m_val + p_val / 2.0);
            auto q_over_2s = SymbolicExpr::divide(q_expr, SymbolicExpr::multiply(num(2), s_expr));
            
            auto quad1_c_expr = sub(m_plus_p_half, q_over_2s)->simplify();
            auto quad2_c_expr = SymbolicExpr::add(m_plus_p_half, q_over_2s)->simplify();
            
            auto y_roots1 = solve_quadratic_internal(num(1), s_expr, quad1_c_expr);
            auto neg_s_expr = negate(s_expr);
            auto y_roots2 = solve_quadratic_internal(num(1), neg_s_expr, quad2_c_expr);
            
            std::vector<std::shared_ptr<SymbolicExpr>> results;
            for (const auto& y : y_roots1) {
                results.push_back(sub(y, shift_expr)->simplify());
            }
            for (const auto& y : y_roots2) {
                results.push_back(sub(y, shift_expr)->simplify());
            }
            return results;
        }
        
        double s_val = std::sqrt(2.0 * m_val);
        
        // quad1: y^2 + s*y + (m + p/2 - q/(2s)) = 0
        double quad1_c_val = m_val + p_val / 2.0 - q_val / (2.0 * s_val);
        // quad2: y^2 - s*y + (m + p/2 + q/(2s)) = 0
        double quad2_c_val = m_val + p_val / 2.0 + q_val / (2.0 * s_val);
        
        // Step 5: Solve both quadratics
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        
        // Quadratic 1: y^2 + s*y + quad1_c = 0
        double disc1 = s_val * s_val - 4.0 * quad1_c_val;
        if (disc1 >= -eps) {
            double sqrt_disc1 = std::sqrt(std::max(0.0, disc1));
            double y1 = (-s_val + sqrt_disc1) / 2.0;
            double y2 = (-s_val - sqrt_disc1) / 2.0;
            results.push_back(SymbolicExpr::number(y1 - shift_val));
            results.push_back(SymbolicExpr::number(y2 - shift_val));
        } else {
            // Complex roots
            double real_part = -s_val / 2.0 - shift_val;
            double imag_part = std::sqrt(-disc1) / 2.0;
            auto i_unit = SymbolicExpr::sqrt(num(-1));
            results.push_back(SymbolicExpr::add(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part), i_unit))->simplify());
            results.push_back(sub(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part), i_unit))->simplify());
        }
        
        // Quadratic 2: y^2 - s*y + quad2_c = 0
        double disc2 = s_val * s_val - 4.0 * quad2_c_val;
        if (disc2 >= -eps) {
            double sqrt_disc2 = std::sqrt(std::max(0.0, disc2));
            double y3 = (s_val + sqrt_disc2) / 2.0;
            double y4 = (s_val - sqrt_disc2) / 2.0;
            results.push_back(SymbolicExpr::number(y3 - shift_val));
            results.push_back(SymbolicExpr::number(y4 - shift_val));
        } else {
            // Complex roots
            double real_part = s_val / 2.0 - shift_val;
            double imag_part = std::sqrt(-disc2) / 2.0;
            auto i_unit = SymbolicExpr::sqrt(num(-1));
            results.push_back(SymbolicExpr::add(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part), i_unit))->simplify());
            results.push_back(sub(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part), i_unit))->simplify());
        }
        
        return results;
    }
    
    // ===== Symbolic path: coefficients contain parameters =====
    // Return expressions with sqrt/cbrt without evaluating discriminant sign
    
    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto a4 = SymbolicExpr::power(a, num(4));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    auto b4 = SymbolicExpr::power(b, num(4));
    
    // p = (8ac - 3b^2) / (8a^2)
    auto eight_ac = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a, c));
    auto three_b2 = SymbolicExpr::multiply(num(3), b2);
    auto p = SymbolicExpr::divide(sub(eight_ac, three_b2), SymbolicExpr::multiply(num(8), a2))->simplify();
    
    // q = (b^3 - 4abc + 8a^2*d) / (8a^3)
    auto four_abc = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto eight_a2d = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(b3, four_abc), eight_a2d);
    auto q = SymbolicExpr::divide(q_num, SymbolicExpr::multiply(num(8), a3))->simplify();
    
    // r = (-3b^4 + 256a^3*e - 64a^2*b*d + 16a*b^2*c) / (256a^4)
    auto neg3_b4 = SymbolicExpr::multiply(num(-3), b4);
    auto t256_a3e = SymbolicExpr::multiply(num(256), SymbolicExpr::multiply(a3, e));
    auto neg64_a2bd = SymbolicExpr::multiply(num(-64), SymbolicExpr::multiply(a2, SymbolicExpr::multiply(b, d)));
    auto t16_ab2c = SymbolicExpr::multiply(num(16), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b2, c)));
    auto r_num = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::add(neg3_b4, t256_a3e), neg64_a2bd), t16_ab2c);
    auto r = SymbolicExpr::divide(r_num, SymbolicExpr::multiply(num(256), a4))->simplify();
    
    // Shift for back-substitution: x = y - b/(4a)
    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(4), a))->simplify();
    
    // Step 2: If q == 0, degenerate to biquadratic in y: y^4 + py^2 + r = 0
    auto q_simp = q->simplify();
    if (q_simp->get_number_value_is_zero()) {
        // Solve u^2 + p*u + r = 0 where u = y^2
        auto u_roots = solve_quadratic_internal(num(1), p, r);
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (const auto& u : u_roots) {
            // y = ±sqrt(u), then x = y - shift
            auto pos_y = SymbolicExpr::sqrt(u)->simplify();
            auto neg_y = negate(SymbolicExpr::sqrt(u))->simplify();
            results.push_back(sub(pos_y, shift)->simplify());
            results.push_back(sub(neg_y, shift)->simplify());
        }
        return results;
    }
    
    // Step 3: Solve resolvent cubic 8m^3 + 8pm^2 + (2p^2 - 8r)m - q^2 = 0
    auto p2 = SymbolicExpr::power(p, num(2));
    auto q2 = SymbolicExpr::power(q, num(2));
    auto eight_p = SymbolicExpr::multiply(num(8), p);
    auto two_p2_minus_8r = sub(SymbolicExpr::multiply(num(2), p2), SymbolicExpr::multiply(num(8), r))->simplify();
    auto neg_q2 = negate(q2)->simplify();
    
    auto cubic_roots = solve_cubic(num(8), eight_p, two_p2_minus_8r, neg_q2, var);
    
    // Pick a non-zero root m from the resolvent cubic
    // For symbolic coefficients, just use the first root
    std::shared_ptr<SymbolicExpr> m = nullptr;
    if (!cubic_roots.empty()) {
        m = cubic_roots[0];
    }
    
    if (!m) return {};
    
    // Step 4: Factor the depressed quartic into two quadratics
    // s = sqrt(2m)
    auto two_m = SymbolicExpr::multiply(num(2), m);
    auto s = SymbolicExpr::sqrt(two_m)->simplify();
    
    // p_half = p/2
    auto p_half = SymbolicExpr::divide(p, num(2))->simplify();
    
    // m + p/2
    auto m_plus_p_half = SymbolicExpr::add(m, p_half)->simplify();
    
    // q/(2s)
    auto q_over_2s = SymbolicExpr::divide(q, SymbolicExpr::multiply(num(2), s))->simplify();
    
    // quad1: y^2 + s*y + (m + p/2 - q/(2s)) = 0
    auto quad1_c = sub(m_plus_p_half, q_over_2s)->simplify();
    
    // quad2: y^2 - s*y + (m + p/2 + q/(2s)) = 0
    auto neg_s = negate(s)->simplify();
    auto quad2_c = SymbolicExpr::add(m_plus_p_half, q_over_2s)->simplify();
    
    // Step 5: Solve both quadratics
    auto y_roots1 = solve_quadratic_internal(num(1), s, quad1_c);
    auto y_roots2 = solve_quadratic_internal(num(1), neg_s, quad2_c);
    
    // Step 6: Back-substitute x_k = y_k - b/(4a)
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        results.push_back(sub(y, shift)->simplify());
    }
    for (const auto& y : y_roots2) {
        results.push_back(sub(y, shift)->simplify());
    }
    
    return results;
}

// ============================================================
// find_rational_roots - Rational Root Theorem with multiplicity
// ============================================================

// Helper: find all positive divisors of a BigInt (absolute value).
// For very large numbers, we cap the search to avoid combinatorial explosion.
static std::vector<BigInt> positive_divisors(const BigInt& n) {
    std::vector<BigInt> divs;
    if (n.is_zero()) {
        // Zero has no meaningful divisors for our purposes
        return divs;
    }
    BigInt abs_n = n.Abs();
    
    // Trial division up to sqrt(abs_n).
    // For practical polynomial solving, coefficients are usually small.
    // We cap at 1000 to avoid extremely long loops for huge coefficients.
    BigInt i(1);
    BigInt limit(1000);
    // We'll collect divisors by trial division
    // If abs_n is larger than limit^2, we still try up to limit
    while (i * i <= abs_n && i <= limit) {
        BigInt rem = abs_n % i;
        if (rem.is_zero()) {
            divs.push_back(i);
            BigInt counterpart = abs_n / i;
            if (counterpart != i) {
                divs.push_back(counterpart);
            }
        }
        i = i + BigInt(1);
    }
    
    // Sort divisors for deterministic candidate ordering
    std::sort(divs.begin(), divs.end());
    return divs;
}

// Square-free factorization: decompose poly into coprime factors f₁·f₂²·f₃³·…
// Uses the classical Yun's algorithm over Q.
// Returns pairs of (factor, multiplicity) where each factor is square-free and monic.
std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly) {
    
    std::vector<std::pair<Polynomial<Rational>, int>> factors;
    
    if (poly.is_zero() || poly.degree() <= 0) {
        return factors;
    }
    
    // Work with a monic copy for cleaner results
    Polynomial<Rational> f = poly.make_monic();
    
    // If degree is 1, it's trivially square-free
    if (f.degree() == 1) {
        factors.push_back({f, 1});
        return factors;
    }
    
    // Step 1: Compute g = gcd(f, f')
    Polynomial<Rational> f_prime = f.differentiate();
    Polynomial<Rational> g = Polynomial<Rational>::gcd(f, f_prime);
    
    // Step 2: If gcd is constant (degree 0), f is already square-free
    if (g.degree() <= 0) {
        factors.push_back({f, 1});
        return factors;
    }
    
    // Step 3: Yun's algorithm
    // w = f / g (the "square-free part" but may share factors with g)
    auto [w, w_rem] = f.div_mod(g);
    
    int multiplicity = 1;
    
    while (w.degree() > 0) {
        // y = gcd(g, w)
        Polynomial<Rational> y = Polynomial<Rational>::gcd(g, w);
        
        // z = w / y — this is the square-free factor with current multiplicity
        auto [z, z_rem] = w.div_mod(y);
        
        if (z.degree() > 0) {
            // z is a non-trivial square-free factor with this multiplicity
            factors.push_back({z.make_monic(), multiplicity});
        }
        
        // Update for next iteration
        // g = g / y, w = y, multiplicity++
        auto [g_next, g_rem] = g.div_mod(y);
        g = g_next;
        w = y;
        multiplicity++;
    }
    
    // If g still has degree > 0 after the loop, it contributes a factor
    // with the current multiplicity (this handles the case where all factors
    // have the same multiplicity)
    if (g.degree() > 0) {
        factors.push_back({g.make_monic(), multiplicity});
    }
    
    return factors;
}

// ============================================================
// solve_by_factoring - Orchestrator for polynomial preprocessing
// Pipeline: extract rational roots → square-free factorization →
// for each factor, dispatch to closed-form solver if degree ≤ 4,
// otherwise emit RootOf placeholders via make_rootof_solutions.
// Preserves multiplicities so the returned vector size matches
// the original polynomial degree.
// ============================================================

// Helper: solve a Polynomial<SymbolicPolyCoeff> of degree 1-4 using closed-form
static std::vector<std::shared_ptr<SymbolicExpr>> solve_closed_form_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var)
{
    int deg = poly.degree();
    if (deg <= 0) return {};

    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {
        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        return { SymbolicExpr::divide(neg_b, a)->simplify() };
    } else if (deg == 2) {
        return solve_quadratic_internal(get_coeff(2), get_coeff(1), get_coeff(0));
    } else if (deg == 3) {
        return solve_cubic(get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var);
    } else if (deg == 4) {
        return solve_quartic(get_coeff(4), get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var);
    }

    return {};
}

// Helper: try to convert a SymbolicPolyCoeff polynomial to Rational polynomial
// Returns true if all coefficients are numeric (no free variables)
static bool convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly)
{
    out_poly = Polynomial<Rational>(sym_poly.variable_name);
    out_poly.coeffs.resize(sym_poly.coeffs.size(), Rational(0));

    for (size_t i = 0; i < sym_poly.coeffs.size(); ++i) {
        auto coeff_expr = sym_poly.coeffs[i].val;
        if (!coeff_expr) {
            out_poly.coeffs[i] = Rational(0);
            continue;
        }

        auto simplified = coeff_expr->simplify();

        // Check if the coefficient has free variables (parametric)
        if (simplified->root) {
            struct VarCheck : public SymbolicVisitor {
                bool found = false;
                void visit(NumberNode&) override {}
                void visit(VariableNode&) override { found = true; }
                void visit(AddNode& n) override { for (auto& op : n.operands) { if (found) return; op->accept(*this); } }
                void visit(MultiplyNode& n) override { for (auto& op : n.operands) { if (found) return; op->accept(*this); } }
                void visit(PowerNode& n) override { n.base->accept(*this); if (!found) n.exponent->accept(*this); }
                void visit(FunctionNode& n) override { for (auto& arg : n.arguments) { if (found) return; arg->accept(*this); } }
                void visit(MatrixNode&) override {}
            } checker;
            simplified->root->accept(checker);
            if (checker.found) return false;
        }

        // Try to extract a Rational value
        if (auto num_node = std::dynamic_pointer_cast<NumberNode>(simplified->root)) {
            if (std::holds_alternative<Rational>(num_node->value)) {
                out_poly.coeffs[i] = std::get<Rational>(num_node->value);
            } else if (std::holds_alternative<BigInt>(num_node->value)) {
                out_poly.coeffs[i] = Rational(std::get<BigInt>(num_node->value));
            } else {
                out_poly.coeffs[i] = Rational::from_double(std::get<double>(num_node->value));
            }
        } else if (simplified->is_zero()) {
            out_poly.coeffs[i] = Rational(0);
        } else if (simplified->is_one()) {
            out_poly.coeffs[i] = Rational(1);
        } else {
            // Try numeric evaluation as last resort
            double val = simplified->to_numeric();
            if (std::isnan(val) || std::isinf(val)) {
                return false;
            }
            out_poly.coeffs[i] = Rational::from_double(val);
        }
    }

    out_poly.trim();
    return true;
}

// Helper: convert a Rational polynomial factor to a SymbolicPolyCoeff polynomial
static Polynomial<SymbolicPolyCoeff> rational_to_symbolic_poly(
    const Polynomial<Rational>& rat_poly)
{
    std::vector<SymbolicPolyCoeff> sym_coeffs;
    sym_coeffs.reserve(rat_poly.coeffs.size());
    for (const auto& c : rat_poly.coeffs) {
        sym_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(c)));
    }
    return Polynomial<SymbolicPolyCoeff>(sym_coeffs, rat_poly.variable_name);
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    if (poly.is_zero() || poly.degree() < 1) return {};

    // If degree <= 4, just use closed-form directly
    if (poly.degree() <= 4) {
        return solve_closed_form_poly(poly, var);
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;

    // Step 1: Try to convert to Rational polynomial
    Polynomial<Rational> rat_poly;
    if (!convert_to_rational_poly(poly, rat_poly)) {
        // Cannot convert to rational - generate RootOf for the whole polynomial
        return make_rootof_solutions(poly, var);
    }

    // Step 2: Extract rational roots with multiplicity
    auto rational_roots = find_rational_roots(rat_poly);

    // Add rational roots to results
    for (const auto& r : rational_roots) {
        results.push_back(SymbolicExpr::number(r));
    }

    // Step 3: Divide out the rational roots to get the remaining quotient
    Polynomial<Rational> quotient = rat_poly;
    for (const auto& r : rational_roots) {
        Polynomial<Rational> linear_factor({-r, Rational(1)}, var);
        auto [q, rem] = quotient.div_mod(linear_factor);
        quotient = q;
    }

    // If quotient is constant or zero, we're done
    if (quotient.degree() <= 0) {
        return results;
    }

    // If quotient degree <= 4, solve directly with closed-form
    if (quotient.degree() <= 4) {
        auto sym_quotient = rational_to_symbolic_poly(quotient);
        auto factor_roots = solve_closed_form_poly(sym_quotient, var);
        results.insert(results.end(), factor_roots.begin(), factor_roots.end());
        return results;
    }

    // Step 4: Apply square-free factorization on the quotient
    auto sqfree_factors = square_free_factorization(quotient);

    // Step 5: For each (factor, multiplicity) pair, solve or generate RootOf
    for (const auto& [factor, multiplicity] : sqfree_factors) {
        if (factor.degree() <= 0) continue;

        std::vector<std::shared_ptr<SymbolicExpr>> factor_solutions;

        if (factor.degree() <= 4) {
            // Solve using closed-form
            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = solve_closed_form_poly(sym_factor, var);
        } else {
            // Generate RootOf placeholders
            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = make_rootof_solutions(sym_factor, var);
        }

        // Repeat each solution `multiplicity` times to preserve total root count
        for (int m = 0; m < multiplicity; ++m) {
            results.insert(results.end(), factor_solutions.begin(), factor_solutions.end());
        }
    }

    return results;
}

std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly) {
    std::vector<Rational> roots;
    
    if (poly.is_zero() || poly.degree() < 1) {
        return roots;
    }
    
    // Work on a mutable copy of the polynomial
    Polynomial<Rational> current = poly;
    
    // Special case: factor out maximum power of x when constant term is zero
    // (Requirement 4.7)
    while (current.degree() >= 1 && current.coeffs[0] == Rational(0)) {
        roots.push_back(Rational(0));
        // Divide by x: shift coefficients down
        std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
        current = Polynomial<Rational>(new_coeffs, current.variable_name);
    }
    
    // Now search for non-zero rational roots using the Rational Root Theorem
    // Continue until no more rational roots or degree drops to <= 4
    while (current.degree() >= 1) {
        // Get constant term (a₀) and leading coefficient (aₙ)
        Rational a0 = current.coeffs[0];
        Rational an = current.lead_coeff();
        
        if (a0 == Rational(0)) {
            // Constant term became zero again after division - factor out x
            roots.push_back(Rational(0));
            std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
            current = Polynomial<Rational>(new_coeffs, current.variable_name);
            continue;
        }
        
        // For the rational root theorem, we need integer numerators/denominators.
        // a0 = a0_num / a0_den, an = an_num / an_den
        // If p/q is a root of the polynomial with integer coefficients,
        // then p | a0 and q | an (for the integer-coefficient version).
        // For rational coefficients, we first convert to integer coefficients
        // by multiplying through by the LCM of all denominators.
        
        // Compute LCM of all coefficient denominators to get integer polynomial
        BigInt lcm_den(1);
        for (const auto& c : current.coeffs) {
            BigInt d = c.get_denominator();
            lcm_den = BigInt::lcm(lcm_den, d);
        }
        
        // Integer polynomial coefficients: coeff[i] * lcm_den
        // The constant term of the integer poly = a0 * lcm_den
        // The leading coeff of the integer poly = an * lcm_den
        BigInt int_a0 = (a0 * Rational(lcm_den)).get_numerator();
        BigInt int_an = (an * Rational(lcm_den)).get_numerator();
        
        // Get positive divisors of |int_a0| and |int_an|
        std::vector<BigInt> p_divs = positive_divisors(int_a0);
        std::vector<BigInt> q_divs = positive_divisors(int_an);
        
        if (p_divs.empty() || q_divs.empty()) {
            break; // Cannot enumerate candidates
        }
        
        // Enumerate candidates ±p/q
        bool found_root = false;
        Rational found_r;
        
        for (const auto& p : p_divs) {
            for (const auto& q : q_divs) {
                // Try +p/q and -p/q
                Rational candidate_pos(p, q);
                Rational candidate_neg(-p, q);
                
                if (current.eval(candidate_pos) == Rational(0)) {
                    found_root = true;
                    found_r = candidate_pos;
                    break;
                }
                if (current.eval(candidate_neg) == Rational(0)) {
                    found_root = true;
                    found_r = candidate_neg;
                    break;
                }
            }
            if (found_root) break;
        }
        
        if (!found_root) {
            break; // No more rational roots in the current polynomial
        }
        
        // Found a root. Determine its full multiplicity by repeated division.
        // Divide by (x - r) repeatedly while the quotient still vanishes at r.
        Polynomial<Rational> linear_factor({-found_r, Rational(1)}, current.variable_name);
        
        while (current.degree() >= 1) {
            auto [quotient, remainder] = current.div_mod(linear_factor);
            // Check if division was exact (remainder is zero)
            if (!remainder.is_zero()) {
                break;
            }
            roots.push_back(found_r);
            current = quotient;
            
            // Check if the quotient still vanishes at r
            if (current.degree() < 1 || current.eval(found_r) != Rational(0)) {
                break;
            }
        }
        
        // If degree dropped to <= 4, stop searching for rational roots
        // (the caller can use closed-form solvers for the remainder)
        if (current.degree() <= 4) {
            break;
        }
    }
    
    return roots;
}

} // namespace lamina
