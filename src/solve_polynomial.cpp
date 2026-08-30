#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"
#include "root_of_utils.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include <cmath>
#include <algorithm>
#include <optional>
#include <stdexcept>

namespace lamina {

static bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return true;

    struct VarDetector : public lamina::detail::RecursiveSymbolicVisitor {
        bool has_var = false;
        void visit(const NumberNode&) override {}
        void visit(const VariableNode&) override { has_var = true; }
        void visit(const AddNode& n) override { for (auto& op : n.operands()) { if (has_var) return; op->accept(*this); } }
        void visit(const MultiplyNode& n) override { for (auto& op : n.operands()) { if (has_var) return; op->accept(*this); } }
        void visit(const PowerNode& n) override { n.base()->accept(*this); if (!has_var) n.exponent()->accept(*this); }
        void visit(const FunctionNode& n) override { for (auto& arg : n.arguments()) { if (has_var) return; arg->accept(*this); } }
        void visit(const MatrixNode&) override {}
        void visit(const RelationalNode& n) override { n.left()->accept(*this); if (!has_var) n.right()->accept(*this); }
        void visit(const LogicalNode& n) override { n.left()->accept(*this); if (!has_var && n.right()) n.right()->accept(*this); }
        void visit(const PiecewiseNode& n) override {
            for (const auto& branch : n.branches()) {
                if (has_var) return;
                branch.expression->accept(*this);
                if (!has_var) branch.condition->accept(*this);
            }
            if (!has_var && n.default_expr()) n.default_expr()->accept(*this);
        }
        void visit(const SummationNode& n) override { n.body()->accept(*this); if (!has_var) n.lower_bound()->accept(*this); if (!has_var) n.upper_bound()->accept(*this); }
        void visit(const ProductNode& n) override { n.body()->accept(*this); if (!has_var) n.lower_bound()->accept(*this); if (!has_var) n.upper_bound()->accept(*this); }
        void visit(const TransformNode& n) override { n.body()->accept(*this); }
        void visit(const QuantifierNode& n) override { n.domain()->accept(*this); if (!has_var) n.predicate()->accept(*this); }
        void visit(const SetBuilderNode& n) override { n.domain()->accept(*this); if (!has_var) n.predicate()->accept(*this); }
        void visit(const ComplexNode& n) override { n.real()->accept(*this); if (!has_var) n.imag()->accept(*this); }
        void visit(const FiniteSetNode& n) override { for (const auto& e : n.elements()) { if (has_var) return; e->accept(*this); } }
        void visit(const IntervalNode& n) override { n.lower()->accept(*this); if (!has_var) n.upper()->accept(*this); }
        void visit(const MembershipNode& n) override { n.element()->accept(*this); if (!has_var) n.set()->accept(*this); }
        void visit(const QuantityNode& n) override { n.value()->accept(*this); }
    } detector;

    lamina::detail::node(expr)->accept(detector);
    return !detector.has_var;
}

static std::shared_ptr<SymbolicExpr> cbrt_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3)));
}

[[maybe_unused]] static std::shared_ptr<SymbolicExpr> acos_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(x)}));
}

static std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
}

static std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    return SymbolicExpr::add(a, negate(b));
}

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
[[maybe_unused]] static std::shared_ptr<SymbolicExpr> num_r(int p, int q) { return SymbolicExpr::number(Rational(p, q)); }

static std::optional<double> finite_numeric_value(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return std::nullopt;
    ComputationContext context;
    auto evaluated = evaluate_numeric(*expr, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    auto simplified = expr->simplify();
    if (simplified && lamina::detail::node(simplified)) {
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
            if (std::holds_alternative<Rational>(num->value())) {
                const auto& value = std::get<Rational>(num->value());
                if (!value.get_numerator().is_zero() && evaluated.value().value == 0.0) {
                    return std::nullopt;
                }
            } else if (std::holds_alternative<BigInt>(num->value())) {
                const auto& value = std::get<BigInt>(num->value());
                if (!value.is_zero() && evaluated.value().value == 0.0) {
                    return std::nullopt;
                }
            }
        }
    }
    return evaluated.value().value;
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c);

static bool convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b) {

    auto neg_b = negate(b);
    return { SymbolicExpr::divide(neg_b, a)->simplify() };
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {

    auto b2 = SymbolicExpr::power(b, num(2));
    auto four_ac = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, c));
    auto delta = sub(b2, four_ac)->simplify();

    auto neg_b = negate(b);
    auto two_a = SymbolicExpr::multiply(num(2), a);
    if (delta && delta->is_zero()) {
        return { SymbolicExpr::divide(neg_b, two_a)->simplify() };
    }

    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
    auto x2 = SymbolicExpr::divide(sub(neg_b, sqrt_delta), two_a)->simplify();

    return { x1, x2 };
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string&) {

    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {

        auto b_simp = b->simplify();
        if (b_simp->get_number_value_is_zero()) {

            return solve_linear_internal(c, d);
        }

        return solve_quadratic_internal(b, c, d);
    }

    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));

    auto three_ac = SymbolicExpr::multiply(num(3), SymbolicExpr::multiply(a, c));
    auto p_num = sub(three_ac, b2);
    auto p_den = SymbolicExpr::multiply(num(3), a2);
    auto p = SymbolicExpr::divide(p_num, p_den)->simplify();

    auto two_b3 = SymbolicExpr::multiply(num(2), b3);
    auto nine_abc = SymbolicExpr::multiply(num(9), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto twentyseven_a2d = SymbolicExpr::multiply(num(27), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(two_b3, nine_abc), twentyseven_a2d);
    auto q_den = SymbolicExpr::multiply(num(27), a3);
    auto q = SymbolicExpr::divide(q_num, q_den)->simplify();

    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(3), a))->simplify();

    bool all_numeric = is_purely_numeric(p) && is_purely_numeric(q);

    std::vector<std::shared_ptr<SymbolicExpr>> roots;

    auto p_checked = finite_numeric_value(p);
    auto q_checked = finite_numeric_value(q);
    auto shift_checked = finite_numeric_value(shift);

    if (all_numeric && p_checked && q_checked && shift_checked) {

        double p_val = *p_checked;
        double q_val = *q_checked;
        double shift_val = *shift_checked;

        double D = (q_val / 2.0) * (q_val / 2.0) + (p_val / 3.0) * (p_val / 3.0) * (p_val / 3.0);

        const double eps = 1e-12;

        if (std::abs(p_val) < eps && std::abs(q_val) < eps) {

            auto root = SymbolicExpr::number(-shift_val);
            roots = { root, root, root };
        } else if (std::abs(D) < eps) {

            double q_half_val = q_val / 2.0;
            double cbrt_q_half_val = std::cbrt(q_half_val);

            double t1_val = -2.0 * cbrt_q_half_val;

            double t2_val = cbrt_q_half_val;

            auto x1 = SymbolicExpr::number(t1_val - shift_val);
            auto x2 = SymbolicExpr::number(t2_val - shift_val);

            roots = { x1, x2, x2 };
        } else if (D > eps) {

            double neg_q_half_val = -q_val / 2.0;
            double sqrt_D_val = std::sqrt(D);

            double u_arg_val = neg_q_half_val + sqrt_D_val;
            double v_arg_val = neg_q_half_val - sqrt_D_val;

            double u_val = std::cbrt(u_arg_val);
            double v_val = std::cbrt(v_arg_val);

            double t1_val = u_val + v_val;
            double x1_val = t1_val - shift_val;

            double real_part = -(u_val + v_val) / 2.0 - shift_val;
            double imag_part_val = std::sqrt(3.0) * (u_val - v_val) / 2.0;

            auto x1 = SymbolicExpr::number(x1_val);

            auto i_unit = SymbolicExpr::sqrt(num(-1));
            auto x2 = SymbolicExpr::add(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part_val), i_unit))->simplify();
            auto x3 = sub(
                SymbolicExpr::number(real_part),
                SymbolicExpr::multiply(SymbolicExpr::number(imag_part_val), i_unit))->simplify();

            roots = { x1, x2, x3 };
        } else {

            double neg_p3_over_27 = -(p_val * p_val * p_val) / 27.0;
            double r_val = std::sqrt(neg_p3_over_27);
            double cos_arg = -q_val / (2.0 * r_val);

            if (cos_arg > 1.0) cos_arg = 1.0;
            if (cos_arg < -1.0) cos_arg = -1.0;
            double theta_val = std::acos(cos_arg);
            double two_cbrt_r = 2.0 * std::cbrt(r_val);

            for (int k = 0; k < 3; ++k) {
                double angle = (theta_val + 2.0 * k * LMMC_CONST_PI) / 3.0;
                double t_k_val = two_cbrt_r * std::cos(angle);
                double x_k_val = t_k_val - shift_val;
                roots.push_back(SymbolicExpr::number(x_k_val));
            }
        }
    } else {

        auto q_half = SymbolicExpr::divide(q, num(2));
        auto p_third = SymbolicExpr::divide(p, num(3));
        auto D_expr = SymbolicExpr::add(
            SymbolicExpr::power(q_half, num(2)),
            SymbolicExpr::power(p_third, num(3)))->simplify();

        auto neg_q_half = negate(q_half)->simplify();
        auto sqrt_D = SymbolicExpr::sqrt(D_expr);

        auto u = cbrt_expr(SymbolicExpr::add(neg_q_half, sqrt_D))->simplify();
        auto v = cbrt_expr(sub(neg_q_half, sqrt_D))->simplify();

        auto t1 = SymbolicExpr::add(u, v);
        auto x1 = sub(t1, shift)->simplify();

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
    const std::string&) {

    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) && is_purely_numeric(c);

    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);

    if (all_numeric && a_checked && b_checked && c_checked && *a_checked != 0.0) {

        double av = *a_checked;
        double bv = *b_checked;
        double cv = *c_checked;

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

                    auto u_expr = SymbolicExpr::number(u_vals[i]);
                    results.push_back(SymbolicExpr::sqrt(u_expr)->simplify());
                    results.push_back(negate(SymbolicExpr::sqrt(u_expr))->simplify());
                }
            }
        } else {

            auto u_roots = solve_quadratic_internal(a, b, c);
            for (const auto& u : u_roots) {
                results.push_back(SymbolicExpr::sqrt(u)->simplify());
                results.push_back(negate(SymbolicExpr::sqrt(u))->simplify());
            }
        }

        return results;
    }

    auto u_roots = solve_quadratic_internal(a, b, c);

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

    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {

        return solve_cubic(b, c, d, e, var);
    }

    {
        Polynomial<SymbolicPolyCoeff> symbolic_poly(
            {SymbolicPolyCoeff(e), SymbolicPolyCoeff(d), SymbolicPolyCoeff(c),
             SymbolicPolyCoeff(b), SymbolicPolyCoeff(a)},
            var);
        Polynomial<Rational> rational_poly;
        if (convert_to_rational_poly(symbolic_poly, rational_poly) &&
            rational_poly.degree() == 4) {
            auto rational_roots = find_rational_roots(rational_poly);
            if (rational_roots.size() == 4) {
                std::vector<std::shared_ptr<SymbolicExpr>> proven_roots;
                proven_roots.reserve(rational_roots.size());
                for (const auto& root : rational_roots) {
                    proven_roots.push_back(SymbolicExpr::number(root));
                }
                return proven_roots;
            }
        }
    }

    auto b_simp = b->simplify();
    auto d_simp = d->simplify();
    if (b_simp->get_number_value_is_zero() && d_simp->get_number_value_is_zero()) {
        return solve_biquadratic(a, c, e, var);
    }

    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) &&
                       is_purely_numeric(c) && is_purely_numeric(d) && is_purely_numeric(e);

    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);
    auto d_checked = finite_numeric_value(d);
    auto e_checked = finite_numeric_value(e);

    if (all_numeric && a_checked && b_checked && c_checked && d_checked && e_checked &&
        *a_checked != 0.0) {

        double av = *a_checked;
        double bv = *b_checked;
        double cv = *c_checked;
        double dv = *d_checked;
        double ev = *e_checked;

        double shift_val = bv / (4.0 * av);
        double p_val = (8.0*av*cv - 3.0*bv*bv) / (8.0*av*av);
        double q_val = (bv*bv*bv - 4.0*av*bv*cv + 8.0*av*av*dv) / (8.0*av*av*av);
        double r_val = (-3.0*bv*bv*bv*bv + 256.0*av*av*av*ev - 64.0*av*av*bv*dv + 16.0*av*bv*bv*cv) / (256.0*av*av*av*av);

        if (!std::isfinite(shift_val) || !std::isfinite(p_val) ||
            !std::isfinite(q_val) || !std::isfinite(r_val)) {
            all_numeric = false;
        }

        if (all_numeric) {

        const double eps = 1e-12;

        if (std::abs(q_val) < eps) {

            double disc = p_val * p_val - 4.0 * r_val;
            double sqrt_disc = std::sqrt(std::abs(disc));
            double u1, u2;
            if (disc >= 0) {
                u1 = (-p_val + sqrt_disc) / 2.0;
                u2 = (-p_val - sqrt_disc) / 2.0;
            } else {

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

                    auto u_expr = SymbolicExpr::number(u_vals[i]);
                    auto shift_expr = SymbolicExpr::number(shift_val);
                    results.push_back(sub(SymbolicExpr::sqrt(u_expr), shift_expr)->simplify());
                    results.push_back(sub(negate(SymbolicExpr::sqrt(u_expr)), shift_expr)->simplify());
                }
            }
            return results;
        }

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

        double m_val = 0.0;
        bool found_m = false;

        for (const auto& root : cubic_roots) {
            auto maybe_val = finite_numeric_value(root);
            if (maybe_val && *maybe_val > eps) {
                if (!found_m || *maybe_val > m_val) {
                    m_val = *maybe_val;
                    found_m = true;
                }
            }
        }

        if (!found_m) {
            for (const auto& root : cubic_roots) {
                auto maybe_val = finite_numeric_value(root);
                if (maybe_val && std::abs(*maybe_val) > eps) {
                    if (!found_m || std::abs(*maybe_val) > std::abs(m_val)) {
                        m_val = *maybe_val;
                        found_m = true;
                    }
                }
            }
        }

        if (!found_m && !cubic_roots.empty()) {
            auto maybe_val = finite_numeric_value(cubic_roots[0]);
            if (maybe_val) {
                m_val = *maybe_val;
                found_m = true;
            }
        }

        if (!found_m) return {};

        if (m_val < 0) {

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

        double quad1_c_val = m_val + p_val / 2.0 - q_val / (2.0 * s_val);

        double quad2_c_val = m_val + p_val / 2.0 + q_val / (2.0 * s_val);

        std::vector<std::shared_ptr<SymbolicExpr>> results;

        double disc1 = s_val * s_val - 4.0 * quad1_c_val;
        if (disc1 >= -eps) {
            double sqrt_disc1 = std::sqrt(std::max(0.0, disc1));
            double y1 = (-s_val + sqrt_disc1) / 2.0;
            double y2 = (-s_val - sqrt_disc1) / 2.0;
            results.push_back(SymbolicExpr::number(y1 - shift_val));
            results.push_back(SymbolicExpr::number(y2 - shift_val));
        } else {

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

        double disc2 = s_val * s_val - 4.0 * quad2_c_val;
        if (disc2 >= -eps) {
            double sqrt_disc2 = std::sqrt(std::max(0.0, disc2));
            double y3 = (s_val + sqrt_disc2) / 2.0;
            double y4 = (s_val - sqrt_disc2) / 2.0;
            results.push_back(SymbolicExpr::number(y3 - shift_val));
            results.push_back(SymbolicExpr::number(y4 - shift_val));
        } else {

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
    }

    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto a4 = SymbolicExpr::power(a, num(4));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    auto b4 = SymbolicExpr::power(b, num(4));

    auto eight_ac = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a, c));
    auto three_b2 = SymbolicExpr::multiply(num(3), b2);
    auto p = SymbolicExpr::divide(sub(eight_ac, three_b2), SymbolicExpr::multiply(num(8), a2))->simplify();

    auto four_abc = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto eight_a2d = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(b3, four_abc), eight_a2d);
    auto q = SymbolicExpr::divide(q_num, SymbolicExpr::multiply(num(8), a3))->simplify();

    auto neg3_b4 = SymbolicExpr::multiply(num(-3), b4);
    auto t256_a3e = SymbolicExpr::multiply(num(256), SymbolicExpr::multiply(a3, e));
    auto neg64_a2bd = SymbolicExpr::multiply(num(-64), SymbolicExpr::multiply(a2, SymbolicExpr::multiply(b, d)));
    auto t16_ab2c = SymbolicExpr::multiply(num(16), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b2, c)));
    auto r_num = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::add(neg3_b4, t256_a3e), neg64_a2bd), t16_ab2c);
    auto r = SymbolicExpr::divide(r_num, SymbolicExpr::multiply(num(256), a4))->simplify();

    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(4), a))->simplify();

    auto q_simp = q->simplify();
    if (q_simp->get_number_value_is_zero()) {

        auto u_roots = solve_quadratic_internal(num(1), p, r);
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (const auto& u : u_roots) {

            auto pos_y = SymbolicExpr::sqrt(u)->simplify();
            auto neg_y = negate(SymbolicExpr::sqrt(u))->simplify();
            results.push_back(sub(pos_y, shift)->simplify());
            results.push_back(sub(neg_y, shift)->simplify());
        }
        return results;
    }

    auto p2 = SymbolicExpr::power(p, num(2));
    auto q2 = SymbolicExpr::power(q, num(2));
    auto eight_p = SymbolicExpr::multiply(num(8), p);
    auto two_p2_minus_8r = sub(SymbolicExpr::multiply(num(2), p2), SymbolicExpr::multiply(num(8), r))->simplify();
    auto neg_q2 = negate(q2)->simplify();

    auto cubic_roots = solve_cubic(num(8), eight_p, two_p2_minus_8r, neg_q2, var);

    std::shared_ptr<SymbolicExpr> m = nullptr;
    if (!cubic_roots.empty()) {
        m = cubic_roots[0];
    }

    if (!m) return {};

    auto two_m = SymbolicExpr::multiply(num(2), m);
    auto s = SymbolicExpr::sqrt(two_m)->simplify();

    auto p_half = SymbolicExpr::divide(p, num(2))->simplify();

    auto m_plus_p_half = SymbolicExpr::add(m, p_half)->simplify();

    auto q_over_2s = SymbolicExpr::divide(q, SymbolicExpr::multiply(num(2), s))->simplify();

    auto quad1_c = sub(m_plus_p_half, q_over_2s)->simplify();

    auto neg_s = negate(s)->simplify();
    auto quad2_c = SymbolicExpr::add(m_plus_p_half, q_over_2s)->simplify();

    auto y_roots1 = solve_quadratic_internal(num(1), s, quad1_c);
    auto y_roots2 = solve_quadratic_internal(num(1), neg_s, quad2_c);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        results.push_back(sub(y, shift)->simplify());
    }
    for (const auto& y : y_roots2) {
        results.push_back(sub(y, shift)->simplify());
    }

    return results;
}

static std::vector<BigInt> positive_divisors(const BigInt& n) {
    std::vector<BigInt> divs;
    if (n.is_zero()) {

        return divs;
    }
    BigInt abs_n = n.Abs();

    BigInt i(1);
    BigInt limit(1000);

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

    std::sort(divs.begin(), divs.end());
    return divs;
}

std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly) {

    std::vector<std::pair<Polynomial<Rational>, int>> factors;

    if (poly.is_zero() || poly.degree() <= 0) {
        return factors;
    }

    Polynomial<Rational> f = poly.make_monic();

    if (f.degree() == 1) {
        factors.push_back({f, 1});
        return factors;
    }

    Polynomial<Rational> f_prime = f.differentiate();
    Polynomial<Rational> g = Polynomial<Rational>::gcd(f, f_prime);

    if (g.degree() <= 0) {
        factors.push_back({f, 1});
        return factors;
    }

    auto [w, w_rem] = f.div_mod(g);

    int multiplicity = 1;

    while (w.degree() > 0) {

        Polynomial<Rational> y = Polynomial<Rational>::gcd(g, w);

        auto [z, z_rem] = w.div_mod(y);

        if (z.degree() > 0) {

            factors.push_back({z.make_monic(), multiplicity});
        }

        auto [g_next, g_rem] = g.div_mod(y);
        g = g_next;
        w = y;
        multiplicity++;
    }

    if (g.degree() > 0) {
        factors.push_back({g.make_monic(), multiplicity});
    }

    return factors;
}

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

        if (lamina::detail::node(simplified)) {
            struct VarCheck : public lamina::detail::RecursiveSymbolicVisitor {
                bool found = false;
                void visit(const NumberNode&) override {}
                void visit(const VariableNode&) override { found = true; }
                void visit(const AddNode& n) override { for (auto& op : n.operands()) { if (found) return; op->accept(*this); } }
                void visit(const MultiplyNode& n) override { for (auto& op : n.operands()) { if (found) return; op->accept(*this); } }
                void visit(const PowerNode& n) override { n.base()->accept(*this); if (!found) n.exponent()->accept(*this); }
                void visit(const FunctionNode& n) override { for (auto& arg : n.arguments()) { if (found) return; arg->accept(*this); } }
                void visit(const MatrixNode&) override {}
                void visit(const RelationalNode& n) override { n.left()->accept(*this); if (!found) n.right()->accept(*this); }
                void visit(const LogicalNode& n) override { n.left()->accept(*this); if (!found && n.right()) n.right()->accept(*this); }
                void visit(const PiecewiseNode& n) override {
                    for (const auto& branch : n.branches()) {
                        if (found) return;
                        branch.expression->accept(*this);
                        if (!found) branch.condition->accept(*this);
                    }
                    if (!found && n.default_expr()) n.default_expr()->accept(*this);
                }
                void visit(const SummationNode& n) override { n.body()->accept(*this); if (!found) n.lower_bound()->accept(*this); if (!found) n.upper_bound()->accept(*this); }
                void visit(const ProductNode& n) override { n.body()->accept(*this); if (!found) n.lower_bound()->accept(*this); if (!found) n.upper_bound()->accept(*this); }
                void visit(const TransformNode& n) override { n.body()->accept(*this); }
                void visit(const QuantifierNode& n) override { n.domain()->accept(*this); if (!found) n.predicate()->accept(*this); }
                void visit(const SetBuilderNode& n) override { n.domain()->accept(*this); if (!found) n.predicate()->accept(*this); }
                void visit(const ComplexNode& n) override { n.real()->accept(*this); if (!found) n.imag()->accept(*this); }
                void visit(const FiniteSetNode& n) override { for (const auto& e : n.elements()) { if (found) return; e->accept(*this); } }
                void visit(const IntervalNode& n) override { n.lower()->accept(*this); if (!found) n.upper()->accept(*this); }
                void visit(const MembershipNode& n) override { n.element()->accept(*this); if (!found) n.set()->accept(*this); }
                void visit(const QuantityNode& n) override { n.value()->accept(*this); }
            } checker;
            lamina::detail::node(simplified)->accept(checker);
            if (checker.found) return false;
        }

        if (auto num_node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
            if (std::holds_alternative<Rational>(num_node->value())) {
                out_poly.coeffs[i] = std::get<Rational>(num_node->value());
            } else if (std::holds_alternative<BigInt>(num_node->value())) {
                out_poly.coeffs[i] = Rational(std::get<BigInt>(num_node->value()));
            } else {
                out_poly.coeffs[i] = Rational::from_double(std::get<double>(num_node->value()));
            }
        } else if (simplified->is_zero()) {
            out_poly.coeffs[i] = Rational(0);
        } else if (simplified->is_one()) {
            out_poly.coeffs[i] = Rational(1);
        } else {
            ComputationContext context;
            auto evaluated = evaluate_numeric(*simplified, NumericBindings{}, context);
            if (!evaluated || !evaluated.value().is_finite() ||
                !std::isfinite(evaluated.value().value)) {
                return false;
            }
            out_poly.coeffs[i] = Rational::from_double(evaluated.value().value);
        }
    }

    out_poly.trim();
    return true;
}

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

    if (poly.degree() <= 4) {
        return solve_closed_form_poly(poly, var);
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;

    Polynomial<Rational> rat_poly;
    if (!convert_to_rational_poly(poly, rat_poly)) {

        return make_rootof_solutions(poly, var);
    }

    auto rational_roots = find_rational_roots(rat_poly);

    for (const auto& r : rational_roots) {
        results.push_back(SymbolicExpr::number(r));
    }

    Polynomial<Rational> quotient = rat_poly;
    for (const auto& r : rational_roots) {
        Polynomial<Rational> linear_factor({-r, Rational(1)}, var);
        auto [q, rem] = quotient.div_mod(linear_factor);
        quotient = q;
    }

    if (quotient.degree() <= 0) {
        return results;
    }

    if (quotient.degree() <= 4) {
        auto sym_quotient = rational_to_symbolic_poly(quotient);
        auto factor_roots = solve_closed_form_poly(sym_quotient, var);
        results.insert(results.end(), factor_roots.begin(), factor_roots.end());
        return results;
    }

    auto sqfree_factors = square_free_factorization(quotient);

    for (const auto& [factor, multiplicity] : sqfree_factors) {
        if (factor.degree() <= 0) continue;

        std::vector<std::shared_ptr<SymbolicExpr>> factor_solutions;

        if (factor.degree() <= 4) {

            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = solve_closed_form_poly(sym_factor, var);
        } else {

            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = make_rootof_solutions(sym_factor, var);
        }

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

    Polynomial<Rational> current = poly;

    while (current.degree() >= 1 && current.coeffs[0] == Rational(0)) {
        roots.push_back(Rational(0));

        std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
        current = Polynomial<Rational>(new_coeffs, current.variable_name);
    }

    while (current.degree() >= 1) {

        Rational a0 = current.coeffs[0];
        Rational an = current.lead_coeff();

        if (a0 == Rational(0)) {

            roots.push_back(Rational(0));
            std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
            current = Polynomial<Rational>(new_coeffs, current.variable_name);
            continue;
        }

        BigInt lcm_den(1);
        for (const auto& c : current.coeffs) {
            BigInt d = c.get_denominator();
            lcm_den = BigInt::lcm(lcm_den, d);
        }

        BigInt int_a0 = (a0 * Rational(lcm_den)).get_numerator();
        BigInt int_an = (an * Rational(lcm_den)).get_numerator();

        std::vector<BigInt> p_divs = positive_divisors(int_a0);
        std::vector<BigInt> q_divs = positive_divisors(int_an);

        if (p_divs.empty() || q_divs.empty()) {
            break;
        }

        bool found_root = false;
        Rational found_r;

        for (const auto& p : p_divs) {
            for (const auto& q : q_divs) {

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
            break;
        }

        Polynomial<Rational> linear_factor({-found_r, Rational(1)}, current.variable_name);

        while (current.degree() >= 1) {
            auto [quotient, remainder] = current.div_mod(linear_factor);

            if (!remainder.is_zero()) {
                break;
            }
            roots.push_back(found_r);
            current = quotient;

            if (current.degree() < 1 || current.eval(found_r) != Rational(0)) {
                break;
            }
        }
    }

    return roots;
}

}
