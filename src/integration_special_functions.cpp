#include "internal/integration_support.hpp"

namespace lamina {

namespace {

// Build a single-argument FunctionNode wrapped in SymbolicExpr.
inline std::shared_ptr<SymbolicExpr> sf_make_fn(
    FunctionNode::FuncType t,
    const std::shared_ptr<SymbolicExpr>& arg) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            t, std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(arg)}));
}

// Build sqrt(pi).
inline std::shared_ptr<SymbolicExpr> sf_sqrt_pi() {
    return SymbolicExpr::sqrt(SymbolicExpr::variable("pi"));
}

// Test whether `node` is a single-argument FunctionNode of the given type
// whose argument is exactly the integration variable.
bool sf_is_fn_of_var(const std::shared_ptr<const SymbolicNode>& node,
                     FunctionNode::FuncType t,
                     const std::string& var) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn || fn->type() != t) return false;
    if (fn->arguments().size() != 1) return false;
    auto v = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
    return v && v->name() == var;
}

// Detect a 1/x factor: a PowerNode whose base is the integration variable and
// exponent is the integer -1.
bool sf_is_inv_var(const std::shared_ptr<const SymbolicNode>& node,
                   const std::string& var) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) return false;
    auto b = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (!b || b->name() != var) return false;
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) return false;
    if (std::holds_alternative<BigInt>(en->value())) {
        return std::get<BigInt>(en->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(en->value())) {
        return std::get<Rational>(en->value()) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(en->value())) {
        lmmc_real_t d = std::get<lmmc_real_t>(en->value());
        int eq = 0;
        lmmc_double_nearly_equal_tol(d, -1.0, 1e-12, 1e-12, &eq);
        return eq != 0;
    }
    return false;
}

// Detect a 1/ln(var) factor, i.e. PowerNode(ln(var), -1).
bool sf_is_inv_ln_var(const std::shared_ptr<const SymbolicNode>& node,
                      const std::string& var) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) return false;
    if (!sf_is_fn_of_var(pw->base(), FunctionNode::FuncType::Ln, var)) return false;
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) return false;
    if (std::holds_alternative<BigInt>(en->value())) {
        return std::get<BigInt>(en->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(en->value())) {
        return std::get<Rational>(en->value()) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(en->value())) {
        lmmc_real_t d = std::get<lmmc_real_t>(en->value());
        int eq = 0;
        lmmc_double_nearly_equal_tol(d, -1.0, 1e-12, 1e-12, &eq);
        return eq != 0;
    }
    return false;
}

// Split a node into (factors, has_inv_var) where the inv-var factor is removed
// from `factors` if present. Returns false if there is more than one inv-var
// factor (which would be 1/x^2 and is not the form we handle here).
bool sf_split_inv_var(const std::shared_ptr<const SymbolicNode>& node,
                      const std::string& var,
                      std::vector<std::shared_ptr<const SymbolicNode>>& other_factors,
                      bool& has_inv_var) {
    other_factors.clear();
    has_inv_var = false;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        factors = mul->operands();
    } else {
        factors.push_back(node);
    }
    for (const auto& f : factors) {
        if (sf_is_inv_var(f, var)) {
            if (has_inv_var) return false; // multiple 1/x factors (1/x^2 form)
            has_inv_var = true;
        } else {
            other_factors.push_back(f);
        }
    }
    return true;
}

// Detect exp(-c*var^2) where c is a constant w.r.t. var that is rational and
// strictly positive. On success, sets `c_out` to that rational coefficient.
//
// We accept exp(arg) where the argument is a polynomial in var of degree
// exactly 2, no constant or linear term, and the leading rational coefficient
// is strictly negative. (We then return its absolute value as c.)
bool sf_match_exp_neg_quad(const std::shared_ptr<const SymbolicNode>& node,
                           const std::string& var,
                           Rational& c_out) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn || fn->type() != FunctionNode::FuncType::Exp) return false;
    if (fn->arguments().size() != 1) return false;
    auto arg = lamina::detail::expression_from_node(fn->arguments()[0]);
    Polynomial<Rational> poly;
    try {
        poly = symbolic_to_poly<Rational>(lamina::detail::make_expression_ptr(arg), var);
    } catch (...) {
        return false;
    }
    if (poly.degree() != 2) return false;
    if (poly.coeffs.size() < 3) return false;
    // Constant and linear terms must be zero.
    if (!(poly.coeffs[0] == Rational(0))) return false;
    if (!(poly.coeffs[1] == Rational(0))) return false;
    Rational a2 = poly.coeffs[2];
    // Need a strictly negative coefficient: exp(-c*x^2) with c > 0.
    if (!(a2 < Rational(0))) return false;
    c_out = Rational(0) - a2; // c = -a2 > 0
    return true;
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> SpecialFunctionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext&, int depth) {
    (void)ctx;
    (void)depth;

    using FT = FunctionNode::FuncType;

    auto v = SymbolicExpr::variable(var);

    if (sf_is_inv_ln_var(lamina::detail::node(expr), var)) {
        auto li = sf_make_fn(FT::Li, v);
        auto simp = li->simplify();
        return simp ? simp : li;
    }

    {
        Rational c_rat;
        if (sf_match_exp_neg_quad(lamina::detail::node(expr), var, c_rat)) {
            // Result: (sqrt(pi) / (2 * sqrt(c))) * erf(sqrt(c) * x)
            auto sqrt_pi = sf_sqrt_pi();
            auto sqrt_c = SymbolicExpr::sqrt(SymbolicExpr::number(c_rat));
            auto scaled_x = SymbolicExpr::multiply(sqrt_c, v);
            auto erf_node = sf_make_fn(FT::Erf, scaled_x);

            auto two = SymbolicExpr::number(2);
            auto two_sqrt_c = SymbolicExpr::multiply(two, sqrt_c);
            auto coeff = SymbolicExpr::multiply(
                sqrt_pi,
                SymbolicExpr::power(two_sqrt_c, SymbolicExpr::number(-1)));
            auto result = SymbolicExpr::multiply(coeff, erf_node);
            auto simp = result->simplify();
            return simp ? simp : result;
        }
    }

    // These all have shape (something)*1/x, where the "something" is a
    // FunctionNode of var. Other 1/x patterns (e.g. 1/x alone, x*1/x) are
    // not our concern.
    {
        std::vector<std::shared_ptr<const SymbolicNode>> others;
        bool has_inv = false;
        if (sf_split_inv_var(lamina::detail::node(expr), var, others, has_inv) && has_inv && others.size() == 1) {
            const auto& other = others[0];
            if (sf_is_fn_of_var(other, FT::Exp, var)) {
                auto ei = sf_make_fn(FT::Ei, v);
                auto simp = ei->simplify();
                return simp ? simp : ei;
            }
            if (sf_is_fn_of_var(other, FT::Sin, var)) {
                auto si = sf_make_fn(FT::Si, v);
                auto simp = si->simplify();
                return simp ? simp : si;
            }
            if (sf_is_fn_of_var(other, FT::Cos, var)) {
                auto ci = sf_make_fn(FT::Ci, v);
                auto simp = ci->simplify();
                return simp ? simp : ci;
            }
        }
    }

    return nullptr;
}


} // namespace lamina
