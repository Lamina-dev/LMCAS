#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"

namespace LMCAS {

namespace {


// Build sqrt(pi).
inline std::shared_ptr<SymbolicExpr> sf_sqrt_pi() {
    return SymbolicExpr::sqrt(LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<VariableNode>("pi", true)));
}

// Test whether `node` is a single-argument FunctionNode of the given type
// whose argument is exactly the integration variable.
bool sf_is_fn_of_var(const std::shared_ptr<const SymbolicNode>& node,
                     FunctionNode::FuncType t,
                     const std::string& var) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn || fn->type() != t) { return false; }
    if (fn->arguments().size() != 1) { return false; }
    auto v = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
    return v && !v->is_constant() && v->name() == var;
}

// Detect a 1/x factor: a PowerNode whose base is the integration variable and
// exponent is the integer -1.
bool sf_is_inv_var(const std::shared_ptr<const SymbolicNode>& node,
                   const std::string& var) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) { return false; }
    auto b = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (!b || b->is_constant() || b->name() != var) { return false; }
    const auto exponent = detail::exact_rational_value(pw->exponent());
    return exponent && *exponent == Rational(-1);
}

// Detect a 1/ln(var) factor, i.e. PowerNode(ln(var), -1).
bool sf_is_inv_ln_var(const std::shared_ptr<const SymbolicNode>& node,
                      const std::string& var) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) { return false; }
    if (!sf_is_fn_of_var(pw->base(), FunctionNode::FuncType::Ln, var)) { return false; }
    const auto exponent = detail::exact_rational_value(pw->exponent());
    return exponent && *exponent == Rational(-1);
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
            if (has_inv_var) { return false; } /** @brief 多个 1/x 因子形成 1/x^2 类幂次，超出此路径的处理范围。 */
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
Result<bool> sf_match_exp_neg_quad(const std::shared_ptr<const SymbolicNode>& node,
                           const std::string& var,
                           Rational& c_out, ComputationContext& context) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn || fn->type() != FunctionNode::FuncType::Exp) { return false; }
    if (fn->arguments().size() != 1) { return false; }
    auto arg = LMCAS::detail::expression_from_node(fn->arguments()[0]);
    auto converted = recognize_rational_polynomial(arg, var, context);
    if (!converted) { return Result<bool>::failure(converted.error()); }
    if (!converted.value()) { return false; }
    const auto& poly = *converted.value();
    if (poly.degree() != 2) { return false; }
    if (poly.coeffs.size() < 3) { return false; }
    // Constant and linear terms must be zero.
    if (!(poly.coeffs[0] == Rational(0))) { return false; }
    if (!(poly.coeffs[1] == Rational(0))) { return false; }
    Rational a2 = poly.coeffs[2];
    // Need a strictly negative coefficient: exp(-c*x^2) with c > 0.
    if (!(a2 < Rational(0))) { return false; }
    c_out = Rational(0) - a2; // c = -a2 > 0
    return true;
}

} // anonymous namespace

Result<std::optional<detail::IntegralParts>> Integrator::try_special_function(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int) try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");

    using FT = FunctionNode::FuncType;

    auto v = SymbolicExpr::variable(var);

    if (sf_is_inv_ln_var(LMCAS::detail::node(expr), var)) {
        auto li = make_unary_function(FT::Li, v);
        auto simp = detail::make_expression_ptr(detail::normalize_node(detail::node(li), transform_budget));
        return std::optional<detail::IntegralParts>{{*(simp ? simp : li), zero_}};
    }

    {
        Rational c_rat;
        auto matched = sf_match_exp_neg_quad(LMCAS::detail::node(expr), var, c_rat, context);
        if (!matched) {
            return Result<std::optional<detail::IntegralParts>>::failure(matched.error());
        }
        if (matched.value()) {
            // Result: (sqrt(pi) / (2 * sqrt(c))) * erf(sqrt(c) * x)
            auto sqrt_pi = sf_sqrt_pi();
            auto sqrt_c = SymbolicExpr::sqrt(SymbolicExpr::number(c_rat));
            auto scaled_x = SymbolicExpr::multiply(sqrt_c, v);
            auto erf_node = make_unary_function(FT::Erf, scaled_x);

            auto two = SymbolicExpr::number(2);
            auto two_sqrt_c = SymbolicExpr::multiply(two, sqrt_c);
            auto coeff = SymbolicExpr::multiply(
                sqrt_pi,
                SymbolicExpr::power(two_sqrt_c, SymbolicExpr::number(-1)));
            auto result = SymbolicExpr::multiply(coeff, erf_node);
            auto simp = detail::make_expression_ptr(detail::normalize_node(detail::node(result), transform_budget));
            return std::optional<detail::IntegralParts>{{*(simp ? simp : result), zero_}};
        }
    }

    // These all have shape (something)*1/x, where the "something" is a
    // FunctionNode of var. Other 1/x patterns (e.g. 1/x alone, x*1/x) are
    // not our concern.
    {
        std::vector<std::shared_ptr<const SymbolicNode>> others;
        bool has_inv = false;
        if (sf_split_inv_var(LMCAS::detail::node(expr), var, others, has_inv) && has_inv && others.size() == 1) {
            const auto& other = others[0];
            if (sf_is_fn_of_var(other, FT::Exp, var)) {
                auto ei = make_unary_function(FT::Ei, v);
                auto simp = detail::make_expression_ptr(detail::normalize_node(detail::node(ei), transform_budget));
                return std::optional<detail::IntegralParts>{{*(simp ? simp : ei), zero_}};
            }
            if (sf_is_fn_of_var(other, FT::Sin, var)) {
                auto si = make_unary_function(FT::Si, v);
                auto simp = detail::make_expression_ptr(detail::normalize_node(detail::node(si), transform_budget));
                return std::optional<detail::IntegralParts>{{*(simp ? simp : si), zero_}};
            }
            if (sf_is_fn_of_var(other, FT::Cos, var)) {
                auto ci = make_unary_function(FT::Ci, v);
                auto simp = detail::make_expression_ptr(detail::normalize_node(detail::node(ci), transform_budget));
                return std::optional<detail::IntegralParts>{{*(simp ? simp : ci), zero_}};
            }
        }
    }

    return std::optional<detail::IntegralParts>{};
} catch (const CasError& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(error);
}


} // namespace LMCAS
