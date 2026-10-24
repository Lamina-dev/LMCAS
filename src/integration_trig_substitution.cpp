#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "polynomial_conversion.hpp"

namespace LMCAS {

namespace {

// 标准二次根式使用精确有理系数，保留原被积式的指数和尺度。
struct QuadRadical {
    int pattern = 0;       // 1: a²-x², 2: a²+x², 3: x²-a²
    Rational exponent{0};
    Rational a_sq{0};
};

bool trigsub_half_exponent(const NumberNode& number, Rational& exponent) {
    const auto exact = detail::exact_rational_value(number);
    if (!exact || (*exact != Rational(1, 2) && *exact != Rational(-1, 2))) {
        return false;
    }
    exponent = *exact;
    return true;
}

bool trigsub_quadratic_coefficients(const detail::SymbolicNodePtr& radicand,
    const std::string& var, ComputationContext& context, Rational& c0, Rational& c2) {
    auto polynomial = recognize_rational_polynomial(
        detail::expression_from_node(radicand), var, context);
    if (!polynomial) { throw polynomial.error(); }
    if (!polynomial.value() || polynomial.value()->degree() != 2) { return false; }
    const auto& coefficients = polynomial.value()->coeffs;
    if (coefficients[1] != Rational(0) ||
        (coefficients[2] != Rational(1) && coefficients[2] != Rational(-1))) {
        return false;
    }
    c0 = coefficients[0];
    c2 = coefficients[2];
    return true;
}

bool trigsub_match_radical(const std::shared_ptr<const SymbolicNode>& node,
                           const std::string& var,
                           ComputationContext& context,
                           QuadRadical& out) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) { return false; }
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) { return false; }
    Rational exponent;
    auto radicand = pw->base();
    const auto outer = detail::exact_rational_value(*en);
    if (outer && *outer == Rational(-1)) {
        if (const auto* root = dynamic_cast<const PowerNode*>(radicand.get())) {
            const auto inner = detail::exact_rational_value(root->exponent());
            if (!inner || *inner != Rational(1, 2)) { return false; }
            radicand = root->base();
        } else if (const auto* root = dynamic_cast<const FunctionNode*>(radicand.get());
                   root && root->type() == FunctionNode::FuncType::Sqrt &&
                   root->arguments().size() == 1) {
            radicand = root->arguments()[0];
        } else {
            return false;
        }
        exponent = Rational(-1, 2);
    } else if (!trigsub_half_exponent(*en, exponent)) {
        return false;
    }
    Rational c0, c2;
    if (!trigsub_quadratic_coefficients(radicand, var, context, c0, c2)) { return false; }
    out.exponent = exponent;
    if (c2 < Rational(0) && c0 > Rational(0)) { out.pattern = 1; out.a_sq = c0; }
    else if (c2 > Rational(0) && c0 > Rational(0)) { out.pattern = 2; out.a_sq = c0; }
    else if (c2 > Rational(0) && c0 < Rational(0)) { out.pattern = 3; out.a_sq = Rational(0) - c0; }
    else return false;
    return true;
}

} // anonymous namespace

Result<std::optional<detail::IntegralParts>> Integrator::try_trig_substitution(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!LMCAS::detail::node(expr)) { return std::optional<detail::IntegralParts>{}; }

    auto x = SymbolicExpr::variable(var);

    /// 仅处理被积函数恰为单个二次根式幂的常见闭式情形：
    ///   ∫ (a²-x²)^(-1/2) dx = arcsin(x/a)
    ///   ∫ (a²+x²)^(-1/2) dx = arcsinh(x/a) = ln(x + √(x²+a²))
    ///   ∫ (x²-a²)^(-1/2) dx = arccosh(x/a) = ln(x + √(x²-a²))
    ///   ∫ (a²-x²)^( 1/2) dx = (x/2)√(a²-x²) + (a²/2)arcsin(x/a)
    QuadRadical qr;
    if (!trigsub_match_radical(
            LMCAS::detail::node(expr), var, context, qr)) {
        return std::optional<detail::IntegralParts>{};
    }

    auto a_sq_expr = SymbolicExpr::number(qr.a_sq);
    auto a_expr = SymbolicExpr::sqrt(a_sq_expr); // a = √(a²)

    /// 构造 √(模式) 表达式
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto make_radicand = [&]() -> std::shared_ptr<SymbolicExpr> {
        if (qr.pattern == 1) // a²-x²
            return SymbolicExpr::add(a_sq_expr, SymbolicExpr::multiply(SymbolicExpr::number(-1), x_sq));
        if (qr.pattern == 2) // a²+x²
            return SymbolicExpr::add(a_sq_expr, x_sq);
        return SymbolicExpr::add(x_sq, SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)); // x²-a²
    };
    auto x_over_a = SymbolicExpr::divide(x, a_expr);

    auto arcsin = [&](const std::shared_ptr<SymbolicExpr>& u) {
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcSin, std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(u)}));
    };
    auto ln = [&](const std::shared_ptr<SymbolicExpr>& u) { return SymbolicExpr::ln(u); };

    if (qr.exponent < Rational(0)) {
        /// (...)^(-1/2)
        if (qr.pattern == 1) {
            /// arcsin(x/a)
            return std::optional<detail::IntegralParts>{{detail::expression_from_node(
                detail::normalize_node(detail::node(arcsin(x_over_a)), transform_budget)), zero_}};
        }
        if (qr.pattern == 2) {
            /// ln(x + √(x²+a²))
            auto rad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq, a_sq_expr));
            return std::optional<detail::IntegralParts>{{detail::expression_from_node(
                detail::normalize_node(detail::node(ln(SymbolicExpr::add(x, rad))), transform_budget)), zero_}};
        }
        /// pattern 3: ln(x + √(x²-a²))
        auto rad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)));
        return std::optional<detail::IntegralParts>{{detail::expression_from_node(
            detail::normalize_node(detail::node(ln(SymbolicExpr::add(x, rad))), transform_budget)), zero_}};
    } else {
        /// (...)^(+1/2)
        auto rad = SymbolicExpr::sqrt(make_radicand());
        auto half = SymbolicExpr::number(Rational(1, 2));
        if (qr.pattern == 1) {
            /// (x/2)√(a²-x²) + (a²/2)arcsin(x/a)
            auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
            auto term2 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, a_sq_expr), arcsin(x_over_a));
            return std::optional<detail::IntegralParts>{{detail::expression_from_node(
                detail::normalize_node(detail::node(SymbolicExpr::add(term1, term2)), transform_budget)), zero_}};
        }
        if (qr.pattern == 2) {
            /// (x/2)√(a²+x²) + (a²/2)ln(x+√(x²+a²))
            auto lrad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq, a_sq_expr));
            auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
            auto term2 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, a_sq_expr),
                ln(SymbolicExpr::add(x, lrad)));
            return std::optional<detail::IntegralParts>{{detail::expression_from_node(
                detail::normalize_node(detail::node(SymbolicExpr::add(term1, term2)), transform_budget)), zero_}};
        }
        /// pattern 3: (x/2)√(x²-a²) - (a²/2)ln(x+√(x²-a²))
        auto lrad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)));
        auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
        auto term2 = SymbolicExpr::multiply(
            SymbolicExpr::multiply(SymbolicExpr::multiply(SymbolicExpr::number(-1), half), a_sq_expr),
            ln(SymbolicExpr::add(x, lrad)));
        return std::optional<detail::IntegralParts>{{detail::expression_from_node(
            detail::normalize_node(detail::node(SymbolicExpr::add(term1, term2)), transform_budget)), zero_}};
    }
} catch (const CasError& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(error);
}

}
