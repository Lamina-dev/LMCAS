#include "internal/visitors/limit_visitor.hpp"
#include "polynomial_conversion.hpp"
#include <algorithm>

namespace LMCAS {

bool LimitVisitor::can_rationalize_sum(const AddNode& node,
    const std::vector<std::shared_ptr<const SymbolicNode>>& values) const {
    if (!is_limit_at_infinity() || node.operands().size() != 2 || values.size() != 2) {
        return false;
    }
    if (!is_inf(values[0]) || !is_inf(values[1]) ||
        is_neg_inf(values[0]) == is_neg_inf(values[1])) {
        return false;
    }
    return defined_near_point(node.operands()[0]) && defined_near_point(node.operands()[1]);
}

std::optional<Polynomial<Rational>> LimitVisitor::square_algebraic_term(
    const std::shared_ptr<const SymbolicNode>& term) const {
    charge();
    std::shared_ptr<const SymbolicNode> radicand;
    if (auto root = std::dynamic_pointer_cast<const PowerNode>(term)) {
        auto exponent = extract_coeff_value<Rational>(detail::make_expression_ptr(root->exponent()));
        if (exponent && exponent.value() == Rational(BigInt(1), BigInt(2)))
            { radicand = root->base(); }
    } else if (auto root = std::dynamic_pointer_cast<const FunctionNode>(term)) {
        if (root->type() == FunctionNode::FuncType::Sqrt && root->arguments().size() == 1)
            { radicand = root->arguments()[0]; }
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(term)) {
        Polynomial<Rational> result(std::vector<Rational>{Rational(1)}, var);
        for (const auto& factor : product->operands()) {
            auto factor_square = square_algebraic_term(factor);
            if (!factor_square) { return std::nullopt; }
            result = result * *factor_square;
        }
        return result;
    }
    auto polynomial = symbolic_to_poly<Rational>(
        detail::make_expression_ptr(radicand ? radicand : term), var);
    if (!polynomial) {
        if (polynomial.error().code != CasErrc::UnsupportedExpression) { throw polynomial.error(); }
        return std::nullopt;
    }
    return radicand ? polynomial.value() : polynomial.value() * polynomial.value();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::rationalized_sum_limit(
    const AddNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    if (!can_rationalize_sum(node, values)) { return nullptr; }
    /**
     * @brief 用 (a+b)(a-b)=a²-b² 有理化两个代数项。
     * 两项趋向异号无穷，保证共轭式首项系数非零；截断后全零的级数不足以证明该条件。
     */
    auto a = square_algebraic_term(node.operands()[0]), b = square_algebraic_term(node.operands()[1]);
    if (!a || !b) { return nullptr; }
    if (a->is_zero() || b->is_zero()) { return nullptr; }
    if (a->lead_coeff() <= Rational(0) || b->lead_coeff() <= Rational(0)) { return nullptr; }
    auto numerator = *a - *b;
    if (numerator.is_zero()) { return detail::make_node<NumberNode>(BigInt(0)); }
    const int degree = std::max(a->degree(), b->degree());
    const int denominator_sign = is_neg_inf(values[0]) ? -1 : 1;
    if (2 * numerator.degree() < degree) { return detail::make_node<NumberNode>(BigInt(0)); }
    if (2 * numerator.degree() > degree) {
        return signed_infinity((numerator.lead_coeff() > Rational(0) ? 1 : -1) * denominator_sign);
    }
    auto leading = [&](const Polynomial<Rational>& polynomial) {
        return polynomial.degree() == degree ? SymbolicExpr::sqrt(SymbolicExpr::number(polynomial.lead_coeff()))
            : SymbolicExpr::number(0);
    };
    auto denominator = SymbolicExpr::multiply(SymbolicExpr::number(denominator_sign),
        SymbolicExpr::add(leading(*a), leading(*b)));
    return detail::node(SymbolicExpr::divide(SymbolicExpr::number(numerator.lead_coeff()), denominator)->simplify());
}

std::pair<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const SymbolicNode>>
LimitVisitor::combine_sum_ratio(
    const std::vector<std::shared_ptr<const SymbolicNode>>& numerators,
    const std::vector<std::shared_ptr<const SymbolicNode>>& denominators) {
    auto common_denominator = make_product_or_one(denominators);
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    for (std::size_t i = 0; i < numerators.size(); ++i) {
        std::vector<std::shared_ptr<const SymbolicNode>> factors{numerators[i]};
        for (std::size_t j = 0; j < denominators.size(); ++j) {
            if (j != i) { factors.push_back(denominators[j]); }
        }
        terms.push_back(make_product_or_one(factors));
    }
    auto numerator = terms.size() == 1 ? terms[0] :
        std::static_pointer_cast<const SymbolicNode>(LMCAS::detail::make_node<AddNode>(terms));
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    numerator->accept(norm);
    numerator = norm.get_result();
    common_denominator->accept(norm);
    common_denominator = norm.get_result();
    return {std::move(numerator), std::move(common_denominator)};
}



std::shared_ptr<const SymbolicNode> LimitVisitor::denominator_factor(
    const std::shared_ptr<const SymbolicNode>& node) const {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return nullptr;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    if (!exponent) {
        return nullptr;
    }
    auto value = extract_coeff_value<Rational>(detail::make_expression_ptr(power->exponent()));
    if (!value || value.value() >= Rational(0)) return nullptr;
    if (value.value() == Rational(-1)) return power->base();
    return LMCAS::detail::make_node<PowerNode>(
        power->base(), LMCAS::detail::make_node<NumberNode>(-value.value()));
}

std::pair<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const SymbolicNode>>
LimitVisitor::extract_num_den(const std::shared_ptr<const SymbolicNode>& expr) {
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        std::vector<std::shared_ptr<const SymbolicNode>> numerator, denominator;
        for (const auto& operand : product->operands()) {
            auto factor = denominator_factor(operand);
            if (factor) denominator.push_back(factor);
            else numerator.push_back(operand);
        }
        if (!denominator.empty()) {
            return {make_product_or_one(numerator), make_product_or_one(denominator)};
        }
    }
    auto denominator = denominator_factor(expr);
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    if (denominator) {
        return {one, denominator};
    }
    return {expr, one};
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_zero_infinity_pair(
    const std::shared_ptr<const SymbolicNode>& zero,
    const std::shared_ptr<const SymbolicNode>& infinity) {
    bool prefer_infinity = false;
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        auto zero_growth = classify_growth(zero);
        auto infinity_growth = classify_growth(infinity);
        prefer_infinity = zero_growth == GrowthClass::Exponential &&
                          infinity_growth == GrowthClass::Polynomial;
    }
    const auto& first = prefer_infinity ? infinity : zero;
    const auto& second = prefer_infinity ? zero : infinity;
    auto inverse = LMCAS::detail::make_node<PowerNode>(
        second, LMCAS::detail::make_node<NumberNode>(BigInt(-1)));
    auto value = apply_lhopital(first, inverse);
    if (value) {
        return value;
    }
    inverse = LMCAS::detail::make_node<PowerNode>(
        first, LMCAS::detail::make_node<NumberNode>(BigInt(-1)));
    return apply_lhopital(second, inverse);
}

std::shared_ptr<const SymbolicNode> LimitVisitor::multiply_remaining_limits(
    const std::shared_ptr<const SymbolicNode>& value,
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors) {
    if (!value || factors.empty()) {
        return value;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> values{value};
    for (const auto& factor : factors) {
        auto limit = eval_limit(factor);
        if (!limit) return nullptr;
        values.push_back(limit);
    }
    if (values.size() == 1) {
        return values[0];
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    LMCAS::detail::make_node<MultiplyNode>(values)->accept(norm);
    return norm.get_result();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_zero_times_inf(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors,
    const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    std::vector<std::shared_ptr<const SymbolicNode>> zero, infinity, other;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (is_zero_limit(values[i])) zero.push_back(factors[i]);
        else if (is_inf(values[i])) infinity.push_back(factors[i]);
        else other.push_back(factors[i]);
    }
    if (zero.empty() || infinity.empty()) {
        return nullptr;
    }
    auto first = make_product_or_one(zero);
    auto second = make_product_or_one(infinity);
    return multiply_remaining_limits(resolve_zero_infinity_pair(first, second), other);
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_combined_ratio(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    auto numerator = eval_limit(num);
    auto denominator = eval_limit(den);
    bool numerator_zero = is_zero_limit(numerator);
    bool denominator_zero = is_zero_limit(denominator);
    if ((numerator_zero && denominator_zero) || (is_inf(numerator) && is_inf(denominator))) {
        return apply_lhopital(num, den);
    }
    if (denominator_zero || is_inf(denominator) || !denominator) {
        return nullptr;
    }
    auto ratio = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            num, LMCAS::detail::make_node<PowerNode>(
                den, LMCAS::detail::make_node<NumberNode>(BigInt(-1)))});
    return eval_limit(ratio);
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_inf_minus_inf(
    const AddNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    auto rationalized = rationalized_sum_limit(node, values);
    if (rationalized) { return rationalized; }
    std::vector<std::shared_ptr<const SymbolicNode>> numerators, denominators;
    bool has_denominator = false;
    for (const auto& operand : node.operands()) {
        auto fraction = extract_num_den(operand);
        numerators.push_back(fraction.first);
        denominators.push_back(fraction.second);
        if (!fraction.second->is_one()) { has_denominator = true; }
    }
    if (!has_denominator) {
        auto whole = node.operands().size() == 1 ? node.operands()[0] :
            std::static_pointer_cast<const SymbolicNode>(LMCAS::detail::make_node<AddNode>(node.operands()));
        return taylor_fallback(whole, LMCAS::detail::make_node<NumberNode>(BigInt(1)));
    }
    auto combined = combine_sum_ratio(numerators, denominators);
    return resolve_combined_ratio(combined.first, combined.second);
}

void LimitVisitor::visit(const AddNode& node) {
    if (is_limit_at_neg_infinity()) {
        auto value = handle_neg_infinity_limit(LMCAS::detail::make_node<AddNode>(node.operands()));
        if (value) {
            result = value;
            return;
        }
    }
    if (is_limit_at_infinity()) {
        auto polynomial = LMCAS::detail::make_node<AddNode>(node.operands());
        if (get_polynomial_degree(polynomial) > 0) {
            const int sign = get_sign(get_leading_coefficient(polynomial));
            if (sign != 0) {
                result = signed_infinity(sign);
                return;
            }
        }
    }
    auto squeezed = try_squeeze(LMCAS::detail::make_node<AddNode>(node.operands()));
    if (squeezed) {
        result = squeezed;
        return;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> values;
    for (const auto& operand : node.operands()) {
        operand->accept(*this);
        if (!result) return;
        values.push_back(result);
    }
    if (classify_add_form(values) == IndeterminateForm::InfMinusInf) {
        auto resolved = resolve_inf_minus_inf(node, values);
        if (resolved) {
            result = resolved;
            return;
        }
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    LMCAS::detail::make_node<AddNode>(values)->accept(norm);
    result = norm.get_result();
}

}
