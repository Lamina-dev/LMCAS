#include "internal/visitors/limit_visitor.hpp"
#include "internal/facts_query.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/normalization_utils.hpp"
#include "polynomial_conversion.hpp"

namespace LMCAS {

std::shared_ptr<const SymbolicNode> LimitVisitor::infinite_base_power_limit(
    const std::shared_ptr<const SymbolicNode>& base,
    const std::shared_ptr<const NumberNode>& number, bool integral, const BigInt& integer) {
    if (is_neg_inf(base) && !integral) { return nullptr; }
    if (get_sign(number) < 0) { return detail::make_node<NumberNode>(BigInt(0)); }
    if (number->is_positive()) { return signed_infinity(is_neg_inf(base) && integer.is_odd() ? -1 : 1); }
    return nullptr;
}

namespace {
std::shared_ptr<const SymbolicNode> collect_log_power_factors(
    const std::shared_ptr<const SymbolicNode>& exponent,
    std::vector<std::shared_ptr<const SymbolicNode>>& coefficients) {
    std::shared_ptr<const SymbolicNode> logarithm;
    auto collect = [&](const std::shared_ptr<const SymbolicNode>& op) {
        auto f = std::dynamic_pointer_cast<const FunctionNode>(op);
        if (!logarithm && f && f->type() == FunctionNode::FuncType::Ln && f->arguments().size() == 1)
            { logarithm = op; }
        else { coefficients.push_back(op); }
    };
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(exponent))
        { for (const auto& op : product->operands()) { collect(op); } }
    else { collect(exponent); }
    return logarithm;
}
}



std::shared_ptr<const SymbolicNode> LimitVisitor::linear_exponent_coefficient(
    const std::shared_ptr<const SymbolicNode>& node) const {
    if (is_limit_variable(node)) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return nullptr;
    }
    bool has_variable = false;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    for (const auto& operand : product->operands()) {
        if (!has_variable && is_limit_variable(operand)) has_variable = true;
        else factors.push_back(operand);
    }
    if (!has_variable) {
        return nullptr;
    }
    return make_product_or_one(factors);
}

bool LimitVisitor::exponential_power_limit(const PowerNode& node) {
    auto sum = std::dynamic_pointer_cast<const AddNode>(node.base());
    if (!sum) {
        return false;
    }
    auto coefficient = linear_exponent_coefficient(node.exponent());
    if (!coefficient) {
        return false;
    }
    auto term = one_plus_term(*sum);
    if (!term) {
        return false;
    }
    auto numerator = variable_fraction_numerator(term);
    if (!numerator) {
        return false;
    }
    if (!is_zero_limit(eval_limit(term))) return false;
    auto inner = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{numerator, coefficient});
    auto value = eval_limit(inner);
    if (!value) return false;
    auto exponential = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Exp,
        std::vector<std::shared_ptr<const SymbolicNode>>{value});
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    exponential->accept(norm);
    result = norm.get_result();
    return true;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_exponential_form(
    const std::shared_ptr<const SymbolicNode>& base,
    const std::shared_ptr<const SymbolicNode>& exponent) {
    auto logarithm = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<const SymbolicNode>>{base});
    auto product = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{exponent, logarithm});
    auto value = eval_limit(product);
    if (!value) {
        return nullptr;
    }
    if (!is_inf(value)) {
        auto exponential = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp, std::vector<std::shared_ptr<const SymbolicNode>>{value});
        NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
        exponential->accept(norm);
        return norm.get_result();
    }
    if (!is_neg_inf(value)) {
        return signed_infinity(1);
    }
    return LMCAS::detail::make_node<NumberNode>(BigInt(0));
}

std::shared_ptr<const SymbolicNode> LimitVisitor::singular_power_limit(
    const PowerNode& node, const std::shared_ptr<const SymbolicNode>& base,
    const std::shared_ptr<const SymbolicNode>& exponent) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(exponent);
    if (!number) {
        return nullptr;
    }
    BigInt integer;
    const bool integral = try_get_integer_value(number, integer);
    double value = get_numeric_value(number);
    if (is_inf(base)) {
        return infinite_base_power_limit(base, number, integral, integer);
    }
    if (!is_zero_limit(base) || !(value < 0) || direction.empty()) {
        return nullptr;
    }
    auto sign = determine_sign_near_point(node.base(), direction);
    if (!sign || *sign == 0) { return nullptr; }
    if (!integral && *sign < 0) { return nullptr; }
    return signed_infinity(integral && integer.is_odd() ? *sign : 1);
}

void LimitVisitor::visit(const PowerNode& node) {
    if (is_limit_at_neg_infinity()) {
        auto value = handle_neg_infinity_limit(
            LMCAS::detail::make_node<PowerNode>(node.base(), node.exponent()));
        if (value) {
            result = value;
            return;
        }
    }
    if (is_limit_at_infinity() && exponential_power_limit(node)) {
        return;
    }
    if (node.base()->is_one() && defined_near_point(node.exponent())) {
        result = detail::make_node<NumberNode>(BigInt(1));
        return;
    }
    if (node.exponent()->is_zero()) {
        auto sign = determine_sign_near_point(node.base(), direction);
        if (sign && *sign != 0) {
            result = detail::make_node<NumberNode>(BigInt(1));
            return;
        }
    }
    node.base()->accept(*this);
    auto base = result;
    node.exponent()->accept(*this);
    auto exponent = result;
    if (!base || !exponent) { result = nullptr; return; }
    if (classify_power_form(base, exponent) != IndeterminateForm::None) {
        auto value = resolve_exponential_form(node.base(), node.exponent());
        if (value) {
            result = value;
            return;
        }
    }
    auto singular = singular_power_limit(node, base, exponent);
    if (singular) {
        result = singular;
        return;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    LMCAS::detail::make_node<PowerNode>(base, exponent)->accept(norm);
    result = norm.get_result();
}

bool LimitVisitor::infinite_function_limit(FunctionNode::FuncType type, bool negative) {
    switch (type) {
        case FunctionNode::FuncType::Exp:
            result = negative ? std::static_pointer_cast<const SymbolicNode>(
                LMCAS::detail::make_node<NumberNode>(BigInt(0))) : signed_infinity(1);
            return true;
        case FunctionNode::FuncType::Ln:
            result = negative ? nullptr : signed_infinity(1);
            return !negative;
        case FunctionNode::FuncType::ArcTan:
            result = LMCAS::detail::make_node<NumberNode>(Rational(negative ? -1 : 1, 2));
            result = LMCAS::detail::make_node<MultiplyNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    result, LMCAS::detail::make_node<VariableNode>("pi")});
            return true;
        default:
            return false;
    }
}

bool LimitVisitor::composed_function_limit(const FunctionNode& node) {
    if (node.arguments().size() != 1 ||
        node.type() == FunctionNode::FuncType::Infinity) {
        return false;
    }
    auto inner = eval_limit(node.arguments()[0]);
    if (!inner) {
        return false;
    }
    if (!is_inf(inner)) {
        auto candidate = LMCAS::detail::make_node<FunctionNode>(node.type(),
            std::vector<std::shared_ptr<const SymbolicNode>>{inner});
        std::optional<detail::AssumptionFacts> assumed;
        if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
        const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
        auto defined = detail::query_definedness(candidate, facts, Domain::Real, *context_);
        if (!defined) { throw defined.error(); }
        if (defined && defined.value() == Tribool::True) {
            NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
            candidate->accept(norm);
            result = norm.get_result();
            return true;
        }
        if (node.type() == FunctionNode::FuncType::Ln && is_zero_limit(inner)) {
            auto sign = determine_sign_near_point(node.arguments()[0], direction);
            if (sign && *sign > 0) {
                result = signed_infinity(-1);
                return true;
            }
        }
    }
    if (is_inf(inner)) {
        return infinite_function_limit(node.type(), is_neg_inf(inner));
    }
    return false;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::positive_tail_log_power(const FunctionNode& node) const {
    if (!is_limit_at_infinity() || node.type() != FunctionNode::FuncType::Exp ||
        node.arguments().size() != 1) { return nullptr; }
    std::vector<std::shared_ptr<const SymbolicNode>> coefficients;
    auto logarithm = collect_log_power_factors(node.arguments()[0], coefficients);
    if (!logarithm) { return nullptr; }
    auto f = std::dynamic_pointer_cast<const FunctionNode>(logarithm);
    std::optional<detail::AssumptionFacts> assumed;
    if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto exponent = make_product_or_one(coefficients);
    auto real_exponent = detail::query_real_value(exponent, facts, *context_);
    if (!real_exponent) { throw real_exponent.error(); }
    if (real_exponent.value() != Tribool::True || !defined_near_point(f->arguments()[0])) { return nullptr; }
    auto positive = detail::query_positive_value(f->arguments()[0], facts, *context_);
    if (!positive) { throw positive.error(); }
    if (positive.value() != Tribool::True) {
        LimitVisitor witness(var, point, direction, assumption_ctx_, context_);
        auto sign = witness.determine_sign_near_point(f->arguments()[0], direction);
        if (!sign || *sign <= 0) { return nullptr; }
    }
    return detail::make_node<PowerNode>(f->arguments()[0], exponent);
}

void LimitVisitor::visit(const FunctionNode& node) {
    if (is_limit_at_neg_infinity() && node.type() != FunctionNode::FuncType::Infinity) {
        auto value = handle_neg_infinity_limit(
            LMCAS::detail::make_node<FunctionNode>(node.type(), node.arguments()));
        if (value) {
            result = value;
            return;
        }
    }
    if (auto rewritten = positive_tail_log_power(node)) {
        result = eval_limit(rewritten);
        return;
    }
    if (node.type() == FunctionNode::FuncType::Sgn && node.arguments().size() == 1) {
        auto value = evaluate_sgn_limit(node.arguments()[0]);
        result = value ? *value : nullptr;
        return;
    }
    if (node.type() == FunctionNode::FuncType::Abs && node.arguments().size() == 1) {
        auto value = evaluate_abs_limit(node.arguments()[0]);
        if (value) {
            result = *value;
            return;
        }
    }
    if (composed_function_limit(node)) {
        return;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    for (const auto& argument : node.arguments()) {
        argument->accept(*this);
        if (!result) return;
        arguments.push_back(result);
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    LMCAS::detail::make_node<FunctionNode>(node.type(), arguments)->accept(norm);
    result = norm.get_result();
}

}
