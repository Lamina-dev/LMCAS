#include "internal/visitors/limit_visitor.hpp"

namespace LMCAS {

bool LimitVisitor::is_limit_variable(const std::shared_ptr<const SymbolicNode>& node) const {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(node);
    return variable && !variable->is_constant() && variable->name() == var;
}

bool LimitVisitor::is_positive_variable_power(const std::shared_ptr<const SymbolicNode>& node) const {
    if (is_limit_variable(node)) {
        return true;
    }
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    return is_limit_variable(power->base()) && exponent && exponent->is_positive();
}

bool LimitVisitor::is_variable_reciprocal(const std::shared_ptr<const SymbolicNode>& node) const {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power || !is_limit_variable(power->base())) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    return exponent && get_numeric_value(exponent) == -1.0;
}

bool LimitVisitor::is_log_variable(const std::shared_ptr<const SymbolicNode>& node) const {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!function || function->arguments().size() != 1) {
        return false;
    }
    if (function->type() != FunctionNode::FuncType::Ln &&
        function->type() != FunctionNode::FuncType::Log) {
        return false;
    }
    return is_limit_variable(function->arguments()[0]);
}

bool LimitVisitor::logarithmic_zero_product(const MultiplyNode& node) const {
    auto number = std::dynamic_pointer_cast<const NumberNode>(point);
    if (!number) {
        return false;
    }
    if (!number->is_zero() || direction != "+") {
        return false;
    }
    bool has_power = false;
    bool has_logarithm = false;
    for (const auto& operand : node.operands()) {
        if (!is_positive_variable_power(operand) && !is_log_variable(operand) &&
            !std::dynamic_pointer_cast<const NumberNode>(operand)) { return false; }
        has_power = has_power || is_positive_variable_power(operand);
        has_logarithm = has_logarithm || is_log_variable(operand);
    }
    return has_power && has_logarithm;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::one_plus_term(const AddNode& node) const {
    bool has_one = false;
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    for (const auto& operand : node.operands()) {
        if (!has_one && operand->is_one()) has_one = true;
        else terms.push_back(operand);
    }
    if (!has_one || terms.empty()) return nullptr;
    return terms.size() == 1 ? terms[0] : detail::make_node<AddNode>(std::move(terms));
}

std::shared_ptr<const SymbolicNode> LimitVisitor::variable_fraction_numerator(
    const std::shared_ptr<const SymbolicNode>& node) const {
    if (is_variable_reciprocal(node)) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return nullptr;
    }
    bool has_denominator = false;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    for (const auto& operand : product->operands()) {
        if (!has_denominator && is_variable_reciprocal(operand)) has_denominator = true;
        else factors.push_back(operand);
    }
    if (!has_denominator) {
        return nullptr;
    }
    return make_product_or_one(factors);
}

std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::logarithmic_remaining_limit(
    const MultiplyNode& node, const SymbolicNode* logarithm,
    const std::shared_ptr<const SymbolicNode>& numerator) {
    std::vector<std::shared_ptr<const SymbolicNode>> values{numerator};
    bool has_variable = false;
    for (const auto& operand : node.operands()) {
        if (operand.get() == logarithm) {
            continue;
        }
        if (!has_variable && is_limit_variable(operand)) {
            has_variable = true;
            continue;
        }
        auto value = eval_limit(operand);
        if (!value) {
            return std::nullopt;
        }
        values.push_back(value);
    }
    if (!has_variable) {
        return std::nullopt;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    make_product_or_one(values)->accept(norm);
    return norm.get_result();
}

bool LimitVisitor::logarithmic_product_limit(const MultiplyNode& node) {
    for (const auto& operand : node.operands()) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(operand);
        if (!function || function->arguments().size() != 1) {
            continue;
        }
        if (function->type() != FunctionNode::FuncType::Ln &&
            function->type() != FunctionNode::FuncType::Log) {
            continue;
        }
        auto sum = std::dynamic_pointer_cast<const AddNode>(function->arguments()[0]);
        if (!sum) {
            continue;
        }
        auto small_term = one_plus_term(*sum);
        if (!small_term) {
            continue;
        }
        if (!is_zero_limit(eval_limit(small_term))) continue;
        auto numerator = variable_fraction_numerator(small_term);
        if (!numerator) {
            continue;
        }
        auto value = logarithmic_remaining_limit(node, operand.get(), numerator);
        if (value) {
            result = *value;
            return true;
        }
    }
    return false;
}

bool LimitVisitor::has_slower_factors(const MultiplyNode& node, const SymbolicNode* exponential) const {
    for (const auto& operand : node.operands()) {
        if (operand.get() == exponential) {
            continue;
        }
        if (!polynomial_logarithmic_bound(operand)) {
            return false;
        }
    }
    return true;
}

bool LimitVisitor::decaying_exponential_product(const MultiplyNode& node) {
    for (const auto& operand : node.operands()) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(operand);
        if (!function || function->arguments().size() != 1 ||
            function->type() != FunctionNode::FuncType::Exp) {
            continue;
        }
        auto exponent_limit = eval_limit(function->arguments()[0]);
        if (!is_neg_inf(exponent_limit)) {
            continue;
        }
        if (proves_superlogarithmic_magnitude(function->arguments()[0]) &&
            has_slower_factors(node, operand.get())) {
            return true;
        }
    }
    return false;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::matching_reciprocal(
    const MultiplyNode& node, const std::shared_ptr<const SymbolicNode>& argument) const {
    for (const auto& operand : node.operands()) {
        auto power = std::dynamic_pointer_cast<const PowerNode>(operand);
        if (!power) {
            continue;
        }
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        if (!exponent || !power->base() || power->base()->compare(*argument) != 0) {
            continue;
        }
        if (get_numeric_value(exponent) == -1.0) {
            return operand;
        }
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::trigonometric_remaining_limit(
    const MultiplyNode& node, const SymbolicNode* function, const SymbolicNode* reciprocal) {
    std::vector<std::shared_ptr<const SymbolicNode>> values;
    for (const auto& operand : node.operands()) {
        if (operand.get() == function || operand.get() == reciprocal) {
            continue;
        }
        auto value = eval_limit(operand);
        if (!value) {
            return nullptr;
        }
        values.push_back(value);
    }
    values.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    make_product_or_one(values)->accept(norm);
    return norm.get_result();
}

bool LimitVisitor::trigonometric_product_limit(const MultiplyNode& node) {
    for (const auto& operand : node.operands()) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(operand);
        if (!function || function->arguments().size() != 1) {
            continue;
        }
        if (function->type() != FunctionNode::FuncType::Sin &&
            function->type() != FunctionNode::FuncType::Tan) {
            continue;
        }
        const auto& argument = function->arguments()[0];
        auto reciprocal = matching_reciprocal(node, argument);
        if (!reciprocal) {
            continue;
        }
        auto argument_limit = eval_limit(argument);
        if (is_zero_limit(argument_limit)) {
            result = trigonometric_remaining_limit(node, operand.get(), reciprocal.get());
            return result != nullptr;
        }
    }
    return false;
}

}
