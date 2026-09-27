#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

bool NormalizationVisitor::normalize_logarithm(
    FunctionNode::FuncType type,
    const std::vector<std::shared_ptr<const SymbolicNode>>& arguments) {
    if (type == FunctionNode::FuncType::Ln && arguments.size() == 1) {
        return normalize_natural_logarithm(arguments.front());
    }
    if (type != FunctionNode::FuncType::Log || arguments.size() != 2) {
        return false;
    }
    if (const auto power = std::dynamic_pointer_cast<const PowerNode>(arguments.front())) {
        normalization_check_children(rewrite_budget(), 1, arguments[1]);
        auto logarithm = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{arguments[1]});
        if (power->base()->compare(*arguments[1]) == 0 &&
            normalization_proved(detail::query_positive_value(arguments[1], facts_, context_)) &&
            normalization_proved(detail::query_real_value(power->exponent(), facts_, context_)) &&
            normalization_proved(detail::query_nonzero_value(logarithm, facts_, Domain::Real, context_))) {
            set_result(power->exponent());
            return true;
        }
    }
    normalization_check_children(rewrite_budget(), 1, arguments[0]);
    auto numerator = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{arguments[0]});
    normalization_check_children(rewrite_budget(), 1, arguments[1]);
    auto denominator = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{arguments[1]});
    normalization_check_children(rewrite_budget(), 2, denominator);
    auto reciprocal = detail::make_node<PowerNode>(denominator, detail::make_node<NumberNode>(BigInt(-1)));
    auto quotient = make_normalized_multiply_node({numerator, reciprocal}, rewrite_budget());
    NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
    quotient->accept(visitor);
    set_result(visitor.get_result());
    return true;
}

bool NormalizationVisitor::normalize_natural_logarithm(
    const std::shared_ptr<const SymbolicNode>& argument) {
    if (const auto power = std::dynamic_pointer_cast<const PowerNode>(argument)) {
        if (!normalization_proved(detail::query_positive_value(power->base(), facts_, context_)) ||
            !normalization_proved(detail::query_real_value(power->exponent(), facts_, context_))) {
            return false;
        }
        normalization_check_children(rewrite_budget(), 1, power->base());
        auto logarithm = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{power->base()});
        auto product = make_normalized_multiply_node({power->exponent(), logarithm}, rewrite_budget());
        NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
        product->accept(visitor);
        set_result(visitor.get_result());
        return true;
    }
    const auto function = std::dynamic_pointer_cast<const FunctionNode>(argument);
    if (!function || function->type() != FunctionNode::FuncType::Exp) {
        return false;
    }
    if (function->arguments().size() != 1 ||
        !normalization_proved(detail::query_real_value(function->arguments().front(), facts_, context_))) {
        return false;
    }
    set_result(function->arguments().front());
    return true;
}

bool NormalizationVisitor::normalize_function_parity(
    FunctionNode::FuncType type, const std::shared_ptr<const SymbolicNode>& argument) {
    std::shared_ptr<const SymbolicNode> positive;
    if (!check_negative_arg(argument, positive, rewrite_budget())) {
        return false;
    }
    switch (type) {
        case FunctionNode::FuncType::Sin:
        case FunctionNode::FuncType::Tan:
        case FunctionNode::FuncType::ArcTan:
        case FunctionNode::FuncType::ArcSin: {
            normalization_check_children(rewrite_budget(), 1, positive);
            auto function = detail::make_node<FunctionNode>(type,
                std::vector<std::shared_ptr<const SymbolicNode>>{positive});
            set_result(make_normalized_multiply_node(
                {detail::make_node<NumberNode>(BigInt(-1)), function}, rewrite_budget()));
            return true;
        }
        case FunctionNode::FuncType::Cos:
            normalization_check_children(rewrite_budget(), 1, positive);
            set_result(detail::make_node<FunctionNode>(type,
                std::vector<std::shared_ptr<const SymbolicNode>>{positive}));
            return true;
        default:
            return false;
    }
}

bool NormalizationVisitor::normalize_unary_function(
    FunctionNode::FuncType type, const std::shared_ptr<const SymbolicNode>& argument) {
    Rational coefficient;
    if (get_pi_coeff(argument, coefficient)) {
        if (auto value = normalization_pi_function(type, coefficient, rewrite_budget())) {
            set_result(std::move(value));
            return true;
        }
    }
    return normalize_function_parity(type, argument);
}

void NormalizationVisitor::visit(const FunctionNode& node) {
    std::shared_ptr<const SymbolicNode> argument;
    if (try_normalize_squared_norm(node, argument)) {
        return;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    std::size_t nodes = 0;
    if (argument) {
        normalization_append(rewrite_budget(), nodes, arguments, argument);
    } else {
        for (const auto& operand : node.arguments()) {
            operand->accept(*this);
            normalization_append(rewrite_budget(), nodes, arguments, result);
        }
    }
    if (arguments.size() == 1) {
        if (const auto number = std::dynamic_pointer_cast<const NumberNode>(arguments.front())) {
            if (auto value = normalization_numeric_function(node.type(), number, rewrite_budget())) {
                set_result(std::move(value));
                return;
            }
        }
    }
    if (normalize_logarithm(node.type(), arguments)) {
        return;
    }
    if (arguments.size() == 1 && normalize_unary_function(node.type(), arguments.front())) {
        return;
    }
    if (rewrite_budget()) { rewrite_budget()->append_size(1, nodes); }
    auto function = detail::make_node<FunctionNode>(node.type(), arguments);
    if (auto simplified = try_assumption_simplify(function)) {
        simplified->accept(*this);
        return;
    }
    set_result(function);
}

void NormalizationVisitor::visit(const UninterpretedFunctionNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    std::size_t nodes = 1;
    normalization_check_count(rewrite_budget(), node.arguments().size());
    arguments.reserve(node.arguments().size());
    for (const auto& argument : node.arguments()) {
        argument->accept(*this);
        normalization_append(rewrite_budget(), nodes, arguments, result);
    }
    set_result(LMCAS::detail::make_node<UninterpretedFunctionNode>(
        node.name(), std::move(arguments)));
}
std::shared_ptr<const SymbolicNode> NormalizationVisitor::simplify_assumed_square_root(
    const std::shared_ptr<const SymbolicNode>& argument) const {
    const auto power = std::dynamic_pointer_cast<const PowerNode>(argument);
    if (!power) {
        return nullptr;
    }
    const auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    BigInt value;
    if (!try_get_integer_value(exponent, value)) {
        return nullptr;
    }
    if (value != BigInt(2)) {
        return nullptr;
    }
    auto nonnegative = detail::query_nonnegative_value(power->base(), facts_, context_);
    if (normalization_proved(nonnegative)) {
        return power->base();
    }
    auto real = detail::query_real_value(power->base(), facts_, context_);
    if (!real) {
        throw real.error();
    }
    if (real.value() != Tribool::True) {
        return nullptr;
    }
    normalization_check_children(rewrite_budget(), 1, power->base());
    return detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{power->base()});
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::simplify_assumed_absolute(
    const std::shared_ptr<const SymbolicNode>& argument) const {
    auto positive = detail::query_nonnegative_value(argument, facts_, context_);
    if (normalization_proved(positive)) {
        return argument;
    }
    if (!normalization_proved(detail::query_real_value(argument, facts_, context_)) ||
        !normalization_proved(facts_.is_negative(argument, context_))) {
        return nullptr;
    }
    return make_normalized_multiply_node({
        detail::make_node<NumberNode>(BigInt(-1)), argument}, rewrite_budget());
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::try_assumption_simplify(
    const std::shared_ptr<const SymbolicNode>& node) {
    const auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!function) {
        return nullptr;
    }
    if (function->arguments().size() != 1) {
        return nullptr;
    }
    const auto& argument = function->arguments().front();
    switch (function->type()) {
        case FunctionNode::FuncType::Sqrt:
            return simplify_assumed_square_root(argument);
        case FunctionNode::FuncType::Abs:
            return simplify_assumed_absolute(argument);
        default:
            return nullptr;
    }
}

}
