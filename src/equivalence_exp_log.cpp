#include "internal/equivalence_support.hpp"

#include "assumption_context.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace LMCAS::equivalence_detail {
namespace {

std::vector<detail::SymbolicNodePtr> rewrite_operands(
    const std::vector<detail::SymbolicNodePtr>& operands,
    const AssumptionContext* assumptions, detail::RewriteBudget& budget) {
    std::vector<detail::SymbolicNodePtr> rewritten;
    budget.require_nodes(operands.size());
    std::size_t size = 0;
    rewritten.reserve(operands.size());
    for (const auto& operand : operands) {
        auto child = rewrite_exp_log_basic_identity(operand, assumptions, budget);
        size = budget.append_size(size, budget.measure(detail::node(child)));
        rewritten.push_back(detail::node(child));
    }
    return rewritten;
}

ExprPtr positive_log_argument(const detail::SymbolicNodePtr& node,
                              const AssumptionContext* assumptions) {
    auto inner_ln = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!inner_ln || inner_ln->type() != FunctionNode::FuncType::Ln ||
        inner_ln->arguments().size() != 1) {
        return nullptr;
    }
    auto argument = detail::make_expression_ptr(inner_ln->arguments()[0]);
    auto positive = assumptions
        ? assumptions->is_positive(*argument)
        : Result<Tribool>::success(Tribool::Unknown);
    const bool known_positive = inner_ln->arguments()[0]->is_positive() ||
        (positive && positive.value() == Tribool::True);
    return known_positive ? argument : nullptr;
}
std::optional<ExprPtr> rewrite_negative_lambert_exponential(
    const detail::SymbolicNodePtr& node, detail::RewriteBudget& budget) {
    const auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product || product->operands().size() != 2) {
        return std::nullopt;
    }
    detail::SymbolicNodePtr lambert_node;
    if (exact_integer_node(product->operands()[0], -1)) {
        lambert_node = product->operands()[1];
    } else if (exact_integer_node(product->operands()[1], -1)) {
        lambert_node = product->operands()[0];
    } else {
        return std::nullopt;
    }
    const auto lambert = std::dynamic_pointer_cast<const FunctionNode>(lambert_node);
    if (!lambert || lambert->type() != FunctionNode::FuncType::LambertW ||
        lambert->arguments().size() != 1) {
        return std::nullopt;
    }
    const auto argument =
        std::dynamic_pointer_cast<const NumberNode>(lambert->arguments()[0]);
    if (!argument || std::holds_alternative<lmmc_real_t>(argument->value()) ||
        !argument->is_positive()) {
        return std::nullopt;
    }
    if (argument->is_one()) {
        budget.measure(lambert_node);
        return detail::make_expression_ptr(lambert_node);
    }
    budget.require_nodes(3);
    auto reciprocal = detail::make_node<PowerNode>(
        lambert->arguments()[0], detail::make_node<NumberNode>(BigInt(-1)));
    auto rewritten = SymbolicFactory::create_multiply(
        {lambert_node, std::move(reciprocal)}, &budget);
    return normalize_equivalence(
        detail::make_expression_ptr(std::move(rewritten)), budget);
}

ExprPtr rewrite_exp_log_function(const FunctionNode& function,
                                 const AssumptionContext* assumptions,
                                 detail::RewriteBudget& budget) {
    auto arguments = rewrite_operands(function.arguments(), assumptions, budget);
    if (arguments.size() == 1) {
        if (function.type() == FunctionNode::FuncType::Exp) {
            if (exact_integer_node(arguments[0], 0)) {
                budget.require_nodes(1);
                return SymbolicExpr::number(1);
            }
            auto argument = positive_log_argument(arguments[0], assumptions);
            if (auto lambert =
                    rewrite_negative_lambert_exponential(arguments[0], budget)) {
                return std::move(*lambert);
            }
            if (argument) { return argument; }
        }
        if (function.type() == FunctionNode::FuncType::Ln &&
            exact_integer_node(arguments[0], 1)) {
            budget.require_nodes(1);
            return SymbolicExpr::number(0);
        }
    }
    normalization_check_children(&budget, 1, arguments);
    return detail::make_expression_ptr(detail::make_node<FunctionNode>(
        function.type(), std::move(arguments)));
}

}

ExprPtr rewrite_exp_log_basic_identity(const detail::SymbolicNodePtr& node,
                                       const AssumptionContext* assumptions,
                                       detail::RewriteBudget& budget) {
    if (!node) { return nullptr; }
    detail::RewriteScope scope(&budget);
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return detail::make_expression_ptr(SymbolicFactory::create_add(
            rewrite_operands(add->operands(), assumptions, budget), &budget));
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return detail::make_expression_ptr(SymbolicFactory::create_multiply(
            rewrite_operands(multiply->operands(), assumptions, budget), &budget));
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_exp_log_basic_identity(power->base(), assumptions, budget);
        auto exponent = rewrite_exp_log_basic_identity(power->exponent(), assumptions, budget);
        budget.append_size(budget.append_size(1, budget.measure(detail::node(base))),
                           budget.measure(detail::node(exponent)));
        return SymbolicExpr::power(base, exponent);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return rewrite_exp_log_function(*function, assumptions, budget);
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = rewrite_exp_log_basic_identity(complex_node->real(), assumptions, budget);
        auto imag = rewrite_exp_log_basic_identity(complex_node->imag(), assumptions, budget);
        return rewritten_complex(real, imag, node, budget);
    }
    budget.measure(node);
    return detail::make_expression_ptr(node);
}

}
