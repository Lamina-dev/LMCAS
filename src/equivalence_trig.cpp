#include "internal/equivalence_support.hpp"

#include <utility>
#include <vector>

namespace LMCAS::equivalence_detail {
namespace {

struct TrigSquare {
    FunctionNode::FuncType type;
    detail::SymbolicNodePtr argument;
    detail::SymbolicNodePtr coefficient;
};

std::optional<TrigSquare> trig_square(
    const detail::SymbolicNodePtr& node) {
    auto candidate = node;
    detail::SymbolicNodePtr coefficient =
        detail::make_node<NumberNode>(BigInt(1));
    if (const auto multiply =
            std::dynamic_pointer_cast<const MultiplyNode>(node);
        multiply && multiply->operands().size() == 2) {
        const auto first_number =
            std::dynamic_pointer_cast<const NumberNode>(
                multiply->operands()[0]);
        const auto second_number =
            std::dynamic_pointer_cast<const NumberNode>(
                multiply->operands()[1]);
        if (first_number) {
            coefficient = first_number;
            candidate = multiply->operands()[1];
        } else if (second_number) {
            coefficient = second_number;
            candidate = multiply->operands()[0];
        }
    }
    const auto power =
        std::dynamic_pointer_cast<const PowerNode>(candidate);
    if (!power || !exact_integer_node(power->exponent(), 2)) {
        return std::nullopt;
    }
    const auto function =
        std::dynamic_pointer_cast<const FunctionNode>(power->base());
    if (!function || function->arguments().size() != 1 ||
        (function->type() != FunctionNode::FuncType::Sin &&
         function->type() != FunctionNode::FuncType::Cos)) {
        return std::nullopt;
    }
    return TrigSquare{
        function->type(), function->arguments()[0], coefficient};
}

std::optional<ExprPtr> unwrap_trig_negated_argument(
    const detail::SymbolicNodePtr& node, detail::RewriteBudget& budget) {
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) { return std::nullopt; }

    bool found_negative_one = false;
    std::vector<detail::SymbolicNodePtr> remaining;
    budget.require_nodes(multiply->operands().size());
    std::size_t size = 0;
    remaining.reserve(multiply->operands().size());
    for (const auto& operand : multiply->operands()) {
        if (!found_negative_one && exact_integer_node(operand, -1)) {
            found_negative_one = true;
            continue;
        }
        size = budget.append_size(size, budget.measure(operand));
        remaining.push_back(operand);
    }
    if (!found_negative_one || remaining.empty()) { return std::nullopt; }
    if (remaining.size() == 1) {
        return detail::make_expression_ptr(remaining.front());
    }
    return detail::make_expression_ptr(
        SymbolicFactory::create_multiply(std::move(remaining), &budget));
}

std::vector<detail::SymbolicNodePtr> rewrite_operands(
    const std::vector<detail::SymbolicNodePtr>& operands, detail::RewriteBudget& budget) {
    std::vector<detail::SymbolicNodePtr> rewritten;
    budget.require_nodes(operands.size());
    std::size_t size = 0;
    rewritten.reserve(operands.size());
    for (const auto& operand : operands) {
        auto child = rewrite_trig_basic_identity(operand, budget);
        size = budget.append_size(size, budget.measure(detail::node(child)));
        rewritten.push_back(detail::node(child));
    }
    return rewritten;
}

std::optional<std::pair<std::size_t, detail::SymbolicNodePtr>>
complementary_square(
    const std::vector<detail::SymbolicNodePtr>& nodes,
    const std::vector<bool>& used, std::size_t index) {
    const auto square = trig_square(nodes[index]);
    if (!square) { return std::nullopt; }
    const auto complement =
        square->type == FunctionNode::FuncType::Sin
        ? FunctionNode::FuncType::Cos
        : FunctionNode::FuncType::Sin;
    for (std::size_t other = index + 1; other < nodes.size(); ++other) {
        if (used[other]) { continue; }
        const auto candidate = trig_square(nodes[other]);
        if (candidate && candidate->type == complement &&
            square->argument->equals(*candidate->argument) &&
            square->coefficient->equals(*candidate->coefficient)) {
            return std::pair<std::size_t, detail::SymbolicNodePtr>{
                other, square->coefficient};
        }
    }
    return std::nullopt;
}

ExprPtr rewrite_trig_sum(const AddNode& add, detail::RewriteBudget& budget) {
    auto rewritten = rewrite_operands(add.operands(), budget);
    std::vector<bool> used(rewritten.size(), false);
    std::vector<detail::SymbolicNodePtr> result_nodes;
    std::size_t size = 0;
    for (std::size_t index = 0; index < rewritten.size(); ++index) {
        if (used[index]) { continue; }
        auto partner = complementary_square(rewritten, used, index);
        if (partner) {
            used[index] = true;
            used[partner->first] = true;
            size = budget.append_size(
                size, budget.measure(partner->second));
            result_nodes.push_back(partner->second);
        } else {
            size = budget.append_size(size, budget.measure(rewritten[index]));
            result_nodes.push_back(rewritten[index]);
        }
    }
    if (result_nodes.empty()) {
        budget.require_nodes(1);
        return SymbolicExpr::number(0);
    }
    return detail::make_expression_ptr(
        SymbolicFactory::create_add(std::move(result_nodes), &budget));
}

ExprPtr rewrite_trig_function(const FunctionNode& function, detail::RewriteBudget& budget) {
    auto arguments = rewrite_operands(function.arguments(), budget);
    if (arguments.size() == 1) {
        auto positive_argument = unwrap_trig_negated_argument(arguments[0], budget);
        if (positive_argument && *positive_argument) {
            if (function.type() == FunctionNode::FuncType::Sin) {
                budget.append_size(3, budget.measure(detail::node(*positive_argument)));
                return SymbolicExpr::multiply(
                    SymbolicExpr::number(-1),
                    SymbolicExpr::sin(*positive_argument));
            }
            if (function.type() == FunctionNode::FuncType::Cos) {
                budget.append_size(1, budget.measure(detail::node(*positive_argument)));
                return SymbolicExpr::cos(*positive_argument);
            }
        }
    }
    normalization_check_children(&budget, 1, arguments);
    return detail::make_expression_ptr(
        detail::make_node<FunctionNode>(function.type(), std::move(arguments)));
}

}

ExprPtr rewrite_trig_basic_identity(const detail::SymbolicNodePtr& node,
                                     detail::RewriteBudget& budget) {
    if (!node) { return nullptr; }
    detail::RewriteScope scope(&budget);
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return rewrite_trig_sum(*add, budget);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return detail::make_expression_ptr(SymbolicFactory::create_multiply(
            rewrite_operands(multiply->operands(), budget), &budget));
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_trig_basic_identity(power->base(), budget);
        auto exponent = rewrite_trig_basic_identity(power->exponent(), budget);
        budget.append_size(budget.append_size(1, budget.measure(detail::node(base))),
                           budget.measure(detail::node(exponent)));
        return SymbolicExpr::power(base, exponent);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return rewrite_trig_function(*function, budget);
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = rewrite_trig_basic_identity(complex_node->real(), budget);
        auto imag = rewrite_trig_basic_identity(complex_node->imag(), budget);
        return rewritten_complex(real, imag, node, budget);
    }
    budget.measure(node);
    return detail::make_expression_ptr(node);
}

}
