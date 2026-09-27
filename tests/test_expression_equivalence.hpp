#pragma once

#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>
#include <vector>

using namespace LMCAS;
inline bool test_is_minus_one_number(const std::shared_ptr<const SymbolicNode> &node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) {
        return false;
    }
    if (std::holds_alternative<BigInt>(num->value())) {
        return std::get<BigInt>(num->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return std::get<Rational>(num->value()) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        return std::get<lmmc_real_t>(num->value()) == -1.0;
    }
    return false;
}

inline void test_mark_inverse_pair(
    const std::vector<std::shared_ptr<const SymbolicNode>> &operands,
    std::vector<bool> &used, std::size_t index) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(operands[index]);
    if (!power || !test_is_minus_one_number(power->exponent()))
        return;
    for (std::size_t other = 0; other < operands.size(); ++other) {
        if (index == other || used[other])
            continue;
        if (power->base()->compare(*operands[other]) == 0) {
            used[index] = true;
            used[other] = true;
            return;
        }
    }
}

inline std::shared_ptr<const SymbolicNode> test_cancel_inverse_factors(
    const std::vector<std::shared_ptr<const SymbolicNode>> &operands) {
    std::vector<bool> used(operands.size(), false);
    for (std::size_t index = 0; index < operands.size(); ++index) {
        if (!used[index])
            test_mark_inverse_pair(operands, used, index);
    }
    std::vector<std::shared_ptr<const SymbolicNode>> kept;
    kept.reserve(operands.size());
    for (std::size_t index = 0; index < operands.size(); ++index) {
        if (!used[index])
            kept.push_back(operands[index]);
    }
    if (kept.empty()) {
        return LMCAS::detail::node(SymbolicExpr::number(1));
    }
    if (kept.size() == 1) {
        return kept[0];
    }
    return LMCAS::detail::make_node<MultiplyNode>(kept);
}

inline std::shared_ptr<const SymbolicNode> test_cancel_inverse_products_node(
    const std::shared_ptr<const SymbolicNode> &node) {
    if (!node) {
        return node;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.reserve(add->operands().size());
        for (const auto &op : add->operands()) {
            ops.push_back(test_cancel_inverse_products_node(op));
        }
        return LMCAS::detail::make_node<AddNode>(ops);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.reserve(mul->operands().size());
        for (const auto &op : mul->operands()) {
            ops.push_back(test_cancel_inverse_products_node(op));
        }

        return test_cancel_inverse_factors(ops);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return LMCAS::detail::make_node<PowerNode>(
            test_cancel_inverse_products_node(pow->base()),
            test_cancel_inverse_products_node(pow->exponent()));
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(func->arguments().size());
        for (const auto &arg : func->arguments()) {
            args.push_back(test_cancel_inverse_products_node(arg));
        }
        return LMCAS::detail::make_node<FunctionNode>(func->type(), args);
    }

    return node;
}

inline std::shared_ptr<SymbolicExpr> test_cancel_inverse_products(
    const std::shared_ptr<SymbolicExpr> &expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return expr;
    }
    return LMCAS::detail::make_expression_ptr(
        test_cancel_inverse_products_node(LMCAS::detail::node(expr)));
}

inline std::shared_ptr<SymbolicExpr> test_simplified_zero(
    const std::shared_ptr<SymbolicExpr> &expression) {
    if (!expression) {
        return nullptr;
    }
    auto simplified = expression->simplify();
    return simplified && simplified->is_zero() ? simplified : nullptr;
}

inline std::shared_ptr<SymbolicExpr> test_expanded_zero(
    const std::shared_ptr<SymbolicExpr> &expression) {
    auto expanded = expression->expand();
    if (!expanded) {
        return nullptr;
    }
    auto simplified = expanded->simplify();
    if (simplified && simplified->is_zero()) {
        return simplified;
    }
    const auto &candidate = simplified ? simplified : expanded;
    auto inverse_zero = test_simplified_zero(test_cancel_inverse_products(candidate));
    if (inverse_zero) {
        return inverse_zero;
    }
    return test_simplified_zero(candidate->cancel());
}

inline std::shared_ptr<SymbolicExpr> test_normalized_delta(
    const std::shared_ptr<SymbolicExpr> &actual,
    const std::shared_ptr<SymbolicExpr> &expected) {
    if (!actual || !expected) {
        return nullptr;
    }

    auto neg_expected = SymbolicExpr::multiply(SymbolicExpr::number(-1), expected);
    auto delta = SymbolicExpr::add(actual, neg_expected);
    if (!delta) {
        return nullptr;
    }

    auto simplified = delta->simplify();
    if (simplified && simplified->is_zero()) {
        return simplified;
    }

    auto inverse_zero = test_simplified_zero(
        test_cancel_inverse_products(simplified ? simplified : delta));
    if (inverse_zero) {
        return inverse_zero;
    }
    auto expanded_zero = test_expanded_zero(delta);
    if (expanded_zero) {
        return expanded_zero;
    }
    auto cancelled_zero = test_simplified_zero(
        (simplified ? simplified : delta)->cancel());
    if (cancelled_zero) {
        return cancelled_zero;
    }
    return simplified;
}

inline bool test_expr_equivalent(
    const std::shared_ptr<SymbolicExpr> &actual,
    const std::shared_ptr<SymbolicExpr> &expected) {
    auto delta = test_normalized_delta(actual, expected);
    return delta && delta->is_zero();
}
