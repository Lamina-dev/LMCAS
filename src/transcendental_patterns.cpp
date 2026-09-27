#include "internal/transcendental_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"

#include <cmath>

namespace LMCAS::detail::transcendental {
namespace {

bool invertible_function(const FunctionNode& function, const std::string&) {
    if (function.arguments().size() != 1) {
        return false;
    }
    using Type = FunctionNode::FuncType;
    switch (function.type()) {
    case Type::Sin: {
        return true;
    }
    case Type::Cos: {
        return true;
    }
    case Type::Tan: {
        return true;
    }
    case Type::Exp: {
        return true;
    }
    case Type::Ln: {
        return true;
    }
    default:
        return false;
    }
}

std::shared_ptr<SymbolicExpr> constant_sum(
    std::vector<std::shared_ptr<const SymbolicNode>> terms) {
    if (terms.empty()) {
        return SymbolicExpr::number(0);
    }
    if (terms.size() == 1) {
        return make_expression_ptr(terms[0]);
    }
    return make_expression_ptr(make_node<AddNode>(std::move(terms)));
}

struct ScaledFunction {
    std::shared_ptr<const FunctionNode> function;
    std::shared_ptr<const SymbolicNode> coefficient;
};

std::optional<ScaledFunction> scaled_function(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (function && invertible_function(*function, var) &&
        expression_depends_on_variable(node, var)) {
        return ScaledFunction{function, make_node<NumberNode>(BigInt(1))};
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return std::nullopt;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> coefficients;
    std::shared_ptr<const FunctionNode> found;
    for (const auto& operand : product->operands()) {
        auto candidate = std::dynamic_pointer_cast<const FunctionNode>(operand);
        if (candidate && invertible_function(*candidate, var) && !found &&
            expression_depends_on_variable(operand, var)) {
            found = candidate;
        } else {
            coefficients.push_back(operand);
        }
    }
    if (!found) {
        return std::nullopt;
    }
    auto coefficient = SymbolicFactory::create_multiply(std::move(coefficients));
    if (expression_depends_on_variable(coefficient, var)) {
        return std::nullopt;
    }
    return ScaledFunction{found, coefficient};
}

std::optional<InversePattern> additive_function_pattern(
    const AddNode& sum, const std::string& var) {
    std::optional<ScaledFunction> function;
    std::vector<std::shared_ptr<const SymbolicNode>> constants;
    for (const auto& operand : sum.operands()) {
        auto candidate = scaled_function(operand, var);
        if (candidate && !function) {
            function = std::move(candidate);
            continue;
        }
        if (expression_depends_on_variable(operand, var)) {
            return std::nullopt;
        }
        constants.push_back(operand);
    }
    if (!function) {
        return std::nullopt;
    }
    auto negative = SymbolicExpr::multiply(
        constant_sum(std::move(constants)), SymbolicExpr::number(-1));
    return InversePattern{function->function->type(),
        make_expression_ptr(function->function->arguments()[0]), negative,
        make_expression_ptr(function->coefficient)};
}

std::shared_ptr<const SymbolicNode> product_without(
    const MultiplyNode& product, std::size_t index) {
    std::vector<std::shared_ptr<const SymbolicNode>> remaining;
    for (std::size_t i = 0; i < product.operands().size(); ++i) {
        if (i != index) {
            remaining.push_back(product.operands()[i]);
        }
    }
    if (remaining.size() == 1) {
        return remaining[0];
    }
    return make_node<MultiplyNode>(std::move(remaining));
}

std::optional<std::shared_ptr<SymbolicExpr>> lambert_inner(
    const MultiplyNode& product, const std::string& var) {
    for (std::size_t i = 0; i < product.operands().size(); ++i) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(product.operands()[i]);
        if (!function || function->type() != FunctionNode::FuncType::Exp ||
            function->arguments().size() != 1) {
            continue;
        }
        const auto& argument = function->arguments()[0];
        auto remainder = product_without(product, i);
        const bool equal = !remainder || !argument
            ? remainder == argument : remainder->equals(*argument);
        if (equal && expression_depends_on_variable(remainder, var)) {
            return make_expression_ptr(remainder);
        }
    }
    return std::nullopt;
}

std::optional<LambertWPattern> opposite_exponential_lambert_pattern(
    const AddNode& sum, const std::string& var) {
    if (sum.operands().size() != 2) {
        return std::nullopt;
    }
    for (std::size_t index = 0; index < 2; ++index) {
        auto exponential = std::dynamic_pointer_cast<const FunctionNode>(
            sum.operands()[index]);
        if (!exponential ||
            exponential->type() != FunctionNode::FuncType::Exp ||
            exponential->arguments().size() != 1) {
            continue;
        }
        const auto& argument = exponential->arguments()[0];
        const auto& other = sum.operands()[1 - index];
        if (!argument || !other || !argument->equals(*other) ||
            !expression_depends_on_variable(argument, var)) {
            continue;
        }
        auto inner = SymbolicExpr::multiply(
            SymbolicExpr::number(-1), make_expression_ptr(argument))->simplify();
        return LambertWPattern{std::move(inner), SymbolicExpr::number(1)};
    }
    return std::nullopt;
}

std::optional<LambertWPattern> additive_lambert_pattern(
    const AddNode& sum, const std::string& var) {
    if (auto direct = opposite_exponential_lambert_pattern(sum, var)) {
        return direct;
    }
    std::optional<std::shared_ptr<SymbolicExpr>> found;
    std::vector<std::shared_ptr<const SymbolicNode>> constants;
    bool other_variables = false;
    for (const auto& operand : sum.operands()) {
        auto product = std::dynamic_pointer_cast<const MultiplyNode>(operand);
        if (product) {
            auto candidate = lambert_inner(*product, var);
            if (candidate && !found) {
                found = candidate;
                continue;
            }
        }
        if (!expression_depends_on_variable(operand, var)) {
            constants.push_back(operand);
        } else {
            other_variables = true;
        }
    }
    if (!found || other_variables) {
        return std::nullopt;
    }
    auto rhs = SymbolicExpr::multiply(
        constant_sum(std::move(constants)), SymbolicExpr::number(-1))->simplify();
    return LambertWPattern{*found, rhs};
}

std::optional<ExpBasePattern> power_pattern(const PowerNode& power, const std::string& var) {
    if (expression_depends_on_variable(power.base(), var) ||
        !expression_depends_on_variable(power.exponent(), var)) {
        return std::nullopt;
    }
    auto base = make_expression_ptr(power.base());
    auto exponent = make_expression_ptr(power.exponent());
    lmmc_real_t value;
    if (try_evaluate_numeric(base, value) && value > 0 &&
        std::abs(value - 1.0) > LMMC_REAL_EPSILON) {
        return ExpBasePattern{base, exponent, nullptr};
    }
    return std::nullopt;
}

std::optional<ExpBasePattern> additive_power_pattern(const AddNode& sum, const std::string& var) {
    std::optional<ExpBasePattern> found;
    std::vector<std::shared_ptr<const SymbolicNode>> constants;
    bool other_variables = false;
    for (const auto& operand : sum.operands()) {
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(operand)) {
            auto candidate = power_pattern(*power, var);
            if (candidate && !found) {
                found = std::move(candidate);
                continue;
            }
        }
        if (!expression_depends_on_variable(operand, var)) {
            constants.push_back(operand);
        } else {
            other_variables = true;
        }
    }
    if (!found || other_variables) {
        return std::nullopt;
    }
    found->rhs = SymbolicExpr::multiply(
        constant_sum(std::move(constants)), SymbolicExpr::number(-1))->simplify();
    return found;
}

}

std::optional<InversePattern> decompose_trig_exp_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !node(expr)) {
        return std::nullopt;
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node(expr))) {
        if (invertible_function(*function, var)) {
            return InversePattern{function->type(),
                make_expression_ptr(function->arguments()[0]), SymbolicExpr::number(0), nullptr};
        }
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node(expr))) {
        return additive_function_pattern(*sum, var);
    }
    return std::nullopt;
}

std::optional<LambertWPattern> decompose_lambert_w_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !node(expr)) {
        return std::nullopt;
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node(expr))) {
        auto inner = lambert_inner(*product, var);
        if (inner) {
            return LambertWPattern{*inner, SymbolicExpr::number(0)};
        }
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node(expr))) {
        return additive_lambert_pattern(*sum, var);
    }
    return std::nullopt;
}

std::optional<ExpBasePattern> decompose_exp_base_pattern(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !node(expr)) {
        return std::nullopt;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node(expr))) {
        auto pattern = power_pattern(*power, var);
        if (pattern) {
            pattern->rhs = SymbolicExpr::number(0);
            return pattern;
        }
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node(expr))) {
        return additive_power_pattern(*sum, var);
    }
    return std::nullopt;
}

}
