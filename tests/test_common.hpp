#pragma once
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "symbolic.hpp"
#include "assumption.hpp"
#include "internal/symbolic_ast.hpp"
#include "solve_strategies.hpp"
#include "residual_verification.hpp"
using namespace LMCAS;

inline std::shared_ptr<const SymbolicNode> test_variable_node(
    const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

inline std::shared_ptr<const SymbolicNode> test_integer_node(int value) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(value));
}

inline SymbolicExpr test_expression_from_node(
    std::shared_ptr<const SymbolicNode> node) {
    return LMCAS::detail::expression_from_node(std::move(node));
}

inline bool test_proved_equivalent(const std::shared_ptr<SymbolicExpr> &actual,
                                   const std::shared_ptr<SymbolicExpr> &expected) {
    ComputationContext context;
    auto proof = check_equivalent(actual, expected, context);
    return proof && std::holds_alternative<ProvedZeroResidual>(proof.value());
}

inline ::testing::AssertionResult test_expression_text(
    const std::shared_ptr<SymbolicExpr> &actual,
    const std::string &expected) {
    const std::string actual_text = actual ? actual->to_string() : "null";
    if (actual_text == expected) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "Actual: " << actual_text << "\nExpected: " << expected;
}

inline ::testing::AssertionResult test_expression_text(
    const std::shared_ptr<SymbolicExpr> &actual,
    const std::shared_ptr<SymbolicExpr> &expected) {
    return test_expression_text(
        actual, expected ? expected->to_string() : "null");
}

inline ::testing::AssertionResult test_same_expression(
    const std::shared_ptr<SymbolicExpr>& actual,
    const std::shared_ptr<SymbolicExpr>& expected) {
    if (!actual || !expected) {
        return ::testing::AssertionFailure() << "Actual and expected expressions require values";
    }
    const auto& actual_node = LMCAS::detail::node(*actual);
    const auto& expected_node = LMCAS::detail::node(*expected);
    if (!actual_node || !expected_node) {
        return ::testing::AssertionFailure() << "Actual and expected expressions require AST roots";
    }
    if (actual_node->compare(*expected_node) == 0) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "Actual: " << actual->to_string() << "\nExpected: " << expected->to_string();
}

/**
 * @brief 通过公开谓词查询符号，Zero 为 NonZero 的逻辑补集。
 * 模板将推导实现头文件的依赖限定在使用方。
 */
template <class Engine>
inline Result<Tribool> test_sign_query(
    const Engine &engine, const SymbolicExpr &expression, Sign target) {
    switch (target) {
    case Sign::Positive: {
        return engine.query_positive_checked(expression);
    }
    case Sign::Negative: {
        return engine.query_negative_checked(expression);
    }
    case Sign::NonNegative: {
        return engine.query_nonnegative_checked(expression);
    }
    case Sign::NonPositive: {
        return engine.query_nonpositive_checked(expression);
    }
    case Sign::NonZero: {
        return engine.query_nonzero_checked(expression);
    }
    case Sign::Zero: {
        auto nonzero = engine.query_nonzero_checked(expression);
        if (!nonzero) {
            return nonzero;
        }
        return nonzero.value() == Tribool::Unknown ? Tribool::Unknown : nonzero.value() == Tribool::True ? Tribool::False
                                                                                                         : Tribool::True;
    }
    }
    return Tribool::Unknown;
}

inline std::vector<std::shared_ptr<SymbolicExpr>> solve_vector_for_test(
    const std::shared_ptr<SymbolicExpr> &expression,
    const std::string &variable,
    const LMCAS::SolveOptions &options = {}) {
    auto result = LMCAS::solve_equation(expression, variable, options);
    if (!result) {
        throw std::runtime_error(
            "checked solve failed: " + result.error().message);
    }
    if (std::holds_alternative<LMCAS::EmptySolutions>(result.value())) {
        return {};
    }
    const auto *finite =
        std::get_if<LMCAS::FiniteSolutions>(&result.value());
    if (!finite) {
        throw std::runtime_error("solve result is not finitely enumerable");
    }
    std::vector<std::shared_ptr<SymbolicExpr>> values;
    for (const auto &solution : finite->values) {
        for (std::size_t copy = 0; copy < solution.multiplicity; ++copy) {
            values.push_back(solution.value);
        }
    }
    return values;
}

// Recursive numeric evaluator for symbolic expressions after substitution.
// Finite mode uses nullopt for unbound variables and invalid powers/functions.
// Algebra mode preserves NaN/infinity and evaluates negative odd roots
// for algebraic solution checks.
#include <optional>
#include <limits>

enum class TestNumericMode { Finite, Algebra };
inline std::optional<double> test_numeric_eval(
    const std::shared_ptr<SymbolicExpr> &expression,
    TestNumericMode mode = TestNumericMode::Finite);

inline double test_numeric_basic_function(FunctionNode::FuncType type, double argument) {
    using Type = FunctionNode::FuncType;
    switch (type) {
    case Type::Exp: {
        return std::exp(argument);
    }
    case Type::Ln: {
        return std::log(argument);
    }
    case Type::Sqrt: {
        return std::sqrt(argument);
    }
    case Type::Abs: {
        return std::abs(argument);
    }
    case Type::Sinh: {
        return std::sinh(argument);
    }
    case Type::Cosh: {
        return std::cosh(argument);
    }
    case Type::Tanh: {
        return std::tanh(argument);
    }
    case Type::ArcSin: {
        return std::asin(argument);
    }
    case Type::ArcCos: {
        return std::acos(argument);
    }
    case Type::ArcTan: {
        return std::atan(argument);
    }
    default: {
        return std::numeric_limits<double>::quiet_NaN();
    }
    }
}

inline double test_numeric_quotient(double numerator, double denominator) {
    if (denominator == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return numerator / denominator;
}

inline double test_numeric_function(FunctionNode::FuncType type, double argument) {
    using Type = FunctionNode::FuncType;
    switch (type) {
    case Type::Sin: {
        return std::sin(argument);
    }
    case Type::Cos: {
        return std::cos(argument);
    }
    case Type::Tan: {
        return std::tan(argument);
    }
    case Type::Sec: {
        return test_numeric_quotient(1.0, std::cos(argument));
    }
    case Type::Csc: {
        return test_numeric_quotient(1.0, std::sin(argument));
    }
    case Type::Cot: {
        return test_numeric_quotient(std::cos(argument), std::sin(argument));
    }
    default: {
        return test_numeric_basic_function(type, argument);
    }
    }
}

inline bool test_numeric_algebra_unsupported(FunctionNode::FuncType type) {
    using Type = FunctionNode::FuncType;
    return type == Type::Sinh || type == Type::Cosh || type == Type::Tanh ||
           type == Type::Sec || type == Type::Csc || type == Type::Cot;
}

inline std::optional<double> test_numeric_function_node(
    const FunctionNode &function, TestNumericMode mode) {
    if (function.arguments().size() != 1) {
        return mode == TestNumericMode::Algebra
                   ? std::optional<double>(std::nan("")) : std::nullopt;
    }
    auto argument = test_numeric_eval(
        LMCAS::detail::make_expression_ptr(function.arguments()[0]), mode);
    if (!argument) {
        return std::nullopt;
    }
    if (mode == TestNumericMode::Algebra &&
        test_numeric_algebra_unsupported(function.type())) {
        return std::nan("");
    }
    const double value = test_numeric_function(function.type(), *argument);
    if (mode == TestNumericMode::Finite && !std::isfinite(value)) {
        return std::nullopt;
    }
    return value;
}

template <bool Product>
inline std::optional<double> test_numeric_fold(
    const std::vector<std::shared_ptr<const SymbolicNode>> &operands,
    TestNumericMode mode) {
    double value = Product ? 1.0 : 0.0;
    for (const auto &operand : operands) {
        auto number = test_numeric_eval(
            LMCAS::detail::make_expression_ptr(operand), mode);
        if (!number) {
            return std::nullopt;
        }
        if constexpr (Product)
            value *= *number;
        else
            value += *number;
    }
    return value;
}

inline std::optional<double> test_numeric_power(
    const PowerNode &power, TestNumericMode mode) {
    auto base = test_numeric_eval(
        LMCAS::detail::make_expression_ptr(power.base()), mode);
    auto exponent = test_numeric_eval(
        LMCAS::detail::make_expression_ptr(power.exponent()), mode);
    if (!base || !exponent) {
        return std::nullopt;
    }
    if (mode == TestNumericMode::Algebra) {
        if (*base < 0.0 && std::abs(*exponent - std::round(*exponent)) > 1e-15) {
            double denominator = std::round(1.0 / *exponent);
            if (std::abs(*exponent * denominator - 1.0) < 1e-12 &&
                (static_cast<int>(denominator) % 2 == 1)) {
                return -std::pow(-*base, *exponent);
            }
            return std::nan("");
        }
        return std::pow(*base, *exponent);
    }
    if (*base == 0.0 && *exponent < 0.0) {
        return std::nullopt;
    }
    const double value = std::pow(*base, *exponent);
    if (!std::isfinite(value)) {
        return std::nullopt;
    }
    return value;
}

inline std::optional<double> test_numeric_number(const NumberNode &number) {
    if (std::holds_alternative<lmmc_real_t>(number.value())) {
        return std::get<lmmc_real_t>(number.value());
    }
    if (std::holds_alternative<BigInt>(number.value())) {
        return std::get<BigInt>(number.value()).to_double();
    }
    if (std::holds_alternative<Rational>(number.value())) {
        return std::get<Rational>(number.value()).to_double();
    }
    return std::nullopt;
}

inline std::optional<double> test_numeric_eval(
    const std::shared_ptr<SymbolicExpr> &expression, TestNumericMode mode) {
    if (!expression || !LMCAS::detail::node(expression)) {
        return mode == TestNumericMode::Algebra
                   ? std::optional<double>(0.0) : std::nullopt;
    }
    auto root = LMCAS::detail::node(expression);
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(root)) {
        return test_numeric_number(*number);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(root)) {
        return test_numeric_fold<false>(add->operands(), mode);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(root)) {
        return test_numeric_fold<true>(multiply->operands(), mode);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(root)) {
        return test_numeric_power(*power, mode);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(root)) {
        return test_numeric_function_node(*function, mode);
    }
    return mode == TestNumericMode::Algebra
               ? std::optional<double>(std::nan("")) : std::nullopt;
}

inline double test_numeric_value(const std::shared_ptr<SymbolicExpr> &expression) {
    return test_numeric_eval(expression, TestNumericMode::Algebra).value();
}
