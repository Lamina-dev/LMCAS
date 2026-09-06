#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "solve_strategies.hpp"
using namespace LMCAS;

inline int g_failures = 0;
inline int g_passes = 0;

inline std::vector<std::shared_ptr<SymbolicExpr>> solve_vector_for_test(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const LMCAS::SolveOptions& options = {}) {
    auto result = LMCAS::solve_equation(expression, variable, options);
    if (!result) {
        if (result.error().code == LMCAS::CasErrc::Inconclusive) return {};
        throw std::runtime_error(
            "checked solve failed: " + result.error().message);
    }
    if (std::holds_alternative<LMCAS::EmptySolutions>(result.value())) return {};
    const auto* finite =
        std::get_if<LMCAS::FiniteSolutions>(&result.value());
    if (!finite) {
        throw std::runtime_error("solve result is not finitely enumerable");
    }
    std::vector<std::shared_ptr<SymbolicExpr>> values;
    for (const auto& solution : finite->values) {
        for (std::size_t copy = 0; copy < solution.multiplicity; ++copy) {
            values.push_back(solution.value);
        }
    }
    return values;
}

inline void TEST_CASE(const std::string& name) {
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Test Case: " << name << std::endl;
}

inline void EXPECT_EQ_STR(const std::string& actual, const std::string& expected, const std::string& msg) {
    if (actual != expected) {
        std::cerr << "[FAIL] " << msg << "\n  Expected: " << expected << "\n  Got:      " << actual << std::endl;
        g_failures++;
    } else {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    }
}

inline void EXPECT_EQ_EXPR(const std::shared_ptr<SymbolicExpr>& actual, const std::shared_ptr<SymbolicExpr>& expected, const std::string& msg) {
    std::string s_actual = actual ? actual->to_string() : "null";
    std::string s_expected = expected ? expected->to_string() : "null";
    EXPECT_EQ_STR(s_actual, s_expected, msg);
}

inline void EXPECT_EQ_EXPR_STR(const std::shared_ptr<SymbolicExpr>& actual, const std::string& expected_str, const std::string& msg) {
    std::string s_actual = actual ? actual->to_string() : "null";
    EXPECT_EQ_STR(s_actual, expected_str, msg);
}

inline void EXPECT_TRUE(bool condition, const std::string& msg) {
    if (!condition) {
        std::cerr << "[FAIL] " << msg << " is expected to be TRUE but is FALSE" << std::endl;
        g_failures++;
    } else {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    }
}

inline void EXPECT_FALSE(bool condition, const std::string& msg) {
    if (condition) {
        std::cerr << "[FAIL] " << msg << " is expected to be FALSE but is TRUE" << std::endl;
        g_failures++;
    } else {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    }
}

inline void EXPECT_NEAR(double actual, double expected, double tolerance, const std::string& msg) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "[FAIL] " << msg << "\n  Expected: " << expected
                  << " (+/- " << tolerance << ")\n  Got:      " << actual << std::endl;
        g_failures++;
    } else {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    }
}

inline void EXPECT_CONTAINS(const std::string& actual, const std::vector<std::string>& tokens, const std::string& msg) {
    for(const auto& t : tokens) {
        if (actual.find(t) == std::string::npos) {
            std::cerr << "[FAIL] " << msg << " | Missing token: " << t << "\n  In: " << actual << std::endl;
            g_failures++;
            return;
        }
    }
    std::cout << "[PASS] " << msg << std::endl;
    g_passes++;
}

inline bool test_is_minus_one_number(const std::shared_ptr<const SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        return std::get<BigInt>(num->value()).to_int() == -1;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return std::get<Rational>(num->value()).to_double() == -1.0;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        return std::get<lmmc_real_t>(num->value()) == -1.0;
    }
    return false;
}

inline std::shared_ptr<const SymbolicNode> test_cancel_inverse_products_node(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return node;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.reserve(add->operands().size());
        for (const auto& op : add->operands()) {
            ops.push_back(test_cancel_inverse_products_node(op));
        }
        return LMCAS::detail::make_node<AddNode>(ops);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.reserve(mul->operands().size());
        for (const auto& op : mul->operands()) {
            ops.push_back(test_cancel_inverse_products_node(op));
        }

        std::vector<bool> used(ops.size(), false);
        for (size_t i = 0; i < ops.size(); ++i) {
            if (used[i]) continue;
            std::shared_ptr<const SymbolicNode> inverse_base;
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(ops[i])) {
                if (test_is_minus_one_number(pow->exponent())) {
                    inverse_base = pow->base();
                }
            }
            if (!inverse_base) continue;

            for (size_t j = 0; j < ops.size(); ++j) {
                if (i == j || used[j]) continue;
                if (inverse_base->compare(*ops[j]) == 0) {
                    used[i] = true;
                    used[j] = true;
                    break;
                }
            }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> kept;
        for (size_t i = 0; i < ops.size(); ++i) {
            if (!used[i]) kept.push_back(ops[i]);
        }
        if (kept.empty()) return LMCAS::detail::node(SymbolicExpr::number(1));
        if (kept.size() == 1) return kept[0];
        return LMCAS::detail::make_node<MultiplyNode>(kept);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return LMCAS::detail::make_node<PowerNode>(
            test_cancel_inverse_products_node(pow->base()),
            test_cancel_inverse_products_node(pow->exponent()));
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(func->arguments().size());
        for (const auto& arg : func->arguments()) {
            args.push_back(test_cancel_inverse_products_node(arg));
        }
        return LMCAS::detail::make_node<FunctionNode>(func->type(), args);
    }

    return node;
}

inline std::shared_ptr<SymbolicExpr> test_cancel_inverse_products(
    const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return expr;
    return LMCAS::detail::make_expression_ptr(
        test_cancel_inverse_products_node(LMCAS::detail::node(expr)));
}

inline std::shared_ptr<SymbolicExpr> test_normalized_delta(
    const std::shared_ptr<SymbolicExpr>& actual,
    const std::shared_ptr<SymbolicExpr>& expected) {
    if (!actual || !expected) return nullptr;

    auto neg_expected = SymbolicExpr::multiply(SymbolicExpr::number(-1), expected);
    auto delta = SymbolicExpr::add(actual, neg_expected);
    if (!delta) return nullptr;

    auto simplified = delta->simplify();
    if (simplified && simplified->is_zero()) return simplified;

    auto inverse_cancelled = test_cancel_inverse_products(simplified ? simplified : delta);
    if (inverse_cancelled) {
        auto inverse_simplified = inverse_cancelled->simplify();
        if (inverse_simplified && inverse_simplified->is_zero()) {
            return inverse_simplified;
        }
    }

    auto expanded = delta->expand();
    if (expanded) {
        auto expanded_simplified = expanded->simplify();
        if (expanded_simplified && expanded_simplified->is_zero()) {
            return expanded_simplified;
        }

        auto expanded_inverse = test_cancel_inverse_products(expanded_simplified ? expanded_simplified : expanded);
        if (expanded_inverse) {
            auto expanded_inverse_simplified = expanded_inverse->simplify();
            if (expanded_inverse_simplified && expanded_inverse_simplified->is_zero()) {
                return expanded_inverse_simplified;
            }
        }

        auto cancelled = expanded_simplified ? expanded_simplified->cancel() : expanded->cancel();
        if (cancelled) {
            auto cancelled_simplified = cancelled->simplify();
            if (cancelled_simplified && cancelled_simplified->is_zero()) {
                return cancelled_simplified;
            }
        }
    }

    auto cancelled = simplified ? simplified->cancel() : delta->cancel();
    if (cancelled) {
        auto cancelled_simplified = cancelled->simplify();
        if (cancelled_simplified && cancelled_simplified->is_zero()) {
            return cancelled_simplified;
        }
    }
    return simplified;
}

inline bool test_expr_equivalent(
    const std::shared_ptr<SymbolicExpr>& actual,
    const std::shared_ptr<SymbolicExpr>& expected) {
    auto delta = test_normalized_delta(actual, expected);
    return delta && delta->is_zero();
}

/// Returns 0 if all tests passed, 1 otherwise.
/// Use as: return TEST_REPORT();
inline int TEST_REPORT() {
    std::cout << "\n===================================================" << std::endl;
    std::cout << "Results: " << g_passes << " passed, " << g_failures << " failed." << std::endl;
    if (g_failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    }
    std::cerr << g_failures << " failure(s) encountered." << std::endl;
    return 1;
}

// Recursive numeric evaluator for symbolic expressions that may contain
// AddNode, MultiplyNode, PowerNode after variable substitution. Returns
// std::nullopt when the expression still contains free variables or when
// a math error occurs (division by zero, etc.).
#include "symbolic_ast.hpp"
#include <optional>
#include <cmath>
#include <limits>


inline std::optional<double> test_numeric_eval(const std::shared_ptr<SymbolicExpr>& e) {
    if (!e || !LMCAS::detail::node(e)) return std::nullopt;
    auto root = LMCAS::detail::node(e);

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) return std::get<lmmc_real_t>(num->value());
        if (std::holds_alternative<BigInt>(num->value())) return (double)std::get<BigInt>(num->value()).to_double();
        if (std::holds_alternative<Rational>(num->value())) return (double)std::get<Rational>(num->value()).to_double();
        return std::nullopt;
    }
    if (std::dynamic_pointer_cast<const VariableNode>(root)) {
        return std::nullopt;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(root)) {
        double s = 0.0;
        for (auto& op : add->operands()) {
            auto v = test_numeric_eval(LMCAS::detail::make_expression_ptr(op));
            if (!v) return std::nullopt;
            s += *v;
        }
        return s;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(root)) {
        double s = 1.0;
        for (auto& op : mul->operands()) {
            auto v = test_numeric_eval(LMCAS::detail::make_expression_ptr(op));
            if (!v) return std::nullopt;
            s *= *v;
        }
        return s;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(root)) {
        auto b = test_numeric_eval(LMCAS::detail::make_expression_ptr(pow->base()));
        auto x = test_numeric_eval(LMCAS::detail::make_expression_ptr(pow->exponent()));
        if (!b || !x) return std::nullopt;
        if (*b == 0.0 && *x < 0.0) return std::nullopt;
        double v = std::pow(*b, *x);
        if (!std::isfinite(v)) return std::nullopt;
        return v;
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(root)) {
        if (func->arguments().size() == 1) {
            auto arg = test_numeric_eval(LMCAS::detail::make_expression_ptr(func->arguments()[0]));
            if (!arg) return std::nullopt;
            double v = std::numeric_limits<double>::quiet_NaN();
            switch (func->type()) {
                case FunctionNode::FuncType::Sin: v = std::sin(*arg); break;
                case FunctionNode::FuncType::Cos: v = std::cos(*arg); break;
                case FunctionNode::FuncType::Tan: v = std::tan(*arg); break;
                case FunctionNode::FuncType::Exp: v = std::exp(*arg); break;
                case FunctionNode::FuncType::Ln:  v = std::log(*arg); break;
                case FunctionNode::FuncType::Sqrt: v = std::sqrt(*arg); break;
                case FunctionNode::FuncType::Abs: v = std::abs(*arg); break;
                // Reciprocal trig: sec=1/cos, csc=1/sin, cot=cos/sin.
                case FunctionNode::FuncType::Sec: {
                    double c = std::cos(*arg);
                    if (c == 0.0) return std::nullopt;
                    v = 1.0 / c;
                    break;
                }
                case FunctionNode::FuncType::Csc: {
                    double s = std::sin(*arg);
                    if (s == 0.0) return std::nullopt;
                    v = 1.0 / s;
                    break;
                }
                case FunctionNode::FuncType::Cot: {
                    double s = std::sin(*arg);
                    if (s == 0.0) return std::nullopt;
                    v = std::cos(*arg) / s;
                    break;
                }
                // Hyperbolic trig.
                case FunctionNode::FuncType::Sinh: v = std::sinh(*arg); break;
                case FunctionNode::FuncType::Cosh: v = std::cosh(*arg); break;
                case FunctionNode::FuncType::Tanh: v = std::tanh(*arg); break;
                // Inverse trig.
                case FunctionNode::FuncType::ArcSin:
                    if (*arg < -1.0 || *arg > 1.0) return std::nullopt;
                    v = std::asin(*arg);
                    break;
                case FunctionNode::FuncType::ArcCos:
                    if (*arg < -1.0 || *arg > 1.0) return std::nullopt;
                    v = std::acos(*arg);
                    break;
                case FunctionNode::FuncType::ArcTan: v = std::atan(*arg); break;
                // Special functions have no elementary closed form; treat as
                // non-numeric for now so round-trip checks fall back to AST equality.
                case FunctionNode::FuncType::Erf:
                case FunctionNode::FuncType::Ei:
                case FunctionNode::FuncType::Si:
                case FunctionNode::FuncType::Ci:
                case FunctionNode::FuncType::Li:
                    return std::nullopt;
                default: return std::nullopt;
            }
            if (!std::isfinite(v)) return std::nullopt;
            return v;
        }
    }
    return std::nullopt;
}
