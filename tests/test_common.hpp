#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "symbolic.hpp"

inline int g_failures = 0;
inline int g_passes = 0;

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

inline std::optional<double> test_numeric_eval(const std::shared_ptr<SymbolicExpr>& e) {
    if (!e || !e->root) return 0.0;
    auto root = e->root;

    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value)) return std::get<lmmc_real_t>(num->value);
        if (std::holds_alternative<BigInt>(num->value)) return (double)std::get<BigInt>(num->value).to_double();
        if (std::holds_alternative<Rational>(num->value)) return (double)std::get<Rational>(num->value).to_double();
        return 0.0;
    }
    if (std::dynamic_pointer_cast<VariableNode>(root)) {
        return std::nullopt;
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(root)) {
        double s = 0.0;
        for (auto& op : add->operands) {
            auto v = test_numeric_eval(std::make_shared<SymbolicExpr>(op));
            if (!v) return std::nullopt;
            s += *v;
        }
        return s;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(root)) {
        double s = 1.0;
        for (auto& op : mul->operands) {
            auto v = test_numeric_eval(std::make_shared<SymbolicExpr>(op));
            if (!v) return std::nullopt;
            s *= *v;
        }
        return s;
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(root)) {
        auto b = test_numeric_eval(std::make_shared<SymbolicExpr>(pow->base));
        auto x = test_numeric_eval(std::make_shared<SymbolicExpr>(pow->exponent));
        if (!b || !x) return std::nullopt;
        if (*b == 0.0 && *x < 0.0) return std::nullopt;
        double v = std::pow(*b, *x);
        if (!std::isfinite(v)) return std::nullopt;
        return v;
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(root)) {
        if (func->arguments.size() == 1) {
            auto arg = test_numeric_eval(std::make_shared<SymbolicExpr>(func->arguments[0]));
            if (!arg) return std::nullopt;
            switch (func->type) {
                case FunctionNode::FuncType::Sin: return std::sin(*arg);
                case FunctionNode::FuncType::Cos: return std::cos(*arg);
                case FunctionNode::FuncType::Tan: return std::tan(*arg);
                case FunctionNode::FuncType::Exp: return std::exp(*arg);
                case FunctionNode::FuncType::Ln:  return std::log(*arg);
                case FunctionNode::FuncType::Sqrt: return std::sqrt(*arg);
                case FunctionNode::FuncType::Abs: return std::abs(*arg);
                default: break;
            }
        }
    }
    return std::nullopt;
}
