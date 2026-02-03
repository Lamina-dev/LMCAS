#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include "../symbolic.hpp"

inline int g_failures = 0;

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
    }
}

inline void EXPECT_CONTAINS(const std::string& actual, const std::vector<std::string>& tokens, const std::string& msg) {
    bool ok = true;
    for(const auto& t : tokens) {
        if (actual.find(t) == std::string::npos) {
            std::cerr << "[FAIL] " << msg << " | Missing token: " << t << "\n  In: " << actual << std::endl;
            g_failures++;
            return;
        }
    }
    std::cout << "[PASS] " << msg << std::endl;
}

inline int TEST_REPORT() {
    if (g_failures == 0) {
        std::cout << "\n===================================================" << std::endl;
        std::cout << "All tests passed!" << std::endl;
        return 0;
    }
    std::cerr << "\n===================================================" << std::endl;
    std::cerr << g_failures << " failures encountered." << std::endl;
    return 1;
}
