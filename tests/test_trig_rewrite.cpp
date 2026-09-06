#include "../include/symbolic.hpp"
#include "../include/matcher.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "test_common.hpp"
#include <iostream>
#include <cassert>

using namespace LMCAS;

bool test_sin_sq_add_cos_sq() {
    std::cout << "Testing sin(x)^2 + cos(x)^2 -> 1" << std::endl;
    auto x = SymbolicExpr::variable("x");
    auto sinx = SymbolicExpr::sin(x);
    auto cosx = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sinx, SymbolicExpr::number(2));
    auto cos2 = SymbolicExpr::power(cosx, SymbolicExpr::number(2));

    auto expr = SymbolicExpr::add(sin2, cos2);
    auto simplified = expr->simplify_trig();

    std::cout << "Original: " << expr->to_string() << std::endl;
    std::cout << "Simplified: " << simplified->to_string() << std::endl;

    return simplified->to_string() == "1";
}

bool test_sin_2x() {
    std::cout << "Testing sin(2*x) -> 2*sin(x)*cos(x)" << std::endl;
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto sin2x = SymbolicExpr::sin(two_x);

    auto simplified = sin2x->simplify_trig();

    std::cout << "Original: " << sin2x->to_string() << std::endl;
    std::cout << "Simplified: " << simplified->to_string() << std::endl;

    return simplified->to_string() == "2*sin(x)*cos(x)" ||
           simplified->to_string() == "2*cos(x)*sin(x)" ||
           simplified->to_string() == "(sin(x) * cos(x) * 2)";
}

bool test_cos_2x() {
    std::cout << "Testing cos(2*x) -> cos(x)^2 - sin(x)^2" << std::endl;
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto cos2x = SymbolicExpr::cos(two_x);

    auto simplified = cos2x->simplify_trig();

    std::cout << "Original: " << cos2x->to_string() << std::endl;
    std::cout << "Simplified: " << simplified->to_string() << std::endl;

    std::string s = simplified->to_string();
    std::cout << "DEBUG: " << s << std::endl;

    bool result = (s.find("cos(x)^2") != std::string::npos || s.find("cos(x) ^ 2") != std::string::npos) &&
           (s.find("sin(x)^2") != std::string::npos || s.find("sin(x) ^ 2") != std::string::npos) &&
           (s.find("-1") != std::string::npos || s.find("- 1") != std::string::npos || s.find("- sin") != std::string::npos || s.find("-sin") != std::string::npos);
    if (!result) {

        if (s == "cos(x)^2 + -1*sin(x)^2" || s == "cos(x)^2 - 1*(sin(x)^2)") result = true;
    }
    return result;
}

int main() {
    std::cout << "Starting tests..." << std::endl;
    EXPECT_TRUE(test_sin_sq_add_cos_sq(), "sin^2 + cos^2 simplifies to 1");
    EXPECT_TRUE(test_sin_2x(), "sin(2x) rewrites to 2*sin(x)*cos(x)");
    EXPECT_TRUE(test_cos_2x(), "cos(2x) rewrites to cos(x)^2 - sin(x)^2");

    return TEST_REPORT();
}
