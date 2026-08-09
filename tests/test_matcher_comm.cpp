#include "../include/symbolic.hpp"
#include "../include/matcher.hpp"
#include "test_common.hpp"
#include <iostream>
#include <vector>
#include <string>

using namespace lamina;

std::shared_ptr<SymbolicExpr> mk_wildcard(const std::string& name) {
    return lamina::detail::make_expression_ptr(wildcard(name));
}

bool test_commutative_add() {
    std::cout << "Testing commutative addition..." << std::endl;

    auto x = mk_wildcard("x");
    auto y = mk_wildcard("y");
    auto pattern = SymbolicExpr::add(x, y);

    auto n2 = SymbolicExpr::number(2);
    auto n3 = SymbolicExpr::number(3);
    auto target = SymbolicExpr::add(n3, n2);

    MatchMap results;
    std::unordered_set<std::string> w = {"x", "y"};

    if (!Matcher::match(*pattern, *target, w, results)) {
        std::cout << "Failed to match x+y with 3+2" << std::endl;
        return false;
    }

    std::cout << "Matched: x=" << results.at("x").to_string() << ", y=" << results.at("y").to_string() << std::endl;
    return true;
}

bool test_commutative_mul() {
    std::cout << "Testing commutative multiplication..." << std::endl;

    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B");
    auto pattern = SymbolicExpr::multiply(A, B);

    auto vx = SymbolicExpr::variable("x");
    auto vy = SymbolicExpr::variable("y");
    auto target = SymbolicExpr::multiply(vy, vx);

    MatchMap results;
    std::unordered_set<std::string> w = {"A", "B"};

    if (!Matcher::match(*pattern, *target, w, results)) {
        std::cout << "Failed to match A*B with y*x" << std::endl;
        return false;
    }
     std::cout << "Matched: A=" << results.at("A").to_string() << ", B=" << results.at("B").to_string() << std::endl;
    return true;
}

bool test_commutative_nested() {
    std::cout << "Testing nested commutative..." << std::endl;

    auto x = mk_wildcard("x");
    auto sinx = SymbolicExpr::sin(x);
    auto cosx = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sinx, SymbolicExpr::number(2));
    auto cos2 = SymbolicExpr::power(cosx, SymbolicExpr::number(2));
    auto pattern = SymbolicExpr::add(sin2, cos2);

    auto y = SymbolicExpr::variable("y");
    auto siny = SymbolicExpr::sin(y);
    auto cosy = SymbolicExpr::cos(y);
    auto sin2y = SymbolicExpr::power(siny, SymbolicExpr::number(2));
    auto cos2y = SymbolicExpr::power(cosy, SymbolicExpr::number(2));

    auto target = SymbolicExpr::add(cos2y, sin2y);

    MatchMap results;
    std::unordered_set<std::string> w = {"x"};

    if (!Matcher::match(*pattern, *target, w, results)) {
         std::cout << "Failed to match sin(x)^2 + cos(x)^2 with cos(y)^2 + sin(y)^2" << std::endl;
         return false;
    }

    std::cout << "Matched: x=" << results.at("x").to_string() << std::endl;
    return results.at("x").to_string() == "y";
}

bool test_subset_match();

int main() {
    std::cerr << "Starting matcher tests..." << std::endl;
    EXPECT_TRUE(test_commutative_add(), "commutative addition matcher");
    EXPECT_TRUE(test_commutative_mul(), "commutative multiplication matcher");
    EXPECT_TRUE(test_commutative_nested(), "nested commutative matcher");
    EXPECT_TRUE(test_subset_match(), "subset AddNode matcher");

    return TEST_REPORT();
}

bool test_subset_match() {
    std::cout << "Testing subset matching (AddNode)..." << std::endl;
    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B");
    auto pattern = SymbolicExpr::add(A, B);

    auto va = SymbolicExpr::variable("a");
    auto vb = SymbolicExpr::variable("b");
    auto vc = SymbolicExpr::variable("c");
    auto vd = SymbolicExpr::variable("d");

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(lamina::detail::node(va));
    ops.push_back(lamina::detail::node(vb));
    ops.push_back(lamina::detail::node(vc));
    ops.push_back(lamina::detail::node(vd));
    auto target = lamina::detail::expression_from_node(SymbolicFactory::create_add(ops));

    MatchMap results;
    std::unordered_set<std::string> w = {"A", "B"};

    if (!Matcher::match(*pattern, target, w, results)) {
        std::cout << "Match failed for subset (A+B inside a+b+c+d)" << std::endl;
        return false;
    }

    if (results.find("__Add_REST__") == results.end()) {
        std::cout << "Missing __Add_REST__ binding" << std::endl;
        return false;
    }

    std::cout << "Remainder: " << results.at("__Add_REST__").to_string() << std::endl;
    return true;
}
