#include "../include/symbolic.hpp"
#include "../include/matcher.hpp"
#include <iostream>
#include <vector>
#include <string>

using namespace lamina;

std::shared_ptr<SymbolicExpr> mk_wildcard(const std::string& name) {
    return std::make_shared<SymbolicExpr>(wildcard(name));
}

bool test_commutative_add() {
    std::cout << "Testing commutative addition..." << std::endl;
    // Pattern: x + y
    auto x = mk_wildcard("x");
    auto y = mk_wildcard("y");
    auto pattern = SymbolicExpr::add(x, y);

    // Target: 2 + 3
    auto n2 = SymbolicExpr::number(2);
    auto n3 = SymbolicExpr::number(3);
    auto target = SymbolicExpr::add(n3, n2); // Order reversed relative to "2+3" if we imagine pattern is "first + second"

    MatchMap results;
    std::unordered_set<std::string> w = {"x", "y"};
    
    // Normal match should work (2 matches x, 3 matches y Or 3 matches x, 2 matches y)
    if (!Matcher::match(*pattern, *target, w, results)) {
        std::cout << "Failed to match x+y with 3+2" << std::endl;
        return false;
    }
    
    std::cout << "Matched: x=" << results["x"].to_string() << ", y=" << results["y"].to_string() << std::endl;
    return true;
}

bool test_commutative_mul() {
    std::cout << "Testing commutative multiplication..." << std::endl;
    // Pattern: A * B
    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B"); // Note: "B" is wildcard name
    auto pattern = SymbolicExpr::multiply(A, B);

    // Target: y * x
    auto vx = SymbolicExpr::variable("x");
    auto vy = SymbolicExpr::variable("y");
    auto target = SymbolicExpr::multiply(vy, vx);

    MatchMap results;
    std::unordered_set<std::string> w = {"A", "B"};
    
    if (!Matcher::match(*pattern, *target, w, results)) {
        std::cout << "Failed to match A*B with y*x" << std::endl;
        return false;
    }
     std::cout << "Matched: A=" << results["A"].to_string() << ", B=" << results["B"].to_string() << std::endl;
    return true;
}

bool test_commutative_nested() {
    std::cout << "Testing nested commutative..." << std::endl;
    // Pattern: sin(x)^2 + cos(x)^2 -> 1
    // The pattern is Add(Pow(sin(x), 2), Pow(cos(x), 2))
    // We want to match Add(Pow(cos(y), 2), Pow(sin(y), 2))
    
    auto x = mk_wildcard("x");
    auto sinx = SymbolicExpr::sin(x);
    auto cosx = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sinx, SymbolicExpr::number(2));
    auto cos2 = SymbolicExpr::power(cosx, SymbolicExpr::number(2));
    auto pattern = SymbolicExpr::add(sin2, cos2);
    
    // Target
    auto y = SymbolicExpr::variable("y");
    auto siny = SymbolicExpr::sin(y);
    auto cosy = SymbolicExpr::cos(y);
    auto sin2y = SymbolicExpr::power(siny, SymbolicExpr::number(2));
    auto cos2y = SymbolicExpr::power(cosy, SymbolicExpr::number(2));
    // Target is cos^2 + sin^2 (reversed order)
    auto target = SymbolicExpr::add(cos2y, sin2y);
    
    MatchMap results;
    std::unordered_set<std::string> w = {"x"};
    
    if (!Matcher::match(*pattern, *target, w, results)) {
         std::cout << "Failed to match sin(x)^2 + cos(x)^2 with cos(y)^2 + sin(y)^2" << std::endl;
         return false;
    }
    
    std::cout << "Matched: x=" << results["x"].to_string() << std::endl;
    return results["x"].to_string() == "y";
}


bool test_subset_match();

int main() {
    std::cerr << "Starting matcher tests..." << std::endl;
    bool pass = true;
    if (!test_commutative_add()) {
        std::cerr << "test_commutative_add failed" << std::endl;
        pass = false;
    }
    if (!test_commutative_mul()) {
        std::cerr << "test_commutative_mul failed" << std::endl;
        pass = false;
    }
    if (!test_commutative_nested()) {
        std::cerr << "test_commutative_nested failed" << std::endl;
        pass = false;
    }
    if (!test_subset_match()) {
        std::cerr << "test_subset_match failed" << std::endl;
        pass = false;
    }
    
    if (pass) std::cerr << "All matcher tests passed!" << std::endl;
    else std::cerr << "Some matcher tests failed." << std::endl;
    
    return pass ? 0 : 1;
}

bool test_subset_match() {
    std::cout << "Testing subset matching (AddNode)..." << std::endl;
    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B");
    auto pattern = SymbolicExpr::add(A, B);
    
    // Target: a + b + c + d
    auto va = SymbolicExpr::variable("a");
    auto vb = SymbolicExpr::variable("b");
    auto vc = SymbolicExpr::variable("c");
    auto vd = SymbolicExpr::variable("d");
    
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(va->root);
    ops.push_back(vb->root);
    ops.push_back(vc->root);
    ops.push_back(vd->root);
    auto target = SymbolicExpr(SymbolicFactory::create_add(ops));
    
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
    
    std::cout << "Remainder: " << results["__Add_REST__"].to_string() << std::endl;
    return true;
}

