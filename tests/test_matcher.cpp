#include <iostream>
#include <cassert>
#include "symbolic.hpp"
#include "matcher.hpp"

using namespace lamina;

SymbolicExpr var(const std::string& name) {
    return SymbolicExpr(SymbolicFactory::create_variable(name));
}

SymbolicExpr num(int n) {
    return SymbolicExpr(SymbolicFactory::create_number(BigInt(n)));
}

SymbolicExpr add(SymbolicExpr a, SymbolicExpr b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(a.root);
    ops.push_back(b.root);
    return SymbolicExpr(SymbolicFactory::create_add(ops));
}

SymbolicExpr mul(SymbolicExpr a, SymbolicExpr b) {
    std::vector<std::shared_ptr<SymbolicNode>> ops;
    ops.push_back(a.root);
    ops.push_back(b.root);
    return SymbolicExpr(SymbolicFactory::create_multiply(ops));
}

SymbolicExpr pow(SymbolicExpr b, SymbolicExpr e) {
    return SymbolicExpr(std::make_shared<PowerNode>(b.root, e.root));
}

SymbolicExpr func(FunctionNode::FuncType type, SymbolicExpr arg) {
    std::vector<std::shared_ptr<SymbolicNode>> args;
    args.push_back(arg.root);
    return SymbolicExpr(std::make_shared<FunctionNode>(type, args));
}

void test_basic_match() {
    std::cout << "Testing basic match..." << std::endl;

    auto x = var("x");
    auto one = num(1);
    auto p = add(x, one);

    auto y = var("y");
    auto t = add(y, one);

    MatchMap bindings;
    std::unordered_set<std::string> wildcards = {"x"};

    bool matched = Matcher::match(p, t, wildcards, bindings);
    if (!matched) {
        std::cerr << "Basic match failed!" << std::endl;
        std::cerr << "Pattern: x + 1" << std::endl;
        std::cerr << "Target: y + 1" << std::endl;
    }
    assert(matched);
    assert(bindings.count("x"));

    assert(bindings["x"].root->equals(*y.root));

    std::cout << "Basic match passed." << std::endl;
}

void test_trig_identity() {
    std::cout << "Testing trig identity sin(x)^2 + cos(x)^2 -> 1..." << std::endl;

    auto A = var("A");
    auto two = num(2);
    auto sinA = func(FunctionNode::FuncType::Sin, A);
    auto cosA = func(FunctionNode::FuncType::Cos, A);
    auto term1 = pow(sinA, two);
    auto term2 = pow(cosA, two);
    auto pattern = add(term1, term2);

    auto y = var("y");
    auto sinY = func(FunctionNode::FuncType::Sin, y);
    auto cosY = func(FunctionNode::FuncType::Cos, y);
    auto target = add(pow(sinY, two), pow(cosY, two));

    MatchMap bindings;
    std::unordered_set<std::string> wildcards = {"A"};

    bool matched = Matcher::match(pattern, target, wildcards, bindings);
    if (!matched) {
        std::cerr << "Trig identity match failed!" << std::endl;
    }
    assert(matched);
    assert(bindings.count("A"));
    assert(bindings["A"].root->equals(*y.root));

    auto replacement = num(1);
    auto result = Matcher::replace(replacement, bindings);
    assert(result.root->is_one());

    std::cout << "Trig identity passed." << std::endl;
}

void test_rewrite_engine() {
    std::cout << "Testing rewrite engine..." << std::endl;

    RewriteEngine engine;

    auto x = var("x");
    auto zero = num(0);
    auto pattern = add(x, zero);
    auto replacement = x;

    engine.add_rule(Rule(pattern, replacement, {"x"}));

    auto y = var("y");
    auto one = num(1);
    auto inner = add(y, one);
    auto target = add(inner, zero);

    auto result = engine.apply(target);

    assert(result.root->equals(*inner.root));

    std::cout << "Rewrite engine passed." << std::endl;
}

int main() {
    try {
        test_basic_match();
        test_trig_identity();
        test_rewrite_engine();
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
