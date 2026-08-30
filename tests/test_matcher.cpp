#include <iostream>
#include "symbolic.hpp"
#include "matcher.hpp"
#include "test_common.hpp"

using namespace lamina;

SymbolicExpr var(const std::string& name) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

SymbolicExpr num(int n) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(n)));
}

SymbolicExpr add(SymbolicExpr a, SymbolicExpr b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(lamina::detail::node(a));
    ops.push_back(lamina::detail::node(b));
    return lamina::detail::expression_from_node(SymbolicFactory::create_add(ops));
}

SymbolicExpr mul(SymbolicExpr a, SymbolicExpr b) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(lamina::detail::node(a));
    ops.push_back(lamina::detail::node(b));
    return lamina::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
}

SymbolicExpr pow(SymbolicExpr b, SymbolicExpr e) {
    return lamina::detail::expression_from_node(lamina::detail::make_node<PowerNode>(lamina::detail::node(b), lamina::detail::node(e)));
}

SymbolicExpr func(FunctionNode::FuncType type, SymbolicExpr arg) {
    std::vector<std::shared_ptr<const SymbolicNode>> args;
    args.push_back(lamina::detail::node(arg));
    return lamina::detail::expression_from_node(lamina::detail::make_node<FunctionNode>(type, args));
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
    EXPECT_TRUE(matched, "basic wildcard pattern matches target");
    EXPECT_TRUE(bindings.count("x") == 1, "basic match binds wildcard x");
    EXPECT_TRUE(bindings.count("x") == 1 &&
                    lamina::detail::node(bindings.at("x"))->equals(*lamina::detail::node(y)),
                "basic match binds x to y");

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
    EXPECT_TRUE(matched, "trig identity pattern matches target");
    EXPECT_TRUE(bindings.count("A") == 1, "trig identity binds wildcard A");
    EXPECT_TRUE(bindings.count("A") == 1 &&
                    lamina::detail::node(bindings.at("A"))->equals(*lamina::detail::node(y)),
                "trig identity binds A to y");

    auto replacement = num(1);
    auto result = Matcher::replace(replacement, bindings);
    EXPECT_TRUE(lamina::detail::node(result)->is_one(), "trig identity replacement returns one");

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
    EXPECT_TRUE(result.has_value(), "rewrite engine application succeeds");
    if (!result) return;

    EXPECT_TRUE(lamina::detail::node(result.value())->equals(
                    *lamina::detail::node(inner)),
                "rewrite engine removes additive zero");

    std::cout << "Rewrite engine passed." << std::endl;
}

int main() {
    try {
        test_basic_match();
        test_trig_identity();
        test_rewrite_engine();
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }

    std::cout << "All tests passed!" << std::endl;
    return TEST_REPORT();
}
