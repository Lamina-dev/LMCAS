#include "symbolic.hpp"
#include "matcher.hpp"
#include "test_common.hpp"
#include <vector>
#include <string>

using namespace LMCAS;

std::shared_ptr<SymbolicExpr> mk_wildcard(const std::string &name) {
    return LMCAS::detail::make_expression_ptr(wildcard(name));
}

TEST(LmcasMatcherComm, CommutativeAdd) {
    auto x = mk_wildcard("x");
    auto y = mk_wildcard("y");
    auto pattern = SymbolicExpr::add(x, y);

    auto n2 = SymbolicExpr::number(2);
    auto n3 = SymbolicExpr::number(3);
    auto target = SymbolicExpr::add(n3, n2);

    MatchMap results;
    std::unordered_set<std::string> w = {"x", "y"};

    ASSERT_TRUE(Matcher::match(*pattern, *target, w, results));
    EXPECT_EQ(results.count("x"), 1u);
    EXPECT_EQ(results.count("y"), 1u);
}

TEST(LmcasMatcherComm, CommutativeMul) {
    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B");
    auto pattern = SymbolicExpr::multiply(A, B);

    auto vx = SymbolicExpr::variable("x");
    auto vy = SymbolicExpr::variable("y");
    auto target = SymbolicExpr::multiply(vy, vx);

    MatchMap results;
    std::unordered_set<std::string> w = {"A", "B"};

    ASSERT_TRUE(Matcher::match(*pattern, *target, w, results));
    EXPECT_EQ(results.count("A"), 1u);
    EXPECT_EQ(results.count("B"), 1u);
}

TEST(LmcasMatcherComm, CommutativeNested) {
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

    ASSERT_TRUE(Matcher::match(*pattern, *target, w, results));
    EXPECT_EQ(results.at("x").to_string(), "y");
}

TEST(LmcasMatcherComm, SubsetMatch) {
    auto A = mk_wildcard("A");
    auto B = mk_wildcard("B");
    auto pattern = SymbolicExpr::add(A, B);

    auto va = SymbolicExpr::variable("a");
    auto vb = SymbolicExpr::variable("b");
    auto vc = SymbolicExpr::variable("c");
    auto vd = SymbolicExpr::variable("d");

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    ops.push_back(LMCAS::detail::node(va));
    ops.push_back(LMCAS::detail::node(vb));
    ops.push_back(LMCAS::detail::node(vc));
    ops.push_back(LMCAS::detail::node(vd));
    auto target = LMCAS::detail::expression_from_node(SymbolicFactory::create_add(ops));

    MatchMap results;
    std::unordered_set<std::string> w = {"A", "B"};

    ASSERT_TRUE(Matcher::match(*pattern, target, w, results));
    EXPECT_NE(results.find("__Add_REST__"), results.end());
}
