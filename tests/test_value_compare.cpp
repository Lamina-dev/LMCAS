#include "test_common.hpp"
#include <iostream>
#include <set>
#include <map>
#include "value.hpp"
#include "symbolic.hpp"

int main() {
    Value v1(10);
    Value v2(20);
    Value v3(10);
    Value v4 = nullptr;
    Value v5 = nullptr;
    Value v6(3.14);

    EXPECT_TRUE(v1 == v3, "v1 should equal v3");
    EXPECT_TRUE(!(v1 == v2), "v1 should not equal v2");
    EXPECT_TRUE(v1 < v2, "v1 should be less than v2");
    EXPECT_TRUE(!(v2 < v1), "v2 should not be less than v1");
    EXPECT_TRUE(v4 == v5, "null should equal null");
    EXPECT_TRUE(!(v4 < v5), "null should not be less than null");

    EXPECT_TRUE(v1 < v6 || v6 < v1, "Different types should be comparable");

    auto x = SymbolicExpr::variable("x");
    auto expr1 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    auto y = SymbolicExpr::variable("x");
    auto expr2 = SymbolicExpr::add(y, SymbolicExpr::number(1));

    Value sym1(expr1);
    Value sym2(expr2);

    EXPECT_TRUE(sym1 == sym2, "SymbolicExpr AST structure should compare equal");
    EXPECT_TRUE(!(sym1 < sym2) && !(sym2 < sym1), "SymbolicExpr AST structure < should be false for equal trees");

    std::set<Value> s;
    s.insert(v1);
    s.insert(v2);
    s.insert(v3);
    s.insert(v4);
    s.insert(v5);
    s.insert(v6);
    s.insert(sym1);
    s.insert(sym2);

    EXPECT_TRUE(s.size() == 5, "Set should have 5 unique elements (10, 20, null, 3.14, x+1)");

    return TEST_REPORT();
}
