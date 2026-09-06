#include "../include/expr.hpp"
#include "test_common.hpp"
#include <cmath>
#include <vector>

using namespace LMCAS;

int main() {
    TEST_CASE("Printed arithmetic preserves power grouping when parsed");
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto sum = SymbolicExpr::add(x, SymbolicExpr::number(1));
    struct Case {
        LMCAS::ExprPtr expression;
        double expected;
    };
    const std::vector<Case> cases = {
        {sum, 5},
        {SymbolicExpr::multiply(sum, two), 10},
        {SymbolicExpr::sqrt(SymbolicExpr::add(
             SymbolicExpr::power(x, two), SymbolicExpr::number(9))), 5},
        {SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3))), std::cbrt(4.0)},
        {SymbolicExpr::power(x, SymbolicExpr::number(Rational(-1, 2))), 0.5},
        {SymbolicExpr::power(SymbolicExpr::number(Rational(3, 2)), x), 5.0625},
        {SymbolicExpr::power(SymbolicExpr::number(-2), x), 16},
        {SymbolicExpr::power(SymbolicExpr::number(-2.0), x), 16},
        {SymbolicExpr::power(SymbolicExpr::power(x, two), three), 4096},
        {SymbolicExpr::power(x, SymbolicExpr::power(two, three)), 65536},
    };
    for (const auto& test : cases) {
        auto parsed = LMCAS::parse_expr(test.expression->to_string());
        EXPECT_TRUE(parsed.has_value(), "printed arithmetic can be parsed");
        if (!parsed) continue;
        auto value = LMCAS::evalf(*parsed.value(), {{"x", 4}});
        EXPECT_TRUE(value.has_value(), "parsed arithmetic has a finite value");
        if (value) {
            EXPECT_TRUE(std::abs(value.value().value - test.expected) < 1e-12,
                        "printed arithmetic preserves its numeric meaning");
        }
    }
    return TEST_REPORT();
}
