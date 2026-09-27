#include "expr.hpp"
#include "test_common.hpp"
#include <cmath>
#include <vector>

using namespace LMCAS;

TEST(Print, ArithmeticPreservesPowerGroupingWhenParsed) {
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
             SymbolicExpr::power(x, two), SymbolicExpr::number(9))),
         5},
        {SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3))), std::cbrt(4.0)},
        {SymbolicExpr::power(x, SymbolicExpr::number(Rational(-1, 2))), 0.5},
        {SymbolicExpr::power(SymbolicExpr::number(Rational(3, 2)), x), 5.0625},
        {SymbolicExpr::power(SymbolicExpr::number(-2), x), 16},
        {SymbolicExpr::power(SymbolicExpr::number(-2.0), x), 16},
        {SymbolicExpr::power(SymbolicExpr::power(x, two), three), 4096},
        {SymbolicExpr::power(x, SymbolicExpr::power(two, three)), 65536},
        {SymbolicExpr::multiply(x, SymbolicExpr::number(Rational(3, 2))), 6},
        {SymbolicExpr::multiply(x, SymbolicExpr::multiply(two, three)), 24},
        {SymbolicExpr::power(SymbolicExpr::number(1.25), two), 1.5625},
    };
    for (const auto &test : cases) {
        SCOPED_TRACE(test.expression->to_string());
        auto parsed = LMCAS::parse_expr(test.expression->to_string());
        EXPECT_TRUE(parsed.has_value()) << "printed arithmetic can be parsed";
        if (!parsed)
            continue;
        auto value = LMCAS::evalf(*parsed.value(), {{"x", 4}});
        EXPECT_TRUE(value.has_value()) << "parsed arithmetic has a finite value";
        if (value) {
            EXPECT_TRUE(std::isfinite(value.value().value));
            EXPECT_TRUE(std::isfinite(test.expected));
            EXPECT_NEAR(value.value().value, test.expected, 1e-12) << "printed arithmetic preserves its numeric meaning";
        }
    }
}

TEST(Print, RelationalAndMembershipTreesRetainGrouping) {
    for (const char *source : {"(x < 1) < 2", "x < (1 < 2)",
                               "(x == 1) == 2", "x == (1 == 2)",
                               "(x in {1}) in {2}", "not (x < 1 and x > 0)"}) {
        SCOPED_TRACE(source);
        auto original = parse_expr(source);
        EXPECT_TRUE(original.has_value()) << "nested relation syntax parses";
        if (!original)
            continue;
        auto parsed = parse_expr(original.value()->to_string());
        EXPECT_TRUE(parsed && structurally_equal(*original.value(), *parsed.value())) << "printing preserves the relational expression tree";
    }
}
