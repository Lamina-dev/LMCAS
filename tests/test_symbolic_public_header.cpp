#include "symbolic.hpp"
#include <gtest/gtest.h>

using namespace LMCAS;

#ifdef LMCAS_INTERNAL_AST_INCLUDED
#error "symbolic.hpp must not include the AST implementation"
#endif

static_assert(sizeof(SymbolicExpr) > 0, "SymbolicExpr must be a complete value type");

TEST(LmcasSymbolicPublicHeader, ConstructAndSubstitute) {
    auto one = SymbolicExpr::number(1);
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(x, one);
    ASSERT_NE(one, nullptr);
    ASSERT_NE(x, nullptr);
    ASSERT_NE(expression, nullptr);
    EXPECT_TRUE(one->is_number());
    EXPECT_EQ(one->get_int(), 1);
    auto substituted = expression->substitute("x", SymbolicExpr::number(2));
    ASSERT_NE(substituted, nullptr);
    auto result = substituted->simplify();
    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(result->is_number());
    EXPECT_EQ(result->get_int(), 3);
}
