/**
 * @file test_squeeze_theorem.cpp
 * @brief 夹逼定理极限测试：验证 LimitVisitor 的 try_squeeze() 功能。
 *
 * 测试用例覆盖：
 * - x·sin(1/x) as x→0 = 0 (基本夹逼)
 * - x²·cos(1/x) as x→0 = 0 (高次零因子)
 * - x·cos(x) as x→0 = 0 (bounded × zero)
 * - sin(1/x)/x as x→∞ = 0 (无穷处夹逼)
 * - 一般夹逼定理：f = g + bounded×zero → lim f = lim g
 * - arctan 作为有界函数的识别
 * - 非夹逼情况不误触发
 */
#include "test_common.hpp"
#include "symbolic_ast.hpp"
#include "visitors/limit_visitor.hpp"
#include "visitors/differentiation_visitor.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto five = SymbolicExpr::number(5);
    auto neg_one = SymbolicExpr::number(-1);

    // --- Test 1: lim(x→0) x·sin(1/x) = 0 ---
    TEST_CASE("Squeeze: x*sin(1/x) as x->0");
    {
        // Build: x * sin(1/x) = x * sin(x^(-1))
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto sin_inv_x = SymbolicExpr::sin(inv_x);
        auto expr = SymbolicExpr::multiply(x, sin_inv_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x*sin(1/x), x->0) = 0");
    }

    // --- Test 2: lim(x→0) x^2·cos(1/x) = 0 ---
    TEST_CASE("Squeeze: x^2*cos(1/x) as x->0");
    {
        auto x_sq = SymbolicExpr::power(x, two);
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto cos_inv_x = SymbolicExpr::cos(inv_x);
        auto expr = SymbolicExpr::multiply(x_sq, cos_inv_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x^2*cos(1/x), x->0) = 0");
    }

    // --- Test 3: lim(x→0) x·cos(x) = 0 ---
    TEST_CASE("Squeeze: x*cos(x) as x->0");
    {
        auto cos_x = SymbolicExpr::cos(x);
        auto expr = SymbolicExpr::multiply(x, cos_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x*cos(x), x->0) = 0");
    }

    // --- Test 4: lim(x→0) x·sin(x) = 0 ---
    TEST_CASE("Squeeze: x*sin(x) as x->0");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto expr = SymbolicExpr::multiply(x, sin_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x*sin(x), x->0) = 0");
    }

    // --- Test 5: lim(x→0) sin(x) should NOT be affected (no zero-tending factor) ---
    TEST_CASE("Non-squeeze: sin(x) as x->0 = 0 (direct substitution)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto lim = sin_x->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(sin(x), x->0) = 0");
    }

    // --- Test 6: lim(x→2) x·sin(x) should NOT use squeeze (x doesn't tend to 0) ---
    TEST_CASE("Non-squeeze: x*sin(x) as x->2 (normal substitution)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto expr = SymbolicExpr::multiply(x, sin_x);
        auto lim = expr->limit("x", two);
        // Should be 2*sin(2), not 0
        EXPECT_TRUE(lim != nullptr, "limit(x*sin(x), x->2) is not null");
        auto lim_str = lim->to_string();
        EXPECT_TRUE(lim_str != "0", "limit(x*sin(x), x->2) != 0");
    }

    // --- Test 7: lim(x→0) x·sin(1/x)·cos(1/x) = 0 (multiple bounded factors) ---
    TEST_CASE("Squeeze: x*sin(1/x)*cos(1/x) as x->0");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto sin_inv_x = SymbolicExpr::sin(inv_x);
        auto cos_inv_x = SymbolicExpr::cos(inv_x);
        auto expr = SymbolicExpr::multiply(SymbolicExpr::multiply(x, sin_inv_x), cos_inv_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x*sin(1/x)*cos(1/x), x->0) = 0");
    }

    // --- Test 8: lim(x→0) x·arctan(1/x) = 0 (arctan is bounded) ---
    TEST_CASE("Squeeze: x*arctan(1/x) as x->0");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        // Create arctan(1/x) using FunctionNode directly
        auto atan_inv_x = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcTan,
                std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(inv_x)}));
        auto expr = SymbolicExpr::multiply(x, atan_inv_x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(x*arctan(1/x), x->0) = 0");
    }

    // --- Test 9: lim(x→0) (5 + x·sin(1/x)) = 5 (general squeeze: constant + squeeze-to-zero) ---
    TEST_CASE("General squeeze: 5 + x*sin(1/x) as x->0 = 5");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto sin_inv_x = SymbolicExpr::sin(inv_x);
        auto squeeze_term = SymbolicExpr::multiply(x, sin_inv_x);
        auto expr = SymbolicExpr::add(five, squeeze_term);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "5", "limit(5 + x*sin(1/x), x->0) = 5");
    }

    // --- Test 10: lim(x→0) (3 + x^2·cos(1/x)) = 3 ---
    TEST_CASE("General squeeze: 3 + x^2*cos(1/x) as x->0 = 3");
    {
        auto x_sq = SymbolicExpr::power(x, two);
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto cos_inv_x = SymbolicExpr::cos(inv_x);
        auto squeeze_term = SymbolicExpr::multiply(x_sq, cos_inv_x);
        auto expr = SymbolicExpr::add(three, squeeze_term);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "3", "limit(3 + x^2*cos(1/x), x->0) = 3");
    }

    // --- Test 11: lim(x→0) sin(1/x)^2 * x = 0 (bounded^2 × zero) ---
    TEST_CASE("Squeeze: sin(1/x)^2 * x as x->0");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto sin_inv_x = SymbolicExpr::sin(inv_x);
        auto sin_sq = SymbolicExpr::power(sin_inv_x, two);
        auto expr = SymbolicExpr::multiply(sin_sq, x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit(sin(1/x)^2 * x, x->0) = 0");
    }

    // --- Test 12: lim(x→0) (2+sin(1/x))*x = 0 (bounded_expression × zero) ---
    TEST_CASE("Squeeze: (2+sin(1/x))*x as x->0");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto sin_inv_x = SymbolicExpr::sin(inv_x);
        auto bounded_expr = SymbolicExpr::add(two, sin_inv_x);
        auto expr = SymbolicExpr::multiply(bounded_expr, x);

        auto lim = expr->limit("x", zero);
        EXPECT_EQ_EXPR_STR(lim, "0", "limit((2+sin(1/x))*x, x->0) = 0");
    }

    return TEST_REPORT();
}
