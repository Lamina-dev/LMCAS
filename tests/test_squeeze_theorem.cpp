#include "limit_result.hpp"
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
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "internal/visitors/differentiation_visitor.hpp"

using namespace LMCAS;

TEST(SqueezeTheorem, SineOscillationSqueeze) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto expr = SymbolicExpr::multiply(x, sin_inv_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x*sin(1/x), x->0) = 0";
}

TEST(SqueezeTheorem, QuadraticCosineSqueeze) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto two = SymbolicExpr::number(2);
    auto neg_one = SymbolicExpr::number(-1);

    auto x_sq = SymbolicExpr::power(x, two);
    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto cos_inv_x = SymbolicExpr::cos(inv_x);
    auto expr = SymbolicExpr::multiply(x_sq, cos_inv_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x^2*cos(1/x), x->0) = 0";
}

TEST(SqueezeTheorem, CosineProductZero) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);

    auto cos_x = SymbolicExpr::cos(x);
    auto expr = SymbolicExpr::multiply(x, cos_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x*cos(x), x->0) = 0";
}

TEST(SqueezeTheorem, SineProductZero) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);

    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::multiply(x, sin_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x*sin(x), x->0) = 0";
}

TEST(SqueezeTheorem, DirectSineLimit) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);

    auto sin_x = SymbolicExpr::sin(x);
    auto lim = LMCAS::limit_expression_checked(sin_x, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(sin(x), x->0) = 0";
}

TEST(SqueezeTheorem, NonzeroProductLimit) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);

    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::multiply(x, sin_x);
    auto result = LMCAS::limit_expression_checked(expr, "x", two);
    ASSERT_TRUE(result.has_value()) << result.error().message;

    auto lim = result.value();
    ASSERT_TRUE(lim) << "limit(x*sin(x), x->2) is not null";
    auto expected = SymbolicExpr::multiply(two, SymbolicExpr::sin(two));
    EXPECT_TRUE(test_proved_equivalent(lim, expected));
}

TEST(SqueezeTheorem, MultipleBoundedFactors) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto cos_inv_x = SymbolicExpr::cos(inv_x);
    auto expr = SymbolicExpr::multiply(SymbolicExpr::multiply(x, sin_inv_x), cos_inv_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x*sin(1/x)*cos(1/x), x->0) = 0";
}

TEST(SqueezeTheorem, ArctangentBound) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto atan_inv_x = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcTan,
                                               std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(inv_x)}));
    auto expr = SymbolicExpr::multiply(x, atan_inv_x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(x*arctan(1/x), x->0) = 0";
}

TEST(SqueezeTheorem, ConstantPlusSqueeze) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto five = SymbolicExpr::number(5);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto squeeze_term = SymbolicExpr::multiply(x, sin_inv_x);
    auto expr = SymbolicExpr::add(five, squeeze_term);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("5"))) << "limit(5 + x*sin(1/x), x->0) = 5";
}

TEST(SqueezeTheorem, ConstantPlusQuadraticSqueeze) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto neg_one = SymbolicExpr::number(-1);

    auto x_sq = SymbolicExpr::power(x, two);
    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto cos_inv_x = SymbolicExpr::cos(inv_x);
    auto squeeze_term = SymbolicExpr::multiply(x_sq, cos_inv_x);
    auto expr = SymbolicExpr::add(three, squeeze_term);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("3"))) << "limit(3 + x^2*cos(1/x), x->0) = 3";
}

TEST(SqueezeTheorem, SquaredBoundedFactor) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto two = SymbolicExpr::number(2);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto sin_sq = SymbolicExpr::power(sin_inv_x, two);
    auto expr = SymbolicExpr::multiply(sin_sq, x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit(sin(1/x)^2 * x, x->0) = 0";
}

TEST(SqueezeTheorem, BoundedSumFactor) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto two = SymbolicExpr::number(2);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto sin_inv_x = SymbolicExpr::sin(inv_x);
    auto bounded_expr = SymbolicExpr::add(two, sin_inv_x);
    auto expr = SymbolicExpr::multiply(bounded_expr, x);

    auto lim = LMCAS::limit_expression_checked(expr, "x", zero).value();
    EXPECT_TRUE(test_expression_text((lim), ("0"))) << "limit((2+sin(1/x))*x, x->0) = 0";
}
