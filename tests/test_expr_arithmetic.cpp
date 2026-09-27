#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprArithmetic, ArithmeticEvaluation) {
    auto expr_two = LMCAS::integer(2);
    auto expr_three = LMCAS::integer(3);
    auto expr_sum = LMCAS::add(expr_two.value(), expr_three.value());
    auto expr_product =
        LMCAS::mul(expr_sum.value(), SymbolicExpr::number(4));
    auto expr_quotient =
        LMCAS::div(expr_product.value(), SymbolicExpr::number(2));
    auto expr_difference =
        LMCAS::sub(expr_quotient.value(), SymbolicExpr::number(5));
    auto expr_negated = LMCAS::neg(expr_difference.value());
    EXPECT_TRUE((expr_sum && expr_product && expr_quotient &&
                 expr_difference && expr_negated))
        << "Expr arithmetic wrappers construct symbolic expressions";

    auto expr_negated_value =
        expr_negated ? LMCAS::evalf(*expr_negated.value())
                     : LMCAS::Result<LMCAS::ApproxReal>::failure(
                           LMCAS::CasErrc::InternalInvariant,
                           "arithmetic construction failed", "test");
    EXPECT_NEAR(expr_negated_value.value().value, -5.0, 0.0) << "Expr arithmetic wrappers evaluate explicitly";
}

TEST(ExprArithmetic, ArithmeticSolutionSets) {
    auto expr_x = LMCAS::sym("expr_x");
    auto expr_five = LMCAS::integer(5);
    auto expr_polynomial =
        LMCAS::add(expr_x.value(), expr_five.value());
    auto expr_equation =
        LMCAS::eq(expr_polynomial.value(), SymbolicExpr::number(0));
    auto expr_solved =
        LMCAS::solve_expr_set(expr_polynomial.value(), "expr_x");
    EXPECT_TRUE((expr_equation && LMCAS::structurally_equal(
                                      *expr_equation.value(),
                                      *LMCAS::parse_expr("expr_x + 5 == 0").value())))
        << "Expr eq wrapper constructs a relational expression";
    EXPECT_TRUE((expr_solved &&
                 expr_solved.value().contains(*SymbolicExpr::number(-5))))
        << "Expr arithmetic wrappers feed the LMCAS solve set facade";
}

TEST(ExprArithmetic, ArithmeticInvalidInputs) {
    auto expr_two = LMCAS::integer(2);
    auto null_add = LMCAS::add(nullptr, expr_two.value());
    auto null_div = LMCAS::div(expr_two.value(), nullptr);
    auto null_neg = LMCAS::neg(nullptr);
    auto null_eq = LMCAS::eq(expr_two.value(), nullptr);
    EXPECT_TRUE((!null_add &&
                 null_add.error().code == LMCAS::CasErrc::InvalidArgument))
        << "Expr add rejects null input";
    EXPECT_TRUE((!null_div &&
                 null_div.error().code == LMCAS::CasErrc::InvalidArgument))
        << "Expr div rejects null input";
    EXPECT_TRUE((!null_neg &&
                 null_neg.error().code == LMCAS::CasErrc::InvalidArgument))
        << "Expr neg rejects null input";
    EXPECT_TRUE((!null_eq &&
                 null_eq.error().code == LMCAS::CasErrc::InvalidArgument))
        << "Expr eq rejects null input";
}

TEST(ExprArithmetic, ArithmeticBudget) {
    auto expr_two = LMCAS::integer(2);
    auto expr_three = LMCAS::integer(3);
    LMCAS::ResourceLimits exhausted_expr_limits;
    exhausted_expr_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_expr_context(exhausted_expr_limits);
    auto exhausted_expr_add =
        LMCAS::add(expr_two.value(), expr_three.value(),
                   exhausted_expr_context);
    EXPECT_TRUE((!exhausted_expr_add &&
                 exhausted_expr_add.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "Expr arithmetic wrappers observe the computation budget";
}

TEST(ExprArithmetic, SimplifyExpression) {
    auto transform_x = LMCAS::sym("transform_x");
    auto transform_zero = LMCAS::integer(0);
    auto transform_x_plus_zero =
        LMCAS::add(transform_x.value(), transform_zero.value());
    auto transform_simplified =
        LMCAS::simplify(transform_x_plus_zero.value());
    auto transform_simplified_value =
        transform_simplified
            ? LMCAS::evalf(*transform_simplified.value(),
                           LMCAS::NumericBindings{{"transform_x", 7.0}})
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "simplify construction failed", "test");
    EXPECT_TRUE((transform_simplified &&
                 transform_simplified_value &&
                 transform_simplified_value.value().value == 7.0))
        << "simplify lowers through the LMCAS Result facade";
}

TEST(ExprArithmetic, ExpandExpression) {
    auto transform_x = LMCAS::sym("transform_x");
    auto transform_one = LMCAS::integer(1);
    auto transform_two = LMCAS::integer(2);
    auto transform_left =
        LMCAS::add(transform_x.value(), transform_one.value());
    auto transform_right =
        LMCAS::add(transform_x.value(), transform_two.value());
    auto transform_product =
        LMCAS::mul(transform_left.value(), transform_right.value());
    auto transform_expanded =
        LMCAS::expand(transform_product.value());
    auto transform_expanded_value =
        transform_expanded
            ? LMCAS::evalf(*transform_expanded.value(),
                           LMCAS::NumericBindings{{"transform_x", 3.0}})
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "expand construction failed", "test");
    EXPECT_TRUE((transform_expanded && transform_expanded_value &&
                 transform_expanded_value.value().value == 20.0))
        << "expand lowers through the LMCAS Result facade";
}

TEST(ExprArithmetic, DifferentiateExpression) {
    auto transform_x = LMCAS::sym("transform_x");
    auto transform_three = LMCAS::integer(3);
    auto transform_x_cubed =
        LMCAS::pow(transform_x.value(), transform_three.value());
    auto transform_derivative =
        LMCAS::differentiate(transform_x_cubed.value(), "transform_x");
    auto transform_derivative_value =
        transform_derivative
            ? LMCAS::evalf(*transform_derivative.value(),
                           LMCAS::NumericBindings{{"transform_x", 2.0}})
            : LMCAS::Result<LMCAS::ApproxReal>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "differentiate construction failed", "test");
    EXPECT_TRUE((transform_derivative && transform_derivative_value &&
                 transform_derivative_value.value().value == 12.0))
        << "differentiate lowers through the LMCAS Result facade";
}

TEST(ExprArithmetic, TransformInvalidInputs) {
    auto transform_x = LMCAS::sym("transform_x");
    auto null_simplify = LMCAS::simplify(nullptr);
    auto null_expand = LMCAS::expand(nullptr);
    auto null_differentiate = LMCAS::differentiate(nullptr, "x");
    auto empty_variable =
        LMCAS::differentiate(transform_x.value(), "");
    EXPECT_TRUE((!null_simplify &&
                 null_simplify.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "simplify rejects null input";
    EXPECT_TRUE((!null_expand &&
                 null_expand.error().code == LMCAS::CasErrc::InvalidArgument))
        << "expand rejects null input";
    EXPECT_TRUE((!null_differentiate &&
                 null_differentiate.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "differentiate rejects null input";
    EXPECT_TRUE((!empty_variable &&
                 empty_variable.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "differentiate rejects empty variable names";
}

TEST(ExprArithmetic, TransformBudget) {
    auto transform_x = LMCAS::sym("transform_x");
    LMCAS::ResourceLimits exhausted_transform_limits;
    exhausted_transform_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_transform_context(
        exhausted_transform_limits);
    auto exhausted_transform =
        LMCAS::simplify(transform_x.value(), exhausted_transform_context);
    EXPECT_TRUE((!exhausted_transform &&
                 exhausted_transform.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "Expr transform wrappers observe the computation budget";
}

} // namespace
