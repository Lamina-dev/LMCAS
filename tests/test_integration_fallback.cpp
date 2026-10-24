#include "test_common.hpp"
#include "integration.hpp"
#include <string>
#include "expr.hpp"
#include "numeric_evaluation.hpp"
#include "internal/normalization_utils.hpp"

using namespace LMCAS;

TEST(LmcasIntegrationFallback, IndependentVariablePrimitive) {
    auto y = SymbolicExpr::variable("y");
    Integrator integrator;
    ComputationContext context;
    auto res = integrator.integrate_checked(*y, "x", context);
    ASSERT_TRUE(res) << res.error().message;
    std::string s = res.value().to_string();
    bool has_xy = (s.find("y*x") != std::string::npos) || (s.find("x*y") != std::string::npos);
    bool has_unevaluated = (s.find("integral(") != std::string::npos) || (s.find("Integral(") != std::string::npos);
    EXPECT_TRUE((has_xy)) << ("Integration of y w.r.t x should produce a product of x and y");
    EXPECT_TRUE((!has_unevaluated)) << ("Result should not be an unevaluated integral");
}

static void expect_partial_structure(const ExprPtr& actual, const ExprPtr& expected) {
    class SplitIntegral final : public detail::SymbolicRewriter {
    public:
        detail::SymbolicNodePtr remainder;
        std::size_t count = 0;
        void visit(const IntegralNode& integral) override {
            EXPECT_EQ(integral.variable(), "x");
            remainder = integral.body();
            ++count;
            set_result(detail::node(SymbolicExpr::number(0)));
        }
    } actual_split, expected_split;
    auto primitive = actual_split.rewrite(detail::node(actual));
    auto expected_primitive = expected_split.rewrite(detail::node(expected));
    ASSERT_EQ(actual_split.count, 1u);
    ASSERT_EQ(expected_split.count, 1u);
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(primitive),
        detail::make_expression_ptr(expected_primitive))) << actual->to_string();
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(actual_split.remainder),
        detail::make_expression_ptr(expected_split.remainder))) << actual->to_string();
}

static void expect_partial_identity(const ExprPtr &source, const ExprPtr &expected) {
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    expect_partial_structure(detail::make_expression_ptr(result.value()), expected);
    auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
    ASSERT_TRUE(derivative) << derivative.error().message;
    auto proof = check_equivalent(derivative.value(), source, context);
    ASSERT_TRUE(proof) << proof.error().message;
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
}

TEST(LmcasIntegrationFallback, PartialSum) {
    auto x = SymbolicExpr::variable("x");
    auto unsolved = SymbolicExpr::sin(SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto primitive = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    expect_partial_identity(SymbolicExpr::add(x, unsolved),
        SymbolicExpr::add(primitive, SymbolicExpr::make_integral(unsolved, "x")));
}

TEST(LmcasIntegrationFallback, ScaledPartialSum) {
    auto x = SymbolicExpr::variable("x");
    auto unsolved = SymbolicExpr::sin(SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto primitive = SymbolicExpr::multiply(SymbolicExpr::number(Rational(3, 2)),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    expect_partial_identity(SymbolicExpr::multiply(SymbolicExpr::number(3),
        SymbolicExpr::add(x, unsolved)),
        SymbolicExpr::add(primitive, SymbolicExpr::make_integral(
            SymbolicExpr::multiply(SymbolicExpr::number(3), unsolved), "x")));
}

TEST(LmcasIntegrationFallback, UnsolvedLeafRetainsOriginal) {
    auto x = SymbolicExpr::variable("x");
    auto source = SymbolicExpr::sin(SymbolicExpr::power(x, SymbolicExpr::number(3)));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    auto integral = std::dynamic_pointer_cast<const IntegralNode>(detail::node(result.value()));
    ASSERT_NE(integral, nullptr);
    EXPECT_EQ(integral->body(), detail::node(source));
    EXPECT_EQ(integral->variable(), "x");
}

TEST(LmcasIntegrationFallback, RepeatedQuadraticRetainsProperFraction) {
    auto x = SymbolicExpr::variable("x");
    auto denominator = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto fraction = SymbolicExpr::power(denominator, SymbolicExpr::number(-2));
    auto primitive = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    expect_partial_identity(SymbolicExpr::add(x, fraction),
        SymbolicExpr::add(primitive, SymbolicExpr::make_integral(fraction, "x")));
}

TEST(LmcasIntegrationFallback, DecimalPartialSumRetainsOriginalLeaf) {
    auto x = SymbolicExpr::variable("x");
    auto unsolved = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(0.3),
        SymbolicExpr::power(x, SymbolicExpr::number(3))));
    auto source = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(0.1), x), unsolved);
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    class IntegralBody final : public detail::RecursiveSymbolicVisitor {
      public:
        detail::SymbolicNodePtr body;
        void visit(const IntegralNode &node) override { body = node.body(); }
    } visitor;
    detail::node(result.value())->accept(visitor);
    ASSERT_NE(visitor.body, nullptr);
    EXPECT_EQ(visitor.body, detail::node(unsolved));
    auto exact_leaf = SymbolicExpr::sin(SymbolicExpr::multiply(
        SymbolicExpr::number(Rational(3, 10)),
        SymbolicExpr::power(x, SymbolicExpr::number(3))));
    EXPECT_TRUE(structurally_equal(detail::expression_from_node(visitor.body), *exact_leaf));
    auto expected = SymbolicExpr::add(SymbolicExpr::multiply(
        SymbolicExpr::number(Rational(1, 20)),
        SymbolicExpr::power(x, SymbolicExpr::number(2))),
        SymbolicExpr::make_integral(exact_leaf, "x"));
    expect_partial_structure(detail::make_expression_ptr(result.value()), expected);
    auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
    ASSERT_TRUE(derivative) << derivative.error().message;
    auto value = evaluate_numeric(*derivative.value(), {{"x", 2.0}}, context);
    ASSERT_TRUE(value) << value.error().message;
    EXPECT_NEAR(value.value().value, 0.2 + std::sin(2.4), 1e-12);
}

TEST(LmcasIntegrationFallback, ImproperFractionPreservesIntegratedQuotient) {
    auto x = SymbolicExpr::variable("x");
    auto denominator = SymbolicExpr::power(SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(1)),
        SymbolicExpr::number(2));
    auto source = SymbolicExpr::divide(SymbolicExpr::add(
        SymbolicExpr::multiply(x, denominator), SymbolicExpr::number(1)), denominator);
    auto primitive = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    auto remainder = SymbolicExpr::divide(SymbolicExpr::number(1), denominator);
    expect_partial_identity(source,
        SymbolicExpr::add(primitive, SymbolicExpr::make_integral(remainder, "x")));
}

TEST(LmcasIntegrationFallback, PiecewisePreservesPartialBranches) {
    auto x = SymbolicExpr::variable("x");
    auto unsolved = SymbolicExpr::sin(SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto condition = detail::make_node<RelationalNode>(
        detail::node(x), detail::node(SymbolicExpr::number(0)), RelationOp::LT);
    auto source = detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{detail::node(SymbolicExpr::add(x, unsolved)), condition}},
        detail::node(SymbolicExpr::multiply(SymbolicExpr::number(2), x))));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    auto piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(detail::node(result.value()));
    ASSERT_NE(piecewise, nullptr);
    ASSERT_EQ(piecewise->branches().size(), 1u);
    EXPECT_TRUE(piecewise->branches()[0].condition->equals(*condition));
    auto first = detail::make_expression_ptr(piecewise->branches()[0].expression);
    auto expected_first = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
            SymbolicExpr::power(x, SymbolicExpr::number(2))),
        SymbolicExpr::make_integral(unsolved, "x"));
    expect_partial_structure(first, expected_first);
    ASSERT_NE(piecewise->default_expr(), nullptr);
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(piecewise->default_expr()),
        SymbolicExpr::power(x, SymbolicExpr::number(2))));
    auto derivative = differentiate(first, "x", context);
    ASSERT_TRUE(derivative) << derivative.error().message;
    auto proof = check_equivalent(derivative.value(), SymbolicExpr::add(x, unsolved), context);
    ASSERT_TRUE(proof) << proof.error().message;
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
}

TEST(LmcasIntegrationFallback, IndefiniteIntegralDependsOnItsVariable) {
    auto x = SymbolicExpr::variable("x");
    auto inner = SymbolicExpr::make_integral(x, "x");
    Integrator integrator;
    for (const auto& source : {inner, SymbolicExpr::multiply(SymbolicExpr::number(2), inner)}) {
        ComputationContext context;
        auto result = integrator.integrate_checked(*source, "x", context);
        ASSERT_TRUE(result) << result.error().message;
        auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
        ASSERT_TRUE(derivative) << derivative.error().message;
        EXPECT_TRUE(test_proved_equivalent(derivative.value(), source));
    }
}

TEST(LmcasIntegrationFallback, CancellingRemaindersComplete) {
    auto term = SymbolicExpr::sin(SymbolicExpr::power(
        SymbolicExpr::variable("x"), SymbolicExpr::number(3)));
    auto source = SymbolicExpr::add(term,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(result.value().is_zero()) << result.value().to_string();
}

TEST(LmcasIntegrationFallback, ScaledPiecewiseCompletes) {
    auto x = SymbolicExpr::variable("x");
    auto condition = detail::make_node<RelationalNode>(
        detail::node(x), detail::node(SymbolicExpr::number(0)), RelationOp::LT);
    auto piecewise = detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{detail::node(x), condition}},
        detail::node(SymbolicExpr::multiply(SymbolicExpr::number(2), x))));
    auto source = SymbolicExpr::multiply(SymbolicExpr::number(2), piecewise);
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    const auto* completed = dynamic_cast<const PiecewiseNode*>(detail::node(result.value()).get());
    ASSERT_NE(completed, nullptr);
    ASSERT_EQ(completed->branches().size(), 1u);
    EXPECT_TRUE(completed->branches()[0].condition->equals(*condition));
    const std::pair<detail::SymbolicNodePtr, double> branches[] = {
        {completed->branches()[0].expression, -2.0}, {completed->default_expr(), 2.0}};
    for (const auto& [branch, value] : branches) {
        auto derivative = differentiate(detail::make_expression_ptr(branch), "x", context);
        ASSERT_TRUE(derivative) << derivative.error().message;
        auto actual = evaluate_numeric(*derivative.value(), {{"x", value}}, context);
        ASSERT_TRUE(actual) << actual.error().message;
        EXPECT_DOUBLE_EQ(actual.value().value, value < 0 ? -4.0 : 8.0);
    }
}
