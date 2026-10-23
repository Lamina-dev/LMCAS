#include <memory>
#include <vector>

#include "expr.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/equivalence_engine.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include "internal/visitors/expand_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

using namespace LMCAS;
using LMCAS::detail::SymbolicNodePtr;

class LmcasExpansionRegressions : public ::testing::Test {
  protected:
    SymbolicNodePtr a = detail::make_node<VariableNode>("a");
    SymbolicNodePtr b = detail::make_node<VariableNode>("b");
    SymbolicNodePtr c = detail::make_node<VariableNode>("c");
    SymbolicNodePtr d = detail::make_node<VariableNode>("d");
    SymbolicNodePtr left = detail::make_node<AddNode>(std::vector<SymbolicNodePtr>{a, b});
    SymbolicNodePtr right = detail::make_node<AddNode>(std::vector<SymbolicNodePtr>{c, d});
    SymbolicNodePtr input = detail::make_node<MultiplyNode>(std::vector<SymbolicNodePtr>{left, right});

    SymbolicNodePtr expanded_product() const {
        std::vector<SymbolicNodePtr> terms;
        for (const auto &first : {a, b}) {
            for (const auto &second : {c, d}) {
                terms.push_back(detail::make_node<MultiplyNode>(
                    std::vector<SymbolicNodePtr>{first, second}));
            }
        }
        return detail::make_node<AddNode>(std::move(terms));
    }
};

TEST_F(LmcasExpansionRegressions, ExactNodeBudgetBoundary) {
    auto expected = expanded_product();
    ComputationContext context;
    detail::RewriteBudget sufficient_budget(context, 64, 13, kEquivalentOperation);
    ExpandVisitor sufficient(sufficient_budget);
    input->accept(sufficient);
    EXPECT_TRUE(sufficient.get_result() && sufficient.get_result()->equals(*expected))
        << "thirteen nodes admit exactly a*c + a*d + b*c + b*d";

    detail::RewriteBudget limited_budget(context, 64, 12, kEquivalentOperation);
    ExpandVisitor limited(limited_budget);
    try {
        input->accept(limited);
        FAIL() << "twelve nodes must reject the expanded candidate";
    } catch (const CasError &error) {
        EXPECT_EQ(error.code, CasErrc::ResourceLimit);
    }
}

TEST_F(LmcasExpansionRegressions, ExpandedChildSumsShareOutputRoot) {
    auto e = detail::make_node<VariableNode>("e");
    auto nested_sum = detail::make_node<AddNode>(std::vector<SymbolicNodePtr>{input, e});
    auto expected_sum = detail::make_node<AddNode>(std::vector<SymbolicNodePtr>{expanded_product(), e});
    ComputationContext context;
    detail::RewriteBudget flattened_budget(context, 64, 14, kEquivalentOperation);
    ExpandVisitor flattened(flattened_budget);
    nested_sum->accept(flattened);
    EXPECT_TRUE(flattened.get_result() && flattened.get_result()->equals(*expected_sum))
        << "expanding (a+b)*(c+d)+e fits fourteen retained nodes";
}

TEST_F(LmcasExpansionRegressions, DistributesOneSumAcrossFactor) {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    NormalizationVisitor visitor(transform_budget, detail::no_facts(), Domain::Real);
    auto expanded = visitor.expand_product(left, c);
    auto sum = std::dynamic_pointer_cast<const AddNode>(expanded);
    ASSERT_NE(sum, nullptr);
    EXPECT_EQ(sum->operands().size(), 2u) << "(a+b)*c expands to two terms";
}

TEST_F(LmcasExpansionRegressions, DistributesTwoSums) {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    NormalizationVisitor visitor(transform_budget, detail::no_facts(), Domain::Real);
    auto expanded = visitor.expand_product(left, right);
    auto sum = std::dynamic_pointer_cast<const AddNode>(expanded);
    ASSERT_NE(sum, nullptr);
    EXPECT_EQ(sum->operands().size(), 4u) << "(a+b)*(c+d) expands to four terms";
}

TEST_F(LmcasExpansionRegressions, CubicBinomialCoefficients) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        SymbolicExpr::number(3));
    auto expected = parse_expr("x^3+3*x^2+3*x+1");
    ASSERT_TRUE(expected);

    auto expanded = LMCAS::expand(expression).value();
    ASSERT_TRUE(expanded);
    EXPECT_TRUE(test_proved_equivalent(expanded, expected.value()));
}

TEST_F(LmcasExpansionRegressions, SharedContextKeepsConsumedSteps) {
    auto expression = detail::make_expression_ptr(input);
    ComputationContext measured;
    auto first = LMCAS::expand(expression, measured);
    ASSERT_TRUE(first);
    ASSERT_GT(measured.steps_used(), 1u);

    ResourceLimits limits;
    limits.max_steps = measured.steps_used() + 1;
    ComputationContext context(limits);
    ASSERT_TRUE(LMCAS::expand(expression, context));
    const auto used = context.steps_used();
    auto second = LMCAS::expand(expression, context);
    ASSERT_FALSE(second);
    EXPECT_EQ(second.error().code, CasErrc::ResourceLimit);
    EXPECT_GT(context.steps_used(), used);
}

TEST_F(LmcasExpansionRegressions, DifferentiationConsumesInnerSteps) {
    ResourceLimits limits;
    limits.max_steps = 8;
    ComputationContext context(limits);
    auto expression = SymbolicExpr::power(
        SymbolicExpr::variable("x"), SymbolicExpr::number(4));
    auto derivative = LMCAS::differentiate(expression, "x", context);
    ASSERT_FALSE(derivative);
    EXPECT_EQ(derivative.error().code, CasErrc::ResourceLimit);
    EXPECT_GT(context.steps_used(), 1u);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST_F(LmcasExpansionRegressions, EmptySubstitutionsObserveCancellation) {
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext context({}, cancellation);
    auto result = LMCAS::substitute(detail::make_expression_ptr(input),
                                    std::vector<Binding>{}, context);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Cancelled);
}

TEST_F(LmcasExpansionRegressions, SubstitutionAvoidsIntegralCapture) {
    auto x = detail::make_node<VariableNode>("x");
    auto y = detail::make_node<VariableNode>("y");
    auto body = detail::make_node<AddNode>(std::vector<SymbolicNodePtr>{x, y});
    auto integral = detail::make_expression_ptr(
        detail::make_node<IntegralNode>(body, "y"));
    ComputationContext context;
    auto result = LMCAS::substitute(integral, "x",
                                    detail::make_expression_ptr(y), context);
    ASSERT_TRUE(result);
    auto node = std::dynamic_pointer_cast<const IntegralNode>(
        detail::node(result.value()));
    ASSERT_NE(node, nullptr);
    EXPECT_NE(node->variable(), "y");
    auto rebound = LMCAS::substitute(result.value(), "y", SymbolicExpr::number(7),
                                     context);
    ASSERT_TRUE(rebound);
    auto rebound_integral = std::dynamic_pointer_cast<const IntegralNode>(
        detail::node(rebound.value()));
    ASSERT_NE(rebound_integral, nullptr);
    auto expected_body = SymbolicExpr::add(
        SymbolicExpr::number(7), SymbolicExpr::variable(node->variable()));
    EXPECT_TRUE(test_proved_equivalent(
        detail::make_expression_ptr(rebound_integral->body()), expected_body));
}

TEST_F(LmcasExpansionRegressions, MatrixExpansionSharesUnexpandedElements) {
    auto matrix = detail::make_node<MatrixNode>(1, 1, MatrixNode::DenseStorage{input});
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 100, "test.expand");
    ExpandVisitor visitor(budget);
    matrix->accept(visitor);
    auto expanded = std::dynamic_pointer_cast<const MatrixNode>(visitor.get_result());
    ASSERT_NE(expanded, nullptr);
    EXPECT_EQ(expanded->get(0, 0).get(), input.get());
}

TEST_F(LmcasExpansionRegressions, IntegralDerivativeSharesImmutableBody) {
    auto integral = detail::make_node<IntegralNode>(input, "x");
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 100, "test.differentiate");
    DifferentiationVisitor visitor("x", budget);
    integral->accept(visitor);
    EXPECT_EQ(visitor.get_result().get(), input.get());
}
