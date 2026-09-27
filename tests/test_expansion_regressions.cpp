#include <memory>
#include <vector>

#include "expr.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/equivalence_engine.hpp"
#include "internal/rewrite_budget.hpp"
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
    ExpandVisitor sufficient(&sufficient_budget);
    input->accept(sufficient);
    EXPECT_TRUE(sufficient.get_result() && sufficient.get_result()->equals(*expected))
        << "thirteen nodes admit exactly a*c + a*d + b*c + b*d";

    detail::RewriteBudget limited_budget(context, 64, 12, kEquivalentOperation);
    ExpandVisitor limited(&limited_budget);
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
    ExpandVisitor flattened(&flattened_budget);
    nested_sum->accept(flattened);
    EXPECT_TRUE(flattened.get_result() && flattened.get_result()->equals(*expected_sum))
        << "expanding (a+b)*(c+d)+e fits fourteen retained nodes";
}

TEST_F(LmcasExpansionRegressions, DistributesOneSumAcrossFactor) {
    NormalizationVisitor visitor;
    auto expanded = visitor.expand_product(left, c);
    auto sum = std::dynamic_pointer_cast<const AddNode>(expanded);
    ASSERT_NE(sum, nullptr);
    EXPECT_EQ(sum->operands().size(), 2u) << "(a+b)*c expands to two terms";
}

TEST_F(LmcasExpansionRegressions, DistributesTwoSums) {
    NormalizationVisitor visitor;
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

    auto expanded = expression->expand();
    ASSERT_TRUE(expanded);
    EXPECT_TRUE(test_proved_equivalent(expanded, expected.value()));
}
