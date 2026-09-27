#include "test_common.hpp"
#include "internal/equivalence_engine.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/visitors/normalization_visitor.hpp"

using namespace LMCAS;

TEST(AstDepth, FixedAstDepthLimit) {
    auto expr = SymbolicExpr::variable("x");
    for (int i = 0; i < 510; ++i) {
        expr = SymbolicExpr::sin(expr);
    }

    bool caught = false;
    try {
        expr->simplify();
    } catch (const std::runtime_error &) {
        caught = true;
    }

    EXPECT_TRUE((caught)) << "AST depth limit should trigger an exception on deep trees";
}

TEST(AstDepth, RuntimeNormalizationDepth) {
    auto x = detail::make_node<VariableNode>("x");
    auto y = detail::make_node<VariableNode>("y");
    auto implication = detail::make_node<LogicalNode>(x, y, LogicalNode::Op::Implies);
    ComputationContext context;
    detail::RewriteBudget budget(context, 3, 100, kEquivalentOperation);
    NormalizationVisitor limited(&budget);
    bool exhausted = false;
    try {
        implication->accept(limited);
    } catch (const CasError &error) {
        exhausted = error.code == CasErrc::ResourceLimit;
    }
    EXPECT_TRUE((exhausted)) << "generated Or/Not recursion exceeds depth three";
    x->accept(limited);
    EXPECT_TRUE((limited.get_result() && limited.get_result()->equals(*x))) << "the same visitor, budget and context accept x after exception unwinding";

    detail::RewriteBudget sufficient_budget(context, 4, 100, kEquivalentOperation);
    NormalizationVisitor sufficient(&sufficient_budget);
    implication->accept(sufficient);
    auto not_x = detail::make_node<LogicalNode>(x, nullptr, LogicalNode::Op::Not);
    auto expected = detail::make_node<LogicalNode>(not_x, y, LogicalNode::Op::Or);
    EXPECT_TRUE((sufficient.get_result() && sufficient.get_result()->equals(*expected))) << "depth four admits exactly Or(Not(x), y)";
}
