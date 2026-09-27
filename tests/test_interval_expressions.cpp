#include "test_interval_support.hpp"

namespace {

struct ExpectedRelation {
    RelationalNode::Op op;
    int endpoint;
};

void collect_relations(
    const std::shared_ptr<const SymbolicNode> &node,
    std::vector<std::shared_ptr<const RelationalNode>> &relations) {
    if (auto relation =
            std::dynamic_pointer_cast<const RelationalNode>(node)) {
        relations.push_back(relation);
        return;
    }
    auto logical = std::dynamic_pointer_cast<const LogicalNode>(node);
    ASSERT_NE(logical, nullptr);
    collect_relations(logical->left(), relations);
    collect_relations(logical->right(), relations);
}

void expect_relations(
    const std::vector<std::shared_ptr<const RelationalNode>> &relations,
    const std::vector<ExpectedRelation> &expected) {
    ASSERT_EQ(relations.size(), expected.size());
    std::vector<bool> matched(expected.size(), false);
    for (const auto &relation : relations) {
        auto variable =
            std::dynamic_pointer_cast<const VariableNode>(relation->left());
        ASSERT_NE(variable, nullptr);
        EXPECT_EQ(variable->name(), "x");
        auto endpoint = detail::make_expression_ptr(relation->right());
        ASSERT_NE(endpoint, nullptr);

        bool found = false;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            if (!matched[i] && relation->op() == expected[i].op &&
                endpoint->compare(SymbolicExpr::number(
                    expected[i].endpoint)) == 0) {
                matched[i] = true;
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "unexpected interval relation";
    }
    for (bool relation_matched : matched) {
        EXPECT_TRUE(relation_matched) << "missing interval relation";
    }
}

void expect_relation_tree(
    const std::shared_ptr<SymbolicExpr> &expression,
    LogicalNode::Op logical_op,
    const std::vector<ExpectedRelation> &expected) {
    ASSERT_NE(expression, nullptr);
    auto logical = std::dynamic_pointer_cast<const LogicalNode>(
        detail::node(expression));
    ASSERT_NE(logical, nullptr);
    EXPECT_EQ(logical->op(), logical_op);
    std::vector<std::shared_ptr<const RelationalNode>> relations;
    collect_relations(detail::node(expression), relations);
    expect_relations(relations, expected);
}

void expect_single_relation(
    const std::shared_ptr<SymbolicExpr> &expression,
    RelationalNode::Op op,
    int endpoint) {
    ASSERT_NE(expression, nullptr);
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(
        detail::node(expression));
    ASSERT_NE(relation, nullptr);
    auto variable =
        std::dynamic_pointer_cast<const VariableNode>(relation->left());
    ASSERT_NE(variable, nullptr);
    EXPECT_EQ(variable->name(), "x");
    EXPECT_EQ(relation->op(), op);
    EXPECT_EQ(detail::make_expression_ptr(relation->right())->compare(
                  SymbolicExpr::number(endpoint)),
              0);
}

} // namespace

TEST(IntervalExpressions, ToExprEmptySetReturnsCanonicalFiniteSet) {
    auto expression = IntervalUnion::empty().to_expr("x");

    ASSERT_NE(expression, nullptr);
    auto empty = std::dynamic_pointer_cast<const FiniteSetNode>(
        detail::node(expression));
    ASSERT_NE(empty, nullptr);
    EXPECT_TRUE(empty->elements().empty());
}

TEST(IntervalExpressions, ToExprEntireLineReturnsCanonicalTrueRelation) {
    auto expression = IntervalUnion::entire_line().to_expr("x");

    ASSERT_NE(expression, nullptr);
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(
        detail::node(expression));
    ASSERT_NE(relation, nullptr);
    EXPECT_EQ(relation->op(), RelationalNode::Op::EQ);
    auto left = std::dynamic_pointer_cast<const NumberNode>(relation->left());
    auto right =
        std::dynamic_pointer_cast<const NumberNode>(relation->right());
    ASSERT_NE(left, nullptr);
    ASSERT_NE(right, nullptr);
    EXPECT_TRUE(left->is_zero());
    EXPECT_TRUE(right->is_zero());
}

TEST(IntervalExpressions, OpenFiniteIntervalUsesStrictRelations) {
    auto interval = IntervalUnion::from_single(Interval{
        Endpoint::open(SymbolicExpr::number(2)),
        Endpoint::open(SymbolicExpr::number(5))});
    expect_relation_tree(
        interval.to_expr("x"), LogicalNode::Op::And,
        {{RelationalNode::Op::GT, 2}, {RelationalNode::Op::LT, 5}});
}

TEST(IntervalExpressions, UnboundedIntervalsUseOneRelation) {
    expect_single_relation(
        IntervalUnion::from_single(Interval{
            Endpoint::open(SymbolicExpr::number(3)),
            Endpoint::pos_inf()})
            .to_expr("x"),
        RelationalNode::Op::GT, 3);
    expect_single_relation(
        IntervalUnion::from_single(Interval{
            Endpoint::neg_inf(),
            Endpoint::closed(SymbolicExpr::number(4))})
            .to_expr("x"),
        RelationalNode::Op::LEQ, 4);
}

TEST(IntervalExpressions, ClosedFiniteIntervalUsesNonstrictRelations) {
    auto interval = IntervalUnion::from_single(Interval{
        Endpoint::closed(SymbolicExpr::number(1)),
        Endpoint::closed(SymbolicExpr::number(3))});
    expect_relation_tree(
        interval.to_expr("x"), LogicalNode::Op::And,
        {{RelationalNode::Op::GEQ, 1},
         {RelationalNode::Op::LEQ, 3}});
}

TEST(IntervalExpressions, DisjointIntervalsUseRelationDisjunction) {
    auto two_intervals = IntervalUnion({
        Interval{Endpoint::neg_inf(),
                 Endpoint::open(SymbolicExpr::number(-2))},
        Interval{Endpoint::open(SymbolicExpr::number(3)),
                 Endpoint::pos_inf()}});
    expect_relation_tree(
        two_intervals.to_expr("x"), LogicalNode::Op::Or,
        {{RelationalNode::Op::LT, -2},
         {RelationalNode::Op::GT, 3}});

    auto three_intervals = IntervalUnion({
        Interval{Endpoint::neg_inf(),
                 Endpoint::open(SymbolicExpr::number(-5))},
        Interval{Endpoint::open(SymbolicExpr::number(-1)),
                 Endpoint::open(SymbolicExpr::number(1))},
        Interval{Endpoint::open(SymbolicExpr::number(5)),
                 Endpoint::pos_inf()}});
    auto expression = three_intervals.to_expr("x");
    ASSERT_NE(expression, nullptr);
    auto outer = std::dynamic_pointer_cast<const LogicalNode>(
        detail::node(expression));
    ASSERT_NE(outer, nullptr);
    EXPECT_EQ(outer->op(), LogicalNode::Op::Or);
    std::vector<std::shared_ptr<const RelationalNode>> relations;
    collect_relations(detail::node(expression), relations);
    expect_relations(
        relations,
        {{RelationalNode::Op::LT, -5},
         {RelationalNode::Op::GT, -1},
         {RelationalNode::Op::LT, 1},
         {RelationalNode::Op::GT, 5}});
}
