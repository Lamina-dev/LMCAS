#include "limit_result.hpp"
#include "test_common.hpp"
#include "internal/symbolic_ast.hpp"

using namespace LMCAS;

namespace {

std::shared_ptr<SymbolicExpr> make_piecewise(
    const std::shared_ptr<SymbolicExpr> &left_expression,
    RelationalNode::Op left_op,
    const std::shared_ptr<SymbolicExpr> &right_expression,
    RelationalNode::Op right_op,
    const std::shared_ptr<SymbolicExpr> &point) {
    auto x = detail::node(SymbolicExpr::variable("x"));
    auto boundary = detail::node(point);
    return detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{
            {detail::node(left_expression),
             detail::make_node<RelationalNode>(x, boundary, left_op)},
            {detail::node(right_expression),
             detail::make_node<RelationalNode>(x, boundary, right_op)}}));
}

std::shared_ptr<SymbolicExpr> make_unary(
    FunctionNode::FuncType type,
    const std::shared_ptr<SymbolicExpr> &argument) {
    return detail::make_expression_ptr(detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{
                  detail::node(argument)}));
}

void expect_finite_limit(const LimitResult &result,
                         const std::shared_ptr<SymbolicExpr> &expected) {
    ASSERT_TRUE(result.has_value()) << result.error().message;
    const auto *finite = std::get_if<FiniteLimit>(&result.value().value);
    ASSERT_NE(finite, nullptr);
    EXPECT_TRUE(test_proved_equivalent(finite->value, expected));
}

} // namespace

TEST(PiecewiseLimit, DirectionalBranchesUseCheckedOutcomes) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto piecewise = make_piecewise(
        SymbolicExpr::add(x, SymbolicExpr::number(1)),
        RelationalNode::Op::GT,
        SymbolicExpr::add(x, SymbolicExpr::number(-1)),
        RelationalNode::Op::LT,
        zero);

    expect_finite_limit(
        limit_checked(piecewise, "x", zero, LimitDirection::FromAbove),
        SymbolicExpr::number(1));
    expect_finite_limit(
        limit_checked(piecewise, "x", zero, LimitDirection::FromBelow),
        SymbolicExpr::number(-1));

    auto two_sided = limit_checked(
        piecewise, "x", zero, LimitDirection::Both);
    ASSERT_TRUE(two_sided.has_value()) << two_sided.error().message;
    EXPECT_TRUE(std::holds_alternative<LimitDoesNotExist>(
        two_sided.value().value));
}

TEST(PiecewiseLimit, EqualBranchesHaveFiniteTwoSidedLimit) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto square = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto piecewise = make_piecewise(
        square, RelationalNode::Op::GT,
        square, RelationalNode::Op::LT, zero);

    expect_finite_limit(
        limit_checked(piecewise, "x", zero, LimitDirection::Both),
        zero);
}

TEST(PiecewiseLimit, LargeFiniteBoundaryKeepsRelationOrientation) {
    auto x = detail::node(SymbolicExpr::variable("x"));
    auto point = detail::node(SymbolicExpr::number(1.0e308));
    auto one = detail::node(SymbolicExpr::number(1));
    auto negative_one = detail::node(SymbolicExpr::number(-1));
    auto boundary = detail::make_expression_ptr(point);

    auto variable_left = detail::make_expression_ptr(
        detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{
                {one, detail::make_node<RelationalNode>(
                          x, point, RelationalNode::Op::GT)},
                {negative_one, detail::make_node<RelationalNode>(
                                   x, point, RelationalNode::Op::LEQ)}}));
    expect_finite_limit(
        limit_checked(variable_left, "x", boundary,
                      LimitDirection::FromAbove),
        SymbolicExpr::number(1));

    auto variable_right = detail::make_expression_ptr(
        detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{
                {one, detail::make_node<RelationalNode>(
                          point, x, RelationalNode::Op::LT)},
                {negative_one, detail::make_node<RelationalNode>(
                                   point, x, RelationalNode::Op::GEQ)}}));
    expect_finite_limit(
        limit_checked(variable_right, "x", boundary,
                      LimitDirection::FromAbove),
        SymbolicExpr::number(1));
}

TEST(PiecewiseLimit, SignAndAbsoluteDirectionalLimitsAreFinite) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto sign = make_unary(FunctionNode::FuncType::Sgn, x);
    auto absolute = make_unary(FunctionNode::FuncType::Abs, x);

    expect_finite_limit(
        limit_checked(sign, "x", zero, LimitDirection::FromAbove),
        SymbolicExpr::number(1));
    expect_finite_limit(
        limit_checked(sign, "x", zero, LimitDirection::FromBelow),
        SymbolicExpr::number(-1));
    expect_finite_limit(
        limit_checked(absolute, "x", zero, LimitDirection::FromAbove),
        zero);
    expect_finite_limit(
        limit_checked(absolute, "x", zero, LimitDirection::FromBelow),
        zero);
}

TEST(PiecewiseLimit, PoleDirectionsUseTypedInfinityOutcomes) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto reciprocal = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto inverse_square = SymbolicExpr::power(x, SymbolicExpr::number(-2));

    auto right = limit_checked(
        reciprocal, "x", zero, LimitDirection::FromAbove);
    ASSERT_TRUE(right.has_value()) << right.error().message;
    EXPECT_TRUE(std::holds_alternative<PositiveInfinityLimit>(
        right.value().value));

    auto left = limit_checked(
        reciprocal, "x", zero, LimitDirection::FromBelow);
    ASSERT_TRUE(left.has_value()) << left.error().message;
    EXPECT_TRUE(std::holds_alternative<NegativeInfinityLimit>(
        left.value().value));

    for (auto direction :
         {LimitDirection::FromAbove, LimitDirection::FromBelow}) {
        auto even = limit_checked(inverse_square, "x", zero, direction);
        ASSERT_TRUE(even.has_value()) << even.error().message;
        EXPECT_TRUE(std::holds_alternative<PositiveInfinityLimit>(
            even.value().value));
    }
}

TEST(PiecewiseLimit, NonstrictBoundaryHasMatchingDirectionalLimits) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto piecewise = make_piecewise(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        RelationalNode::Op::GEQ,
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(-1)),
        RelationalNode::Op::LT,
        one);

    for (auto direction : {LimitDirection::FromAbove,
                           LimitDirection::FromBelow,
                           LimitDirection::Both}) {
        expect_finite_limit(
            limit_checked(piecewise, "x", one, direction), one);
    }
}

TEST(PiecewiseLimit, ZeroLimitIsNotEventualEquality) {
    auto reciprocal = SymbolicExpr::power(
        SymbolicExpr::variable("x"), SymbolicExpr::number(-1));
    auto condition = detail::make_node<RelationalNode>(
        detail::node(reciprocal), detail::node(SymbolicExpr::number(0)),
        RelationOp::EQ);
    auto piecewise = detail::make_expression_ptr(
        detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{
                {detail::node(SymbolicExpr::number(1)), condition}},
            detail::node(SymbolicExpr::number(2))));

    expect_finite_limit(
        limit_checked(piecewise, "x", SymbolicExpr::infinity(),
                      LimitDirection::Both),
        SymbolicExpr::number(2));
}
