#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>
#include <initializer_list>
#include <type_traits>
#include <utility>

using namespace LMCAS;

namespace {

using ElementReference = decltype(*std::declval<const ExprSet &>().elements()[0]);
using ExpressionReference = decltype(*std::declval<const ExprSet &>().expression());
static_assert(std::is_const_v<std::remove_reference_t<ElementReference>>);
static_assert(std::is_const_v<std::remove_reference_t<ExpressionReference>>);
static_assert(!std::is_assignable_v<ElementReference, SymbolicExpr>);
static_assert(!std::is_assignable_v<ExpressionReference, SymbolicExpr>);

void expect_integer_set(const ExprSet &set, std::initializer_list<int> values) {
    EXPECT_TRUE((set.size() == values.size() && set.empty() == (values.size() == 0))) << "set size and emptiness agree with its mathematical members";
    std::vector<ExprPtr> expected_elements;
    for (int value : values) {
        auto expected = SymbolicExpr::number(value);
        EXPECT_TRUE((set.contains(*expected))) << "expected member remains in the set";
        expected_elements.push_back(std::move(expected));
    }
    ComputationContext context;
    auto expected = make_finite_set(std::move(expected_elements), context);
    EXPECT_TRUE((expected && set.expression() &&
                 structurally_equal(*expected.value(), *set.expression())))
        << "expression view represents exactly the expected finite set";
    for (const auto &element : set.elements()) {
        bool expected_member = false;
        for (int value : values) {
            expected_member = expected_member ||
                              structurally_equal(*element, *SymbolicExpr::number(value));
        }
        EXPECT_TRUE((expected_member && set.contains(*element))) << "read-only projection agrees with set membership";
    }
}

TEST(ExprSets, ImmutableWrapperViews) {
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);
    auto left = expr_set({two, one, one});
    auto right = expr_set({three, two});
    EXPECT_TRUE((left && right)) << "valid finite sets construct successfully";
    if (!left || !right)
        return;
    auto copied = left.value();
    auto united = left.value().set_union(right.value());
    auto common = left.value().intersection(right.value());
    auto different = left.value().difference(right.value());
    auto symmetric = left.value().symmetric_difference(right.value());
    auto element_copy = std::make_shared<SymbolicExpr>(*copied.elements()[0]);
    auto expression_copy = std::make_shared<SymbolicExpr>(*copied.expression());
    *one = *SymbolicExpr::number(10);
    *two = *SymbolicExpr::number(20);
    *three = *SymbolicExpr::number(30);
    *element_copy = *SymbolicExpr::number(40);
    *expression_copy = *SymbolicExpr::number(50);
    expect_integer_set(left.value(), {1, 2});
    expect_integer_set(copied, {1, 2});
    expect_integer_set(right.value(), {2, 3});
    expect_integer_set(united, {1, 2, 3});
    expect_integer_set(common, {2});
    expect_integer_set(different, {1});
    expect_integer_set(symmetric, {1, 3});
    for (const ExprSet *set : {&united, &common, &different, &symmetric}) {
        auto mutable_copy = std::make_shared<SymbolicExpr>(*set->elements()[0]);
        *mutable_copy = *SymbolicExpr::number(99);
        EXPECT_TRUE((!set->contains(*mutable_copy))) << "operation results do not expose assignable member aliases";
    }
    expect_integer_set(united, {1, 2, 3});
    expect_integer_set(common, {2});
    expect_integer_set(different, {1});
    expect_integer_set(symmetric, {1, 3});
}

TEST(ExprSets, CanonicalDefaultEmpty) {
    ExprSet empty;
    expect_integer_set(empty, {});
    auto populated = expr_set({SymbolicExpr::number(1)});
    ASSERT_TRUE((populated.has_value())) << "singleton construction succeeds";
    if (!populated)
        return;
    EXPECT_TRUE((!empty.contains(*SymbolicExpr::number(1)) &&
                 empty.subset_of(populated.value())))
        << "default empty membership and subset agree";
    expect_integer_set(empty.set_union(populated.value()), {1});
    expect_integer_set(empty.intersection(populated.value()), {});
    expect_integer_set(empty.difference(populated.value()), {});
    expect_integer_set(empty.symmetric_difference(populated.value()), {1});
}

TEST(ExprSets, ExactApproximateDistinction) {
    auto exact = rational(Rational(BigInt(1), BigInt(2)));
    auto approximate = approx_real(0.5);
    EXPECT_TRUE((exact && approximate)) << "exact and approximate literals construct";
    if (!exact || !approximate) {
        return;
    }
    auto set = expr_set({exact.value(), approximate.value(), exact.value()});
    EXPECT_TRUE((set && set.value().size() == 2 &&
                 set.value().contains(*exact.value()) &&
                 set.value().contains(*approximate.value())))
        << "structural deduplication keeps distinct exact and approximate values";
    if (!set) {
        return;
    }
    auto exact_set = expr_set({exact.value()});
    auto remainder = set.value().difference(exact_set.value());
    EXPECT_TRUE((remainder.size() == 1 &&
                 remainder.contains(*approximate.value()) &&
                 !remainder.contains(*exact.value())))
        << "set difference removes only the exact representation";
}

TEST(ExprSets, FiniteCollectionMembership) {
    auto one_a = SymbolicExpr::number(1);
    auto one_b = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto base_set = LMCAS::expr_set({one_a, one_b, two});
    EXPECT_TRUE((base_set && base_set.value().size() == 2)) << "set<Expr> removes structurally equal duplicates";
    EXPECT_TRUE((base_set && base_set.value().contains(*SymbolicExpr::number(1)))) << "set<Expr> membership uses structural equality";
    auto facade_contains_one = base_set
                                   ? LMCAS::expr_set_contains(base_set.value(),
                                                              SymbolicExpr::number(1))
                                   : LMCAS::Result<bool>::failure(LMCAS::CasErrc::InternalInvariant,
                                                                  "set construction failed", "test");
    EXPECT_TRUE((facade_contains_one && facade_contains_one.value())) << "set<Expr> in operator facade reports membership";
    auto facade_not_contains_three = base_set
                                         ? LMCAS::expr_set_not_contains(base_set.value(),
                                                                        SymbolicExpr::number(3))
                                         : LMCAS::Result<bool>::failure(LMCAS::CasErrc::InternalInvariant,
                                                                        "set construction failed", "test");
    EXPECT_TRUE((facade_not_contains_three && facade_not_contains_three.value())) << "set<Expr> not in operator facade reports non-membership";
    auto empty_expr_set = LMCAS::expr_set({});
    EXPECT_TRUE((empty_expr_set && empty_expr_set.value().empty())) << "set<Expr> can represent the empty finite set";
}

TEST(ExprSets, FiniteUnion) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto rhs_set = LMCAS::expr_set({SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto union_set = base_set.value().set_union(rhs_set.value());
    auto facade_union =
        LMCAS::expr_set_union(base_set.value(), rhs_set.value());
    EXPECT_TRUE((facade_union && facade_union.value().size() == 3)) << "set<Expr> union facade returns deduplicated elements";
    EXPECT_TRUE((union_set.size() == 3)) << "set<Expr> union returns deduplicated elements";
}

TEST(ExprSets, FiniteIntersection) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto rhs_set = LMCAS::expr_set({SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto intersection_set = base_set.value().intersection(rhs_set.value());
    auto facade_intersection = LMCAS::expr_set_intersection(
        base_set.value(), rhs_set.value());
    EXPECT_TRUE((facade_intersection &&
                 facade_intersection.value().size() == 1 &&
                 facade_intersection.value().contains(
                     *SymbolicExpr::number(2))))
        << "set<Expr> intersection facade keeps common elements";
    EXPECT_TRUE((intersection_set.size() == 1 &&
                 intersection_set.contains(*SymbolicExpr::number(2))))
        << "set<Expr> intersection keeps common elements";
}

TEST(ExprSets, FiniteDifference) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto rhs_set = LMCAS::expr_set({SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto difference_set = base_set.value().difference(rhs_set.value());
    auto facade_difference = LMCAS::expr_set_difference(
        base_set.value(), rhs_set.value());
    EXPECT_TRUE((facade_difference &&
                 facade_difference.value().size() == 1 &&
                 facade_difference.value().contains(
                     *SymbolicExpr::number(1))))
        << "set<Expr> difference facade removes right-hand elements";
    EXPECT_TRUE((difference_set.size() == 1 &&
                 difference_set.contains(*SymbolicExpr::number(1))))
        << "set<Expr> difference removes right-hand elements";
}

TEST(ExprSets, FiniteSymmetricDifference) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto rhs_set = LMCAS::expr_set({SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto symmetric = base_set.value().symmetric_difference(rhs_set.value());
    auto facade_symmetric = LMCAS::expr_set_symmetric_difference(
        base_set.value(), rhs_set.value());
    EXPECT_TRUE((facade_symmetric &&
                 facade_symmetric.value().size() == 2 &&
                 facade_symmetric.value().contains(
                     *SymbolicExpr::number(1)) &&
                 facade_symmetric.value().contains(
                     *SymbolicExpr::number(3))))
        << "set<Expr> xor facade follows symmetric difference semantics";
    EXPECT_TRUE((symmetric.size() == 2 &&
                 symmetric.contains(*SymbolicExpr::number(1)) &&
                 symmetric.contains(*SymbolicExpr::number(3))))
        << "set<Expr> symmetric difference follows xor semantics";
}

TEST(ExprSets, FiniteSubset) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto rhs_set = LMCAS::expr_set({SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto union_set = base_set.value().set_union(rhs_set.value());
    auto intersection_set = base_set.value().intersection(rhs_set.value());
    auto facade_subset =
        LMCAS::expr_set_subset(intersection_set, union_set);
    EXPECT_TRUE((facade_subset && facade_subset.value())) << "set<Expr> subset facade checks membership of every element";
    EXPECT_TRUE((intersection_set.subset_of(union_set))) << "set<Expr> subset checks membership of every element";
}

TEST(ExprSets, EmptySetIdentities) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto empty_expr_set = LMCAS::expr_set({});
    EXPECT_TRUE((empty_expr_set.value().subset_of(base_set.value()))) << "empty set<Expr> is a subset of every set<Expr>";
    EXPECT_TRUE((base_set.value().set_union(empty_expr_set.value()).size() ==
                 base_set.value().size()))
        << "set<Expr> union with empty preserves the left set";
    EXPECT_TRUE((base_set.value().intersection(empty_expr_set.value()).empty())) << "set<Expr> intersection with empty is empty";
    EXPECT_TRUE((base_set.value().difference(empty_expr_set.value()).size() ==
                 base_set.value().size()))
        << "set<Expr> difference by empty preserves the left set";
}

TEST(ExprSets, InvalidSetElements) {
    auto base_set = LMCAS::expr_set({SymbolicExpr::number(1), SymbolicExpr::number(2)});
    auto null_set = LMCAS::expr_set({nullptr});
    EXPECT_TRUE((!null_set &&
                 null_set.error().code == LMCAS::CasErrc::InvalidArgument))
        << "set<Expr> rejects null elements";
    EXPECT_TRUE((!null_set &&
                 std::string(LMCAS::error_name(null_set.error())) ==
                     "SetElementTypeMismatch"))
        << "set<Expr> construction exposes the LMCAS element type diagnostic";
    if (base_set) {
        auto null_membership =
            LMCAS::expr_set_contains(base_set.value(), nullptr);
        EXPECT_TRUE((!null_membership &&
                     null_membership.error().code ==
                         LMCAS::CasErrc::InvalidArgument))
            << "set<Expr> membership rejects null elements";
        EXPECT_TRUE((!null_membership &&
                     std::string(LMCAS::error_name(
                         null_membership.error())) == "SetElementTypeMismatch"))
            << "set<Expr> membership exposes the element type diagnostic";
    }
}

} // namespace
