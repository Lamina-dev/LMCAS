
#include "test_common.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "computation_context.hpp"
#include <string>

using namespace LMCAS;

static Interval exact_closed_interval(int lower, int upper) {
    return Interval{Endpoint::closed(SymbolicExpr::number(lower)),
                    Endpoint::closed(SymbolicExpr::number(upper))};
}

static void expect_exact_bounds(const PropertyStore &store,
                                const std::string &symbol,
                                int lower, int upper) {
    const auto bounds = store.get_bounds(symbol);
    ASSERT_TRUE(bounds);
    ASSERT_TRUE(bounds->lower.value);
    ASSERT_TRUE(bounds->upper.value);
    EXPECT_FALSE(bounds->lower.is_open);
    EXPECT_FALSE(bounds->upper.is_open);
    EXPECT_DOUBLE_EQ(bounds->lower.value->to_numeric(),
                     static_cast<double>(lower));
    EXPECT_DOUBLE_EQ(bounds->upper.value->to_numeric(),
                     static_cast<double>(upper));
}

TEST(PropertyStoreParityBounded, DeclareParityEven) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_parity("x") == Parity::Even)) << "x has Even parity";
}

TEST(PropertyStoreParityBounded, DeclareParityOdd) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_parity("x") == Parity::Odd)) << "x has Odd parity";
}

TEST(PropertyStoreParityBounded, ParityDefaultUnknown) {
    PropertyStore store;

    EXPECT_TRUE((store.get_parity("undeclared") == Parity::Unknown)) << "Undeclared symbol has Unknown parity";
}

TEST(PropertyStoreParityBounded, ParityAutoPromotesToIntegerFromComplex) {
    PropertyStore store;
    // Default domain is Complex
    EXPECT_TRUE((store.get_domain("x") == Domain::Complex)) << "x starts with Complex domain";

    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_domain("x") == Domain::Integer)) << "x domain promoted to Integer after Even parity declaration";
}

TEST(PropertyStoreParityBounded, ParityOddAutoPromotesToInteger) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_domain("x") == Domain::Integer)) << "x domain promoted to Integer after Odd parity declaration";
}

TEST(PropertyStoreParityBounded, ParityAutoPromotesFromReal) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Real).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_domain("x") == Domain::Integer)) << "x domain promoted from Real to Integer after Even parity";
}

TEST(PropertyStoreParityBounded, ParityDoesNotDemoteMoreSpecificDomain) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Natural).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    // Natural is more specific than Integer, so domain should remain Natural
    EXPECT_TRUE((store.get_domain("x") == Domain::Natural)) << "x domain remains Natural (more specific than Integer)";
}

TEST(PropertyStoreParityBounded, ParityDoesNotDemotePositiveint) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::PositiveInt).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_domain("x") == Domain::PositiveInt)) << "x domain remains PositiveInt (more specific than Integer)";
}

TEST(PropertyStoreParityBounded, ParityIdempotentEven) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_parity("x") == Parity::Even)) << "x still has Even parity after re-declaration";
}

TEST(PropertyStoreParityBounded, ParityIdempotentOdd) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_parity("x") == Parity::Odd)) << "x still has Odd parity after re-declaration";
}

TEST(PropertyStoreParityBounded, ParityContradictionEvenThenOdd) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    auto failure_110 = store.declare_parity("x", Parity::Odd);
    EXPECT_TRUE((!failure_110.has_value())) << "Declaring Odd after Even returns InvalidArgument";
    EXPECT_TRUE((store.get_parity("x") == Parity::Even)) << "x parity remains Even after failed Odd declaration";
}

TEST(PropertyStoreParityBounded, ParityContradictionOddThenEven) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    auto failure_126 = store.declare_parity("x", Parity::Even);
    EXPECT_TRUE((!failure_126.has_value())) << "Declaring Even after Odd returns InvalidArgument";
    EXPECT_TRUE((store.get_parity("x") == Parity::Odd)) << "x parity remains Odd after failed Even declaration";
}

TEST(PropertyStoreParityBounded, ParityUnknownCanBeSet) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Unknown).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_parity("x") == Parity::Unknown)) << "x parity set to Unknown";
}

TEST(PropertyStoreParityBounded, DeclareBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "x has Bounded";
}

TEST(PropertyStoreParityBounded, DeclareUnbounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Unbounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Unbounded)) << "x has Unbounded";
}

TEST(PropertyStoreParityBounded, BoundednessDefaultUnknown) {
    PropertyStore store;

    EXPECT_TRUE((store.get_boundedness("undeclared") == Boundedness::Unknown)) << "Undeclared symbol has Unknown boundedness";
}

TEST(PropertyStoreParityBounded, DeclareBoundedWithInterval) {
    PropertyStore store;

    auto lower_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded, bounds).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "x has Bounded";
    EXPECT_TRUE((store.get_bounds("x").has_value())) << "x has bounds stored";
}

TEST(PropertyStoreParityBounded, DeclareBoundedWithoutInterval) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "x has Bounded";
    EXPECT_FALSE((store.get_bounds("x").has_value())) << "x has no bounds stored (none provided)";
}

TEST(PropertyStoreParityBounded, RepeatedBoundsIntersectTransactionally) {
    PropertyStore store;
    const auto initial_revision = store.revision();
    ASSERT_TRUE(store.declare_bounded("x", Boundedness::Bounded));
    EXPECT_EQ(store.revision(), initial_revision + 1);
    EXPECT_FALSE(store.get_bounds("x"));

    auto outer = exact_closed_interval(0, 10);
    const auto outer_revision = store.revision();
    ASSERT_TRUE(store.declare_bounded_checked(
        "x", Boundedness::Bounded, outer));
    EXPECT_EQ(store.revision(), outer_revision + 1);
    expect_exact_bounds(store, "x", 0, 10);

    auto inner = exact_closed_interval(2, 8);
    const auto inner_revision = store.revision();
    ASSERT_TRUE(store.declare_bounded_checked(
        "x", Boundedness::Bounded, inner));
    EXPECT_EQ(store.revision(), inner_revision + 1);
    expect_exact_bounds(store, "x", 2, 8);

    auto disjoint = exact_closed_interval(20, 30);
    const auto rejected_revision = store.revision();
    auto rejected = store.declare_bounded_checked(
        "x", Boundedness::Bounded, disjoint);
    ASSERT_FALSE(rejected);
    EXPECT_EQ(rejected.error().code, CasErrc::InvalidArgument);
    EXPECT_EQ(store.revision(), rejected_revision);
    EXPECT_EQ(store.get_boundedness("x"), Boundedness::Bounded);
    expect_exact_bounds(store, "x", 2, 8);
}

TEST(PropertyStoreParityBounded, FiniteClassificationAcceptsLaterBounds) {
    PropertyStore store;
    ASSERT_TRUE(store.declare_finiteness("x", Finiteness::Finite));
    ASSERT_EQ(store.get_boundedness("x"), Boundedness::Bounded);
    ASSERT_FALSE(store.get_bounds("x"));

    auto bounds = exact_closed_interval(-3, 4);
    const auto revision = store.revision();
    ASSERT_TRUE(store.declare_bounded_checked(
        "x", Boundedness::Bounded, bounds));
    EXPECT_EQ(store.revision(), revision + 1);
    EXPECT_EQ(store.get_finiteness("x"), Finiteness::Finite);
    expect_exact_bounds(store, "x", -3, 4);
}

TEST(PropertyStoreParityBounded, CheckedIntersectionFailureRollsBackState) {
    PropertyStore store;
    auto original = exact_closed_interval(2, 8);
    ASSERT_TRUE(store.declare_bounded_checked(
        "x", Boundedness::Bounded, original));

    auto narrower = exact_closed_interval(3, 7);
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled_context({}, token);
    const auto cancelled_revision = store.revision();
    auto cancelled = store.declare_bounded_checked(
        "x", Boundedness::Bounded, narrower, cancelled_context);
    ASSERT_FALSE(cancelled);
    EXPECT_EQ(cancelled.error().code, CasErrc::Cancelled);
    EXPECT_EQ(store.revision(), cancelled_revision);
    expect_exact_bounds(store, "x", 2, 8);

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted_context(limits);
    const auto exhausted_revision = store.revision();
    auto exhausted = store.declare_bounded_checked(
        "x", Boundedness::Bounded, narrower, exhausted_context);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(store.revision(), exhausted_revision);
    expect_exact_bounds(store, "x", 2, 8);
}

TEST(PropertyStoreParityBounded, BoundednessIdempotentBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "x still has Bounded after re-declaration";
}

TEST(PropertyStoreParityBounded, BoundednessIdempotentUnbounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Unbounded).has_value())) << "boundedness declaration succeeds";
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Unbounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Unbounded)) << "x still has Unbounded after re-declaration";
}

TEST(PropertyStoreParityBounded, BoundednessContradictionBoundedThenUnbounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    auto failure_231 = store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE((!failure_231.has_value())) << "Declaring Unbounded after Bounded returns InvalidArgument";
    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Bounded)) << "x boundedness remains Bounded after failed Unbounded declaration";
}

TEST(PropertyStoreParityBounded, BoundednessContradictionUnboundedThenBounded) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Unbounded).has_value())) << "boundedness declaration succeeds";

    auto failure_247 = store.declare_bounded("x", Boundedness::Bounded);
    EXPECT_TRUE((!failure_247.has_value())) << "Declaring Bounded after Unbounded returns InvalidArgument";
    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Unbounded)) << "x boundedness remains Unbounded after failed Bounded declaration";
}

TEST(PropertyStoreParityBounded, BoundednessUnknownCanBeSet) {
    PropertyStore store;

    auto lower_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Bounded, bounds).has_value())) << "boundedness declaration succeeds";
    EXPECT_TRUE((store.declare_bounded("x", Boundedness::Unknown).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((store.get_boundedness("x") == Boundedness::Unknown)) << "x boundedness set to Unknown";
    EXPECT_FALSE((store.get_bounds("x").has_value())) << "x bounds cleared when set to Unknown";
}

TEST(PropertyStoreParityBounded, UndeclaredSymbolHasNoBounds) {
    PropertyStore store;

    EXPECT_FALSE((store.get_bounds("undeclared").has_value())) << "Undeclared symbol has no bounds";
}

TEST(PropertyStoreParityBounded, IntervalQueriesPreserveExactLargeEndpoints) {
    PropertyStore store;
    const BigInt two_to_53("9007199254740992");
    const BigInt next_integer = two_to_53 + BigInt(1);

    Interval first_point = Interval::point(SymbolicExpr::number(two_to_53));
    Interval second_point = Interval::point(SymbolicExpr::number(next_integer));
    EXPECT_TRUE((store.declare_differentiable("f", first_point).has_value())) << "differentiability declaration succeeds";
    auto success_278 = store.declare_continuous("f", second_point);
    EXPECT_TRUE((success_278.has_value())) << "adjacent large integer points do not falsely overlap after exact comparison";

    Interval closed_span{
        Endpoint::closed(SymbolicExpr::number(two_to_53)),
        Endpoint::closed(SymbolicExpr::number(next_integer))};
    EXPECT_TRUE((store.declare_continuous("g", closed_span).has_value())) << "continuity declaration succeeds";
    EXPECT_TRUE((store.is_continuous("g", second_point).value())) << "closed span covers its exact large upper endpoint";

    Interval open_upper_span{
        Endpoint::closed(SymbolicExpr::number(two_to_53)),
        Endpoint::open(SymbolicExpr::number(next_integer))};
    EXPECT_TRUE((store.declare_continuous("h", open_upper_span).has_value())) << "continuity declaration succeeds";
    EXPECT_TRUE((!store.is_continuous("h", second_point).value())) << "open upper endpoint does not cover the exact large boundary point";
}

TEST(PropertyStoreParityBounded, ParityEvenWithIntegerDomainAlreadySet) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Integer).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((store.get_domain("x") == Domain::Integer)) << "x domain remains Integer";
    EXPECT_TRUE((store.get_parity("x") == Parity::Even)) << "x has Even parity";
}
