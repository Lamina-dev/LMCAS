
#include "test_common.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include <stdexcept>
#include <string>

using namespace lamina;


void test_declare_parity_even() {
    TEST_CASE("Declare parity Even stores Even");
    PropertyStore store;
    store.declare_parity("x", Parity::Even);

    EXPECT_TRUE(store.get_parity("x") == Parity::Even, "x has Even parity");
}

void test_declare_parity_odd() {
    TEST_CASE("Declare parity Odd stores Odd");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);

    EXPECT_TRUE(store.get_parity("x") == Parity::Odd, "x has Odd parity");
}

void test_parity_default_unknown() {
    TEST_CASE("Default parity is Unknown for undeclared symbol");
    PropertyStore store;

    EXPECT_TRUE(store.get_parity("undeclared") == Parity::Unknown,
                "Undeclared symbol has Unknown parity");
}

void test_parity_auto_promotes_to_integer_from_complex() {
    TEST_CASE("Even parity auto-promotes domain from Complex to Integer");
    PropertyStore store;
    // Default domain is Complex
    EXPECT_TRUE(store.get_domain("x") == Domain::Complex, "x starts with Complex domain");

    store.declare_parity("x", Parity::Even);

    EXPECT_TRUE(store.get_domain("x") == Domain::Integer,
                "x domain promoted to Integer after Even parity declaration");
}

void test_parity_odd_auto_promotes_to_integer() {
    TEST_CASE("Odd parity auto-promotes domain from Complex to Integer");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);

    EXPECT_TRUE(store.get_domain("x") == Domain::Integer,
                "x domain promoted to Integer after Odd parity declaration");
}

void test_parity_auto_promotes_from_real() {
    TEST_CASE("Even parity auto-promotes domain from Real to Integer");
    PropertyStore store;
    store.declare_domain("x", Domain::Real);
    store.declare_parity("x", Parity::Even);

    EXPECT_TRUE(store.get_domain("x") == Domain::Integer,
                "x domain promoted from Real to Integer after Even parity");
}

void test_parity_does_not_demote_more_specific_domain() {
    TEST_CASE("Even parity does not demote Natural domain to Integer");
    PropertyStore store;
    store.declare_domain("x", Domain::Natural);
    store.declare_parity("x", Parity::Even);

    // Natural is more specific than Integer, so domain should remain Natural
    EXPECT_TRUE(store.get_domain("x") == Domain::Natural,
                "x domain remains Natural (more specific than Integer)");
}

void test_parity_does_not_demote_positiveint() {
    TEST_CASE("Odd parity does not demote PositiveInt domain to Integer");
    PropertyStore store;
    store.declare_domain("x", Domain::PositiveInt);
    store.declare_parity("x", Parity::Odd);

    EXPECT_TRUE(store.get_domain("x") == Domain::PositiveInt,
                "x domain remains PositiveInt (more specific than Integer)");
}

void test_parity_idempotent_even() {
    TEST_CASE("Idempotent re-declaration of Even parity");
    PropertyStore store;
    store.declare_parity("x", Parity::Even);
    store.declare_parity("x", Parity::Even);  // Should be no-op

    EXPECT_TRUE(store.get_parity("x") == Parity::Even,
                "x still has Even parity after re-declaration");
}

void test_parity_idempotent_odd() {
    TEST_CASE("Idempotent re-declaration of Odd parity");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);
    store.declare_parity("x", Parity::Odd);  // Should be no-op

    EXPECT_TRUE(store.get_parity("x") == Parity::Odd,
                "x still has Odd parity after re-declaration");
}

void test_parity_contradiction_even_then_odd() {
    TEST_CASE("Contradiction: Even then Odd returns failure");
    PropertyStore store;
    store.declare_parity("x", Parity::Even);

    auto failure_110 = store.declare_parity("x", Parity::Odd);
    EXPECT_TRUE(!failure_110.has_value(), "Declaring Odd after Even returns InvalidArgument");
    EXPECT_TRUE(store.get_parity("x") == Parity::Even,
                "x parity remains Even after failed Odd declaration");
}

void test_parity_contradiction_odd_then_even() {
    TEST_CASE("Contradiction: Odd then Even returns failure");
    PropertyStore store;
    store.declare_parity("x", Parity::Odd);

    auto failure_126 = store.declare_parity("x", Parity::Even);
    EXPECT_TRUE(!failure_126.has_value(), "Declaring Even after Odd returns InvalidArgument");
    EXPECT_TRUE(store.get_parity("x") == Parity::Odd,
                "x parity remains Odd after failed Even declaration");
}

void test_parity_unknown_can_be_set() {
    TEST_CASE("Setting parity to Unknown is allowed");
    PropertyStore store;
    store.declare_parity("x", Parity::Even);
    store.declare_parity("x", Parity::Unknown);

    EXPECT_TRUE(store.get_parity("x") == Parity::Unknown,
                "x parity set to Unknown");
}


void test_declare_bounded() {
    TEST_CASE("Declare Bounded stores Bounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
                "x has Bounded");
}

void test_declare_unbounded() {
    TEST_CASE("Declare Unbounded stores Unbounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
                "x has Unbounded");
}

void test_boundedness_default_unknown() {
    TEST_CASE("Default boundedness is Unknown for undeclared symbol");
    PropertyStore store;

    EXPECT_TRUE(store.get_boundedness("undeclared") == Boundedness::Unknown,
                "Undeclared symbol has Unknown boundedness");
}

void test_declare_bounded_with_interval() {
    TEST_CASE("Declare Bounded with Interval stores bounds");
    PropertyStore store;

    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    store.declare_bounded("x", Boundedness::Bounded, bounds);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
                "x has Bounded");
    EXPECT_TRUE(store.get_bounds("x").has_value(),
                "x has bounds stored");
}

void test_declare_bounded_without_interval() {
    TEST_CASE("Declare Bounded without Interval stores no bounds");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
                "x has Bounded");
    EXPECT_FALSE(store.get_bounds("x").has_value(),
                 "x has no bounds stored (none provided)");
}

void test_boundedness_idempotent_bounded() {
    TEST_CASE("Idempotent re-declaration of Bounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);
    store.declare_bounded("x", Boundedness::Bounded);  // Should be no-op

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
                "x still has Bounded after re-declaration");
}

void test_boundedness_idempotent_unbounded() {
    TEST_CASE("Idempotent re-declaration of Unbounded");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);
    store.declare_bounded("x", Boundedness::Unbounded);  // Should be no-op

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
                "x still has Unbounded after re-declaration");
}

void test_boundedness_contradiction_bounded_then_unbounded() {
    TEST_CASE("Contradiction: Bounded then Unbounded returns failure");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Bounded);

    auto failure_231 = store.declare_bounded("x", Boundedness::Unbounded);
    EXPECT_TRUE(!failure_231.has_value(), "Declaring Unbounded after Bounded returns InvalidArgument");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Bounded,
                "x boundedness remains Bounded after failed Unbounded declaration");
}

void test_boundedness_contradiction_unbounded_then_bounded() {
    TEST_CASE("Contradiction: Unbounded then Bounded returns failure");
    PropertyStore store;
    store.declare_bounded("x", Boundedness::Unbounded);

    auto failure_247 = store.declare_bounded("x", Boundedness::Bounded);
    EXPECT_TRUE(!failure_247.has_value(), "Declaring Bounded after Unbounded returns InvalidArgument");
    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unbounded,
                "x boundedness remains Unbounded after failed Bounded declaration");
}

void test_boundedness_unknown_can_be_set() {
    TEST_CASE("Setting boundedness to Unknown is allowed and clears bounds");
    PropertyStore store;

    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(10)));

    Interval bounds;
    bounds.lower = Endpoint::closed(lower_val);
    bounds.upper = Endpoint::closed(upper_val);

    store.declare_bounded("x", Boundedness::Bounded, bounds);
    store.declare_bounded("x", Boundedness::Unknown);

    EXPECT_TRUE(store.get_boundedness("x") == Boundedness::Unknown,
                "x boundedness set to Unknown");
    EXPECT_FALSE(store.get_bounds("x").has_value(),
                 "x bounds cleared when set to Unknown");
}

void test_undeclared_symbol_has_no_bounds() {
    TEST_CASE("Undeclared symbol has no bounds");
    PropertyStore store;

    EXPECT_FALSE(store.get_bounds("undeclared").has_value(),
                 "Undeclared symbol has no bounds");
}

void test_interval_queries_preserve_exact_large_endpoints() {
    TEST_CASE("PropertyStore interval queries preserve exact large endpoints");
    PropertyStore store;
    const BigInt two_to_53("9007199254740992");
    const BigInt next_integer = two_to_53 + BigInt(1);

    Interval first_point = Interval::point(SymbolicExpr::number(two_to_53));
    Interval second_point = Interval::point(SymbolicExpr::number(next_integer));
    store.declare_differentiable("f", first_point);
    auto success_278 = store.declare_continuous("f", second_point);
    EXPECT_TRUE(success_278.has_value(), "adjacent large integer points do not falsely overlap after exact comparison");

    Interval closed_span{
        Endpoint::closed(SymbolicExpr::number(two_to_53)),
        Endpoint::closed(SymbolicExpr::number(next_integer))
    };
    store.declare_continuous("g", closed_span);
    EXPECT_TRUE(store.is_continuous("g", second_point).value(),
                "closed span covers its exact large upper endpoint");

    Interval open_upper_span{
        Endpoint::closed(SymbolicExpr::number(two_to_53)),
        Endpoint::open(SymbolicExpr::number(next_integer))
    };
    store.declare_continuous("h", open_upper_span);
    EXPECT_TRUE(!store.is_continuous("h", second_point).value(),
                "open upper endpoint does not cover the exact large boundary point");
}


void test_parity_even_with_integer_domain_already_set() {
    TEST_CASE("Even parity with Integer domain already set is fine");
    PropertyStore store;
    store.declare_domain("x", Domain::Integer);
    store.declare_parity("x", Parity::Even);

    EXPECT_TRUE(store.get_domain("x") == Domain::Integer,
                "x domain remains Integer");
    EXPECT_TRUE(store.get_parity("x") == Parity::Even,
                "x has Even parity");
}

int main() {
    // Parity tests
    test_declare_parity_even();
    test_declare_parity_odd();
    test_parity_default_unknown();
    test_parity_auto_promotes_to_integer_from_complex();
    test_parity_odd_auto_promotes_to_integer();
    test_parity_auto_promotes_from_real();
    test_parity_does_not_demote_more_specific_domain();
    test_parity_does_not_demote_positiveint();
    test_parity_idempotent_even();
    test_parity_idempotent_odd();
    test_parity_contradiction_even_then_odd();
    test_parity_contradiction_odd_then_even();
    test_parity_unknown_can_be_set();

    // Boundedness tests
    test_declare_bounded();
    test_declare_unbounded();
    test_boundedness_default_unknown();
    test_declare_bounded_with_interval();
    test_declare_bounded_without_interval();
    test_boundedness_idempotent_bounded();
    test_boundedness_idempotent_unbounded();
    test_boundedness_contradiction_bounded_then_unbounded();
    test_boundedness_contradiction_unbounded_then_bounded();
    test_boundedness_unknown_can_be_set();
    test_undeclared_symbol_has_no_bounds();
    test_interval_queries_preserve_exact_large_endpoints();

    // Combined tests
    test_parity_even_with_integer_domain_already_set();

    return TEST_REPORT();
}
