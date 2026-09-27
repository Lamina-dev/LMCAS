#include "test_common.hpp"
#include "interval.hpp"
#include "expr.hpp"
#include <cstdint>
#include <cstring>
#include <limits>
#include <cmath>

using namespace LMCAS;

TEST(IntervalSerialize, EmptySetToString) {
    auto empty = IntervalUnion::empty();
    std::string s = empty.to_string();
    EXPECT_EQ((s), ("\xe2\x88\x85")) << "Empty set should be ∅";
}

TEST(IntervalSerialize, EntireRealLineToString) {
    auto entire = IntervalUnion::entire_line();
    std::string s = entire.to_string();
    EXPECT_EQ((s), ("(-\xe2\x88\x9e, +\xe2\x88\x9e)")) << "Entire line should be (-∞, +∞)";
}

TEST(IntervalSerialize, SingleFiniteClosedInterval35) {
    auto iv = IntervalUnion::from_single(Interval{
        Endpoint::closed(SymbolicExpr::number(3)),
        Endpoint::closed(SymbolicExpr::number(5))});
    std::string s = iv.to_string();
    EXPECT_EQ((s), ("[3, 5]")) << "Closed interval [3, 5]";
}

TEST(IntervalSerialize, OpenInterval14) {
    auto iv = IntervalUnion::from_single(Interval{
        Endpoint::open(SymbolicExpr::number(1)),
        Endpoint::open(SymbolicExpr::number(4))});
    std::string s = iv.to_string();
    EXPECT_EQ((s), ("(1, 4)")) << "Open interval (1, 4)";
}

TEST(IntervalSerialize, HalfOpenInterval27) {
    auto iv = IntervalUnion::from_single(Interval{
        Endpoint::closed(SymbolicExpr::number(2)),
        Endpoint::open(SymbolicExpr::number(7))});
    std::string s = iv.to_string();
    EXPECT_EQ((s), ("[2, 7)")) << "Half-open interval [2, 7)";
}

TEST(IntervalSerialize, Interval2) {
    auto iv = IntervalUnion::from_single(Interval{
        Endpoint::neg_inf(),
        Endpoint::open(SymbolicExpr::number(-2))});
    std::string s = iv.to_string();
    EXPECT_EQ((s), ("(-\xe2\x88\x9e, -2)")) << "Interval (-∞, -2)";
}

TEST(IntervalSerialize, Interval3) {
    auto iv = IntervalUnion::from_single(Interval{
        Endpoint::closed(SymbolicExpr::number(3)),
        Endpoint::pos_inf()});
    std::string s = iv.to_string();
    EXPECT_EQ((s), ("[3, +\xe2\x88\x9e)")) << "Interval [3, +∞)";
}

TEST(IntervalSerialize, Union23) {
    std::vector<Interval> ivs = {
        Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-2))},
        Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::pos_inf()}};
    auto u = IntervalUnion(ivs);
    std::string s = u.to_string();
    EXPECT_EQ((s), ("(-\xe2\x88\x9e, -2) \xe2\x88\xaa [3, +\xe2\x88\x9e)")) << "Union (-∞, -2) ∪ [3, +∞)";
}

TEST(IntervalSerialize, Parse) {
    auto result = IntervalUnion::parse("\xe2\x88\x85");
    EXPECT_TRUE((result.has_value())) << "Parse ∅ should succeed";
    EXPECT_TRUE((result->is_empty())) << "Parsed ∅ should be empty";
}

TEST(IntervalSerialize, ParseContract) {
    auto result = IntervalUnion::parse("(-\xe2\x88\x9e, +\xe2\x88\x9e)");
    EXPECT_TRUE((result.has_value())) << "Parse (-∞, +∞) should succeed";
    EXPECT_TRUE((result->is_entire_line())) << "Parsed (-∞, +∞) should be entire line";
}

TEST(IntervalSerialize, Parse35) {
    auto result = IntervalUnion::parse("[3, 5]");
    EXPECT_TRUE((result.has_value())) << "Parse [3, 5] should succeed";
    EXPECT_TRUE((!result->is_empty())) << "Parsed [3, 5] should not be empty";
    EXPECT_TRUE((result->contains(4.0))) << "Parsed [3, 5] should contain 4";
    EXPECT_TRUE((result->contains(3.0))) << "Parsed [3, 5] should contain 3";
    EXPECT_TRUE((result->contains(5.0))) << "Parsed [3, 5] should contain 5";
    EXPECT_TRUE((!result->contains(2.9))) << "Parsed [3, 5] should not contain 2.9";
    EXPECT_TRUE((!result->contains(5.1))) << "Parsed [3, 5] should not contain 5.1";
}

TEST(IntervalSerialize, Parse23) {
    auto result = IntervalUnion::parse("(-\xe2\x88\x9e, -2) \xe2\x88\xaa [3, +\xe2\x88\x9e)");
    EXPECT_TRUE((result.has_value())) << "Parse union should succeed";
    EXPECT_TRUE((result->contains(-5.0))) << "Should contain -5";
    EXPECT_TRUE((!result->contains(-2.0))) << "Should not contain -2 (open)";
    EXPECT_TRUE((!result->contains(0.0))) << "Should not contain 0";
    EXPECT_TRUE((result->contains(3.0))) << "Should contain 3 (closed)";
    EXPECT_TRUE((result->contains(100.0))) << "Should contain 100";
}

TEST(IntervalSerialize, ParseInvalidInput) {
    auto r1 = IntervalUnion::parse("");
    EXPECT_TRUE((!r1.has_value())) << "Empty string should return nullopt";

    auto r2 = IntervalUnion::parse("invalid");
    EXPECT_TRUE((!r2.has_value())) << "Invalid string should return nullopt";

    auto r3 = IntervalUnion::parse("[3, ]");
    EXPECT_TRUE((!r3.has_value())) << "Missing value should return nullopt";
}

TEST(IntervalSerialize, RoundTripToStringThenParse) {
    std::vector<Interval> ivs = {
        Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-2))},
        Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::pos_inf()}};
    auto original = IntervalUnion(ivs);
    std::string s = original.to_string();
    auto parsed = IntervalUnion::parse(s);
    EXPECT_TRUE((parsed.has_value())) << "Round-trip parse should succeed";

    EXPECT_TRUE((parsed->contains(-10.0) == original.contains(-10.0))) << "Round-trip: -10";
    EXPECT_TRUE((parsed->contains(-2.0) == original.contains(-2.0))) << "Round-trip: -2";
    EXPECT_TRUE((parsed->contains(0.0) == original.contains(0.0))) << "Round-trip: 0";
    EXPECT_TRUE((parsed->contains(3.0) == original.contains(3.0))) << "Round-trip: 3";
    EXPECT_TRUE((parsed->contains(10.0) == original.contains(10.0))) << "Round-trip: 10";
}

TEST(IntervalSerialize, ParseIntervalWithNegativeNumbers51) {
    auto result = IntervalUnion::parse("[-5, -1]");
    EXPECT_TRUE((result.has_value())) << "Parse [-5, -1] should succeed";
    EXPECT_TRUE((result->contains(-3.0))) << "Should contain -3";
    EXPECT_TRUE((result->contains(-5.0))) << "Should contain -5";
    EXPECT_TRUE((result->contains(-1.0))) << "Should contain -1";
    EXPECT_TRUE((!result->contains(0.0))) << "Should not contain 0";
}

TEST(IntervalSerialize, ParseIntervalWithDecimals1537) {
    auto result = IntervalUnion::parse("(1.5, 3.7)");
    EXPECT_TRUE((result.has_value())) << "Parse (1.5, 3.7) should succeed";
    EXPECT_TRUE((result->contains(2.0))) << "Should contain 2.0";
    EXPECT_TRUE((!result->contains(1.5))) << "Should not contain 1.5 (open)";
    EXPECT_TRUE((!result->contains(3.7))) << "Should not contain 3.7 (open)";
}

TEST(IntervalSerialize, ParsePreservesExactRationalAndArbitraryPrecisionEndpoints) {
    auto rational = IntervalUnion::parse("[1/3, 2/3]");
    ASSERT_TRUE((rational.has_value())) << "Exact rational endpoints should parse";
    if (rational) {
        EXPECT_EQ((rational->to_string()), ("[1/3, 2/3]")) << "Exact rational endpoints should round-trip without approximation";
    }

    const std::string large = "1234567890123456789012345678901234567890";
    auto arbitrary_precision = IntervalUnion::parse("[" + large + ", " + large + "]");
    ASSERT_TRUE((arbitrary_precision.has_value())) << "Arbitrary-precision integer endpoints should parse";
    if (arbitrary_precision) {
        EXPECT_EQ((arbitrary_precision->to_string()), ("[" + large + ", " + large + "]")) << "Arbitrary-precision integer endpoints should round-trip exactly";
    }
}

TEST(IntervalSerialize, ParseRejectsInvalidNumericIntervalsWithoutLeakingExceptions) {
    std::optional<IntervalUnion> invalid;
    EXPECT_NO_THROW(invalid = IntervalUnion::parse("[1/0, 2]"))
        << "Malformed external input remains inside the optional parser contract";
    EXPECT_TRUE((!invalid.has_value())) << "Zero rational denominator should be rejected";
    EXPECT_TRUE((!IntervalUnion::parse("[+\xe2\x88\x9e, 1]").has_value())) << "Positive infinity cannot be a lower endpoint";
}

static void expect_approximate_endpoint_values() {
    for (double expected : {0.5, 1.0, 1.2345678901234567,
                            std::numeric_limits<double>::denorm_min(),
                            std::nextafter(std::numeric_limits<double>::min(), 0.0),
                            std::numeric_limits<double>::min(), std::numeric_limits<double>::max(),
                            std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0), 0.0, -0.0}) {
        auto expression = approx_real(expected).value();
        auto original = IntervalUnion::from_single(Interval::point(expression));
        auto parsed = IntervalUnion::parse(original.to_string());
        EXPECT_TRUE((parsed && parsed->intervals().size() == 1)) << "approximate point interval roundtrip preserves its sole interval";
        if (!parsed || parsed->intervals().size() != 1) {
            continue;
        }
        for (const auto *endpoint : {&parsed->intervals()[0].lower,
                                     &parsed->intervals()[0].upper}) {
            auto value = evalf(*endpoint->value);
            ASSERT_TRUE((value.has_value())) << "parsed approximate endpoint evaluates";
            if (!value) {
                continue;
            }
            std::uint64_t actual_bits, expected_bits;
            std::memcpy(&actual_bits, &value.value().value, sizeof actual_bits);
            std::memcpy(&expected_bits, &expected, sizeof expected_bits);
            EXPECT_TRUE((actual_bits == expected_bits)) << "endpoint preserves all binary64 bits";
            EXPECT_TRUE((structurally_equal(*endpoint->value, *expression))) << "endpoint remains approximate rather than exact";
        }
    }
}

static void expect_approximate_interval_syntax() {
    auto open = IntervalUnion::parse("(approx( -0 ), approx( +1 ))");
    EXPECT_TRUE((open && !open->contains(0.0) && open->contains(0.5) && !open->contains(1.0))) << "atom and open interval closing delimiters remain distinct";
    auto legacy = IntervalUnion::parse("[0.5, 1.0]");
    EXPECT_TRUE((legacy && structurally_equal(*legacy->intervals()[0].lower.value, *approx_real(0.5).value()) &&
                 structurally_equal(*legacy->intervals()[0].upper.value,
                                    *approx_real(1.0).value())))
        << "legacy bare decimals retain approximate semantics";
}

TEST(IntervalSerialize, ApproximateEndpointRoundtrip) {
    expect_approximate_endpoint_values();
    expect_approximate_interval_syntax();
    for (const char *literal : {"approx()", "approx(x)", "approx(1/2)", "approx(1,2)",
                                "approx(1e)", "approx(nan)", "approx(inf)", "approx(0x1p0)",
                                "approx(1e400)", "approx(1e-4000)", "approx(1 2)", "approx(1)junk"}) {
        EXPECT_TRUE((!IntervalUnion::parse(std::string("[") + literal + ", 2]"))) << "interval parser rejects malformed approximate atoms";
    }
}
