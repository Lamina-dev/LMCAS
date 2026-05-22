#include "test_common.hpp"
#include "interval.hpp"

using namespace lamina;

int main() {
    // Test 1: Empty set to_string
    TEST_CASE("Empty set to_string");
    {
        auto empty = IntervalUnion::empty();
        std::string s = empty.to_string();
        EXPECT_EQ_STR(s, "\xe2\x88\x85", "Empty set should be ∅");
    }

    // Test 2: Entire real line to_string
    TEST_CASE("Entire real line to_string");
    {
        auto entire = IntervalUnion::entire_line();
        std::string s = entire.to_string();
        EXPECT_EQ_STR(s, "(-\xe2\x88\x9e, +\xe2\x88\x9e)", "Entire line should be (-∞, +∞)");
    }

    // Test 3: Single finite interval [3, 5]
    TEST_CASE("Single finite closed interval [3, 5]");
    {
        auto iv = IntervalUnion::from_single(Interval{
            Endpoint::closed(SymbolicExpr::number(3)),
            Endpoint::closed(SymbolicExpr::number(5))
        });
        std::string s = iv.to_string();
        EXPECT_EQ_STR(s, "[3, 5]", "Closed interval [3, 5]");
    }

    // Test 4: Open interval (1, 4)
    TEST_CASE("Open interval (1, 4)");
    {
        auto iv = IntervalUnion::from_single(Interval{
            Endpoint::open(SymbolicExpr::number(1)),
            Endpoint::open(SymbolicExpr::number(4))
        });
        std::string s = iv.to_string();
        EXPECT_EQ_STR(s, "(1, 4)", "Open interval (1, 4)");
    }

    // Test 5: Half-open interval [2, 7)
    TEST_CASE("Half-open interval [2, 7)");
    {
        auto iv = IntervalUnion::from_single(Interval{
            Endpoint::closed(SymbolicExpr::number(2)),
            Endpoint::open(SymbolicExpr::number(7))
        });
        std::string s = iv.to_string();
        EXPECT_EQ_STR(s, "[2, 7)", "Half-open interval [2, 7)");
    }

    // Test 6: Interval with infinity (-∞, -2)
    TEST_CASE("Interval (-∞, -2)");
    {
        auto iv = IntervalUnion::from_single(Interval{
            Endpoint::neg_inf(),
            Endpoint::open(SymbolicExpr::number(-2))
        });
        std::string s = iv.to_string();
        EXPECT_EQ_STR(s, "(-\xe2\x88\x9e, -2)", "Interval (-∞, -2)");
    }

    // Test 7: Interval with infinity [3, +∞)
    TEST_CASE("Interval [3, +∞)");
    {
        auto iv = IntervalUnion::from_single(Interval{
            Endpoint::closed(SymbolicExpr::number(3)),
            Endpoint::pos_inf()
        });
        std::string s = iv.to_string();
        EXPECT_EQ_STR(s, "[3, +\xe2\x88\x9e)", "Interval [3, +∞)");
    }

    // Test 8: Union of two intervals (-∞, -2) ∪ [3, +∞)
    TEST_CASE("Union (-∞, -2) ∪ [3, +∞)");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-2))},
            Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::pos_inf()}
        };
        auto u = IntervalUnion(ivs);
        std::string s = u.to_string();
        EXPECT_EQ_STR(s, "(-\xe2\x88\x9e, -2) \xe2\x88\xaa [3, +\xe2\x88\x9e)",
                      "Union (-∞, -2) ∪ [3, +∞)");
    }

    // Test 9: Parse empty set
    TEST_CASE("Parse ∅");
    {
        auto result = IntervalUnion::parse("\xe2\x88\x85");
        EXPECT_TRUE(result.has_value(), "Parse ∅ should succeed");
        EXPECT_TRUE(result->is_empty(), "Parsed ∅ should be empty");
    }

    // Test 10: Parse entire line
    TEST_CASE("Parse (-∞, +∞)");
    {
        auto result = IntervalUnion::parse("(-\xe2\x88\x9e, +\xe2\x88\x9e)");
        EXPECT_TRUE(result.has_value(), "Parse (-∞, +∞) should succeed");
        EXPECT_TRUE(result->is_entire_line(), "Parsed (-∞, +∞) should be entire line");
    }

    // Test 11: Parse single interval [3, 5]
    TEST_CASE("Parse [3, 5]");
    {
        auto result = IntervalUnion::parse("[3, 5]");
        EXPECT_TRUE(result.has_value(), "Parse [3, 5] should succeed");
        EXPECT_TRUE(!result->is_empty(), "Parsed [3, 5] should not be empty");
        EXPECT_TRUE(result->contains(4.0), "Parsed [3, 5] should contain 4");
        EXPECT_TRUE(result->contains(3.0), "Parsed [3, 5] should contain 3");
        EXPECT_TRUE(result->contains(5.0), "Parsed [3, 5] should contain 5");
        EXPECT_TRUE(!result->contains(2.9), "Parsed [3, 5] should not contain 2.9");
        EXPECT_TRUE(!result->contains(5.1), "Parsed [3, 5] should not contain 5.1");
    }

    // Test 12: Parse union (-∞, -2) ∪ [3, +∞)
    TEST_CASE("Parse (-∞, -2) ∪ [3, +∞)");
    {
        auto result = IntervalUnion::parse("(-\xe2\x88\x9e, -2) \xe2\x88\xaa [3, +\xe2\x88\x9e)");
        EXPECT_TRUE(result.has_value(), "Parse union should succeed");
        EXPECT_TRUE(result->contains(-5.0), "Should contain -5");
        EXPECT_TRUE(!result->contains(-2.0), "Should not contain -2 (open)");
        EXPECT_TRUE(!result->contains(0.0), "Should not contain 0");
        EXPECT_TRUE(result->contains(3.0), "Should contain 3 (closed)");
        EXPECT_TRUE(result->contains(100.0), "Should contain 100");
    }

    // Test 13: Parse invalid input
    TEST_CASE("Parse invalid input");
    {
        auto r1 = IntervalUnion::parse("");
        EXPECT_TRUE(!r1.has_value(), "Empty string should return nullopt");

        auto r2 = IntervalUnion::parse("invalid");
        EXPECT_TRUE(!r2.has_value(), "Invalid string should return nullopt");

        auto r3 = IntervalUnion::parse("[3, ]");
        EXPECT_TRUE(!r3.has_value(), "Missing value should return nullopt");
    }

    // Test 14: Round-trip test
    TEST_CASE("Round-trip: to_string then parse");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-2))},
            Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::pos_inf()}
        };
        auto original = IntervalUnion(ivs);
        std::string s = original.to_string();
        auto parsed = IntervalUnion::parse(s);
        EXPECT_TRUE(parsed.has_value(), "Round-trip parse should succeed");
        // Verify equivalence by checking containment at key points
        EXPECT_TRUE(parsed->contains(-10.0) == original.contains(-10.0), "Round-trip: -10");
        EXPECT_TRUE(parsed->contains(-2.0) == original.contains(-2.0), "Round-trip: -2");
        EXPECT_TRUE(parsed->contains(0.0) == original.contains(0.0), "Round-trip: 0");
        EXPECT_TRUE(parsed->contains(3.0) == original.contains(3.0), "Round-trip: 3");
        EXPECT_TRUE(parsed->contains(10.0) == original.contains(10.0), "Round-trip: 10");
    }

    // Test 15: Parse negative numbers
    TEST_CASE("Parse interval with negative numbers [-5, -1]");
    {
        auto result = IntervalUnion::parse("[-5, -1]");
        EXPECT_TRUE(result.has_value(), "Parse [-5, -1] should succeed");
        EXPECT_TRUE(result->contains(-3.0), "Should contain -3");
        EXPECT_TRUE(result->contains(-5.0), "Should contain -5");
        EXPECT_TRUE(result->contains(-1.0), "Should contain -1");
        EXPECT_TRUE(!result->contains(0.0), "Should not contain 0");
    }

    // Test 16: Parse decimal numbers
    TEST_CASE("Parse interval with decimals (1.5, 3.7)");
    {
        auto result = IntervalUnion::parse("(1.5, 3.7)");
        EXPECT_TRUE(result.has_value(), "Parse (1.5, 3.7) should succeed");
        EXPECT_TRUE(result->contains(2.0), "Should contain 2.0");
        EXPECT_TRUE(!result->contains(1.5), "Should not contain 1.5 (open)");
        EXPECT_TRUE(!result->contains(3.7), "Should not contain 3.7 (open)");
    }

    return TEST_REPORT();
}
