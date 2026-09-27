#include "test_interval_support.hpp"

TEST(Interval, EmptyInterval) {
    auto empty_iv = Interval::empty();
    EXPECT_TRUE((empty_iv.is_empty())) << "Interval::empty().is_empty() == true";
    EXPECT_TRUE((!empty_iv.contains(0))) << "Interval::empty().contains(0) == false";
    EXPECT_TRUE((!empty_iv.contains(-100))) << "Interval::empty().contains(-100) == false";
    EXPECT_TRUE((!empty_iv.contains(100))) << "Interval::empty().contains(100) == false";
}

TEST(Interval, SinglePointInterval33) {
    auto pt = Interval::point(SymbolicExpr::number(3));
    EXPECT_TRUE((!pt.is_empty())) << "Point interval is not empty";
    EXPECT_TRUE((pt.contains(3.0))) << "Interval::point(3).contains(3) == true";
    EXPECT_TRUE((!pt.contains(3.1))) << "Interval::point(3).contains(3.1) == false";
    EXPECT_TRUE((!pt.contains(2.9))) << "Interval::point(3).contains(2.9) == false";
    EXPECT_TRUE((!pt.contains(0))) << "Interval::point(3).contains(0) == false";
}

TEST(Interval, EntireLineInfInf) {
    auto entire = Interval::entire_line();
    EXPECT_TRUE((entire.is_entire_line())) << "Interval::entire_line().is_entire_line() == true";
    EXPECT_TRUE((!entire.is_empty())) << "Entire line is not empty";
    EXPECT_TRUE((entire.contains(0))) << "Entire line contains 0";
    EXPECT_TRUE((entire.contains(-1e15))) << "Entire line contains -1e15";
    EXPECT_TRUE((entire.contains(1e15))) << "Entire line contains 1e15";
    EXPECT_TRUE((entire.contains(-999.5))) << "Entire line contains -999.5";
}

TEST(Interval, InfinityEndpointsAlwaysOpen) {
    auto neg = Endpoint::neg_inf();
    EXPECT_TRUE((neg.is_open)) << "neg_inf().is_open == true";
    EXPECT_TRUE((neg.is_neg_infinity)) << "neg_inf().is_neg_infinity == true";

    auto pos = Endpoint::pos_inf();
    EXPECT_TRUE((pos.is_open)) << "pos_inf().is_open == true";
    EXPECT_TRUE((pos.is_pos_infinity)) << "pos_inf().is_pos_infinity == true";
}

TEST(Interval, ContainsAtOpenEndpoints15) {
    auto iv = Interval{
        Endpoint::open(SymbolicExpr::number(1)),
        Endpoint::open(SymbolicExpr::number(5))};
    EXPECT_TRUE((!iv.contains(1.0))) << "(1, 5) should NOT contain 1";
    EXPECT_TRUE((!iv.contains(5.0))) << "(1, 5) should NOT contain 5";
    EXPECT_TRUE((iv.contains(3.0))) << "(1, 5) should contain 3";
    EXPECT_TRUE((iv.contains(1.001))) << "(1, 5) should contain 1.001";
    EXPECT_TRUE((iv.contains(4.999))) << "(1, 5) should contain 4.999";
    EXPECT_TRUE((!iv.contains(0.999))) << "(1, 5) should NOT contain 0.999";
    EXPECT_TRUE((!iv.contains(5.001))) << "(1, 5) should NOT contain 5.001";
}

TEST(Interval, ContainsAtClosedEndpoints15) {
    auto iv = Interval{
        Endpoint::closed(SymbolicExpr::number(1)),
        Endpoint::closed(SymbolicExpr::number(5))};
    EXPECT_TRUE((iv.contains(1.0))) << "[1, 5] should contain 1";
    EXPECT_TRUE((iv.contains(5.0))) << "[1, 5] should contain 5";
    EXPECT_TRUE((iv.contains(3.0))) << "[1, 5] should contain 3";
    EXPECT_TRUE((!iv.contains(0.999))) << "[1, 5] should NOT contain 0.999";
    EXPECT_TRUE((!iv.contains(5.001))) << "[1, 5] should NOT contain 5.001";
}

TEST(Interval, NormalizeMergesOverlapping153717) {
    std::vector<Interval> ivs = {
        Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(5))},
        Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::closed(SymbolicExpr::number(7))}};
    auto u = IntervalUnion(ivs);
    EXPECT_TRUE((u.intervals().size() == 1)) << "Overlapping intervals should merge into 1";
    EXPECT_TRUE((u.contains(1.0))) << "Merged [1,7] contains 1";
    EXPECT_TRUE((u.contains(7.0))) << "Merged [1,7] contains 7";
    EXPECT_TRUE((u.contains(4.0))) << "Merged [1,7] contains 4";
    EXPECT_TRUE((!u.contains(0.5))) << "Merged [1,7] does not contain 0.5";
    EXPECT_TRUE((!u.contains(7.5))) << "Merged [1,7] does not contain 7.5";
}

TEST(Interval, NormalizeMergesAdjacent133515) {
    std::vector<Interval> ivs = {
        Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(3))},
        Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::closed(SymbolicExpr::number(5))}};
    auto u = IntervalUnion(ivs);
    EXPECT_TRUE((u.intervals().size() == 1)) << "Adjacent intervals [1,3]∪[3,5] should merge into 1";
    EXPECT_TRUE((u.contains(1.0))) << "Merged [1,5] contains 1";
    EXPECT_TRUE((u.contains(3.0))) << "Merged [1,5] contains 3";
    EXPECT_TRUE((u.contains(5.0))) << "Merged [1,5] contains 5";
}

TEST(Interval, NormalizeDoesnTMergeDisjoint1245) {
    std::vector<Interval> ivs = {
        Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(2))},
        Interval{Endpoint::closed(SymbolicExpr::number(4)), Endpoint::closed(SymbolicExpr::number(5))}};
    auto u = IntervalUnion(ivs);
    EXPECT_TRUE((u.intervals().size() == 2)) << "Disjoint intervals [1,2]∪[4,5] should stay as 2";
    EXPECT_TRUE((u.contains(1.5))) << "Contains 1.5 (in first interval)";
    EXPECT_TRUE((u.contains(4.5))) << "Contains 4.5 (in second interval)";
    EXPECT_TRUE((!u.contains(3.0))) << "Does not contain 3.0 (in gap)";
}

TEST(Interval, ParseInvalidStringsReturnNullopt) {
    auto r1 = IntervalUnion::parse("");
    EXPECT_TRUE((!r1.has_value())) << "Empty string \"\" → nullopt";

    auto r2 = IntervalUnion::parse("invalid");
    EXPECT_TRUE((!r2.has_value())) << "\"invalid\" → nullopt";

    auto r3 = IntervalUnion::parse("[3, ]");
    EXPECT_TRUE((!r3.has_value())) << "\"[3, ]\" → nullopt";

    auto r4 = IntervalUnion::parse("abc");
    EXPECT_TRUE((!r4.has_value())) << "\"abc\" → nullopt";
}

TEST(Interval, ParseEmptySet) {
    auto result = IntervalUnion::parse("\xe2\x88\x85");
    EXPECT_TRUE((result.has_value())) << "Parse ∅ should succeed";
    EXPECT_TRUE((result->is_empty())) << "Parsed ∅ should be empty IntervalUnion";
    EXPECT_TRUE((!result->contains(0))) << "Parsed ∅ should not contain 0";
}

TEST(Interval, ParseEntireLine) {
    auto result = IntervalUnion::parse("(-\xe2\x88\x9e, +\xe2\x88\x9e)");
    EXPECT_TRUE((result.has_value())) << "Parse (-∞, +∞) should succeed";
    EXPECT_TRUE((result->is_entire_line())) << "Parsed (-∞, +∞) should be entire line";
    EXPECT_TRUE((result->contains(0))) << "Entire line contains 0";
    EXPECT_TRUE((result->contains(-1e10))) << "Entire line contains -1e10";
    EXPECT_TRUE((result->contains(1e10))) << "Entire line contains 1e10";
}
