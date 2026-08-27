#include "test_common.hpp"
#include "interval.hpp"
#include <random>
#include <sstream>
#include <cmath>
#include <limits>

using namespace lamina;

int main() {

    TEST_CASE("Empty interval");
    {
        auto empty_iv = Interval::empty();
        EXPECT_TRUE(empty_iv.is_empty(), "Interval::empty().is_empty() == true");
        EXPECT_TRUE(!empty_iv.contains(0), "Interval::empty().contains(0) == false");
        EXPECT_TRUE(!empty_iv.contains(-100), "Interval::empty().contains(-100) == false");
        EXPECT_TRUE(!empty_iv.contains(100), "Interval::empty().contains(100) == false");
    }

    TEST_CASE("Single point interval [3, 3]");
    {
        auto pt = Interval::point(SymbolicExpr::number(3));
        EXPECT_TRUE(!pt.is_empty(), "Point interval is not empty");
        EXPECT_TRUE(pt.contains(3.0), "Interval::point(3).contains(3) == true");
        EXPECT_TRUE(!pt.contains(3.1), "Interval::point(3).contains(3.1) == false");
        EXPECT_TRUE(!pt.contains(2.9), "Interval::point(3).contains(2.9) == false");
        EXPECT_TRUE(!pt.contains(0), "Interval::point(3).contains(0) == false");
    }

    TEST_CASE("Entire line (-inf, +inf)");
    {
        auto entire = Interval::entire_line();
        EXPECT_TRUE(entire.is_entire_line(), "Interval::entire_line().is_entire_line() == true");
        EXPECT_TRUE(!entire.is_empty(), "Entire line is not empty");
        EXPECT_TRUE(entire.contains(0), "Entire line contains 0");
        EXPECT_TRUE(entire.contains(-1e15), "Entire line contains -1e15");
        EXPECT_TRUE(entire.contains(1e15), "Entire line contains 1e15");
        EXPECT_TRUE(entire.contains(-999.5), "Entire line contains -999.5");
    }

    TEST_CASE("Infinity endpoints always open");
    {
        auto neg = Endpoint::neg_inf();
        EXPECT_TRUE(neg.is_open, "neg_inf().is_open == true");
        EXPECT_TRUE(neg.is_neg_infinity, "neg_inf().is_neg_infinity == true");

        auto pos = Endpoint::pos_inf();
        EXPECT_TRUE(pos.is_open, "pos_inf().is_open == true");
        EXPECT_TRUE(pos.is_pos_infinity, "pos_inf().is_pos_infinity == true");
    }

    TEST_CASE("Contains at open endpoints: (1, 5)");
    {
        auto iv = Interval{
            Endpoint::open(SymbolicExpr::number(1)),
            Endpoint::open(SymbolicExpr::number(5))
        };
        EXPECT_TRUE(!iv.contains(1.0), "(1, 5) should NOT contain 1");
        EXPECT_TRUE(!iv.contains(5.0), "(1, 5) should NOT contain 5");
        EXPECT_TRUE(iv.contains(3.0), "(1, 5) should contain 3");
        EXPECT_TRUE(iv.contains(1.001), "(1, 5) should contain 1.001");
        EXPECT_TRUE(iv.contains(4.999), "(1, 5) should contain 4.999");
        EXPECT_TRUE(!iv.contains(0.999), "(1, 5) should NOT contain 0.999");
        EXPECT_TRUE(!iv.contains(5.001), "(1, 5) should NOT contain 5.001");
    }

    TEST_CASE("Contains at closed endpoints: [1, 5]");
    {
        auto iv = Interval{
            Endpoint::closed(SymbolicExpr::number(1)),
            Endpoint::closed(SymbolicExpr::number(5))
        };
        EXPECT_TRUE(iv.contains(1.0), "[1, 5] should contain 1");
        EXPECT_TRUE(iv.contains(5.0), "[1, 5] should contain 5");
        EXPECT_TRUE(iv.contains(3.0), "[1, 5] should contain 3");
        EXPECT_TRUE(!iv.contains(0.999), "[1, 5] should NOT contain 0.999");
        EXPECT_TRUE(!iv.contains(5.001), "[1, 5] should NOT contain 5.001");
    }

    TEST_CASE("Normalize merges overlapping: [1,5] ∪ [3,7] → [1,7]");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(5))},
            Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::closed(SymbolicExpr::number(7))}
        };
        auto u = IntervalUnion(ivs);
        EXPECT_TRUE(u.intervals().size() == 1, "Overlapping intervals should merge into 1");
        EXPECT_TRUE(u.contains(1.0), "Merged [1,7] contains 1");
        EXPECT_TRUE(u.contains(7.0), "Merged [1,7] contains 7");
        EXPECT_TRUE(u.contains(4.0), "Merged [1,7] contains 4");
        EXPECT_TRUE(!u.contains(0.5), "Merged [1,7] does not contain 0.5");
        EXPECT_TRUE(!u.contains(7.5), "Merged [1,7] does not contain 7.5");
    }

    TEST_CASE("Normalize merges adjacent: [1,3] ∪ [3,5] → [1,5]");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(3))},
            Interval{Endpoint::closed(SymbolicExpr::number(3)), Endpoint::closed(SymbolicExpr::number(5))}
        };
        auto u = IntervalUnion(ivs);
        EXPECT_TRUE(u.intervals().size() == 1, "Adjacent intervals [1,3]∪[3,5] should merge into 1");
        EXPECT_TRUE(u.contains(1.0), "Merged [1,5] contains 1");
        EXPECT_TRUE(u.contains(3.0), "Merged [1,5] contains 3");
        EXPECT_TRUE(u.contains(5.0), "Merged [1,5] contains 5");
    }

    TEST_CASE("Normalize doesn't merge disjoint: [1,2] ∪ [4,5]");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::closed(SymbolicExpr::number(1)), Endpoint::closed(SymbolicExpr::number(2))},
            Interval{Endpoint::closed(SymbolicExpr::number(4)), Endpoint::closed(SymbolicExpr::number(5))}
        };
        auto u = IntervalUnion(ivs);
        EXPECT_TRUE(u.intervals().size() == 2, "Disjoint intervals [1,2]∪[4,5] should stay as 2");
        EXPECT_TRUE(u.contains(1.5), "Contains 1.5 (in first interval)");
        EXPECT_TRUE(u.contains(4.5), "Contains 4.5 (in second interval)");
        EXPECT_TRUE(!u.contains(3.0), "Does not contain 3.0 (in gap)");
    }

    TEST_CASE("Parse invalid strings return nullopt");
    {
        auto r1 = IntervalUnion::parse("");
        EXPECT_TRUE(!r1.has_value(), "Empty string \"\" → nullopt");

        auto r2 = IntervalUnion::parse("invalid");
        EXPECT_TRUE(!r2.has_value(), "\"invalid\" → nullopt");

        auto r3 = IntervalUnion::parse("[3, ]");
        EXPECT_TRUE(!r3.has_value(), "\"[3, ]\" → nullopt");

        auto r4 = IntervalUnion::parse("abc");
        EXPECT_TRUE(!r4.has_value(), "\"abc\" → nullopt");
    }

    TEST_CASE("Parse empty set ∅");
    {
        auto result = IntervalUnion::parse("\xe2\x88\x85");
        EXPECT_TRUE(result.has_value(), "Parse ∅ should succeed");
        EXPECT_TRUE(result->is_empty(), "Parsed ∅ should be empty IntervalUnion");
        EXPECT_TRUE(!result->contains(0), "Parsed ∅ should not contain 0");
    }

    TEST_CASE("Parse entire line (-∞, +∞)");
    {
        auto result = IntervalUnion::parse("(-\xe2\x88\x9e, +\xe2\x88\x9e)");
        EXPECT_TRUE(result.has_value(), "Parse (-∞, +∞) should succeed");
        EXPECT_TRUE(result->is_entire_line(), "Parsed (-∞, +∞) should be entire line");
        EXPECT_TRUE(result->contains(0), "Entire line contains 0");
        EXPECT_TRUE(result->contains(-1e10), "Entire line contains -1e10");
        EXPECT_TRUE(result->contains(1e10), "Entire line contains 1e10");
    }

    TEST_CASE("IntervalUnion Serialization Round-trip");
    {
        std::mt19937 rng(42);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLE_POINTS = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> val_dist(-20, 20);
        std::uniform_int_distribution<int> bool_dist(0, 1);
        std::uniform_int_distribution<int> inf_dist(0, 5);
        std::uniform_int_distribution<int> count_dist(1, 3);
        std::uniform_real_distribution<double> sample_dist(-100.0, 100.0);

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int count = count_dist(rng);
            std::vector<Interval> intervals;

            for (int c = 0; c < count; ++c) {
                Endpoint lower, upper;

                if (inf_dist(rng) == 0) {
                    lower = Endpoint::neg_inf();
                } else {
                    int val = val_dist(rng);
                    lower = bool_dist(rng) ? Endpoint::open(SymbolicExpr::number(val))
                                           : Endpoint::closed(SymbolicExpr::number(val));
                }

                if (inf_dist(rng) == 0) {
                    upper = Endpoint::pos_inf();
                } else {
                    int val = val_dist(rng);

                    if (!lower.is_neg_infinity && lower.value) {
                        int lower_val = static_cast<int>(lower.value->to_numeric());
                        if (val < lower_val) val = lower_val;
                    }
                    upper = bool_dist(rng) ? Endpoint::open(SymbolicExpr::number(val))
                                           : Endpoint::closed(SymbolicExpr::number(val));
                }

                intervals.push_back(Interval{lower, upper});
            }

            IntervalUnion original(intervals);

            std::string str = original.to_string();

            auto parsed = IntervalUnion::parse(str);

            if (!parsed.has_value()) {
                std::ostringstream oss;
                oss << "Iter " << iter << ": parse() failed for: " << str;
                EXPECT_TRUE(false, oss.str());
                continue;
            }

            bool equivalent = true;
            for (int s = 0; s < NUM_SAMPLE_POINTS; ++s) {
                double point = sample_dist(rng);
                bool orig_contains = original.contains(point);
                bool parsed_contains = parsed->contains(point);

                if (orig_contains != parsed_contains) {
                    std::ostringstream oss;
                    oss << "Iter " << iter << ": mismatch at " << point
                        << " (orig=" << orig_contains << ", parsed=" << parsed_contains
                        << ") str: " << str;
                    EXPECT_TRUE(false, oss.str());
                    equivalent = false;
                    break;
                }
            }

            if (equivalent) {
                ++pass_count;
            }
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed serialization round-trip";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Set Operations Correctness (intersect/unite/complement sampling)");
    {
        std::mt19937 rng(100);
        std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
        std::uniform_int_distribution<int> open_dist(0, 1);
        std::uniform_int_distribution<int> count_dist(1, 4);
        std::uniform_real_distribution<double> sample_dist(-100.0, 100.0);

        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 100;
        int pass_count = 0;

        auto gen_interval = [&]() -> Interval {
            double a = std::round(val_dist(rng) * 10.0) / 10.0;
            double b = std::round(val_dist(rng) * 10.0) / 10.0;
            if (a > b) std::swap(a, b);
            Endpoint lower = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(a))
                                            : Endpoint::closed(SymbolicExpr::number(a));
            Endpoint upper = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(b))
                                            : Endpoint::closed(SymbolicExpr::number(b));
            return Interval{lower, upper};
        };

        auto gen_union = [&]() -> IntervalUnion {
            int count = count_dist(rng);
            std::vector<Interval> intervals;
            for (int i = 0; i < count; ++i) {
                intervals.push_back(gen_interval());
            }
            return IntervalUnion(std::move(intervals));
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
            IntervalUnion A = gen_union();
            IntervalUnion B = gen_union();

            IntervalUnion A_intersect_B = A.intersect(B);
            IntervalUnion A_unite_B = A.unite(B);
            IntervalUnion A_complement = A.complement();

            bool iter_passed = true;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double x = sample_dist(rng);

                bool in_intersect = A_intersect_B.contains(x);
                bool expected_intersect = A.contains(x) && B.contains(x);
                if (in_intersect != expected_intersect) {
                    std::ostringstream oss;
                    oss << "(intersect) failed: iter=" << iter << " x=" << x
                        << " in_intersect=" << in_intersect << " expected=" << expected_intersect;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }

                bool in_unite = A_unite_B.contains(x);
                bool expected_unite = A.contains(x) || B.contains(x);
                if (in_unite != expected_unite) {
                    std::ostringstream oss;
                    oss << "(unite) failed: iter=" << iter << " x=" << x
                        << " in_unite=" << in_unite << " expected=" << expected_unite;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }

                bool in_complement = A_complement.contains(x);
                bool expected_complement = !A.contains(x);
                if (in_complement != expected_complement) {
                    std::ostringstream oss;
                    oss << "(complement) failed: iter=" << iter << " x=" << x
                        << " in_complement=" << in_complement << " expected=" << expected_complement;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed set operations correctness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("De Morgan's Law (A ∩ complement(A) = ∅, A ∪ complement(A) = ℝ)");
    {
        std::mt19937 rng(200);
        std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
        std::uniform_int_distribution<int> open_dist(0, 1);
        std::uniform_int_distribution<int> count_dist(1, 4);

        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        auto gen_interval = [&]() -> Interval {
            double a = std::round(val_dist(rng) * 10.0) / 10.0;
            double b = std::round(val_dist(rng) * 10.0) / 10.0;
            if (a > b) std::swap(a, b);
            Endpoint lower = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(a))
                                            : Endpoint::closed(SymbolicExpr::number(a));
            Endpoint upper = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(b))
                                            : Endpoint::closed(SymbolicExpr::number(b));
            return Interval{lower, upper};
        };

        auto gen_union = [&]() -> IntervalUnion {
            int count = count_dist(rng);
            std::vector<Interval> intervals;
            for (int i = 0; i < count; ++i) {
                intervals.push_back(gen_interval());
            }
            return IntervalUnion(std::move(intervals));
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
            IntervalUnion A = gen_union();
            IntervalUnion A_comp = A.complement();

            IntervalUnion intersection = A.intersect(A_comp);
            if (!intersection.is_empty()) {
                std::ostringstream oss;
                oss << "(A ∩ complement(A) != ∅): iter=" << iter
                    << " got " << intersection.intervals().size() << " intervals: "
                    << intersection.to_string();
                EXPECT_TRUE(false, oss.str());
                continue;
            }

            IntervalUnion union_result = A.unite(A_comp);
            if (!union_result.is_entire_line()) {
                std::ostringstream oss;
                oss << "(A ∪ complement(A) != ℝ): iter=" << iter
                    << " got: " << union_result.to_string();
                EXPECT_TRUE(false, oss.str());
                continue;
            }

            ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed De Morgan's Law";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("IntervalUnion Invariant (sorted by lower bound, pairwise disjoint)");
    {
        std::mt19937 rng(300);
        std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
        std::uniform_int_distribution<int> open_dist(0, 1);
        std::uniform_int_distribution<int> count_dist(1, 4);

        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        auto get_ep_value = [](const Endpoint& ep) -> double {
            if (ep.is_neg_infinity) return -std::numeric_limits<double>::infinity();
            if (ep.is_pos_infinity) return std::numeric_limits<double>::infinity();
            if (ep.value) return ep.value->to_numeric();
            return 0.0;
        };

        auto gen_interval = [&]() -> Interval {
            double a = std::round(val_dist(rng) * 10.0) / 10.0;
            double b = std::round(val_dist(rng) * 10.0) / 10.0;
            if (a > b) std::swap(a, b);
            Endpoint lower = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(a))
                                            : Endpoint::closed(SymbolicExpr::number(a));
            Endpoint upper = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(b))
                                            : Endpoint::closed(SymbolicExpr::number(b));
            return Interval{lower, upper};
        };

        auto gen_union = [&]() -> IntervalUnion {
            int count = count_dist(rng);
            std::vector<Interval> intervals;
            for (int i = 0; i < count; ++i) {
                intervals.push_back(gen_interval());
            }
            return IntervalUnion(std::move(intervals));
        };

        auto check_invariant = [&](const IntervalUnion& u, const std::string& label, int iter) -> bool {
            const auto& ivs = u.intervals();
            for (size_t i = 1; i < ivs.size(); ++i) {
                double prev_lower = get_ep_value(ivs[i - 1].lower);
                double curr_lower = get_ep_value(ivs[i].lower);
                if (prev_lower > curr_lower) {
                    std::ostringstream oss;
                    oss << "(" << label << " not sorted): iter=" << iter
                        << " ivs[" << (i-1) << "].lower=" << prev_lower
                        << " > ivs[" << i << "].lower=" << curr_lower;
                    EXPECT_TRUE(false, oss.str());
                    return false;
                }

                double prev_upper = get_ep_value(ivs[i - 1].upper);
                if (prev_upper > curr_lower) {
                    std::ostringstream oss;
                    oss << "(" << label << " not disjoint): iter=" << iter
                        << " ivs[" << (i-1) << "].upper=" << prev_upper
                        << " > ivs[" << i << "].lower=" << curr_lower;
                    EXPECT_TRUE(false, oss.str());
                    return false;
                }

                if (prev_upper == curr_lower &&
                    !ivs[i - 1].upper.is_open && !ivs[i].lower.is_open) {
                    std::ostringstream oss;
                    oss << "(" << label << " adjacent not merged): iter=" << iter
                        << " both closed at " << prev_upper;
                    EXPECT_TRUE(false, oss.str());
                    return false;
                }
            }
            return true;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
            IntervalUnion A = gen_union();
            IntervalUnion B = gen_union();

            bool ok = true;

            ok = ok && check_invariant(A, "A", iter);

            ok = ok && check_invariant(B, "B", iter);

            ok = ok && check_invariant(A.intersect(B), "A∩B", iter);

            ok = ok && check_invariant(A.unite(B), "A∪B", iter);

            ok = ok && check_invariant(A.complement(), "complement(A)", iter);

            if (ok) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed IntervalUnion invariant";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("to_expr: empty set returns nullptr");
    {
        auto empty_u = IntervalUnion::empty();
        auto expr = empty_u.to_expr("x");
        EXPECT_TRUE(expr == nullptr, "Empty IntervalUnion to_expr() returns nullptr");
    }

    TEST_CASE("to_expr: entire line returns nullptr");
    {
        auto entire = IntervalUnion::entire_line();
        auto expr = entire.to_expr("x");
        EXPECT_TRUE(expr == nullptr, "Entire line to_expr() returns nullptr (no constraint)");
    }

    TEST_CASE("to_expr: single interval (2, 5) -> x > 2 And x < 5");
    {
        auto iv = Interval{
            Endpoint::open(SymbolicExpr::number(2)),
            Endpoint::open(SymbolicExpr::number(5))
        };
        auto u = IntervalUnion::from_single(iv);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for (2, 5) should not be nullptr");
        std::string str = expr->to_string();

        EXPECT_TRUE(str.find("and") != std::string::npos, "to_expr (2,5) should contain 'and'");
        EXPECT_TRUE(str.find(">") != std::string::npos, "to_expr (2,5) should contain '>'");
        EXPECT_TRUE(str.find("<") != std::string::npos, "to_expr (2,5) should contain '<'");
    }

    TEST_CASE("to_expr: single interval (3, +inf) -> x > 3");
    {
        auto iv = Interval{
            Endpoint::open(SymbolicExpr::number(3)),
            Endpoint::pos_inf()
        };
        auto u = IntervalUnion::from_single(iv);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for (3, +inf) should not be nullptr");
        std::string str = expr->to_string();
        EXPECT_TRUE(str.find(">") != std::string::npos, "to_expr (3,+inf) should contain '>'");

        EXPECT_TRUE(str.find("and") == std::string::npos, "to_expr (3,+inf) should NOT contain 'and'");
    }

    TEST_CASE("to_expr: single interval (-inf, 4] -> x <= 4");
    {
        auto iv = Interval{
            Endpoint::neg_inf(),
            Endpoint::closed(SymbolicExpr::number(4))
        };
        auto u = IntervalUnion::from_single(iv);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for (-inf, 4] should not be nullptr");
        std::string str = expr->to_string();
        EXPECT_TRUE(str.find("<=") != std::string::npos, "to_expr (-inf,4] should contain '<='");
        EXPECT_TRUE(str.find("and") == std::string::npos, "to_expr (-inf,4] should NOT contain 'and'");
    }

    TEST_CASE("to_expr: single interval [1, 3] -> x >= 1 And x <= 3");
    {
        auto iv = Interval{
            Endpoint::closed(SymbolicExpr::number(1)),
            Endpoint::closed(SymbolicExpr::number(3))
        };
        auto u = IntervalUnion::from_single(iv);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for [1, 3] should not be nullptr");
        std::string str = expr->to_string();
        EXPECT_TRUE(str.find(">=") != std::string::npos, "to_expr [1,3] should contain '>='");
        EXPECT_TRUE(str.find("<=") != std::string::npos, "to_expr [1,3] should contain '<='");
        EXPECT_TRUE(str.find("and") != std::string::npos, "to_expr [1,3] should contain 'and'");
    }

    TEST_CASE("to_expr: two intervals (-inf, -2) ∪ (3, +inf) -> Or");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-2))},
            Interval{Endpoint::open(SymbolicExpr::number(3)), Endpoint::pos_inf()}
        };
        auto u = IntervalUnion(ivs);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for (-inf,-2)∪(3,+inf) should not be nullptr");
        std::string str = expr->to_string();
        EXPECT_TRUE(str.find("or") != std::string::npos, "to_expr multi-interval should contain 'or'");
        EXPECT_TRUE(str.find("<") != std::string::npos, "to_expr should contain '<'");
        EXPECT_TRUE(str.find(">") != std::string::npos, "to_expr should contain '>'");
    }

    TEST_CASE("to_expr: three intervals uses nested Or");
    {
        std::vector<Interval> ivs = {
            Interval{Endpoint::neg_inf(), Endpoint::open(SymbolicExpr::number(-5))},
            Interval{Endpoint::open(SymbolicExpr::number(-1)), Endpoint::open(SymbolicExpr::number(1))},
            Interval{Endpoint::open(SymbolicExpr::number(5)), Endpoint::pos_inf()}
        };
        auto u = IntervalUnion(ivs);
        auto expr = u.to_expr("x");
        EXPECT_TRUE(expr != nullptr, "to_expr() for 3 intervals should not be nullptr");
        std::string str = expr->to_string();
        EXPECT_TRUE(str.find("or") != std::string::npos, "to_expr 3 intervals should contain 'or'");
    }

    return TEST_REPORT();
}
