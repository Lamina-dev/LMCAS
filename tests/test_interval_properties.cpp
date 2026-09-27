#include "test_interval_support.hpp"

static Interval generate_sample_interval(std::mt19937 &rng,
                                         std::uniform_real_distribution<double> &val_dist,
                                         std::uniform_int_distribution<int> &open_dist) {

    double a = std::round(val_dist(rng) * 10.0) / 10.0;
    double b = std::round(val_dist(rng) * 10.0) / 10.0;
    if (a > b)
        std::swap(a, b);
    Endpoint lower = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(a))
                                    : Endpoint::closed(SymbolicExpr::number(a));
    Endpoint upper = open_dist(rng) ? Endpoint::open(SymbolicExpr::number(b))
                                    : Endpoint::closed(SymbolicExpr::number(b));
    return Interval{lower, upper};
}

static IntervalUnion generate_sample_union(std::mt19937 &rng,
                                           std::uniform_real_distribution<double> &val_dist,
                                           std::uniform_int_distribution<int> &open_dist,
                                           std::uniform_int_distribution<int> &count_dist) {
    int count = count_dist(rng);
    std::vector<Interval> intervals;
    for (int i = 0; i < count; ++i) {
        intervals.push_back(generate_sample_interval(rng, val_dist, open_dist));
    }
    return IntervalUnion(std::move(intervals));
}

static double interval_endpoint_value(const Endpoint &ep) {

    if (ep.is_neg_infinity) {
        return -std::numeric_limits<double>::infinity();
    }
    if (ep.is_pos_infinity) {
        return std::numeric_limits<double>::infinity();
    }
    if (ep.value) {
        return ep.value->to_numeric();
    }
    return 0.0;
}

static bool check_union_invariant(const IntervalUnion &u, const std::string &label, int iter) {

    const auto &ivs = u.intervals();
    for (size_t i = 1; i < ivs.size(); ++i) {
        double prev_lower = interval_endpoint_value(ivs[i - 1].lower);
        double curr_lower = interval_endpoint_value(ivs[i].lower);
        if (prev_lower > curr_lower) {
            std::ostringstream oss;
            oss << "(" << label << " not sorted): iter=" << iter
                << " ivs[" << (i - 1) << "].lower=" << prev_lower
                << " > ivs[" << i << "].lower=" << curr_lower;
            ADD_FAILURE() << oss.str();
            return false;
        }

        double prev_upper = interval_endpoint_value(ivs[i - 1].upper);
        if (prev_upper > curr_lower) {
            std::ostringstream oss;
            oss << "(" << label << " not disjoint): iter=" << iter
                << " ivs[" << (i - 1) << "].upper=" << prev_upper
                << " > ivs[" << i << "].lower=" << curr_lower;
            ADD_FAILURE() << oss.str();
            return false;
        }

        if (prev_upper == curr_lower &&
            !ivs[i - 1].upper.is_open && !ivs[i].lower.is_open) {
            std::ostringstream oss;
            oss << "(" << label << " adjacent not merged): iter=" << iter
                << " both closed at " << prev_upper;
            ADD_FAILURE() << oss.str();
            return false;
        }
    }
    return true;
}

static IntervalUnion generate_serialization_union(std::mt19937 &rng,
                                                  std::uniform_int_distribution<int> &val_dist,
                                                  std::uniform_int_distribution<int> &bool_dist,
                                                  std::uniform_int_distribution<int> &inf_dist,
                                                  std::uniform_int_distribution<int> &count_dist) {
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
                if (val < lower_val) {
                    val = lower_val;
                }
            }
            upper = bool_dist(rng) ? Endpoint::open(SymbolicExpr::number(val))
                                   : Endpoint::closed(SymbolicExpr::number(val));
        }

        intervals.push_back(Interval{lower, upper});
    }

    return IntervalUnion(intervals);
}

static bool serialization_samples_match(const IntervalUnion &original,
                                        const IntervalUnion &parsed, const std::string &str, int iter,
                                        std::mt19937 &rng, std::uniform_real_distribution<double> &sample_dist,
                                        int NUM_SAMPLE_POINTS) {
    bool equivalent = true;
    for (int s = 0; s < NUM_SAMPLE_POINTS; ++s) {
        double point = sample_dist(rng);
        bool orig_contains = original.contains(point);
        bool parsed_contains = parsed.contains(point);

        if (orig_contains != parsed_contains) {
            std::ostringstream oss;
            oss << "Iter " << iter << ": mismatch at " << point
                << " (orig=" << orig_contains << ", parsed=" << parsed_contains
                << ") str: " << str;
            ADD_FAILURE() << oss.str();
            equivalent = false;
            break;
        }
    }

    return equivalent;
}

static bool set_operation_samples_match(const IntervalUnion &A,
                                        const IntervalUnion &B, const IntervalUnion &A_intersect_B,
                                        const IntervalUnion &A_unite_B, const IntervalUnion &A_complement,
                                        int iter, std::mt19937 &rng,
                                        std::uniform_real_distribution<double> &sample_dist, int NUM_SAMPLES) {
    bool iter_passed = true;
    for (int s = 0; s < NUM_SAMPLES; ++s) {
        double x = sample_dist(rng);

        bool in_intersect = A_intersect_B.contains(x);
        bool expected_intersect = A.contains(x) && B.contains(x);
        if (in_intersect != expected_intersect) {
            std::ostringstream oss;
            oss << "(intersect) failed: iter=" << iter << " x=" << x
                << " in_intersect=" << in_intersect << " expected=" << expected_intersect;
            ADD_FAILURE() << oss.str();
            iter_passed = false;
            break;
        }

        bool in_unite = A_unite_B.contains(x);
        bool expected_unite = A.contains(x) || B.contains(x);
        if (in_unite != expected_unite) {
            std::ostringstream oss;
            oss << "(unite) failed: iter=" << iter << " x=" << x
                << " in_unite=" << in_unite << " expected=" << expected_unite;
            ADD_FAILURE() << oss.str();
            iter_passed = false;
            break;
        }

        bool in_complement = A_complement.contains(x);
        bool expected_complement = !A.contains(x);
        if (in_complement != expected_complement) {
            std::ostringstream oss;
            oss << "(complement) failed: iter=" << iter << " x=" << x
                << " in_complement=" << in_complement << " expected=" << expected_complement;
            ADD_FAILURE() << oss.str();
            iter_passed = false;
            break;
        }
    }

    return iter_passed;
}

TEST(IntervalProperties, IntervalunionSerializationRoundTrip) {
    std::mt19937 rng(42);
    const int NUM_ITERATIONS = 100;
    const int NUM_SAMPLE_POINTS = 100;

    std::uniform_int_distribution<int> val_dist(-20, 20);
    std::uniform_int_distribution<int> bool_dist(0, 1);
    std::uniform_int_distribution<int> inf_dist(0, 5);
    std::uniform_int_distribution<int> count_dist(1, 3);
    std::uniform_real_distribution<double> sample_dist(-100.0, 100.0);

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        IntervalUnion original = generate_serialization_union(
            rng, val_dist, bool_dist, inf_dist, count_dist);

        std::string str = original.to_string();

        auto parsed = IntervalUnion::parse(str);

        if (!parsed.has_value()) {
            std::ostringstream oss;
            oss << "Iter " << iter << ": parse() failed for: " << str;
            ADD_FAILURE() << oss.str();
            continue;
        }

        EXPECT_TRUE(serialization_samples_match(
            original, *parsed, str, iter, rng, sample_dist, NUM_SAMPLE_POINTS));
    }
}

TEST(IntervalProperties, SetOperationsCorrectnessIntersectUniteComplementSampling) {
    std::mt19937 rng(100);
    std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
    std::uniform_int_distribution<int> open_dist(0, 1);
    std::uniform_int_distribution<int> count_dist(1, 4);
    std::uniform_real_distribution<double> sample_dist(-100.0, 100.0);

    const int NUM_ITERATIONS = 100;
    const int NUM_SAMPLES = 100;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        IntervalUnion A = generate_sample_union(rng, val_dist, open_dist, count_dist);
        IntervalUnion B = generate_sample_union(rng, val_dist, open_dist, count_dist);

        IntervalUnion A_intersect_B = A.intersect(B);
        IntervalUnion A_unite_B = A.unite(B);
        IntervalUnion A_complement = A.complement();

        EXPECT_TRUE(set_operation_samples_match(A, B,
                                                A_intersect_B, A_unite_B, A_complement, iter, rng, sample_dist, NUM_SAMPLES));
    }
}

TEST(IntervalProperties, DeMorganSLawAComplementAAComplementA) {
    std::mt19937 rng(200);
    std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
    std::uniform_int_distribution<int> open_dist(0, 1);
    std::uniform_int_distribution<int> count_dist(1, 4);

    const int NUM_ITERATIONS = 100;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        IntervalUnion A = generate_sample_union(rng, val_dist, open_dist, count_dist);
        IntervalUnion A_comp = A.complement();

        IntervalUnion intersection = A.intersect(A_comp);
        EXPECT_TRUE(intersection.is_empty())
            << "(A ∩ complement(A) != ∅): iter=" << iter
            << " got " << intersection.intervals().size() << " intervals: "
            << intersection.to_string();
        if (!intersection.is_empty()) {
            continue;
        }

        IntervalUnion union_result = A.unite(A_comp);
        EXPECT_TRUE(union_result.is_entire_line())
            << "(A ∪ complement(A) != ℝ): iter=" << iter
            << " got: " << union_result.to_string();
        if (!union_result.is_entire_line()) {
            continue;
        }
    }
}

TEST(IntervalProperties, IntervalunionInvariantSortedByLowerBoundPairwiseDisjoint) {
    std::mt19937 rng(300);
    std::uniform_real_distribution<double> val_dist(-50.0, 50.0);
    std::uniform_int_distribution<int> open_dist(0, 1);
    std::uniform_int_distribution<int> count_dist(1, 4);

    const int NUM_ITERATIONS = 100;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        IntervalUnion A = generate_sample_union(rng, val_dist, open_dist, count_dist);
        IntervalUnion B = generate_sample_union(rng, val_dist, open_dist, count_dist);

        bool ok = check_union_invariant(A, "A", iter);

        if (ok) {
            ok = check_union_invariant(B, "B", iter);
        }

        if (ok) {
            ok = check_union_invariant(A.intersect(B), "A∩B", iter);
        }

        if (ok) {
            ok = check_union_invariant(A.unite(B), "A∪B", iter);
        }

        if (ok) {
            ok = check_union_invariant(A.complement(), "complement(A)", iter);
        }

        EXPECT_TRUE(ok) << "IntervalUnion invariant";
    }
}
