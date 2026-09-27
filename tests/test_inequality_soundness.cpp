#include "test_inequality_support.hpp"

bool sample_inside(const IntervalUnion &iu, std::mt19937 &rng, double &out) {
    const auto &intervals = iu.intervals();
    if (intervals.empty()) {
        return false;
    }

    std::uniform_int_distribution<size_t> idx_dist(0, intervals.size() - 1);
    size_t idx = idx_dist(rng);
    const auto &iv = intervals[idx];

    double lo = -1000.0;
    double hi = 1000.0;

    if (!iv.lower.is_neg_infinity && iv.lower.value) {
        lo = iv.lower.value->to_numeric();
    }
    if (!iv.upper.is_pos_infinity && iv.upper.value) {
        hi = iv.upper.value->to_numeric();
    }

    double epsilon = 1e-6;
    if (iv.lower.is_open && !iv.lower.is_neg_infinity) {
        lo += epsilon;
    }
    if (iv.upper.is_open && !iv.upper.is_pos_infinity) {
        hi -= epsilon;
    }

    if (lo >= hi) {
        out = (lo + hi) / 2.0;
        return true;
    }

    std::uniform_real_distribution<double> point_dist(lo, hi);
    out = point_dist(rng);
    return true;
}

bool check_soundness_samples(
    const std::shared_ptr<SymbolicExpr> &poly,
    const std::vector<int> &coeffs,
    InequalityType type,
    const IntervalUnion &solution,
    std::mt19937 &rng,
    int iter) {
    for (int sample = 0; sample < 100; ++sample) {
        double point;
        if (!sample_inside(solution, rng, point)) {
            continue;
        }
        double value = eval_poly(poly, point);
        if (satisfies(value, type)) {
            continue;
        }
        std::ostringstream message;
        message << "FAIL: iter=" << iter
                << " point=" << point << " poly_value=" << value
                << " type=" << static_cast<int>(type)
                << " poly=" << format_coefficients(coeffs);
        ADD_FAILURE() << message.str();
        return false;
    }
    return true;
}

TEST(InequalitySoundness, SolutionSoundness) {
    std::mt19937 rng(42);
    const int num_iterations = 100;
    std::uniform_int_distribution<int> type_dist(0, 3);
    for (int iter = 0; iter < num_iterations; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        auto coeffs = generate_coefficients(rng);
        auto poly = build_poly(coeffs);
        InequalityType type = inequality_types[type_dist(rng)];
        auto solution = InequalitySolver::solve_inequality(poly, type, "x");
        if (solution.is_empty()) {
            continue;
        }
        EXPECT_TRUE(check_soundness_samples(poly, coeffs, type, solution, rng, iter));
    }
}
