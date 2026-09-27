#include "test_inequality_support.hpp"

bool near_solution_boundary(const IntervalUnion &solution, double point) {
    for (const auto &interval : solution.intervals()) {
        if (!interval.lower.is_neg_infinity && interval.lower.value) {
            double boundary = interval.lower.value->to_numeric();
            if (std::abs(point - boundary) < 1e-8)
                return true;
        }
        if (!interval.upper.is_pos_infinity && interval.upper.value) {
            double boundary = interval.upper.value->to_numeric();
            if (std::abs(point - boundary) < 1e-8)
                return true;
        }
    }
    return false;
}

bool sample_outside(const IntervalUnion &solution, std::mt19937 &rng, double &out) {
    std::uniform_real_distribution<double> range_dist(-100.0, 100.0);
    for (int attempt = 0; attempt < 1000; ++attempt) {
        double candidate = range_dist(rng);
        if (solution.contains(candidate))
            continue;
        if (near_solution_boundary(solution, candidate))
            continue;
        out = candidate;
        return true;
    }
    return false;
}

bool check_completeness_samples(
    const std::shared_ptr<SymbolicExpr> &poly,
    const std::vector<int> &coeffs,
    InequalityType type,
    const IntervalUnion &solution,
    std::mt19937 &rng,
    int iter) {
    for (int sample = 0; sample < 100; ++sample) {
        double point;
        if (!sample_outside(solution, rng, point))
            continue;
        double value = eval_poly(poly, point);
        if (!satisfies(value, type))
            continue;
        std::ostringstream message;
        message << "FAIL: iter=" << iter
                << " point=" << point << " poly_value=" << value
                << " type=" << static_cast<int>(type)
                << " solution=" << solution.to_string()
                << " poly=" << format_coefficients(coeffs);
        ADD_FAILURE() << message.str();
        return false;
    }
    return true;
}

TEST(InequalityCompleteness, SolutionCompleteness) {
    std::mt19937 rng(123);
    const int num_iterations = 100;
    std::uniform_int_distribution<int> type_dist(0, 3);
    for (int iter = 0; iter < num_iterations; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        auto coeffs = generate_coefficients(rng);
        auto poly = build_poly(coeffs);
        InequalityType type = inequality_types[type_dist(rng)];
        auto solution = InequalitySolver::solve_inequality(poly, type, "x");
        if (solution.is_entire_line()) {
            continue;
        }
        EXPECT_TRUE(check_completeness_samples(poly, coeffs, type, solution, rng, iter));
    }
}
