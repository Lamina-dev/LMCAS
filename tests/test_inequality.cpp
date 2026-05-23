#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include <cmath>
#include <random>
#include <set>
#include <sstream>
#include <vector>

using namespace lamina;

static std::shared_ptr<SymbolicExpr> linear(int a, int b) {
    auto x = SymbolicExpr::variable("x");
    auto ax = SymbolicExpr::multiply(SymbolicExpr::number(a), x);
    return SymbolicExpr::add(ax, SymbolicExpr::number(b));
}

static std::shared_ptr<SymbolicExpr> negate(std::shared_ptr<SymbolicExpr> e) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), e);
}

int main() {

    auto x = SymbolicExpr::variable("x");

    TEST_CASE("Linear inequality: 2x - 3 > 0");
    {

        auto expr = linear(2, -3);
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "2x - 3 > 0 should not be empty");
        EXPECT_TRUE(!result.contains(0.0), "2x - 3 > 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "2x - 3 > 0: x=1 not in solution");
        EXPECT_TRUE(!result.contains(1.5), "2x - 3 > 0: x=1.5 (root) not in solution (strict)");
        EXPECT_TRUE(result.contains(2.0), "2x - 3 > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(100.0), "2x - 3 > 0: x=100 in solution");
        EXPECT_TRUE(!result.contains(-5.0), "2x - 3 > 0: x=-5 not in solution");
    }

    TEST_CASE("Quadratic inequality: x^2 - 4 >= 0");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(!result.is_empty(), "x^2 - 4 >= 0 should not be empty");

        EXPECT_TRUE(result.contains(-10.0), "x^2 - 4 >= 0: x=-10 in solution");
        EXPECT_TRUE(result.contains(-3.0), "x^2 - 4 >= 0: x=-3 in solution");
        EXPECT_TRUE(result.contains(3.0), "x^2 - 4 >= 0: x=3 in solution");
        EXPECT_TRUE(result.contains(10.0), "x^2 - 4 >= 0: x=10 in solution");

        EXPECT_TRUE(!result.contains(0.0), "x^2 - 4 >= 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "x^2 - 4 >= 0: x=1 not in solution");
        EXPECT_TRUE(!result.contains(-1.0), "x^2 - 4 >= 0: x=-1 not in solution");
        EXPECT_TRUE(!result.contains(1.9), "x^2 - 4 >= 0: x=1.9 not in solution");
        EXPECT_TRUE(!result.contains(-1.9), "x^2 - 4 >= 0: x=-1.9 not in solution");

        EXPECT_TRUE(result.contains(2.0), "x^2 - 4 >= 0: x=2 in solution (non-strict, root)");
        EXPECT_TRUE(result.contains(-2.0), "x^2 - 4 >= 0: x=-2 in solution (non-strict, root)");
    }

    TEST_CASE("Cubic inequality: x^3 - x < 0");
    {

        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto expr = SymbolicExpr::add(x3, negate(x));
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::LessThan, "x");

        EXPECT_TRUE(!result.is_empty(), "x^3 - x < 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "x^3 - x < 0: x=-5 in solution");
        EXPECT_TRUE(result.contains(-2.0), "x^3 - x < 0: x=-2 in solution");
        EXPECT_TRUE(!result.contains(-1.0), "x^3 - x < 0: x=-1 not in solution (strict, root)");
        EXPECT_TRUE(!result.contains(1.0), "x^3 - x < 0: x=1 not in solution (strict, root)");
        EXPECT_TRUE(result.contains(0.5), "x^3 - x < 0: x=0.5 in solution");
        EXPECT_TRUE(!result.contains(2.0), "x^3 - x < 0: x=2 not in solution");
        EXPECT_TRUE(!result.contains(-0.5), "x^3 - x < 0: x=-0.5 not in solution (between -1 and 0)");

        EXPECT_TRUE(result.contains(0.001), "x^3 - x < 0: x=0.001 in solution (just above 0)");
    }

    TEST_CASE("Repeated root: (x-1)^2*(x+2) > 0");
    {

        auto x_minus_1 = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto x_minus_1_sq = SymbolicExpr::power(x_minus_1, SymbolicExpr::number(2));
        auto x_plus_2 = SymbolicExpr::add(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::multiply(x_minus_1_sq, x_plus_2)->expand();

        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)^2*(x+2) > 0 should not be empty");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)^2*(x+2) > 0: x=-2 not in solution (root, strict)");
        EXPECT_TRUE(!result.contains(-3.0), "(x-1)^2*(x+2) > 0: x=-3 not in solution");
        EXPECT_TRUE(result.contains(0.0), "(x-1)^2*(x+2) > 0: x=0 in solution");
        EXPECT_TRUE(!result.contains(1.0), "(x-1)^2*(x+2) > 0: x=1 not in solution (root, strict)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)^2*(x+2) > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(10.0), "(x-1)^2*(x+2) > 0: x=10 in solution");
        EXPECT_TRUE(result.contains(-1.0), "(x-1)^2*(x+2) > 0: x=-1 in solution");
    }

    TEST_CASE("Rational inequality: (x-1)/(x+2) > 0");
    {
        auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

        auto result = InequalitySolver::solve_rational_inequality(
            numerator, denominator, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)/(x+2) > 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "(x-1)/(x+2) > 0: x=-5 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)/(x+2) > 0: x=-2 not in solution (den root)");
        EXPECT_TRUE(!result.contains(0.0), "(x-1)/(x+2) > 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "(x-1)/(x+2) > 0: x=1 not in solution (strict, num root)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)/(x+2) > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(100.0), "(x-1)/(x+2) > 0: x=100 in solution");
    }

    TEST_CASE("Rational inequality non-strict: (x-1)/(x+2) >= 0");
    {
        auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

        auto result = InequalitySolver::solve_rational_inequality(
            numerator, denominator, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)/(x+2) >= 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "(x-1)/(x+2) >= 0: x=-5 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)/(x+2) >= 0: x=-2 not in solution (den root excluded)");
        EXPECT_TRUE(!result.contains(0.0), "(x-1)/(x+2) >= 0: x=0 not in solution");
        EXPECT_TRUE(result.contains(1.0), "(x-1)/(x+2) >= 0: x=1 in solution (num root, non-strict)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)/(x+2) >= 0: x=2 in solution");
    }

    TEST_CASE("System of inequalities: {x^2 - 4 > 0, x < 5}");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr1 = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

        auto expr2 = SymbolicExpr::add(x, SymbolicExpr::number(-5));

        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> system = {
            {expr1, InequalityType::GreaterThan},
            {expr2, InequalityType::LessThan}
        };

        auto result = InequalitySolver::solve_inequalities(system, "x");

        EXPECT_TRUE(!result.is_empty(), "{x^2-4>0, x<5} should not be empty");
        EXPECT_TRUE(result.contains(-10.0), "{x^2-4>0, x<5}: x=-10 in solution");
        EXPECT_TRUE(result.contains(-3.0), "{x^2-4>0, x<5}: x=-3 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "{x^2-4>0, x<5}: x=-2 not in solution (strict root)");
        EXPECT_TRUE(!result.contains(0.0), "{x^2-4>0, x<5}: x=0 not in solution");
        EXPECT_TRUE(!result.contains(2.0), "{x^2-4>0, x<5}: x=2 not in solution (strict root)");
        EXPECT_TRUE(result.contains(3.0), "{x^2-4>0, x<5}: x=3 in solution");
        EXPECT_TRUE(result.contains(4.0), "{x^2-4>0, x<5}: x=4 in solution");
        EXPECT_TRUE(!result.contains(5.0), "{x^2-4>0, x<5}: x=5 not in solution (strict)");
        EXPECT_TRUE(!result.contains(6.0), "{x^2-4>0, x<5}: x=6 not in solution");
    }

    TEST_CASE("Zero polynomial: 0 > 0 -> empty");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(result.is_empty(), "0 > 0 should be empty");
        EXPECT_TRUE(!result.contains(0.0), "0 > 0: contains nothing");
        EXPECT_TRUE(!result.contains(1.0), "0 > 0: contains nothing");
    }

    TEST_CASE("Zero polynomial: 0 >= 0 -> entire line");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(result.is_entire_line(), "0 >= 0 should be entire line");
        EXPECT_TRUE(result.contains(0.0), "0 >= 0: contains 0");
        EXPECT_TRUE(result.contains(-1000.0), "0 >= 0: contains -1000");
        EXPECT_TRUE(result.contains(1000.0), "0 >= 0: contains 1000");
    }

    TEST_CASE("Zero polynomial: 0 < 0 -> empty");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessThan, "x");

        EXPECT_TRUE(result.is_empty(), "0 < 0 should be empty");
    }

    TEST_CASE("Zero polynomial: 0 <= 0 -> entire line");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessEqual, "x");

        EXPECT_TRUE(result.is_entire_line(), "0 <= 0 should be entire line");
    }

    TEST_CASE("Non-polynomial: sin(x) > 0 -> empty (cannot solve)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto result = InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(result.is_empty(), "sin(x) > 0 should return empty (non-polynomial)");
    }

    TEST_CASE("Non-polynomial: sin(x) >= 0 -> empty (cannot solve)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto result = InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(result.is_empty(), "sin(x) >= 0 should return empty (non-polynomial)");
    }

    TEST_CASE("Property 2: Solution Soundness - points in solution set satisfy inequality");
    {
        std::mt19937 rng(42);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> degree_dist(2, 2);
        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto build_poly = [&](const std::vector<int>& coeffs) -> std::shared_ptr<SymbolicExpr> {
            auto xv = SymbolicExpr::variable("x");
            std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(coeffs[0]);
            for (size_t i = 1; i < coeffs.size(); ++i) {
                if (coeffs[i] == 0) continue;
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::number(coeffs[i]),
                    SymbolicExpr::power(xv, SymbolicExpr::number(static_cast<int>(i)))
                );
                result = SymbolicExpr::add(result, term);
            }
            return result;
        };

        auto eval_poly = [](const std::shared_ptr<SymbolicExpr>& poly, double point) -> double {
            auto val_expr = SymbolicExpr::number(point);
            auto substituted = poly->substitute("x", val_expr);
            return substituted->to_numeric();
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        auto sample_inside = [&](const IntervalUnion& iu, double& out) -> bool {
            const auto& intervals = iu.intervals();
            if (intervals.empty()) return false;

            std::uniform_int_distribution<size_t> idx_dist(0, intervals.size() - 1);
            size_t idx = idx_dist(rng);
            const auto& iv = intervals[idx];

            double lo = -1000.0;
            double hi = 1000.0;

            if (!iv.lower.is_neg_infinity && iv.lower.value) {
                lo = iv.lower.value->to_numeric();
            }
            if (!iv.upper.is_pos_infinity && iv.upper.value) {
                hi = iv.upper.value->to_numeric();
            }

            double epsilon = 1e-6;
            if (iv.lower.is_open && !iv.lower.is_neg_infinity) lo += epsilon;
            if (iv.upper.is_open && !iv.upper.is_pos_infinity) hi -= epsilon;

            if (lo >= hi) {
                out = (lo + hi) / 2.0;
                return true;
            }

            std::uniform_real_distribution<double> point_dist(lo, hi);
            out = point_dist(rng);
            return true;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int degree = degree_dist(rng);
            std::vector<int> coeffs(degree + 1);
            for (int i = 0; i <= degree; ++i) {
                coeffs[i] = coeff_dist(rng);
            }

            while (coeffs[degree] == 0) {
                coeffs[degree] = coeff_dist(rng);
            }

            auto poly = build_poly(coeffs);

            InequalityType type = types[type_dist(rng)];

            auto solution = InequalitySolver::solve_inequality(poly, type, "x");

            if (solution.is_empty()) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double point;
                if (!sample_inside(solution, point)) continue;

                double value = eval_poly(poly, point);

                if (!satisfies(value, type)) {
                    std::ostringstream oss;
                    oss << "Property 2 FAIL: iter=" << iter
                        << " point=" << point << " poly_value=" << value
                        << " type=" << static_cast<int>(type)
                        << " poly=[";
                    for (size_t i = 0; i < coeffs.size(); ++i) {
                        if (i > 0) oss << ",";
                        oss << coeffs[i];
                    }
                    oss << "]";
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "Property 2: " << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed solution soundness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Property 3: Solution Completeness - points outside solution set violate inequality");
    {
        std::mt19937 rng(123);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> degree_dist(2, 2);
        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto build_poly = [&](const std::vector<int>& coeffs) -> std::shared_ptr<SymbolicExpr> {
            auto xv = SymbolicExpr::variable("x");
            std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(coeffs[0]);
            for (size_t i = 1; i < coeffs.size(); ++i) {
                if (coeffs[i] == 0) continue;
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::number(coeffs[i]),
                    SymbolicExpr::power(xv, SymbolicExpr::number(static_cast<int>(i)))
                );
                result = SymbolicExpr::add(result, term);
            }
            return result;
        };

        auto eval_poly = [](const std::shared_ptr<SymbolicExpr>& poly, double point) -> double {
            auto val_expr = SymbolicExpr::number(point);
            auto substituted = poly->substitute("x", val_expr);
            return substituted->to_numeric();
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        auto sample_outside = [&](const IntervalUnion& iu, double& out) -> bool {
            std::uniform_real_distribution<double> range_dist(-100.0, 100.0);

            for (int attempt = 0; attempt < 1000; ++attempt) {
                double candidate = range_dist(rng);

                if (iu.contains(candidate)) continue;

                bool near_boundary = false;
                for (const auto& iv : iu.intervals()) {
                    if (!iv.lower.is_neg_infinity && iv.lower.value) {
                        double boundary = iv.lower.value->to_numeric();
                        if (std::abs(candidate - boundary) < 1e-8) {
                            near_boundary = true;
                            break;
                        }
                    }
                    if (!iv.upper.is_pos_infinity && iv.upper.value) {
                        double boundary = iv.upper.value->to_numeric();
                        if (std::abs(candidate - boundary) < 1e-8) {
                            near_boundary = true;
                            break;
                        }
                    }
                }

                if (!near_boundary) {
                    out = candidate;
                    return true;
                }
            }
            return false;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int degree = degree_dist(rng);
            std::vector<int> coeffs(degree + 1);
            for (int i = 0; i <= degree; ++i) {
                coeffs[i] = coeff_dist(rng);
            }

            while (coeffs[degree] == 0) {
                coeffs[degree] = coeff_dist(rng);
            }

            auto poly = build_poly(coeffs);

            InequalityType type = types[type_dist(rng)];

            auto solution = InequalitySolver::solve_inequality(poly, type, "x");

            if (solution.is_entire_line()) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double point;
                if (!sample_outside(solution, point)) continue;

                double value = eval_poly(poly, point);

                if (satisfies(value, type)) {
                    std::ostringstream oss;
                    oss << "Property 3 FAIL: iter=" << iter
                        << " point=" << point << " poly_value=" << value
                        << " type=" << static_cast<int>(type)
                        << " solution=" << solution.to_string()
                        << " poly=[";
                    for (size_t i = 0; i < coeffs.size(); ++i) {
                        if (i > 0) oss << ",";
                        oss << coeffs[i];
                    }
                    oss << "]";
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "Property 3: " << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed solution completeness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Property 7: Endpoint Correctness (Strict vs Non-strict)");
    {
        std::mt19937 rng(777);
        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> root_count_dist(1, 2);
        std::uniform_int_distribution<int> root_val_dist(-10, 10);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int num_roots = root_count_dist(rng);
            std::set<int> root_set;
            while ((int)root_set.size() < num_roots) {
                root_set.insert(root_val_dist(rng));
            }
            std::vector<int> roots(root_set.begin(), root_set.end());

            auto x_var = SymbolicExpr::variable("x");
            auto factor0 = SymbolicExpr::add(x_var, SymbolicExpr::number(-roots[0]));
            std::shared_ptr<SymbolicExpr> poly = factor0;
            for (size_t i = 1; i < roots.size(); ++i) {
                auto factor = SymbolicExpr::add(x_var, SymbolicExpr::number(-roots[i]));
                poly = SymbolicExpr::multiply(poly, factor);
            }
            poly = poly->expand();

            InequalityType ineq_type = types[type_dist(rng)];
            bool is_strict = (ineq_type == InequalityType::GreaterThan ||
                              ineq_type == InequalityType::LessThan);

            auto solution = InequalitySolver::solve_inequality(poly, ineq_type, "x");

            bool iter_passed = true;
            for (int r : roots) {
                bool root_in_solution = solution.contains((double)r);

                if (is_strict) {

                    if (root_in_solution) {
                        std::ostringstream oss;
                        oss << "Property 7 FAIL (strict): iter=" << iter
                            << " root=" << r << " is in solution but shouldn't be"
                            << " poly=" << poly->to_string()
                            << " type=" << (int)ineq_type;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                } else {

                    if (!root_in_solution) {
                        std::ostringstream oss;
                        oss << "Property 7 FAIL (non-strict): iter=" << iter
                            << " root=" << r << " is NOT in solution but should be"
                            << " poly=" << poly->to_string()
                            << " type=" << (int)ineq_type;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "Property 7: " << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed endpoint correctness (strict vs non-strict)";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Property 8: Multiplicity Sign Change Correctness");
    {
        std::mt19937 rng(888);
        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> root_val_dist(-5, 5);
        std::uniform_int_distribution<int> even_mult_dist(0, 1);
        std::uniform_int_distribution<int> odd_mult_dist(0, 1);
        std::uniform_int_distribution<int> parity_dist(0, 1);

        const double epsilon = 0.5;

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int root = root_val_dist(rng);
            bool use_even = (parity_dist(rng) == 0);
            int multiplicity;
            if (use_even) {
                multiplicity = (even_mult_dist(rng) == 0) ? 2 : 4;
            } else {
                multiplicity = (odd_mult_dist(rng) == 0) ? 1 : 3;
            }

            auto x_var = SymbolicExpr::variable("x");
            auto factor = SymbolicExpr::add(x_var, SymbolicExpr::number(-root));
            auto poly = factor;
            for (int i = 1; i < multiplicity; ++i) {
                poly = SymbolicExpr::multiply(poly, factor);
            }
            poly = poly->expand();

            double left_val, right_val;
            {
                auto sub_left = poly->substitute("x", SymbolicExpr::number((double)root - epsilon));
                auto simp_left = sub_left->simplify();
                left_val = simp_left->to_numeric();

                auto sub_right = poly->substitute("x", SymbolicExpr::number((double)root + epsilon));
                auto simp_right = sub_right->simplify();
                right_val = simp_right->to_numeric();
            }

            const double sign_tol = 1e-10;
            int left_sign = (left_val > sign_tol) ? 1 : ((left_val < -sign_tol) ? -1 : 0);
            int right_sign = (right_val > sign_tol) ? 1 : ((right_val < -sign_tol) ? -1 : 0);

            if (left_sign == 0 || right_sign == 0) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            if (use_even) {

                if (left_sign != right_sign) {
                    std::ostringstream oss;
                    oss << "Property 8 FAIL (even mult): iter=" << iter
                        << " root=" << root << " mult=" << multiplicity
                        << " left_sign=" << left_sign << " right_sign=" << right_sign
                        << " left_val=" << left_val << " right_val=" << right_val;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                }
            } else {

                if (left_sign == right_sign) {
                    std::ostringstream oss;
                    oss << "Property 8 FAIL (odd mult): iter=" << iter
                        << " root=" << root << " mult=" << multiplicity
                        << " left_sign=" << left_sign << " right_sign=" << right_sign
                        << " left_val=" << left_val << " right_val=" << right_val;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "Property 8: " << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed multiplicity sign change correctness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Property 11: Parametric Inequality Consistency");
    {
        std::mt19937 rng(1111);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 50;
        int pass_count = 0;

        std::uniform_int_distribution<int> leading_dist(1, 4);
        std::uniform_int_distribution<int> sign_dist(0, 1);
        std::uniform_int_distribution<int> param_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int a_val = leading_dist(rng);
            if (sign_dist(rng)) a_val = -a_val;

            int p_val = param_dist(rng);
            int q_val = param_dist(rng);

            int disc = p_val * p_val - 4 * a_val * q_val;
            if (disc < 0) {

                if (a_val > 0) {
                    q_val = (p_val * p_val) / (4 * a_val) - 1;
                } else {
                    q_val = (p_val * p_val) / (4 * a_val) + 1;
                }
                disc = p_val * p_val - 4 * a_val * q_val;
                if (disc < 0) {
                    ++pass_count;
                    continue;
                }
            }

            InequalityType ineq_type = types[type_dist(rng)];

            auto x_var = SymbolicExpr::variable("x");
            auto p_param = SymbolicExpr::variable("p");
            auto q_param = SymbolicExpr::variable("q");

            auto ax2 = SymbolicExpr::multiply(
                SymbolicExpr::number(a_val),
                SymbolicExpr::power(x_var, SymbolicExpr::number(2)));
            auto px = SymbolicExpr::multiply(p_param, x_var);
            auto parametric_expr = SymbolicExpr::add(SymbolicExpr::add(ax2, px), q_param);

            auto parametric_result = InequalitySolver::solve_parametric_inequality(
                parametric_expr, ineq_type, "x", {"p", "q"});

            auto concrete_expr = SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::number(a_val),
                        SymbolicExpr::power(x_var, SymbolicExpr::number(2))),
                    SymbolicExpr::multiply(SymbolicExpr::number(p_val), x_var)),
                SymbolicExpr::number(q_val));

            auto direct_solution = InequalitySolver::solve_inequality(
                concrete_expr, ineq_type, "x");

            IntervalUnion parametric_solution = IntervalUnion::empty();
            bool found_case = false;

            if (parametric_result.is_empty()) {
                ++pass_count;
                continue;
            }

            if (parametric_result.is_single()) {
                parametric_solution = parametric_result.single_solution();
                found_case = true;
            } else {

                for (const auto& pcase : parametric_result.cases) {
                    if (!pcase.condition) {
                        parametric_solution = pcase.solution;
                        found_case = true;
                        break;
                    }
                }
                if (!found_case && !parametric_result.cases.empty()) {
                    parametric_solution = parametric_result.cases[0].solution;
                    found_case = true;
                }
            }

            if (!found_case) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            std::uniform_real_distribution<double> sample_dist(-20.0, 20.0);

            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double test_point = sample_dist(rng);

                double poly_val = (double)a_val * test_point * test_point
                                + (double)p_val * test_point
                                + (double)q_val;

                if (std::abs(poly_val) < 1e-6) continue;

                bool expected_result = satisfies(poly_val, ineq_type);

                bool in_direct = direct_solution.contains(test_point);

                bool in_parametric = false;
                const auto& param_intervals = parametric_solution.intervals();

                if (parametric_solution.is_empty()) {
                    in_parametric = false;
                } else if (parametric_solution.is_entire_line()) {
                    in_parametric = true;
                } else {
                    for (const auto& iv : param_intervals) {
                        double lo = -1e18;
                        double hi = 1e18;
                        bool lo_open = true;
                        bool hi_open = true;
                        bool lo_valid = true;
                        bool hi_valid = true;

                        if (!iv.lower.is_neg_infinity && iv.lower.value) {
                            auto lo_expr = iv.lower.value->substitute("p", SymbolicExpr::number(p_val));
                            lo_expr = lo_expr->substitute("q", SymbolicExpr::number(q_val));
                            lo_expr = lo_expr->simplify();
                            auto lo_opt = test_numeric_eval(lo_expr);
                            if (lo_opt && std::isfinite(*lo_opt)) {
                                lo = *lo_opt;
                                lo_open = iv.lower.is_open;
                            } else {
                                // Cannot evaluate lower bound; skip this interval.
                                lo_valid = false;
                            }
                        }

                        if (!iv.upper.is_pos_infinity && iv.upper.value) {
                            auto hi_expr = iv.upper.value->substitute("p", SymbolicExpr::number(p_val));
                            hi_expr = hi_expr->substitute("q", SymbolicExpr::number(q_val));
                            hi_expr = hi_expr->simplify();
                            auto hi_opt = test_numeric_eval(hi_expr);
                            if (hi_opt && std::isfinite(*hi_opt)) {
                                hi = *hi_opt;
                                hi_open = iv.upper.is_open;
                            } else {
                                hi_valid = false;
                            }
                        }

                        // If we couldn't evaluate an endpoint, try the
                        // non-simplified version as well (simplify may have
                        // introduced a form that numeric_eval can't handle).
                        if (!lo_valid && !iv.lower.is_neg_infinity && iv.lower.value) {
                            auto lo_expr = iv.lower.value->substitute("p", SymbolicExpr::number(p_val));
                            lo_expr = lo_expr->substitute("q", SymbolicExpr::number(q_val));
                            auto lo_opt = test_numeric_eval(lo_expr);
                            if (lo_opt && std::isfinite(*lo_opt)) {
                                lo = *lo_opt;
                                lo_open = iv.lower.is_open;
                                lo_valid = true;
                            }
                        }
                        if (!hi_valid && !iv.upper.is_pos_infinity && iv.upper.value) {
                            auto hi_expr = iv.upper.value->substitute("p", SymbolicExpr::number(p_val));
                            hi_expr = hi_expr->substitute("q", SymbolicExpr::number(q_val));
                            auto hi_opt = test_numeric_eval(hi_expr);
                            if (hi_opt && std::isfinite(*hi_opt)) {
                                hi = *hi_opt;
                                hi_open = iv.upper.is_open;
                                hi_valid = true;
                            }
                        }

                        if (!lo_valid || !hi_valid) continue;

                        bool in_interval = false;
                        if (lo_open) {
                            in_interval = (test_point > lo + 1e-10);
                        } else {
                            in_interval = (test_point >= lo - 1e-10);
                        }
                        if (in_interval) {
                            if (hi_open) {
                                in_interval = (test_point < hi - 1e-10);
                            } else {
                                in_interval = (test_point <= hi + 1e-10);
                            }
                        }

                        if (in_interval) {
                            in_parametric = true;
                            break;
                        }
                    }
                }

                if (in_parametric != expected_result) {

                    // If we couldn't evaluate any interval endpoint (all were
                    // skipped due to lo_valid/hi_valid failures), we cannot
                    // reliably compare. Skip this sample point.
                    if (!in_parametric && !parametric_solution.is_empty() &&
                        !parametric_solution.is_entire_line()) {
                        // Check if we actually evaluated at least one interval.
                        bool any_evaluated = false;
                        for (const auto& iv2 : param_intervals) {
                            bool l_ok = iv2.lower.is_neg_infinity;
                            bool u_ok = iv2.upper.is_pos_infinity;
                            if (!l_ok && iv2.lower.value) {
                                auto le = iv2.lower.value->substitute("p", SymbolicExpr::number(p_val));
                                le = le->substitute("q", SymbolicExpr::number(q_val));
                                auto lv = test_numeric_eval(le);
                                l_ok = lv.has_value();
                            } else { l_ok = true; }
                            if (!u_ok && iv2.upper.value) {
                                auto ue = iv2.upper.value->substitute("p", SymbolicExpr::number(p_val));
                                ue = ue->substitute("q", SymbolicExpr::number(q_val));
                                auto uv = test_numeric_eval(ue);
                                u_ok = uv.has_value();
                            } else { u_ok = true; }
                            if (l_ok && u_ok) { any_evaluated = true; break; }
                        }
                        if (!any_evaluated) continue;
                    }

                    bool near_root = false;
                    double root1 = (-p_val + std::sqrt(std::abs(disc))) / (2.0 * a_val);
                    double root2 = (-p_val - std::sqrt(std::abs(disc))) / (2.0 * a_val);
                    if (std::abs(test_point - root1) < 1e-4 ||
                        std::abs(test_point - root2) < 1e-4) {
                        near_root = true;
                    }

                    if (!near_root) {
                        std::ostringstream oss;
                        oss << "Property 11 FAIL: iter=" << iter
                            << " a=" << a_val << " p=" << p_val << " q=" << q_val
                            << " type=" << static_cast<int>(ineq_type)
                            << " point=" << test_point
                            << " expected=" << expected_result
                            << " in_parametric=" << in_parametric
                            << " in_direct=" << in_direct
                            << " poly_val=" << poly_val;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "Property 11: " << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed parametric inequality consistency";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    return TEST_REPORT();
}
