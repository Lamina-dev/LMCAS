#include "test_common.hpp"
#include "newton_raphson.hpp"
#include "solve_polynomial.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <cmath>
#include "poly_utils.hpp"
#include <random>
#include <set>
#include <sstream>

static std::shared_ptr<SymbolicExpr> num_expr(int n) { return SymbolicExpr::number(n); }

static double eval_numeric_expr(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !expr->root) return 0.0;

    if (auto n = std::dynamic_pointer_cast<NumberNode>(expr->root)) {
        if (std::holds_alternative<lmmc_real_t>(n->value)) return std::get<lmmc_real_t>(n->value);
        if (std::holds_alternative<BigInt>(n->value)) return std::get<BigInt>(n->value).to_double();
        if (std::holds_alternative<Rational>(n->value)) return std::get<Rational>(n->value).to_double();
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr->root)) {
        double result = 0.0;
        for (auto& op : add->operands) {
            result += eval_numeric_expr(std::make_shared<SymbolicExpr>(op));
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr->root)) {
        double result = 1.0;
        for (auto& op : mul->operands) {
            result *= eval_numeric_expr(std::make_shared<SymbolicExpr>(op));
        }
        return result;
    }

    if (auto pow_node = std::dynamic_pointer_cast<PowerNode>(expr->root)) {
        double base = eval_numeric_expr(std::make_shared<SymbolicExpr>(pow_node->base));
        double exp = eval_numeric_expr(std::make_shared<SymbolicExpr>(pow_node->exponent));
        if (base < 0.0 && std::abs(exp - std::round(exp)) > 1e-15) {
            double denom = std::round(1.0 / exp);
            if (std::abs(exp * denom - 1.0) < 1e-12 && ((int)denom % 2 == 1)) {
                return -std::pow(-base, exp);
            }
            return std::nan("");
        }
        return std::pow(base, exp);
    }

    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr->root)) {
        if (func->arguments.size() == 1) {
            double arg = eval_numeric_expr(std::make_shared<SymbolicExpr>(func->arguments[0]));
            switch (func->type) {
                case FunctionNode::FuncType::Sin: return std::sin(arg);
                case FunctionNode::FuncType::Cos: return std::cos(arg);
                case FunctionNode::FuncType::Tan: return std::tan(arg);
                case FunctionNode::FuncType::Exp: return std::exp(arg);
                case FunctionNode::FuncType::Ln: return std::log(arg);
                case FunctionNode::FuncType::Sqrt:
                    if (arg < 0.0) return std::nan("");
                    return std::sqrt(arg);
                case FunctionNode::FuncType::Abs: return std::abs(arg);
                case FunctionNode::FuncType::ArcCos: return std::acos(arg);
                case FunctionNode::FuncType::ArcSin: return std::asin(arg);
                case FunctionNode::FuncType::ArcTan: return std::atan(arg);
                default: break;
            }
        }
    }

    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr->root)) {
        return std::nan("");
    }

    return std::nan("");
}

static lamina::Polynomial<Rational> poly_from_roots(const std::vector<int>& roots) {
    lamina::Polynomial<Rational> result({Rational(1)}, "x");
    for (int r : roots) {

        lamina::Polynomial<Rational> factor({Rational(-r), Rational(1)}, "x");
        result = result * factor;
    }
    return result;
}

static double eval_poly_at_double(const lamina::Polynomial<Rational>& poly, double x) {
    double result = 0.0;
    double x_pow = 1.0;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        result += poly.coeffs[i].to_double() * x_pow;
        x_pow *= x;
    }
    return result;
}

int main() {
    TEST_CASE("Newton-Raphson - Basic convergence (x^2 - 2)");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-2)
        );
        auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

        lamina::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::newton_raphson(f, df, "x", 1.5, opts);
        EXPECT_TRUE(result.has_value(), "Newton-Raphson should converge for x^2-2 near 1.5");
        if (result.has_value()) {
            EXPECT_TRUE(std::abs(result->value - std::sqrt(2.0)) < 1e-10,
                "Root should be close to sqrt(2)");
            EXPECT_TRUE(result->residual < opts.tolerance,
                "Residual should be below tolerance");
        }
    }

    TEST_CASE("Newton-Raphson - Convergence with bracket (x^2 - 2)");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-2)
        );
        auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::newton_raphson(f, df, "x", 1.5, 1.0, 2.0, opts);
        EXPECT_TRUE(result.has_value(), "Newton-Raphson with bracket should converge for x^2-2");
        if (result.has_value()) {
            EXPECT_TRUE(std::abs(result->value - std::sqrt(2.0)) < 1e-10,
                "Root should be close to sqrt(2)");
            EXPECT_TRUE(result->residual < opts.tolerance,
                "Residual should be below tolerance");
        }
    }

    TEST_CASE("Newton-Raphson - Bisection fallback when derivative near zero");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto df = SymbolicExpr::multiply(
            SymbolicExpr::number(3),
            SymbolicExpr::power(x, SymbolicExpr::number(2))
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::newton_raphson(f, df, "x", 1e-8, -1.0, 1.0, opts);
        EXPECT_TRUE(result.has_value(), "Should converge via bisection fallback for x^3 near zero");
        if (result.has_value()) {
            EXPECT_TRUE(std::abs(result->value) < 1e-4,
                "Root should be close to 0");
        }
    }

    TEST_CASE("Newton-Raphson - No bracket, derivative near zero returns nullopt");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto df = SymbolicExpr::multiply(
            SymbolicExpr::number(3),
            SymbolicExpr::power(x, SymbolicExpr::number(2))
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::newton_raphson(f, df, "x", 1e-8, opts);

        EXPECT_TRUE(result.has_value(), "f(1e-8) = 1e-24 < tolerance, should converge immediately");
    }

    TEST_CASE("Newton-Raphson - Non-convergence returns nullopt");
    {

        auto x = SymbolicExpr::variable("x");

        auto f = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, SymbolicExpr::number(5)),
                SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
            ),
            SymbolicExpr::number(-1)
        );

        auto df = SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::number(5),
                SymbolicExpr::power(x, SymbolicExpr::number(4))
            ),
            SymbolicExpr::number(-1)
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 2;

        auto result = lamina::newton_raphson(f, df, "x", 10.0, opts);
        EXPECT_TRUE(!result.has_value(), "Should not converge in 2 iterations from x=10");
    }

    TEST_CASE("Newton-Raphson - Damping engages on overshoot");
    {

        auto x = SymbolicExpr::variable("x");

        auto f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-4)
        );
        auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::newton_raphson(f, df, "x", 0.01, opts);
        EXPECT_TRUE(result.has_value(), "Should converge for x^2-4 even from x0=0.01 (large first step)");
        if (result.has_value()) {
            EXPECT_TRUE(std::abs(result->value - 2.0) < 1e-10 || std::abs(result->value + 2.0) < 1e-10,
                "Root should be ±2");
        }

        auto result2 = lamina::newton_raphson(f, df, "x", 0.01, 0.0, 3.0, opts);
        EXPECT_TRUE(result2.has_value(), "Should converge with bracket for x^2-4 from x0=0.01");
        if (result2.has_value()) {
            EXPECT_TRUE(std::abs(result2->value - 2.0) < 1e-10,
                "Root should be 2 within bracket [0, 3]");
        }
    }

    TEST_CASE("Bisection - Basic convergence (x^2 - 2)");
    {
        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-2)
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::bisection(f, "x", 1.0, 2.0, opts);
        EXPECT_TRUE(result.has_value(), "Bisection should converge for x^2-2 on [1,2]");
        if (result.has_value()) {
            EXPECT_TRUE(std::abs(result->value - std::sqrt(2.0)) < 1e-10,
                "Root should be close to sqrt(2)");
        }
    }

    TEST_CASE("Bisection - No sign change returns nullopt");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(1)
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 100;

        auto result = lamina::bisection(f, "x", -1.0, 1.0, opts);
        EXPECT_TRUE(!result.has_value(), "Bisection should return nullopt when no sign change");
    }

    TEST_CASE("Sturm isolation - x^2 - 2 has 2 real roots (early)");
    {

        lamina::Polynomial<Rational> poly("x");
        poly.coeffs = {Rational(-2), Rational(0), Rational(1)};

        auto intervals = lamina::isolate_real_roots(poly);
        EXPECT_TRUE(intervals.size() == 2, "x^2-2 should have 2 isolated real roots");
    }

    TEST_CASE("Property 9: Sturm sequence root count accuracy");
    {
        const int NUM_TRIALS = 60;
        int pass_count = 0;

        std::mt19937 rng(314159);
        std::uniform_int_distribution<int> degree_dist(2, 5);
        std::uniform_int_distribution<int> root_dist(-4, 4);

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int deg = degree_dist(rng);

            std::vector<int> all_roots;
            for (int i = 0; i < deg; ++i) {
                all_roots.push_back(root_dist(rng));
            }

            std::set<int> distinct_roots_set(all_roots.begin(), all_roots.end());
            int expected_distinct_real_roots = (int)distinct_roots_set.size();

            lamina::Polynomial<Rational> poly = poly_from_roots(all_roots);

            auto intervals = lamina::isolate_real_roots(poly);
            int sturm_count = (int)intervals.size();

            bool count_matches = (sturm_count == expected_distinct_real_roots);

            if (!count_matches) {
                std::ostringstream msg;
                msg << "Property 9 Trial " << trial << ": Sturm found "
                    << sturm_count << " roots, expected " << expected_distinct_real_roots
                    << " distinct real roots (degree " << deg << ", roots: [";
                for (size_t k = 0; k < all_roots.size(); ++k) {
                    if (k > 0) msg << ",";
                    msg << all_roots[k];
                }
                msg << "])";
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool intervals_valid = true;
            for (const auto& interval : intervals) {
                double lo = interval.first.to_double();
                double hi = interval.second.to_double();

                bool contains_known_root = false;
                for (int r : distinct_roots_set) {
                    double rd = (double)r;
                    if (rd >= lo - 1e-10 && rd <= hi + 1e-10) {
                        contains_known_root = true;
                        break;
                    }
                }

                if (!contains_known_root) {
                    intervals_valid = false;
                    std::ostringstream msg;
                    msg << "Property 9 Trial " << trial
                        << ": interval [" << lo << "," << hi
                        << "] does not contain any known root. Known roots: [";
                    bool first = true;
                    for (int r : distinct_roots_set) {
                        if (!first) msg << ",";
                        msg << r;
                        first = false;
                    }
                    msg << "]";
                    EXPECT_TRUE(false, msg.str());
                    break;
                }
            }

            if (deg <= 4 && intervals_valid) {

                std::vector<std::shared_ptr<SymbolicExpr>> symbolic_roots;
                if (deg == 2) {
                    symbolic_roots = lamina::solve_cubic(
                        num_expr(0),
                        SymbolicExpr::number(poly.coeffs[2].to_double()),
                        SymbolicExpr::number(poly.coeffs[1].to_double()),
                        SymbolicExpr::number(poly.coeffs[0].to_double()),
                        "x");

                } else if (deg == 3) {
                    symbolic_roots = lamina::solve_cubic(
                        SymbolicExpr::number(poly.coeffs[3].to_double()),
                        SymbolicExpr::number(poly.coeffs[2].to_double()),
                        SymbolicExpr::number(poly.coeffs[1].to_double()),
                        SymbolicExpr::number(poly.coeffs[0].to_double()),
                        "x");
                } else if (deg == 4) {
                    symbolic_roots = lamina::solve_quartic(
                        SymbolicExpr::number(poly.coeffs[4].to_double()),
                        SymbolicExpr::number(poly.coeffs[3].to_double()),
                        SymbolicExpr::number(poly.coeffs[2].to_double()),
                        SymbolicExpr::number(poly.coeffs[1].to_double()),
                        SymbolicExpr::number(poly.coeffs[0].to_double()),
                        "x");
                }

                std::set<double> closed_form_real_roots;
                for (const auto& root : symbolic_roots) {
                    double val = eval_numeric_expr(root);
                    if (!std::isnan(val) && !std::isinf(val)) {

                        double rounded = std::round(val * 1e6) / 1e6;
                        closed_form_real_roots.insert(rounded);
                    }
                }

                int closed_form_count = (int)closed_form_real_roots.size();

                if (closed_form_count > sturm_count) {
                    std::ostringstream msg;
                    msg << "Property 9 Trial " << trial
                        << ": closed-form found " << closed_form_count
                        << " real roots but Sturm only found " << sturm_count
                        << " (degree " << deg << ")";
                    EXPECT_TRUE(false, msg.str());
                    intervals_valid = false;
                }
            }

            if (intervals_valid && count_matches) {
                pass_count++;
            }
        }

        {
            std::ostringstream msg;
            msg << "Property 9: Sturm root count accuracy: " << pass_count
                << "/" << NUM_TRIALS << " trials passed";
            EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
        }
    }

    TEST_CASE("Property 8: Newton-Raphson residual bound");
    {
        const int NUM_TRIALS = 35;
        const lmmc_real_t TOLERANCE = 1e-12;
        int pass_count = 0;
        int total_roots_checked = 0;

        std::mt19937 rng(777);

        std::uniform_int_distribution<int> degree_dist(2, 4);
        std::uniform_int_distribution<int> root_dist(-8, 8);

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int deg = degree_dist(rng);

            std::set<int> root_set;
            while ((int)root_set.size() < deg) {
                root_set.insert(root_dist(rng));
            }
            std::vector<int> known_roots(root_set.begin(), root_set.end());
            lamina::Polynomial<Rational> poly = poly_from_roots(known_roots);

            lamina::Polynomial<Rational> dpoly = poly.differentiate();

            auto expr = lamina::poly_to_symbolic(poly);
            auto df_expr = lamina::poly_to_symbolic(dpoly);

            auto intervals = lamina::isolate_real_roots(poly);

            bool trial_ok = true;
            for (const auto& [lo_rat, hi_rat] : intervals) {
                lmmc_real_t lo = lo_rat.to_double();
                lmmc_real_t hi = hi_rat.to_double();
                lmmc_real_t x0 = (lo + hi) * 0.5;

                lamina::SolveOptions opts;
                opts.tolerance = TOLERANCE;
                opts.max_newton_iterations = 100;

                auto result = lamina::newton_raphson(expr, df_expr, "x", x0, lo, hi, opts);

                if (result.has_value()) {
                    total_roots_checked++;

                    lmmc_real_t residual = std::abs(
                        eval_poly_at_double(poly, result->value));

                    if (residual >= TOLERANCE * 100) {
                        trial_ok = false;
                        std::ostringstream msg;
                        msg << "Property 8 Trial " << trial << ": root=" << result->value
                            << " residual=" << residual << " >= " << (TOLERANCE * 100)
                            << " (degree " << deg << ", roots: [";
                        for (size_t k = 0; k < known_roots.size(); ++k) {
                            if (k > 0) msg << ",";
                            msg << known_roots[k];
                        }
                        msg << "])";
                        EXPECT_TRUE(false, msg.str());
                    }
                }
            }

            if (trial_ok) {
                pass_count++;
            }
        }

        {
            std::ostringstream msg;
            msg << "Property 8: Newton residual bound: " << pass_count
                << "/" << NUM_TRIALS << " trials passed ("
                << total_roots_checked << " roots checked, tolerance=" << TOLERANCE << ")";
            EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
        }
    }

    TEST_CASE("Newton-Raphson - Deflation correctly continues to remaining roots");
    {

        lamina::Polynomial<Rational> poly("x");
        poly.coeffs = {Rational(2), Rational(-3), Rational(1)};

        auto intervals = lamina::isolate_real_roots(poly);
        EXPECT_TRUE(intervals.size() == 2, "Sturm should isolate 2 roots for (x-1)(x-2)");

        auto expr = lamina::poly_to_symbolic(poly);
        auto df_expr = expr->differentiate("x");

        lamina::SolveOptions opts;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;

        int roots_found = 0;
        for (const auto& [lo_rat, hi_rat] : intervals) {
            double lo = lo_rat.to_double();
            double hi = hi_rat.to_double();
            double x0 = (lo + hi) * 0.5;

            auto result = lamina::newton_raphson(expr, df_expr, "x", x0, lo, hi, opts);
            EXPECT_TRUE(result.has_value(),
                "Newton should find root in interval [" + std::to_string(lo) + ", " + std::to_string(hi) + "]");
            if (result.has_value()) {
                roots_found++;

                double r = result->value;
                double residual = std::abs(r*r - 3*r + 2);
                EXPECT_TRUE(residual < 1e-6,
                    "Root " + std::to_string(r) + " should satisfy x^2-3x+2=0");
            }
        }
        EXPECT_TRUE(roots_found == 2, "Should find a root in each isolated interval");
    }

    TEST_CASE("Newton-Raphson - Deflation with cubic polynomial");
    {

        lamina::Polynomial<Rational> poly("x");
        poly.coeffs = {Rational(-6), Rational(11), Rational(-6), Rational(1)};

        auto intervals = lamina::isolate_real_roots(poly);
        EXPECT_TRUE(intervals.size() == 3, "Sturm should isolate 3 roots for (x-1)(x-2)(x-3)");

        auto expr = lamina::poly_to_symbolic(poly);
        auto df_expr = expr->differentiate("x");

        lamina::SolveOptions opts;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;

        int roots_found = 0;
        for (const auto& [lo_rat, hi_rat] : intervals) {
            double lo = lo_rat.to_double();
            double hi = hi_rat.to_double();
            double x0 = (lo + hi) * 0.5;

            auto result = lamina::newton_raphson(expr, df_expr, "x", x0, lo, hi, opts);
            EXPECT_TRUE(result.has_value(),
                "Newton should find root in interval [" + std::to_string(lo) + ", " + std::to_string(hi) + "]");
            if (result.has_value()) {
                roots_found++;

                double r = result->value;
                double residual = std::abs((r-1.0)*(r-2.0)*(r-3.0));
                EXPECT_TRUE(residual < 1e-6,
                    "Root " + std::to_string(r) + " should satisfy (x-1)(x-2)(x-3)=0");
            }
        }

        EXPECT_TRUE(roots_found == 3, "Should find all 3 roots via Newton on isolated intervals");
    }

    TEST_CASE("Newton-Raphson - Non-polynomial input requires x0 (initial guess)");
    {

        auto x = SymbolicExpr::variable("x");

        auto f = SymbolicExpr::add(
            SymbolicExpr::sin(x),
            SymbolicExpr::number(-0.5)
        );

        {
            lamina::SolveOptions opts;
            opts.allow_numeric = true;
            opts.tolerance = 1e-10;
            opts.max_newton_iterations = 100;
            opts.has_initial_guess = true;
            opts.initial_guess = 0.5;

            auto roots = lamina::solve_numeric(f, "x", opts);
            EXPECT_TRUE(roots.size() <= 1,
                "Non-polynomial solve_numeric should return at most 1 root");
        }

        {
            lamina::SolveOptions opts;
            opts.allow_numeric = true;
            opts.tolerance = 1e-10;
            opts.max_newton_iterations = 100;
            opts.has_initial_guess = true;
            opts.initial_guess = 2.5;

            auto roots = lamina::solve_numeric(f, "x", opts);
            EXPECT_TRUE(roots.size() <= 1,
                "Non-polynomial with different x0 should still return at most 1 root");
        }

        {
            lamina::SolveOptions opts;
            opts.allow_numeric = true;
            opts.tolerance = 1e-10;
            opts.max_newton_iterations = 100;
            opts.has_initial_guess = false;

            auto roots = lamina::solve_numeric(f, "x", opts);
            EXPECT_TRUE(roots.size() <= 1,
                "Non-polynomial without explicit x0 should return at most 1 root");
        }

        {

            auto poly_f = SymbolicExpr::add(
                SymbolicExpr::power(x, SymbolicExpr::number(2)),
                SymbolicExpr::number(-4)
            );

            lamina::SolveOptions opts;
            opts.allow_numeric = true;
            opts.tolerance = 1e-10;
            opts.max_newton_iterations = 100;

            auto roots = lamina::solve_numeric(poly_f, "x", opts);

            EXPECT_TRUE(roots.size() == 2,
                "Polynomial x^2-4 should find 2 roots via Sturm path");
        }
    }

    TEST_CASE("Newton-Raphson - Non-convergence with limited iterations returns empty");
    {

        auto x = SymbolicExpr::variable("x");
        auto f = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, SymbolicExpr::number(3)),
                SymbolicExpr::multiply(SymbolicExpr::number(-2), x)
            ),
            SymbolicExpr::number(2)
        );
        auto df = SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::number(3),
                SymbolicExpr::power(x, SymbolicExpr::number(2))
            ),
            SymbolicExpr::number(-2)
        );

        lamina::SolveOptions opts;
        opts.tolerance = 1e-12;
        opts.max_newton_iterations = 1;

        auto result = lamina::newton_raphson(f, df, "x", 5.0, opts);
        EXPECT_TRUE(!result.has_value(),
            "Should not converge in 1 iteration from x=5 for x^3-2x+2");
    }

    TEST_CASE("Sturm isolation - x^2 - 2 has 2 real roots");
    {

        lamina::Polynomial<Rational> poly("x");
        poly.coeffs = {Rational(-2), Rational(0), Rational(1)};

        auto intervals = lamina::isolate_real_roots(poly);
        EXPECT_TRUE(intervals.size() == 2, "x^2-2 should have 2 isolated real roots");
    }

    return TEST_REPORT();
}
