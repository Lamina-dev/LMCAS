#include "test_newton_raphson_support.hpp"

static bool newton_interval_residuals_are_bounded(
    const LMCAS::Polynomial<Rational> &poly,
    const std::shared_ptr<SymbolicExpr> &expr,
    const std::shared_ptr<SymbolicExpr> &df_expr,
    const std::vector<std::pair<Rational, Rational>> &intervals,
    const std::vector<int> &known_roots, int trial, int deg,
    lmmc_real_t TOLERANCE) {
    bool trial_ok = true;
    for (const auto &[lo_rat, hi_rat] : intervals) {
        lmmc_real_t lo = lo_rat.to_double();
        lmmc_real_t hi = hi_rat.to_double();
        lmmc_real_t x0 = (lo + hi) * 0.5;

        LMCAS::SolveOptions opts;
        opts.tolerance = TOLERANCE;
        opts.max_newton_iterations = 100;

        auto result = LMCAS::newton_raphson_checked(expr, df_expr, "x", x0, lo, hi, opts).value();

        if (result.has_value()) {
            lmmc_real_t residual = std::abs(
                eval_poly_at_double(poly, result->value));

            if (residual >= TOLERANCE * 100) {
                trial_ok = false;
                std::ostringstream msg;
                msg << "Trial " << trial << ": root=" << result->value
                    << " residual=" << residual << " >= " << (TOLERANCE * 100)
                    << " (degree " << deg << ", roots: [";
                for (size_t k = 0; k < known_roots.size(); ++k) {
                    if (k > 0)
                        msg << ",";
                    msg << known_roots[k];
                }
                msg << "])";
                ADD_FAILURE() << msg.str();
            }
        }
    }

    return trial_ok;
}

TEST(NewtonRaphsonResidual, NewtonRaphsonResidualBound) {
    const int NUM_TRIALS = 35;
    const lmmc_real_t TOLERANCE = 1e-12;

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
        LMCAS::Polynomial<Rational> poly = poly_from_roots(known_roots);

        LMCAS::Polynomial<Rational> dpoly = poly.differentiate();

        auto expr = LMCAS::poly_to_symbolic(poly);
        auto df_expr = LMCAS::poly_to_symbolic(dpoly);

        auto intervals = LMCAS::isolate_real_roots_checked(poly).value();

        bool trial_ok = newton_interval_residuals_are_bounded(
            poly, expr, df_expr, intervals, known_roots, trial, deg,
            TOLERANCE);

        EXPECT_TRUE(trial_ok) << "Trial " << trial << ", degree=" << deg
                              << ", tolerance=" << TOLERANCE;
    }
}
