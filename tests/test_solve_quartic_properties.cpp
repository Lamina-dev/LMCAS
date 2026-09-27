#include "test_solve_quartic_support.hpp"

static void verify_general_quartic_trials(std::mt19937 &rng_quartic,
                                          std::uniform_int_distribution<int> &root_dist,
                                          int NUM_GENERAL_TRIALS, double RESIDUAL_TOL) {
    for (int trial = 0; trial < NUM_GENERAL_TRIALS; ++trial) {
        int r1 = root_dist(rng_quartic);
        int r2 = root_dist(rng_quartic);
        int r3 = root_dist(rng_quartic);
        int r4 = root_dist(rng_quartic);

        int b_val = -(r1 + r2 + r3 + r4);
        int c_val = r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4;
        int d_val = -(r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4);
        int e_val = r1 * r2 * r3 * r4;

        auto roots = LMCAS::solve_quartic(
            SymbolicExpr::number(1),
            SymbolicExpr::number(b_val),
            SymbolicExpr::number(c_val),
            SymbolicExpr::number(d_val),
            SymbolicExpr::number(e_val),
            "x");

        EXPECT_EQ(roots.size(), 4u)
            << "General Trial " << trial
            << " (roots=" << r1 << "," << r2 << "," << r3 << "," << r4 << ")";
        if (roots.size() != 4) {
            continue;
        }

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            if (std::isnan(val) || std::isinf(val)) {
                continue;
            }
            double residual = eval_quartic(1.0, (double)b_val, (double)c_val,
                                           (double)d_val, (double)e_val, val);
            EXPECT_FALSE(std::abs(residual) >= RESIDUAL_TOL)
                << "General Trial " << trial << " root " << i
                << " (roots=" << r1 << "," << r2 << "," << r3 << "," << r4
                << "): |f(r)| = " << std::abs(residual)
                << ", tolerance=" << RESIDUAL_TOL << " (r = " << val << ")";
        }
    }
}

static void verify_biquadratic_trials(std::mt19937 &rng_quartic,
                                      std::uniform_int_distribution<int> &root_dist,
                                      int NUM_BIQUADRATIC_TRIALS, double RESIDUAL_TOL) {
    for (int trial = 0; trial < NUM_BIQUADRATIC_TRIALS; ++trial) {
        int r1 = root_dist(rng_quartic);
        int r2 = root_dist(rng_quartic);

        while (r1 == 0 && r2 == 0) {
            r1 = root_dist(rng_quartic);
            r2 = root_dist(rng_quartic);
        }

        int a_coeff = 1;
        int b_coeff = 0;
        int c_coeff = -(r1 * r1 + r2 * r2);
        int d_coeff = 0;
        int e_coeff = r1 * r1 * r2 * r2;

        auto roots = LMCAS::solve_quartic(
            SymbolicExpr::number(a_coeff),
            SymbolicExpr::number(b_coeff),
            SymbolicExpr::number(c_coeff),
            SymbolicExpr::number(d_coeff),
            SymbolicExpr::number(e_coeff),
            "x");

        EXPECT_EQ(roots.size(), 4u)
            << "Biquadratic Trial " << trial
            << " (r1=" << r1 << ", r2=" << r2 << ")";
        if (roots.size() != 4) {
            continue;
        }

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            if (std::isnan(val) || std::isinf(val)) {
                continue;
            }
            double residual = eval_quartic((double)a_coeff, (double)b_coeff,
                                           (double)c_coeff, (double)d_coeff,
                                           (double)e_coeff, val);
            EXPECT_FALSE(std::abs(residual) >= RESIDUAL_TOL)
                << "Biquadratic Trial " << trial << " root " << i
                << " (r1=" << r1 << ", r2=" << r2
                << "): |f(r)| = " << std::abs(residual)
                << ", tolerance=" << RESIDUAL_TOL << " (r = " << val << ")";
        }
    }
}

static void verify_depressed_even_trials(std::mt19937 &rng_quartic,
                                         std::uniform_int_distribution<int> &root_dist,
                                         int NUM_Q0_TRIALS, double RESIDUAL_TOL) {
    for (int trial = 0; trial < NUM_Q0_TRIALS; ++trial) {
        int center = root_dist(rng_quartic);
        int k = root_dist(rng_quartic);
        while (k == 0)
            k = root_dist(rng_quartic);

        int r1 = center, r2 = center, r3 = center + k, r4 = center - k;

        int b_val = -(r1 + r2 + r3 + r4);
        int c_val = r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4;
        int d_val = -(r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4);
        int e_val = r1 * r2 * r3 * r4;

        auto roots = LMCAS::solve_quartic(
            SymbolicExpr::number(1),
            SymbolicExpr::number(b_val),
            SymbolicExpr::number(c_val),
            SymbolicExpr::number(d_val),
            SymbolicExpr::number(e_val),
            "x");

        EXPECT_EQ(roots.size(), 4u)
            << "q=0 Trial " << trial
            << " (center=" << center << ", k=" << k << ")";
        if (roots.size() != 4) {
            continue;
        }

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            if (std::isnan(val) || std::isinf(val)) {
                continue;
            }
            double residual = eval_quartic(1.0, (double)b_val, (double)c_val,
                                           (double)d_val, (double)e_val, val);
            EXPECT_FALSE(std::abs(residual) >= RESIDUAL_TOL)
                << "q=0 Trial " << trial << " root " << i
                << " (center=" << center << ", k=" << k
                << "): |f(r)| = " << std::abs(residual)
                << ", tolerance=" << RESIDUAL_TOL << " (r = " << val << ")";
        }
    }
}

TEST(SolveQuarticProperties, QuarticRootVerificationRandomFromKnownRoots) {
    const double RESIDUAL_TOL = 1e-10;
    const int NUM_GENERAL_TRIALS = 40;
    const int NUM_BIQUADRATIC_TRIALS = 10;
    const int NUM_Q0_TRIALS = 10;

    std::mt19937 rng_quartic(314159);
    std::uniform_int_distribution<int> root_dist(-5, 5);

    verify_general_quartic_trials(rng_quartic, root_dist,
                                  NUM_GENERAL_TRIALS, RESIDUAL_TOL);

    verify_biquadratic_trials(rng_quartic, root_dist,
                              NUM_BIQUADRATIC_TRIALS, RESIDUAL_TOL);

    verify_depressed_even_trials(rng_quartic, root_dist,
                                 NUM_Q0_TRIALS, RESIDUAL_TOL);
}
