#include "test_common.hpp"
#include "solve_polynomial.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <string>

static double eval_quartic(double a, double b, double c, double d, double e, double x) {
    return a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
}

int main() {
    TEST_CASE("Solve Quartic - Biquadratic");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(0);
        auto c = SymbolicExpr::number(-5);
        auto d = SymbolicExpr::number(0);
        auto e = SymbolicExpr::number(4);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Biquadratic should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(1, 0, -5, 0, 4, val);
            EXPECT_TRUE(std::abs(residual) < 1e-8,
                "Biquadratic root " + std::to_string(i) + " residual < 1e-8 (got " + std::to_string(residual) + ")");
        }
    }

    TEST_CASE("Solve Quartic - General (Ferrari)");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(-10);
        auto c = SymbolicExpr::number(35);
        auto d = SymbolicExpr::number(-50);
        auto e = SymbolicExpr::number(24);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "General quartic should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(1, -10, 35, -50, 24, val);
            EXPECT_TRUE(std::abs(residual) < 1e-8,
                "General quartic root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }
    }

    TEST_CASE("Solve Quartic - Depressed q=0 case");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(2);
        auto c = SymbolicExpr::number(-7);
        auto d = SymbolicExpr::number(-8);
        auto e = SymbolicExpr::number(12);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Quartic (x-1)(x+2)(x-2)(x+3) should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(1, 2, -7, -8, 12, val);
            EXPECT_TRUE(std::abs(residual) < 1e-8,
                "Root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }
    }

    TEST_CASE("Solve Quartic - Leading coefficient != 1");

    {
        auto a = SymbolicExpr::number(2);
        auto b = SymbolicExpr::number(-20);
        auto c = SymbolicExpr::number(70);
        auto d = SymbolicExpr::number(-100);
        auto e = SymbolicExpr::number(48);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Quartic with a=2 should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(2, -20, 70, -100, 48, val);
            EXPECT_TRUE(std::abs(residual) < 1e-8,
                "Root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }
    }

    TEST_CASE("Solve Quartic - a=0 delegates to cubic");

    {
        auto a = SymbolicExpr::number(0);
        auto b = SymbolicExpr::number(1);
        auto c = SymbolicExpr::number(-6);
        auto d = SymbolicExpr::number(11);
        auto e = SymbolicExpr::number(-6);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 3, "Quartic with a=0 should delegate to cubic and return 3 roots");
    }

    TEST_CASE("Solve Quartic - Biquadratic shortcut verifies root values ±1, ±2");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(0);
        auto c = SymbolicExpr::number(-5);
        auto d = SymbolicExpr::number(0);
        auto e = SymbolicExpr::number(4);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Biquadratic x^4-5x^2+4 should return 4 roots");

        std::vector<double> vals;
        for (auto& r : roots) vals.push_back(r->to_numeric());
        std::sort(vals.begin(), vals.end());

        double expected[] = {-2.0, -1.0, 1.0, 2.0};
        bool all_match = true;
        for (int i = 0; i < 4; ++i) {
            if (std::abs(vals[i] - expected[i]) > 1e-8) {
                all_match = false;
                break;
            }
        }
        EXPECT_TRUE(all_match,
            "Biquadratic roots should be {-2, -1, 1, 2} (got {"
            + std::to_string(vals[0]) + ", " + std::to_string(vals[1]) + ", "
            + std::to_string(vals[2]) + ", " + std::to_string(vals[3]) + "})");
    }

    TEST_CASE("Solve Biquadratic - exact huge coefficient avoids unsafe numeric conversion");
    {
        std::string huge_digits = "1" + std::string(400, '0');
        auto roots = lamina::solve_biquadratic(
            SymbolicExpr::number(BigInt(huge_digits)),
            SymbolicExpr::number(0),
            SymbolicExpr::number(-1),
            "x");

        EXPECT_TRUE(roots.size() == 4,
            "huge exact biquadratic should return symbolic square-root roots");

        int zero_roots = 0;
        for (const auto& root : roots) {
            std::string text = root ? root->to_string() : "";
            EXPECT_TRUE(text.find("inf") == std::string::npos &&
                        text.find("nan") == std::string::npos,
                "huge exact biquadratic roots should not contain fabricated inf/nan");
            if (root && root->simplify()->is_zero()) {
                ++zero_roots;
            }
        }
        EXPECT_TRUE(zero_roots < 4,
            "huge exact biquadratic must not collapse nonzero roots to zero");
    }

    TEST_CASE("Solve Quartic - q=0 after depression (non-biquadratic)");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(8);
        auto c = SymbolicExpr::number(22);
        auto d = SymbolicExpr::number(24);
        auto e = SymbolicExpr::number(9);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "q=0 quartic (x+1)^2(x+3)^2 should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(1, 8, 22, 24, 9, val);
            EXPECT_TRUE(std::abs(residual) < 1e-8,
                "q=0 quartic root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }

        std::vector<double> vals;
        for (auto& r : roots) vals.push_back(r->to_numeric());
        std::sort(vals.begin(), vals.end());

        int count_neg3 = 0, count_neg1 = 0;
        for (double v : vals) {
            if (std::abs(v - (-3.0)) < 1e-8) count_neg3++;
            else if (std::abs(v - (-1.0)) < 1e-8) count_neg1++;
        }
        EXPECT_TRUE(count_neg3 == 2 && count_neg1 == 2,
            "q=0 quartic roots should be {-3, -3, -1, -1}");
    }

    TEST_CASE("Solve Quartic - Repeated root (x-1)^4");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(-4);
        auto c = SymbolicExpr::number(6);
        auto d = SymbolicExpr::number(-4);
        auto e = SymbolicExpr::number(1);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "(x-1)^4 should return 4 roots");

        bool all_one = true;
        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            if (std::abs(val - 1.0) > 1e-6) {
                all_one = false;
            }
            double residual = eval_quartic(1, -4, 6, -4, 1, val);
            EXPECT_TRUE(std::abs(residual) < 1e-6,
                "(x-1)^4 root " + std::to_string(i) + " residual < 1e-6 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }
        EXPECT_TRUE(all_one, "(x-1)^4 all roots should equal 1");
    }

    TEST_CASE("Solve Quartic - Repeated roots (x-1)^2*(x-2)^2");

    {
        auto a = SymbolicExpr::number(1);
        auto b = SymbolicExpr::number(-6);
        auto c = SymbolicExpr::number(13);
        auto d = SymbolicExpr::number(-12);
        auto e = SymbolicExpr::number(4);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "(x-1)^2(x-2)^2 should return 4 roots");

        for (size_t i = 0; i < roots.size(); ++i) {
            double val = roots[i]->to_numeric();
            double residual = eval_quartic(1, -6, 13, -12, 4, val);
            EXPECT_TRUE(std::abs(residual) < 1e-6,
                "(x-1)^2(x-2)^2 root " + std::to_string(i) + " residual < 1e-6 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")");
        }

        std::vector<double> vals;
        for (auto& r : roots) vals.push_back(r->to_numeric());
        std::sort(vals.begin(), vals.end());

        int count_1 = 0, count_2 = 0;
        for (double v : vals) {
            if (std::abs(v - 1.0) < 1e-6) count_1++;
            else if (std::abs(v - 2.0) < 1e-6) count_2++;
        }
        EXPECT_TRUE(count_1 == 2 && count_2 == 2,
            "(x-1)^2(x-2)^2 roots should be {1, 1, 2, 2} (got {"
            + std::to_string(vals[0]) + ", " + std::to_string(vals[1]) + ", "
            + std::to_string(vals[2]) + ", " + std::to_string(vals[3]) + "})");
    }

    TEST_CASE("Solve Quartic - Vieta's formulas: sum and product of roots");

    {
        double ca = 2.0, cb = 0.0, cc = -20.0, cd = 0.0, ce = 18.0;
        auto a = SymbolicExpr::number(ca);
        auto b = SymbolicExpr::number(cb);
        auto c = SymbolicExpr::number(cc);
        auto d = SymbolicExpr::number(cd);
        auto e = SymbolicExpr::number(ce);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Vieta quartic should return 4 roots");

        double sum = 0.0;
        double product = 1.0;
        for (auto& r : roots) {
            double val = r->to_numeric();
            sum += val;
            product *= val;
        }

        double expected_sum = -cb / ca;
        double expected_product = ce / ca;

        EXPECT_TRUE(std::abs(sum - expected_sum) < 1e-8,
            "Vieta sum: Σrᵢ = -b/a = " + std::to_string(expected_sum) + " (got " + std::to_string(sum) + ")");
        EXPECT_TRUE(std::abs(product - expected_product) < 1e-8,
            "Vieta product: ∏rᵢ = e/a = " + std::to_string(expected_product) + " (got " + std::to_string(product) + ")");
    }

    {
        double ca = 1.0, cb = -10.0, cc = 35.0, cd = -50.0, ce = 24.0;
        auto a = SymbolicExpr::number(ca);
        auto b = SymbolicExpr::number(cb);
        auto c = SymbolicExpr::number(cc);
        auto d = SymbolicExpr::number(cd);
        auto e = SymbolicExpr::number(ce);

        auto roots = lamina::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE(roots.size() == 4, "Vieta quartic (1,2,3,4) should return 4 roots");

        double sum = 0.0;
        double product = 1.0;
        for (auto& r : roots) {
            double val = r->to_numeric();
            sum += val;
            product *= val;
        }

        double expected_sum = -cb / ca;
        double expected_product = ce / ca;

        EXPECT_TRUE(std::abs(sum - expected_sum) < 1e-8,
            "Vieta sum (1+2+3+4=10): got " + std::to_string(sum));
        EXPECT_TRUE(std::abs(product - expected_product) < 1e-8,
            "Vieta product (1*2*3*4=24): got " + std::to_string(product));
    }

    TEST_CASE("Quartic root verification (random from known roots)");

    {
        const double RESIDUAL_TOL = 1e-10;
        const int NUM_GENERAL_TRIALS = 40;
        const int NUM_BIQUADRATIC_TRIALS = 10;
        const int NUM_Q0_TRIALS = 10;
        int quartic_verify_pass_count = 0;
        int quartic_verify_total_roots = 0;
        int total_trials = NUM_GENERAL_TRIALS + NUM_BIQUADRATIC_TRIALS + NUM_Q0_TRIALS;

        std::mt19937 rng_quartic(314159);
        std::uniform_int_distribution<int> root_dist(-5, 5);

        for (int trial = 0; trial < NUM_GENERAL_TRIALS; ++trial) {
            int r1 = root_dist(rng_quartic);
            int r2 = root_dist(rng_quartic);
            int r3 = root_dist(rng_quartic);
            int r4 = root_dist(rng_quartic);

            int b_val = -(r1 + r2 + r3 + r4);
            int c_val = r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4;
            int d_val = -(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4);
            int e_val = r1*r2*r3*r4;

            auto roots = lamina::solve_quartic(
                SymbolicExpr::number(1),
                SymbolicExpr::number(b_val),
                SymbolicExpr::number(c_val),
                SymbolicExpr::number(d_val),
                SymbolicExpr::number(e_val),
                "x");

            if (roots.size() != 4) {
                std::ostringstream msg;
                msg << "General Trial " << trial
                    << " (roots=" << r1 << "," << r2 << "," << r3 << "," << r4
                    << "): expected 4 roots, got " << roots.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool trial_ok = true;
            for (size_t i = 0; i < roots.size(); ++i) {
                double val = roots[i]->to_numeric();
                if (std::isnan(val) || std::isinf(val)) {
                    continue;
                }
                quartic_verify_total_roots++;
                double residual = eval_quartic(1.0, (double)b_val, (double)c_val,
                                              (double)d_val, (double)e_val, val);
                if (std::abs(residual) >= RESIDUAL_TOL) {
                    trial_ok = false;
                    std::ostringstream msg;
                    msg << "General Trial " << trial << " root " << i
                        << " (roots=" << r1 << "," << r2 << "," << r3 << "," << r4
                        << "): |f(r)| = " << std::abs(residual) << " >= 1e-10"
                        << " (r = " << val << ")";
                    EXPECT_TRUE(false, msg.str());
                }
            }
            if (trial_ok) quartic_verify_pass_count++;
        }

        for (int trial = 0; trial < NUM_BIQUADRATIC_TRIALS; ++trial) {
            int r1 = root_dist(rng_quartic);
            int r2 = root_dist(rng_quartic);

            while (r1 == 0 && r2 == 0) {
                r1 = root_dist(rng_quartic);
                r2 = root_dist(rng_quartic);
            }

            int a_coeff = 1;
            int b_coeff = 0;
            int c_coeff = -(r1*r1 + r2*r2);
            int d_coeff = 0;
            int e_coeff = r1*r1 * r2*r2;

            auto roots = lamina::solve_quartic(
                SymbolicExpr::number(a_coeff),
                SymbolicExpr::number(b_coeff),
                SymbolicExpr::number(c_coeff),
                SymbolicExpr::number(d_coeff),
                SymbolicExpr::number(e_coeff),
                "x");

            if (roots.size() != 4) {
                std::ostringstream msg;
                msg << "Biquadratic Trial " << trial
                    << " (r1=" << r1 << ", r2=" << r2
                    << "): expected 4 roots, got " << roots.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool trial_ok = true;
            for (size_t i = 0; i < roots.size(); ++i) {
                double val = roots[i]->to_numeric();
                if (std::isnan(val) || std::isinf(val)) {
                    continue;
                }
                quartic_verify_total_roots++;
                double residual = eval_quartic((double)a_coeff, (double)b_coeff,
                                              (double)c_coeff, (double)d_coeff,
                                              (double)e_coeff, val);
                if (std::abs(residual) >= RESIDUAL_TOL) {
                    trial_ok = false;
                    std::ostringstream msg;
                    msg << "Biquadratic Trial " << trial << " root " << i
                        << " (r1=" << r1 << ", r2=" << r2
                        << "): |f(r)| = " << std::abs(residual) << " >= 1e-10"
                        << " (r = " << val << ")";
                    EXPECT_TRUE(false, msg.str());
                }
            }
            if (trial_ok) quartic_verify_pass_count++;
        }

        for (int trial = 0; trial < NUM_Q0_TRIALS; ++trial) {
            int center = root_dist(rng_quartic);
            int k = root_dist(rng_quartic);
            while (k == 0) k = root_dist(rng_quartic);

            int r1 = center, r2 = center, r3 = center + k, r4 = center - k;

            int b_val = -(r1 + r2 + r3 + r4);
            int c_val = r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4;
            int d_val = -(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4);
            int e_val = r1*r2*r3*r4;

            auto roots = lamina::solve_quartic(
                SymbolicExpr::number(1),
                SymbolicExpr::number(b_val),
                SymbolicExpr::number(c_val),
                SymbolicExpr::number(d_val),
                SymbolicExpr::number(e_val),
                "x");

            if (roots.size() != 4) {
                std::ostringstream msg;
                msg << "q=0 Trial " << trial
                    << " (center=" << center << ", k=" << k
                    << "): expected 4 roots, got " << roots.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool trial_ok = true;
            for (size_t i = 0; i < roots.size(); ++i) {
                double val = roots[i]->to_numeric();
                if (std::isnan(val) || std::isinf(val)) {
                    continue;
                }
                quartic_verify_total_roots++;
                double residual = eval_quartic(1.0, (double)b_val, (double)c_val,
                                              (double)d_val, (double)e_val, val);
                if (std::abs(residual) >= RESIDUAL_TOL) {
                    trial_ok = false;
                    std::ostringstream msg;
                    msg << "q=0 Trial " << trial << " root " << i
                        << " (center=" << center << ", k=" << k
                        << "): |f(r)| = " << std::abs(residual) << " >= 1e-10"
                        << " (r = " << val << ")";
                    EXPECT_TRUE(false, msg.str());
                }
            }
            if (trial_ok) quartic_verify_pass_count++;
        }

        {
            std::ostringstream msg;
            msg << "Quartic root verification: " << quartic_verify_pass_count
                << "/" << total_trials << " trials passed ("
                << quartic_verify_total_roots << " real roots verified)";
            EXPECT_TRUE(quartic_verify_pass_count == total_trials, msg.str());
        }
    }

    return TEST_REPORT();
}
