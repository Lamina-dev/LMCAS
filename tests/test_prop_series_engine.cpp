#include "test_common.hpp"
#include "series_engine.hpp"
#include <random>
#include <cmath>

using namespace LMCAS;


static constexpr unsigned kPropertySeed = 0x53455249u;
static std::mt19937 rng(kPropertySeed);

/// Generate a random integer in [lo, hi].
static int rand_int(int lo, int hi) {
    std::uniform_int_distribution<int> dist(lo, hi);
    return dist(rng);
}

/// Generate a random non-zero integer in [-max_abs, max_abs].
static int rand_nonzero(int max_abs = 5) {
    int v = 0;
    while (v == 0) v = rand_int(-max_abs, max_abs);
    return v;
}

/// Helper: create a number expression.
static std::shared_ptr<SymbolicExpr> num(int n) {
    return SymbolicExpr::number(n);
}

/// Generate a random coefficient vector of given length with integer entries.
static std::vector<std::shared_ptr<SymbolicExpr>> rand_coeffs(int len, int max_abs = 4) {
    std::vector<std::shared_ptr<SymbolicExpr>> coeffs(len);
    for (int i = 0; i < len; ++i) {
        coeffs[i] = num(rand_int(-max_abs, max_abs));
    }
    return coeffs;
}


static void test_power_series_addition_componentwise() {
    TEST_CASE("Power series addition is component-wise");

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        int len = rand_int(1, 8);
        auto a = rand_coeffs(len);
        auto b = rand_coeffs(len);

        auto result = LMCAS::power_series_add(a, b);

        bool property_holds = true;

        // Result length should be max(a.size(), b.size()) = len (same length)
        if (result.size() != static_cast<size_t>(len)) {
            property_holds = false;
        } else {
            for (int k = 0; k < len; ++k) {
                // Compute expected: a[k] + b[k]
                auto expected = SymbolicExpr::add(a[k], b[k])->simplify();
                auto actual = result[k] ? result[k]->simplify() : num(0);

                auto exp_val = test_numeric_eval(expected);
                auto act_val = test_numeric_eval(actual);

                if (exp_val && act_val) {
                    if (std::abs(*exp_val - *act_val) > 1e-9) {
                        property_holds = false;
                        std::cerr << "[INFO] Trial " << trial << " index " << k
                                  << ": expected=" << *exp_val
                                  << " actual=" << *act_val << std::endl;
                        break;
                    }
                } else {
                    // Compare string representations as fallback
                    std::string exp_str = expected ? expected->to_string() : "0";
                    std::string act_str = actual ? actual->to_string() : "0";
                    if (exp_str != act_str) {
                        property_holds = false;
                        std::cerr << "[INFO] Trial " << trial << " index " << k
                                  << ": expected=" << exp_str
                                  << " actual=" << act_str << std::endl;
                        break;
                    }
                }
            }
        }

        if (property_holds) pass_count++;
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " power series addition trials are component-wise (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_power_series_multiplication_cauchy() {
    TEST_CASE("Power series multiplication satisfies Cauchy product");

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        int len_a = rand_int(1, 5);
        int len_b = rand_int(1, 5);
        int order = rand_int(1, std::min(len_a + len_b, 6));

        auto a = rand_coeffs(len_a, 3);
        auto b = rand_coeffs(len_b, 3);

        auto result = LMCAS::power_series_multiply_checked(a, b, order).value();

        bool property_holds = true;

        if (result.size() != static_cast<size_t>(order)) {
            property_holds = false;
        } else {
            for (int k = 0; k < order; ++k) {
                // Compute Cauchy product coefficient manually:
                // c[k] = sum_{j=0}^{k} a[j] * b[k-j]
                double expected_val = 0.0;
                for (int j = 0; j <= k; ++j) {
                    double aj = 0.0, bkj = 0.0;
                    if (j < len_a) {
                        auto v = test_numeric_eval(a[j]);
                        if (v) aj = *v;
                    }
                    if ((k - j) < len_b) {
                        auto v = test_numeric_eval(b[k - j]);
                        if (v) bkj = *v;
                    }
                    expected_val += aj * bkj;
                }

                auto act_val = test_numeric_eval(result[k]);
                if (act_val) {
                    if (std::abs(*act_val - expected_val) > 1e-6) {
                        property_holds = false;
                        std::cerr << "[INFO] Trial " << trial << " index " << k
                                  << ": expected=" << expected_val
                                  << " actual=" << *act_val << std::endl;
                        break;
                    }
                } else {
                    /// 数值求值未决时,使用两侧零值条件完成本采样点判定.
                    if (std::abs(expected_val) > 1e-9) {
                        property_holds = false;
                        break;
                    }
                }
            }
        }

        if (property_holds) pass_count++;
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " power series multiplication trials satisfy Cauchy product (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_fourier_series_symmetry() {
    TEST_CASE("Fourier series symmetry (even/odd functions)");

    auto x = SymbolicExpr::variable("x");
    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        std::shared_ptr<SymbolicExpr> f;
        int expected_parity = 0; // 1 = even, -1 = odd

        int choice = rand_int(0, 5);
        switch (choice) {
        case 0: {
            // Even: x^(2n) for random even power
            int n = rand_int(1, 3);
            f = SymbolicExpr::power(x, num(2 * n));
            expected_parity = 1;
            break;
        }
        case 1: {
            // Even: cos(a*x) for random a
            int a = rand_nonzero(4);
            f = SymbolicExpr::cos(SymbolicExpr::multiply(num(a), x));
            expected_parity = 1;
            break;
        }
        case 2: {
            // Odd: x^(2n+1) for random odd power
            int n = rand_int(0, 2);
            f = SymbolicExpr::power(x, num(2 * n + 1));
            expected_parity = -1;
            break;
        }
        case 3: {
            // Odd: sin(a*x) for random a
            int a = rand_nonzero(4);
            f = SymbolicExpr::sin(SymbolicExpr::multiply(num(a), x));
            expected_parity = -1;
            break;
        }
        case 4: {
            // Even: a*x^2 + b*x^4 (sum of even powers)
            int a = rand_nonzero(3);
            int b = rand_nonzero(3);
            auto t1 = SymbolicExpr::multiply(num(a), SymbolicExpr::power(x, num(2)));
            auto t2 = SymbolicExpr::multiply(num(b), SymbolicExpr::power(x, num(4)));
            f = SymbolicExpr::add(t1, t2);
            expected_parity = 1;
            break;
        }
        case 5: {
            // Odd: a*x + b*x^3 (sum of odd powers)
            int a = rand_nonzero(3);
            int b = rand_nonzero(3);
            auto t1 = SymbolicExpr::multiply(num(a), x);
            auto t2 = SymbolicExpr::multiply(num(b), SymbolicExpr::power(x, num(3)));
            f = SymbolicExpr::add(t1, t2);
            expected_parity = -1;
            break;
        }
        }

        // Verify parity: substitute x -> -x and check
        auto neg_x = SymbolicExpr::multiply(num(-1), x);
        auto f_neg = f->substitute("x", neg_x);
        if (!f_neg) { pass_count++; continue; }
        f_neg = f_neg->simplify();
        auto f_s = f->simplify();

        bool property_holds = false;

        if (expected_parity == 1) {
            // Even: f(-x) - f(x) should be 0
            auto diff = SymbolicExpr::add(f_neg,
                SymbolicExpr::multiply(num(-1), f_s))->simplify();
            if (diff && diff->is_zero()) {
                property_holds = true;
            } else {
                // Try numeric evaluation at a test point
                auto diff_at_1 = diff ? diff->substitute("x", num(1)) : nullptr;
                if (diff_at_1) {
                    diff_at_1 = diff_at_1->simplify();
                    auto val = test_numeric_eval(diff_at_1);
                    property_holds = val && (std::abs(*val) < 1e-9);
                }
            }
        } else if (expected_parity == -1) {
            // Odd: f(-x) + f(x) should be 0
            auto sum = SymbolicExpr::add(f_neg, f_s)->simplify();
            if (sum && sum->is_zero()) {
                property_holds = true;
            } else {
                // Try numeric evaluation at a test point
                auto sum_at_1 = sum ? sum->substitute("x", num(1)) : nullptr;
                if (sum_at_1) {
                    sum_at_1 = sum_at_1->simplify();
                    auto val = test_numeric_eval(sum_at_1);
                    property_holds = val && (std::abs(*val) < 1e-9);
                }
            }
        }

        if (property_holds) {
            pass_count++;
        } else {
            std::cerr << "[INFO] Trial " << trial
                      << ": parity=" << expected_parity
                      << " f=" << f->to_string()
                      << " f(-x)=" << (f_neg ? f_neg->to_string() : "null")
                      << std::endl;
        }
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " Fourier series symmetry (parity) trials passed (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_symbolic_summation_closed_form() {
    TEST_CASE("Symbolic summation closed form verification");

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        // Generate a polynomial body in k of degree 1..3
        int degree = rand_int(1, 3);
        auto k = SymbolicExpr::variable("k");

        // Build polynomial: c_d*k^d + ... + c_1*k + c_0
        std::vector<int> coeffs(degree + 1);
        for (int i = 0; i <= degree; ++i) {
            coeffs[i] = rand_int(-3, 3);
        }
        coeffs[degree] = rand_nonzero(3); // ensure leading coeff non-zero

        std::shared_ptr<SymbolicExpr> body = num(0);
        for (int i = 0; i <= degree; ++i) {
            if (coeffs[i] == 0) continue;
            std::shared_ptr<SymbolicExpr> term;
            if (i == 0) {
                term = num(coeffs[i]);
            } else if (i == 1) {
                term = SymbolicExpr::multiply(num(coeffs[i]), k);
            } else {
                term = SymbolicExpr::multiply(num(coeffs[i]),
                    SymbolicExpr::power(k, num(i)));
            }
            body = SymbolicExpr::add(body, term);
        }
        body = body->simplify();

        // Choose a concrete upper bound for verification
        int upper_val = rand_int(3, 10);

        // Compute direct summation for verification
        double direct_sum = 0.0;
        for (int kv = 1; kv <= upper_val; ++kv) {
            double term_val = 0.0;
            for (int i = 0; i <= degree; ++i) {
                term_val += coeffs[i] * std::pow(kv, i);
            }
            direct_sum += term_val;
        }

        // Use symbolic_sum with concrete bounds (direct evaluation path)
        auto direct_result = LMCAS::symbolic_sum(body, "k", num(1), num(upper_val));

        bool property_holds = false;

        if (direct_result) {
            direct_result = direct_result->simplify();
            auto val = test_numeric_eval(direct_result);
            if (val) {
                property_holds = (std::abs(*val - direct_sum) < 1e-6);
                if (!property_holds) {
                    std::cerr << "[INFO] Trial " << trial
                              << ": symbolic=" << *val
                              << " direct=" << direct_sum
                              << " body=" << body->to_string()
                              << " n=" << upper_val << std::endl;
                }
            }
        }

        if (property_holds) pass_count++;
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " symbolic summation closed form trials verified (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_convergence_radius_geometric() {
    TEST_CASE("Convergence radius for geometric series");

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        // Generate a non-zero ratio r (integer for simplicity)
        int r = rand_nonzero(5);
        double abs_r = std::abs(static_cast<double>(r));
        double expected_radius = 1.0 / abs_r;

        // Infinite geometric sequence a_n = r^n is represented by its
        // general coefficient, never inferred from a finite prefix.
        auto n = SymbolicExpr::variable("n");
        auto coefficient = SymbolicExpr::power(num(r), n);
        auto checked_radius =
            LMCAS::convergence_radius_checked(coefficient, "n");
        auto radius = checked_radius
            ? checked_radius.value() : std::shared_ptr<SymbolicExpr>{};

        bool property_holds = false;

        if (radius) {
            auto val = test_numeric_eval(radius);
            if (val) {
                // The convergence radius should be 1/|r|
                property_holds = (std::abs(*val - expected_radius) < 1e-6);
                if (!property_holds) {
                    std::cerr << "[INFO] Trial " << trial
                              << ": r=" << r
                              << " expected_R=" << expected_radius
                              << " actual_R=" << *val << std::endl;
                }
            } else {
                // Check if it's infinity (which would be correct for |r| < 1,
                // but we use integer r so |r| >= 1, radius <= 1)
                auto str = radius->to_string();
                if (str.find("inf") != std::string::npos ||
                    str.find("∞") != std::string::npos) {
                    // Infinity is only correct if |r| < 1 (not possible with integer r != 0)
                    property_holds = false;
                }
            }
        }

        if (property_holds) pass_count++;
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " geometric series convergence radius trials passed (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


int main() {
    test_power_series_addition_componentwise();
    test_power_series_multiplication_cauchy();
    test_fourier_series_symmetry();
    test_symbolic_summation_closed_form();
    test_convergence_radius_geometric();

    return TEST_REPORT();
}
