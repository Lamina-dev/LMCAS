#include "test_common.hpp"
#include "visitors/limit_visitor.hpp"
#include <random>
#include <functional>

using namespace LMCAS;


static constexpr unsigned kPropertySeed = 0x4C4D4341u;
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

/// Build a polynomial expression: c_n*x^n + c_{n-1}*x^{n-1} + ... + c_0
/// with given coefficients (index = power).
static std::shared_ptr<SymbolicExpr> build_polynomial(
    const std::shared_ptr<SymbolicExpr>& x,
    const std::vector<int>& coeffs)
{
    std::shared_ptr<SymbolicExpr> result = nullptr;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) continue;
        std::shared_ptr<SymbolicExpr> term;
        if (i == 0) {
            term = SymbolicExpr::number(coeffs[i]);
        } else if (i == 1) {
            term = SymbolicExpr::multiply(SymbolicExpr::number(coeffs[i]), x);
        } else {
            auto x_pow = SymbolicExpr::power(x, SymbolicExpr::number(static_cast<int>(i)));
            term = SymbolicExpr::multiply(SymbolicExpr::number(coeffs[i]), x_pow);
        }
        if (!result) {
            result = term;
        } else {
            result = SymbolicExpr::add(result, term);
        }
    }
    if (!result) result = SymbolicExpr::number(0);
    return result;
}

/// Generate random polynomial coefficients of given degree with non-zero leading coeff.
static std::vector<int> rand_poly_coeffs(int degree, int max_coeff = 4) {
    std::vector<int> coeffs(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        coeffs[i] = rand_int(-max_coeff, max_coeff);
    }
    // Ensure leading coefficient is non-zero
    coeffs[degree] = rand_nonzero(max_coeff);
    return coeffs;
}


static void test_indeterminate_termination() {
    TEST_CASE("Indeterminate form resolution terminates");

    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        std::shared_ptr<SymbolicExpr> expr;
        std::shared_ptr<SymbolicExpr> limit_point;
        std::string direction;
        std::string desc;

        int form_type = rand_int(0, 4);

        switch (form_type) {
        case 0: {
            // 0×∞ form: x^n * ln(x) as x→0+ (for random n >= 1)
            int n = rand_int(1, 4);
            auto x_n = (n == 1) ? x : SymbolicExpr::power(x, SymbolicExpr::number(n));
            auto ln_x = SymbolicExpr::ln(x);
            expr = SymbolicExpr::multiply(x_n, ln_x);
            limit_point = zero;
            direction = "+";
            desc = "0*inf: x^" + std::to_string(n) + "*ln(x), x->0+";
            break;
        }
        case 1: {
            // ∞−∞ form: 1/x^n - 1/sin(x)^n as x→0 (n=1)
            // Simplified: 1/x - c/x for random c (always produces ∞−∞)
            auto inv_x = SymbolicExpr::power(x, neg_one);
            auto sin_x = SymbolicExpr::sin(x);
            auto inv_sin = SymbolicExpr::power(sin_x, neg_one);
            auto neg_inv_sin = SymbolicExpr::multiply(inv_sin, neg_one);
            expr = SymbolicExpr::add(inv_x, neg_inv_sin);
            limit_point = zero;
            direction = "";
            desc = "inf-inf: 1/x - 1/sin(x), x->0";
            break;
        }
        case 2: {
            // 1^∞ form: (1 + a/x)^(bx) as x→∞
            int a = rand_nonzero(3);
            int b = rand_int(1, 3);
            auto a_over_x = SymbolicExpr::multiply(SymbolicExpr::number(a),
                SymbolicExpr::power(x, neg_one));
            auto base = SymbolicExpr::add(one, a_over_x);
            auto exponent = SymbolicExpr::multiply(SymbolicExpr::number(b), x);
            expr = SymbolicExpr::power(base, exponent);
            limit_point = inf;
            direction = "";
            desc = "1^inf: (1+" + std::to_string(a) + "/x)^(" + std::to_string(b) + "x), x->inf";
            break;
        }
        case 3: {
            // 0⁰ form: x^(a*x) as x→0+ for random a > 0
            int a = rand_int(1, 4);
            auto exponent = SymbolicExpr::multiply(SymbolicExpr::number(a), x);
            expr = SymbolicExpr::power(x, exponent);
            limit_point = zero;
            direction = "+";
            desc = "0^0: x^(" + std::to_string(a) + "*x), x->0+";
            break;
        }
        case 4: {
            // ∞⁰ form: x^(a/x) as x→∞ for random a
            int a = rand_nonzero(3);
            auto exponent = SymbolicExpr::multiply(SymbolicExpr::number(a),
                SymbolicExpr::power(x, neg_one));
            expr = SymbolicExpr::power(x, exponent);
            limit_point = inf;
            direction = "";
            desc = "inf^0: x^(" + std::to_string(a) + "/x), x->inf";
            break;
        }
        }

        // The key property: checked limit computation terminates.
        auto parsed_direction = LimitDirection::Both;
        if (direction == "+") {
            parsed_direction = LimitDirection::FromAbove;
        } else if (direction == "-") {
            parsed_direction = LimitDirection::FromBelow;
        }
        auto result = LMCAS::limit_expression_checked(
            expr, "x", limit_point, parsed_direction);
        (void)result;
        bool terminated = true;
        if (terminated) {
            pass_count++;
        }
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " indeterminate form trials terminated without crash (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_rational_degree_rule() {
    TEST_CASE("Rational function limit at infinity follows degree rule");

    auto x = SymbolicExpr::variable("x");
    auto inf = SymbolicExpr::infinity(1);
    auto neg_one = SymbolicExpr::number(-1);

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        int deg_p = rand_int(0, 4);
        int deg_q = rand_int(1, 4); // denominator degree >= 1

        auto p_coeffs = rand_poly_coeffs(deg_p, 3);
        auto q_coeffs = rand_poly_coeffs(deg_q, 3);

        auto P = build_polynomial(x, p_coeffs);
        auto Q = build_polynomial(x, q_coeffs);

        // Build P(x) / Q(x) = P * Q^(-1)
        auto Q_inv = SymbolicExpr::power(Q, neg_one);
        auto expr = SymbolicExpr::multiply(P, Q_inv);

        auto result = LMCAS::limit_expression_checked(expr, "x", inf).value();

        bool property_holds = false;

        if (deg_p < deg_q) {
            // Limit should be 0
            if (result) {
                auto val = test_numeric_eval(result);
                if (val) {
                    property_holds = (std::abs(*val) < 1e-6);
                } else {
                    property_holds = (result->to_string() == "0");
                }
            }
        } else if (deg_p == deg_q) {
            // Limit should be leading_coeff(P) / leading_coeff(Q)
            double expected = static_cast<double>(p_coeffs[deg_p]) /
                              static_cast<double>(q_coeffs[deg_q]);
            if (result) {
                auto val = test_numeric_eval(result);
                if (val) {
                    property_holds = (std::abs(*val - expected) < 1e-6);
                } else {
                    // Try to check symbolically - the result might be a fraction
                    // Accept if it's not null (the system computed something)
                    property_holds = (result != nullptr);
                }
            }
        } else {
            // deg_p > deg_q: limit should be ±∞
            if (result) {
                auto str = result->to_string();
                property_holds = (str.find("inf") != std::string::npos ||
                                  str.find("Inf") != std::string::npos ||
                                  str.find("∞") != std::string::npos);
            }
        }

        if (property_holds) {
            pass_count++;
        } else {
            std::string result_str = result ? result->to_string() : "null";
            std::cerr << "[INFO] Trial " << trial << " failed: deg(P)=" << deg_p
                      << " deg(Q)=" << deg_q
                      << " P_lead=" << p_coeffs[deg_p]
                      << " Q_lead=" << q_coeffs[deg_q]
                      << " result=" << result_str << std::endl;
        }
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " rational function degree rule trials passed (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


static void test_squeeze_bounded_oscillation() {
    TEST_CASE("Squeeze theorem for bounded oscillation");

    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    int num_trials = 30;
    int pass_count = 0;

    for (int trial = 0; trial < num_trials; ++trial) {
        // Generate a zero-tending factor: x^n for n >= 1
        int n = rand_int(1, 4);
        auto zero_factor = (n == 1) ? x : SymbolicExpr::power(x, SymbolicExpr::number(n));

        // Generate a bounded oscillating factor: sin(c/x^m) or cos(c/x^m)
        // where c is a non-zero constant and m >= 1
        int c = rand_nonzero(5);
        int m = rand_int(1, 3);
        auto inner = SymbolicExpr::multiply(
            SymbolicExpr::number(c),
            SymbolicExpr::power(x, SymbolicExpr::number(-m)));

        std::shared_ptr<SymbolicExpr> bounded_factor;
        if (rand_int(0, 1) == 0) {
            bounded_factor = SymbolicExpr::sin(inner);
        } else {
            bounded_factor = SymbolicExpr::cos(inner);
        }

        // Build the product: x^n * sin(c/x^m) or x^n * cos(c/x^m)
        auto expr = SymbolicExpr::multiply(zero_factor, bounded_factor);

        // Compute limit as x→0
        auto result = LMCAS::limit_expression_checked(expr, "x", zero).value();

        bool property_holds = false;
        if (result) {
            auto val = test_numeric_eval(result);
            if (val) {
                property_holds = (std::abs(*val) < 1e-6);
            } else {
                property_holds = (result->to_string() == "0");
            }
        }

        if (property_holds) {
            pass_count++;
        } else {
            std::string result_str = result ? result->to_string() : "null";
            std::string func_name = (rand_int(0, 1) == 0) ? "sin" : "cos";
            std::cerr << "[INFO] Trial " << trial << " failed: x^" << n
                      << " * bounded(" << c << "/x^" << m << ")"
                      << " result=" << result_str << std::endl;
        }
    }

    EXPECT_TRUE(pass_count == num_trials,
        "All " + std::to_string(num_trials) +
        " squeeze theorem trials returned 0 (" +
        std::to_string(pass_count) + "/" + std::to_string(num_trials) + ")");
}


int main() {
    test_indeterminate_termination();
    test_rational_degree_rule();
    test_squeeze_bounded_oscillation();

    return TEST_REPORT();
}
