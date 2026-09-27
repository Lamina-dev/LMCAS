#include "limit_result.hpp"
#include "test_common.hpp"
#include "internal/visitors/limit_visitor.hpp"
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
    while (v == 0)
        v = rand_int(-max_abs, max_abs);
    return v;
}

/// Build a polynomial expression: c_n*x^n + c_{n-1}*x^{n-1} + ... + c_0
/// with given coefficients (index = power).
static std::shared_ptr<SymbolicExpr> build_polynomial(
    const std::shared_ptr<SymbolicExpr> &x,
    const std::vector<int> &coeffs) {
    std::shared_ptr<SymbolicExpr> result = nullptr;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) {
            continue;
        }
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
    if (!result) {
        result = SymbolicExpr::number(0);
    }
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

static bool matches_rational_degree_limit(
    const std::shared_ptr<SymbolicExpr> &result,
    const std::vector<int> &p_coeffs, const std::vector<int> &q_coeffs) {
    if (!result) {
        return false;
    }
    if (p_coeffs.size() > q_coeffs.size()) {
        auto text = result->to_string();
        return text.find("inf") != std::string::npos ||
               text.find("Inf") != std::string::npos ||
               text.find("∞") != std::string::npos;
    }
    const double expected = p_coeffs.size() < q_coeffs.size()
                                ? 0.0
                                : static_cast<double>(p_coeffs.back()) / q_coeffs.back();
    auto value = test_numeric_eval(result);
    if (value) {
        return std::abs(*value - expected) < 1e-6;
    }
    auto exact = p_coeffs.size() < q_coeffs.size()
                     ? SymbolicExpr::number(0)
                     : SymbolicExpr::number(Rational(p_coeffs.back(), q_coeffs.back()));
    return SymbolicExpr::add(result, SymbolicExpr::multiply(
                                         SymbolicExpr::number(-1), exact))
        ->simplify()
        ->is_zero();
}

TEST(PropLimits, RationalDegreeRule) {
    auto x = SymbolicExpr::variable("x");
    auto inf = SymbolicExpr::infinity(1);
    auto neg_one = SymbolicExpr::number(-1);

    int num_trials = 30;

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

        auto evaluated = LMCAS::limit_expression_checked(expr, "x", inf);
        ASSERT_TRUE(evaluated) << "Trial " << trial << ": " << expr->to_string()
                               << "; " << evaluated.error().operation
                               << ": " << evaluated.error().message;
        auto result = evaluated.value();

        bool property_holds = matches_rational_degree_limit(result, p_coeffs, q_coeffs);

        EXPECT_TRUE(property_holds)
            << "Trial " << trial << ": deg(P)=" << deg_p
            << " deg(Q)=" << deg_q
            << " P_lead=" << p_coeffs[deg_p]
            << " Q_lead=" << q_coeffs[deg_q]
            << " result=" << (result ? result->to_string() : "null");
    }
}

TEST(PropLimits, SqueezeBoundedOscillation) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto neg_one = SymbolicExpr::number(-1);

    int num_trials = 30;

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

        EXPECT_TRUE(property_holds)
            << "Trial " << trial << ": x^" << n
            << " * bounded(" << c << "/x^" << m << ")"
            << " result=" << (result ? result->to_string() : "null");
        if (!property_holds) {
            std::string result_str = result ? result->to_string() : "null";
            std::string func_name = (rand_int(0, 1) == 0) ? "sin" : "cos";
            std::cerr << "[INFO] Trial " << trial << " failed: x^" << n
                      << " * bounded(" << c << "/x^" << m << ")"
                      << " result=" << result_str << std::endl;
        }
    }
}
