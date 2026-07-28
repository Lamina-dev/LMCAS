// Feature: integration-enhancements, Property 4: Trigonometric sin^m·cos^n round-trip
//
// Validates: Requirements 3.1, 3.9
//
// Property 4: For all integer pairs (m, n) with m >= 0, n >= 0, and
// m + n <= 8, integrating sin^m(x)·cos^n(x) and differentiating the
// result SHALL yield an expression numerically equal to sin^m(x)·cos^n(x)
// at sample points x in {0.5, 1.0, 1.5, 2.0, 2.5} within tolerance 1e-10.
//
// Approach
// --------
//   * Iterate (m, n) over all 45 integer pairs with m >= 0, n >= 0,
//     m + n <= 8 (1 + 2 + ... + 9 = 45 pairs).
//
//   * For each pair, build the integrand sin^m(x)·cos^n(x) using the
//     simplest AST form available:
//       - m = n = 0  -> 1
//       - m = k, n = 0 -> sin(x) for k=1, sin(x)^k for k>=2
//       - m = 0, n = k -> cos(x) for k=1, cos(x)^k for k>=2
//       - m >= 1, n >= 1 -> sin(x)^m * cos(x)^n
//     (If the exponent is exactly 1 we skip the PowerNode wrapper so the
//     expression stays in the canonical shape the trig-combination
//     strategy expects.)
//
//   * Call Integrator::integrate(integrand, "x") to obtain a closed form.
//     If the integrator returns an unevaluated integral node for any
//     pair within the property's scope, that's a failure of Requirement
//     3.1 and the test reports it as such.
//
//   * Differentiate the result symbolically.
//
//   * Numerically evaluate both the integrand and the derivative at the
//     sample points x in {0.5, 1.0, 1.5, 2.0, 2.5} and compare within
//     tolerance 1e-10. test_numeric_eval supports Sin/Cos/Exp/Ln/Sqrt/Abs
//     plus the algebraic node kinds, which is sufficient for the
//     antiderivatives produced by TrigCombinationStrategy (multiples of
//     sin(k*x) / cos(k*x) and rational combinations thereof).
//
//   * The test FAILS if any pair produces a numeric mismatch above
//     tolerance, an unevaluated integral, or a sample point whose values
//     fail to evaluate (which would indicate the result contains a
//     function the evaluator cannot handle and therefore the round-trip
//     cannot be verified).

#include "test_common.hpp"
#include "integration.hpp"

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using lamina::Integrator;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-10;

// ----- AST construction ---------------------------------------------------

std::shared_ptr<SymbolicExpr> num_int(long long n) {
    return SymbolicExpr::number(n);
}

// Build f(x)^k, but if k == 1 return f(x) unwrapped so the integrand stays
// in the exact shape the TrigCombinationStrategy detects.
std::shared_ptr<SymbolicExpr> pow_or_self(std::shared_ptr<SymbolicExpr> base, int k) {
    if (k == 1) return base;
    return SymbolicExpr::power(base, num_int(k));
}

// Build the integrand sin(x)^m * cos(x)^n in the simplest shape.
std::shared_ptr<SymbolicExpr> build_integrand(int m, int n) {
    auto x_var = SymbolicExpr::variable(kVarName);
    if (m == 0 && n == 0) {
        return num_int(1);
    }
    if (n == 0) {
        return pow_or_self(SymbolicExpr::sin(x_var), m);
    }
    if (m == 0) {
        return pow_or_self(SymbolicExpr::cos(x_var), n);
    }
    auto sin_part = pow_or_self(SymbolicExpr::sin(x_var), m);
    auto cos_part = pow_or_self(SymbolicExpr::cos(x_var), n);
    return SymbolicExpr::multiply(sin_part, cos_part);
}

bool has_integral_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == FunctionNode::FuncType::Calculus_Integral) return true;
        for (auto& a : fn->arguments())
            if (has_integral_node(a)) return true;
    } else if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands())
            if (has_integral_node(op)) return true;
    } else if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands())
            if (has_integral_node(op)) return true;
    } else if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (has_integral_node(pow->base())) return true;
        if (has_integral_node(pow->exponent())) return true;
    }
    return false;
}

const std::vector<double>& sample_points() {
    static const std::vector<double> S = {0.5, 1.0, 1.5, 2.0, 2.5};
    return S;
}

// ----- Per-pair check -----------------------------------------------------

struct PairReport {
    bool failed = false;
    int  matches = 0;
    std::string detail;
};

PairReport verify_pair(int m, int n) {
    PairReport rep;

    auto integrand = build_integrand(m, n);
    if (!integrand) {
        rep.failed = true;
        rep.detail = "integrand build returned null";
        return rep;
    }

    Integrator integ;
    std::shared_ptr<SymbolicExpr> result;
    try {
        result = lamina::detail::make_expression_ptr(integ.integrate(*integrand, kVarName));
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception during integration: ") + e.what();
        return rep;
    }

    if (has_integral_node(lamina::detail::node(result))) {
        rep.failed = true;
        rep.detail = "integrator returned unevaluated integral: " + result->to_string();
        return rep;
    }

    auto deriv = result->differentiate(kVarName);
    if (!deriv) {
        rep.failed = true;
        rep.detail = "differentiation returned null";
        return rep;
    }
    auto deriv_simp = deriv->simplify();
    if (!deriv_simp) deriv_simp = deriv;

    auto integrand_simp = integrand->simplify();
    if (!integrand_simp) integrand_simp = integrand;

    for (double xv : sample_points()) {
        auto x_val = SymbolicExpr::number(xv);
        auto integrand_at = integrand_simp->substitute(kVarName, x_val);
        auto deriv_at = deriv_simp->substitute(kVarName, x_val);
        if (!integrand_at || !deriv_at) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv << ": substitute returned null"
                << " | result=" << result->to_string();
            rep.detail = oss.str();
            return rep;
        }
        integrand_at = integrand_at->simplify();
        deriv_at = deriv_at->simplify();

        auto pv = test_numeric_eval(integrand_at);
        auto dv = test_numeric_eval(deriv_at);
        if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv << ": numeric evaluation failed"
                << " | result=" << result->to_string()
                << " | derivative=" << deriv_simp->to_string();
            rep.detail = oss.str();
            return rep;
        }
        double delta = std::abs(*pv - *dv);
        if (delta > kTolerance) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv
                << ": integrand=" << *pv
                << " vs d/dx(result)=" << *dv
                << " |delta|=" << delta
                << " | result=" << result->to_string();
            rep.detail = oss.str();
            return rep;
        }
        ++rep.matches;
    }
    return rep;
}

}// anonymous namespace

int main() {
    TEST_CASE("Property 4: Trigonometric sin^m * cos^n round-trip");

    int total_pairs = 0;
    int verified_pairs = 0;
    int failed_pairs = 0;
    int total_matches = 0;

    for (int m = 0; m <= 8; ++m) {
        for (int n = 0; n + m <= 8; ++n) {
            ++total_pairs;
            PairReport rep = verify_pair(m, n);

            std::ostringstream label;
            label << "sin^" << m << "(x) * cos^" << n << "(x)";

            if (rep.failed) {
                ++failed_pairs;
                std::cerr << "[FAIL] " << label.str() << " : " << rep.detail << std::endl;
                EXPECT_TRUE(false, label.str() + ": numeric round-trip mismatch");
            } else {
                ++verified_pairs;
                total_matches += rep.matches;
                std::ostringstream oss;
                oss << label.str() << ": " << rep.matches << " sample point(s) matched";
                EXPECT_TRUE(!rep.failed && rep.matches > 0, oss.str());
            }
        }
    }

    std::cout << "\nSummary:"
              << " total_pairs="     << total_pairs
              << " verified_pairs="  << verified_pairs
              << " failed_pairs="    << failed_pairs
              << " total_matches="   << total_matches
              << std::endl;

    EXPECT_TRUE(verified_pairs == 45,
                "all 45 (m, n) pairs with m+n<=8 verified by numeric round-trip "
                "(verified=" + std::to_string(verified_pairs) + ", expected=45)");

    return TEST_REPORT();
}
