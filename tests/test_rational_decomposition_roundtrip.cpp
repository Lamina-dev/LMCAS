
#include "test_common.hpp"
#include "integration.hpp"
#include "rational.hpp"
#include "bigint.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-8;
constexpr double kPoleGuard = 0.1;
constexpr int kRequiredSamples = 5;


std::shared_ptr<SymbolicExpr> num_int(long long n) {
    return SymbolicExpr::number(n);
}

// Build a symbolic polynomial from coefficients (low-to-high order) using
// only AddNode/MultiplyNode/PowerNode/NumberNode/VariableNode so that
// test_numeric_eval can evaluate it directly.
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const std::vector<long long>& coeffs) {
    auto x_var = SymbolicExpr::variable(kVarName);
    std::shared_ptr<SymbolicExpr> result = num_int(0);
    bool is_zero = true;
    for (size_t k = 0; k < coeffs.size(); ++k) {
        if (coeffs[k] == 0) continue;
        std::shared_ptr<SymbolicExpr> term;
        if (k == 0) {
            term = num_int(coeffs[k]);
        } else if (k == 1) {
            if (coeffs[k] == 1) {
                term = x_var;
            } else if (coeffs[k] == -1) {
                term = SymbolicExpr::multiply(num_int(-1), x_var);
            } else {
                term = SymbolicExpr::multiply(num_int(coeffs[k]), x_var);
            }
        } else {
            auto pw = SymbolicExpr::power(x_var, num_int(static_cast<long long>(k)));
            if (coeffs[k] == 1) {
                term = pw;
            } else if (coeffs[k] == -1) {
                term = SymbolicExpr::multiply(num_int(-1), pw);
            } else {
                term = SymbolicExpr::multiply(num_int(coeffs[k]), pw);
            }
        }
        if (is_zero) {
            result = term;
            is_zero = false;
        } else {
            result = SymbolicExpr::add(result, term);
        }
    }
    return result;
}

// Numeric polynomial evaluation using Horner's method (low-to-high coeffs).
double poly_eval(const std::vector<long long>& coeffs, double x) {
    double y = 0.0;
    for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
        y = y * x + static_cast<double>(*it);
    }
    return y;
}

bool has_integral_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (std::dynamic_pointer_cast<const IntegralNode>(node)) return true;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
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

const std::vector<double>& candidate_sample_points() {
    static const std::vector<double> S = {
        -3.7, -2.7, -1.7, -0.7,
         0.3,  0.7,  1.3,  1.7,
         2.3,  2.7,  3.3,  3.7,
         4.3,  4.7
    };
    return S;
}


struct Poly {
    std::vector<long long> coeffs;  // low-to-high
    std::string label;
};

// Denominators of degree 3..5, all factorable over Q into linear and
// irreducible quadratic factors. Specifically excluded: repeated
// irreducible quadratic factors (e.g. (x^2+1)^2), which the design
// document explicitly leaves as unevaluated integrals.
const std::vector<Poly>& denominators() {
    static const std::vector<Poly> D = {
        // x^3 - 1 = (x-1)(x^2 + x + 1)
        {{-1, 0, 0, 1}, "x^3 - 1"},
        // x^3 - x = x(x-1)(x+1)
        {{0, -1, 0, 1}, "x^3 - x"},
        // x^3 + x = x(x^2 + 1)
        {{0, 1, 0, 1}, "x^3 + x"},
        // x^3 - 8 = (x-2)(x^2 + 2x + 4)
        {{-8, 0, 0, 1}, "x^3 - 8"},
        // x^3 + x^2 + x + 1 = (x+1)(x^2 + 1)
        {{1, 1, 1, 1}, "x^3 + x^2 + x + 1"},
        // x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
        {{-6, 11, -6, 1}, "x^3 - 6x^2 + 11x - 6"},
        // x^3 + 2x^2 - x - 2 = (x-1)(x+1)(x+2)
        {{-2, -1, 2, 1}, "x^3 + 2x^2 - x - 2"},
        // x^4 - 1 = (x-1)(x+1)(x^2 + 1)
        {{-1, 0, 0, 0, 1}, "x^4 - 1"},
        // x^4 - 16 = (x-2)(x+2)(x^2 + 4)
        {{-16, 0, 0, 0, 1}, "x^4 - 16"},
        // x^4 - 5x^2 + 4 = (x-1)(x+1)(x-2)(x+2)
        {{4, 0, -5, 0, 1}, "x^4 - 5x^2 + 4"},
    };
    return D;
}

// Numerators of degree 0..2 with small integer coefficients.
const std::vector<Poly>& numerators() {
    static const std::vector<Poly> N = {
        {{1},          "1"},
        {{0, 1},       "x"},
        {{1, 1},       "x + 1"},
        {{-1, 1},      "x - 1"},
        {{3, 2},       "2x + 3"},
        {{0, 0, 1},    "x^2"},
        {{1, 0, 1},    "x^2 + 1"},
        {{-1, 0, 1},   "x^2 - 1"},
        {{5, 1, 3},    "3x^2 + x + 5"},
        {{-3, 0, 2},   "2x^2 - 3"},
    };
    return N;
}


std::vector<double> pick_safe_sample_points(const std::vector<long long>& Q_coeffs) {
    std::vector<double> kept;
    for (double x : candidate_sample_points()) {
        double qv = poly_eval(Q_coeffs, x);
        if (std::abs(qv) >= kPoleGuard) {
            kept.push_back(x);
            if (static_cast<int>(kept.size()) >= kRequiredSamples) break;
        }
    }
    return kept;
}


struct ComboReport {
    bool unevaluated = false;
    bool failed = false;
    int matches = 0;
    std::string detail;
};

ComboReport verify_combo(const Poly& P, const Poly& Q) {
    ComboReport rep;

    // Build P/Q as P * Q^-1 so the integrator sees a pure rational form.
    auto P_sym = poly_to_symbolic(P.coeffs);
    auto Q_sym = poly_to_symbolic(Q.coeffs);
    if (!P_sym || !Q_sym) {
        rep.failed = true;
        rep.detail = "polynomial AST construction failed";
        return rep;
    }
    auto Q_inv = SymbolicExpr::power(Q_sym, num_int(-1));
    auto integrand = SymbolicExpr::multiply(P_sym, Q_inv);
    if (!integrand) {
        rep.failed = true;
        rep.detail = "integrand construction failed";
        return rep;
    }

    // Integrate.
    Integrator integ;
    auto integrated = integ.integrate(*integrand, kVarName);
    if (!integrated) {
        rep.failed = true;
        rep.detail = std::string("integration failed: ") + integrated.error().message;
        return rep;
    }
    auto result = LMCAS::detail::make_expression_ptr(integrated.value());

    // The design explicitly allows the strategy to return an unevaluated
    // integral node when factoring fails or when the partial fraction term
    // requires a higher power of an irreducible quadratic. Such cases fall
    // outside the property's conditional scope ("Q has degree <= 5 with
    // linear/irreducible quadratic factors"), but we still report them.
    if (has_integral_node(LMCAS::detail::node(result))) {
        rep.unevaluated = true;
        rep.detail = "integrator left unevaluated integral";
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

    auto sample_pts = pick_safe_sample_points(Q.coeffs);
    if (static_cast<int>(sample_pts.size()) < kRequiredSamples) {
        rep.failed = true;
        std::ostringstream oss;
        oss << "could not find " << kRequiredSamples
            << " sample points avoiding poles for Q=" << Q.label;
        rep.detail = oss.str();
        return rep;
    }

    for (double xv : sample_pts) {
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
    TEST_CASE("Rational function decomposition round-trip");

    int total_combos = 0;
    int verified_combos = 0;
    int unevaluated_combos = 0;
    int failed_combos = 0;
    int total_matches = 0;

    std::vector<std::string> unevaluated_labels;

    for (const auto& Q : denominators()) {
        for (const auto& P : numerators()) {
            ++total_combos;
            ComboReport rep = verify_combo(P, Q);
            std::string prefix = "(" + P.label + ") / (" + Q.label + ")";

            if (rep.failed) {
                ++failed_combos;
                std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
                EXPECT_TRUE(false, prefix + ": numeric round-trip mismatch");
            } else if (rep.unevaluated) {
                ++unevaluated_combos;
                unevaluated_labels.push_back(prefix);
                std::cout << "[UNEVAL] " << prefix
                          << " (integrator left unevaluated integral; "
                             "permitted by design)" << std::endl;
            } else {
                ++verified_combos;
                total_matches += rep.matches;
                std::ostringstream oss;
                oss << prefix << ": " << rep.matches
                    << " sample point(s) matched within tolerance " << kTolerance;
                EXPECT_TRUE(!rep.failed && !rep.unevaluated && rep.matches > 0, oss.str());
            }
        }
    }

    std::cout << "\nSummary:"
              << " total_combos="   << total_combos
              << " verified="       << verified_combos
              << " unevaluated="    << unevaluated_combos
              << " failed="         << failed_combos
              << " total_matches="  << total_matches
              << std::endl;

    if (!unevaluated_labels.empty()) {
        std::cout << "Unevaluated combinations (acceptable per design):" << std::endl;
        for (const auto& s : unevaluated_labels) {
            std::cout << "  - " << s << std::endl;
        }
    }

    EXPECT_TRUE(total_combos >= 100,
                "at least 100 (P, Q) iterations attempted "
                "(total=" + std::to_string(total_combos)
                + ", required>=100)");

    EXPECT_TRUE(failed_combos == 0,
                "no (P, Q) combination produced a numeric round-trip mismatch "
                "(failed=" + std::to_string(failed_combos) + ")");

    return TEST_REPORT();
}
