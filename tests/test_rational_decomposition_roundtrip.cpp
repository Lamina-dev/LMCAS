#include "expr.hpp"
#include "test_common.hpp"
#include "integration.hpp"
#include "residual_verification.hpp"
#include "rational.hpp"
#include "bigint.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char *kVarName = "x";
constexpr double kTolerance = 1e-8;
constexpr double kPoleGuard = 0.1;
constexpr std::size_t kRequiredSamples = 5;

std::shared_ptr<SymbolicExpr> num_int(const BigInt &n) {
    return SymbolicExpr::number(n);
}

// Build a symbolic polynomial from coefficients (low-to-high order) using
// only AddNode/MultiplyNode/PowerNode/NumberNode/VariableNode so that
// test_numeric_eval can evaluate it directly.
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const std::vector<BigInt> &coeffs) {
    auto x_var = SymbolicExpr::variable(kVarName);
    std::shared_ptr<SymbolicExpr> result = num_int(0);
    bool is_zero = true;
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        if (coeffs[k] == 0) {
            continue;
        }
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
            auto pw = SymbolicExpr::power(x_var, num_int(BigInt(k)));
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
double poly_eval(const std::vector<BigInt> &coeffs, double x) {
    double y = 0.0;
    for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
        y = y * x + it->to_double();
    }
    return y;
}

const std::vector<double> &candidate_sample_points() {
    static const std::vector<double> S = {
        -3.7, -2.7, -1.7, -0.7,
        0.3, 0.7, 1.3, 1.7,
        2.3, 2.7, 3.3, 3.7,
        4.3, 4.7};
    return S;
}

struct Poly {
    std::vector<BigInt> coeffs; /**< 按次数升序排列的系数。 */
    std::string label;
};

// Denominators of degree 3..5, all factorable over Q into linear and
// irreducible quadratic factors. Specifically excluded: repeated
// irreducible quadratic factors (e.g. (x^2+1)^2), which the design
// document explicitly leaves as unevaluated integrals.
const std::vector<Poly> &denominators() {
    static const std::vector<Poly> D = {
        {{1, 0, 2}, "2x^2 + 1"},
        {{-2, 0, 1}, "x^2 - 2"},
        {{1, 2, 1}, "(x + 1)^2"},
        {{-2, 0, 3}, "3x^2 - 2"},
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
const std::vector<Poly> &numerators() {
    static const std::vector<Poly> N = {
        {{1}, "1"},
        {{0, 1}, "x"},
        {{1, 1}, "x + 1"},
        {{-1, 1}, "x - 1"},
        {{3, 2}, "2x + 3"},
        {{0, 0, 1}, "x^2"},
        {{1, 0, 1}, "x^2 + 1"},
        {{-1, 0, 1}, "x^2 - 1"},
        {{5, 1, 3}, "3x^2 + x + 5"},
        {{-3, 0, 2}, "2x^2 - 3"},
    };
    return N;
}

std::vector<double> pick_safe_sample_points(const std::vector<BigInt> &Q_coeffs) {
    std::vector<double> kept;
    for (double x : candidate_sample_points()) {
        double qv = poly_eval(Q_coeffs, x);
        if (std::abs(qv) >= kPoleGuard) {
            kept.push_back(x);
            if (kept.size() >= kRequiredSamples) {
                break;
            }
        }
    }
    return kept;
}

struct ComboReport {
    bool unevaluated = false;
    bool failed = false;
    std::size_t matches = 0;
    std::string detail;
};

static void compare_roundtrip_samples(ComboReport &rep,
                                      const std::shared_ptr<SymbolicExpr> &integrand_simp,
                                      const std::shared_ptr<SymbolicExpr> &deriv_simp,
                                      const std::shared_ptr<SymbolicExpr> &result,
                                      const std::vector<double> &sample_pts) {
    for (double xv : sample_pts) {
        auto x_val = SymbolicExpr::number(xv);
        auto integrand_at = LMCAS::substitute(integrand_simp, kVarName, x_val).value();
        auto deriv_at = LMCAS::substitute(deriv_simp, kVarName, x_val).value();
        if (!integrand_at || !deriv_at) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv << ": substitute returned null"
                << " | result=" << result->to_string();
            rep.detail = oss.str();
            return;
        }
        integrand_at = LMCAS::simplify(integrand_at).value();
        deriv_at = LMCAS::simplify(deriv_at).value();

        auto pv = test_numeric_eval(integrand_at);
        auto dv = test_numeric_eval(deriv_at);
        if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv << ": numeric evaluation failed"
                << " | result=" << result->to_string()
                << " | derivative=" << deriv_simp->to_string();
            rep.detail = oss.str();
            return;
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
            return;
        }
        ++rep.matches;
    }
}

ComboReport verify_combo(const Poly &P, const Poly &Q) {
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
    ComputationContext context;
    auto integrated = integ.integrate_checked(*integrand, kVarName, context);
    if (!integrated) {
        rep.failed = true;
        rep.detail = std::string("integration failed: ") + integrated.error().message;
        return rep;
    }
    auto result = LMCAS::detail::make_expression_ptr(integrated.value());

    if (LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(result))) {
        rep.unevaluated = true;
        rep.detail = "integrator left unevaluated integral";
        return rep;
    }

    auto deriv = LMCAS::differentiate(result, kVarName).value();
    if (!deriv) {
        rep.failed = true;
        rep.detail = "differentiation returned null";
        return rep;
    }
    auto proof = check_equivalent(deriv, integrand, context);
    if (!proof || !std::holds_alternative<ProvedZeroResidual>(proof.value())) {
        rep.failed = true;
        rep.detail = proof ? "derivative residual was not proved zero" : proof.error().message;
        return rep;
    }
    auto deriv_simp = LMCAS::simplify(deriv).value();
    if (!deriv_simp) {
        deriv_simp = deriv;
    }

    auto integrand_simp = LMCAS::simplify(integrand).value();
    if (!integrand_simp) {
        integrand_simp = integrand;
    }

    auto sample_pts = pick_safe_sample_points(Q.coeffs);
    if (sample_pts.size() < kRequiredSamples) {
        rep.failed = true;
        std::ostringstream oss;
        oss << "could not find " << kRequiredSamples
            << " sample points avoiding poles for Q=" << Q.label;
        rep.detail = oss.str();
        return rep;
    }

    compare_roundtrip_samples(rep, integrand_simp, deriv_simp, result, sample_pts);
    return rep;
}

} // anonymous namespace

TEST(RationalDecompositionRoundtrip, PolynomialQuotientsDifferentiateToTheirIntegrands) {
    std::size_t verified_combos = 0;
    for (const auto &Q : denominators()) {
        for (const auto &P : numerators()) {
            SCOPED_TRACE("(" + P.label + ") / (" + Q.label + ")");
            const ComboReport rep = verify_combo(P, Q);
            ASSERT_FALSE(rep.failed) << rep.detail;
            ASSERT_FALSE(rep.unevaluated) << rep.detail;
            EXPECT_EQ(rep.matches, kRequiredSamples) << rep.detail;
            ++verified_combos;
        }
    }
    EXPECT_EQ(verified_combos, denominators().size() * numerators().size());
}
