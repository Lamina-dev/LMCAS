#include "test_common.hpp"
#include "integration.hpp"
#include "rational.hpp"
#include "bigint.hpp"

#include <cstddef>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <functional>
#include <sstream>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char *kVarName = "x";
constexpr double kTolerance = 1e-10;

std::shared_ptr<SymbolicExpr> num_rat(const BigInt &n, const BigInt &d) {
    return SymbolicExpr::number(Rational(n, d));
}

std::shared_ptr<SymbolicExpr> num_int(const BigInt &n) {
    return SymbolicExpr::number(n);
}

// A pattern is a builder that takes a sub-expression `arg` and produces
// f(arg) as a fresh SymbolicExpr. We deliberately restrict ourselves to
// patterns whose antiderivative-and-derivative chain reduces to function
// nodes that test_numeric_eval supports.
struct Pattern {
    std::string name;
    std::function<std::shared_ptr<SymbolicExpr>(std::shared_ptr<SymbolicExpr>)> build;
};

const std::vector<Pattern> &patterns() {
    static const std::vector<Pattern> P = {
        // Trigonometric / exponential: derivative chain stays in {sin, cos, exp}.
        {"sin(arg)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::sin(a); }},
        {"cos(arg)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::cos(a); }},
        {"exp(arg)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::exp(a); }},

        // Polynomial powers via the table's x^n / 1/x rules.
        {"arg^2", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(a, num_int(2)); }},
        {"arg^3", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(a, num_int(3)); }},
        {"1/arg", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(a, num_int(-1)); }},

        // Squared trig: round-trip flows through half-angle identities, but
        // the derivative stays in {sin, cos} and is therefore evaluable.
        {"sin(arg)^2", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(SymbolicExpr::sin(a), num_int(2)); }},
        {"cos(arg)^2", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(SymbolicExpr::cos(a), num_int(2)); }},
    };
    return P;
}

struct ABPair {
    std::shared_ptr<SymbolicExpr> a;
    std::shared_ptr<SymbolicExpr> b;
    std::string label;
};

const std::vector<ABPair> &ab_pairs() {
    static const std::vector<ABPair> pairs = []() {
        // 5 values of a (all non-zero) x 5 values of b = 25 pairs.
        struct Frac {
            BigInt n, d;
        };
        const std::vector<Frac> a_vals = {
            {1, 1}, {2, 1}, {-1, 1}, {3, 1}, {1, 2}};
        const std::vector<Frac> b_vals = {
            {0, 1}, {1, 1}, {-1, 1}, {2, 1}, {1, 2}};
        std::vector<ABPair> out;
        out.reserve(a_vals.size() * b_vals.size());
        for (const auto &fa : a_vals) {
            for (const auto &fb : b_vals) {
                ABPair p;
                p.a = num_rat(fa.n, fa.d);
                p.b = num_rat(fb.n, fb.d);
                std::ostringstream oss;
                oss << "a=" << fa.n.to_string() << "/" << fa.d.to_string()
                    << ",b=" << fb.n.to_string() << "/" << fb.d.to_string();
                p.label = oss.str();
                out.push_back(std::move(p));
            }
        }
        return out;
    }();
    return pairs;
}

const std::vector<double> &sample_points() {
    static const std::vector<double> S = {0.5, 1.0, 1.5, 2.0, 2.5};
    return S;
}

struct ComboReport {
    bool unevaluated = false; // integrator returned an unevaluated integral
    bool failed = false;      // numeric mismatch above tolerance
    std::size_t matches = 0;
    std::size_t skipped = 0;
    std::string detail; // failure / status message
};

static void compare_roundtrip_samples(ComboReport &rep,
                                      const std::shared_ptr<SymbolicExpr> &integrand_simp,
                                      const std::shared_ptr<SymbolicExpr> &deriv_simp,
                                      const std::shared_ptr<SymbolicExpr> &result) {
    for (double xv : sample_points()) {
        auto x_val = SymbolicExpr::number(xv);
        auto integrand_at = integrand_simp->substitute(kVarName, x_val);
        auto deriv_at = deriv_simp->substitute(kVarName, x_val);
        if (!integrand_at || !deriv_at) {
            ++rep.skipped;
            continue;
        }
        integrand_at = integrand_at->simplify();
        deriv_at = deriv_at->simplify();

        auto pv = test_numeric_eval(integrand_at);
        auto dv = test_numeric_eval(deriv_at);
        if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
            ++rep.skipped;
            continue;
        }
        double delta = std::abs(*pv - *dv);
        if (delta <= kTolerance) {
            ++rep.matches;
        } else {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv
                << ": integrand=" << *pv
                << " vs d/dx(result)=" << *dv
                << " |delta|=" << delta
                << " | result=" << result->to_string();
            rep.detail = oss.str();
            break;
        }
    }
}

ComboReport verify_combo(const Pattern &pat, const ABPair &ab) {
    ComboReport rep;

    // Build integrand: f(a*x + b)
    auto x_var = SymbolicExpr::variable(kVarName);
    auto ax = SymbolicExpr::multiply(ab.a, x_var);
    auto arg = SymbolicExpr::add(ax, ab.b);
    auto integrand_ptr = pat.build(arg);
    if (!integrand_ptr) {
        rep.detail = "integrand build returned null";
        rep.failed = true;
        return rep;
    }

    // Integrate.
    Integrator integ;
    auto integrated = integ.integrate(*integrand_ptr, kVarName);
    if (!integrated) {
        rep.failed = true;
        rep.detail = std::string("integration failed: ") + integrated.error().message;
        return rep;
    }
    auto result = LMCAS::detail::make_expression_ptr(integrated.value());

    // Combinations that left an unevaluated integral are not evidence for
    // the round-trip property; they fall outside its conditional scope.
    if (LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(result))) {
        rep.unevaluated = true;
        rep.detail = "unevaluated integral in result";
        return rep;
    }

    auto deriv = result->differentiate(kVarName);
    if (!deriv) {
        rep.failed = true;
        rep.detail = "differentiation returned null";
        return rep;
    }
    auto deriv_simp = deriv->simplify();
    if (!deriv_simp) {
        deriv_simp = deriv;
    }

    auto integrand_simp = integrand_ptr->simplify();
    if (!integrand_simp) {
        integrand_simp = integrand_ptr;
    }

    compare_roundtrip_samples(rep, integrand_simp, deriv_simp, result);
    return rep;
}

} // anonymous namespace

TEST(LinearSubstitutionRoundtrip, AffineArgumentsDifferentiateToTheirIntegrands) {
    std::size_t verified_combos = 0;
    for (const auto &pat : patterns()) {
        for (const auto &ab : ab_pairs()) {
            SCOPED_TRACE(pat.name + " [" + ab.label + "]");
            const ComboReport rep = verify_combo(pat, ab);
            ASSERT_FALSE(rep.failed) << rep.detail;
            ASSERT_FALSE(rep.unevaluated) << rep.detail;
            ASSERT_GT(rep.matches, 0u) << rep.detail;
            ++verified_combos;
        }
    }
    EXPECT_EQ(verified_combos, patterns().size() * ab_pairs().size());
}
