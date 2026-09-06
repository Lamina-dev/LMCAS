
#include "test_common.hpp"
#include "integration.hpp"
#include "rational.hpp"
#include "bigint.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <functional>
#include <sstream>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-10;


std::shared_ptr<SymbolicExpr> num_rat(long long n, long long d) {
    return SymbolicExpr::number(Rational(BigInt(n), BigInt(d)));
}

std::shared_ptr<SymbolicExpr> num_int(long long n) {
    return SymbolicExpr::number(n);
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


// A pattern is a builder that takes a sub-expression `arg` and produces
// f(arg) as a fresh SymbolicExpr. We deliberately restrict ourselves to
// patterns whose antiderivative-and-derivative chain reduces to function
// nodes that test_numeric_eval supports.
struct Pattern {
    std::string name;
    std::function<std::shared_ptr<SymbolicExpr>(std::shared_ptr<SymbolicExpr>)> build;
};

const std::vector<Pattern>& patterns() {
    static const std::vector<Pattern> P = {
        // Trigonometric / exponential: derivative chain stays in {sin, cos, exp}.
        {"sin(arg)",     [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::sin(a); }},
        {"cos(arg)",     [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::cos(a); }},
        {"exp(arg)",     [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::exp(a); }},

        // Polynomial powers via the table's x^n / 1/x rules.
        {"arg^2",        [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(a, num_int(2)); }},
        {"arg^3",        [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(a, num_int(3)); }},
        {"1/arg",        [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(a, num_int(-1)); }},

        // Squared trig: round-trip flows through half-angle identities, but
        // the derivative stays in {sin, cos} and is therefore evaluable.
        {"sin(arg)^2",   [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(SymbolicExpr::sin(a), num_int(2)); }},
        {"cos(arg)^2",   [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(SymbolicExpr::cos(a), num_int(2)); }},
    };
    return P;
}


struct ABPair {
    std::shared_ptr<SymbolicExpr> a;
    std::shared_ptr<SymbolicExpr> b;
    std::string label;
};

const std::vector<ABPair>& ab_pairs() {
    static const std::vector<ABPair> pairs = []() {
        // 5 values of a (all non-zero) x 5 values of b = 25 pairs.
        struct Frac { long long n, d; };
        const std::vector<Frac> a_vals = {
            {1, 1}, {2, 1}, {-1, 1}, {3, 1}, {1, 2}
        };
        const std::vector<Frac> b_vals = {
            {0, 1}, {1, 1}, {-1, 1}, {2, 1}, {1, 2}
        };
        std::vector<ABPair> out;
        out.reserve(a_vals.size() * b_vals.size());
        for (const auto& fa : a_vals) {
            for (const auto& fb : b_vals) {
                ABPair p;
                p.a = num_rat(fa.n, fa.d);
                p.b = num_rat(fb.n, fb.d);
                std::ostringstream oss;
                oss << "a=" << fa.n << "/" << fa.d
                    << ",b=" << fb.n << "/" << fb.d;
                p.label = oss.str();
                out.push_back(std::move(p));
            }
        }
        return out;
    }();
    return pairs;
}

const std::vector<double>& sample_points() {
    static const std::vector<double> S = {0.5, 1.0, 1.5, 2.0, 2.5};
    return S;
}


struct ComboReport {
    bool unevaluated = false;   // integrator returned an unevaluated integral
    bool failed = false;        // numeric mismatch above tolerance
    int  matches = 0;
    int  skipped = 0;
    std::string detail;         // failure / status message
};

ComboReport verify_combo(const Pattern& pat, const ABPair& ab) {
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
    if (has_integral_node(LMCAS::detail::node(result))) {
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
    if (!deriv_simp) deriv_simp = deriv;

    auto integrand_simp = integrand_ptr->simplify();
    if (!integrand_simp) integrand_simp = integrand_ptr;

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
    return rep;
}

}// anonymous namespace

int main() {
    TEST_CASE("Linear substitution round-trip");

    int total_combos = 0;
    int verified_combos = 0;
    int unevaluated_combos = 0;
    int skipped_combos = 0;
    int failed_combos = 0;
    int total_matches = 0;
    int total_skipped_pts = 0;

    for (const auto& pat : patterns()) {
        for (const auto& ab : ab_pairs()) {
            ++total_combos;

            ComboReport rep = verify_combo(pat, ab);
            std::string prefix = pat.name + " [" + ab.label + "]";

            if (rep.failed) {
                ++failed_combos;
                std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
                EXPECT_TRUE(false, prefix + ": numeric round-trip mismatch");
            } else if (rep.unevaluated) {
                ++unevaluated_combos;
                std::cout << "[UNEVAL] " << prefix
                          << " (integrator left unevaluated integral; outside property scope)"
                          << std::endl;
            } else if (rep.matches > 0) {
                ++verified_combos;
                total_matches += rep.matches;
                total_skipped_pts += rep.skipped;
            } else {
                ++skipped_combos;
                total_skipped_pts += rep.skipped;
                std::cout << "[SKIP] " << prefix
                          << " (no sample point produced a numerically evaluable comparison)"
                          << std::endl;
            }
        }
    }

    std::cout << "\nSummary:"
              << " total_combos="     << total_combos
              << " verified_combos="  << verified_combos
              << " unevaluated="      << unevaluated_combos
              << " skipped_combos="   << skipped_combos
              << " failed_combos="    << failed_combos
              << " total_matches="    << total_matches
              << " skipped_points="   << total_skipped_pts
              << std::endl;

    EXPECT_TRUE(verified_combos >= 100,
                "at least 100 (table entry, a, b) combinations verified by round-trip "
                "(verified=" + std::to_string(verified_combos)
                + ", required>=100)");

    return TEST_REPORT();
}
