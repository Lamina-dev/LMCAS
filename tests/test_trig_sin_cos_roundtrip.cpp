
#include "test_common.hpp"
#include "integration.hpp"
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
constexpr double kTolerance = 1e-10;

std::shared_ptr<SymbolicExpr> num_int(const BigInt &n) {
    return SymbolicExpr::number(n);
}

// Build f(x)^k, but if k == 1 return f(x) unwrapped so the integrand stays
// in the exact shape the TrigCombinationStrategy detects.
std::shared_ptr<SymbolicExpr> pow_or_self(std::shared_ptr<SymbolicExpr> base, const BigInt &k) {
    if (k == 1) {
        return base;
    }
    return SymbolicExpr::power(base, num_int(k));
}

// Build the integrand sin(x)^m * cos(x)^n in the simplest shape.
std::shared_ptr<SymbolicExpr> build_integrand(const BigInt &m, const BigInt &n) {
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

const std::vector<double> &sample_points() {
    static const std::vector<double> S = {0.5, 1.0, 1.5, 2.0, 2.5};
    return S;
}

struct PairReport {
    bool failed = false;
    std::size_t matches = 0;
    std::string detail;
};

void compare_roundtrip_samples(
    const std::shared_ptr<SymbolicExpr> &integrand_simp,
    const std::shared_ptr<SymbolicExpr> &deriv_simp,
    const std::shared_ptr<SymbolicExpr> &result, PairReport &rep) {
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
            return;
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

PairReport verify_pair(const BigInt &m, const BigInt &n) {
    PairReport rep;

    auto integrand = build_integrand(m, n);
    if (!integrand) {
        rep.failed = true;
        rep.detail = "integrand build returned null";
        return rep;
    }

    Integrator integ;
    auto integrated = integ.integrate(*integrand, kVarName);
    if (!integrated) {
        rep.failed = true;
        rep.detail = std::string("integration failed: ") + integrated.error().message;
        return rep;
    }
    auto result = LMCAS::detail::make_expression_ptr(integrated.value());

    if (LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(result))) {
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
    if (!deriv_simp) {
        deriv_simp = deriv;
    }

    auto integrand_simp = integrand->simplify();
    if (!integrand_simp) {
        integrand_simp = integrand;
    }

    compare_roundtrip_samples(integrand_simp, deriv_simp, result, rep);
    return rep;
}

} // anonymous namespace

TEST(TrigSinCosRoundtrip, AllPowersWithTotalDegreeAtMostEight) {
    for (std::size_t m = 0; m <= 8; ++m) {
        for (std::size_t n = 0; n + m <= 8; ++n) {
            SCOPED_TRACE("sin^" + std::to_string(m) + "(x) * cos^" + std::to_string(n) + "(x)");
            const PairReport rep = verify_pair(BigInt(m), BigInt(n));
            EXPECT_FALSE(rep.failed) << rep.detail;
            EXPECT_GT(rep.matches, 0u);
        }
    }
}
