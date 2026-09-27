
#include "test_common.hpp"
#include "integration.hpp"
#include "internal/symbolic_ast.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char *kVarName = "x";
constexpr double kTolerance = 1e-10;

std::shared_ptr<SymbolicExpr> sec_of(std::shared_ptr<SymbolicExpr> arg) {
    using FT = FunctionNode::FuncType;
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FT::Sec,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(arg)}));
}

std::shared_ptr<SymbolicExpr> sec_pow(std::shared_ptr<SymbolicExpr> arg, int n) {
    auto s = sec_of(arg);
    if (n == 1) {
        return s;
    }
    return SymbolicExpr::power(s, SymbolicExpr::number(n));
}

const std::vector<double> &sample_points() {
    static const std::vector<double> S = {0.3, 0.5, 0.7, 0.9, 1.1};
    return S;
}

struct NReport {
    bool unevaluated = false;
    bool failed = false;
    int matches = 0;
    int skipped = 0;
    std::string detail;
};

void compare_roundtrip_samples(
    const std::shared_ptr<SymbolicExpr> &integrand_simp,
    const std::shared_ptr<SymbolicExpr> &deriv_simp,
    const std::shared_ptr<SymbolicExpr> &result, NReport &rep) {
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

NReport verify_n(int n) {
    NReport rep;

    auto x_var = SymbolicExpr::variable(kVarName);
    auto integrand = sec_pow(x_var, n);

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
        rep.unevaluated = true;
        rep.detail = "unevaluated integral in result: " + result->to_string();
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

TEST(TrigSecRoundtrip, EvenPowersThroughEight) {
    for (int n : {2, 4, 6, 8}) {
        SCOPED_TRACE("sec(x)^" + std::to_string(n));
        const NReport rep = verify_n(n);
        EXPECT_FALSE(rep.failed) << rep.detail;
        EXPECT_FALSE(rep.unevaluated) << rep.detail;
        EXPECT_GT(rep.matches, 0) << "no evaluable sample point";
    }
}
