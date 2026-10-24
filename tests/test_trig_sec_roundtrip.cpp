#include "expr.hpp"
#include "test_common.hpp"
#include "integration.hpp"
#include "residual_verification.hpp"
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
        auto integrand_at = LMCAS::substitute(integrand_simp, kVarName, x_val).value();
        auto deriv_at = LMCAS::substitute(deriv_simp, kVarName, x_val).value();
        if (!integrand_at || !deriv_at) {
            ++rep.skipped;
            continue;
        }
        integrand_at = LMCAS::simplify(integrand_at).value();
        deriv_at = LMCAS::simplify(deriv_at).value();

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
        rep.detail = "unevaluated integral in result: " + result->to_string();
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
