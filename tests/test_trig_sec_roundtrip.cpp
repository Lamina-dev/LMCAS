// Feature: integration-enhancements, Property 6: Trigonometric sec^n round-trip
//
// Validates: Requirements 3.8
//
// Property 6: For all even integers n with 2 <= n <= 8, integrating sec^n(x)
// and differentiating the result SHALL yield an expression numerically equal
// to sec^n(x) at sample points x in {0.3, 0.5, 0.7, 0.9, 1.1} within
// tolerance 1e-10.
//
// Approach
// --------
//   * For every even integer n in {2, 4, 6, 8}, build the integrand
//     sec(x)^n by direct AST construction (FunctionNode of FuncType::Sec
//     wrapped in a PowerNode).
//   * Run the integrator and verify the result is closed form (no leftover
//     unevaluated integral nodes).
//   * Symbolically differentiate the result via SymbolicExpr::differentiate.
//   * For each sample point x in {0.3, 0.5, 0.7, 0.9, 1.1}, substitute,
//     simplify, and numerically evaluate both the integrand sec^n(x) and
//     the derivative of the antiderivative. test_numeric_eval has been
//     extended to evaluate FunctionNode::FuncType::Sec via 1/cos(x).
//   * Sample points where either side returns std::nullopt or non-finite
//     values (e.g. due to cos(x) crossing zero) are skipped.
//   * The test FAILS if any (n, x) pair produces a numeric mismatch above
//     tolerance, or if any n in {2, 4, 6, 8} cannot be verified at any
//     sample point.

#include "test_common.hpp"
#include "integration.hpp"
#include "symbolic_ast.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>

using lamina::Integrator;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-10;

// ----- AST helpers --------------------------------------------------------

std::shared_ptr<SymbolicExpr> sec_of(std::shared_ptr<SymbolicExpr> arg) {
    using FT = FunctionNode::FuncType;
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            FT::Sec,
            std::vector<std::shared_ptr<SymbolicNode>>{arg->root}));
}

std::shared_ptr<SymbolicExpr> sec_pow(std::shared_ptr<SymbolicExpr> arg, int n) {
    auto s = sec_of(arg);
    if (n == 1) return s;
    return SymbolicExpr::power(s, SymbolicExpr::number(n));
}

bool has_integral_node(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Calculus_Integral) return true;
        for (auto& a : fn->arguments)
            if (has_integral_node(a)) return true;
    } else if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands)
            if (has_integral_node(op)) return true;
    } else if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands)
            if (has_integral_node(op)) return true;
    } else if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (has_integral_node(pow->base)) return true;
        if (has_integral_node(pow->exponent)) return true;
    }
    return false;
}

const std::vector<double>& sample_points() {
    static const std::vector<double> S = {0.3, 0.5, 0.7, 0.9, 1.1};
    return S;
}

// ----- Per-n verification -------------------------------------------------

struct NReport {
    bool unevaluated = false;
    bool failed = false;
    int matches = 0;
    int skipped = 0;
    std::string detail;
};

NReport verify_n(int n) {
    NReport rep;

    auto x_var = SymbolicExpr::variable(kVarName);
    auto integrand = sec_pow(x_var, n);

    Integrator integ;
    SymbolicExpr result;
    try {
        result = integ.integrate(*integrand, kVarName);
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception during integration: ") + e.what();
        return rep;
    }

    if (has_integral_node(result.root)) {
        rep.unevaluated = true;
        rep.detail = "unevaluated integral in result: " + result.to_string();
        return rep;
    }

    auto deriv = result.differentiate(kVarName);
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
                << " | result=" << result.to_string();
            rep.detail = oss.str();
            break;
        }
    }
    return rep;
}

}// anonymous namespace

int main() {
    TEST_CASE("Property 6: Trigonometric sec^n round-trip");

    const std::vector<int> ns = {2, 4, 6, 8};

    int total_n = 0;
    int verified_n = 0;
    int unevaluated_n = 0;
    int failed_n = 0;
    int total_matches = 0;
    int total_skipped = 0;

    for (int n : ns) {
        ++total_n;
        NReport rep = verify_n(n);
        std::string prefix = "sec(x)^" + std::to_string(n);

        if (rep.failed) {
            ++failed_n;
            std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
            EXPECT_TRUE(false, prefix + ": numeric round-trip mismatch");
        } else if (rep.unevaluated) {
            ++unevaluated_n;
            std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
            EXPECT_TRUE(false, prefix + ": integrator left unevaluated integral");
        } else if (rep.matches > 0) {
            ++verified_n;
            total_matches += rep.matches;
            total_skipped += rep.skipped;
            std::ostringstream oss;
            oss << prefix << ": " << rep.matches << " match(es), "
                << rep.skipped << " skipped point(s)";
            EXPECT_TRUE(true, oss.str());
        } else {
            // No sample point produced an evaluable comparison. For sec^n
            // with n in {2,4,6,8} on the chosen sample points, sec is well
            // defined everywhere, so this would indicate a numeric-eval
            // gap rather than a domain issue. Treat it as a failure to
            // surface the problem.
            ++failed_n;
            std::cerr << "[FAIL] " << prefix
                      << " : no sample point produced an evaluable comparison" << std::endl;
            EXPECT_TRUE(false, prefix + ": no evaluable sample point");
        }
    }

    std::cout << "\nSummary:"
              << " total_n="     << total_n
              << " verified_n="  << verified_n
              << " unevaluated=" << unevaluated_n
              << " failed_n="    << failed_n
              << " total_matches=" << total_matches
              << " total_skipped=" << total_skipped
              << std::endl;

    EXPECT_TRUE(verified_n == 4,
                "all 4 even values of n in {2,4,6,8} verified by round-trip "
                "(verified=" + std::to_string(verified_n) + ", required=4)");

    return TEST_REPORT();
}
