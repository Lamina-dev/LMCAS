
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

const std::vector<double>& sample_points() {
    static const std::vector<double> S = {0.5, 1.0, 1.5, 2.0, 2.5};
    return S;
}


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
    auto integrated = integ.integrate(*integrand, kVarName);
    if (!integrated) {
        rep.failed = true;
        rep.detail = std::string("integration failed: ") + integrated.error().message;
        return rep;
    }
    auto result = lamina::detail::make_expression_ptr(integrated.value());

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
    TEST_CASE("Trigonometric sin^m * cos^n round-trip");

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
