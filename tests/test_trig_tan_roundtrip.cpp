
#include "test_common.hpp"
#include "integration.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <sstream>

using lamina::Integrator;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-10;


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
    static const std::vector<double> S = {0.3, 0.5, 0.7, 0.9, 1.1};
    return S;
}


struct NReport {
    bool unevaluated = false;   // integrator returned an unevaluated integral
    bool failed = false;        // numeric mismatch above tolerance
    int  matches = 0;
    int  skipped = 0;
    std::string detail;
};

NReport verify_n(int n) {
    NReport rep;

    // Build integrand: tan(x)^n
    auto x_var = SymbolicExpr::variable(kVarName);
    auto tan_x = SymbolicExpr::tan(x_var);
    auto integrand_ptr =
        SymbolicExpr::power(tan_x, SymbolicExpr::number(n));
    if (!integrand_ptr) {
        rep.failed = true;
        rep.detail = "integrand build returned null";
        return rep;
    }

    // Integrate.
    Integrator integ;
    std::shared_ptr<SymbolicExpr> result;
    try {
        result = lamina::detail::make_expression_ptr(integ.integrate(*integrand_ptr, kVarName));
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception during integration: ") + e.what();
        return rep;
    }

    if (has_integral_node(lamina::detail::node(result))) {
        rep.unevaluated = true;
        rep.detail = "integrator left an unevaluated integral: "
                     + result->to_string();
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
    TEST_CASE("Trigonometric tan^n round-trip");

    int total = 0;
    int verified = 0;
    int unevaluated = 0;
    int failed = 0;

    for (int n = 2; n <= 8; ++n) {
        ++total;
        NReport rep = verify_n(n);
        std::string prefix = "tan(x)^" + std::to_string(n);

        if (rep.failed) {
            ++failed;
            std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
            EXPECT_TRUE(false, prefix + ": numeric round-trip mismatch");
        } else if (rep.unevaluated) {
            ++unevaluated;
            std::cerr << "[FAIL] " << prefix << " : " << rep.detail << std::endl;
            EXPECT_TRUE(false, prefix + ": integrator left unevaluated integral "
                                        "(requires success on n in [2,8])");
        } else if (rep.matches > 0) {
            ++verified;
            std::ostringstream oss;
            oss << prefix << ": " << rep.matches << " match(es), "
                << rep.skipped << " skipped point(s)";
            EXPECT_TRUE(!rep.failed && !rep.unevaluated && rep.matches > 0, oss.str());
        } else {
            ++failed;
            std::cerr << "[FAIL] " << prefix
                      << " : no sample point produced a numerically evaluable "
                         "comparison (skipped=" << rep.skipped << ")" << std::endl;
            EXPECT_TRUE(false, prefix + ": no numerically verified sample");
        }
    }

    std::cout << "\nSummary:"
              << " total="       << total
              << " verified="    << verified
              << " unevaluated=" << unevaluated
              << " failed="      << failed
              << std::endl;

    EXPECT_TRUE(verified == 7,
                "all 7 values of n in [2, 8] verified by round-trip "
                "(verified=" + std::to_string(verified) + ", required=7)");

    return TEST_REPORT();
}
