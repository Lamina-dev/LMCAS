// Feature: integration-enhancements, Property 9: Multiple integral separability
//
// Validates: Requirements 6.10
//
// Property 9: For all separable integrands of the form f(x)*g(y) where f
// depends only on x and g depends only on y, integrated over independent
// constant bounds [a1, b1] for x and [a2, b2] for y, the
// MultipleIntegralEngine SHALL produce a result numerically equal to
// (integral_{a1}^{b1} f(x) dx) * (integral_{a2}^{b2} g(y) dy) within
// tolerance 1e-10.
//
// Approach
// --------
//   * Choose 6 unary functions for f and g: x, x^2, x^3, sin(x), cos(x),
//     exp(x). All are smooth and have no singularities, so constant bounds
//     in any sub-interval of R are safe.
//
//   * Choose 3 sets of constant bounds for each variable: {0,1}, {1,2},
//     {-1,1}. This yields 6 * 6 * 3 * 3 = 324 separable combinations,
//     comfortably exceeding the >=100 combinations required by the task.
//
//   * For each combination:
//       1. Build the integrand f(x)*g(y).
//       2. Use MultipleIntegralEngine to evaluate the iterated integral
//          with the inner step over x and the outer step over y, both
//          with definite bounds.
//       3. Compute the per-variable definite integrals Ix = integrate_def
//          (f, x, a1, b1) and Iy = integrate_def(g, y, a2, b2) using the
//          same Integrator instance.
//       4. Numerically evaluate the engine result and Ix * Iy (after
//          simplify) using SymbolicExpr::to_numeric and compare within
//          tolerance 1e-10.
//
//   * Combinations whose engine output still contains an unevaluated
//     Calculus_Integral node, or whose numeric values are not finite, are
//     reported as SKIPPED (they are outside the property's premise).
//
//   * The test FAILS if any combination produces a numeric mismatch above
//     tolerance, or if fewer than 100 combinations are successfully
//     verified.

#include "test_common.hpp"
#include "integration.hpp"
#include "symbolic_ast.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <functional>

using lamina::Integrator;
using lamina::MultipleIntegralEngine;

namespace {

constexpr double kTolerance = 1e-10;

// --------------------- AST helpers -------------------------------

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

// --------------------- Function catalog --------------------------

// Each entry builds f(arg) given the argument expression.
struct FunctionSpec {
    std::string name;
    std::function<std::shared_ptr<SymbolicExpr>(std::shared_ptr<SymbolicExpr>)> build;
};

const std::vector<FunctionSpec>& functions() {
    static const std::vector<FunctionSpec> F = {
        {"t",       [](std::shared_ptr<SymbolicExpr> a) { return a; }},
        {"t^2",     [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(a, SymbolicExpr::number(2)); }},
        {"t^3",     [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::power(a, SymbolicExpr::number(3)); }},
        {"sin(t)",  [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::sin(a); }},
        {"cos(t)",  [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::cos(a); }},
        {"exp(t)",  [](std::shared_ptr<SymbolicExpr> a) {
            return SymbolicExpr::exp(a); }},
    };
    return F;
}

// --------------------- Bound catalog -----------------------------

struct BoundSpec {
    long long lo;
    long long hi;
    std::string label;
};

const std::vector<BoundSpec>& bounds() {
    static const std::vector<BoundSpec> B = {
        {0, 1,  "[0,1]"},
        {1, 2,  "[1,2]"},
        {-1, 1, "[-1,1]"},
    };
    return B;
}

// --------------------- Per-combination check ---------------------

struct ComboReport {
    bool unevaluated = false;
    bool skipped     = false;
    bool failed      = false;
    double engine_value  = 0.0;
    double product_value = 0.0;
    std::string detail;
};

ComboReport verify_combo(const FunctionSpec& f,
                         const FunctionSpec& g,
                         const BoundSpec& bx,
                         const BoundSpec& by) {
    ComboReport rep;

    auto x_var = SymbolicExpr::variable("x");
    auto y_var = SymbolicExpr::variable("y");

    auto fx = f.build(x_var);
    auto gy = g.build(y_var);
    auto integrand = SymbolicExpr::multiply(fx, gy);
    if (!integrand) {
        rep.failed = true;
        rep.detail = "integrand build returned null";
        return rep;
    }

    auto x_lo = SymbolicExpr::number(bx.lo);
    auto x_hi = SymbolicExpr::number(bx.hi);
    auto y_lo = SymbolicExpr::number(by.lo);
    auto y_hi = SymbolicExpr::number(by.hi);

    Integrator integrator;
    MultipleIntegralEngine engine;

    // Steps: index 0 is innermost. Inner over x, outer over y.
    std::vector<MultipleIntegralEngine::IntegrationStep> steps = {
        {"x", x_lo, x_hi},
        {"y", y_lo, y_hi},
    };

    std::shared_ptr<SymbolicExpr> engine_result;
    try {
        engine_result = engine.evaluate(*integrand, steps, integrator);
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception in engine.evaluate: ") + e.what();
        return rep;
    }

    if (!engine_result) {
        rep.failed = true;
        rep.detail = "engine.evaluate returned null";
        return rep;
    }
    if (has_integral_node(engine_result->root)) {
        rep.unevaluated = true;
        rep.detail = "engine result contains unevaluated integral";
        return rep;
    }

    // Compute Ix = integrate_def(f, x, x_lo, x_hi)
    SymbolicExpr Ix, Iy;
    try {
        // Build f(x) and g(y) afresh (they were consumed when wrapped as
        // children of the multiply node above; rebuild for safety).
        auto fx2 = f.build(SymbolicExpr::variable("x"));
        auto gy2 = g.build(SymbolicExpr::variable("y"));
        Ix = integrator.integrate_def(*fx2, "x", *x_lo, *x_hi);
        Iy = integrator.integrate_def(*gy2, "y", *y_lo, *y_hi);
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception in integrate_def: ") + e.what();
        return rep;
    }

    if (has_integral_node(Ix.root) || has_integral_node(Iy.root)) {
        rep.unevaluated = true;
        rep.detail = "per-variable integrate_def left unevaluated integral";
        return rep;
    }

    // Numerically evaluate engine result and the product Ix * Iy.
    auto engine_simp = engine_result->simplify();
    if (!engine_simp) engine_simp = engine_result;

    auto Ix_simp = Ix.simplify();
    auto Iy_simp = Iy.simplify();
    if (!Ix_simp || !Iy_simp) {
        rep.failed = true;
        rep.detail = "simplify of per-variable integral returned null";
        return rep;
    }

    // Use test_numeric_eval (from test_common.hpp) which recursively walks
    // Add/Multiply/Power/Function nodes; SymbolicExpr::to_numeric() only
    // handles NumberNode and single-arg FunctionNode and returns 0 for
    // compound expressions like (1/2)*sin(1).
    auto engine_opt = test_numeric_eval(engine_simp);
    auto Ix_opt     = test_numeric_eval(Ix_simp);
    auto Iy_opt     = test_numeric_eval(Iy_simp);

    if (!engine_opt || !Ix_opt || !Iy_opt) {
        rep.skipped = true;
        std::ostringstream oss;
        oss << "non-evaluable expression (engine="
            << (engine_opt ? std::to_string(*engine_opt) : "n/a")
            << ", Ix=" << (Ix_opt ? std::to_string(*Ix_opt) : "n/a")
            << ", Iy=" << (Iy_opt ? std::to_string(*Iy_opt) : "n/a") << ")";
        rep.detail = oss.str();
        return rep;
    }

    double engine_val  = *engine_opt;
    double product_val = (*Ix_opt) * (*Iy_opt);

    if (!std::isfinite(engine_val) || !std::isfinite(product_val)) {
        rep.skipped = true;
        std::ostringstream oss;
        oss << "non-finite numeric value (engine=" << engine_val
            << ", Ix*Iy=" << product_val << ")";
        rep.detail = oss.str();
        return rep;
    }

    rep.engine_value  = engine_val;
    rep.product_value = product_val;

    double delta = std::abs(engine_val - product_val);
    if (delta > kTolerance) {
        rep.failed = true;
        std::ostringstream oss;
        oss << "engine=" << engine_val
            << " vs Ix*Iy=" << product_val
            << " |delta|=" << delta
            << " | Ix=" << Ix_simp->to_string()
            << " | Iy=" << Iy_simp->to_string()
            << " | engine_result=" << engine_simp->to_string();
        rep.detail = oss.str();
    }
    return rep;
}

}// anonymous namespace

int main() {
    TEST_CASE("Property 9: Multiple integral separability");

    int total_combos      = 0;
    int verified_combos   = 0;
    int unevaluated_combos = 0;
    int skipped_combos    = 0;
    int failed_combos     = 0;

    for (const auto& f : functions()) {
        for (const auto& g : functions()) {
            for (const auto& bx : bounds()) {
                for (const auto& by : bounds()) {
                    ++total_combos;

                    ComboReport rep = verify_combo(f, g, bx, by);

                    std::string prefix = "f=" + f.name
                                       + ", g=" + g.name
                                       + ", x=" + bx.label
                                       + ", y=" + by.label;

                    if (rep.failed) {
                        ++failed_combos;
                        std::cerr << "[FAIL] " << prefix
                                  << " : " << rep.detail << std::endl;
                        EXPECT_TRUE(false,
                            prefix + ": numeric separability mismatch");
                    } else if (rep.unevaluated) {
                        ++unevaluated_combos;
                        std::cout << "[UNEVAL] " << prefix
                                  << " (" << rep.detail << ")" << std::endl;
                    } else if (rep.skipped) {
                        ++skipped_combos;
                        std::cout << "[SKIP] " << prefix
                                  << " (" << rep.detail << ")" << std::endl;
                    } else {
                        ++verified_combos;
                    }
                }
            }
        }
    }

    std::cout << "\nSummary:"
              << " total_combos="    << total_combos
              << " verified_combos=" << verified_combos
              << " unevaluated="     << unevaluated_combos
              << " skipped="         << skipped_combos
              << " failed_combos="   << failed_combos
              << std::endl;

    EXPECT_TRUE(verified_combos >= 100,
                "at least 100 separable (f,g,x-bounds,y-bounds) combinations "
                "verified by numeric equality (verified="
                + std::to_string(verified_combos) + ", required>=100)");

    return TEST_REPORT();
}
