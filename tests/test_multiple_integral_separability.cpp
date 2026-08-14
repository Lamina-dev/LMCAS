
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
using lamina::IntegrationStep;

namespace {

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
    lamina::ComputationContext context;

    std::vector<IntegrationStep> steps = {
        {"x", x_lo, x_hi},
        {"y", y_lo, y_hi},
    };

    std::shared_ptr<SymbolicExpr> engine_result;
    auto integrated = lamina::integrate_multiple_checked(
        *integrand, steps, integrator, context);
    if (!integrated) {
        rep.failed = true;
        rep.detail = "multiple integration failed: " + integrated.error().message;
        return rep;
    }
    engine_result = lamina::detail::make_expression_ptr(integrated.value());

    if (!engine_result) {
        rep.failed = true;
        rep.detail = "engine.evaluate returned null";
        return rep;
    }
    if (has_integral_node(lamina::detail::node(engine_result))) {
        rep.unevaluated = true;
        rep.detail = "engine result contains unevaluated integral";
        return rep;
    }

    // Compute Ix = integrate_def(f, x, x_lo, x_hi)
    std::shared_ptr<SymbolicExpr> Ix;
    std::shared_ptr<SymbolicExpr> Iy;
    try {
        // Build f(x) and g(y) afresh (they were consumed when wrapped as
        // children of the multiply node above; rebuild for safety).
        auto fx2 = f.build(SymbolicExpr::variable("x"));
        auto gy2 = g.build(SymbolicExpr::variable("y"));
        Ix = lamina::detail::make_expression_ptr(
            integrator.integrate_def(*fx2, "x", *x_lo, *x_hi));
        Iy = lamina::detail::make_expression_ptr(
            integrator.integrate_def(*gy2, "y", *y_lo, *y_hi));
    } catch (const std::exception& e) {
        rep.failed = true;
        rep.detail = std::string("exception in integrate_def: ") + e.what();
        return rep;
    }

    if (has_integral_node(lamina::detail::node(Ix)) || has_integral_node(lamina::detail::node(Iy))) {
        rep.unevaluated = true;
        rep.detail = "per-variable integrate_def left unevaluated integral";
        return rep;
    }

    // Numerically evaluate engine result and the product Ix * Iy.
    auto engine_simp = engine_result->simplify();
    if (!engine_simp) engine_simp = engine_result;

    auto Ix_simp = Ix->simplify();
    auto Iy_simp = Iy->simplify();
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
    TEST_CASE("Multiple integral separability");

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
