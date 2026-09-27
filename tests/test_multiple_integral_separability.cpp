
#include "test_common.hpp"
#include "integration.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <functional>

using namespace LMCAS;

using LMCAS::IntegrationStep;
using LMCAS::Integrator;

namespace {

constexpr double kTolerance = 1e-10;

// Each entry builds f(arg) given the argument expression.
struct FunctionSpec {
    std::string name;
    std::function<std::shared_ptr<SymbolicExpr>(std::shared_ptr<SymbolicExpr>)> build;
};

const std::vector<FunctionSpec> &functions() {
    static const std::vector<FunctionSpec> F = {
        {"t", [](std::shared_ptr<SymbolicExpr> a) { return a; }},
        {"t^2", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(a, SymbolicExpr::number(2)); }},
        {"t^3", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::power(a, SymbolicExpr::number(3)); }},
        {"sin(t)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::sin(a); }},
        {"cos(t)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::cos(a); }},
        {"exp(t)", [](std::shared_ptr<SymbolicExpr> a) { return SymbolicExpr::exp(a); }},
    };
    return F;
}

struct BoundSpec {
    BigInt lo;
    BigInt hi;
    std::string label;
};

const std::vector<BoundSpec> &bounds() {
    static const std::vector<BoundSpec> B = {
        {0, 1, "[0,1]"},
        {1, 2, "[1,2]"},
        {-1, 1, "[-1,1]"},
    };
    return B;
}

struct ComboReport {
    bool unevaluated = false;
    bool skipped = false;
    bool failed = false;
    double engine_value = 0.0;
    double product_value = 0.0;
    std::string detail;
};

void compare_numeric_integrals(
    const std::shared_ptr<SymbolicExpr> &engine_simp,
    const std::shared_ptr<SymbolicExpr> &Ix_simp,
    const std::shared_ptr<SymbolicExpr> &Iy_simp, ComboReport &rep) {
    /**
     * @brief 使用 test_common.hpp 的 test_numeric_eval 递归求值 Add、Multiply、Power、Function 节点。
     * SymbolicExpr::to_numeric() 仅处理 NumberNode 和单参数 FunctionNode，
     * 对 (1/2)*sin(1) 等复合表达式返回 0。
     */
    auto engine_opt = test_numeric_eval(engine_simp);
    auto Ix_opt = test_numeric_eval(Ix_simp);
    auto Iy_opt = test_numeric_eval(Iy_simp);

    if (!engine_opt || !Ix_opt || !Iy_opt) {
        rep.skipped = true;
        std::ostringstream oss;
        oss << "non-evaluable expression (engine="
            << (engine_opt ? std::to_string(*engine_opt) : "n/a")
            << ", Ix=" << (Ix_opt ? std::to_string(*Ix_opt) : "n/a")
            << ", Iy=" << (Iy_opt ? std::to_string(*Iy_opt) : "n/a") << ")";
        rep.detail = oss.str();
        return;
    }

    double engine_val = *engine_opt;
    double product_val = (*Ix_opt) * (*Iy_opt);

    if (!std::isfinite(engine_val) || !std::isfinite(product_val)) {
        rep.skipped = true;
        std::ostringstream oss;
        oss << "non-finite numeric value (engine=" << engine_val
            << ", Ix*Iy=" << product_val << ")";
        rep.detail = oss.str();
        return;
    }

    rep.engine_value = engine_val;
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
}

void compare_integral_results(
    const std::shared_ptr<SymbolicExpr> &engine_result,
    const std::shared_ptr<SymbolicExpr> &Ix,
    const std::shared_ptr<SymbolicExpr> &Iy, ComboReport &rep) {
    auto engine_simp = engine_result->simplify();
    if (!engine_simp) {
        engine_simp = engine_result;
    }

    auto Ix_simp = Ix->simplify();
    auto Iy_simp = Iy->simplify();
    if (!Ix_simp || !Iy_simp) {
        rep.failed = true;
        rep.detail = "simplify of per-variable integral returned null";
        return;
    }
    compare_numeric_integrals(engine_simp, Ix_simp, Iy_simp, rep);
}

void verify_separated_product(
    const FunctionSpec &f, const FunctionSpec &g,
    const std::shared_ptr<SymbolicExpr> &x_lo,
    const std::shared_ptr<SymbolicExpr> &x_hi,
    const std::shared_ptr<SymbolicExpr> &y_lo,
    const std::shared_ptr<SymbolicExpr> &y_hi, Integrator &integrator,
    const std::shared_ptr<SymbolicExpr> &engine_result, ComboReport &rep) {
    auto fx2 = f.build(SymbolicExpr::variable("x"));
    auto gy2 = g.build(SymbolicExpr::variable("y"));
    auto ix_result = integrator.integrate_def(*fx2, "x", *x_lo, *x_hi);
    if (!ix_result) {
        rep.failed = true;
        rep.detail = std::string("x integration failed: ") + ix_result.error().message;
        return;
    }
    auto iy_result = integrator.integrate_def(*gy2, "y", *y_lo, *y_hi);
    if (!iy_result) {
        rep.failed = true;
        rep.detail = std::string("y integration failed: ") + iy_result.error().message;
        return;
    }
    auto Ix = LMCAS::detail::make_expression_ptr(ix_result.value());
    auto Iy = LMCAS::detail::make_expression_ptr(iy_result.value());

    if (LMCAS::detail::contains_node_type<IntegralNode>(LMCAS::detail::node(Ix)) ||
        LMCAS::detail::contains_node_type<IntegralNode>(LMCAS::detail::node(Iy))) {
        rep.unevaluated = true;
        rep.detail = "per-variable integrate_def left unevaluated integral";
        return;
    }
    compare_integral_results(engine_result, Ix, Iy, rep);
}

ComboReport verify_combo(const FunctionSpec &f,
                         const FunctionSpec &g,
                         const BoundSpec &bx,
                         const BoundSpec &by) {
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
    LMCAS::ComputationContext context;

    std::vector<IntegrationStep> steps = {
        {"x", x_lo, x_hi},
        {"y", y_lo, y_hi},
    };

    std::shared_ptr<SymbolicExpr> engine_result;
    auto integrated = LMCAS::integrate_multiple_checked(
        *integrand, steps, integrator, context);
    if (!integrated) {
        rep.failed = true;
        rep.detail = "multiple integration failed: " + integrated.error().message;
        return rep;
    }
    engine_result = LMCAS::detail::make_expression_ptr(integrated.value());

    if (!engine_result) {
        rep.failed = true;
        rep.detail = "engine.evaluate returned null";
        return rep;
    }
    if (LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(engine_result))) {
        rep.unevaluated = true;
        rep.detail = "engine result contains unevaluated integral";
        return rep;
    }

    verify_separated_product(f, g, x_lo, x_hi, y_lo, y_hi,
                             integrator, engine_result, rep);
    return rep;
}

} // anonymous namespace

TEST(MultipleIntegralSeparability, ConstantIntegralsObeyComputationContext) {
    Integrator integrator;
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    std::vector<IntegrationStep> steps{{"x", zero, one}};

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto result = integrate_multiple_checked(*one, steps, integrator, cancelled);
    EXPECT_TRUE(!result && result.error().code == CasErrc::Cancelled) << "constant definite integral respects cancellation";

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    result = integrate_multiple_checked(*one, steps, integrator, exhausted);
    EXPECT_TRUE(!result && result.error().code == CasErrc::ResourceLimit) << "constant definite integral respects an exhausted step budget";

    limits.max_steps = 1;
    ComputationContext shared_budget(limits);
    steps.push_back({"y", zero, one});
    result = integrate_multiple_checked(*one, steps, integrator, shared_budget);
    EXPECT_TRUE(!result && result.error().code == CasErrc::ResourceLimit) << "successive constant integrals share the step budget";
}

TEST(MultipleIntegralSeparability, ConstantIntegrandUsesDefiniteEndpointValidation) {
    Integrator integrator;
    auto zero = SymbolicExpr::number(0);
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    ComputationContext direct_context;
    auto direct = integrator.integrate_def_checked(
        *y, "x", *zero, *x, direct_context);
    ASSERT_FALSE(direct.has_value());
    EXPECT_EQ(direct.error().code, CasErrc::InvalidArgument);

    ComputationContext multiple_context;
    std::vector<IntegrationStep> invalid_steps{{"x", zero, x}};
    auto multiple = integrate_multiple_checked(
        *y, invalid_steps, integrator, multiple_context);
    ASSERT_FALSE(multiple.has_value());
    EXPECT_EQ(multiple.error().code, direct.error().code)
        << "multiple integration applies the same endpoint contract as the definite API";
}

TEST(MultipleIntegralSeparability, ConstantIntegrandPreservesUnknownEqualEndpointDomain) {
    Integrator integrator;
    auto a = SymbolicExpr::variable("a");
    auto integrand = SymbolicExpr::sqrt(
        SymbolicExpr::variable("parameter"));

    ComputationContext direct_context;
    auto direct = integrator.integrate_def_checked(
        *integrand, "x", *a, *a, direct_context);
    ASSERT_TRUE(direct.has_value()) << direct.error().message;
    EXPECT_TRUE(LMCAS::detail::contains_node_type<IntegralNode>(
        LMCAS::detail::node(direct.value())));

    ComputationContext multiple_context;
    std::vector<IntegrationStep> steps{{"x", a, a}};
    auto multiple = integrate_multiple_checked(
        *integrand, steps, integrator, multiple_context);
    ASSERT_TRUE(multiple.has_value()) << multiple.error().message;
    EXPECT_TRUE(LMCAS::detail::contains_node_type<IntegralNode>(
        LMCAS::detail::node(multiple.value())));
}

TEST(MultipleIntegralSeparability, ConstantIntegrandComputesRectangleArea) {
    Integrator integrator;
    ComputationContext context;
    auto zero = SymbolicExpr::number(0);
    std::vector<IntegrationStep> steps{
        {"x", zero, SymbolicExpr::number(2)},
        {"y", zero, SymbolicExpr::number(3)}};
    auto result = integrate_multiple_checked(
        *SymbolicExpr::number(1), steps, integrator, context);
    ASSERT_TRUE(result.has_value()) << result.error().message;
    auto value = test_numeric_eval(
        LMCAS::detail::make_expression_ptr(result.value()));
    ASSERT_TRUE(value.has_value());
    EXPECT_NEAR(*value, 6.0, 1e-12);
}

TEST(MultipleIntegralSeparability, SupportsMoreThanThreeBinders) {
    Integrator integrator;
    LMCAS::ComputationContext context;
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    std::vector<IntegrationStep> steps{
        {"w", zero, one},
        {"z", zero, one},
        {"y", zero, one},
        {"x", zero, one}};
    auto fourfold = integrate_multiple_checked(
        *SymbolicExpr::number(1), steps, integrator, context);
    EXPECT_TRUE(fourfold.has_value()) << "fourfold exact integral succeeds";
    if (fourfold) {
        auto value = test_numeric_eval(
            LMCAS::detail::make_expression_ptr(fourfold.value()));
        EXPECT_TRUE(value && std::abs(*value - 1.0) < 1e-12) << "unit four-cube integrates to one";
    }
}

TEST(MultipleIntegralSeparability, SeparatedProductsMatchIndependentIntegrals) {
    std::size_t verified_combos = 0;
    for (const auto &f : functions()) {
        for (const auto &g : functions()) {
            for (const auto &bx : bounds()) {
                for (const auto &by : bounds()) {
                    SCOPED_TRACE("f=" + f.name + ", g=" + g.name + ", x=" + bx.label + ", y=" + by.label);
                    const ComboReport rep = verify_combo(f, g, bx, by);
                    EXPECT_FALSE(rep.failed) << rep.detail;
                    if (!rep.failed && !rep.unevaluated && !rep.skipped) {
                        ++verified_combos;
                    }
                }
            }
        }
    }
    EXPECT_GE(verified_combos, 100u)
        << "at least 100 separable (f,g,x-bounds,y-bounds) combinations "
           "verified by numeric equality";
}
