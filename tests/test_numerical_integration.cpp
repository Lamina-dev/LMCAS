#include "numerical_integration.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"
#include <cmath>
#include <limits>

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static bool close(const std::shared_ptr<SymbolicExpr>& e, double expected, double tol = 1e-6) {
    if (!e) return false;
    auto s = e->simplify();
    if (s->is_number()) return std::abs(s->to_numeric() - expected) < tol;
    // Fall back to recursive numeric evaluation (handles sqrt, etc.).
    auto v = test_numeric_eval(s);
    if (v.has_value()) return std::abs(*v - expected) < tol;
    return false;
}

static bool close(const ApproxReal& value, double expected,
                  double tolerance = 1e-6) {
    return value.status == NumericStatus::Finite &&
           std::abs(value.value - expected) < tolerance;
}

int main() {
    auto x = SymbolicExpr::variable("x");

    // ∫₀¹ x dx = 1/2 (Simpson exact for linear)
    {
        auto r = quadrature_simpson_numeric(
            x, "x", num(0), num(1), 10).value();
        EXPECT_TRUE(close(r, 0.5), "Simpson ∫₀¹ x dx = 0.5");
    }

    // ∫₀¹ x^2 dx = 1/3 (Simpson exact for quadratic)
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = quadrature_simpson_numeric(
            f, "x", num(0), num(1), 10).value();
        EXPECT_TRUE(close(r, 1.0/3.0), "Simpson ∫₀¹ x² dx = 1/3");
    }

    // Fixed Simpson derives its error estimate from an actual refinement.
    {
        auto f = SymbolicExpr::power(x, num(4));
        auto r = quadrature_simpson_numeric(f, "x", num(0), num(1), 10);
        const double observed_error =
            r ? std::abs(r.value().value - 0.2) : 0.0;
        EXPECT_TRUE(r && observed_error < 2e-5,
                    "checked fixed Simpson integrates x^4 accurately");
        EXPECT_TRUE(r && r.value().status == NumericStatus::Finite &&
                        r.value().absolute_error > 1e-8 &&
                        observed_error <= r.value().absolute_error + 1e-14,
                    "fixed Simpson reports its step-doubling error estimate");
    }

    // ∫₀¹ x^3 dx = 1/4 via Gauss-Legendre (n=2 exact up to degree 3)
    {
        auto f = SymbolicExpr::power(x, num(3));
        auto r = quadrature_gaussian_numeric(
            f, "x", num(0), num(1), 2).value();
        EXPECT_TRUE(close(r, 0.25, 1e-6), "Gauss ∫₀¹ x³ dx = 1/4");
    }

    // checked Gauss-Legendre reports the same value through Result<ApproxReal>
    {
        auto f = SymbolicExpr::power(x, num(3));
        auto r = quadrature_gaussian_numeric(f, "x", num(0), num(1), 2);
        EXPECT_TRUE(r && std::abs(r.value().value - 0.25) < 1e-12,
                    "checked Gauss ∫₀¹ x³ dx = 1/4");
    }

    // Orders above three must use LMMC Gauss-Legendre, not a Simpson fallback.
    {
        auto eighth_power = SymbolicExpr::power(x, num(8));
        auto order_five = quadrature_gaussian_numeric(
            eighth_power, "x", num(0), num(1), 5);
        EXPECT_TRUE(order_five &&
                        std::abs(order_five.value().value - 1.0 / 9.0) < 1e-13,
                    "five-point Gauss integrates degree-eight polynomials exactly");

        auto order_four = quadrature_gaussian_numeric(
            eighth_power, "x", num(0), num(1), 4);
        const double observed_error = order_four
            ? std::abs(order_four.value().value - 1.0 / 9.0)
            : 0.0;
        EXPECT_TRUE(order_four && order_four.value().absolute_error > 1e-8 &&
                        observed_error <= order_four.value().absolute_error + 1e-14,
                    "Gaussian error metadata comes from an adjacent-order rule");

        auto order_twenty = quadrature_gaussian_numeric(
            eighth_power, "x", num(0), num(1), 20);
        EXPECT_TRUE(order_twenty &&
                        std::abs(order_twenty.value().value - 1.0 / 9.0) < 1e-13,
                    "twenty-point Gauss is accepted and evaluated directly");

        auto unsupported = quadrature_gaussian_numeric(
            eighth_power, "x", num(0), num(1), 21);
        EXPECT_TRUE(!unsupported &&
                        unsupported.error().code == CasErrc::InvalidArgument,
                    "Gaussian orders above twenty are rejected explicitly");
    }

    // adaptive_simpson on ∫₀¹ x^2 dx = 1/3
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = adaptive_simpson_numeric(
            f, "x", num(0), num(1), 1e-10).value();
        EXPECT_TRUE(close(r, 1.0/3.0, 1e-8), "adaptive Simpson ∫₀¹ x² dx = 1/3");
    }

    // checked adaptive_simpson reports success with explicit error metadata
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = adaptive_simpson_numeric(f, "x", num(0), num(1), 1e-10);
        EXPECT_TRUE(r.has_value(), "checked adaptive Simpson returns Result success");
        EXPECT_TRUE(r.has_value() && std::abs(r.value().value - 1.0/3.0) < 1e-8,
                    "checked adaptive Simpson ∫₀¹ x² dx = 1/3");
        EXPECT_TRUE(r.has_value() && r.value().absolute_error >= 0.0 &&
                        std::abs(r.value().value - 1.0/3.0) <= r.value().absolute_error,
                    "reported Simpson error bounds the observed quadratic error");
    }

    // Crossing a singularity is a domain failure, not a principal-value success.
    {
        auto reciprocal = SymbolicExpr::divide(num(1), x);
        auto r = adaptive_simpson_numeric(
            reciprocal, "x", num(-1), num(1), 1e-10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::DomainError,
                    "adaptive Simpson rejects an interior pole");
    }

    // A caller-imposed recursion limit must not return an unverified estimate.
    {
        auto absolute = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(x)}));
        ComputationContext context;
        auto r = adaptive_simpson_numeric(
            absolute, "x", num(-1), num(1), context, 1e-14, 0);
        EXPECT_TRUE(!r && r.error().code == CasErrc::ResourceLimit,
                    "insufficient refinement depth returns ResourceLimit");
    }

    // Cancellation is checked through the shared evaluation context.
    {
        CancellationToken cancellation;
        cancellation.cancel();
        ComputationContext context({}, cancellation);
        auto r = adaptive_simpson_numeric(x, "x", num(0), num(1), context, 1e-10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::Cancelled,
                    "adaptive Simpson propagates cancellation");
    }

    // Finite endpoints can still define a width that overflows IEEE double.
    {
        const double largest = std::numeric_limits<double>::max();
        auto r = adaptive_simpson_numeric(
            x, "x", SymbolicExpr::number(-largest), SymbolicExpr::number(largest), 1e-10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::NumericFailure,
                    "unrepresentable interval width returns NumericFailure");
    }

    // Reversed limits preserve orientation.
    {
        auto square = SymbolicExpr::multiply(x, x);
        auto r = adaptive_simpson_numeric(square, "x", num(1), num(0), 1e-10);
        EXPECT_TRUE(r && std::abs(r.value().value + 1.0/3.0) <=
                            r.value().absolute_error + 1e-15,
                    "reversed limits negate the integral");
    }

    // checked path must not hide missing variables as zero/null.
    {
        auto y = SymbolicExpr::variable("y");
        auto f = SymbolicExpr::add(x, y);
        auto r = adaptive_simpson_numeric(f, "x", num(0), num(1), 1e-10);
        EXPECT_TRUE(!r.has_value(), "checked adaptive Simpson rejects unbound variables");
        EXPECT_TRUE(!r.has_value() && r.error().code == CasErrc::UnboundSymbol,
                    "checked adaptive Simpson reports UnboundSymbol");
    }

    // checked endpoint evaluation must also use explicit numeric errors.
    {
        auto r = adaptive_simpson_numeric(x, "x", SymbolicExpr::variable("a"), num(1), 1e-10);
        EXPECT_TRUE(!r.has_value(), "checked adaptive Simpson rejects symbolic bounds");
        EXPECT_TRUE(!r.has_value() && r.error().code == CasErrc::UnboundSymbol,
                    "checked adaptive Simpson reports symbolic bound as UnboundSymbol");
    }

    // Resource budgets must be honored before deep recursion starts.
    {
        ResourceLimits limits;
        limits.max_steps = 0;
        ComputationContext context(limits);
        auto r = adaptive_simpson_numeric(x, "x", num(0), num(1), context, 1e-10);
        EXPECT_TRUE(!r.has_value(), "checked adaptive Simpson honors exhausted step budget");
        EXPECT_TRUE(!r.has_value() && r.error().code == CasErrc::ResourceLimit,
                    "checked adaptive Simpson reports ResourceLimit");
    }

    // checked fixed Simpson propagates unbound variables.
    {
        auto y = SymbolicExpr::variable("y");
        auto f = SymbolicExpr::add(x, y);
        auto r = quadrature_simpson_numeric(f, "x", num(0), num(1), 10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::UnboundSymbol,
                    "checked fixed Simpson reports UnboundSymbol");
    }

    // checked fixed Simpson rejects domain errors at samples.
    {
        auto f = SymbolicExpr::ln(x);
        auto r = quadrature_simpson_numeric(f, "x", num(-1), num(1), 10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::DomainError,
                    "checked fixed Simpson reports sample DomainError");
    }

    // checked fixed quadrature validates arguments and observes cancellation.
    {
        auto invalid_n = quadrature_simpson_numeric(x, "x", num(0), num(1), 0);
        EXPECT_TRUE(!invalid_n && invalid_n.error().code == CasErrc::InvalidArgument,
                    "checked fixed Simpson rejects non-positive sample counts");

        auto odd_n = quadrature_simpson_numeric(x, "x", num(0), num(1), 9);
        EXPECT_TRUE(!odd_n && odd_n.error().code == CasErrc::InvalidArgument,
                    "checked fixed Simpson rejects odd subinterval counts");

        CancellationToken cancellation;
        cancellation.cancel();
        ComputationContext context({}, cancellation);
        auto cancelled = quadrature_gaussian_numeric(x, "x", num(0), num(1), context, 2);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked Gaussian observes cancellation");
    }


    // checked fixed quadrature consumes the shared step budget.
    {
        ResourceLimits limits;
        limits.max_steps = 1;
        ComputationContext context(limits);
        auto r = quadrature_simpson_numeric(x, "x", num(0), num(1), context, 10);
        EXPECT_TRUE(!r && r.error().code == CasErrc::ResourceLimit,
                    "checked fixed Simpson reports ResourceLimit");
    }

    // adaptive_simpson on a transcendental: ∫₀^π sin(x) dx = 2
    {
        auto f = SymbolicExpr::sin(x);
        auto pi = SymbolicExpr::number(LMMC_CONST_PI);
        auto r = adaptive_simpson_numeric(
            f, "x", num(0), pi, 1e-10).value();
        EXPECT_TRUE(close(r, 2.0, 1e-6), "adaptive Simpson ∫₀^π sin x dx = 2");
    }

    // numerical_integrate wrapper
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = numerical_integrate_numeric(
            f, "x", num(0), num(1), 100).value();
        EXPECT_TRUE(close(r, 1.0/3.0, 1e-6), "numerical_integrate ∫₀¹ x² dx = 1/3");
    }

    // checked numerical_integrate wrapper uses the explicit fixed Simpson contract.
    {
        auto f = SymbolicExpr::multiply(x, x);
        auto r = numerical_integrate_numeric(f, "x", num(0), num(1), 100);
        EXPECT_TRUE(r && std::abs(r.value().value - 1.0 / 3.0) < 1e-10,
                    "checked numerical_integrate ∫₀¹ x² dx = 1/3");
    }

    return TEST_REPORT();
}
