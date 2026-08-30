#include "../include/symbolic_ode.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/assumption_context.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <cmath>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

namespace lamina {

namespace {

constexpr const char* kSolveLinear2OdeOperation = "solve_linear2_ode";

} // namespace

/// Check if the dependent variable is known Positive in the given context.
static bool dep_var_is_positive(const std::string& y, const AssumptionContext* ctx) {
    if (!ctx) return false;
    auto y_expr = SymbolicExpr::variable(y);
    return detail::propagate_result(ctx->is_positive(*y_expr)) == Tribool::True;
}

/// Check if a symbolic expression is known NonZero in the given context.
static bool expr_is_nonzero(const std::shared_ptr<SymbolicExpr>& expr, const AssumptionContext* ctx) {
    if (!ctx || !expr) return false;
    return detail::propagate_result(ctx->is_nonzero(*expr)) == Tribool::True;
}

/// Wrap an expression in abs() to signal positive-branch preference.
static std::shared_ptr<SymbolicExpr> make_abs(std::shared_ptr<SymbolicExpr> expr) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(expr)}));
}

// solve_separable_ode

std::shared_ptr<SymbolicExpr> solve_separable_ode(
    std::shared_ptr<SymbolicExpr> rhs,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {

    auto inv_y = SymbolicExpr::divide(SymbolicExpr::number(1), rhs);
    auto int_y = inv_y->integrate(y);
    auto int_x = SymbolicExpr::number(1)->integrate(x);

    auto result = SymbolicExpr::add(int_y, SymbolicExpr::multiply(SymbolicExpr::number(-1), int_x));

    // When the dependent variable is known Positive, prefer the positive
    // solution branch by wrapping in abs() (which simplifies to identity for
    // positive expressions, signaling downstream that only positive values apply).
    if (dep_var_is_positive(y, ctx)) {
        result = make_abs(result);
    }

    return result;
}

// solve_linear1_ode

std::shared_ptr<SymbolicExpr> solve_linear1_ode(
    std::shared_ptr<SymbolicExpr> Px,
    std::shared_ptr<SymbolicExpr> Qx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {

    auto intP = Px->integrate(x);
    auto mu = SymbolicExpr::exp(intP);

    auto Qmu = SymbolicExpr::multiply(Qx, mu);
    auto intQmu = Qmu->integrate(x);

    auto C = SymbolicExpr::variable("C");
    auto num = SymbolicExpr::add(intQmu, C);
    auto result = SymbolicExpr::divide(num, mu);

    // When the dependent variable is known Positive, prefer positive branch.
    if (dep_var_is_positive(y, ctx)) {
        result = make_abs(result);
    }

    return result;
}

// solve_linear2_ode

std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {

    // Determine whether the leading coefficient 'a' is known NonZero via
    // assumptions. If so, skip the zero-coefficient degenerate case check.
    bool a_known_nonzero = false;
    if (ctx) {
        auto a_expr = SymbolicExpr::number(a);
        a_known_nonzero = expr_is_nonzero(a_expr, ctx);
    }

    if (!a_known_nonzero) {
        // Guard a == 0: equation degenerates to first-order (or constant) form.
        int a_is_zero = 0;
        lmmc_double_nearly_equal_tol(a, 0.0, 1e-12, 1e-12, &a_is_zero);
        if (a_is_zero) {
            int b_is_zero = 0;
            lmmc_double_nearly_equal_tol(b, 0.0, 1e-12, 1e-12, &b_is_zero);
            if (b_is_zero) {
                // a == 0 and b == 0: not a proper second-order ODE.
                throw std::invalid_argument(
                    "solve_linear2_ode: leading and first-derivative coefficients are both zero");
            }
            // a == 0, b != 0: degenerates to b*y' + c*y = f(x), i.e. y' + (c/b)*y = f/b.
            auto Px = SymbolicExpr::number(c / b);
            auto Qx = fx->is_zero()
                          ? SymbolicExpr::number(0)
                          : SymbolicExpr::divide(fx, SymbolicExpr::number(b));
            return solve_linear1_ode(Px, Qx, x, y, ctx);
        }
    }

    double D = b*b - 4*a*c;
    auto C1 = SymbolicExpr::variable("C1");
    auto C2 = SymbolicExpr::variable("C2");
    std::shared_ptr<SymbolicExpr> yh;
    int eq;
    lmmc_double_nearly_equal_tol(D, 0.0, 1e-12, 1e-12, &eq);
    if (!eq && D > 0) {
        double r1 = (-b + std::sqrt(D)) / (2*a);
        double r2 = (-b - std::sqrt(D)) / (2*a);
        yh = SymbolicExpr::add(
            SymbolicExpr::multiply(C1, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r1), SymbolicExpr::variable(x)))),
            SymbolicExpr::multiply(C2, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r2), SymbolicExpr::variable(x))))
        );
    } else if (eq) {
        double r = -b / (2*a);
        yh = SymbolicExpr::add(
            SymbolicExpr::multiply(C1, SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r), SymbolicExpr::variable(x)))),
            SymbolicExpr::multiply(C2, SymbolicExpr::multiply(SymbolicExpr::variable(x), SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(r), SymbolicExpr::variable(x)))))
        );
    } else {
        double real = -b / (2*a);
        double imag = std::sqrt(-D) / (2*a);

        auto exp_part = SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(real), SymbolicExpr::variable(x)));
        auto cos_part = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        auto sin_part = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(imag), SymbolicExpr::variable(x)));
        yh = SymbolicExpr::multiply(exp_part, SymbolicExpr::add(SymbolicExpr::multiply(C1, cos_part), SymbolicExpr::multiply(C2, sin_part)));
    }

    if (!fx->is_zero()) {
        /// 非齐次特解当前位于支持域之外;显式诊断保留齐次解与通解的语义边界.
        throw std::logic_error(
            "solve_linear2_ode: non-homogeneous case is outside the current support domain");
    }

    // When the dependent variable is known Positive, prefer positive branch.
    if (dep_var_is_positive(y, ctx)) {
        yh = make_abs(yh);
    }

    return yh;
}

Result<std::shared_ptr<SymbolicExpr>> solve_linear2_ode_checked(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const AssumptionContext* ctx
) {
    auto step = context.consume_steps(1, kSolveLinear2OdeOperation);
    if (!step) return Result<std::shared_ptr<SymbolicExpr>>::failure(step.error());

    if (!fx || !lamina::detail::node(fx)) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InvalidArgument,
            "forcing expression must not be null",
            kSolveLinear2OdeOperation);
    }
    if (x.empty() || y.empty()) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InvalidArgument,
            "ODE variables must not be empty",
            kSolveLinear2OdeOperation);
    }

    try {
        if (!fx->is_zero()) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(
                CasErrc::Inconclusive,
                "non-homogeneous second-order constant-coefficient ODEs are outside the current support domain",
                kSolveLinear2OdeOperation);
        }

        return Result<std::shared_ptr<SymbolicExpr>>::success(
            solve_linear2_ode(a, b, c, std::move(fx), x, y, ctx));
    } catch (const std::invalid_argument& ex) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InvalidArgument, ex.what(), kSolveLinear2OdeOperation);
    } catch (const std::bad_alloc&) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::ResourceLimit,
            "ODE solving allocation failed",
            kSolveLinear2OdeOperation);
    } catch (const std::logic_error& ex) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::Inconclusive, ex.what(), kSolveLinear2OdeOperation);
    } catch (const std::exception& ex) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InternalInvariant, ex.what(), kSolveLinear2OdeOperation);
    }
}

Result<std::shared_ptr<SymbolicExpr>> solve_linear2_ode_checked(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {
    ComputationContext context;
    return solve_linear2_ode_checked(a, b, c, std::move(fx), x, y, context, ctx);
}

}
