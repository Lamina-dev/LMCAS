#include "../include/symbolic_ode.hpp"
#include "internal/ode_characteristic_roots.hpp"
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

namespace LMCAS {

namespace {

constexpr const char* kSolveLinear2OdeOperation = "solve_linear2_ode";

} // namespace

/// Check if the dependent variable is known Positive in the given context.
static bool dep_var_is_positive(const std::string& y, const AssumptionContext* ctx) {
    if (!ctx) return false;
    auto y_expr = SymbolicExpr::variable(y);
    auto positive = ctx->is_positive(*y_expr);
    return positive && positive.value() == Tribool::True;
}

/// Wrap an expression in abs() to signal positive-branch preference.
static std::shared_ptr<SymbolicExpr> make_abs(std::shared_ptr<SymbolicExpr> expr) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(expr)}));
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

static Result<std::shared_ptr<SymbolicExpr>>
solve_linear2_homogeneous(
    double a, double b, double c,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx)
{
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InvalidArgument,
            "solve_linear2_ode: coefficients must be finite",
            kSolveLinear2OdeOperation);
    }
    if (a == 0.0) {
        if (b == 0.0) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(
                CasErrc::InvalidArgument,
                "solve_linear2_ode: leading and first-derivative coefficients are both zero",
                kSolveLinear2OdeOperation);
        }
        const double normalized_c = c / b;
        if (!std::isfinite(normalized_c)) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(
                CasErrc::NumericFailure,
                "solve_linear2_ode: degenerate coefficient normalization is non-finite",
                kSolveLinear2OdeOperation);
        }
        auto solution = solve_linear1_ode(
            SymbolicExpr::number(normalized_c),
            SymbolicExpr::number(0), x, y, ctx);
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            std::move(solution));
    }

    auto roots = ode_root_detail::find_characteristic_roots(
        {a, b, c}, kSolveLinear2OdeOperation);
    if (!roots) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(roots.error());
    }

    auto variable = SymbolicExpr::variable(x);
    auto solution = SymbolicExpr::number(0);
    int constant_index = 1;
    for (const auto& root : roots.value()) {
        auto exponential = SymbolicExpr::exp(
            SymbolicExpr::multiply(
                SymbolicExpr::number(root.real_part), variable));
        if (root.is_complex) {
            auto argument = SymbolicExpr::multiply(
                SymbolicExpr::number(root.imag_part), variable);
            auto cosine_basis = SymbolicExpr::multiply(
                exponential, SymbolicExpr::cos(argument));
            auto sine_basis = SymbolicExpr::multiply(
                exponential, SymbolicExpr::sin(argument));
            auto cosine_term = SymbolicExpr::multiply(
                SymbolicExpr::variable(
                    "C" + std::to_string(constant_index++)),
                cosine_basis);
            auto sine_term = SymbolicExpr::multiply(
                SymbolicExpr::variable(
                    "C" + std::to_string(constant_index++)),
                sine_basis);
            solution = SymbolicExpr::add(
                solution,
                SymbolicExpr::add(cosine_term, sine_term));
            continue;
        }
        for (int power = 0; power < root.multiplicity; ++power) {
            auto basis = exponential;
            if (power > 0) {
                basis = SymbolicExpr::multiply(
                    SymbolicExpr::power(
                        variable, SymbolicExpr::number(power)),
                    exponential);
            }
            solution = SymbolicExpr::add(
                solution,
                SymbolicExpr::multiply(
                    SymbolicExpr::variable(
                        "C" + std::to_string(constant_index++)),
                    basis));
        }
    }
    solution = solution->simplify();
    if (dep_var_is_positive(y, ctx)) solution = make_abs(solution);
    return Result<std::shared_ptr<SymbolicExpr>>::success(
        std::move(solution));
}

std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {
    if (!fx || !LMCAS::detail::node(fx)) {
        throw std::invalid_argument(
            "solve_linear2_ode: forcing expression must not be null");
    }
    if (!fx->is_zero()) {
        throw std::logic_error(
            "solve_linear2_ode: non-homogeneous case is outside the current support domain");
    }
    auto solution = solve_linear2_homogeneous(a, b, c, x, y, ctx);
    if (!solution) {
        if (solution.error().code == CasErrc::InvalidArgument) {
            throw std::invalid_argument(solution.error().message);
        }
        if (solution.error().code == CasErrc::NumericFailure) {
            throw std::overflow_error(solution.error().message);
        }
        throw std::logic_error(solution.error().message);
    }
    return std::move(solution.value());
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

    if (!fx || !LMCAS::detail::node(fx)) {
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

        return solve_linear2_homogeneous(a, b, c, x, y, ctx);
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
