#include "symbolic_ode.hpp"
#include "internal/ode_characteristic_roots.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/normalization_utils.hpp"
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "assumption_context.hpp"
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

// Separation only describes regions where g(y) is finite and nonzero.
// In that region reciprocal integer powers can be combined without extending
// the domain of the original ODE, unlike unconditional normalization.
std::shared_ptr<const SymbolicNode> reciprocal_on_regular_region(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (const auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> factors;
        factors.reserve(product->operands().size());
        for (const auto& factor : product->operands()) {
            factors.push_back(reciprocal_on_regular_region(factor));
        }
        return detail::make_node<MultiplyNode>(std::move(factors));
    }
    if (const auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        BigInt exponent;
        if (try_get_integer_value(
                std::dynamic_pointer_cast<const NumberNode>(power->exponent()), exponent)) {
            return detail::make_node<PowerNode>(
                power->base(), detail::make_node<NumberNode>(-exponent));
        }
    }
    return detail::make_node<PowerNode>(node, detail::make_node<NumberNode>(BigInt(-1)));
}

}

/// Check if the dependent variable is known Positive in the given context.
static bool dep_var_is_positive(const std::string& y, const AssumptionContext* ctx) {
    if (!ctx) {
        return false;
    }
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

    if (!is_separable(rhs, x, y)) {
        return nullptr;
    }

    std::shared_ptr<SymbolicExpr> f;
    std::shared_ptr<SymbolicExpr> g;
    const auto& root = detail::node(rhs);
    if (!expression_depends_on_variable(root, y)) {
        f = rhs;
        g = SymbolicExpr::number(1);
    } else if (!expression_depends_on_variable(root, x)) {
        f = SymbolicExpr::number(1);
        g = rhs;
    } else {
        // The shared classifier guarantees a product of uncoupled factors here.
        const auto& product = static_cast<const MultiplyNode&>(*root);
        std::vector<std::shared_ptr<const SymbolicNode>> x_factors;
        std::vector<std::shared_ptr<const SymbolicNode>> y_factors;
        for (const auto& factor : product.operands()) {
            auto& factors = expression_depends_on_variable(factor, y) ? y_factors : x_factors;
            factors.push_back(factor);
        }
        f = detail::make_expression_ptr(detail::make_node<MultiplyNode>(std::move(x_factors)));
        g = detail::make_expression_ptr(detail::make_node<MultiplyNode>(std::move(y_factors)));
    }

    auto inv_y = detail::make_expression_ptr(reciprocal_on_regular_region(detail::node(g)))->simplify();
    auto int_y = inv_y->integrate(y);
    auto int_x = f->simplify()->integrate(x);
    auto result = SymbolicExpr::add(
        int_y, SymbolicExpr::multiply(SymbolicExpr::number(-1), int_x))->simplify();

    // Preserve the existing positive-context branch marker.
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

static std::shared_ptr<SymbolicExpr> linear2_root_basis(
    const std::vector<ode_root_detail::CharRoot>& roots,
    const std::string& x) {
    auto variable = SymbolicExpr::variable(x);
    auto solution = SymbolicExpr::number(0);
    int constant_index = 1;
    for (const auto& root : roots) {
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
                solution, SymbolicExpr::add(cosine_term, sine_term));
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
    return solution;
}

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

    auto solution = linear2_root_basis(roots.value(), x);
    solution = solution->simplify();
    if (dep_var_is_positive(y, ctx)) {
        solution = make_abs(solution);
    }
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
    if (!step) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(step.error());
    }

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
