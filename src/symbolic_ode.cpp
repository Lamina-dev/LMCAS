#include "internal/expression_transform.hpp"
#include "symbolic_ode.hpp"
#include "internal/ode_characteristic_roots.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/normalization_utils.hpp"
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "assumption_context.hpp"
#include "integrator.hpp"
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
static Result<bool> dep_var_is_positive(
    const std::string& y, const AssumptionContext* ctx, ComputationContext& context) {
    if (!ctx) {
        return false;
    }
    auto y_expr = SymbolicExpr::variable(y);
    auto positive = ctx->is_positive_checked(*y_expr, context);
    if (!positive) { return Result<bool>::failure(positive.error()); }
    return positive.value() == Tribool::True;
}

/// Wrap an expression in abs() to signal positive-branch preference.
static std::shared_ptr<SymbolicExpr> make_abs(std::shared_ptr<SymbolicExpr> expr) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(expr)}));
}

// solve_separable_ode

ExprPtr solve_separable_ode(ExprPtr rhs, const std::string& x,
                           const std::string& y, const AssumptionContext* ctx) {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");

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

    auto inv_y = detail::make_expression_ptr(detail::normalize_node(reciprocal_on_regular_region(detail::node(g)), transform_budget));
    Integrator integrator;
    auto integrated_y = integrator.integrate_checked(*inv_y, y, context);
    if (!integrated_y) { throw integrated_y.error(); }
    auto int_y = std::make_shared<SymbolicExpr>(std::move(integrated_y.value()));
    auto normalized_f = detail::make_expression_ptr(
        detail::normalize_node(detail::node(f), transform_budget));
    auto integrated_x = integrator.integrate_checked(*normalized_f, x, context);
    if (!integrated_x) { throw integrated_x.error(); }
    auto int_x = std::make_shared<SymbolicExpr>(std::move(integrated_x.value()));
    auto result = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::add(
        int_y, SymbolicExpr::multiply(SymbolicExpr::number(-1), int_x))), transform_budget));

    // Preserve the existing positive-context branch marker.
    auto positive = dep_var_is_positive(y, ctx, context);
    if (!positive) { throw positive.error(); }
    if (positive.value()) {
        result = make_abs(result);
    }

    return result;
}


// solve_linear1_ode

ExprPtr solve_linear1_ode(ExprPtr Px, ExprPtr Qx, const std::string& x,
                         const std::string& y, const AssumptionContext* ctx) {
    ComputationContext context;
    Integrator integrator;
    auto integrated_p = integrator.integrate_checked(*Px, x, context);
    if (!integrated_p) { throw integrated_p.error(); }
    auto intP = std::make_shared<SymbolicExpr>(std::move(integrated_p.value()));
    auto mu = SymbolicExpr::exp(intP);

    auto Qmu = SymbolicExpr::multiply(Qx, mu);
    auto integrated_q = integrator.integrate_checked(*Qmu, x, context);
    if (!integrated_q) { throw integrated_q.error(); }
    auto intQmu = std::make_shared<SymbolicExpr>(std::move(integrated_q.value()));

    auto C = SymbolicExpr::variable("C");
    auto num = SymbolicExpr::add(intQmu, C);
    auto result = SymbolicExpr::divide(num, mu);

    // When the dependent variable is known Positive, prefer positive branch.
    auto positive = dep_var_is_positive(y, ctx, context);
    if (!positive) { throw positive.error(); }
    if (positive.value()) {
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
    const AssumptionContext* ctx, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            CasErrc::InvalidArgument,
            "solve_linear2_ode: coefficients must be finite",
            kSolveLinear2OdeOperation);
    }
    ExprPtr solution;
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
        // 常系数退化方程以特征根 -c/b 构造指数基，沿用数值系数的求解语义。
        solution = SymbolicExpr::multiply(SymbolicExpr::variable("C"),
            SymbolicExpr::exp(SymbolicExpr::multiply(
                SymbolicExpr::number(-normalized_c), SymbolicExpr::variable(x))));
    } else {
        auto roots = ode_root_detail::find_characteristic_roots(
            {a, b, c}, kSolveLinear2OdeOperation);
        if (!roots) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(roots.error());
        }
        solution = linear2_root_basis(roots.value(), x);
    }
    solution = detail::make_expression_ptr(detail::normalize_node(detail::node(solution), transform_budget));
    auto positive = dep_var_is_positive(y, ctx, context);
    if (!positive) { return ExpressionResult::failure(positive.error()); }
    if (positive.value()) {
        solution = make_abs(solution);
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(
        std::move(solution));
} catch (const CasError& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, "allocation failed", "solve_linear2_homogeneous");
} catch (const std::length_error& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, error.what(), "solve_linear2_homogeneous");
} catch (const std::exception& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::InternalInvariant, error.what(), "solve_linear2_homogeneous");
}

std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx
) {
    ComputationContext context;
    if (!fx || !LMCAS::detail::node(fx)) {
        throw std::invalid_argument(
            "solve_linear2_ode: forcing expression must not be null");
    }
    if (!fx->is_zero()) {
        throw std::logic_error(
            "solve_linear2_ode: non-homogeneous case is outside the current support domain");
    }
    auto solution = solve_linear2_homogeneous(a, b, c, x, y, ctx, context);
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

        return solve_linear2_homogeneous(a, b, c, x, y, ctx, context);
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
