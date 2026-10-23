#include "internal/expression_transform.hpp"
#include "inequality_solver.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/numeric_probe.hpp"
#include "internal/inequality_solver_support.hpp"
#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace LMCAS::detail::inequality_support {

std::optional<double> try_checked_numeric_constant(
    const SymbolicExpr& expr, ComputationContext& context) {
    return LMCAS::detail::try_finite_numeric(expr, context);
}

int exact_numeric_sign(const std::shared_ptr<SymbolicExpr>& expr,
                       ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!expr || !LMCAS::detail::node(expr)) {
        return 0;
    }
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(expr), transform_budget));
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return 0;
    }
    auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(simplified));
    if (!num) {
        return 0;
    }
    if (std::holds_alternative<BigInt>(num->value())) {
        const auto& value = std::get<BigInt>(num->value());
        if (value.is_zero()) {
            return 0;
        }
        return value.is_negative() ? -1 : 1;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        const auto& value = std::get<Rational>(num->value());
        if (value.get_numerator().is_zero()) {
            return 0;
        }
        return value.get_numerator().is_negative() ? -1 : 1;
    }
    const auto value = std::get<lmmc_real_t>(num->value());
    if (!std::isfinite(value) || value == 0.0) {
        return 0;
    }
    return value < 0 ? -1 : 1;
}

int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& poly,
                           ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (poly.is_zero()) {
        return 0;
    }
    auto lc = poly.lead_coeff().val;
    if (!lc) {
        return 1;
    }
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(lc), transform_budget));
    if (!simplified) {
        return 1;
    }

    if (auto val = try_checked_numeric_constant(*simplified, context)) {
        if (*val > 0) {
            return 1;
        }
        if (*val < 0) {
            return -1;
        }
    }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(num->value())) {
            return std::get<BigInt>(num->value()).is_negative() ? -1 : 1;
        }
        if (std::holds_alternative<Rational>(num->value())) {
            return std::get<Rational>(num->value()).get_numerator().is_negative() ? -1 : 1;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            return std::get<lmmc_real_t>(num->value()) < 0 ? -1 : 1;
        }
    }
    return 1;
}

static std::shared_ptr<SymbolicExpr> snap_verified_integer_root(
    const Polynomial<Rational>& poly,
    const std::shared_ptr<SymbolicExpr>& root, ComputationContext& context) {
    if (!root || poly.is_zero()) {
        return root;
    }
    auto numeric = try_checked_numeric_constant(*root, context);
    if (!numeric) {
        return root;
    }

    const Rational candidate = Rational::from_double(*numeric);
    if (!candidate.is_integer()) {
        return root;
    }
    if (poly.eval(candidate) == Rational(0)) {
        return SymbolicExpr::number(candidate);
    }
    return root;
}

using RootMultiplicityList = std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>>;

static void append_verified_roots(
    RootMultiplicityList& result,
    const Polynomial<Rational>& polynomial,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    int multiplicity, ComputationContext& context) {
    for (const auto& root : roots) {
        if (!root) {
            continue;
        }
        auto verified = snap_verified_integer_root(polynomial, root, context);
        if (try_checked_numeric_constant(*verified, context)) {
            result.push_back({verified, multiplicity});
        }
    }
}

static void append_affine_root(RootMultiplicityList& result,
                               const Polynomial<Rational>& factor, int mult,
                               ComputationContext& context) {
    Rational a = factor.coeffs[1];
    Rational b = factor.coeffs[0];
    if (a != Rational(0)) {
        Rational root_val = Rational(0) - b / a;
        auto root_expr = SymbolicExpr::number(root_val);
        if (try_checked_numeric_constant(*root_expr, context)) {
            result.push_back({root_expr, mult});
        }
    }
}

static void append_quadratic_roots(RootMultiplicityList& result,
                                  const Polynomial<Rational>& factor, int mult,
                                  ComputationContext& context) {
    Rational a = factor.coeffs[2];
    Rational b = factor.coeffs[1];
    Rational c = factor.coeffs[0];

    Rational disc = b * b - Rational(4) * a * c;
    if (disc < Rational(0)) {
        return;
    }
    double disc_val = disc.to_double();

    double a_val = a.to_double();
    double b_val = b.to_double();
    double sqrt_disc = std::sqrt(disc_val);

    double r1;
    double r2;
    if (sqrt_disc == 0.0) {
        r1 = -b_val / (2.0 * a_val);
        r2 = r1;
    } else {
        const double q =
            -0.5 * b_val - 0.5 * std::copysign(sqrt_disc, b_val);
        r1 = q / a_val;
        r2 = c.to_double() / q;
    }

    if (std::isfinite(r1)) {
        auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r1), context);
        result.push_back({root, mult});
    }
    if (std::isfinite(r2) && r2 != r1) {
        auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r2), context);
        result.push_back({root, mult});
    }
}

static void append_factor_roots(RootMultiplicityList& result,
                                const Polynomial<Rational>& factor, int mult,
                                const std::string& variable, ComputationContext& context) {
    std::vector<SymbolicPolyCoeff> spc_coeffs;
    for (int i = 0; i <= factor.degree(); ++i) {
        spc_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(factor.coeffs[i])));
    }
    Polynomial<SymbolicPolyCoeff> factor_spc(spc_coeffs, variable);

    auto factor_roots = polynomial_solver_detail::solve_by_factoring_impl(
        factor_spc, variable, context);
    append_verified_roots(result, factor, factor_roots, mult, context);
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> find_roots_with_multiplicity(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable, ComputationContext& context) {

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> result;

    auto converted = symbolic_to_poly<Rational>(expr, variable);
    if (!converted) throw std::invalid_argument(converted.error().message);
    const auto& poly_rat = converted.value();
    if (poly_rat.is_zero() || poly_rat.degree() <= 0) {
        return result;
    }

    auto factors = square_free_factorization(poly_rat);
    if (factors.empty()) {
        auto symbolic = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
        if (!symbolic) throw std::invalid_argument(symbolic.error().message);
        const auto& poly_spc = symbolic.value();
        if (!poly_spc.is_zero() && poly_spc.degree() >= 1) {
            auto roots = polynomial_solver_detail::solve_by_factoring_impl(
                poly_spc, variable, context);
            append_verified_roots(result, poly_rat, roots, 1, context);
        }
        return result;
    }
    for (const auto& [factor, mult] : factors) {
        if (factor.degree() <= 0) {
            continue;
        }
        if (factor.degree() == 1) {
            append_affine_root(result, factor, mult, context);
        } else if (factor.degree() == 2) {
            append_quadratic_roots(result, factor, mult, context);
        } else {
            append_factor_roots(result, factor, mult, variable, context);
        }
    }
    return result;
}

bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expr,
                                  const std::vector<std::string>& parameters) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    for (const auto& param : parameters) {
        if (expression_depends_on_variable(LMCAS::detail::node(expr), param)) {
            return true;
        }
    }
    return false;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& variable, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");

    if (poly.is_zero() || poly.degree() < 1) {
        return {};
    }

    int deg = poly.degree();
    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) {
            return SymbolicExpr::number(0);
        }
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {

        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto root = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(neg_b, a)), transform_budget));
        return { root };
    }

    auto results = polynomial_solver_detail::solve_by_factoring_impl(poly, variable, context);
    return results;
}
}
