#include "root_of_utils.hpp"
#include "newton_raphson.hpp"
#include "numeric_evaluation.hpp"
#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/exact_root.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace lamina {

RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index,
    ComputationContext& context) {
    constexpr const char* operation = "rootof.construct";
    if (!polynomial || !lamina::detail::node(polynomial) ||
        variable.empty()) {
        return RootOfConstructionResult::failure(
            CasErrc::InvalidArgument,
            "RootOf requires a polynomial and named variable", operation);
    }
    auto recognized = recognize_rational_polynomial(
        *polynomial, variable, context);
    if (!recognized) {
        return RootOfConstructionResult::failure(recognized.error());
    }
    if (!recognized.value()) {
        return RootOfConstructionResult::failure(
            CasErrc::Inconclusive,
            "RootOf coefficients must be exact rational constants",
            operation);
    }
    auto identity = detail::make_exact_root_id(
        *recognized.value(), index, context, operation);
    if (!identity) {
        return RootOfConstructionResult::failure(identity.error());
    }
    return RootOfConstructionResult::success(
        lamina::detail::make_expression_ptr(
            lamina::detail::make_node<RootOfNode>(
                std::move(identity.value()), variable)));
}

RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index) {
    ComputationContext context;
    return make_rootof_checked(polynomial, variable, index, context);
}

static std::shared_ptr<SymbolicExpr> symbolic_poly_to_expr(
    const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return SymbolicExpr::number(0);

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for (int i = poly.degree(); i >= 0; --i) {
        const auto& coeff = poly.coeffs[i];
        if (coeff == SymbolicPolyCoeff(0)) continue;

        auto coeff_expr = coeff.val;
        if (coeff_expr) {
            if (auto simplified = coeff_expr->simplify()) {
                coeff_expr = simplified;
            }
            if (lamina::detail::node(coeff_expr) && coeff_expr->is_zero()) {
                continue;
            }
        }
        if (i == 0) {
            terms.push_back(coeff_expr);
        } else {
            auto var_expr = SymbolicExpr::variable(poly.variable_name);
            std::shared_ptr<SymbolicExpr> var_part;
            if (i == 1) {
                var_part = var_expr;
            } else {
                var_part = SymbolicExpr::power(var_expr, SymbolicExpr::number(i));
            }

            if (coeff_expr->is_one()) {
                terms.push_back(var_part);
            } else {
                terms.push_back(SymbolicExpr::multiply(coeff_expr, var_part));
            }
        }
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];

    auto result = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) {
        result = SymbolicExpr::add(result, terms[i]);
    }
    return result;
}

namespace {

struct RootOfRequest {
    Polynomial<Rational> polynomial;
    std::string variable;
    std::size_t index = 0;
};

Result<RootOfRequest> parse_rootof_request(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context,
    const std::string& operation) {
    if (!rootof_expr || !lamina::detail::node(rootof_expr)) {
        return Result<RootOfRequest>::failure(
            CasErrc::InvalidArgument, "RootOf expression cannot be null", operation);
    }
    auto root = std::dynamic_pointer_cast<const RootOfNode>(
        lamina::detail::node(rootof_expr));
    if (!root) {
        return Result<RootOfRequest>::failure(
            CasErrc::InvalidArgument,
            "expression is not a valid RootOf", operation);
    }
    auto access = context.consume_steps(1, operation);
    if (!access) return Result<RootOfRequest>::failure(access.error());
    return Result<RootOfRequest>::success(RootOfRequest{
        root->exact_id().polynomial, root->variable(), root->index()});
}


} // namespace

std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {
    if (poly.degree() <= 0) return {};
    auto polynomial_expression = symbolic_poly_to_expr(poly);
    ComputationContext context;
    auto recognized = recognize_rational_polynomial(
        *polynomial_expression, var, context);
    if (!recognized || !recognized.value()) return {};
    auto canonical = recognized.value()->square_free_part().make_monic();
    canonical.variable_name = "_root";
    std::vector<std::shared_ptr<SymbolicExpr>> solutions;
    solutions.reserve(static_cast<std::size_t>(canonical.degree()));
    for (int index = 0; index < canonical.degree(); ++index) {
        solutions.push_back(lamina::detail::make_expression_ptr(
            lamina::detail::make_node<RootOfNode>(
                detail::ExactRootId{
                    canonical, static_cast<std::size_t>(index)},
                var)));
    }
    return solutions;
}


RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context)
{
    constexpr const char* operation = "rootof_evaluate";
    auto request = parse_rootof_request(rootof_expr, context, operation);
    if (!request) return RootOfEvaluationResult::failure(request.error());
    auto evaluated = detail::evaluate_root_real(
        detail::ExactRootId{
            request.value().polynomial, request.value().index},
        detail::NumericEvaluationOptions{}, context);
    if (!evaluated) {
        return RootOfEvaluationResult::failure(evaluated.error());
    }
    return RootOfEvaluationResult::success(evaluated.value().value);
}

RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr) {
    ComputationContext context;
    return rootof_evaluate_checked(rootof_expr, context);
}

RootOfComplexEvaluationResult rootof_evaluate_complex_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context) {
    constexpr const char* operation = "rootof_evaluate_complex";
    auto request = parse_rootof_request(rootof_expr, context, operation);
    if (!request) {
        return RootOfComplexEvaluationResult::failure(request.error());
    }
    return detail::evaluate_root_complex(
        detail::ExactRootId{
            request.value().polynomial, request.value().index},
        detail::NumericEvaluationOptions{}, context);
}

RootOfComplexEvaluationResult rootof_evaluate_complex_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr) {
    ComputationContext context;
    return rootof_evaluate_complex_checked(rootof_expr, context);
}


static std::vector<std::shared_ptr<SymbolicExpr>> solve_closed_form_from_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var)
{
    int deg = poly.degree();
    if (deg <= 0) return {};

    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {

        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        return { SymbolicExpr::divide(neg_b, a)->simplify() };
    } else if (deg == 2) {
        auto a = get_coeff(2);
        auto b = get_coeff(1);
        auto c = get_coeff(0);

        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto four_ac = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(b2, SymbolicExpr::multiply(four_ac, SymbolicExpr::number(-1)));
        auto sqrt_delta = SymbolicExpr::sqrt(delta);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);

        auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
        auto x2 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, SymbolicExpr::multiply(sqrt_delta, SymbolicExpr::number(-1))), two_a)->simplify();
        return { x1, x2 };
    } else if (deg == 3) {
        auto a = get_coeff(3);
        auto b = get_coeff(2);
        auto c = get_coeff(1);
        auto d = get_coeff(0);
        return solve_cubic(a, b, c, d, var);
    } else if (deg == 4) {
        auto a = get_coeff(4);
        auto b = get_coeff(3);
        auto c = get_coeff(2);
        auto d = get_coeff(1);
        auto e = get_coeff(0);
        return solve_quartic(a, b, c, d, e, var);
    }

    return {};
}




std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    if (!rootof_expr || !lamina::detail::node(rootof_expr)) return rootof_expr;

    auto root = std::dynamic_pointer_cast<const RootOfNode>(
        lamina::detail::node(rootof_expr));
    if (!root) return rootof_expr;

    const std::string& var = root->variable();
    const int k = static_cast<int>(root->index());
    std::vector<SymbolicPolyCoeff> coefficients;
    coefficients.reserve(root->exact_id().polynomial.coeffs.size());
    for (const auto& coefficient : root->exact_id().polynomial.coeffs) {
        coefficients.emplace_back(SymbolicExpr::number(coefficient));
    }
    Polynomial<SymbolicPolyCoeff> sym_poly(coefficients, var);

    int degree = sym_poly.degree();
    if (degree <= 0) return rootof_expr;

    if (k < 0 || k >= degree) {
        return rootof_expr;
    }

    if (degree == 1) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        return roots.size() == 1 ? roots.front() : rootof_expr;
    }
    if (degree == 2) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        if (roots.size() != 2) return rootof_expr;
        return roots[static_cast<std::size_t>(1 - k)];
    }

    return rootof_expr;
}

}
