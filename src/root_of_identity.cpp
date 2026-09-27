#include "root_of_identity.hpp"
#include "polynomial_conversion.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/exact_root.hpp"

namespace LMCAS {

RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index,
    ComputationContext& context) {
    constexpr const char* operation = "rootof.construct";
    if (!polynomial || !LMCAS::detail::node(polynomial) ||
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
        std::move(*recognized.value()), index, context, operation);
    if (!identity) {
        return RootOfConstructionResult::failure(identity.error());
    }
    return RootOfConstructionResult::success(
        LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<RootOfNode>(
                std::move(identity.value()), variable)));
}

RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index) {
    ComputationContext context;
    return make_rootof_checked(polynomial, variable, index, context);
}
static std::shared_ptr<SymbolicExpr> symbolic_root_term(
    const std::shared_ptr<SymbolicExpr>& coefficient,
    const std::string& variable, int degree) {
    if (degree == 0) return coefficient;
    auto symbol = SymbolicExpr::variable(variable);
    auto power = degree == 1
        ? symbol : SymbolicExpr::power(symbol, SymbolicExpr::number(degree));
    return coefficient->is_one()
        ? power : SymbolicExpr::multiply(coefficient, power);
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
            if (LMCAS::detail::node(coeff_expr) && coeff_expr->is_zero()) {
                continue;
            }
        }
        terms.push_back(symbolic_root_term(coeff_expr, poly.variable_name, i));
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];

    auto result = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) {
        result = SymbolicExpr::add(result, terms[i]);
    }
    return result;
}


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
        solutions.push_back(LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<RootOfNode>(
                detail::ExactRootId{
                    canonical, static_cast<std::size_t>(index)},
                var)));
    }
    return solutions;
}

namespace {


Result<const RootOfNode*> parse_rootof_request(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context,
    const std::string& operation) {
    if (!rootof_expr || !LMCAS::detail::node(rootof_expr)) {
        return Result<const RootOfNode*>::failure(
            CasErrc::InvalidArgument, "RootOf expression cannot be null", operation);
    }
    auto root = std::dynamic_pointer_cast<const RootOfNode>(
        LMCAS::detail::node(rootof_expr));
    if (!root) {
        return Result<const RootOfNode*>::failure(
            CasErrc::InvalidArgument,
            "expression is not a valid RootOf", operation);
    }
    auto access = context.consume_steps(1, operation);
    if (!access) return Result<const RootOfNode*>::failure(access.error());
    return Result<const RootOfNode*>::success(root.get());
}


}

RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context)
{
    constexpr const char* operation = "rootof_evaluate";
    auto request = parse_rootof_request(rootof_expr, context, operation);
    if (!request) return RootOfEvaluationResult::failure(request.error());
    auto evaluated = detail::evaluate_root_real(
        request.value()->exact_id(),
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
        request.value()->exact_id(),
        detail::NumericEvaluationOptions{}, context);
}

RootOfComplexEvaluationResult rootof_evaluate_complex_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr) {
    ComputationContext context;
    return rootof_evaluate_complex_checked(rootof_expr, context);
}

}
