#include "lsr_expr.hpp"

#include <exception>
#include <optional>
#include <utility>
#include <vector>

#include "lsr_expr_internal.hpp"
#include "poly_utils.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

Rational polynomial_coeff_or_zero(const Polynomial<Rational>& polynomial,
                                  std::size_t degree) {
    return degree < polynomial.coeffs.size()
        ? polynomial.coeffs[degree]
        : Rational(0);
}

ExprPtr rational_expression(const Rational& value) {
    return SymbolicExpr::number(value);
}

bool exact_rational_sqrt(const Rational& value, Rational& root) {
    if (value < Rational(0)) return false;
    const BigInt numerator_root = value.get_numerator().sqrt();
    const BigInt denominator_root = value.get_denominator().sqrt();
    if (numerator_root * numerator_root != value.get_numerator() ||
        denominator_root * denominator_root != value.get_denominator()) {
        return false;
    }
    root = Rational(numerator_root, denominator_root);
    return true;
}

ExprPtr sqrt_rational_expression(const Rational& value) {
    Rational root;
    if (exact_rational_sqrt(value, root)) {
        return rational_expression(root);
    }
    // Factor the largest perfect square out of num*den so sqrt(28)
    // canonicalizes to 2*sqrt(7) instead of surviving as sqrt(28).
    BigInt numerator = value.get_numerator();
    BigInt denominator = value.get_denominator();
    BigInt outside_numerator(1);
    BigInt outside_denominator(1);
    for (BigInt factor(2); factor * factor <= numerator; factor = factor + BigInt(1)) {
        while (numerator % (factor * factor) == BigInt(0)) {
            numerator = numerator / (factor * factor);
            outside_numerator = outside_numerator * factor;
        }
    }
    for (BigInt factor(2); factor * factor <= denominator; factor = factor + BigInt(1)) {
        while (denominator % (factor * factor) == BigInt(0)) {
            denominator = denominator / (factor * factor);
            outside_denominator = outside_denominator * factor;
        }
    }
    const Rational inside(numerator, denominator);
    auto radical = SymbolicExpr::sqrt(rational_expression(inside))->simplify();
    if (!radical) {
        return SymbolicExpr::sqrt(rational_expression(value));
    }
    auto scaled = SymbolicExpr::multiply(
        rational_expression(Rational(outside_numerator, outside_denominator)),
        radical)->simplify();
    return scaled ? scaled : radical;
}

ExprResult verified_lsr_complex(ExprPtr real_part, ExprPtr imag_part) {
    auto value = complex(std::move(real_part), std::move(imag_part));
    if (!value) return value;
    auto simplified = value.value()->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return expression_failure(CasErrc::InternalInvariant,
                                  "complex root simplification returned null",
                                  kSolveExprSetOperation);
    }
    return ExprResult::success(std::move(simplified));
}

ExprResult normalize_equation_to_zero(const ExprPtr& expression,
                                      const bool require_equality,
                                      ComputationContext& context) {
    if (!expression) {
        return ExprResult::failure(CasErrc::InvalidArgument,
                                   "equation cannot be null",
                                   kSolveExprSetOperation);
    }
    const auto& node = lamina::detail::node(*expression);
    if (const auto relation =
            std::dynamic_pointer_cast<const RelationalNode>(node)) {
        if (relation->op() == RelationalNode::Op::EQ) {
            auto left = lamina::detail::make_expression_ptr(relation->left());
            auto right = lamina::detail::make_expression_ptr(relation->right());
            if (!left || !right) {
                return ExprResult::failure(
                    CasErrc::InternalInvariant,
                    "equation relation has a null operand",
                    kSolveExprSetOperation);
            }
            auto difference = SymbolicExpr::add(
                std::move(left),
                SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                       std::move(right)));
            if (!difference) {
                return ExprResult::failure(
                    CasErrc::UnsupportedExpression,
                    "equation normalization produced no expression",
                    kSolveExprSetOperation);
            }
            auto simplified = difference->simplify();
            if (!simplified) {
                return ExprResult::failure(
                    CasErrc::InternalInvariant,
                    "equation normalization simplification failed",
                    kSolveExprSetOperation);
            }
            (void)context;
            return ExprResult::success(std::move(simplified));
        }
        if (require_equality) {
            return ExprResult::failure(
                CasErrc::InvalidArgument,
                "solve requires an equality relation",
                kSolveExprSetOperation);
        }
    } else if (require_equality) {
        return ExprResult::failure(CasErrc::InvalidArgument,
                                   "solve requires an equality relation",
                                   kSolveExprSetOperation);
    }
    return ExprResult::success(expression);
}

Result<std::optional<ExprSet>> try_lsr_closed_form_rational_poly_roots(
    const ExprPtr& equation,
    const std::string& variable,
    ComputationContext& context) {
    if (!equation) {
        return Result<std::optional<ExprSet>>::failure(
            CasErrc::InvalidArgument, "equation cannot be null",
            kSolveExprSetOperation);
    }
    if (variable.empty()) {
        return Result<std::optional<ExprSet>>::failure(
            CasErrc::InvalidArgument, "solve variable cannot be empty",
            kSolveExprSetOperation);
    }

    auto recognized = recognize_rational_polynomial(*equation, variable, context);
    if (!recognized) {
        return Result<std::optional<ExprSet>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    const Polynomial<Rational>& polynomial = *recognized.value();
    const int degree = polynomial.degree();
    if (degree < 1 || degree > 2) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    std::vector<ExprPtr> roots;
    if (degree == 1) {
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        if (b.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        roots.push_back(rational_expression((-c) / b)->simplify());
    } else {
        const Rational a = polynomial_coeff_or_zero(polynomial, 2);
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        if (a.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }

        const Rational two_a = Rational(2) * a;
        const Rational discriminant = b * b - Rational(4) * a * c;
        const Rational real_component = (-b) / two_a;

        if (discriminant < Rational(0)) {
            const Rational positive_discriminant = -discriminant;
            const Rational positive_denominator = two_a.abs();
            auto imag_magnitude = SymbolicExpr::divide(
                sqrt_rational_expression(positive_discriminant),
                rational_expression(positive_denominator))->simplify();
            auto positive = verified_lsr_complex(
                rational_expression(real_component), imag_magnitude);
            if (!positive) {
                return Result<std::optional<ExprSet>>::failure(positive.error());
            }
            auto negative_imag = SymbolicExpr::multiply(
                SymbolicExpr::number(-1), imag_magnitude)->simplify();
            auto negative = verified_lsr_complex(
                rational_expression(real_component), negative_imag);
            if (!negative) {
                return Result<std::optional<ExprSet>>::failure(negative.error());
            }
            roots.push_back(negative.value());
            roots.push_back(positive.value());
        } else {
            auto sqrt_discriminant = sqrt_rational_expression(discriminant);
            auto numerator_left = SymbolicExpr::add(
                rational_expression(-b),
                SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                       sqrt_discriminant))->simplify();
            roots.push_back(SymbolicExpr::divide(
                numerator_left, rational_expression(two_a))->simplify());
            if (!discriminant.is_zero()) {
                auto numerator_right = SymbolicExpr::add(
                    rational_expression(-b), sqrt_discriminant)->simplify();
                roots.push_back(SymbolicExpr::divide(
                    numerator_right, rational_expression(two_a))->simplify());
            }
        }
    }

    auto set = ExprSet::make(std::move(roots));
    if (!set) {
        return Result<std::optional<ExprSet>>::failure(set.error());
    }
    return Result<std::optional<ExprSet>>::success(std::move(set.value()));
}

} // namespace

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      ComputationContext& context,
                      const SolveOptions& options) {
    if (!equation) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "equation cannot be null", "lsr.solve_set");
    }
    if (variable.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "solve variable cannot be empty",
                                    "lsr.solve_set");
    }
    return solve_equation(equation, variable, context, options);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      const SolveOptions& options) {
    ComputationContext context;
    return solve_set(equation, variable, context, options);
}

ExprSetResult solve_normalized_expr_set(const ExprPtr& normalized,
                                        const std::string& variable,
                                        ComputationContext& context,
                                        const SolveOptions& options) {
    auto closed_form = try_lsr_closed_form_rational_poly_roots(
        normalized, variable, context);
    if (!closed_form) {
        return ExprSetResult::failure(closed_form.error());
    }
    if (closed_form.value()) {
        return ExprSetResult::success(std::move(*closed_form.value()));
    }

    auto solved = solve_set(normalized, variable, context, options);
    if (!solved) {
        return ExprSetResult::failure(solved.error());
    }

    const auto& solution_set = solved.value();
    if (solution_set.kind() == SolutionSet::Kind::Empty) {
        return expr_set({});
    }
    if (solution_set.kind() != SolutionSet::Kind::Finite) {
        std::string reason = solution_set.reason();
        if (reason.empty()) {
            reason = "solution set is not a finite enumerable set<Expr>";
        }
        return expr_set_failure(CasErrc::Inconclusive, std::move(reason),
                                kSolveExprSetOperation);
    }

    std::vector<ExprPtr> elements;
    elements.reserve(solution_set.finite_solutions().size());
    for (const auto& solution : solution_set.finite_solutions()) {
        if (!solution.conditions.empty()) {
            return expr_set_failure(CasErrc::Inconclusive,
                                    "conditional finite solutions cannot be lowered to set<Expr>",
                                    kSolveExprSetOperation);
        }
        elements.push_back(solution.value);
    }
    return expr_set(std::move(elements));
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             ComputationContext& context,
                             const SolveOptions& options) {
    auto normalized = normalize_equation_to_zero(equation, true, context);
    if (!normalized) {
        return ExprSetResult::failure(normalized.error());
    }
    return solve_normalized_expr_set(normalized.value(), variable, context,
                                     options);
}
ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             const SolveOptions& options) {
    ComputationContext context;
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    auto normalized = normalize_equation_to_zero(expression, false, context);
    if (!normalized) {
        return ExprSetResult::failure(normalized.error());
    }
    return solve_normalized_expr_set(normalized.value(), variable, context,
                                     options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return roots(expression, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return solve(equation, variable, context, options);
}

} // namespace lamina::lsr
