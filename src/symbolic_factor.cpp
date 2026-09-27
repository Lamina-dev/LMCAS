#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <optional>
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/multivariate_conversion.hpp"
#include "transcendental_factor.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "multivariate_factor.hpp"

namespace LMCAS {

namespace {


/// 多元因式分解辅助函数

static std::optional<int> polynomial_exponent(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) { return std::nullopt; }
    BigInt exponent;
    if (std::holds_alternative<BigInt>(number->value())) {
        exponent = std::get<BigInt>(number->value());
    } else if (std::holds_alternative<Rational>(number->value())) {
        const auto& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) { return std::nullopt; }
        exponent = rational.to_bigint();
    } else {
        const auto value = std::get<lmmc_real_t>(number->value());
        if (!std::isfinite(value) || value != std::floor(value)) { return std::nullopt; }
        exponent = Rational::from_double(value).to_bigint();
    }
    if (exponent.is_negative() || exponent > BigInt(100)) { return std::nullopt; }
    return static_cast<int>(exponent.try_to_int64().value());
}

static bool is_poly_expr_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) { return false; }
    if (std::dynamic_pointer_cast<const NumberNode>(node)) { return true; }
    if (std::dynamic_pointer_cast<const VariableNode>(node)) { return true; }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (!is_poly_expr_node(op)) { return false; }
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (!is_poly_expr_node(op)) { return false; }
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return polynomial_exponent(pow->exponent()).has_value() &&
            is_poly_expr_node(pow->base());
    }

    return false;
}

}

static Rational factor_number_coefficient(const NumberNode& number) {
    if (const auto exact = detail::exact_rational_value(number)) {
        return *exact;
    }
    return Rational::from_double(std::get<lmmc_real_t>(number.value()));
}


static LMCAS::MultiPoly symbolic_node_to_multipoly(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::vector<std::string>& vars)
{
    if (!node) { return LMCAS::MultiPoly(Rational(0), vars); }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return LMCAS::MultiPoly(factor_number_coefficient(*num), vars);
    }

    if (auto var_node = std::dynamic_pointer_cast<const VariableNode>(node)) {
        LMCAS::Monomial mono(vars.size(), 0);
        for (size_t i = 0; i < vars.size(); ++i) {
            if (vars[i] == var_node->name()) {
                mono[i] = 1;
                break;
            }
        }
        std::vector<LMCAS::MultiPoly::Term> terms;
        terms.push_back({mono, Rational(1)});
        return LMCAS::MultiPoly(std::move(terms), vars);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        LMCAS::MultiPoly result(Rational(0), vars);
        for (const auto& op : add->operands()) {
            result = result + symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        LMCAS::MultiPoly result(Rational(1), vars);
        for (const auto& op : mul->operands()) {
            result = result * symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base_poly = symbolic_node_to_multipoly(pow->base(), vars);
        const int exp_val = polynomial_exponent(pow->exponent()).value();
        if (exp_val == 0) { return LMCAS::MultiPoly(Rational(1), vars); }
        LMCAS::MultiPoly result(Rational(1), vars);
        for (int i = 0; i < exp_val; ++i) {
            result = result * base_poly;
        }
        return result;
    }

    return LMCAS::MultiPoly(Rational(0), vars);
}

static std::shared_ptr<SymbolicExpr> multipoly_to_symbolic(const LMCAS::MultiPoly& poly) {
    if (poly.is_zero()) { return SymbolicExpr::number(0); }

    const auto& vars = poly.variables();
    const auto& terms = poly.terms();

    std::vector<std::shared_ptr<SymbolicExpr>> term_exprs;
    term_exprs.reserve(terms.size());

    for (const auto& [mono, coeff] : terms) {
        term_exprs.push_back(
            detail::multipoly_term_to_expression(mono, coeff, vars));
    }

    if (term_exprs.empty()) { return SymbolicExpr::number(0); }
    if (term_exprs.size() == 1) { return term_exprs[0]->simplify(); }

    auto result = term_exprs[0];
    for (size_t i = 1; i < term_exprs.size(); ++i) {
        result = SymbolicExpr::add(result, term_exprs[i]);
    }
    return result->simplify();
}

static std::optional<ExpressionResult> factor_numeric_content(
    const std::shared_ptr<SymbolicExpr>& simp,
    const std::shared_ptr<const AddNode>& add_node, ComputationContext& context) {
    auto content_result =
        LMCAS::symbolic_polynomial_content(*simp, context);
    if (!content_result) {
        if (content_result.error().code == CasErrc::UnsupportedExpression) {
            return std::nullopt;
        }
        return LMCAS::ExpressionResult::failure(
            content_result.error());
    }
    const Rational content_value = content_result.value();
    if (content_value != Rational(0) && content_value != Rational(1)) {
        auto content = content_value.is_integer()
            ? SymbolicExpr::number(content_value.to_bigint())
            : SymbolicExpr::number(content_value);
        const Rational inverse_value = Rational(1) / content_value;
        auto inverse_content = inverse_value.is_integer()
            ? SymbolicExpr::number(inverse_value.to_bigint())
            : SymbolicExpr::number(inverse_value);
        std::vector<std::shared_ptr<const SymbolicNode>> primitive_terms;
        primitive_terms.reserve(add_node->operands().size());
        for (const auto& operand : add_node->operands()) {
            auto term = LMCAS::detail::make_expression_ptr(operand);
            auto primitive_term = SymbolicExpr::multiply(term, inverse_content)->simplify();
            primitive_terms.push_back(LMCAS::detail::node(primitive_term));
        }
        auto primitive_sum = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<AddNode>(primitive_terms));
        auto factored_primitive =
            primitive_sum->factor_checked(context);
        if (!factored_primitive) { return factored_primitive; }
        return SymbolicExpr::multiply(
            content, factored_primitive.value())->simplify();
    }
    return std::nullopt;
}

static std::optional<ExpressionResult> factor_common_polynomial(
    const std::shared_ptr<const AddNode>& add_node, ComputationContext& context) {
    std::shared_ptr<SymbolicExpr> common = nullptr;
    for (const auto& op : add_node->operands()) {
         auto expr_op = LMCAS::detail::make_expression_ptr(op);
         if (!common) {
             common = expr_op;
         } else {
             auto gcd = LMCAS::symbolic_polynomial_gcd(*common, *expr_op, context);
             if (!gcd) {
                 if (gcd.error().code == CasErrc::UnsupportedExpression) {
                     return std::nullopt;
                 }
                 return gcd;
             }
             common = std::move(gcd.value());
         }
    }

    if (common && !common->is_one() && !common->is_zero()) {
         std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
         for (const auto& op : add_node->operands()) {
              auto term = LMCAS::detail::make_expression_ptr(op);

              auto inv_common = SymbolicExpr::power(common, SymbolicExpr::number(-1));
              auto quot = SymbolicExpr::multiply(term, inv_common);
              quot = quot->simplify();
              new_ops.push_back(LMCAS::detail::node(quot));
         }
         auto new_sum = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(new_ops));

         auto factored_sum = new_sum->factor_checked(context);
         if (!factored_sum) { return factored_sum; }

         std::vector<std::shared_ptr<const SymbolicNode>> final_ops = {
             LMCAS::detail::node(common),
             LMCAS::detail::node(factored_sum.value())};
         return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(final_ops));
    }
    return std::nullopt;
}

static std::vector<std::shared_ptr<const SymbolicNode>> rational_root_factors(
    const std::vector<Rational>& roots, const Rational& leading,
    const std::shared_ptr<SymbolicExpr>& x_node) {
    std::vector<std::shared_ptr<const SymbolicNode>> factors;

    if (!(leading == Rational(1))) {
        factors.push_back(LMCAS::detail::node(SymbolicExpr::number(leading)));
    }

    for (const auto& r : roots) {
        std::shared_ptr<SymbolicExpr> linear_factor;
        if (r == Rational(0)) {
            linear_factor = x_node;
        } else {
            auto neg_r = SymbolicExpr::number(Rational(0) - r);
            linear_factor = SymbolicExpr::add(x_node, neg_r)->simplify();
        }
        factors.push_back(LMCAS::detail::node(linear_factor));
    }
    return factors;
}

static Polynomial<Rational> rational_root_product(
    const std::vector<Rational>& roots, const Rational& leading,
    const std::string& var) {
    LMCAS::Polynomial<Rational> factored_product({leading}, var);
    for (const auto& r : roots) {
        LMCAS::Polynomial<Rational> lin({Rational(0) - r, Rational(1)}, var);
        factored_product = factored_product * lin;
    }
    return factored_product;
}

struct RationalRootReconstruction {
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    Polynomial<Rational> quotient;
    Polynomial<Rational> remainder;
};

static RationalRootReconstruction reconstruct_rational_roots(
    const Polynomial<Rational>& polynomial,
    const std::vector<Rational>& roots, const std::string& variable) {
    auto symbol =
        detail::make_expression_ptr(detail::make_node<VariableNode>(variable));
    const auto leading = polynomial.lead_coeff();
    auto factors = rational_root_factors(roots, leading, symbol);
    auto product = rational_root_product(roots, leading, variable);
    auto [quotient, remainder] = polynomial.div_mod(product);
    return {
        std::move(factors), std::move(quotient), std::move(remainder)};
}

static std::optional<ExpressionResult> factor_rational_univariate(
    const std::shared_ptr<SymbolicExpr>& simp, const std::string& var,
    ComputationContext& context) {
    try {
        auto converted = symbolic_to_poly<Rational>(simp, var);
        if (!converted) {
            if (converted.error().code == CasErrc::UnsupportedExpression) return std::nullopt;
            return ExpressionResult::failure(converted.error());
        }
        const auto& poly = converted.value();
        if (poly.degree() < 2) { return std::nullopt; }
        auto roots = find_rational_roots(poly);
        if (roots.empty()) { return std::nullopt; }
        auto reconstruction =
            reconstruct_rational_roots(poly, roots, var);
        if (!reconstruction.remainder.is_zero()) return std::nullopt;
        auto& factors = reconstruction.factors;
        auto& quotient = reconstruction.quotient;
        if (quotient.degree() >= 1) {
            auto q_expr = poly_to_symbolic(quotient)->simplify();
            if (quotient.degree() >= 2) {
                auto q_factored = q_expr->factor_checked(context);
                if (!q_factored) { return q_factored; }
                factors.push_back(detail::node(q_factored.value()));
            } else {
                factors.push_back(detail::node(q_expr));
            }
        }
        if (factors.size() > 1) {
            return detail::make_expression_ptr(detail::make_node<MultiplyNode>(factors));
        }
        if (factors.size() == 1) { return detail::make_expression_ptr(factors[0]); }
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
    return std::nullopt;
}

static std::optional<ExpressionResult> factor_symbolic_quadratic(
    const std::shared_ptr<SymbolicExpr>& simp, const std::string& var,
    ComputationContext& context) {
    try {
        auto converted = LMCAS::symbolic_to_poly<LMCAS::SymbolicPolyCoeff>(simp, var);
        if (!converted) {
            if (converted.error().code == CasErrc::UnsupportedExpression) return std::nullopt;
            return ExpressionResult::failure(converted.error());
        }
        const auto& poly = converted.value();
        if (poly.degree() == 2) {
              auto solved =
                 LMCAS::solve_finite_checked(simp, var, context);
             if (!solved) {
                 return LMCAS::ExpressionResult::failure(
                     solved.error());
             }
             auto solutions = std::move(solved.value());
             if (solutions.size() == 2) {
                  auto leading = poly.coeffs[2].val;
                  if (!leading) { leading = SymbolicExpr::number(1); }
                  leading = leading->simplify();

                  auto x_node = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<VariableNode>(var));

                  auto t1 = SymbolicExpr::add(x_node, SymbolicExpr::multiply(solutions[0], SymbolicExpr::number(-1)))->simplify();
                  auto t2 = SymbolicExpr::add(x_node, SymbolicExpr::multiply(solutions[1], SymbolicExpr::number(-1)))->simplify();

                  std::vector<std::shared_ptr<const SymbolicNode>> factors;
                  if (!leading->is_one()) { factors.push_back(LMCAS::detail::node(leading)); }
                  factors.push_back(LMCAS::detail::node(t1));
                  factors.push_back(LMCAS::detail::node(t2));

                  return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(factors));
             }
        }
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
    return std::nullopt;
}

static std::shared_ptr<SymbolicExpr> multivariate_factor_expression(
    const MultiFactorResult& result) {
    std::vector<std::shared_ptr<const SymbolicNode>> factor_nodes;

    if (!(result.constant == Rational(1))) {
        factor_nodes.push_back(
            LMCAS::detail::node(SymbolicExpr::number(result.constant)));
    }

    for (size_t i = 0; i < result.factors.size(); ++i) {
        auto factor_expr = multipoly_to_symbolic(result.factors[i]);
        if (!factor_expr || factor_expr->is_zero()) { continue; }

        int mult = (i < result.multiplicities.size())
            ? result.multiplicities[i] : 1;

        if (mult == 1) {
            factor_nodes.push_back(LMCAS::detail::node(factor_expr));
        } else {
            auto pow_expr = SymbolicExpr::power(
                factor_expr, SymbolicExpr::number(mult));
            factor_nodes.push_back(LMCAS::detail::node(pow_expr));
        }
    }

    if (factor_nodes.size() > 1) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<MultiplyNode>(std::move(factor_nodes)));
    } else if (factor_nodes.size() == 1) {
        return LMCAS::detail::make_expression_ptr(factor_nodes[0]);
    }
    return nullptr;
}

static std::optional<ExpressionResult> factor_multivariate_expression(
    const std::shared_ptr<SymbolicExpr>& simp,
    const std::set<std::string>& factor_variables, ComputationContext& context) {
    if (!is_poly_expr_node(detail::node(simp))) { return std::nullopt; }
    try {
        std::vector<std::string> var_list(
            factor_variables.begin(), factor_variables.end());
        auto mpoly = symbolic_node_to_multipoly(LMCAS::detail::node(simp), var_list);

        if (!mpoly.is_zero() && !mpoly.is_constant()) {
            auto factorization =
                LMCAS::factor_multivariate_checked(
                    mpoly, context);
            if (!factorization) {
                return LMCAS::ExpressionResult::failure(
                    factorization.error());
            }
            const auto& checked = factorization.value();
            if (checked.completeness !=
                LMCAS::Completeness::Complete) {
                return simp;
            }
            const auto& result = checked.value;

            if (!result.factors.empty()) {
                auto expression = multivariate_factor_expression(result);
                if (expression) { return expression; }
            }
        }
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
    return std::nullopt;
}

static std::optional<ExpressionResult> reconstruct_rational_root_factors(
    const Polynomial<Rational>& polynomial, const std::vector<Rational>& roots,
    const std::string& var) {
    auto reconstruction =
        reconstruct_rational_roots(polynomial, roots, var);
    auto& factors = reconstruction.factors;
    auto& quotient = reconstruction.quotient;
    const auto& remainder = reconstruction.remainder;
    if (remainder.is_zero() && quotient.degree() >= 1) {
        auto q_expr = LMCAS::poly_to_symbolic(quotient)->simplify();
        factors.push_back(LMCAS::detail::node(q_expr));
    } else if (remainder.is_zero() && !quotient.is_zero() && quotient.degree() == 0) {
        if (!(quotient.coeffs[0] == Rational(1))) {
            factors.push_back(
                LMCAS::detail::node(SymbolicExpr::number(quotient.coeffs[0])));
        }
    }
    if (factors.size() > 1) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<MultiplyNode>(factors));
    }
    return std::nullopt;
}

static std::optional<ExpressionResult> factor_variable_rational_roots(
    const std::shared_ptr<SymbolicExpr>& simp, const std::string& var) {
    try {
        auto poly = LMCAS::symbolic_to_poly<LMCAS::SymbolicPolyCoeff>(simp, var);
        if (!poly) {
            if (poly.error().code == CasErrc::UnsupportedExpression) { return std::nullopt; }
            return ExpressionResult::failure(poly.error());
        }
        if (poly.value().degree() >= 2) {
            auto converted = LMCAS::symbolic_to_poly<Rational>(simp, var);
            if (!converted) {
                if (converted.error().code == CasErrc::UnsupportedExpression) { return std::nullopt; }
                return ExpressionResult::failure(converted.error());
            }
            const auto& poly_r = converted.value();
            auto roots = LMCAS::find_rational_roots(poly_r);
            if (!roots.empty()) {
                return reconstruct_rational_root_factors(poly_r, roots, var);
            }
        }
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
    return std::nullopt;
}

static ExpressionResult factor_transcendental_expression(
    const std::shared_ptr<SymbolicExpr>& simp,
    ComputationContext& context) {
    const auto transform_variables = LMCAS::free_variables(LMCAS::detail::node(simp));
    if (transform_variables.empty()) { return simp; }

    auto converted = LMCAS::factor_transcendental(
        simp, *transform_variables.begin(), context);
    if (!converted) { return ExpressionResult::failure(converted.error()); }
    const auto& trans_factors = converted.value();
    if (trans_factors.size() > 1) {
        std::vector<std::shared_ptr<const SymbolicNode>> factor_nodes;
        factor_nodes.reserve(trans_factors.size());
        for (const auto& factor : trans_factors) {
            if (factor && LMCAS::detail::node(factor)) {
                factor_nodes.push_back(LMCAS::detail::node(factor));
            }
        }
        if (factor_nodes.size() > 1) {
            return LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<MultiplyNode>(std::move(factor_nodes)));
        }
    }
    return simp;
}

LMCAS::ExpressionResult SymbolicExpr::factor_impl(
    LMCAS::ComputationContext& context) const {
    auto simp = simplify();
    if (!simp || !detail::node(simp)) { return simp; }
    auto add_node = std::dynamic_pointer_cast<const AddNode>(detail::node(simp));
    if (!add_node) { return factor_transcendental_expression(simp, context); }

    auto content = factor_numeric_content(simp, add_node, context);
    if (content) { return std::move(*content); }
    auto common = factor_common_polynomial(add_node, context);
    if (common) { return std::move(*common); }

    auto variables = free_variables(detail::node(simp));
    if (variables.size() == 1) {
        const auto& variable = *variables.begin();
        auto rational = factor_rational_univariate(simp, variable, context);
        if (rational) { return std::move(*rational); }
        auto quadratic = factor_symbolic_quadratic(simp, variable, context);
        if (quadratic) { return std::move(*quadratic); }
    }
    if (variables.size() > 1) {
        auto multivariate = factor_multivariate_expression(simp, variables, context);
        if (multivariate) { return std::move(*multivariate); }
        for (const auto& variable : variables) {
            auto factored = factor_variable_rational_roots(simp, variable);
            if (factored) { return std::move(*factored); }
        }
    }
    return factor_transcendental_expression(simp, context);
}

LMCAS::ExpressionResult SymbolicExpr::factor_checked(
    LMCAS::ComputationContext& context) const
{
    constexpr const char* operation = "factor";
    auto budget = context.consume_steps(1, operation);
    if (!budget) { return LMCAS::ExpressionResult::failure(budget.error()); }
    try {
        return factor_impl(context);
    } catch (const std::bad_alloc&) {
        return LMCAS::ExpressionResult::failure(
            LMCAS::CasErrc::ResourceLimit,
            "allocation failed while factoring expression", operation);
    } catch (const std::exception& ex) {
        return LMCAS::ExpressionResult::failure(
            LMCAS::CasErrc::InternalInvariant, ex.what(), operation);
    }
}

LMCAS::ExpressionResult SymbolicExpr::factor_checked() const
{
    LMCAS::ComputationContext context;
    return factor_checked(context);
}


}
