#include "root_of_utils.hpp"
#include "newton_raphson.hpp"
#include "numeric_evaluation.hpp"
#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace lamina {

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

std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    int n = poly.degree();
    if (n <= 0) return {};

    auto poly_expr = symbolic_poly_to_expr(poly);
    if (auto canonical = poly_expr->simplify()) {
        if (lamina::detail::node(canonical)) {
            poly_expr = canonical;
        }
    }

    std::vector<std::shared_ptr<SymbolicExpr>> solutions;
    solutions.reserve(n);

    for (int k = 0; k < n; ++k) {
        solutions.push_back(SymbolicExpr::root_of(poly_expr, var, k));
    }

    return solutions;
}

std::optional<lmmc_real_t> rootof_evaluate(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    ComputationContext context;
    auto result = rootof_evaluate_checked(rootof_expr, context);
    return result ? std::optional<lmmc_real_t>(result.value()) : std::nullopt;
}

RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context)
{
    constexpr const char* operation = "rootof_evaluate";
    if (!rootof_expr || !lamina::detail::node(rootof_expr)) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf expression cannot be null", operation);
    }

    auto func_node = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(rootof_expr));
    if (!func_node || func_node->type() != FunctionNode::FuncType::RootOf) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "expression is not a RootOf", operation);
    }

    if (func_node->arguments().size() != 3) {
        return RootOfEvaluationResult::failure(
            CasErrc::InternalInvariant, "RootOf must contain three arguments", operation);
    }

    auto poly_node = func_node->arguments()[0];
    auto var_node = func_node->arguments()[1];
    auto index_node = func_node->arguments()[2];

    auto var_sym = std::dynamic_pointer_cast<const VariableNode>(var_node);
    if (!var_sym || var_sym->name().empty()) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf variable must be a named symbol", operation);
    }
    std::string var = var_sym->name();

    auto idx_num = std::dynamic_pointer_cast<const NumberNode>(index_node);
    if (!idx_num || std::holds_alternative<lmmc_real_t>(idx_num->value())) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf index must be an exact integer", operation);
    }

    BigInt index;
    if (std::holds_alternative<BigInt>(idx_num->value())) {
        index = std::get<BigInt>(idx_num->value());
    } else {
        const Rational& rational = std::get<Rational>(idx_num->value());
        if (!rational.is_integer()) {
            return RootOfEvaluationResult::failure(
                CasErrc::InvalidArgument, "RootOf index must be an integer", operation);
        }
        index = rational.to_BigInt();
    }
    auto index_value = index.try_to_int64();
    if (!index_value || *index_value < 0) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf index is negative or too large", operation);
    }

    auto poly_expr = lamina::detail::make_expression_ptr(poly_node);
    auto recognized = recognize_rational_polynomial(*poly_expr, var, context);
    if (!recognized) return RootOfEvaluationResult::failure(recognized.error());
    if (!recognized.value()) {
        return RootOfEvaluationResult::failure(
            CasErrc::Inconclusive,
            "RootOf numeric evaluation supports exact rational polynomials only",
            operation);
    }
    const Polynomial<Rational>& rat_poly = *recognized.value();
    const int degree = rat_poly.degree();
    if (degree <= 0) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf polynomial must have positive degree", operation);
    }
    if (*index_value >= degree) {
        return RootOfEvaluationResult::failure(
            CasErrc::InvalidArgument, "RootOf index exceeds polynomial degree", operation);
    }

    auto isolation_step = context.consume_steps(
        static_cast<std::size_t>(degree), operation);
    if (!isolation_step) return RootOfEvaluationResult::failure(isolation_step.error());
    auto intervals = isolate_real_roots(rat_poly);
    if (static_cast<int>(intervals.size()) != degree) {
        return RootOfEvaluationResult::failure(
            CasErrc::Inconclusive,
            "RootOf double evaluation currently requires every root to be real",
            operation);
    }

    const auto& interval = intervals[static_cast<std::size_t>(*index_value)];
    const lmmc_real_t lo = interval.first.to_double();
    const lmmc_real_t hi = interval.second.to_double();
    const lmmc_real_t x0 = (lo + hi) / 2.0;
    auto df_expr = poly_expr->differentiate(var);

    SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    auto refined = newton_raphson_checked(
        poly_expr, df_expr, var, x0, lo, hi, context, opts);
    if (!refined) return RootOfEvaluationResult::failure(refined.error());
    if (!refined.value()) {
        return RootOfEvaluationResult::failure(
            CasErrc::Inconclusive, "isolated RootOf could not be numerically verified", operation);
    }
    return RootOfEvaluationResult::success(refined.value()->value);
}

static std::pair<std::complex<double>, bool> try_eval_complex(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !lamina::detail::node(expr)) return {{0.0, 0.0}, false};

    ComputationContext context;
    auto evaluated = evaluate_numeric(*expr, {}, context);
    if (evaluated && evaluated.value().is_finite()) {
        return {{evaluated.value().value, 0.0}, true};
    }

    return {{0.0, 0.0}, false};
}

struct RootWithIndex {
    std::shared_ptr<SymbolicExpr> expr;
    double real_part;
    double imag_part;
    bool is_real;
    bool eval_success;
};

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

static std::vector<std::shared_ptr<SymbolicExpr>> sort_roots_by_convention(
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots)
{
    if (roots.empty()) return roots;

    std::vector<RootWithIndex> evaluated;
    evaluated.reserve(roots.size());

    bool all_evaluated = true;
    for (const auto& root : roots) {
        RootWithIndex ri;
        ri.expr = root;
        auto [val, success] = try_eval_complex(root);
        ri.eval_success = success;
        if (success) {
            ri.real_part = val.real();
            ri.imag_part = val.imag();

            ri.is_real = (std::abs(ri.imag_part) < 1e-10);
        } else {
            all_evaluated = false;
            ri.real_part = 0.0;
            ri.imag_part = 0.0;
            ri.is_real = true;
        }
        evaluated.push_back(ri);
    }

    if (!all_evaluated) {
        return roots;
    }

    std::vector<RootWithIndex> real_roots;
    std::vector<RootWithIndex> complex_roots;

    for (auto& ri : evaluated) {
        if (ri.is_real) {
            real_roots.push_back(ri);
        } else {
            complex_roots.push_back(ri);
        }
    }

    std::sort(real_roots.begin(), real_roots.end(),
        [](const RootWithIndex& a, const RootWithIndex& b) {
            return a.real_part < b.real_part;
        });

    std::sort(complex_roots.begin(), complex_roots.end(),
        [](const RootWithIndex& a, const RootWithIndex& b) {
            if (std::abs(a.real_part - b.real_part) > 1e-10) {
                return a.real_part < b.real_part;
            }

            return a.imag_part > b.imag_part;
        });

    std::vector<std::shared_ptr<SymbolicExpr>> sorted;
    sorted.reserve(roots.size());
    for (const auto& ri : real_roots) {
        sorted.push_back(ri.expr);
    }
    for (const auto& ri : complex_roots) {
        sorted.push_back(ri.expr);
    }

    return sorted;
}

static bool is_structural_polynomial(const std::shared_ptr<const SymbolicNode>& node,
                                     const std::string& variable) {
    if (!node) return false;
    if (!expression_depends_on_variable(node, variable)) return true;

    if (auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return symbol->name() == variable;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return std::all_of(add->operands().begin(), add->operands().end(),
            [&](const auto& operand) {
                return is_structural_polynomial(operand, variable);
            });
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return std::all_of(multiply->operands().begin(), multiply->operands().end(),
            [&](const auto& operand) {
                return is_structural_polynomial(operand, variable);
            });
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent_node = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        if (!exponent_node || std::holds_alternative<lmmc_real_t>(exponent_node->value())) {
            return false;
        }
        BigInt exponent;
        if (std::holds_alternative<BigInt>(exponent_node->value())) {
            exponent = std::get<BigInt>(exponent_node->value());
        } else {
            const Rational& rational = std::get<Rational>(exponent_node->value());
            if (!rational.is_integer()) return false;
            exponent = rational.to_BigInt();
        }
        auto value = exponent.try_to_int64();
        return value && *value > 0 && *value < 1000 &&
               is_structural_polynomial(power->base(), variable);
    }
    return false;
}

static std::optional<int> exact_root_index(const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number || std::holds_alternative<lmmc_real_t>(number->value())) {
        return std::nullopt;
    }
    BigInt index;
    if (std::holds_alternative<BigInt>(number->value())) {
        index = std::get<BigInt>(number->value());
    } else {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) return std::nullopt;
        index = rational.to_BigInt();
    }
    auto value = index.try_to_int64();
    if (!value || *value < 0 || *value > std::numeric_limits<int>::max()) {
        return std::nullopt;
    }
    return static_cast<int>(*value);
}

std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    if (!rootof_expr || !lamina::detail::node(rootof_expr)) return rootof_expr;

    auto func_node = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(rootof_expr));
    if (!func_node || func_node->type() != FunctionNode::FuncType::RootOf) {
        return rootof_expr;
    }

    if (func_node->arguments().size() != 3) {
        return rootof_expr;
    }

    auto poly_node = func_node->arguments()[0];
    auto var_node = func_node->arguments()[1];
    auto index_node = func_node->arguments()[2];

    auto var_sym = std::dynamic_pointer_cast<const VariableNode>(var_node);
    if (!var_sym) return rootof_expr;
    std::string var = var_sym->name();

    auto index = exact_root_index(index_node);
    if (!index) return rootof_expr;
    const int k = *index;

    if (!is_structural_polynomial(poly_node, var)) return rootof_expr;

    auto poly_expr = lamina::detail::make_expression_ptr(poly_node);
    auto sym_poly = symbolic_to_poly<SymbolicPolyCoeff>(poly_expr, var);

    int degree = sym_poly.degree();
    if (degree <= 0) return rootof_expr;

    if (k < 0 || k >= degree) {
        return rootof_expr;
    }

    if (degree <= 4) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        if (roots.empty()) {
            return rootof_expr;
        }

        auto sorted_roots = sort_roots_by_convention(roots);

        if (k >= 0 && k < static_cast<int>(sorted_roots.size())) {
            return sorted_roots[k];
        }

        return rootof_expr;
    }

    return rootof_expr;
}

}
