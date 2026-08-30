#include "../include/symbolic.hpp"
#include "../include/parametric_solver.hpp"
#include "../include/solver.hpp"
#include "../include/symbolic_matrix.hpp"
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/poly_utils.hpp"
#include "../include/solve_strategies.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_transcendental.hpp"
#include "../include/solve_mixed_transcendental.hpp"
#include "../include/newton_raphson.hpp"
#include "../include/root_of_utils.hpp"

using namespace lamina;
using SolveVectorResult =
    lamina::Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

std::shared_ptr<SymbolicExpr> get_coeff(const Polynomial<SymbolicPolyCoeff>& p, int deg) {
    if (deg < 0 || deg > p.degree()) return SymbolicExpr::number(0);
    return p.coeffs[deg].val;
}


static bool solve_contains_rootof_node(
    const std::shared_ptr<const SymbolicNode>& node) {
    return lamina::detail::contains_node_type<RootOfNode>(node);
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_filter_rootof_results(
    std::vector<std::shared_ptr<SymbolicExpr>> results,
    const SolveOptions& opts) {
    if (opts.return_rootof) return results;

    std::vector<std::shared_ptr<SymbolicExpr>> filtered;
    filtered.reserve(results.size());
    for (auto& result : results) {
        if (!result || !solve_contains_rootof_node(lamina::detail::node(result))) {
            filtered.push_back(std::move(result));
        }
    }
    return filtered;
}


static std::vector<std::shared_ptr<SymbolicExpr>> solve_polynomial_closed(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    int deg = poly.degree();
    if (deg <= 0) return {};

    if (deg == 1) {
        auto a = get_coeff(poly, 1);
        auto b = get_coeff(poly, 0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto result = SymbolicExpr::divide(neg_b, a);
        return { result->simplify() };
    }

    if (deg == 2) {
        auto a = get_coeff(poly, 2);
        auto b = get_coeff(poly, 1);
        auto c = get_coeff(poly, 0);

        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto ac4 = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(b2, SymbolicExpr::multiply(ac4, SymbolicExpr::number(-1)));

        auto sqrt_delta = SymbolicExpr::power(delta, SymbolicExpr::number(0.5));
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);

        auto numerator1 = SymbolicExpr::add(neg_b, sqrt_delta);
        auto numerator2 = SymbolicExpr::add(neg_b, SymbolicExpr::multiply(sqrt_delta, SymbolicExpr::number(-1)));

        auto x1 = SymbolicExpr::divide(numerator1, two_a);
        auto x2 = SymbolicExpr::divide(numerator2, two_a);

        return { x1->simplify(), x2->simplify() };
    }

    if (deg == 3) {
        auto a = get_coeff(poly, 3);
        auto b = get_coeff(poly, 2);
        auto c = get_coeff(poly, 1);
        auto d = get_coeff(poly, 0);
        return solve_cubic(a, b, c, d, var);
    }

    if (deg == 4) {
        auto a = get_coeff(poly, 4);
        auto b = get_coeff(poly, 3);
        auto c = get_coeff(poly, 2);
        auto d = get_coeff(poly, 1);
        auto e = get_coeff(poly, 0);
        return solve_quartic(a, b, c, d, e, var);
    }

    return {};
}


static std::optional<SolutionSet> solve_direct_function_family(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(
        lamina::detail::node(expression));
    if (!function || function->arguments().size() != 1) {
        return std::nullopt;
    }
    auto argument = std::dynamic_pointer_cast<const VariableNode>(
        function->arguments()[0]);
    if (!argument || argument->name() != variable) return std::nullopt;

    auto integer = SymbolicExpr::variable("_k");
    auto pi = SymbolicExpr::variable("pi");
    auto integer_pi = SymbolicExpr::multiply(integer, pi);
    switch (function->type()) {
        case FunctionNode::FuncType::Sin:
        case FunctionNode::FuncType::Tan:
            return SolutionSet{ParametricSolutions{{
                ParametricSolution{integer_pi, {"_k"}, {}}}}};
        case FunctionNode::FuncType::Cos:
            return SolutionSet{ParametricSolutions{{
                ParametricSolution{
                    SymbolicExpr::add(
                        SymbolicExpr::multiply(
                            SymbolicExpr::number(Rational(1, 2)), pi),
                        integer_pi),
                    {"_k"}, {}}}}};
        case FunctionNode::FuncType::Exp:
            return SolutionSet{EmptySolutions{}};
        case FunctionNode::FuncType::Ln:
            return SolutionSet{FiniteSolutions{{
                FiniteSolution{SymbolicExpr::number(1), 1, {}}}}};
        default:
            return std::nullopt;
    }
}
static SolveVectorResult solve_finite_vector_core(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts);

lamina::SolveResult lamina::solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_equation";
    if (!expr || !lamina::detail::node(expr)) {
        return SolveResult::failure(
            CasErrc::InvalidArgument, "equation expression cannot be null", operation);
    }
    if (var.empty()) {
        return SolveResult::failure(
            CasErrc::InvalidArgument, "solve variable cannot be empty", operation);
    }

    auto step = context.consume_steps(1, operation);
    if (!step) return SolveResult::failure(step.error());

    std::shared_ptr<SymbolicExpr> f_expr = expr;
    if (auto relation = std::dynamic_pointer_cast<const RelationalNode>(lamina::detail::node(expr))) {
        if (relation->op() != RelationalNode::Op::EQ) {
            return SolveResult::failure(
                CasErrc::InvalidArgument,
                "solve_equation accepts equations, not inequalities",
                operation);
        }
        auto left = lamina::detail::make_expression_ptr(relation->left());
        auto right = lamina::detail::make_expression_ptr(relation->right());
        f_expr = SymbolicExpr::add(
            left, SymbolicExpr::multiply(SymbolicExpr::number(-1), right));
    }

    auto polynomial = recognize_rational_polynomial(*f_expr, var, context);
    if (!polynomial) return SolveResult::failure(polynomial.error());
    if (polynomial.value()) {
        const auto& exact = *polynomial.value();
        if (exact.is_zero()) {
            return SolveResult::success(SolutionSet{UniversalSolutions{}});
        }
        if (exact.degree() == 0) {
            return SolveResult::success(SolutionSet{EmptySolutions{}});
        }

        std::vector<SymbolicPolyCoeff> coefficients;
        coefficients.reserve(exact.coeffs.size());
        for (const auto& coefficient : exact.coeffs) {
            coefficients.emplace_back(SymbolicExpr::number(coefficient));
        }
        Polynomial<SymbolicPolyCoeff> symbolic_exact(coefficients, var);

        std::vector<std::shared_ptr<SymbolicExpr>> roots;
        if (exact.degree() <= 4) {
            roots = solve_polynomial_closed(symbolic_exact, var);
        } else {
            roots = solve_by_factoring(symbolic_exact, var);
            roots = solve_filter_rootof_results(std::move(roots), opts);
            if (roots.empty() && opts.return_rootof) {
                roots = make_rootof_solutions(symbolic_exact, var);
            }
        }
        if (roots.empty()) {
            return SolveResult::failure(
                CasErrc::Inconclusive,
                "exact polynomial roots could not be represented", operation);
        }

        std::vector<FiniteSolution> solutions;
        for (auto& root : roots) {
            auto root_step = context.consume_steps(1, operation);
            if (!root_step) return SolveResult::failure(root_step.error());
            if (!root || !lamina::detail::node(root)) {
                return SolveResult::failure(
                    CasErrc::InternalInvariant,
                    "polynomial solver produced a null root", operation);
            }
            auto existing = std::find_if(
                solutions.begin(), solutions.end(),
                [&](const FiniteSolution& solution) {
                    return lamina::detail::node(solution.value)->equals(
                        *lamina::detail::node(root));
                });
            if (existing != solutions.end()) {
                ++existing->multiplicity;
            } else {
                solutions.push_back(
                    FiniteSolution{std::move(root), 1, {}});
            }
        }
        return SolveResult::success(
            SolutionSet{FiniteSolutions{std::move(solutions)}});
    }

    if (auto direct = solve_direct_function_family(f_expr, var)) {
        return SolveResult::success(std::move(*direct));
    }

    auto finite = solve_finite_vector_core(f_expr, var, context, opts);
    if (!finite) return SolveResult::failure(finite.error());
    if (finite.value().empty()) {
        return SolveResult::failure(
            CasErrc::Inconclusive,
            "expression is outside the supported finite exact solve domain",
            operation);
    }

    std::vector<FiniteSolution> solutions;
    solutions.reserve(finite.value().size());
    for (auto& value : finite.value()) {
        solutions.push_back(FiniteSolution{std::move(value), 1, {}});
    }
    return SolveResult::success(
        SolutionSet{FiniteSolutions{std::move(solutions)}});
}

lamina::SolveResult lamina::solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts) {
    ComputationContext context;
    return solve_equation(expr, var, context, opts);
}

static SolveVectorResult solve_finite_vector_core(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_dispatch_vector";
    if (!expr || !lamina::detail::node(expr)) {
        return SolveVectorResult::failure(CasErrc::InvalidArgument,
                                          "expression cannot be null",
                                          operation);
    }
    if (var.empty()) {
        return SolveVectorResult::failure(CasErrc::InvalidArgument,
                                          "solve variable cannot be empty",
                                          operation);
    }

    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return SolveVectorResult::failure(CasErrc::InvalidArgument,
                                          "expression simplification failed",
                                          operation);
    }

    std::shared_ptr<SymbolicExpr> f_expr = simplified;
    if (auto rel = std::dynamic_pointer_cast<const RelationalNode>(lamina::detail::node(simplified))) {
        if (rel->op() == RelationalNode::Op::EQ) {
            auto left = lamina::detail::make_expression_ptr(rel->left());
            auto right = lamina::detail::make_expression_ptr(rel->right());
            f_expr = SymbolicExpr::add(
                left, SymbolicExpr::multiply(right, SymbolicExpr::number(-1)));
        }
    }

    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(f_expr, var);
    if (!polynomial.is_zero() && polynomial.degree() >= 1) {
        if (polynomial.degree() <= 4) {
            auto results = solve_polynomial_closed(polynomial, var);
            results = solve_filter_rootof_results(std::move(results), opts);
            if (!results.empty()) return SolveVectorResult::success(std::move(results));
        }

        if (polynomial.degree() > 4) {
            auto results = solve_by_factoring(polynomial, var);
            results = solve_filter_rootof_results(std::move(results), opts);
            if (!results.empty()) return SolveVectorResult::success(std::move(results));

            if (opts.return_rootof) {
                return SolveVectorResult::success(make_rootof_solutions(polynomial, var));
            }
        }
    }

    auto transcendental_results = solve_transcendental(f_expr, var);
    if (!transcendental_results.empty()) {
        return SolveVectorResult::success(std::move(transcendental_results));
    }

    if (contains_transcendental_of_var(f_expr, var)) {
        auto substitution = detect_trans_substitutions(f_expr, var);
        if (substitution.mappings.empty() || !is_polynomial_after_substitution(substitution)) {
            if (opts.allow_numeric) {
                auto mixed = solve_mixed_transcendental_checked(
                    f_expr, var, opts, context);
                if (!mixed) {
                    return SolveVectorResult::failure(mixed.error());
                }
                if (mixed.value().completeness == Completeness::Complete) {
                    return SolveVectorResult::success(
                        std::move(mixed.value().value));
                }
                return SolveVectorResult::failure(
                    CasErrc::Inconclusive,
                    mixed.value().reason.empty()
                        ? "mixed transcendental search is incomplete"
                        : mixed.value().reason,
                    operation);
            }
            return SolveVectorResult::failure(
                CasErrc::Inconclusive,
                "mixed transcendental equation is outside the finite-vector support domain",
                operation);
        }
    }

    if (opts.allow_numeric) {
        auto numeric_roots = solve_numeric_checked(f_expr, var, context, opts);
        if (!numeric_roots) return SolveVectorResult::failure(numeric_roots.error());
        if (!numeric_roots.value().empty()) {
            std::vector<std::shared_ptr<SymbolicExpr>> results;
            results.reserve(numeric_roots.value().size());
            for (const auto& root : numeric_roots.value()) {
                results.push_back(SymbolicExpr::number(root.value));
            }
            return SolveVectorResult::success(std::move(results));
        }
    }

    return SolveVectorResult::success({});
}

lamina::FiniteSolveResult lamina::solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_finite_projection";
    auto solved = solve_equation(expr, var, context, opts);
    if (!solved) return SolveVectorResult::failure(solved.error());
    if (std::holds_alternative<EmptySolutions>(solved.value())) {
        return SolveVectorResult::success({});
    }
    auto* finite = std::get_if<FiniteSolutions>(&solved.value());
    if (!finite) {
        return SolveVectorResult::failure(
            CasErrc::Inconclusive,
            "solution set is not finitely enumerable", operation);
    }
    std::size_t value_count = 0;
    for (const auto& solution : finite->values) {
        value_count += solution.multiplicity;
    }
    std::vector<std::shared_ptr<SymbolicExpr>> values;
    values.reserve(value_count);
    for (const auto& solution : finite->values) {
        for (std::size_t copy = 0; copy < solution.multiplicity; ++copy) {
            values.push_back(solution.value);
        }
    }
    return SolveVectorResult::success(std::move(values));
}

lamina::FiniteSolveResult lamina::solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const SolveOptions& options) {
    ComputationContext context;
    return solve_finite_checked(
        expression, variable, context, options);
}



std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
SymbolicExpr::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& vars) {
    std::vector<SymbolicExpr> values;
    values.reserve(equations.size());
    for (const auto& equation : equations) {
        if (!equation) return {};
        if (auto relation =
                std::dynamic_pointer_cast<const RelationalNode>(
                    lamina::detail::node(equation))) {
            if (relation->op() != RelationalNode::Op::EQ) return {};
            auto normalized = SymbolicExpr::add(
                lamina::detail::make_expression_ptr(relation->left()),
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-1),
                    lamina::detail::make_expression_ptr(relation->right())));
            values.push_back(*normalized->simplify());
        } else {
            values.push_back(*equation);
        }
    }
    auto solved = lamina::Solver::solve_linear_system(values, vars);
    if (solved.empty()) return {};
    std::map<std::string, std::shared_ptr<SymbolicExpr>> projected;
    for (auto& [variable, value] : solved) {
        projected.emplace(
            variable, std::make_shared<SymbolicExpr>(std::move(value)));
    }
    return {std::move(projected)};
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters) {
    return ParametricSolver::solve_system(equations, unknowns, parameters);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda_name) {
    if (!mat) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!mat_node) return SymbolicExpr::number(0);
    size_t n = mat_node->rows();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    auto lambda = SymbolicExpr::variable(lambda_name);

    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
            auto val = lamina::detail::make_expression_ptr(mat_node->get(i,j));
            if (i == j) {
                data[i][j] = SymbolicExpr::add(val, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
            } else {
                data[i][j] = val;
            }
        }
    }

    auto poly_mat = SymbolicExpr::matrix(data);
    return lamina::matrix_determinant_checked(poly_mat).value();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) {
    auto cp = charpoly(mat, "lambda");
    auto solutions = lamina::solve_finite_checked(cp, "lambda").value();

    std::vector<std::shared_ptr<SymbolicExpr>> distinct_solutions;

    std::set<std::string> seen;
    for(auto& s : solutions) {
        auto str = s->to_string();
        if (seen.find(str) == seen.end()) {
            seen.insert(str);
            distinct_solutions.push_back(s);
        }
    }

    std::vector<std::shared_ptr<const SymbolicNode>> vec_nodes;
    for(auto& s : distinct_solutions) vec_nodes.push_back(lamina::detail::node(s));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_data;
    mat_data.push_back(distinct_solutions);
    return SymbolicExpr::matrix(mat_data);
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> SymbolicExpr::eigenvectors(const std::shared_ptr<SymbolicExpr>& mat) {
    auto evals_expr = eigenvalues(mat);

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> result;

    auto mat_node = evals_expr ? std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(evals_expr)) : nullptr;
    if (!mat_node) {
        return {};
    }

    size_t num_evals = mat_node->cols();

    auto A_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!A_node) return {};
    size_t n = A_node->rows();

    for(size_t i=0; i<num_evals; ++i) {
        auto lambda_node = mat_node->get(0, i);
        auto lambda = lamina::detail::make_expression_ptr(lambda_node);

        std::vector<std::shared_ptr<SymbolicExpr>> equations;
        std::vector<std::string> vars;
        for(size_t k=0; k<n; ++k) vars.push_back("v" + std::to_string(k));

        for(size_t i=0; i<n; ++i) {
            std::vector<std::shared_ptr<SymbolicExpr>> terms;
            for(size_t j=0; j<n; ++j) {
                auto a_ij = lamina::detail::make_expression_ptr(A_node->get(i,j));
                std::shared_ptr<SymbolicExpr> coeff;
                if (i == j) {
                    coeff = SymbolicExpr::add(a_ij, SymbolicExpr::multiply(lambda, SymbolicExpr::number(-1)));
                } else {
                    coeff = a_ij;
                }

                auto var = SymbolicExpr::variable(vars[j]);
                terms.push_back(SymbolicExpr::multiply(coeff, var));
            }
            auto row_eq = terms[0];
            for(size_t k=1; k<terms.size(); ++k) row_eq = SymbolicExpr::add(row_eq, terms[k]);
            equations.push_back(row_eq->simplify());
        }

        auto sols = solve_system(equations, vars);

        std::vector<std::shared_ptr<SymbolicExpr>> eigenvec;
        bool sols_usable = !sols.empty();
        if (sols_usable) {
            for(const auto& v : vars) {
                auto it = sols[0].find(v);
                if (it == sols[0].end()) { sols_usable = false; break; }
                eigenvec.push_back(it->second);
            }
        }

        bool is_non_zero = false;
        if (sols_usable) {
            for(auto& x : eigenvec) if(!x->is_zero()) is_non_zero = true;
        }

        if (sols_usable && is_non_zero) {
             std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data;
             for(auto& val : eigenvec) col_vec_data.push_back({val});
             result.push_back({lambda, {SymbolicExpr::matrix(col_vec_data)}});
        } else {
             if (n > 1) {
                  bool found_vec = false;
                  for(int free_var_idx = n-1; free_var_idx >= 0 && !found_vec; --free_var_idx) {
                       std::vector<std::shared_ptr<SymbolicExpr>> sub_eqs;
                       std::vector<std::string> sub_vars;

                       for(size_t k=0; k<n; ++k) {
                           auto eq_sub = equations[k]->substitute(vars[free_var_idx], SymbolicExpr::number(1))->simplify();
                           if (!eq_sub->is_zero()) {
                               sub_eqs.push_back(eq_sub);
                           }
                       }

                       for(size_t k=0; k<n; ++k) {
                           if (k != (size_t)free_var_idx) sub_vars.push_back(vars[k]);
                       }

                       auto sub_sols = solve_system(sub_eqs, sub_vars);
                       if (!sub_sols.empty()) {
                            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> col_vec_data(n);
                            bool ok = true;
                            for(size_t k=0; k<n; ++k) {
                                if (k == (size_t)free_var_idx) { col_vec_data[k] = {SymbolicExpr::number(1)}; continue; }
                                auto it = sub_sols[0].find(vars[k]);
                                if (it == sub_sols[0].end()) { ok = false; break; }
                                col_vec_data[k] = {it->second};
                            }
                            if (ok) {
                                auto vec_expr = SymbolicExpr::matrix(col_vec_data);
                                result.push_back({lambda, {vec_expr}});
                                found_vec = true;
                            }
                       }
                  }
             }
        }
    }

    return result;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat) return mat;
    auto m_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!m_node) return mat;
    size_t r = m_node->rows();
    size_t c = m_node->cols();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> new_data(c, std::vector<std::shared_ptr<SymbolicExpr>>(r));

    for(size_t i=0; i<r; ++i) {
        for(size_t j=0; j<c; ++j) {
            new_data[j][i] = lamina::detail::make_expression_ptr(m_node->get(i,j));
        }
    }

    return SymbolicExpr::matrix(new_data);
}
