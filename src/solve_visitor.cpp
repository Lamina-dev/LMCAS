#include "../include/symbolic.hpp"
#include "../include/parametric_solver.hpp"
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

std::shared_ptr<SymbolicExpr> get_coeff(const Polynomial<SymbolicPolyCoeff>& p, int deg) {
    if (deg < 0 || deg > p.degree()) return SymbolicExpr::number(0);
    return p.coeffs[deg].val;
}

static bool solve_is_minus_one(const std::shared_ptr<const NumberNode>& number) {
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(-1);
    }
    return std::get<lmmc_real_t>(number->value()) == -1.0;
}

static bool solve_contains_rootof_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::RootOf) return true;
        for (const auto& arg : function->arguments()) {
            if (solve_contains_rootof_node(arg)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (solve_contains_rootof_node(operand)) return true;
        }
        return false;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (solve_contains_rootof_node(operand)) return true;
        }
        return false;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return solve_contains_rootof_node(power->base()) ||
               solve_contains_rootof_node(power->exponent());
    }
    if (auto relation = std::dynamic_pointer_cast<const RelationalNode>(node)) {
        return solve_contains_rootof_node(relation->left()) ||
               solve_contains_rootof_node(relation->right());
    }
    return false;
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

static std::shared_ptr<const SymbolicNode> solve_simplify_with_nonzero_pivots(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return node;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            operands.push_back(solve_simplify_with_nonzero_pivots(operand));
        }
        return lamina::detail::make_node<AddNode>(std::move(operands));
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return lamina::detail::make_node<PowerNode>(
            solve_simplify_with_nonzero_pivots(power->base()),
            solve_simplify_with_nonzero_pivots(power->exponent()));
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto simplified = solve_simplify_with_nonzero_pivots(operand);
            if (auto nested = std::dynamic_pointer_cast<const MultiplyNode>(simplified)) {
                operands.insert(operands.end(), nested->operands().begin(), nested->operands().end());
            } else {
                operands.push_back(simplified);
            }
        }

        std::vector<bool> consumed(operands.size(), false);
        for (size_t i = 0; i < operands.size(); ++i) {
            if (consumed[i]) continue;
            auto reciprocal = std::dynamic_pointer_cast<const PowerNode>(operands[i]);
            auto exponent = reciprocal ? std::dynamic_pointer_cast<const NumberNode>(reciprocal->exponent()) : nullptr;
            if (!reciprocal || !solve_is_minus_one(exponent)) continue;

            for (size_t j = 0; j < operands.size(); ++j) {
                if (i == j || consumed[j]) continue;
                if (reciprocal->base()->compare(*operands[j]) == 0) {
                    consumed[i] = true;
                    consumed[j] = true;
                    break;
                }
            }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> remaining;
        for (size_t i = 0; i < operands.size(); ++i) {
            if (!consumed[i]) remaining.push_back(operands[i]);
        }
        if (remaining.empty()) return lamina::detail::make_node<NumberNode>(BigInt(1));
        if (remaining.size() == 1) return remaining.front();
        return lamina::detail::make_node<MultiplyNode>(std::move(remaining));
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& arg : function->arguments()) {
            args.push_back(solve_simplify_with_nonzero_pivots(arg));
        }
        return lamina::detail::make_node<FunctionNode>(function->type(), std::move(args));
    }

    return node;
}

static std::shared_ptr<SymbolicExpr> solve_simplify_solution_value(
    const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return expr;
    auto node = solve_simplify_with_nonzero_pivots(lamina::detail::node(expr));
    return lamina::detail::make_expression_ptr(node)->simplify();
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

lamina::SolveResult lamina::solve_dispatch_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_dispatch";
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
                "solve_dispatch_checked accepts equations, not inequalities",
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
            return SolveResult::success(SolutionSet::universal());
        }
        if (exact.degree() == 0) {
            return SolveResult::success(SolutionSet::empty());
        }

        std::vector<FiniteSolution> solutions;
        solutions.reserve(static_cast<std::size_t>(exact.degree()));
        for (int index = 0; index < exact.degree(); ++index) {
            auto root_step = context.consume_steps(1, operation);
            if (!root_step) return SolveResult::failure(root_step.error());
            solutions.push_back(FiniteSolution{
                SymbolicExpr::root_of(f_expr, var, index), 1, {}});
        }
        return SolveResult::success(SolutionSet::finite(std::move(solutions)));
    }

    if (!opts.allow_numeric) {
        return SolveResult::success(SolutionSet::inconclusive(
            "expression is outside the exact rational-polynomial support domain"));
    }

    auto numeric = solve_numeric_checked(f_expr, var, context, opts);
    if (!numeric) return SolveResult::failure(numeric.error());
    if (numeric.value().empty()) {
        return SolveResult::success(SolutionSet::inconclusive(
            "numeric iteration produced no verified candidate"));
    }

    std::vector<FiniteSolution> solutions;
    solutions.reserve(numeric.value().size());
    for (const auto& root : numeric.value()) {
        solutions.push_back(FiniteSolution{
            SymbolicExpr::number(root.value), 1, {}});
    }
    return SolveResult::success(SolutionSet::finite(std::move(solutions)));
}

lamina::SolveResult lamina::solve_dispatch_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts) {
    ComputationContext context;
    return solve_dispatch_checked(expr, var, context, opts);
}

lamina::SolveVectorResult lamina::solve_dispatch_vector_checked(
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
            auto L = lamina::detail::make_expression_ptr(rel->left());
            auto R = lamina::detail::make_expression_ptr(rel->right());
            f_expr = SymbolicExpr::add(L, SymbolicExpr::multiply(R, SymbolicExpr::number(-1)));
        }
    }

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(f_expr, var);
    if (!poly.is_zero() && poly.degree() >= 1) {

        if (poly.degree() <= 4) {
            auto results = solve_polynomial_closed(poly, var);
            results = solve_filter_rootof_results(std::move(results), opts);
            if (!results.empty()) return SolveVectorResult::success(std::move(results));
        }

        if (poly.degree() > 4) {
            auto results = solve_by_factoring(poly, var);
            results = solve_filter_rootof_results(std::move(results), opts);
            if (!results.empty()) return SolveVectorResult::success(std::move(results));

            if (opts.return_rootof) {
                return SolveVectorResult::success(make_rootof_solutions(poly, var));
            }
        }
    }

    auto trans_results = solve_transcendental(f_expr, var);
    if (!trans_results.empty()) return SolveVectorResult::success(std::move(trans_results));

    /// 混合超越方程路径：含超越函数且换元无法化为多项式时，委托给混合求解器
    if (contains_transcendental_of_var(f_expr, var)) {
        auto sub_result = detect_trans_substitutions(f_expr, var);
        if (sub_result.mappings.empty() || !is_polynomial_after_substitution(sub_result)) {
            if (opts.allow_numeric) {
                auto mixed_results = solve_mixed_transcendental(f_expr, var, opts);
                if (!mixed_results.empty()) {
                    return SolveVectorResult::success(std::move(mixed_results));
                }
                return SolveVectorResult::failure(
                    CasErrc::Inconclusive,
                    "mixed transcendental numeric search produced no verified vector result",
                    operation);
            }
            return SolveVectorResult::failure(
                CasErrc::Inconclusive,
                "mixed transcendental equation is outside the checked vector support domain",
                operation);
        }
    }

    if (opts.allow_numeric) {
        auto numeric_roots = solve_numeric_checked(f_expr, var, context, opts);
        if (!numeric_roots) return SolveVectorResult::failure(numeric_roots.error());
        if (!numeric_roots.value().empty()) {

            std::vector<std::shared_ptr<SymbolicExpr>> results;
            results.reserve(numeric_roots.value().size());
            for (const auto& nr : numeric_roots.value()) {
                results.push_back(SymbolicExpr::number(nr.value));
            }
            return SolveVectorResult::success(std::move(results));
        }
    }

    return SolveVectorResult::success({});
}

lamina::SolveVectorResult lamina::solve_dispatch_vector_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts) {
    ComputationContext context;
    return solve_dispatch_vector_checked(expr, var, context, opts);
}

std::vector<std::shared_ptr<SymbolicExpr>> lamina::solve_dispatch(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts) {

    ComputationContext context;
    auto result = solve_dispatch_vector_checked(expr, var, context, opts);
    return result ? result.value() : std::vector<std::shared_ptr<SymbolicExpr>>{};
}

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
    if (!eq) return {};
    auto simplified_eq = eq->simplify();

    if (auto rel = std::dynamic_pointer_cast<const RelationalNode>(lamina::detail::node(simplified_eq))) {
        if (rel->op() != RelationalNode::Op::EQ) {
            auto L = lamina::detail::make_expression_ptr(rel->left());
            auto R = lamina::detail::make_expression_ptr(rel->right());
            auto diff = SymbolicExpr::add(L, SymbolicExpr::multiply(R, SymbolicExpr::number(-1)));

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(diff, var_name);

            if (!poly.is_zero()) {
                int max_deg = poly.degree();
                if (max_deg == 1) {
                    auto a = get_coeff(poly, 1);
                    auto b = get_coeff(poly, 0);

                    auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
                    auto a_inv = SymbolicExpr::power(a, SymbolicExpr::number(-1));

                    bool flip = false;

                    try {
                        auto a_simp = a->simplify();
                        if (auto num_a = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(a_simp))) {
                            if (std::holds_alternative<double>(num_a->value())) {
                                 if (std::get<double>(num_a->value()) < 0) flip = true;
                            } else if (std::holds_alternative<BigInt>(num_a->value())) {
                                 if (std::get<BigInt>(num_a->value()).IsNegative()) flip = true;
                            } else if (std::holds_alternative<Rational>(num_a->value())) {
                                 if (std::get<Rational>(num_a->value()).get_numerator().IsNegative()) flip = true;
                            }
                        }
                    } catch(...) {}

                    RelationalNode::Op new_op = rel->op();
                    if (flip) {
                        switch(rel->op()) {
                            case RelationalNode::Op::LT: new_op = RelationalNode::Op::GT; break;
                            case RelationalNode::Op::GT: new_op = RelationalNode::Op::LT; break;
                            case RelationalNode::Op::LEQ: new_op = RelationalNode::Op::GEQ; break;
                            case RelationalNode::Op::GEQ: new_op = RelationalNode::Op::LEQ; break;
                            default: break;
                        }
                    }

                    auto rhs = SymbolicExpr::multiply(neg_b, a_inv)->simplify();
                    auto var_node = lamina::detail::make_node<VariableNode>(var_name);
                    auto res_node = lamina::detail::make_node<RelationalNode>(var_node, lamina::detail::node(rhs), new_op);

                    return { lamina::detail::make_expression_ptr(res_node) };
                }
            }
            return {};
        }
    }

    return lamina::solve_dispatch(simplified_eq, var_name, lamina::SolveOptions{});
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(const std::vector<std::shared_ptr<SymbolicExpr>>& equations, const std::vector<std::string>& vars) {
    size_t n = vars.size();
    size_t m = equations.size();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(m, std::vector<std::shared_ptr<SymbolicExpr>>(n + 1));

    for(size_t i=0; i<m; ++i) {
        auto eq = equations[i]->expand();

        auto C = eq;
        for(const auto& v : vars) C = C->substitute(v, SymbolicExpr::number(0));
        C = C->simplify();
        A[i][n] = SymbolicExpr::multiply(C, SymbolicExpr::number(-1));

        for(size_t j=0; j<n; ++j) {
            auto p = symbolic_to_poly<SymbolicPolyCoeff>(eq, vars[j]);
            if (p.degree() >= 1) {
                 A[i][j] = get_coeff(p, 1);
            } else {
                 A[i][j] = SymbolicExpr::number(0);
            }
        }
    }

    std::vector<size_t> pivot_col_for_row;
    int sign;
    gaussian_eliminate(A, m, n, pivot_col_for_row, sign);

    std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
    for(const auto& v : vars) solution[v] = SymbolicExpr::number(0);

    for(size_t r=0; r<m; ++r) {
        size_t pivot_col = pivot_col_for_row[r];
        if (pivot_col != (size_t)-1) {
            solution[vars[pivot_col]] = solve_simplify_solution_value(A[r][n]);
        }
    }
    return { solution };
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters) {
    return ParametricSolver::solve_system(equations, unknowns, parameters);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::determinant(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat) return SymbolicExpr::number(0);
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!mat_node || mat_node->rows() != mat_node->cols()) return SymbolicExpr::number(0);

    size_t n = mat_node->rows();
    if (n == 1) return lamina::detail::make_expression_ptr(mat_node->get(0,0));
    if (n == 2) {
        auto a = lamina::detail::make_expression_ptr(mat_node->get(0,0));
        auto b = lamina::detail::make_expression_ptr(mat_node->get(0,1));
        auto c = lamina::detail::make_expression_ptr(mat_node->get(1,0));
        auto d = lamina::detail::make_expression_ptr(mat_node->get(1,1));
        return SymbolicExpr::add(SymbolicExpr::multiply(a,d), SymbolicExpr::multiply(SymbolicExpr::multiply(b,c), SymbolicExpr::number(-1)))->simplify();
    }

    if (n > 3) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
        for(size_t i=0; i<n; ++i) {
            for(size_t j=0; j<n; ++j) {
                A[i][j] = lamina::detail::make_expression_ptr(mat_node->get(i,j));
            }
        }

        int sign = 1;
        auto det = SymbolicExpr::number(1);

        for (size_t col = 0; col < n; ++col) {
            size_t pivot_row = col;
            while (pivot_row < n && A[pivot_row][col]->is_zero()) {
                pivot_row++;
            }
            if (pivot_row == n) return SymbolicExpr::number(0);

            if (pivot_row != col) {
                std::swap(A[col], A[pivot_row]);
                sign = -sign;
            }

            auto pivot = A[col][col];
            det = SymbolicExpr::multiply(det, pivot)->simplify();
            auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));

            for (size_t r = col + 1; r < n; ++r) {
                auto factor = A[r][col];
                if (!factor->is_zero()) {
                    auto mult = SymbolicExpr::multiply(factor, pivot_inv)->simplify();
                    auto neg_mult = SymbolicExpr::multiply(mult, SymbolicExpr::number(-1));
                    for (size_t c = col + 1; c < n; ++c) {
                        auto term = SymbolicExpr::multiply(neg_mult, A[col][c]);
                        A[r][c] = SymbolicExpr::add(A[r][c], term)->simplify();
                    }
                    A[r][col] = SymbolicExpr::number(0);
                }
            }
        }

        if (sign == -1) {
            det = SymbolicExpr::multiply(det, SymbolicExpr::number(-1))->simplify();
        }
        return det;
    }

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for(size_t c=0; c<n; ++c) {
        auto elem = lamina::detail::make_expression_ptr(mat_node->get(0,c));
        if (elem->is_zero()) continue;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> minor_data;
        for(size_t r=1; r<n; ++r) {
            std::vector<std::shared_ptr<SymbolicExpr>> row;
            for(size_t k=0; k<n; ++k) {
                if (k == c) continue;
                row.push_back(lamina::detail::make_expression_ptr(mat_node->get(r,k)));
            }
            minor_data.push_back(row);
        }
        auto minor_mat = SymbolicExpr::matrix(minor_data);
        auto minor_det = SymbolicExpr::determinant(minor_mat);

        auto term = SymbolicExpr::multiply(elem, minor_det);
        if (c % 2 == 1) term = SymbolicExpr::multiply(term, SymbolicExpr::number(-1));
        terms.push_back(term);
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    auto result = terms[0];
    for(size_t k=1; k<terms.size(); ++k) result = SymbolicExpr::add(result, terms[k]);
    return result->simplify();
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
    return SymbolicExpr::determinant(poly_mat);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) {
    auto cp = charpoly(mat, "lambda");
    auto solutions = solve(cp, "lambda");

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

std::shared_ptr<SymbolicExpr> SymbolicExpr::rref(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat) return mat;
    auto m_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!m_node) return mat;
    size_t rows = m_node->rows();
    size_t cols = m_node->cols();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> data(rows, std::vector<std::shared_ptr<SymbolicExpr>>(cols));
    for(size_t i=0; i<rows; ++i) {
        for(size_t j=0; j<cols; ++j) {
            data[i][j] = lamina::detail::make_expression_ptr(m_node->get(i,j));
        }
    }

    size_t lead = 0;
    for (size_t r = 0; r < rows; ++r) {
        if (cols <= lead) break;

        size_t i = r;
        while (true) {
             if (i >= rows) {
                 i = r;
                 lead++;
                 if (cols == lead) goto end_loops;
                 continue;
             }
             auto val = data[i][lead]->simplify();
             if (!val->get_number_value_is_zero()) {
                 break;
             }
             i++;
        }

        std::swap(data[i], data[r]);

        auto pivot = data[r][lead];
        auto pivot_inv = SymbolicExpr::power(pivot, SymbolicExpr::number(-1));

        for (size_t j = 0; j < cols; ++j) {
            data[r][j] = SymbolicExpr::multiply(data[r][j], pivot_inv)->simplify();
        }

        for (size_t k = 0; k < rows; ++k) {
            if (k != r) {
                auto factor = data[k][lead];
                if (factor->simplify()->get_number_value_is_zero()) continue;

                auto neg_factor = SymbolicExpr::multiply(factor, SymbolicExpr::number(-1));
                for (size_t j = 0; j < cols; ++j) {
                    auto term = SymbolicExpr::multiply(neg_factor, data[r][j]);
                    data[k][j] = SymbolicExpr::add(data[k][j], term)->simplify();
                }
            }
        }
        lead++;
    }
    end_loops:;

    return SymbolicExpr::matrix(data);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::inverse(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat) return nullptr;
    auto m_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(mat));
    if (!m_node) return nullptr;
    size_t n = m_node->rows();
    if (n != m_node->cols()) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> aug_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(2*n));

    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
            aug_data[i][j] = lamina::detail::make_expression_ptr(m_node->get(i,j));
        }
        for(size_t j=0; j<n; ++j) {
            aug_data[i][n+j] = (i==j ? SymbolicExpr::number(1) : SymbolicExpr::number(0));
        }
    }

    auto aug_mat = SymbolicExpr::matrix(aug_data);
    auto rref_mat = rref(aug_mat);

    if (!rref_mat) return nullptr;
    auto rref_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(rref_mat));

    for(size_t i=0; i<n; ++i) {
        auto diag = lamina::detail::make_expression_ptr(rref_node->get(i,i))->simplify();
        if (!lamina::detail::node(diag)->is_one()) return nullptr;
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> inv_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for(size_t i=0; i<n; ++i) {
        for(size_t j=0; j<n; ++j) {
             inv_data[i][j] = lamina::detail::make_expression_ptr(rref_node->get(i, n+j));
        }
    }

    return SymbolicExpr::matrix(inv_data);
}
