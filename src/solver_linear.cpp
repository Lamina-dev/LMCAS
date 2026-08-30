#include "solver.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/exact_matrix.hpp"
#include "assumption_context.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <queue>
#include <unordered_map>
#include <optional>
#include <limits>
#include "internal/solver_support.hpp"

namespace lamina {
using namespace solver_detail;


std::shared_ptr<SymbolicExpr> solver_detail::to_ptr(const SymbolicExpr& expr) {
    return lamina::detail::make_expression_ptr(expr);
}

static bool get_integer_value(const std::shared_ptr<const SymbolicNode>& node, long long& value) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        value = std::get<BigInt>(num->value()).to_int();
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        const auto& r = std::get<Rational>(num->value());
        if (!r.is_integer()) return false;
        value = r.to_BigInt().to_int();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        lmmc_real_t d = std::get<lmmc_real_t>(num->value());
        int eq;
        lmmc_double_nearly_equal_tol(d, std::round(d), 1e-12, 1e-12, &eq);
        if (!eq) return false;
        value = static_cast<long long>(std::llround(d));
        return true;
    }
    return false;
}

static std::shared_ptr<const NumberNode> add_number_nodes(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
    if (std::holds_alternative<lmmc_real_t>(a->value()) || std::holds_alternative<lmmc_real_t>(b->value())) {
        auto to_real = [](const auto& v) {
            if (std::holds_alternative<lmmc_real_t>(v)) return std::get<lmmc_real_t>(v);
            if (std::holds_alternative<Rational>(v)) return (lmmc_real_t)std::get<Rational>(v).to_double();
            return (lmmc_real_t)std::get<BigInt>(v).to_double();
        };
        lmmc_real_t r1 = to_real(a->value());
        lmmc_real_t r2 = to_real(b->value());
        return lamina::detail::make_node<NumberNode>(r1 + r2);
    }

    if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
        Rational r1 = std::holds_alternative<Rational>(a->value()) ? std::get<Rational>(a->value()) : Rational(std::get<BigInt>(a->value()));
        Rational r2 = std::holds_alternative<Rational>(b->value()) ? std::get<Rational>(b->value()) : Rational(std::get<BigInt>(b->value()));
        return lamina::detail::make_node<NumberNode>(r1 + r2);
    }

    return lamina::detail::make_node<NumberNode>(std::get<BigInt>(a->value()) + std::get<BigInt>(b->value()));
}

static std::shared_ptr<const NumberNode> multiply_number_nodes(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
    if (std::holds_alternative<lmmc_real_t>(a->value()) || std::holds_alternative<lmmc_real_t>(b->value())) {
        auto to_real = [](const auto& v) {
            if (std::holds_alternative<lmmc_real_t>(v)) return std::get<lmmc_real_t>(v);
            if (std::holds_alternative<Rational>(v)) return (lmmc_real_t)std::get<Rational>(v).to_double();
            return (lmmc_real_t)std::get<BigInt>(v).to_double();
        };
        lmmc_real_t r1 = to_real(a->value());
        lmmc_real_t r2 = to_real(b->value());
        return lamina::detail::make_node<NumberNode>(r1 * r2);
    }

    if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
        Rational r1 = std::holds_alternative<Rational>(a->value()) ? std::get<Rational>(a->value()) : Rational(std::get<BigInt>(a->value()));
        Rational r2 = std::holds_alternative<Rational>(b->value()) ? std::get<Rational>(b->value()) : Rational(std::get<BigInt>(b->value()));
        return lamina::detail::make_node<NumberNode>(r1 * r2);
    }

    return lamina::detail::make_node<NumberNode>(std::get<BigInt>(a->value()) * std::get<BigInt>(b->value()));
}

struct NodeLess {
    bool operator()(const std::shared_ptr<const SymbolicNode>& a, const std::shared_ptr<const SymbolicNode>& b) const {
        if (!a || !b) return a < b;
        return a->compare(*b) < 0;
    }
};

std::shared_ptr<SymbolicExpr> solver_detail::multiply_no_expand(
    const std::shared_ptr<const SymbolicNode>& term,
    const std::vector<std::shared_ptr<const SymbolicNode>>& den_factors
) {
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(term)) {
        factors.insert(factors.end(), mul->operands().begin(), mul->operands().end());
    } else if (term) {
        factors.push_back(term);
    }
    factors.insert(factors.end(), den_factors.begin(), den_factors.end());

    std::shared_ptr<const NumberNode> const_acc =
        lamina::detail::make_node<NumberNode>(BigInt(1));
    std::map<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const NumberNode>, NodeLess> bases;

    for (const auto& op : factors) {
        if (!op) continue;
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
            const_acc = multiply_number_nodes(const_acc, num);
            continue;
        }

        std::shared_ptr<const SymbolicNode> base = op;
        std::shared_ptr<const NumberNode> exp = lamina::detail::make_node<NumberNode>(BigInt(1));
        if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
            base = pow->base();
            if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                exp = e_num;
            }
        }

        auto it = bases.find(base);
        if (it == bases.end()) {
            bases[base] = exp;
        } else {
            bases[base] = add_number_nodes(it->second, exp);
        }
    }

    std::vector<std::shared_ptr<const SymbolicNode>> final_ops;
    if (!const_acc->is_one()) final_ops.push_back(const_acc);

    for (const auto& [base, exp] : bases) {
        if (exp->is_zero()) continue;
        if (exp->is_one()) final_ops.push_back(base);
        else final_ops.push_back(lamina::detail::make_node<PowerNode>(base, exp));
    }

    if (final_ops.empty()) return SymbolicExpr::number(1);
    if (final_ops.size() == 1) return lamina::detail::make_expression_ptr(final_ops[0]);
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(final_ops));
}

bool solver_detail::is_polynomial_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (std::dynamic_pointer_cast<const NumberNode>(node)) return true;
    if (std::dynamic_pointer_cast<const VariableNode>(node)) return true;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (!is_polynomial_node(op)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (!is_polynomial_node(op)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        long long exp = 0;
        if (!get_integer_value(pow->exponent(), exp)) return false;
        if (exp < 0) return false;
        return is_polynomial_node(pow->base());
    }

    return false;
}

std::shared_ptr<SymbolicExpr> solver_detail::multiply_factors(const std::vector<std::shared_ptr<const SymbolicNode>>& factors) {
    if (factors.empty()) return SymbolicExpr::number(1);
    auto res = lamina::detail::make_expression_ptr(factors[0]);
    for (size_t i = 1; i < factors.size(); ++i) {
        res = SymbolicExpr::multiply(res, lamina::detail::make_expression_ptr(factors[i]));
    }
    return res->simplify();
}

bool solver_detail::collect_denominator_factors(
    const std::shared_ptr<const SymbolicNode>& node,
    std::vector<std::shared_ptr<const SymbolicNode>>& den_factors,
    std::vector<std::shared_ptr<SymbolicExpr>>& den_constraints
) {
    if (!node) return false;

    if (std::dynamic_pointer_cast<const NumberNode>(node) || std::dynamic_pointer_cast<const VariableNode>(node)) {
        return true;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (!collect_denominator_factors(op, den_factors, den_constraints)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (!collect_denominator_factors(op, den_factors, den_constraints)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        long long exp = 0;
        if (!get_integer_value(pow->exponent(), exp)) return false;

        if (exp < 0) {
            if (!is_polynomial_node(pow->base())) return false;
            long long k = -exp;
            if (k == 1) {
                den_factors.push_back(pow->base());
            } else {
                den_factors.push_back(SymbolicFactory::create_power(pow->base(), SymbolicFactory::create_number(BigInt(k))));
            }
            den_constraints.push_back(lamina::detail::make_expression_ptr(pow->base()));
            return true;
        }

        if (!collect_denominator_factors(pow->base(), den_factors, den_constraints)) return false;
        return true;
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (const auto& arg : func->arguments()) {
            if (!collect_denominator_factors(arg, den_factors, den_constraints)) return false;
        }
        return true;
    }

    return false;
}

static std::pair<SymbolicExpr, SymbolicExpr> isolate_linear_coeff(const SymbolicExpr& expr, const std::string& var) {

    auto expr_ptr = to_ptr(expr);
    auto A_ptr = expr_ptr->differentiate(var);

    std::vector<std::shared_ptr<const SymbolicNode>> mops;
    mops.push_back(lamina::detail::node(A_ptr));
    mops.push_back(SymbolicFactory::create_variable(var));
    auto term_Ax = SymbolicFactory::create_multiply(mops);

    std::vector<std::shared_ptr<const SymbolicNode>> nops;
    nops.push_back(SymbolicFactory::create_number(BigInt(-1)));
    nops.push_back(term_Ax);
    auto neg_term = SymbolicFactory::create_multiply(nops);

    std::vector<std::shared_ptr<const SymbolicNode>> aops;
    aops.push_back(lamina::detail::node(expr));
    aops.push_back(neg_term);
    auto B_node = SymbolicFactory::create_add(aops);

    auto B_expr = lamina::detail::expression_from_node(B_node);
    auto B_simp = to_ptr(B_expr)->simplify();

    auto A_expr = lamina::detail::expression_from_node(lamina::detail::node(A_ptr));
    auto A_simp = to_ptr(A_expr)->simplify();

    return {lamina::detail::expression_from_node(lamina::detail::node(A_simp)),
            lamina::detail::expression_from_node(lamina::detail::node(B_simp))};
}

std::map<std::string, SymbolicExpr> Solver::solve_linear_system(
    const std::vector<SymbolicExpr>& equations_in,
    const std::vector<std::string>& variables)
{

    size_t num_vars = variables.size();
    size_t num_eqs = equations_in.size();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix(num_eqs, std::vector<std::shared_ptr<SymbolicExpr>>(num_vars + 1));

    for (size_t i = 0; i < num_eqs; ++i) {
        SymbolicExpr current_eq = equations_in[i];

        SymbolicExpr constant_part = current_eq;

        for (size_t j = 0; j < num_vars; ++j) {
            auto [coeff, remainder] = isolate_linear_coeff(constant_part, variables[j]);
            matrix[i][j] = lamina::detail::make_expression_ptr(coeff);

            constant_part = remainder;
        }

        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.push_back(SymbolicFactory::create_number(BigInt(-1)));
        ops.push_back(lamina::detail::node(constant_part));
        auto neg_const = SymbolicFactory::create_multiply(ops);
        matrix[i][num_vars] = lamina::detail::make_expression_ptr(neg_const);
    }

    detail::ExactMatrixData augmented{
        num_eqs, num_vars + 1, {}};
    augmented.entries.reserve(num_eqs * (num_vars + 1));
    for (const auto& row : matrix) {
        augmented.entries.insert(
            augmented.entries.end(), row.begin(), row.end());
    }
    ComputationContext context;
    auto solved = detail::solve_linear_exact(
        std::move(augmented), num_vars, context,
        "solve_linear_system");
    if (!solved ||
        std::holds_alternative<detail::InconsistentLinearSolution>(
            solved.value())) {
        return {};
    }

    std::vector<std::shared_ptr<SymbolicExpr>> values;
    if (auto* unique =
            std::get_if<detail::UniqueLinearSolution>(&solved.value())) {
        values = std::move(unique->values);
    } else {
        auto& parametric =
            std::get<detail::ParametricLinearSolution>(solved.value());
        values = parametric.particular;
        for (std::size_t parameter = 0;
             parameter < parametric.free_columns.size(); ++parameter) {
            auto free_variable = SymbolicExpr::variable(
                variables[parametric.free_columns[parameter]]);
            for (std::size_t column = 0; column < num_vars; ++column) {
                auto term = SymbolicExpr::multiply(
                    parametric.nullspace_basis[parameter][column],
                    free_variable);
                values[column] = SymbolicExpr::add(
                    values[column], term)->simplify();
            }
        }
    }

    std::map<std::string, SymbolicExpr> solution;
    for (std::size_t column = 0; column < num_vars; ++column) {
        solution.insert_or_assign(
            variables[column], *values[column]);
    }
    return solution;
}


} // namespace lamina
