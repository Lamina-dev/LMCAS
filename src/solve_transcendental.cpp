#include "solve_transcendental.hpp"
#include "symbolic_ast.hpp"
#include "solve_strategies.hpp"
#include "solve_polynomial.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace lamina {

static std::shared_ptr<SymbolicExpr> make_arcsin(const std::shared_ptr<SymbolicExpr>& c) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcSin,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(c)}));
}

static std::shared_ptr<SymbolicExpr> make_arccos(const std::shared_ptr<SymbolicExpr>& c) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(c)}));
}

static std::shared_ptr<SymbolicExpr> make_arctan(const std::shared_ptr<SymbolicExpr>& c) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcTan,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(c)}));
}

static std::shared_ptr<SymbolicExpr> make_pi() {
    return SymbolicExpr::number(LMMC_CONST_PI);
}

static bool try_evaluate_numeric(const std::shared_ptr<SymbolicExpr>& expr, lmmc_real_t& out) {
    if (!expr || !lamina::detail::node(expr)) return false;
    ComputationContext context;
    auto evaluated = evaluate_numeric(*expr, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return false;
    }
    out = static_cast<lmmc_real_t>(evaluated.value().value);
    return true;
}

struct InversePattern {
    FunctionNode::FuncType func_type;
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
};

static std::optional<InversePattern> decompose_trig_exp_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return std::nullopt;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
        if (func->arguments().size() == 1 && expression_depends_on_variable(func->arguments()[0], var)) {
            auto ft = func->type();
            if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                ft == FunctionNode::FuncType::Ln) {
                return InversePattern{
                    ft,
                    lamina::detail::make_expression_ptr(func->arguments()[0]),
                    SymbolicExpr::number(0)
                };
            }
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {

        std::shared_ptr<const FunctionNode> func_term = nullptr;
        std::shared_ptr<const SymbolicNode> func_coeff = nullptr;
        std::vector<std::shared_ptr<const SymbolicNode>> const_terms;

        for (auto& op : add->operands()) {

            if (auto f = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                if (f->arguments().size() == 1 && expression_depends_on_variable(f->arguments()[0], var)) {
                    auto ft = f->type();
                    if ((ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                         ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                         ft == FunctionNode::FuncType::Ln) && !func_term) {
                        func_term = f;
                        func_coeff = lamina::detail::make_node<NumberNode>(BigInt(1));
                        continue;
                    }
                }
            }

            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
                std::shared_ptr<const FunctionNode> mf = nullptr;
                std::vector<std::shared_ptr<const SymbolicNode>> coeff_parts;
                for (auto& mop : mul->operands()) {
                    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(mop)) {
                        if (f->arguments().size() == 1 && expression_depends_on_variable(f->arguments()[0], var)) {
                            auto ft = f->type();
                            if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                                ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                                ft == FunctionNode::FuncType::Ln) {
                                if (!mf) { mf = f; continue; }
                            }
                        }
                    }
                    coeff_parts.push_back(mop);
                }
                if (mf && !func_term && !expression_depends_on_variable(SymbolicFactory::create_multiply(coeff_parts), var)) {
                    func_term = mf;
                    func_coeff = SymbolicFactory::create_multiply(coeff_parts);
                    continue;
                }
            }

            if (!expression_depends_on_variable(op, var)) {
                const_terms.push_back(op);
            } else {

                return std::nullopt;
            }
        }

        if (func_term) {

            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = lamina::detail::make_expression_ptr(const_terms[0]);
            } else {
                const_sum = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(const_terms));
            }

            auto neg_const = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1));
            auto coeff_expr = lamina::detail::make_expression_ptr(func_coeff);
            auto rhs = SymbolicExpr::divide(neg_const, coeff_expr)->simplify();

            return InversePattern{
                func_term->type(),
                lamina::detail::make_expression_ptr(func_term->arguments()[0]),
                rhs
            };
        }
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        for (auto& op : mul->operands()) {
            if (auto f = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                if (f->arguments().size() == 1 && expression_depends_on_variable(f->arguments()[0], var)) {
                    auto ft = f->type();
                    if (ft == FunctionNode::FuncType::Sin || ft == FunctionNode::FuncType::Cos ||
                        ft == FunctionNode::FuncType::Tan || ft == FunctionNode::FuncType::Exp ||
                        ft == FunctionNode::FuncType::Ln) {
                        return InversePattern{
                            ft,
                            lamina::detail::make_expression_ptr(f->arguments()[0]),
                            SymbolicExpr::number(0)
                        };
                    }
                }
            }
        }
    }

    return std::nullopt;
}

struct LambertWPattern {
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
};

static bool nodes_equal(const std::shared_ptr<const SymbolicNode>& a, const std::shared_ptr<const SymbolicNode>& b) {
    if (!a || !b) return a == b;
    return a->equals(*b);
}

static std::optional<LambertWPattern> decompose_lambert_w_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return std::nullopt;

    auto check_mul_is_g_exp_g = [&](const std::shared_ptr<const MultiplyNode>& mul)
        -> std::optional<std::shared_ptr<SymbolicExpr>> {

        for (size_t i = 0; i < mul->operands().size(); ++i) {
            if (auto f = std::dynamic_pointer_cast<const FunctionNode>(mul->operands()[i])) {
                if (f->type() == FunctionNode::FuncType::Exp && f->arguments().size() == 1) {

                    auto exp_arg = f->arguments()[0];
                    std::vector<std::shared_ptr<const SymbolicNode>> remaining;
                    for (size_t j = 0; j < mul->operands().size(); ++j) {
                        if (j != i) remaining.push_back(mul->operands()[j]);
                    }
                    std::shared_ptr<const SymbolicNode> g_node;
                    if (remaining.size() == 1) {
                        g_node = remaining[0];
                    } else {
                        g_node = lamina::detail::make_node<MultiplyNode>(remaining);
                    }

                    if (nodes_equal(g_node, exp_arg) && expression_depends_on_variable(g_node, var)) {
                        return lamina::detail::make_expression_ptr(g_node);
                    }
                }
            }
        }
        return std::nullopt;
    };

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        auto g = check_mul_is_g_exp_g(mul);
        if (g) {
            return LambertWPattern{*g, SymbolicExpr::number(0)};
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
        std::optional<std::shared_ptr<SymbolicExpr>> g_found;
        std::vector<std::shared_ptr<const SymbolicNode>> const_terms;
        bool has_other_var_terms = false;

        for (auto& op : add->operands()) {
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
                auto g = check_mul_is_g_exp_g(mul);
                if (g && !g_found) {
                    g_found = g;
                    continue;
                }
            }
            if (!expression_depends_on_variable(op, var)) {
                const_terms.push_back(op);
            } else {
                has_other_var_terms = true;
            }
        }

        if (g_found && !has_other_var_terms) {
            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = lamina::detail::make_expression_ptr(const_terms[0]);
            } else {
                const_sum = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(const_terms));
            }
            auto rhs = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1))->simplify();
            return LambertWPattern{*g_found, rhs};
        }
    }

    return std::nullopt;
}

struct ExpBasePattern {
    std::shared_ptr<SymbolicExpr> base;
    std::shared_ptr<SymbolicExpr> inner;
    std::shared_ptr<SymbolicExpr> rhs;
};

static std::optional<ExpBasePattern> decompose_exp_base_pattern(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return std::nullopt;

    auto check_power_is_a_to_g = [&](const std::shared_ptr<const PowerNode>& pow)
        -> std::optional<std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>> {

        if (!expression_depends_on_variable(pow->base(), var) && expression_depends_on_variable(pow->exponent(), var)) {
            auto base_expr = lamina::detail::make_expression_ptr(pow->base());
            auto exp_expr = lamina::detail::make_expression_ptr(pow->exponent());

            lmmc_real_t base_val;
            if (try_evaluate_numeric(base_expr, base_val)) {
                if (base_val > 0 && std::abs(base_val - 1.0) > LMMC_REAL_EPSILON) {
                    return std::make_pair(base_expr, exp_expr);
                }
            }
        }
        return std::nullopt;
    };

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        auto result = check_power_is_a_to_g(pow);
        if (result) {

            return ExpBasePattern{result->first, result->second, SymbolicExpr::number(0)};
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
        std::optional<std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>> found;
        std::vector<std::shared_ptr<const SymbolicNode>> const_terms;
        bool has_other_var_terms = false;

        for (auto& op : add->operands()) {
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
                auto result = check_power_is_a_to_g(pow);
                if (result && !found) {
                    found = result;
                    continue;
                }
            }
            if (!expression_depends_on_variable(op, var)) {
                const_terms.push_back(op);
            } else {
                has_other_var_terms = true;
            }
        }

        if (found && !has_other_var_terms) {
            std::shared_ptr<SymbolicExpr> const_sum;
            if (const_terms.empty()) {
                const_sum = SymbolicExpr::number(0);
            } else if (const_terms.size() == 1) {
                const_sum = lamina::detail::make_expression_ptr(const_terms[0]);
            } else {
                const_sum = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(const_terms));
            }
            auto rhs = SymbolicExpr::multiply(const_sum, SymbolicExpr::number(-1))->simplify();
            return ExpBasePattern{found->first, found->second, rhs};
        }
    }

    return std::nullopt;
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_sin(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {

    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (std::abs(c_val) > 1.0 + LMMC_REAL_EPSILON) {
            return {};
        }
    }

    auto pi = make_pi();
    auto asin_c = make_arcsin(c);

    int num_periods = 1;
    if (max_roots > 0) {

        num_periods = (max_roots + 1) / 2;
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {

        auto two_k_pi = SymbolicExpr::multiply(SymbolicExpr::number(2 * k), pi);
        auto sol1 = SymbolicExpr::add(asin_c, two_k_pi)->simplify();
        results.push_back(sol1);

        auto pi_minus_asin = SymbolicExpr::add(pi, SymbolicExpr::multiply(asin_c, SymbolicExpr::number(-1)));
        auto sol2 = SymbolicExpr::add(pi_minus_asin, two_k_pi)->simplify();
        results.push_back(sol2);
    }

    return results;
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_cos(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {

    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (std::abs(c_val) > 1.0 + LMMC_REAL_EPSILON) {
            return {};
        }
    }

    auto pi = make_pi();
    auto acos_c = make_arccos(c);

    int num_periods = 1;
    if (max_roots > 0) {
        num_periods = (max_roots + 1) / 2;
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {
        auto two_k_pi = SymbolicExpr::multiply(SymbolicExpr::number(2 * k), pi);

        auto sol1 = SymbolicExpr::add(acos_c, two_k_pi)->simplify();
        results.push_back(sol1);

        auto neg_acos = SymbolicExpr::multiply(acos_c, SymbolicExpr::number(-1));
        auto sol2 = SymbolicExpr::add(neg_acos, two_k_pi)->simplify();
        results.push_back(sol2);
    }

    return results;
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_tan(
    const std::shared_ptr<SymbolicExpr>& c, int max_roots) {

    auto pi = make_pi();
    auto atan_c = make_arctan(c);

    int num_periods = 1;
    if (max_roots > 0) {
        num_periods = max_roots;
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (int k = 0; k < num_periods; ++k) {

        auto k_pi = SymbolicExpr::multiply(SymbolicExpr::number(k), pi);
        auto sol = SymbolicExpr::add(atan_c, k_pi)->simplify();
        results.push_back(sol);
    }

    return results;
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_exp(
    const std::shared_ptr<SymbolicExpr>& c) {

    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (c_val <= 0.0) {
            return {};
        }
    }

    auto sol = SymbolicExpr::ln(c)->simplify();
    return {sol};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_ln(
    const std::shared_ptr<SymbolicExpr>& c) {

    auto sol = SymbolicExpr::exp(c)->simplify();
    return {sol};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_lambert_w(
    const std::shared_ptr<SymbolicExpr>& c) {

    auto sol = SymbolicExpr::lambertw(c);
    if (!sol) return {};
    return {sol->simplify()};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_exp_base(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::shared_ptr<SymbolicExpr>& c) {

    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (c_val <= 0.0) {
            return {};
        }
    }

    auto ln_c = SymbolicExpr::ln(c);
    auto ln_a = SymbolicExpr::ln(base);
    auto sol = SymbolicExpr::divide(ln_c, ln_a)->simplify();
    return {sol};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_exp(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string&) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (u_val <= 0.0) return {};
    }
    auto x_val = SymbolicExpr::ln(u_root)->simplify();
    return {x_val};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_sin(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string&) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (std::abs(u_val) > 1.0 + LMMC_REAL_EPSILON) return {};
    }

    auto x_val = make_arcsin(u_root)->simplify();
    return {x_val};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_cos(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string&) {
    lmmc_real_t u_val;
    if (try_evaluate_numeric(u_root, u_val)) {
        if (std::abs(u_val) > 1.0 + LMMC_REAL_EPSILON) return {};
    }
    auto x_val = make_arccos(u_root)->simplify();
    return {x_val};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_tan(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string&) {
    auto x_val = make_arctan(u_root)->simplify();
    return {x_val};
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_h_power(
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string&,
    int n) {
    if (n <= 0) return {};

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    auto exponent = SymbolicExpr::divide(SymbolicExpr::number(1), SymbolicExpr::number(n));

    if (n % 2 == 0) {

        lmmc_real_t u_val;
        if (try_evaluate_numeric(u_root, u_val)) {
            if (u_val < 0.0) return {};
        }

        auto pos_root = SymbolicExpr::power(u_root, exponent)->simplify();
        auto neg_root = SymbolicExpr::multiply(pos_root, SymbolicExpr::number(-1))->simplify();
        results.push_back(pos_root);
        results.push_back(neg_root);
    } else {

        auto root_val = SymbolicExpr::power(u_root, exponent)->simplify();
        results.push_back(root_val);
    }

    return results;
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_polynomial_values(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const std::string& variable) {
    return solve_by_factoring(polynomial, variable);
}

static std::vector<std::shared_ptr<SymbolicExpr>> invert_substitution_h(
    const std::shared_ptr<SymbolicExpr>& h_expr,
    const std::shared_ptr<SymbolicExpr>& u_root,
    const std::string& var) {

    if (!h_expr || !lamina::detail::node(h_expr)) return {};

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(h_expr))) {
        if (func->arguments().size() == 1) {

            if (auto v = std::dynamic_pointer_cast<const VariableNode>(func->arguments()[0])) {
                if (v->name() == var) {
                    switch (func->type()) {
                        case FunctionNode::FuncType::Exp:
                            return invert_h_exp(u_root, var);
                        case FunctionNode::FuncType::Sin:
                            return invert_h_sin(u_root, var);
                        case FunctionNode::FuncType::Cos:
                            return invert_h_cos(u_root, var);
                        case FunctionNode::FuncType::Tan:
                            return invert_h_tan(u_root, var);
                        default:
                            break;
                    }
                }
            }

        }
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(h_expr))) {
        if (v->name() == var) {

            return {u_root};
        }
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(h_expr))) {
        if (auto v = std::dynamic_pointer_cast<const VariableNode>(pow->base())) {
            if (v->name() == var && !expression_depends_on_variable(pow->exponent(), var)) {
                if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                    int n = 0;
                    if (std::holds_alternative<BigInt>(exp_num->value()))
                        n = std::get<BigInt>(exp_num->value()).to_int();
                    else if (std::holds_alternative<double>(exp_num->value()))
                        n = (int)std::get<double>(exp_num->value());
                    else if (std::holds_alternative<Rational>(exp_num->value()))
                        n = (int)std::get<Rational>(exp_num->value()).to_double();
                    if (n >= 2) {
                        return invert_h_power(u_root, var, n);
                    }
                }
            }
        }
    }

    auto eq = SymbolicExpr::add(h_expr,
        SymbolicExpr::multiply(u_root, SymbolicExpr::number(-1)))->simplify();

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(eq, var);
    if (!poly.is_zero() && poly.degree() >= 1) {
        return solve_polynomial_values(poly, var);
    }

    return {};
}

static constexpr int MAX_TRANSCENDENTAL_DEPTH = 5;

static std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental_impl(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    int depth);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_inner_equation(
    const std::shared_ptr<SymbolicExpr>& inner,
    const std::shared_ptr<SymbolicExpr>& value,
    const std::string& var,
    int depth) {

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(inner))) {
        if (v->name() == var) {
            return {value};
        }
    }

    auto inner_eq = SymbolicExpr::add(inner,
        SymbolicExpr::multiply(value, SymbolicExpr::number(-1)))->simplify();

    auto inner_poly = symbolic_to_poly<SymbolicPolyCoeff>(inner_eq, var);
    if (!inner_poly.is_zero() && inner_poly.degree() >= 1) {
        auto inner_solutions = solve_polynomial_values(inner_poly, var);
        if (!inner_solutions.empty()) return inner_solutions;
    }

    return solve_transcendental_impl(inner_eq, var, depth + 1);
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental_impl(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    int depth) {

    if (!expr || !lamina::detail::node(expr)) return {};

    if (depth > MAX_TRANSCENDENTAL_DEPTH) return {};

    SolveOptions opts;
    int max_roots = opts.max_roots;

    auto lambert_pat = decompose_lambert_w_pattern(expr, var);
    if (lambert_pat) {
        auto inverted = invert_lambert_w(lambert_pat->rhs);
        if (inverted.empty()) return {};

        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(lambert_pat->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }

    auto exp_base_pat = decompose_exp_base_pattern(expr, var);
    if (exp_base_pat) {
        auto inverted = invert_exp_base(exp_base_pat->base, exp_base_pat->rhs);
        if (inverted.empty()) return {};

        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(exp_base_pat->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }

    auto pattern = decompose_trig_exp_pattern(expr, var);
    if (pattern) {
        std::vector<std::shared_ptr<SymbolicExpr>> inverted;

        switch (pattern->func_type) {
            case FunctionNode::FuncType::Sin:
                inverted = invert_sin(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Cos:
                inverted = invert_cos(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Tan:
                inverted = invert_tan(pattern->rhs, max_roots);
                break;
            case FunctionNode::FuncType::Exp:
                inverted = invert_exp(pattern->rhs);
                break;
            case FunctionNode::FuncType::Ln:
                inverted = invert_ln(pattern->rhs);
                break;
            default:
                return {};
        }

        if (inverted.empty()) return {};

        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& val : inverted) {
            auto inner_solutions = solve_inner_equation(pattern->inner, val, var, depth);
            results.insert(results.end(), inner_solutions.begin(), inner_solutions.end());
        }
        return results;
    }

    auto subst = detect_substitution(expr, var);
    if (subst) {

        auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(
            subst->poly_in_u, subst->u_var);
        auto u_solutions = solve_polynomial_values(
            polynomial, subst->u_var);
        if (u_solutions.empty()) return {};

        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (auto& u_root : u_solutions) {

            auto x_values = invert_substitution_h(subst->u_expr, u_root, var);
            results.insert(results.end(), x_values.begin(), x_values.end());
        }
        return results;
    }

    return {};
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    return solve_transcendental_impl(expr, var, 0);
}

static void collect_transcendental_subexprs(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var,
    std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {

    if (!node) return;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1 && expression_depends_on_variable(func->arguments()[0], var)) {
            auto ft = func->type();
            if (ft == FunctionNode::FuncType::Exp ||
                ft == FunctionNode::FuncType::Sin ||
                ft == FunctionNode::FuncType::Cos ||
                ft == FunctionNode::FuncType::Tan) {
                candidates.push_back(lamina::detail::make_expression_ptr(node));
            }
        }

        for (auto& arg : func->arguments()) {
            collect_transcendental_subexprs(arg, var, candidates);
        }
        return;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (expression_depends_on_variable(pow->base(), var) && !expression_depends_on_variable(pow->exponent(), var)) {
            if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                int e_val = 0;
                if (std::holds_alternative<BigInt>(exp_num->value()))
                    e_val = std::get<BigInt>(exp_num->value()).to_int();
                else if (std::holds_alternative<double>(exp_num->value()))
                    e_val = (int)std::get<double>(exp_num->value());
                else if (std::holds_alternative<Rational>(exp_num->value()))
                    e_val = (int)std::get<Rational>(exp_num->value()).to_double();

                if (e_val >= 2) {

                    auto base_expr = lamina::detail::make_expression_ptr(pow->base());
                    candidates.push_back(base_expr);
                }
            }
        }

        collect_transcendental_subexprs(pow->base(), var, candidates);
        collect_transcendental_subexprs(pow->exponent(), var, candidates);
        return;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }
}

static void deduplicate_candidates(std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {
    std::vector<std::shared_ptr<SymbolicExpr>> unique;
    std::vector<std::string> seen;
    for (auto& c : candidates) {
        std::string s = c->to_string();
        bool found = false;
        for (auto& existing : seen) {
            if (existing == s) { found = true; break; }
        }
        if (!found) {
            seen.push_back(s);
            unique.push_back(c);
        }
    }
    candidates = std::move(unique);
}

static int extract_exp_multiplier(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!func || func->type() != FunctionNode::FuncType::Exp || func->arguments().size() != 1)
        return 0;

    auto arg = func->arguments()[0];

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(arg)) {
        if (v->name() == var) return 1;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(arg)) {

        std::shared_ptr<const NumberNode> num_part = nullptr;
        bool has_var = false;
        bool has_other = false;

        for (auto& op : mul->operands()) {
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (!num_part) num_part = n;
                else has_other = true;
            } else if (auto v = std::dynamic_pointer_cast<const VariableNode>(op)) {
                if (v->name() == var) has_var = true;
                else has_other = true;
            } else {
                has_other = true;
            }
        }

        if (has_var && num_part && !has_other) {
            int k = 0;
            if (std::holds_alternative<BigInt>(num_part->value()))
                k = std::get<BigInt>(num_part->value()).to_int();
            else if (std::holds_alternative<double>(num_part->value()))
                k = (int)std::get<double>(num_part->value());
            else if (std::holds_alternative<Rational>(num_part->value()))
                k = (int)std::get<Rational>(num_part->value()).to_double();
            if (k > 0) return k;
        }
    }

    return 0;
}

static std::shared_ptr<const SymbolicNode> rewrite_exp_as_u_power(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var,
    const std::string& u_var) {

    if (!node) return node;

    int k = extract_exp_multiplier(node, var);
    if (k > 0) {
        auto u_node = lamina::detail::make_node<VariableNode>(u_var);
        if (k == 1) return u_node;
        return SymbolicFactory::create_power(u_node, lamina::detail::make_node<NumberNode>(BigInt(k)));
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(add->operands().size());
        for (auto& op : add->operands()) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_add(std::move(new_ops));
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (auto& op : mul->operands()) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_multiply(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = rewrite_exp_as_u_power(pow->base(), var, u_var);
        auto new_exp = rewrite_exp_as_u_power(pow->exponent(), var, u_var);
        return SymbolicFactory::create_power(new_base, new_exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (auto& arg : func->arguments()) {
            new_args.push_back(rewrite_exp_as_u_power(arg, var, u_var));
        }
        return lamina::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    return node;
}

static bool has_exp_k_var_terms(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return false;

    if (extract_exp_multiplier(node, var) > 0) return true;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands())
            if (has_exp_k_var_terms(op, var)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands())
            if (has_exp_k_var_terms(op, var)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (has_exp_k_var_terms(pow->base(), var)) return true;
        if (has_exp_k_var_terms(pow->exponent(), var)) return true;
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (auto& arg : func->arguments())
            if (has_exp_k_var_terms(arg, var)) return true;
    }

    return false;
}

std::optional<SubstitutionResult> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return std::nullopt;

    if (!expression_depends_on_variable(lamina::detail::node(expr), var)) return std::nullopt;

    const std::string u_var = "_u_subst";

    std::vector<std::shared_ptr<SymbolicExpr>> candidates;
    collect_transcendental_subexprs(lamina::detail::node(expr), var, candidates);
    deduplicate_candidates(candidates);

    for (auto& h : candidates) {

        if (auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(h))) {
            if (v->name() == var) {

                auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, var);
                if (!poly.is_zero() && poly.degree() >= 2) {

                    continue;
                }
            }
        }

        auto u_expr = SymbolicExpr::variable(u_var);

        struct SubstNodeVisitor {
            std::shared_ptr<const SymbolicNode> target;
            std::string u_name;

            std::shared_ptr<const SymbolicNode> replace(const std::shared_ptr<const SymbolicNode>& node) {
                if (!node) return node;

                if (node->equals(*target)) {
                    return lamina::detail::make_node<VariableNode>(u_name);
                }

                if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
                    new_ops.reserve(add->operands().size());
                    for (auto& op : add->operands()) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_add(std::move(new_ops));
                }
                if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
                    new_ops.reserve(mul->operands().size());
                    for (auto& op : mul->operands()) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_multiply(std::move(new_ops));
                }
                if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
                    auto new_base = replace(pow->base());
                    auto new_exp = replace(pow->exponent());
                    return SymbolicFactory::create_power(new_base, new_exp);
                }
                if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_args;
                    new_args.reserve(func->arguments().size());
                    for (auto& arg : func->arguments()) {
                        new_args.push_back(replace(arg));
                    }
                    return lamina::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
                }

                return node;
            }
        };

        SubstNodeVisitor visitor;
        visitor.target = lamina::detail::node(h);
        visitor.u_name = u_var;

        auto substituted_node = visitor.replace(lamina::detail::node(expr));
        auto substituted_expr = lamina::detail::make_expression_ptr(substituted_node);
        substituted_expr = substituted_expr->simplify();

        if (!expression_depends_on_variable(lamina::detail::node(substituted_expr), var)) {

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(substituted_expr, u_var);
            if (!poly.is_zero() && poly.degree() >= 2) {
                return SubstitutionResult{h, substituted_expr, u_var};
            }
        }
    }

    if (has_exp_k_var_terms(lamina::detail::node(expr), var)) {
        auto rewritten_node = rewrite_exp_as_u_power(lamina::detail::node(expr), var, u_var);
        auto rewritten_expr = lamina::detail::make_expression_ptr(rewritten_node);
        rewritten_expr = rewritten_expr->simplify();

        if (!expression_depends_on_variable(lamina::detail::node(rewritten_expr), var)) {

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(rewritten_expr, u_var);
            if (!poly.is_zero() && poly.degree() >= 2) {

                auto exp_x = SymbolicExpr::exp(SymbolicExpr::variable(var));
                return SubstitutionResult{exp_x, rewritten_expr, u_var};
            }
        }
    }

    return std::nullopt;
}

}
