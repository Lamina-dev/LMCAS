#include "transcendental_factor.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/transcendental_support.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <cstdint>

namespace LMCAS {

static int tf_power_degree(const PowerNode& power, const std::string& var) {
    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    if (!exp_num) return -1;

    BigInt exponent(-1);
    if (std::holds_alternative<BigInt>(exp_num->value())) {
        exponent = std::get<BigInt>(exp_num->value());
    } else if (std::holds_alternative<Rational>(exp_num->value())) {
        const auto& r = std::get<Rational>(exp_num->value());
        if (r.is_integer()) exponent = r.to_bigint();
    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
        lmmc_real_t d = std::get<lmmc_real_t>(exp_num->value());
        if (std::isfinite(d) && d >= 0 && d == std::floor(d) && d < 1000.0) {
            exponent = Rational::from_double(d).to_bigint();
        }
    }

    if (exponent.is_negative() || exponent > BigInt(std::numeric_limits<int>::max())) {
        return -1;
    }
    const int e_val = static_cast<int>(exponent.try_to_int64().value());

    int base_deg = tf_degree_in(power.base(), var);
    if (base_deg < 0 ||
        (e_val != 0 && base_deg > std::numeric_limits<int>::max() / e_val)) return -1;
    return base_deg * e_val;
}

int tf_degree_in(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return 0;

    if (!expression_depends_on_variable(node, var)) return 0;

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? 1 : 0;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        int max_deg = 0;
        for (const auto& op : add->operands()) {
            int d = tf_degree_in(op, var);
            if (d < 0) return -1;
            max_deg = std::max(max_deg, d);
        }
        return max_deg;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_deg = 0;
        for (const auto& op : mul->operands()) {
            int d = tf_degree_in(op, var);
            if (d < 0) return -1;
            if (d > std::numeric_limits<int>::max() - total_deg) return -1;
            total_deg += d;
        }
        return total_deg;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return tf_power_degree(*power, var);
    }
    if (std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return -1;
    }

    return 0;
}

/**
 * @brief 判断换元式是否为候选变量的多项式。
 * 各候选变量的次数须可确定为非负整数，且超越函数不得依赖不定元或原始变量。
 * @param[in] poly_expr 换元后的表达式。
 * @param[in] all_variables 所有候选变量名，包括不定元和原始变量。
 * @return 有效多项式返回 true。
 * @internal
 */
static bool tf_validate_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& all_variables) {

    if (!poly_expr || !LMCAS::detail::node(poly_expr)) return false;

    for (const auto& var : all_variables) {
        if (expression_depends_on_variable(LMCAS::detail::node(poly_expr), var)) {
            int deg = tf_degree_in(LMCAS::detail::node(poly_expr), var);
            if (deg < 0) return false;
        }
    }
    return true;
}

static std::vector<std::string> tf_polynomial_variables(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::vector<std::string>& indeterminates, const std::string& original_var) {
    std::vector<std::string> variables;
    for (const auto& variable : indeterminates) {
        if (expression_depends_on_variable(detail::node(expression), variable)) {
            variables.push_back(variable);
        }
    }
    if (expression_depends_on_variable(detail::node(expression), original_var)) {
        variables.push_back(original_var);
    }
    return variables;
}

static void tf_select_main_variable(TfPolyBuildResult& result,
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& variables, int& max_degree) {
    for (const auto& variable : variables) {
        const int degree = tf_degree_in(detail::node(expanded), variable);
        if (degree > max_degree) {
            max_degree = degree;
            result.main_variable = variable;
        }
    }
    for (const auto& variable : variables) {
        if (variable != result.main_variable) result.param_variables.push_back(variable);
    }
}

static Result<TfPolyBuildResult> tf_convert_main_variable(
    TfPolyBuildResult result, const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& all_variables) {
    const auto main_var = result.main_variable;
    auto converted = symbolic_to_poly<Rational>(expanded, main_var);
    if (converted) {
        result.poly = std::move(converted.value());
        result.success = true;
        return result;
    }
    if (converted.error().code != CasErrc::UnsupportedExpression) {
        return Result<TfPolyBuildResult>::failure(converted.error());
    }
    for (const auto& var : all_variables) {
        if (var == main_var) continue;
        auto trial = symbolic_to_poly<Rational>(expanded, var);
        if (!trial) {
            if (trial.error().code == CasErrc::UnsupportedExpression) continue;
            return Result<TfPolyBuildResult>::failure(trial.error());
        }
        result.main_variable = var;
        result.poly = std::move(trial.value());
        result.param_variables.clear();
        for (const auto& other : all_variables) {
            if (other != var) result.param_variables.push_back(other);
        }
        result.success = true;
        return result;
    }
    result.success = false;
    return result;
}

Result<TfPolyBuildResult> tf_build_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& indeterminates,
    const std::string& original_var) try {

    TfPolyBuildResult result;
    result.success = false;

    if (!poly_expr || !LMCAS::detail::node(poly_expr) || original_var.empty()) {
        return Result<TfPolyBuildResult>::failure(CasErrc::InvalidArgument,
            "expression and variable must be nonempty", "tf_build_polynomial");
    }

    auto before = tf_polynomial_variables(poly_expr, indeterminates, original_var);
    if (!before.empty() && !tf_validate_polynomial(poly_expr, before)) return result;
    auto expanded = poly_expr->expand();
    if (!expanded || !LMCAS::detail::node(expanded)) {
        expanded = poly_expr;
    }
    auto all_variables = tf_polynomial_variables(expanded, indeterminates, original_var);
    if (all_variables.empty()) {
        auto c = extract_coeff_value<Rational>(expanded);
        if (!c) {
            if (c.error().code == CasErrc::UnsupportedExpression) return result;
            return Result<TfPolyBuildResult>::failure(c.error());
        }
        result.poly = Polynomial<Rational>({std::move(c.value())}, "x");
        result.main_variable = "x";
        result.success = true;
        return result;
    }

    if (!tf_validate_polynomial(expanded, all_variables)) {
        return result;
    }
    int max_degree = -1;
    tf_select_main_variable(result, expanded, all_variables, max_degree);
    if (result.main_variable.empty() || max_degree < 0) return result;
    return tf_convert_main_variable(std::move(result), expanded, all_variables);
}
catch (const std::bad_alloc&) {
    return Result<TfPolyBuildResult>::failure(CasErrc::ResourceLimit,
        "allocation failed while constructing polynomial", "tf_build_polynomial");
} catch (const std::exception& ex) {
    return Result<TfPolyBuildResult>::failure(CasErrc::InternalInvariant,
        ex.what(), "tf_build_polynomial");
}

}
