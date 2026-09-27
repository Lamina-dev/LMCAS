#include "internal/inequality_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>

namespace LMCAS {
using namespace detail::inequality_support;
namespace {

int parametric_number_sign(const NumberNode& number) {
    const auto& value = number.value();
    if (const auto* integer = std::get_if<BigInt>(&value)) {
        if (integer->is_negative()) {
            return -1;
        }
        return integer->is_zero() ? 0 : 1;
    }
    if (const auto* rational = std::get_if<Rational>(&value)) {
        if (rational->get_numerator().is_negative()) {
            return -1;
        }
        return rational->get_numerator().is_zero() ? 0 : 1;
    }
    if (const auto* real = std::get_if<lmmc_real_t>(&value)) {
        if (*real < 0) {
            return -1;
        }
        return *real == 0 ? 0 : 1;
    }
    return 1;
}

bool is_parametric_radical(const std::shared_ptr<const SymbolicNode>& node) {
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = LMCAS::detail::make_expression_ptr(power->exponent());
        auto value = try_checked_numeric_constant(*exponent);
        return value && *value > 0 && *value < 1.0;
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return function->type() == FunctionNode::FuncType::Sqrt;
    }
    return false;
}

bool parametric_difference_has_sign(const std::shared_ptr<SymbolicExpr>& expr, int wanted) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    const auto& node = LMCAS::detail::node(expr);
    if (wanted > 0 && is_parametric_radical(node)) {
        return true;
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return false;
    }
    int sign = 1;
    for (const auto& operand : product->operands()) {
        if (auto number = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            int factor_sign = parametric_number_sign(*number);
            if (factor_sign == 0) {
                return false;
            }
            sign *= factor_sign;
        } else if (!is_parametric_radical(operand)) {
            return false;
        }
    }
    return wanted > 0 ? sign > 0 : sign < 0;
}

int parametric_leading_sign(const std::shared_ptr<SymbolicExpr>& coefficient) {
    if (!coefficient) {
        return 0;
    }
    int sign = exact_numeric_sign(coefficient);
    if (sign == 0 && !coefficient->is_zero()) {
        if (auto value = try_checked_numeric_constant(*coefficient)) {
            sign = *value > 0 ? 1 : (*value < 0 ? -1 : 0);
        }
    }
    return sign;
}

void order_parametric_roots(std::vector<std::shared_ptr<SymbolicExpr>>& roots,
                            int leading_sign, int degree) {
    if (roots.size() != 2) {
        std::sort(roots.begin(), roots.end(), root_less_than);
        return;
    }
    bool swapped = false;
    auto difference = SymbolicExpr::add(roots[0],
        SymbolicExpr::multiply(roots[1], SymbolicExpr::number(-1)))->simplify();
    if (!difference) {
        std::sort(roots.begin(), roots.end(), root_less_than);
        swapped = true;
    } else if (auto value = try_checked_numeric_constant(*difference)) {
        if (*value > 0) {
            std::swap(roots[0], roots[1]);
            swapped = true;
        }
    } else if (parametric_difference_has_sign(difference, 1)) {
        std::swap(roots[0], roots[1]);
        swapped = true;
    } else if (parametric_difference_has_sign(difference, -1)) {
        swapped = false;
    }
    if (!swapped && leading_sign > 0 && degree == 2) {
        std::swap(roots[0], roots[1]);
    }
}

std::shared_ptr<SymbolicExpr> leading_condition(
    const std::shared_ptr<SymbolicExpr>& coefficient, RelationalNode::Op relation) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(coefficient),
            LMCAS::detail::node(SymbolicExpr::number(0)), relation));
}

PiecewiseIntervalResult::Case signed_parametric_case(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::shared_ptr<SymbolicExpr>& leading_coeff,
    const std::string& variable, InequalityType type, int sign) {
    PiecewiseIntervalResult::Case branch;
    branch.condition = leading_condition(leading_coeff,
        sign > 0 ? RelationalNode::Op::GT : RelationalNode::Op::LT);
    auto symbolic_roots = solve_symbolic_poly(poly, variable);
    std::vector<int> multiplicities(symbolic_roots.size(), 1);
    std::vector<size_t> indices(symbolic_roots.size());
    for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return root_less_than(symbolic_roots[a], symbolic_roots[b]);
    });
    std::vector<std::shared_ptr<SymbolicExpr>> sorted_roots;
    std::vector<int> sorted_mults;
    for (size_t idx : indices) {
        sorted_roots.push_back(symbolic_roots[idx]);
        sorted_mults.push_back(multiplicities[idx]);
    }
    branch.solution = build_parametric_intervals(sorted_roots, sorted_mults, sign, type);
    return branch;
}

std::shared_ptr<SymbolicExpr> reduced_polynomial_expression(
    const Polynomial<SymbolicPolyCoeff>& reduced_poly, const std::string& variable) {
    auto reduced_expr = SymbolicExpr::number(0);
    auto var_expr = SymbolicExpr::variable(variable);
    for (int i = reduced_poly.degree(); i >= 0; --i) {
        auto coeff_val = reduced_poly.coeffs[i].val;
        if (!coeff_val) {
            continue;
        }
        if (i == 0) {
            reduced_expr = SymbolicExpr::add(reduced_expr, coeff_val);
        } else if (i == 1) {
            reduced_expr = SymbolicExpr::add(reduced_expr,
                SymbolicExpr::multiply(coeff_val, var_expr));
        } else {
            reduced_expr = SymbolicExpr::add(reduced_expr,
                SymbolicExpr::multiply(coeff_val,
                    SymbolicExpr::power(var_expr, SymbolicExpr::number(i))));
        }
    }
    reduced_expr = reduced_expr->simplify();
    return reduced_expr;
}

bool append_degenerate_subcases(
    PiecewiseIntervalResult& result, PiecewiseIntervalResult::Case& parent,
    const Polynomial<SymbolicPolyCoeff>& reduced_poly, InequalityType type,
    const std::string& variable, const std::vector<std::string>& parameters) {
    auto reduced_expr = reduced_polynomial_expression(reduced_poly, variable);
    auto sub_result = InequalitySolver::solve_parametric_inequality(reduced_expr, type, variable, parameters);
    if (sub_result.cases.empty()) {
        parent.solution = IntervalUnion::empty();
        return false;
    }
    for (const auto& sub_case : sub_result.cases) {
        PiecewiseIntervalResult::Case merged;
        if (sub_case.condition) {
            merged.condition = LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<LogicalNode>(
                    LMCAS::detail::node(parent.condition),
                    LMCAS::detail::node(sub_case.condition), LogicalNode::Op::And));
        } else {
            merged.condition = parent.condition;
        }
        merged.solution = sub_case.solution;
        result.cases.push_back(merged);
    }
    return true;
}

IntervalUnion zero_polynomial_solution(InequalityType type) {
    if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
        return IntervalUnion::entire_line();
    }
    return IntervalUnion::empty();
}

void append_degenerate_case(
    PiecewiseIntervalResult& result, const Polynomial<SymbolicPolyCoeff>& poly,
    const std::shared_ptr<SymbolicExpr>& leading_coeff, InequalityType type,
    const std::string& variable, const std::vector<std::string>& parameters) {
    PiecewiseIntervalResult::Case branch;
    branch.condition = leading_condition(leading_coeff, RelationalNode::Op::EQ);
    int degree = poly.degree();
    if (degree < 1) {
        branch.solution = IntervalUnion::empty();
        result.cases.push_back(branch);
        return;
    }
    std::vector<SymbolicPolyCoeff> reduced_coeffs;
    for (int i = 0; i < degree; ++i) reduced_coeffs.push_back(poly.coeffs[i]);
    Polynomial<SymbolicPolyCoeff> reduced_poly(reduced_coeffs, variable);
    if (reduced_poly.is_zero()) {
        branch.solution = zero_polynomial_solution(type);
        result.cases.push_back(branch);
        return;
    }
    auto new_lc = reduced_poly.lead_coeff().val;
    if (new_lc) {
        new_lc = new_lc->simplify();
    }
    if (new_lc && depends_on_any_param(new_lc, parameters)) {
        if (append_degenerate_subcases(result, branch, reduced_poly, type, variable, parameters)) {
            return;
        }
    } else {
        int sign = parametric_leading_sign(new_lc);
        auto roots = solve_symbolic_poly(reduced_poly, variable);
        std::vector<int> multiplicities(roots.size(), 1);
        branch.solution = build_parametric_intervals(roots, multiplicities, sign, type);
    }
    result.cases.push_back(branch);
}

}

PiecewiseIntervalResult InequalitySolver::solve_parametric_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters) {

    PiecewiseIntervalResult result;

    if (!expr) {
        return result;
    }

    if (parameters.empty()) {
        auto solution = solve_inequality(expr, type, variable);
        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;
        single_case.solution = solution;
        result.cases.push_back(single_case);
        return result;
    }

    auto converted = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (!converted) throw std::invalid_argument(converted.error().message);
    const auto& poly = converted.value();

    if (poly.is_zero()) {
        PiecewiseIntervalResult::Case zero_case;
        zero_case.condition = nullptr;
        zero_case.solution = zero_polynomial_solution(type);
        result.cases.push_back(zero_case);
        return result;
    }
    int deg = poly.degree();
    auto leading_coeff = poly.coeffs[deg].val;
    if (!leading_coeff) {
        leading_coeff = SymbolicExpr::number(0);
    }
    leading_coeff = leading_coeff->simplify();

    bool lc_depends_on_params = depends_on_any_param(leading_coeff, parameters);

    if (!lc_depends_on_params) {
        int sign = parametric_leading_sign(leading_coeff);
        auto roots = solve_symbolic_poly(poly, variable);
        order_parametric_roots(roots, sign, deg);
        std::vector<int> multiplicities(roots.size(), 1);
        PiecewiseIntervalResult::Case branch;
        branch.condition = nullptr;
        branch.solution = build_parametric_solution(roots, multiplicities, sign, type);
        result.cases.push_back(branch);
    } else {
        result.cases.push_back(signed_parametric_case(poly, leading_coeff, variable, type, 1));
        result.cases.push_back(signed_parametric_case(poly, leading_coeff, variable, type, -1));
        append_degenerate_case(result, poly, leading_coeff, type, variable, parameters);
    }
    return result;
}

}
