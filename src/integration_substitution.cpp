#include "internal/integration_support.hpp"

namespace LMCAS {

namespace {
std::optional<SymbolicExpr> substitution_argument(const SymbolicExpr& term) {
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(term))) {
        return detail::expression_from_node(power->base());
    }
    auto function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(term));
    if (function && !function->arguments().empty()) {
        return detail::expression_from_node(function->arguments()[0]);
    }
    return std::nullopt;
}

std::shared_ptr<SymbolicExpr> substitution_primitive(
    const SymbolicExpr& candidate_term, const SymbolicExpr& u) {
                if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(candidate_term))) {
                    auto n = LMCAS::detail::expression_from_node(pow->exponent());
                    auto np1 = SymbolicExpr::add(detail::make_expression_ptr(n), SymbolicExpr::number(1))->simplify();
                    if (np1->is_zero()) {
                        return SymbolicExpr::ln(detail::make_expression_ptr(u));
                    } else {
                        return SymbolicExpr::divide(
                            SymbolicExpr::power(detail::make_expression_ptr(u), np1), np1);
                    }
                } else if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(candidate_term))) {
                    if (func->type() == FunctionNode::FuncType::Cos) {
                        return SymbolicExpr::sin(detail::make_expression_ptr(u));
                    } else if (func->type() == FunctionNode::FuncType::Sin) {
                        return SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(detail::make_expression_ptr(u)));
                    } else if (func->type() == FunctionNode::FuncType::Exp) {
                        return SymbolicExpr::exp(detail::make_expression_ptr(u));
                    }
                }
    return nullptr;
}

Result<std::shared_ptr<SymbolicExpr>> substitute_candidate(
    const SymbolicExpr& expr, const SymbolicExpr& candidate_term,
    const std::string& var, ComputationContext& computation) {
    auto u = substitution_argument(candidate_term);
    if (!u || !depends_on_integration_variable(*u, var)) { return nullptr; }
            auto d_ptr = u->differentiate(var);
            if (!d_ptr) { return nullptr; }
            auto du = d_ptr->simplify();
            if (du->is_zero()) { return nullptr; }

            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(detail::make_expression_ptr(f_u), detail::make_expression_ptr(*du));
            auto ratio = SymbolicExpr::divide(detail::make_expression_ptr(expr), term_times_du)->simplify();
            bool ratio_independent = !depends_on_integration_variable(*ratio, var);
            if (!ratio_independent && computation.assumptions()) {
                auto nonzero =
                    computation.assumptions()->is_nonzero(*du);
                if (!nonzero) {
                    return Result<std::shared_ptr<SymbolicExpr>>::failure(
                        nonzero.error());
                }
                if (nonzero.value() == Tribool::True) {
                    ratio_independent = !depends_on_integration_variable(*ratio, var);
                }
            }
    if (!ratio_independent) { return nullptr; }
    auto primitive = substitution_primitive(candidate_term, *u);
    if (!primitive) { return nullptr; }
    return SymbolicExpr::multiply(ratio, primitive);
}
}

Result<std::shared_ptr<SymbolicExpr>> SubstitutionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& computation, int) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(LMCAS::detail::node(expr));
    }
    for (const auto& operand : ops) {
        auto candidate = detail::expression_from_node(operand);
        auto result = substitute_candidate(expr, candidate, var, computation);
        if (!result || result.value()) { return result; }
    }
    return nullptr;
}

namespace {
enum class LinearWrapperKind { None, Function, PowerOfFunction, PowerOfLinear };
struct LinearWrapper {
    LinearWrapperKind kind = LinearWrapperKind::None;
    std::shared_ptr<const FunctionNode> function;
    std::shared_ptr<const PowerNode> power;
    std::optional<SymbolicExpr> argument;
};

bool linear_substitutable_function(const FunctionNode& function) {
    using FT = FunctionNode::FuncType;
    switch (function.type()) {
        case FT::Infinity:
        case FT::Atan2:
        case FT::Log: return false;
        default: return function.arguments().size() == 1;
    }
}

LinearWrapper classify_linear_wrapper(const SymbolicExpr& expr, const std::string& var) {
    LinearWrapper wrapper;
    wrapper.function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(expr));
    if (wrapper.function) {
        if (!linear_substitutable_function(*wrapper.function)) { return wrapper; }
        wrapper.argument = detail::expression_from_node(wrapper.function->arguments()[0]);
        wrapper.kind = LinearWrapperKind::Function;
        return wrapper;
    }
    wrapper.power = std::dynamic_pointer_cast<const PowerNode>(detail::node(expr));
    if (!wrapper.power) { return wrapper; }
    if (expression_depends_on_variable(wrapper.power->exponent(), var)) { return wrapper; }
    wrapper.function = std::dynamic_pointer_cast<const FunctionNode>(wrapper.power->base());
    if (wrapper.function) {
        if (!linear_substitutable_function(*wrapper.function)) { return wrapper; }
        wrapper.argument = detail::expression_from_node(wrapper.function->arguments()[0]);
        wrapper.kind = LinearWrapperKind::PowerOfFunction;
    } else {
        wrapper.argument = detail::expression_from_node(wrapper.power->base());
        wrapper.kind = LinearWrapperKind::PowerOfLinear;
    }
    return wrapper;
}

std::shared_ptr<SymbolicExpr> linear_dummy_expression(
    const LinearWrapper& wrapper, const std::shared_ptr<SymbolicExpr>& dummy_var) {
    if (wrapper.kind == LinearWrapperKind::PowerOfLinear) {
        return detail::make_expression_ptr(detail::make_node<PowerNode>(
            detail::node(dummy_var), wrapper.power->exponent()));
    }
    std::vector<std::shared_ptr<const SymbolicNode>> args{detail::node(dummy_var)};
    auto function = detail::make_node<FunctionNode>(wrapper.function->type(), args);
    if (wrapper.kind == LinearWrapperKind::Function) {
        return detail::make_expression_ptr(function);
    }
    return detail::make_expression_ptr(
        detail::make_node<PowerNode>(function, wrapper.power->exponent()));
}
}

Result<bool> LinearSubstitutionStrategy::extract_linear_arg(
    const SymbolicExpr& arg,
    const std::string& var,
    std::shared_ptr<SymbolicExpr>& a_out,
    std::shared_ptr<SymbolicExpr>& b_out) {
    if (!expression_depends_on_variable(LMCAS::detail::node(arg), var)) { return false; }

    auto converted = symbolic_to_poly<SymbolicPolyCoeff>(detail::make_expression_ptr(arg), var);
    if (!converted) {
        if (converted.error().code == CasErrc::UnsupportedExpression) { return false; }
        return Result<bool>::failure(converted.error());
    }
    auto poly = std::move(converted.value());
    if (poly.degree() != 1) { return false; }
    if (poly.coeffs.size() < 2) { return false; }

    auto a_expr = poly.coeffs[1].val;
    auto b_expr = poly.coeffs[0].val;
    if (!a_expr || !b_expr) { return false; }
    if (expression_depends_on_variable(LMCAS::detail::node(a_expr), var)) { return false; }
    if (expression_depends_on_variable(LMCAS::detail::node(b_expr), var)) { return false; }

    auto a_simp = a_expr->simplify();
    auto b_simp = b_expr->simplify();
    if (!a_simp) a_simp = a_expr;
    if (!b_simp) b_simp = b_expr;
    if (a_simp->is_zero()) { return false; }

    a_out = a_simp;
    b_out = b_simp;
    return true;
}

Result<std::shared_ptr<SymbolicExpr>> LinearSubstitutionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {
    auto wrapper = classify_linear_wrapper(expr, var);
    auto& arg_expr = wrapper.argument;
    std::shared_ptr<SymbolicExpr> a_coeff, b_coeff;
    if (!arg_expr) { return nullptr; }
    auto extracted = extract_linear_arg(*arg_expr, var, a_coeff, b_coeff);
    if (!extracted) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(extracted.error());
    }
    if (!extracted.value()) { return nullptr; }

    /**
     * @brief arg == var（a = 1，b = 0）由前序 TableLookup 或 PowerRule 处理。
     * 此处返回 nullptr，将控制权交给下一策略。
     */
    bool a_is_one  = a_coeff && a_coeff->is_one();
    bool b_is_zero = b_coeff && b_coeff->is_zero();
    if (a_is_one && b_is_zero) { return nullptr; }

    const std::string dummy_name = "__lin_sub_u__";
    auto dummy_var = SymbolicExpr::variable(dummy_name);
    auto test_expr = linear_dummy_expression(wrapper, dummy_var);
    if (!test_expr) { return nullptr; }

    TableLookupStrategy table_only;
    auto table_attempt = table_only.try_integrate(
        *test_expr, dummy_name, ctx, computation, depth);
    if (!table_attempt) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            table_attempt.error());
    }
    auto* table_candidate =
        std::get_if<IntegrationCandidate>(&table_attempt.value());
    if (!table_candidate || !table_candidate->expression) { return nullptr; }
    auto F_substituted = table_candidate->expression->substitute(
        dummy_name, detail::make_expression_ptr(*arg_expr));
    if (!F_substituted) { return nullptr; }

    auto inv_a = SymbolicExpr::power(a_coeff, SymbolicExpr::number(-1));
    auto result = SymbolicExpr::multiply(inv_a, F_substituted);

    auto simplified = result->simplify();
    if (simplified && !contains_unevaluated_integral(LMCAS::detail::node(simplified))) {
        return simplified;
    }
    return result;
}

}
