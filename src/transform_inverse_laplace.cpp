#include "internal/transform_support.hpp"
#include "internal/normalization_utils.hpp"
#include "residual_verification.hpp"
#include <cmath>
#include <exception>

namespace LMCAS {
using namespace transform_detail;

static std::shared_ptr<SymbolicExpr> te_unevaluated_inv_laplace(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s, const std::string& t) {
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::InverseLaplace, LMCAS::detail::node(F)->clone(), s,
        SymbolicFactory::create_variable(t)));
}
static std::shared_ptr<SymbolicExpr> te_inverse_pole(
    const std::shared_ptr<const SymbolicNode>& base, const std::string& variable) {
    auto addition = std::dynamic_pointer_cast<const AddNode>(base);
    if (!addition || addition->operands().size() != 2) return nullptr;
    for (std::size_t i = 0; i < 2; ++i) {
        auto symbol = std::dynamic_pointer_cast<const VariableNode>(addition->operands()[i]);
        if (symbol && !symbol->is_constant() && symbol->name() == variable) {
            return detail::make_expression_ptr(SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), addition->operands()[1 - i]}));
        }
    }
    return nullptr;
}

static std::shared_ptr<SymbolicExpr> te_inv_power(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s,
    const std::string& t, ComputationContext& context) {
    if (!F || !LMCAS::detail::node(F)) {
        return nullptr;
    }
    auto tv = SymbolicExpr::variable(t);
    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(F));
    if (!pw) {
        return nullptr;
    }
    BigInt exponent;
    if (!try_get_integer_value(std::dynamic_pointer_cast<const NumberNode>(pw->exponent()), exponent) ||
        !exponent.is_negative()) {
        return nullptr;
    }
    const BigInt n = -exponent;
    auto make_exp = [&](const std::shared_ptr<SymbolicExpr>& a) {
        auto arg = SymbolicExpr::multiply(a, tv);
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(arg)}));
    };
    auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (bv && !bv->is_constant() && bv->name() == s) {
        if (n == BigInt(1)) {
            return SymbolicExpr::number(1);
        }
        return SymbolicExpr::divide(SymbolicExpr::power(tv, SymbolicExpr::number(n - BigInt(1))),
            SymbolicExpr::number(te_factorial(n - BigInt(1), context)));
    }
    auto pole = te_inverse_pole(pw->base(), s);
        if (pole) {
            if (n == BigInt(1)) {
                if (te_is_zero(pole)) {
                    return SymbolicExpr::number(1);
                }
                return make_exp(pole);
            }
            auto tp = SymbolicExpr::power(tv, SymbolicExpr::number(n - BigInt(1)));
            auto fv = SymbolicExpr::number(te_factorial(n - BigInt(1), context));
            if (te_is_zero(pole)) {
                return SymbolicExpr::divide(tp, fv);
            }
            return SymbolicExpr::divide(SymbolicExpr::multiply(tp, make_exp(pole)), fv);
        }
    return nullptr;
}
static bool te_is_denominator_factor(
    const std::shared_ptr<const SymbolicNode>& expression,
    const std::string& variable) {
    const auto power = std::dynamic_pointer_cast<const PowerNode>(expression);
    if (!power) return false;
    const auto exponent = te_number_value(power->exponent());
    return exponent && *exponent < 0 &&
        expression_depends_on_variable(power->base(), variable);
}

static std::shared_ptr<SymbolicExpr> te_inv_product(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s,
    const std::string& t, ComputationContext& context) {
    if (!F || !LMCAS::detail::node(F)) {
        return nullptr;
    }
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(F));
    if (!mul || mul->operands().size() < 2) {
        return nullptr;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> num_parts;
    std::shared_ptr<const SymbolicNode> den_part = nullptr;
    for (const auto& op : mul->operands()) {
        if (!den_part && te_is_denominator_factor(op, s)) {
            den_part = op;
            continue;
        }
        num_parts.push_back(op);
    }
    if (!den_part) {
        return nullptr;
    }
    std::shared_ptr<SymbolicExpr> num_expr;
    if (num_parts.empty()) {
        num_expr = SymbolicExpr::number(1);
    } else if (num_parts.size() == 1) {
        num_expr = LMCAS::detail::make_expression_ptr(num_parts[0]);
    } else {
        num_expr = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(num_parts));
    }
    if (!te_depends_on(num_expr, s)) {
        auto inv = te_inv_power(LMCAS::detail::make_expression_ptr(den_part), s, t, context);
        if (inv) {
            return SymbolicExpr::multiply(num_expr, inv);
        }
    }
    return nullptr;
}

static std::shared_ptr<SymbolicExpr> inverse_laplace_core(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s,
    const std::string& t, ComputationContext& context) {
    auto step = context.consume_steps(1, "inverse_laplace.recursive");
    if (!step) {
        return nullptr;
    }
    if (!F || !LMCAS::detail::node(F)) {
        return nullptr;
    }
    if (te_is_zero(F)) {
        return SymbolicExpr::number(0);
    }
    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(F));
    if (add) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add->operands()) {
            auto ilt = inverse_laplace_core(
                LMCAS::detail::make_expression_ptr(op), s, t, context);
            if (!ilt) {
                return te_unevaluated_inv_laplace(F, s, t);
            }
            result = result ? SymbolicExpr::add(result, ilt) : ilt;
        }
        return result;
    }
    auto [coeff, body] = te_split_coeff(F, s);
    if (!coeff->is_one()) {
        auto ilt_body = inverse_laplace_core(body, s, t, context);
        if (ilt_body) {
            return SymbolicExpr::multiply(coeff, ilt_body);
        }
    }
    auto pw_res = te_inv_power(F, s, t, context);
    if (pw_res) {
        return pw_res;
    }
    auto prod_res = te_inv_product(F, s, t, context);
    if (prod_res) {
        return prod_res;
    }
    return te_unevaluated_inv_laplace(F, s, t);
}

TransformEngineResult inverse_laplace_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& s,
    const std::string& t,
    ComputationContext& context) {
    const std::string operation = "inverse_laplace";
    auto valid = te_validate_expr_vars(F, s, t, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = inverse_laplace_core(F, s, t, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        auto result = te_wrap_transform_result(expression, operation);
        if (!result) return result;
        auto round_trip = laplace_transform_core(
            expression, t, s, context);
        if (!round_trip ||
            te_contains_transform(LMCAS::detail::node(round_trip))) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Laplace round trip is not proved", operation);
        }
        auto verified = check_equivalent(round_trip, F, context);
        if (!verified) return TransformEngineResult::failure(verified.error());
        if (!std::holds_alternative<ProvedZeroResidual>(verified.value())) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Laplace round trip is not proved", operation);
        }
        result.value().certificate =
            ExactRoundTripProof{SymbolicExpr::number(0)};
        return result;
    } catch (const CasError& error) {
        return TransformEngineResult::failure(error);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult inverse_laplace_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& s,
    const std::string& t) {
    ComputationContext context;
    return inverse_laplace_checked(F, s, t, context);
}

}
