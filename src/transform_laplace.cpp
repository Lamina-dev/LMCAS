#include "internal/transform_support.hpp"
#include "internal/normalization_utils.hpp"
#include <exception>

namespace LMCAS {
using namespace transform_detail;

static std::shared_ptr<SymbolicExpr> te_unevaluated_laplace(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t, const std::string& s) {
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::Laplace, LMCAS::detail::node(f)->clone(), t,
        SymbolicFactory::create_variable(s)));
}
static bool te_is_power_of_var(const std::shared_ptr<SymbolicExpr>& e, const std::string& v, BigInt& n) {
    if (!e || !LMCAS::detail::node(e)) {
        return false;
    }
    auto vn = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(e));
    if (vn && !vn->is_constant() && vn->name() == v) { n = BigInt(1); return true; }
    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(e));
    if (!pw) {
        return false;
    }
    auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (!bv || bv->is_constant() || bv->name() != v) {
        return false;
    }
    return try_get_integer_value(std::dynamic_pointer_cast<const NumberNode>(pw->exponent()), n) &&
           !n.is_negative();
}
static bool te_is_hyp(const std::shared_ptr<SymbolicExpr>& e, const std::string& v,
                      bool& is_sinh, std::shared_ptr<SymbolicExpr>& freq) {
    if (!e || !LMCAS::detail::node(e)) {
        return false;
    }
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(e));
    if (!fn || fn->arguments().empty()) {
        return false;
    }
    if (fn->type() != FunctionNode::FuncType::Sinh && fn->type() != FunctionNode::FuncType::Cosh) {
        return false;
    }
    is_sinh = (fn->type() == FunctionNode::FuncType::Sinh);
    freq = te_linear_coefficient(fn->arguments()[0], v);
    return freq != nullptr;
}

static bool te_extract_exp(const std::shared_ptr<SymbolicExpr>& e, const std::string& t,
                           std::shared_ptr<SymbolicExpr>& a_out, std::shared_ptr<SymbolicExpr>& rem_out) {
    if (!e || !LMCAS::detail::node(e)) {
        return false;
    }
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(e));
    if (fn && fn->type() == FunctionNode::FuncType::Exp && !fn->arguments().empty()) {
        a_out = te_linear_coefficient(fn->arguments()[0], t);
        if (a_out) {
            rem_out = detail::make_expression_ptr(SymbolicFactory::create_number(BigInt(1)));
            return true;
        }
    }
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(e));
    if (mul) {
        for (size_t i = 0; i < mul->operands().size(); ++i) {
            auto fac = LMCAS::detail::make_expression_ptr(mul->operands()[i]);
            std::shared_ptr<SymbolicExpr> at, rt;
            if (te_extract_exp(fac, t, at, rt)) {
                std::vector<std::shared_ptr<const SymbolicNode>> rest = mul->operands();
                rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(i));
                rem_out = detail::make_expression_ptr(
                    SymbolicFactory::create_multiply(std::move(rest)));
                a_out = at; return true;
            }
        }
    }
    return false;
}

static std::shared_ptr<SymbolicExpr> te_laplace_lookup(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& s, ComputationContext& context) {
    if (!f || !LMCAS::detail::node(f)) return nullptr;
    auto sv = SymbolicExpr::variable(s);
    if (!te_depends_on(f, t)) return SymbolicExpr::divide(f, sv);
    BigInt n;
    if (te_is_power_of_var(f, t, n)) {
        return SymbolicExpr::divide(SymbolicExpr::number(te_factorial(n, context)),
            SymbolicExpr::power(sv, SymbolicExpr::number(n + BigInt(1))));
    }
    bool is_sin = false; std::shared_ptr<SymbolicExpr> freq;
    if (te_is_trig(f, t, is_sin, freq)) {
        auto den = SymbolicExpr::add(SymbolicExpr::power(sv, SymbolicExpr::number(2)), SymbolicExpr::power(freq, SymbolicExpr::number(2)));
        return is_sin ? SymbolicExpr::divide(freq, den) : SymbolicExpr::divide(sv, den);
    }
    bool is_sinh = false; std::shared_ptr<SymbolicExpr> hf;
    if (te_is_hyp(f, t, is_sinh, hf)) {
        auto den = SymbolicExpr::add(SymbolicExpr::power(sv, SymbolicExpr::number(2)), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(hf, SymbolicExpr::number(2))));
        return is_sinh ? SymbolicExpr::divide(hf, den) : SymbolicExpr::divide(sv, den);
    }
    std::shared_ptr<SymbolicExpr> a, rem;
    if (te_extract_exp(f, t, a, rem) && !te_depends_on(rem, t))
        return SymbolicExpr::divide(rem, SymbolicExpr::add(sv, SymbolicExpr::multiply(SymbolicExpr::number(-1), a)));
    return nullptr;
}

static std::vector<std::shared_ptr<SymbolicExpr>> te_laplace_roc(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s) {
    if (!f || !LMCAS::detail::node(f)) return {};
    auto sv = SymbolicExpr::variable(s);
    if (!te_depends_on(f, t)) {
        auto condition = te_gt_condition(sv, SymbolicExpr::number(0));
        return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                         : std::vector<std::shared_ptr<SymbolicExpr>>{};
    }
    BigInt n;
    if (te_is_power_of_var(f, t, n)) {
        auto condition = te_gt_condition(sv, SymbolicExpr::number(0));
        return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                         : std::vector<std::shared_ptr<SymbolicExpr>>{};
    }
    bool is_sin = false;
    std::shared_ptr<SymbolicExpr> freq;
    if (te_is_trig(f, t, is_sin, freq)) {
        auto condition = te_gt_condition(sv, SymbolicExpr::number(0));
        return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                         : std::vector<std::shared_ptr<SymbolicExpr>>{};
    }
    std::shared_ptr<SymbolicExpr> a;
    std::shared_ptr<SymbolicExpr> rem;
    if (te_extract_exp(f, t, a, rem) && !te_depends_on(rem, t)) {
        auto condition = te_gt_condition(sv, a);
        return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                         : std::vector<std::shared_ptr<SymbolicExpr>>{};
    }
    return {};
}

static std::shared_ptr<SymbolicExpr> te_laplace_sum(
    const std::shared_ptr<SymbolicExpr>& f, const AddNode& addition,
    const std::string& t, const std::string& s, ComputationContext& context) {
    std::shared_ptr<SymbolicExpr> result;
    for (const auto& operand : addition.operands()) {
        auto transformed = laplace_transform_core(
            detail::make_expression_ptr(operand), t, s, context);
        if (!transformed) return te_unevaluated_laplace(f, t, s);
        result = result ? SymbolicExpr::add(result, transformed) : transformed;
    }
    return result;
}

std::shared_ptr<SymbolicExpr> transform_detail::laplace_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& s, ComputationContext& context) {
    auto step = context.consume_steps(1, "laplace_transform.recursive");
    if (!step) {
        return nullptr;
    }
    if (!f || !LMCAS::detail::node(f)) {
        return nullptr;
    }
    if (te_is_zero(f)) {
        return SymbolicExpr::number(0);
    }
    if (!te_depends_on(f, t)) {
        return SymbolicExpr::divide(f, SymbolicExpr::variable(s));
    }
    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f));
    if (add) {
        return te_laplace_sum(f, *add, t, s, context);
    }
    auto [coeff, body] = te_split_coeff(f, t);
    if (!coeff->is_one()) {
        auto lt_body = laplace_transform_core(body, t, s, context);
        if (lt_body) {
            return SymbolicExpr::multiply(coeff, lt_body);
        }
    }
    std::shared_ptr<SymbolicExpr> a_shift, remainder;
    if (te_extract_exp(f, t, a_shift, remainder)) {
        if (te_depends_on(remainder, t)) {
            auto lt_rem = laplace_transform_core(remainder, t, s, context);
            if (lt_rem) {
                return lt_rem->substitute(s, SymbolicExpr::add(SymbolicExpr::variable(s), SymbolicExpr::multiply(SymbolicExpr::number(-1), a_shift)));
            }
        } else {
            return SymbolicExpr::divide(remainder, SymbolicExpr::add(SymbolicExpr::variable(s), SymbolicExpr::multiply(SymbolicExpr::number(-1), a_shift)));
        }
    }
    auto lookup = te_laplace_lookup(f, t, s, context);
    if (lookup) {
        return lookup;
    }
    return te_unevaluated_laplace(f, t, s);
}

TransformEngineResult laplace_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s,
    ComputationContext& context) {
    const std::string operation = "laplace_transform";
    auto valid = te_validate_expr_vars(f, t, s, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = laplace_transform_core(f, t, s, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        auto result = te_wrap_transform_result(std::move(expression), operation);
        if (result) {
            result.value().value.roc = te_laplace_roc(f, t, s);
        }
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

TransformEngineResult laplace_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& s) {
    ComputationContext context;
    return laplace_transform_checked(f, t, s, context);
}

}
