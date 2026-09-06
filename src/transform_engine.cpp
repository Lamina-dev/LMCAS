/**
 * @file transform_engine.cpp
 * @brief 积分变换引擎实现：Laplace 变换、逆 Laplace 变换、Fourier 变换、Z 变换。
 *
 * 算法来源:
 * - Laplace/Fourier: Schiff, The Laplace Transform: Theory and Applications, Springer.
 * - Z 变换: Oppenheim & Willsky, Signals and Systems, 2nd ed., Chapter 10.
 */
#include "transform_engine.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "residual_verification.hpp"
#include <cmath>
#include <algorithm>
#include <exception>
#include <variant>
#include <map>
#include <functional>

namespace LMCAS {



static bool te_depends_on(const std::shared_ptr<SymbolicExpr>& expression,
                          const std::string& variable) {
    return expression &&
           expression_depends_on_variable(
               LMCAS::detail::node(expression), variable);
}

static Result<void> te_validate_expr_vars(const std::shared_ptr<SymbolicExpr>& expr,
                                          const std::string& input_var,
                                          const std::string& output_var,
                                          ComputationContext& context,
                                          const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "transform expression cannot be null",
                                     operation);
    }
    if (input_var.empty() || output_var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "transform variable names cannot be empty",
                                     operation);
    }
    if (input_var == output_var) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "input and output variables must be distinct",
                                     operation);
    }
    return Result<void>::success();
}

static Result<void> te_validate_convolution_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!f || !LMCAS::detail::node(f) || !g || !LMCAS::detail::node(g)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "convolution expressions cannot be null",
                                     operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "convolution variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

static bool te_contains_transform(
    const std::shared_ptr<const SymbolicNode>& node) {
    return LMCAS::detail::contains_node_type<TransformNode>(node);
}

static TransformEngineResult te_wrap_transform_result(
    std::shared_ptr<SymbolicExpr> expr,
    const std::string& operation) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return TransformEngineResult::failure(
            CasErrc::InternalInvariant,
            "transform construction produced a null expression",
            operation);
    }

    if (te_contains_transform(LMCAS::detail::node(expr))) {
        return TransformEngineResult::failure(
            CasErrc::Inconclusive,
            "transform is outside the current evaluated support domain",
            operation);
    }

    return TransformEngineResult::success(
        Verified<EvaluatedTransform>{
            EvaluatedTransform{std::move(expr), {}, {}},
            ByConstructionProof{}});
}

static std::shared_ptr<SymbolicExpr> te_gt_condition(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs) {
    if (!lhs || !LMCAS::detail::node(lhs) || !rhs || !LMCAS::detail::node(rhs)) return nullptr;
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<RelationalNode>(LMCAS::detail::node(lhs), LMCAS::detail::node(rhs), RelationalNode::Op::GT));
}

static bool te_is_zero(const std::shared_ptr<SymbolicExpr>& e) {
    return e && LMCAS::detail::node(e) && e->is_zero();
}
static double te_numval(const std::shared_ptr<const SymbolicNode>& nd) {
    auto n = std::dynamic_pointer_cast<const NumberNode>(nd);
    if (!n) return 0.0;
    if (std::holds_alternative<BigInt>(n->value())) return std::get<BigInt>(n->value()).to_double();
    if (std::holds_alternative<Rational>(n->value())) return std::get<Rational>(n->value()).to_double();
    return static_cast<double>(std::get<lmmc_real_t>(n->value()));
}
static long long te_factorial(int n) {
    long long r = 1; for (int i = 2; i <= n; ++i) r *= i; return r;
}
static std::shared_ptr<SymbolicExpr> te_unevaluated_laplace(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t, const std::string& s) {
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::Laplace, LMCAS::detail::node(f)->clone(), t,
        SymbolicFactory::create_variable(s)));
}
static std::shared_ptr<SymbolicExpr> te_unevaluated_inv_laplace(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s, const std::string& t) {
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::InverseLaplace, LMCAS::detail::node(F)->clone(), s,
        SymbolicFactory::create_variable(t)));
}
static std::pair<std::shared_ptr<SymbolicExpr>,std::shared_ptr<SymbolicExpr>>
te_split_coeff(const std::shared_ptr<SymbolicExpr>& e, const std::string& v) {
    if (!e || !LMCAS::detail::node(e) || !te_depends_on(e, v)) return {e ? e : SymbolicExpr::number(1), SymbolicExpr::number(1)};
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(e));
    if (!mul) return {SymbolicExpr::number(1), e};
    std::vector<std::shared_ptr<const SymbolicNode>> cp, vp;
    for (auto& op : mul->operands()) {
        if (expression_depends_on_variable(op, v)) vp.push_back(op); else cp.push_back(op);
    }
    if (cp.empty()) return {SymbolicExpr::number(1), e};
    if (vp.empty()) return {e, SymbolicExpr::number(1)};
    auto c = (cp.size()==1) ? LMCAS::detail::make_expression_ptr(cp[0])
        : LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(cp));
    auto b = (vp.size()==1) ? LMCAS::detail::make_expression_ptr(vp[0])
        : LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(vp));
    return {c, b};
}

static bool te_is_power_of_var(const std::shared_ptr<SymbolicExpr>& e, const std::string& v, int& n) {
    if (!e || !LMCAS::detail::node(e)) return false;
    auto vn = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(e));
    if (vn && vn->name() == v) { n = 1; return true; }
    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(e));
    if (!pw) return false;
    auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (!bv || bv->name() != v) return false;
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) return false;
    double val = te_numval(pw->exponent());
    int iv = static_cast<int>(val);
    if (std::abs(val - iv) < 1e-12 && iv >= 0) { n = iv; return true; }
    return false;
}
static bool te_is_trig(const std::shared_ptr<SymbolicExpr>& e, const std::string& v,
                       bool& is_sin, std::shared_ptr<SymbolicExpr>& freq) {
    if (!e || !LMCAS::detail::node(e)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(e));
    if (!fn || fn->arguments().empty()) return false;
    if (fn->type() != FunctionNode::FuncType::Sin && fn->type() != FunctionNode::FuncType::Cos) return false;
    is_sin = (fn->type() == FunctionNode::FuncType::Sin);
    auto arg = fn->arguments()[0];
    auto ml = std::dynamic_pointer_cast<const MultiplyNode>(arg);
    if (ml && ml->operands().size() == 2) {
        for (size_t i = 0; i < 2; ++i) {
            auto vn = std::dynamic_pointer_cast<const VariableNode>(ml->operands()[i]);
            if (vn && vn->name() == v && !expression_depends_on_variable(ml->operands()[1-i], v)) {
                freq = LMCAS::detail::make_expression_ptr(ml->operands()[1-i]); return true;
            }
        }
    }
    auto vn = std::dynamic_pointer_cast<const VariableNode>(arg);
    if (vn && vn->name() == v) { freq = SymbolicExpr::number(1); return true; }
    return false;
}
static bool te_is_hyp(const std::shared_ptr<SymbolicExpr>& e, const std::string& v,
                      bool& is_sinh, std::shared_ptr<SymbolicExpr>& freq) {
    if (!e || !LMCAS::detail::node(e)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(e));
    if (!fn || fn->arguments().empty()) return false;
    if (fn->type() != FunctionNode::FuncType::Sinh && fn->type() != FunctionNode::FuncType::Cosh) return false;
    is_sinh = (fn->type() == FunctionNode::FuncType::Sinh);
    auto arg = fn->arguments()[0];
    auto ml = std::dynamic_pointer_cast<const MultiplyNode>(arg);
    if (ml && ml->operands().size() == 2) {
        for (size_t i = 0; i < 2; ++i) {
            auto vn = std::dynamic_pointer_cast<const VariableNode>(ml->operands()[i]);
            if (vn && vn->name() == v && !expression_depends_on_variable(ml->operands()[1-i], v)) {
                freq = LMCAS::detail::make_expression_ptr(ml->operands()[1-i]); return true;
            }
        }
    }
    auto vn = std::dynamic_pointer_cast<const VariableNode>(arg);
    if (vn && vn->name() == v) { freq = SymbolicExpr::number(1); return true; }
    return false;
}

static bool te_extract_exp(const std::shared_ptr<SymbolicExpr>& e, const std::string& t,
                           std::shared_ptr<SymbolicExpr>& a_out, std::shared_ptr<SymbolicExpr>& rem_out) {
    if (!e || !LMCAS::detail::node(e)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(e));
    if (fn && fn->type() == FunctionNode::FuncType::Exp && !fn->arguments().empty()) {
        auto arg = fn->arguments()[0];
        auto ml = std::dynamic_pointer_cast<const MultiplyNode>(arg);
        if (ml && ml->operands().size() == 2) {
            for (size_t i = 0; i < 2; ++i) {
                auto vn = std::dynamic_pointer_cast<const VariableNode>(ml->operands()[i]);
                if (vn && vn->name() == t && !expression_depends_on_variable(ml->operands()[1-i], t)) {
                    a_out = LMCAS::detail::make_expression_ptr(ml->operands()[1-i]);
                    rem_out = SymbolicExpr::number(1); return true;
                }
            }
        }
        auto vn = std::dynamic_pointer_cast<const VariableNode>(arg);
        if (vn && vn->name() == t) { a_out = SymbolicExpr::number(1); rem_out = SymbolicExpr::number(1); return true; }
    }
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(e));
    if (mul) {
        for (size_t i = 0; i < mul->operands().size(); ++i) {
            auto fac = LMCAS::detail::make_expression_ptr(mul->operands()[i]);
            std::shared_ptr<SymbolicExpr> at, rt;
            if (te_extract_exp(fac, t, at, rt)) {
                std::vector<std::shared_ptr<const SymbolicNode>> rest;
                for (size_t j = 0; j < mul->operands().size(); ++j)
                    if (j != i) rest.push_back(mul->operands()[j]);
                if (rest.empty()) rem_out = SymbolicExpr::number(1);
                else if (rest.size() == 1) rem_out = LMCAS::detail::make_expression_ptr(rest[0]);
                else rem_out = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(std::move(rest)));
                a_out = at; return true;
            }
        }
    }
    return false;
}

TransformTable::TransformTable() { init_laplace_pairs(); }
void TransformTable::add_entry(TransformTableEntry entry) { entries_.push_back(std::move(entry)); }
void TransformTable::init_laplace_pairs() { /* 变换对通过模式匹配硬编码 */ }

static std::shared_ptr<SymbolicExpr> te_laplace_lookup(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t, const std::string& s) {
    if (!f || !LMCAS::detail::node(f)) return nullptr;
    auto sv = SymbolicExpr::variable(s);
    if (!te_depends_on(f, t)) return SymbolicExpr::divide(f, sv);
    int n = 0;
    if (te_is_power_of_var(f, t, n)) {
        return SymbolicExpr::divide(SymbolicExpr::number(static_cast<long long>(te_factorial(n))),
            SymbolicExpr::power(sv, SymbolicExpr::number(n + 1)));
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
    int n = 0;
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

static std::shared_ptr<SymbolicExpr> laplace_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& s, ComputationContext& context) {
    auto step = context.consume_steps(1, "laplace_transform.recursive");
    if (!step) return nullptr;
    if (!f || !LMCAS::detail::node(f)) return nullptr;
    if (te_is_zero(f)) return SymbolicExpr::number(0);
    if (!te_depends_on(f, t)) return SymbolicExpr::divide(f, SymbolicExpr::variable(s));
    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f));
    if (add) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add->operands()) {
            auto lt = laplace_transform_core(
                LMCAS::detail::make_expression_ptr(op), t, s, context);
            if (!lt) return te_unevaluated_laplace(f, t, s);
            result = result ? SymbolicExpr::add(result, lt) : lt;
        }
        return result;
    }
    auto [coeff, body] = te_split_coeff(f, t);
    if (!coeff->is_one()) {
        auto lt_body = laplace_transform_core(body, t, s, context);
        if (lt_body) return SymbolicExpr::multiply(coeff, lt_body);
    }
    std::shared_ptr<SymbolicExpr> a_shift, remainder;
    if (te_extract_exp(f, t, a_shift, remainder)) {
        if (te_depends_on(remainder, t)) {
            auto lt_rem = laplace_transform_core(remainder, t, s, context);
            if (lt_rem) return lt_rem->substitute(s, SymbolicExpr::add(SymbolicExpr::variable(s), SymbolicExpr::multiply(SymbolicExpr::number(-1), a_shift)));
        } else {
            return SymbolicExpr::divide(remainder, SymbolicExpr::add(SymbolicExpr::variable(s), SymbolicExpr::multiply(SymbolicExpr::number(-1), a_shift)));
        }
    }
    auto lookup = te_laplace_lookup(f, t, s);
    if (lookup) return lookup;
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


static std::shared_ptr<SymbolicExpr> te_inv_power(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s, const std::string& t) {
    if (!F || !LMCAS::detail::node(F)) return nullptr;
    auto tv = SymbolicExpr::variable(t);
    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(F));
    if (!pw) return nullptr;
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) return nullptr;
    double ev = te_numval(pw->exponent());
    if (ev >= 0) return nullptr;
    int n = static_cast<int>(-ev);
    if (std::abs(ev + n) > 1e-12 || n < 1) return nullptr;
    auto make_exp = [&](const std::shared_ptr<SymbolicExpr>& a) {
        auto arg = SymbolicExpr::multiply(a, tv);
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(arg)}));
    };
    auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
    if (bv && bv->name() == s) {
        if (n == 1) return SymbolicExpr::number(1);
        return SymbolicExpr::divide(SymbolicExpr::power(tv, SymbolicExpr::number(n-1)),
            SymbolicExpr::number(static_cast<long long>(te_factorial(n-1))));
    }
    auto ba = std::dynamic_pointer_cast<const AddNode>(pw->base());
    if (ba && ba->operands().size() == 2) {
        std::shared_ptr<SymbolicExpr> pole;
        for (size_t i = 0; i < 2; ++i) {
            auto vn = std::dynamic_pointer_cast<const VariableNode>(ba->operands()[i]);
            if (vn && vn->name() == s) {
                auto c = LMCAS::detail::make_expression_ptr(ba->operands()[1-i]);
                pole = SymbolicExpr::multiply(SymbolicExpr::number(-1), c); break;
            }
        }
        if (pole) {
            if (n == 1) {
                if (te_is_zero(pole)) return SymbolicExpr::number(1);
                return make_exp(pole);
            }
            auto tp = SymbolicExpr::power(tv, SymbolicExpr::number(n-1));
            auto fv = SymbolicExpr::number(static_cast<long long>(te_factorial(n-1)));
            if (te_is_zero(pole)) return SymbolicExpr::divide(tp, fv);
            return SymbolicExpr::divide(SymbolicExpr::multiply(tp, make_exp(pole)), fv);
        }
    }
    return nullptr;
}

static std::shared_ptr<SymbolicExpr> te_inv_product(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s, const std::string& t) {
    if (!F || !LMCAS::detail::node(F)) return nullptr;
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(F));
    if (!mul || mul->operands().size() < 2) return nullptr;
    std::vector<std::shared_ptr<const SymbolicNode>> num_parts;
    std::shared_ptr<const SymbolicNode> den_part = nullptr;
    for (auto& op : mul->operands()) {
        auto pw = std::dynamic_pointer_cast<const PowerNode>(op);
        if (pw && !den_part) {
            double ev = te_numval(pw->exponent());
            if (ev < 0 && expression_depends_on_variable(pw->base(), s)) { den_part = op; continue; }
        }
        num_parts.push_back(op);
    }
    if (!den_part) return nullptr;
    std::shared_ptr<SymbolicExpr> num_expr;
    if (num_parts.empty()) num_expr = SymbolicExpr::number(1);
    else if (num_parts.size() == 1) num_expr = LMCAS::detail::make_expression_ptr(num_parts[0]);
    else num_expr = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(num_parts));
    if (!te_depends_on(num_expr, s)) {
        auto inv = te_inv_power(LMCAS::detail::make_expression_ptr(den_part), s, t);
        if (inv) return SymbolicExpr::multiply(num_expr, inv);
    }
    return nullptr;
}

static std::shared_ptr<SymbolicExpr> inverse_laplace_core(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& s,
    const std::string& t, ComputationContext& context) {
    auto step = context.consume_steps(1, "inverse_laplace.recursive");
    if (!step) return nullptr;
    if (!F || !LMCAS::detail::node(F)) return nullptr;
    if (te_is_zero(F)) return SymbolicExpr::number(0);
    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(F));
    if (add) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add->operands()) {
            auto ilt = inverse_laplace_core(
                LMCAS::detail::make_expression_ptr(op), s, t, context);
            if (!ilt) return te_unevaluated_inv_laplace(F, s, t);
            result = result ? SymbolicExpr::add(result, ilt) : ilt;
        }
        return result;
    }
    auto [coeff, body] = te_split_coeff(F, s);
    if (!coeff->is_one()) {
        auto ilt_body = inverse_laplace_core(body, s, t, context);
        if (ilt_body) return SymbolicExpr::multiply(coeff, ilt_body);
    }
    auto pw_res = te_inv_power(F, s, t);
    if (pw_res) return pw_res;
    auto prod_res = te_inv_product(F, s, t);
    if (prod_res) return prod_res;
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



/**
 * @internal
 * @brief 检测 e^(-a*t^2) 高斯形式并提取系数 a > 0。
 */
static bool te_is_gaussian(const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
                            double& a_out) {
    if (!f || !LMCAS::detail::node(f)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(f));
    if (!fn || fn->type() != FunctionNode::FuncType::Exp || fn->arguments().empty()) return false;
    auto arg = fn->arguments()[0];
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(arg);
    if (!mul) return false;
    double coeff_val = 1.0;
    bool has_t_sq = false;
    for (auto& op : mul->operands()) {
        auto nd = std::dynamic_pointer_cast<const NumberNode>(op);
        if (nd) { coeff_val *= te_numval(op); continue; }
        auto pw = std::dynamic_pointer_cast<const PowerNode>(op);
        if (pw) {
            auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
            double ev = te_numval(pw->exponent());
            if (bv && bv->name() == t && std::abs(ev - 2.0) < 1e-12) {
                has_t_sq = true; continue;
            }
        }
        return false; // 含有其他因子
    }
    if (has_t_sq && coeff_val < 0) {
        a_out = -coeff_val;
        return true;
    }
    return false;
}

/**
 * @internal
 * @brief 检测 e^(-a*|t|) 形式并提取系数 a（a 可为符号或数值正数）。
 */
static bool te_is_abs_exp(const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
                           std::shared_ptr<SymbolicExpr>& a_out) {
    if (!f || !LMCAS::detail::node(f)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(f));
    if (!fn || fn->type() != FunctionNode::FuncType::Exp || fn->arguments().empty()) return false;
    auto arg = fn->arguments()[0];
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(arg);
    if (!mul) return false;
    std::shared_ptr<SymbolicExpr> a_candidate;
    bool has_abs_t = false;
    bool has_neg = false;
    for (auto& op : mul->operands()) {
        auto fn2 = std::dynamic_pointer_cast<const FunctionNode>(op);
        if (fn2 && fn2->type() == FunctionNode::FuncType::Abs && fn2->arguments().size() == 1) {
            auto av = std::dynamic_pointer_cast<const VariableNode>(fn2->arguments()[0]);
            if (av && av->name() == t) { has_abs_t = true; continue; }
        }
        auto nd = std::dynamic_pointer_cast<const NumberNode>(op);
        if (nd) {
            double v = te_numval(op);
            if (v < 0) {
                has_neg = true;
                a_candidate = SymbolicExpr::number(-v);
            } else {
                a_candidate = LMCAS::detail::make_expression_ptr(op);
            }
            continue;
        }
        /// 允许符号变量作为系数（如 -a）
        auto vn = std::dynamic_pointer_cast<const VariableNode>(op);
        if (vn) {
            /// 假设符号变量为正 a，前面有负号处理
            a_candidate = LMCAS::detail::make_expression_ptr(op);
            continue;
        }
        return false;
    }
    if (has_abs_t && has_neg && a_candidate) {
        a_out = a_candidate;
        return true;
    }
    return false;
}

static std::shared_ptr<SymbolicExpr> fourier_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& omega, ComputationContext& context) {
    auto step = context.consume_steps(1, "fourier_transform.recursive");
    if (!step) return nullptr;
    if (!f || !LMCAS::detail::node(f)) return nullptr;
    auto wv = SymbolicExpr::variable(omega);

    /// 不依赖 t：无有限 Fourier 变换（δ函数），返回未求值
    if (!te_depends_on(f, t)) {
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Fourier, LMCAS::detail::node(f)->clone(), t,
            SymbolicFactory::create_variable(omega)));
    }

    /// 线性性：处理加法 F{af+bg} = aF{f} + bF{g}
    auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f));
    if (add_node) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add_node->operands()) {
            auto ft = fourier_transform_core(
                LMCAS::detail::make_expression_ptr(op), t, omega, context);
            if (!ft) {
                return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
                    TransformNode::TransformType::Fourier, LMCAS::detail::node(f)->clone(), t,
                    SymbolicFactory::create_variable(omega)));
            }
            result = result ? SymbolicExpr::add(result, ft) : ft;
        }
        return result;
    }

    /// 分离常数系数: c * body
    auto [coeff, body] = te_split_coeff(f, t);
    if (!coeff->is_one()) {
        auto ft_body = fourier_transform_core(body, t, omega, context);
        if (ft_body) return SymbolicExpr::multiply(coeff, ft_body);
    }

    /// F{e^(-a*t^2)} = sqrt(π/a) * e^(-ω²/(4a))
    {
        double a_val = 0.0;
        if (te_is_gaussian(f, t, a_val)) {
            double sqrt_pi_over_a = std::sqrt(3.14159265358979323846 / a_val);
            auto coeff_expr = SymbolicExpr::number(sqrt_pi_over_a);
            auto neg_w_sq_over_4a = SymbolicExpr::multiply(
                SymbolicExpr::number(-1.0 / (4.0 * a_val)),
                SymbolicExpr::power(wv, SymbolicExpr::number(2)));
            auto exp_part = SymbolicExpr::exp(neg_w_sq_over_4a);
            return SymbolicExpr::multiply(coeff_expr, exp_part);
        }
    }

    /// F{e^(-a|t|)} = 2a/(a² + ω²)
    {
        std::shared_ptr<SymbolicExpr> a_expr;
        if (te_is_abs_exp(f, t, a_expr)) {
            auto a_sq = SymbolicExpr::power(a_expr, SymbolicExpr::number(2));
            auto w_sq = SymbolicExpr::power(wv, SymbolicExpr::number(2));
            auto den = SymbolicExpr::add(a_sq, w_sq);
            auto num_expr = SymbolicExpr::multiply(SymbolicExpr::number(2), a_expr);
            return SymbolicExpr::divide(num_expr, den);
        }
    }

    /// F{e^(-a*t)} (因果指数，a > 0): 1/(a + iω)
    {
        /// 提取指数 e^(arg)，其中 arg 关于 t 线性：arg = c * t（c 不依赖 t）。
        auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(f));
        if (fn && fn->type() == FunctionNode::FuncType::Exp && fn->arguments().size() == 1) {
            auto arg = LMCAS::detail::make_expression_ptr(fn->arguments()[0]);
            /// c = d(arg)/dt；若 arg 关于 t 线性则 c 不依赖 t。
            auto c = arg->differentiate(t);
            if (c && !te_depends_on(c, t)) {
                /// arg - c*t 应为 0（纯线性，无常数项）才是标准因果指数。
                auto linear = SymbolicExpr::multiply(c, SymbolicExpr::variable(t));
                auto residual = SymbolicExpr::add(arg,
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), linear))->simplify();
                if (LMCAS::detail::node(residual) && LMCAS::detail::node(residual)->is_zero()) {
                    /// 衰减率 a = -c（要求 a>0，即 c 为负）。ℱ = 1/(a + iω)。
                    auto a = SymbolicExpr::multiply(SymbolicExpr::number(-1), c)->simplify();
                    auto i_unit = LMCAS::detail::make_expression_ptr(
                        SymbolicFactory::create_complex(
                            LMCAS::detail::node(SymbolicExpr::number(0)),
                            LMCAS::detail::node(SymbolicExpr::number(1))));
                    auto iw = SymbolicExpr::multiply(i_unit, wv);
                    auto denom = SymbolicExpr::add(a, iw);
                    return SymbolicExpr::divide(SymbolicExpr::number(1), denom)->simplify();
                }
            }
        }
    }

    /// 未知形式：返回未求值 TransformNode
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::Fourier, LMCAS::detail::node(f)->clone(), t,
        SymbolicFactory::create_variable(omega)));
}

TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega,
    ComputationContext& context) {
    const std::string operation = "fourier_transform";
    auto valid = te_validate_expr_vars(f, t, omega, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = fourier_transform_core(f, t, omega, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        return te_wrap_transform_result(std::move(expression), operation);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega) {
    ComputationContext context;
    return fourier_transform_checked(f, t, omega, context);
}



static std::shared_ptr<SymbolicExpr> inverse_fourier_transform_core(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& omega,
    const std::string& t, ComputationContext& context) {
    auto step = context.consume_steps(1, "inverse_fourier_transform.recursive");
    if (!step) return nullptr;
    if (!F || !LMCAS::detail::node(F)) return nullptr;
    auto tv = SymbolicExpr::variable(t);

    /// 不依赖 ω：返回未求值
    if (!te_depends_on(F, omega)) {
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
            SymbolicFactory::create_variable(t)));
    }

    /// 线性性
    auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(F));
    if (add_node) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add_node->operands()) {
            auto ift = inverse_fourier_transform_core(
                LMCAS::detail::make_expression_ptr(op), omega, t, context);
            if (!ift) {
                return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
                    TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
                    SymbolicFactory::create_variable(t)));
            }
            result = result ? SymbolicExpr::add(result, ift) : ift;
        }
        return result;
    }

    /// F⁻¹{e^(-a*ω²)} = (1/sqrt(4πa)) * e^(-t²/(4a))
    {
        double a_val = 0.0;
        if (te_is_gaussian(F, omega, a_val)) {
            double coeff_val = 1.0 / std::sqrt(4.0 * 3.14159265358979323846 * a_val);
            auto coeff_expr = SymbolicExpr::number(coeff_val);
            auto neg_t_sq_over_4a = SymbolicExpr::multiply(
                SymbolicExpr::number(-1.0 / (4.0 * a_val)),
                SymbolicExpr::power(tv, SymbolicExpr::number(2)));
            auto exp_part = SymbolicExpr::exp(neg_t_sq_over_4a);
            return SymbolicExpr::multiply(coeff_expr, exp_part);
        }
    }

    /// 未知形式：返回未求值
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
        SymbolicFactory::create_variable(t)));
}

TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t,
    ComputationContext& context) {
    const std::string operation = "inverse_fourier_transform";
    auto valid = te_validate_expr_vars(F, omega, t, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = inverse_fourier_transform_core(
            F, omega, t, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        auto result = te_wrap_transform_result(expression, operation);
        if (!result) return result;
        auto round_trip = fourier_transform_core(
            expression, t, omega, context);
        if (!round_trip ||
            te_contains_transform(LMCAS::detail::node(round_trip))) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Fourier round trip is not proved", operation);
        }
        auto verified = check_equivalent(round_trip, F, context);
        if (!verified) return TransformEngineResult::failure(verified.error());
        if (!std::holds_alternative<ProvedZeroResidual>(verified.value())) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Fourier round trip is not proved", operation);
        }
        result.value().certificate =
            ExactRoundTripProof{SymbolicExpr::number(0)};
        return result;
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t) {
    ComputationContext context;
    return inverse_fourier_transform_checked(F, omega, t, context);
}



struct TeGaussian {
    std::shared_ptr<SymbolicExpr> coefficient;
    double rate = 0.0;
};

static bool te_numeric_expr(
    const std::shared_ptr<SymbolicExpr>& expression,
    double& value)
{
    if (!expression || !LMCAS::detail::node(expression)) return false;
    auto number = std::dynamic_pointer_cast<const NumberNode>(
        LMCAS::detail::node(expression));
    if (!number) return false;
    value = te_numval(LMCAS::detail::node(expression));
    return std::isfinite(value);
}

static bool te_var_square(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable)
{
    auto power = expression
        ? std::dynamic_pointer_cast<const PowerNode>(
              LMCAS::detail::node(expression))
        : nullptr;
    if (!power) return false;
    auto base = std::dynamic_pointer_cast<const VariableNode>(power->base());
    auto exponent =
        std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    return base && exponent && base->name() == variable &&
           std::fabs(te_numval(power->exponent()) - 2.0) < 1e-12;
}

static bool te_match_positive_gaussian(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    TeGaussian& gaussian)
{
    auto parts = te_split_coeff(expression, variable);
    double coefficient = 0.0;
    if (!te_numeric_expr(parts.first, coefficient)) return false;

    auto function = std::dynamic_pointer_cast<const FunctionNode>(
        LMCAS::detail::node(parts.second));
    if (!function || function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().size() != 1) {
        return false;
    }
    auto argument =
        LMCAS::detail::make_expression_ptr(function->arguments()[0]);
    auto exponent_parts = te_split_coeff(argument, variable);
    double exponent_coefficient = 0.0;
    if (!te_numeric_expr(exponent_parts.first, exponent_coefficient) ||
        !te_var_square(exponent_parts.second, variable) ||
        !(exponent_coefficient < 0.0)) {
        return false;
    }
    gaussian = {parts.first, -exponent_coefficient};
    return std::isfinite(gaussian.rate) && gaussian.rate > 0.0;
}

static TransformEngineResult te_convolve_core(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& variable,
    ComputationContext& context)
{
    constexpr const char* operation = "convolve";
    auto step = context.consume_steps(1, operation);
    if (!step) return TransformEngineResult::failure(step.error());

    if (te_is_zero(f) || te_is_zero(g)) {
        return te_wrap_transform_result(SymbolicExpr::number(0), operation);
    }

    auto convolve_sum = [&](const std::shared_ptr<SymbolicExpr>& sum,
                            const std::shared_ptr<SymbolicExpr>& other,
                            bool sum_is_left) -> TransformEngineResult {
        auto addition = std::dynamic_pointer_cast<const AddNode>(
            LMCAS::detail::node(sum));
        if (!addition) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "当前双边卷积规则不支持该表达式对", operation);
        }
        std::shared_ptr<SymbolicExpr> accumulated;
        for (const auto& operand : addition->operands()) {
            auto term = LMCAS::detail::make_expression_ptr(operand);
            auto transformed = sum_is_left
                ? te_convolve_core(term, other, variable, context)
                : te_convolve_core(other, term, variable, context);
            if (!transformed) return transformed;
            auto expression = transformed.value().value.expression;
            accumulated = accumulated
                ? SymbolicExpr::add(accumulated, expression)
                : expression;
        }
        return te_wrap_transform_result(
            accumulated ? accumulated->simplify() : SymbolicExpr::number(0),
            operation);
    };

    if (std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f))) {
        return convolve_sum(f, g, true);
    }
    if (std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(g))) {
        return convolve_sum(g, f, false);
    }

    auto left_parts = te_split_coeff(f, variable);
    if (te_depends_on(left_parts.second, variable) &&
        left_parts.second->to_string() != f->to_string()) {
        auto inner =
            te_convolve_core(left_parts.second, g, variable, context);
        if (!inner) return inner;
        auto scaled = SymbolicExpr::multiply(
            left_parts.first, inner.value().value.expression)->simplify();
        return te_wrap_transform_result(scaled, operation);
    }
    auto right_parts = te_split_coeff(g, variable);
    if (te_depends_on(right_parts.second, variable) &&
        right_parts.second->to_string() != g->to_string()) {
        auto inner =
            te_convolve_core(f, right_parts.second, variable, context);
        if (!inner) return inner;
        auto scaled = SymbolicExpr::multiply(
            right_parts.first, inner.value().value.expression)->simplify();
        return te_wrap_transform_result(scaled, operation);
    }

    TeGaussian left;
    TeGaussian right;
    if (!te_match_positive_gaussian(f, variable, left) ||
        !te_match_positive_gaussian(g, variable, right)) {
        return TransformEngineResult::failure(
            CasErrc::Inconclusive,
            "当前双边卷积规则不支持该表达式对", operation);
    }

    const double rate_sum = left.rate + right.rate;
    const double output_rate = left.rate * right.rate / rate_sum;
    auto scale = SymbolicExpr::multiply(
        SymbolicExpr::multiply(left.coefficient, right.coefficient),
        SymbolicExpr::sqrt(SymbolicExpr::divide(
            SymbolicExpr::variable("pi"),
            SymbolicExpr::number(rate_sum))));
    auto square = SymbolicExpr::power(
        SymbolicExpr::variable(variable), SymbolicExpr::number(2));
    auto exponent = SymbolicExpr::multiply(
        SymbolicExpr::number(-output_rate), square);
    auto gaussian = SymbolicExpr::exp(exponent);
    return te_wrap_transform_result(
        SymbolicExpr::multiply(scale, gaussian)->simplify(), operation);
}

TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var,
    ComputationContext& context) {
    const std::string operation = "convolve";
    auto valid =
        te_validate_convolution_inputs(f, g, var, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    try {
        return te_convolve_core(f, g, var, context);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(
            CasErrc::ResourceLimit, "双边卷积分配失败", operation);
    } catch (const std::exception& error) {
        return TransformEngineResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var) {
    ComputationContext context;
    return convolve_checked(f, g, var, context);
}


/**
 * @internal
 * @brief 检测 a^n 形式并提取底数 a（a 不依赖 n）。
 */
static bool zt_is_exp_seq(const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
                           std::shared_ptr<SymbolicExpr>& base_out) {
    if (!f || !LMCAS::detail::node(f)) return false;
    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(f));
    if (!pw) return false;
    auto exp_node = std::dynamic_pointer_cast<const VariableNode>(pw->exponent());
    if (!exp_node || exp_node->name() != n) return false;
    auto base_expr = LMCAS::detail::make_expression_ptr(pw->base());
    if (te_depends_on(base_expr, n)) return false;
    base_out = base_expr;
    return true;
}

/**
 * @internal
 * @brief 检测 sin(w*n) 或 cos(w*n) 形式并提取角频率 w。
 */
static bool zt_is_trig_seq(const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
                            bool& is_sin_out, std::shared_ptr<SymbolicExpr>& omega_out) {
    if (!f || !LMCAS::detail::node(f)) return false;
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(f));
    if (!fn || fn->arguments().empty()) return false;
    if (fn->type() != FunctionNode::FuncType::Sin && fn->type() != FunctionNode::FuncType::Cos) return false;
    is_sin_out = (fn->type() == FunctionNode::FuncType::Sin);
    auto arg = fn->arguments()[0];
    auto ml = std::dynamic_pointer_cast<const MultiplyNode>(arg);
    if (ml && ml->operands().size() == 2) {
        for (size_t i = 0; i < 2; ++i) {
            auto vn = std::dynamic_pointer_cast<const VariableNode>(ml->operands()[i]);
            auto other = LMCAS::detail::make_expression_ptr(ml->operands()[1 - i]);
            if (vn && vn->name() == n && !te_depends_on(other, n)) {
                omega_out = other;
                return true;
            }
        }
    }
    auto vn = std::dynamic_pointer_cast<const VariableNode>(arg);
    if (vn && vn->name() == n) {
        omega_out = SymbolicExpr::number(1);
        return true;
    }
    return false;
}

static std::shared_ptr<SymbolicExpr> zt_abs_expr(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return nullptr;
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(expr)}));
}

static std::vector<std::shared_ptr<SymbolicExpr>> zt_unit_roc(const std::string& z) {
    auto condition = te_gt_condition(zt_abs_expr(SymbolicExpr::variable(z)),
                                     SymbolicExpr::number(1));
    return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                     : std::vector<std::shared_ptr<SymbolicExpr>>{};
}

static std::vector<std::shared_ptr<SymbolicExpr>> zt_base_roc(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::string& z) {
    auto condition = te_gt_condition(zt_abs_expr(SymbolicExpr::variable(z)),
                                     zt_abs_expr(base));
    return condition ? std::vector<std::shared_ptr<SymbolicExpr>>{condition}
                     : std::vector<std::shared_ptr<SymbolicExpr>>{};
}

static std::vector<std::shared_ptr<SymbolicExpr>> zt_roc(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z) {
    if (!f_n || !LMCAS::detail::node(f_n)) return {};

    auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f_n));
    if (add_node) {
        std::vector<std::shared_ptr<SymbolicExpr>> combined;
        for (const auto& op : add_node->operands()) {
            auto term_roc = zt_roc(LMCAS::detail::make_expression_ptr(op), n, z);
            combined.insert(combined.end(), term_roc.begin(), term_roc.end());
        }
        return combined;
    }

    if (!te_depends_on(f_n, n)) return zt_unit_roc(z);

    auto [coeff, body] = te_split_coeff(f_n, n);
    if (!coeff->is_one()) return zt_roc(body, n, z);

    auto vn = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(f_n));
    if (vn && vn->name() == n) return zt_unit_roc(z);

    auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(f_n));
    if (pw) {
        auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
        double ev = te_numval(pw->exponent());
        if (bv && bv->name() == n && std::abs(ev - 2.0) < 1e-12) {
            return zt_unit_roc(z);
        }
    }

    std::shared_ptr<SymbolicExpr> base_expr;
    if (zt_is_exp_seq(f_n, n, base_expr)) return zt_base_roc(base_expr, z);

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(f_n));
    if (mul && mul->operands().size() == 2) {
        for (size_t i = 0; i < 2; ++i) {
            auto term_n = std::dynamic_pointer_cast<const VariableNode>(mul->operands()[i]);
            auto other = LMCAS::detail::make_expression_ptr(mul->operands()[1 - i]);
            if (term_n && term_n->name() == n &&
                zt_is_exp_seq(other, n, base_expr)) {
                return zt_base_roc(base_expr, z);
            }
        }
    }

    bool is_sin = false;
    std::shared_ptr<SymbolicExpr> omega_expr;
    if (zt_is_trig_seq(f_n, n, is_sin, omega_expr)) return zt_unit_roc(z);

    return {};
}

static std::shared_ptr<SymbolicExpr> z_transform_core(
    const std::shared_ptr<SymbolicExpr>& f_n, const std::string& n,
    const std::string& z, ComputationContext& context) {
    auto step = context.consume_steps(1, "z_transform.recursive");
    if (!step) return nullptr;
    if (!f_n || !LMCAS::detail::node(f_n)) return nullptr;
    auto zv = SymbolicExpr::variable(z);

    /// 线性性：处理加法
    auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(f_n));
    if (add_node) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add_node->operands()) {
            auto zt = z_transform_core(
                LMCAS::detail::make_expression_ptr(op), n, z, context);
            if (!zt) {
                return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
                    TransformNode::TransformType::ZTransform, LMCAS::detail::node(f_n)->clone(), n,
                    SymbolicFactory::create_variable(z)));
            }
            result = result ? SymbolicExpr::add(result, zt) : zt;
        }
        return result;
    }

    /// 分离常数系数: c * body → c * Z{body}
    auto [coeff, body] = te_split_coeff(f_n, n);
    if (!coeff->is_one()) {
        auto zt_body = z_transform_core(body, n, z, context);
        if (zt_body) return SymbolicExpr::multiply(coeff, zt_body)->simplify();
    }

    /// Z{c} = c·z/(z-1)（常数序列，即阶跃函数 u[n]）
    if (!te_depends_on(f_n, n)) {
        auto den = SymbolicExpr::add(zv, SymbolicExpr::number(-1));
        return SymbolicExpr::multiply(f_n, SymbolicExpr::divide(zv, den))->simplify();
    }

    /// Z{n} = z/(z-1)^2
    {
        auto vn = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(f_n));
        if (vn && vn->name() == n) {
            auto den = SymbolicExpr::power(
                SymbolicExpr::add(zv, SymbolicExpr::number(-1)),
                SymbolicExpr::number(2));
            return SymbolicExpr::divide(zv, den)->simplify();
        }
    }

    /// Z{n^2} = z*(z+1)/(z-1)^3
    {
        auto pw = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(f_n));
        if (pw) {
            auto bv = std::dynamic_pointer_cast<const VariableNode>(pw->base());
            double ev = te_numval(pw->exponent());
            if (bv && bv->name() == n && std::abs(ev - 2.0) < 1e-12) {
                auto den = SymbolicExpr::power(
                    SymbolicExpr::add(zv, SymbolicExpr::number(-1)),
                    SymbolicExpr::number(3));
                auto num_expr = SymbolicExpr::multiply(zv,
                    SymbolicExpr::add(zv, SymbolicExpr::number(1)));
                return SymbolicExpr::divide(num_expr, den)->simplify();
            }
        }
    }

    /// Z{a^n} = z/(z-a)
    {
        std::shared_ptr<SymbolicExpr> base_expr;
        if (zt_is_exp_seq(f_n, n, base_expr)) {
            auto den = SymbolicExpr::add(zv,
                SymbolicExpr::multiply(SymbolicExpr::number(-1), base_expr));
            return SymbolicExpr::divide(zv, den)->simplify();
        }
    }

    /// Z{n·a^n} = a·z/(z-a)^2
    {
        auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(f_n));
        if (mul && mul->operands().size() == 2) {
            for (size_t i = 0; i < 2; ++i) {
                auto vn = std::dynamic_pointer_cast<const VariableNode>(mul->operands()[i]);
                auto other = LMCAS::detail::make_expression_ptr(mul->operands()[1 - i]);
                if (vn && vn->name() == n) {
                    std::shared_ptr<SymbolicExpr> base_expr;
                    if (zt_is_exp_seq(other, n, base_expr)) {
                        auto den = SymbolicExpr::power(
                            SymbolicExpr::add(zv,
                                SymbolicExpr::multiply(SymbolicExpr::number(-1), base_expr)),
                            SymbolicExpr::number(2));
                        auto num_expr = SymbolicExpr::multiply(base_expr, zv);
                        return SymbolicExpr::divide(num_expr, den)->simplify();
                    }
                }
            }
        }
    }

    /// Z{sin(w·n)} = z·sin(w)/(z^2 - 2z·cos(w) + 1)
    /// Z{cos(w·n)} = z·(z-cos(w))/(z^2 - 2z·cos(w) + 1)
    {
        bool is_sin_s = false;
        std::shared_ptr<SymbolicExpr> omega_expr;
        if (zt_is_trig_seq(f_n, n, is_sin_s, omega_expr)) {
            auto cos_w = SymbolicExpr::cos(omega_expr);
            auto z_sq = SymbolicExpr::power(zv, SymbolicExpr::number(2));
            auto two_z_cos_w = SymbolicExpr::multiply(
                SymbolicExpr::multiply(SymbolicExpr::number(2), zv), cos_w);
            auto den = SymbolicExpr::add(
                SymbolicExpr::add(z_sq,
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), two_z_cos_w)),
                SymbolicExpr::number(1));
            if (is_sin_s) {
                auto num_expr = SymbolicExpr::multiply(zv, SymbolicExpr::sin(omega_expr));
                return SymbolicExpr::divide(num_expr, den)->simplify();
            } else {
                auto z_minus_cos = SymbolicExpr::add(zv,
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), cos_w));
                auto num_expr = SymbolicExpr::multiply(zv, z_minus_cos);
                return SymbolicExpr::divide(num_expr, den)->simplify();
            }
        }
    }

    /// 未知形式：返回未求值节点
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::ZTransform, LMCAS::detail::node(f_n)->clone(), n,
        SymbolicFactory::create_variable(z)));
}

TransformEngineResult z_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z,
    ComputationContext& context) {
    const std::string operation = "z_transform";
    auto valid = te_validate_expr_vars(f_n, n, z, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = z_transform_core(f_n, n, z, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        auto result = te_wrap_transform_result(std::move(expression), operation);
        if (result) {
            result.value().value.roc = zt_roc(f_n, n, z);
        }
        return result;
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult z_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f_n,
    const std::string& n,
    const std::string& z) {
    ComputationContext context;
    return z_transform_checked(f_n, n, z, context);
}


} // namespace LMCAS
