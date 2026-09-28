#include "internal/transform_support.hpp"
#include "internal/normalization_utils.hpp"
#include <exception>

namespace LMCAS {
using namespace transform_detail;

static bool zt_is_exp_seq(const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
                           std::shared_ptr<SymbolicExpr>& base_out) {
    if (!f || !detail::node(f)) {
        return false;
    }
    auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(f));
    if (!power) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const VariableNode>(power->exponent());
    if (!exponent || exponent->is_constant() || exponent->name() != n) {
        return false;
    }
    auto base = detail::make_expression_ptr(power->base());
    if (te_depends_on(base, n)) {
        return false;
    }
    base_out = base;
    return true;
}

static bool zt_is_weighted_exp(const std::shared_ptr<SymbolicExpr>& f,
                              const std::string& n,
                              std::shared_ptr<SymbolicExpr>& base) {
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(detail::node(f));
    if (!product || product->operands().size() != 2) {
        return false;
    }
    for (std::size_t i = 0; i < 2; ++i) {
        auto variable = std::dynamic_pointer_cast<const VariableNode>(product->operands()[i]);
        auto other = detail::make_expression_ptr(product->operands()[1 - i]);
        if (variable && !variable->is_constant() && variable->name() == n &&
            zt_is_exp_seq(other, n, base)) {
            return true;
        }
    }
    return false;
}

enum class ZPolynomialPattern { None, Linear, Quadratic, Cubic };

static ZPolynomialPattern zt_polynomial_pattern(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(detail::node(f));
    if (variable && !variable->is_constant() && variable->name() == n) {
        return ZPolynomialPattern::Linear;
    }
    auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(f));
    if (!power) {
        return ZPolynomialPattern::None;
    }
    variable = std::dynamic_pointer_cast<const VariableNode>(power->base());
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    BigInt integer;
    if (variable && !variable->is_constant() && variable->name() == n && exponent &&
        try_get_integer_value(exponent, integer)) {
        if (integer == BigInt(2)) return ZPolynomialPattern::Quadratic;
        if (integer == BigInt(3)) return ZPolynomialPattern::Cubic;
    }
    return ZPolynomialPattern::None;
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

static std::vector<std::shared_ptr<SymbolicExpr>> zt_lookup_roc(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n, const std::string& z) {
    if (zt_polynomial_pattern(f, n) != ZPolynomialPattern::None) return zt_unit_roc(z);
    std::shared_ptr<SymbolicExpr> base;
    if (zt_is_exp_seq(f, n, base) || zt_is_weighted_exp(f, n, base)) return zt_base_roc(base, z);
    bool is_sin = false;
    std::shared_ptr<SymbolicExpr> frequency;
    if (te_is_trig(f, n, is_sin, frequency)) return zt_unit_roc(z);
    return {};
}

static std::vector<std::shared_ptr<SymbolicExpr>> zt_roc(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n, const std::string& z) {
    if (!f || !detail::node(f)) return {};
    auto addition = std::dynamic_pointer_cast<const AddNode>(detail::node(f));
    if (addition) {
        std::vector<std::shared_ptr<SymbolicExpr>> combined;
        for (const auto& operand : addition->operands()) {
            auto conditions = zt_roc(detail::make_expression_ptr(operand), n, z);
            combined.insert(combined.end(), conditions.begin(), conditions.end());
        }
        return combined;
    }
    if (!te_depends_on(f, n)) return zt_unit_roc(z);
    auto [coefficient, body] = te_split_coeff(f, n);
    if (!coefficient->is_one()) return zt_roc(body, n, z);
    return zt_lookup_roc(f, n, z);
}

static std::shared_ptr<SymbolicExpr> zt_polynomial(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
    const std::shared_ptr<SymbolicExpr>& z, const std::string& z_var) {
    auto pattern = zt_polynomial_pattern(f, n);
    if (pattern == ZPolynomialPattern::None) return nullptr;
    if (pattern == ZPolynomialPattern::Cubic) {
        auto square = SymbolicExpr::power(SymbolicExpr::variable(n), SymbolicExpr::number(2));
        auto quadratic = zt_polynomial(square, n, z, z_var);
        return SymbolicExpr::multiply(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), z),
            quadratic->differentiate(z_var))->simplify();
    }
    auto denominator = SymbolicExpr::power(SymbolicExpr::add(z, SymbolicExpr::number(-1)),
        SymbolicExpr::number(pattern == ZPolynomialPattern::Linear ? 2 : 3));
    auto numerator = pattern == ZPolynomialPattern::Linear ? z : SymbolicExpr::multiply(z,
        SymbolicExpr::add(z, SymbolicExpr::number(1)));
    return SymbolicExpr::divide(numerator, denominator)->simplify();
}

static std::shared_ptr<SymbolicExpr> zt_exponential(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
    const std::shared_ptr<SymbolicExpr>& z) {
    std::shared_ptr<SymbolicExpr> base;
    const bool plain = zt_is_exp_seq(f, n, base);
    if (!plain && !zt_is_weighted_exp(f, n, base)) return nullptr;
    auto difference = SymbolicExpr::add(z, SymbolicExpr::multiply(SymbolicExpr::number(-1), base));
    if (plain) return SymbolicExpr::divide(z, difference)->simplify();
    auto denominator = SymbolicExpr::power(difference, SymbolicExpr::number(2));
    return SymbolicExpr::divide(SymbolicExpr::multiply(base, z), denominator)->simplify();
}

static std::shared_ptr<SymbolicExpr> zt_trigonometric(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
    const std::shared_ptr<SymbolicExpr>& z) {
    bool is_sin = false;
    std::shared_ptr<SymbolicExpr> frequency;
    if (!te_is_trig(f, n, is_sin, frequency)) return nullptr;
    auto cosine = SymbolicExpr::cos(frequency);
    auto cross = SymbolicExpr::multiply(SymbolicExpr::multiply(SymbolicExpr::number(-2), z), cosine);
    auto denominator = SymbolicExpr::add(SymbolicExpr::add(
        SymbolicExpr::power(z, SymbolicExpr::number(2)), cross), SymbolicExpr::number(1));
    auto factor = is_sin ? SymbolicExpr::sin(frequency) : SymbolicExpr::add(z,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), cosine));
    return SymbolicExpr::divide(SymbolicExpr::multiply(z, factor), denominator)->simplify();
}

static std::shared_ptr<SymbolicExpr> zt_unevaluated(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n, const std::string& z) {
    return detail::make_expression_ptr(detail::make_node<TransformNode>(
        TransformNode::TransformType::ZTransform, detail::node(f)->clone(), n,
        SymbolicFactory::create_variable(z)));
}

static std::shared_ptr<SymbolicExpr> z_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
    const std::string& z, ComputationContext& context);

static std::shared_ptr<SymbolicExpr> zt_sum(
    const std::shared_ptr<SymbolicExpr>& f, const AddNode& addition,
    const std::string& n, const std::string& z, ComputationContext& context) {
    std::shared_ptr<SymbolicExpr> result;
    for (const auto& operand : addition.operands()) {
        auto transformed = z_transform_core(detail::make_expression_ptr(operand), n, z, context);
        if (!transformed) return zt_unevaluated(f, n, z);
        result = result ? SymbolicExpr::add(result, transformed) : transformed;
    }
    return result;
}

static std::shared_ptr<SymbolicExpr> z_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& n,
    const std::string& z, ComputationContext& context) {
    auto step = context.consume_steps(1, "z_transform.recursive");
    if (!step) {
        return nullptr;
    }
    if (!f || !detail::node(f)) {
        return nullptr;
    }
    auto zv = SymbolicExpr::variable(z);
    auto addition = std::dynamic_pointer_cast<const AddNode>(detail::node(f));
    if (addition) {
        return zt_sum(f, *addition, n, z, context);
    }
    auto [coefficient, body] = te_split_coeff(f, n);
    if (!coefficient->is_one()) {
        auto transformed = z_transform_core(body, n, z, context);
        if (transformed) {
            return SymbolicExpr::multiply(coefficient, transformed)->simplify();
        }
    }
    if (!te_depends_on(f, n)) {
        auto denominator = SymbolicExpr::add(zv, SymbolicExpr::number(-1));
        return SymbolicExpr::multiply(f, SymbolicExpr::divide(zv, denominator))->simplify();
    }
    auto polynomial = zt_polynomial(f, n, zv, z);
    if (polynomial) {
        return polynomial;
    }
    auto exponential = zt_exponential(f, n, zv);
    if (exponential) {
        return exponential;
    }
    auto trig = zt_trigonometric(f, n, zv);
    return trig ? trig : zt_unevaluated(f, n, z);
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

}
