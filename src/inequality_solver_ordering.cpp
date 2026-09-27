#include "internal/inequality_solver_support.hpp"
#include "internal/symbolic_ast.hpp"
#include <cmath>

namespace LMCAS::detail::inequality_support {
namespace {

int root_number_sign(const NumberNode& number) {
    const auto& value = number.value();
    if (const auto* integer = std::get_if<BigInt>(&value)) {
        if (integer->is_zero()) {
            return 0;
        }
        return integer->is_negative() ? -1 : 1;
    }
    if (const auto* rational = std::get_if<Rational>(&value)) {
        if (rational->get_numerator().is_zero()) {
            return 0;
        }
        return rational->get_numerator().is_negative() ? -1 : 1;
    }
    if (const auto* real = std::get_if<lmmc_real_t>(&value)) {
        if (*real == 0.0) {
            return 0;
        }
        return *real < 0 ? -1 : 1;
    }
    return 0;
}

int root_power_sign(const PowerNode& power) {
    auto exp_expr = LMCAS::detail::make_expression_ptr(power.exponent());
    if (auto ev = try_checked_numeric_constant(*exp_expr)) {
        if (*ev > 0 && *ev < 1.0) {
            return 1;
        }
    }
    auto base_expr = LMCAS::detail::make_expression_ptr(power.base());
    auto bv_checked = try_checked_numeric_constant(*base_expr);
    auto ev_checked = try_checked_numeric_constant(*exp_expr);
    if (!bv_checked || !ev_checked) {
        return 0;
    }
    double bv = *bv_checked;
    double ev = *ev_checked;
    int ei = static_cast<int>(ev);
    if (std::abs(ev - ei) >= 1e-10) {
        return 0;
    }
    if (bv > 0) {
        return 1;
    }
    if (bv < 0) {
        return (ei % 2 == 0) ? 1 : -1;
    }
    return 0;
}

int root_node_sign(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return 0;
    }
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return root_number_sign(*number);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Sqrt) {
            return 1;
        }
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return root_power_sign(*power);
    }
    return 0;
}

int root_expression_sign(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return 0;
    }
    int direct = root_node_sign(LMCAS::detail::node(expr));
    if (direct != 0) {
        return direct;
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr));
    if (!product) {
        return 0;
    }
    int sign = 1;
    for (const auto& operand : product->operands()) {
        int factor_sign = root_node_sign(operand);
        if (factor_sign == 0) {
            return 0;
        }
        sign *= factor_sign;
    }
    return sign;
}

}

bool root_less_than(const std::shared_ptr<SymbolicExpr>& a,
                           const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) {
        return false;
    }
    auto va = try_checked_numeric_constant(*a);
    auto vb = try_checked_numeric_constant(*b);
    if (va && vb) {
        return *va < *vb;
    }
    auto diff = SymbolicExpr::add(a, SymbolicExpr::multiply(b, SymbolicExpr::number(-1)))->simplify();
    if (auto vd = try_checked_numeric_constant(*diff)) {
        return *vd < 0;
    }
    return root_expression_sign(diff) < 0;
}

bool roots_equal(const std::shared_ptr<SymbolicExpr>& a,
                 const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) {
        return false;
    }
    auto numeric_a = try_checked_numeric_constant(*a);
    auto numeric_b = try_checked_numeric_constant(*b);
    if (numeric_a && numeric_b) {
        return *numeric_a == *numeric_b;
    }
    auto difference = SymbolicExpr::add(
        a,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), b));
    auto simplified = difference ? difference->simplify() : nullptr;
    return simplified && simplified->is_zero();
}
}
