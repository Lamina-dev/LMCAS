#include "internal/query_support.hpp"
#include <cmath>

namespace LMCAS {


bool QueryInterface::is_unhandled_type(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    if (!node) return true;
    if (std::dynamic_pointer_cast<const MatrixNode>(node)) return true;
    if (std::dynamic_pointer_cast<const RelationalNode>(node)) return true;
    if (std::dynamic_pointer_cast<const LogicalNode>(node)) return true;
    return false;
}

bool QueryInterface::is_nan_number(const SymbolicExpr& expression) const {
    const auto node = std::dynamic_pointer_cast<const NumberNode>(
        LMCAS::detail::node(expression));
    if (node && std::holds_alternative<lmmc_real_t>(node->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(node->value());
        return std::isnan(v);
    }
    return false;
}

bool QueryInterface::is_infinity_node(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return func->type() == FunctionNode::FuncType::Infinity;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        if (mul->operands().size() == 2) {
            for (const auto& op : mul->operands()) {
                if (auto func = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                    if (func->type() == FunctionNode::FuncType::Infinity) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int QueryInterface::get_infinity_sign(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Infinity) {
            return +1;
        }
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        bool has_infinity = false;
        bool has_neg_one = false;
        for (const auto& op : mul->operands()) {
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                if (func->type() == FunctionNode::FuncType::Infinity) {
                    has_infinity = true;
                }
            }
            const auto number =
                std::dynamic_pointer_cast<const NumberNode>(op);
            if (number && number->is_negative_one()) has_neg_one = true;
        }
        if (has_infinity && has_neg_one) return -1;
        if (has_infinity) return +1; /** @brief 无穷与正因子相乘仍为正无穷。 */
    }
    return 0; /** @brief 符号未定。 */
}

}
