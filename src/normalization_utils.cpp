#include "internal/normalization_utils.hpp"

#include <limits>

namespace LMCAS {

void normalization_append(detail::RewriteBudget* budget, std::size_t& nodes,
    std::vector<std::shared_ptr<const SymbolicNode>>& children,
    std::shared_ptr<const SymbolicNode> child) {
    if (budget) {
        nodes = budget->append_size(nodes, budget->measure(child));
    }
    children.push_back(std::move(child));
}

std::size_t normalization_check_count(detail::RewriteBudget* budget,
    std::size_t rows, std::size_t cols, bool wrapped) {
    if (budget) {
        const auto maximum = std::numeric_limits<std::size_t>::max();
        if (cols != 0 && rows > maximum / cols) {
            budget->append_size(maximum, 1);
        }
        budget->append_size(wrapped ? 1 : 0, rows * cols);
    }
    return rows * cols;
}

lmmc_real_t normalization_number_value(const NumberNode& node) {
    if (const auto* real = std::get_if<lmmc_real_t>(&node.value())) {
        return *real;
    }
    if (const auto* rational = std::get_if<Rational>(&node.value())) {
        return rational->to_double();
    }
    return std::get<BigInt>(node.value()).to_double();
}

bool normalization_factor_less(const std::shared_ptr<const SymbolicNode>& left,
                               const std::shared_ptr<const SymbolicNode>& right) {
    const BigInt a = get_node_degree_helper(left);
    const BigInt b = get_node_degree_helper(right);
    if (a != b) {
        return a > b;
    }
    return left->compare(*right) < 0;
}

bool try_get_integer_value(
    const std::shared_ptr<const NumberNode>& node, BigInt& value) {
    if (!node) {
        return false;
    }
    if (const auto* integer = std::get_if<BigInt>(&node->value())) {
        value = *integer;
        return true;
    }
    if (const auto* rational = std::get_if<Rational>(&node->value())) {
        if (rational->get_denominator() != BigInt(1)) {
            return false;
        }
        value = rational->get_numerator();
        return true;
    }
    const auto real = std::get<lmmc_real_t>(node->value());
    if (!std::isfinite(real) || std::trunc(real) != real) {
        return false;
    }
    value = Rational::from_double(real).get_numerator();
    return true;
}

/**
 * @brief 计算 AST 节点的多项式次数
 * @param node 输入节点
 * @return 节点对应的多项式次数
 */
BigInt get_node_degree_helper(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return BigInt(0);
    }
    if (std::dynamic_pointer_cast<const VariableNode>(node)) {
        return BigInt(1);
    }
    if (auto p = std::dynamic_pointer_cast<const PowerNode>(node)) {
        const auto e = std::dynamic_pointer_cast<const NumberNode>(p->exponent());
        BigInt degree;
        if (try_get_integer_value(e, degree)) {
            return degree;
        }
        return BigInt(1);
    }
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        BigInt d(0);
        for (auto& op : m->operands()) d += get_node_degree_helper(op);
        return d;
    }
    return BigInt(0);
}
std::shared_ptr<const SymbolicNode> make_normalized_multiply_node(
    const std::vector<std::shared_ptr<const SymbolicNode>>& ops,
    detail::RewriteBudget* budget) {
    std::size_t nodes = 0;
    if (budget) {
        for (const auto& op : ops) {
            nodes = budget->append_size(nodes, budget->measure(op));
        }
    }
    normalization_check_arithmetic<MultiplyNode>(budget, ops, nodes, ops.size() != 1);
    if (ops.empty()) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    }
    if (ops.size() == 1) {
        return ops[0];
    }
    return LMCAS::detail::make_node<MultiplyNode>(ops);
}
std::shared_ptr<const SymbolicNode> norm_subst_index(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& index_var,
    const std::shared_ptr<const SymbolicNode>& value, detail::RewriteBudget* budget) {
    normalization_check_children(budget, 0, node);
    return LMCAS::substitute_free(node, index_var, value, budget);
}
bool NodeCompare::operator()(const std::shared_ptr<const SymbolicNode>& lhs, const std::shared_ptr<const SymbolicNode>& rhs) const {
        if (!lhs && !rhs) {
            return false;
        }
        if (!lhs) {
            return true;
        }
        if (!rhs) {
            return false;
        }

        const BigInt d1 = get_node_degree_helper(lhs);
        const BigInt d2 = get_node_degree_helper(rhs);
        if (d1 != d2) {
            return d1 > d2;
        }

        bool isNum1 = std::dynamic_pointer_cast<const NumberNode>(lhs) != nullptr;
        bool isNum2 = std::dynamic_pointer_cast<const NumberNode>(rhs) != nullptr;
        if (isNum1 != isNum2) {
            return isNum1;
        }

        if (lhs->type_priority() != rhs->type_priority()) {
            return lhs->type_priority() < rhs->type_priority();
        }

        return lhs->compare(*rhs) < 0;
    }
bool is_approximate_number(const NumberNode& node) {
    const auto* real = std::get_if<lmmc_real_t>(&node.value());
    return real && (!std::isfinite(*real) || std::trunc(*real) != *real);
}

namespace {
class InexactNumberDetector final : public detail::RecursiveSymbolicVisitor {
public:
    bool found = false;

    void visit(const NumberNode& number) override {
        found = found ||
            std::holds_alternative<lmmc_real_t>(number.value());
    }
};
}

bool contains_inexact_number(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    InexactNumberDetector detector;
    node->accept(detector);
    return detector.found;
}

Rational exact_number_as_rational(const NumberNode& node) {
    if (const auto* rational = std::get_if<Rational>(&node.value())) {
        return *rational;
    }
    if (const auto* integer = std::get_if<BigInt>(&node.value())) {
        return Rational(*integer);
    }
    return Rational::from_double(std::get<lmmc_real_t>(node.value()));
}

/**
 * @brief 两个数值节点相加
 * @param a 加数节点
 * @param b 加数节点
 * @return 和的数值节点
 */
std::shared_ptr<const NumberNode> add_numbers(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
     if (is_approximate_number(*a) || is_approximate_number(*b)) {
         lmmc_real_t v1 = normalization_number_value(*a);
         lmmc_real_t v2 = normalization_number_value(*b);
         lmmc_real_t sum = v1 + v2;
         return LMCAS::detail::make_node<NumberNode>(sum);
     }

     if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
         const Rational r1 = exact_number_as_rational(*a);
         const Rational r2 = exact_number_as_rational(*b);
         return LMCAS::detail::make_node<NumberNode>(r1 + r2);
     }

     BigInt i1;
     BigInt i2;
     try_get_integer_value(a, i1);
     try_get_integer_value(b, i2);
     return LMCAS::detail::make_node<NumberNode>(i1 + i2);
}

/**
 * @brief 两个数值节点相乘
 * @param a 乘数节点
 * @param b 乘数节点
 * @return 积的数值节点
 */
std::shared_ptr<const NumberNode> multiply_numbers(const std::shared_ptr<const NumberNode>& a, const std::shared_ptr<const NumberNode>& b) {
     if (is_approximate_number(*a) || is_approximate_number(*b)) {
         lmmc_real_t v1 = normalization_number_value(*a);
         lmmc_real_t v2 = normalization_number_value(*b);
         lmmc_real_t prod = v1 * v2;
         return LMCAS::detail::make_node<NumberNode>(prod);
     }

     if (std::holds_alternative<Rational>(a->value()) || std::holds_alternative<Rational>(b->value())) {
         const Rational r1 = exact_number_as_rational(*a);
         const Rational r2 = exact_number_as_rational(*b);
         return LMCAS::detail::make_node<NumberNode>(r1 * r2);
     }

     BigInt i1;
     BigInt i2;
     try_get_integer_value(a, i1);
     try_get_integer_value(b, i2);
     return LMCAS::detail::make_node<NumberNode>(i1 * i2);
}
bool get_pi_coeff(const std::shared_ptr<const SymbolicNode>& node, Rational& k) {
    if (const auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (variable->name() == "pi") {
            k = Rational(1);
            return true;
        }
    }
    const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) {
        return false;
    }
    bool has_pi = false;
    k = Rational(1);
    for (const auto& operand : multiply->operands()) {
        if (const auto variable = std::dynamic_pointer_cast<const VariableNode>(operand)) {
            if (variable->name() != "pi" || has_pi) {
                return false;
            }
            has_pi = true;
            continue;
        }
        const auto number = std::dynamic_pointer_cast<const NumberNode>(operand);
        if (!number) {
            return false;
        }
        if (const auto* rational = std::get_if<Rational>(&number->value())) {
            k = k * *rational;
            continue;
        }
        const auto* integer = std::get_if<BigInt>(&number->value());
        if (!integer) {
            return false;
        }
        k = k * Rational(*integer);
    }
    return has_pi;
}

}
