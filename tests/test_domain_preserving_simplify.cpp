#include "test_common.hpp"
#include "../include/visitors/normalization_visitor.hpp"

static std::shared_ptr<const SymbolicNode> normalize_node(const std::shared_ptr<const SymbolicNode>& node) {
    NormalizationVisitor visitor;
    node->accept(visitor);
    return visitor.get_result();
}

static bool is_bigint_value(const std::shared_ptr<const SymbolicNode>& node, long long expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    return number &&
           std::holds_alternative<BigInt>(number->value()) &&
           std::get<BigInt>(number->value()) == BigInt(expected);
}

int main() {
    TEST_CASE("simplify does not distribute products over sums");
    auto x = lamina::detail::make_node<VariableNode>("x");
    auto y = lamina::detail::make_node<VariableNode>("y");
    auto y_plus_one = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{y, lamina::detail::make_node<NumberNode>(BigInt(1))});
    auto product = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, y_plus_one});
    auto simplified_product = normalize_node(product);
    EXPECT_TRUE(std::dynamic_pointer_cast<const MultiplyNode>(simplified_product) != nullptr,
                "simplify keeps x*(y+1) as a product");

    TEST_CASE("zero exponent requires a proved nonzero base");
    auto unknown_zero_power = lamina::detail::make_node<PowerNode>(
        lamina::detail::make_node<VariableNode>("x"), lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto simplified_unknown_zero_power = normalize_node(unknown_zero_power);
    EXPECT_TRUE(std::dynamic_pointer_cast<const PowerNode>(simplified_unknown_zero_power) != nullptr,
                "simplify keeps x^0 without a nonzero assumption");

    auto numeric_zero_power = lamina::detail::make_node<PowerNode>(
        lamina::detail::make_node<NumberNode>(BigInt(2)), lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto simplified_numeric_zero_power = normalize_node(numeric_zero_power);
    EXPECT_TRUE(is_bigint_value(simplified_numeric_zero_power, 1),
                "simplify folds 2^0 because the base is proved nonzero");

    TEST_CASE("multiplicative inverse cancellation is not unconditional");
    auto inv_x = lamina::detail::make_node<PowerNode>(
        lamina::detail::make_node<VariableNode>("x"), lamina::detail::make_node<NumberNode>(BigInt(-1)));
    auto x_times_inv_x = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::make_node<VariableNode>("x"), inv_x});
    auto simplified_inverse = normalize_node(x_times_inv_x);
    EXPECT_FALSE(simplified_inverse->is_one(),
                 "simplify does not turn x*x^-1 into 1 without x != 0");

    TEST_CASE("safe positive integer exponent merging remains available");
    auto x_squared = lamina::detail::make_node<PowerNode>(
        lamina::detail::make_node<VariableNode>("x"), lamina::detail::make_node<NumberNode>(BigInt(2)));
    auto x_times_x_squared = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::make_node<VariableNode>("x"), x_squared});
    auto simplified_positive_merge = normalize_node(x_times_x_squared);
    auto merged_power = std::dynamic_pointer_cast<const PowerNode>(simplified_positive_merge);
    EXPECT_TRUE(merged_power != nullptr, "x*x^2 merges to a power");
    EXPECT_TRUE(merged_power && is_bigint_value(merged_power->exponent(), 3),
                "x*x^2 merges to x^3");

    return TEST_REPORT();
}
