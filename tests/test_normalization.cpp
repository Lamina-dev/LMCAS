#include <memory>
#include <vector>
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

using namespace LMCAS;

std::shared_ptr<const SymbolicNode> normalize(std::shared_ptr<const SymbolicNode> node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

TEST(LmcasNormalization, ZeroProductNormalization) {
    auto x = LMCAS::detail::make_node<VariableNode>("x");
    std::vector<std::shared_ptr<const SymbolicNode>> zero_ops;
    zero_ops.push_back(x);
    zero_ops.push_back(LMCAS::detail::make_node<NumberNode>(0.0));
    auto expr3 = LMCAS::detail::make_node<MultiplyNode>(std::move(zero_ops));

    auto norm3 = normalize(expr3);
    EXPECT_TRUE((norm3 != nullptr && norm3->is_zero())) << ("multiplication by zero normalizes to zero");
}
