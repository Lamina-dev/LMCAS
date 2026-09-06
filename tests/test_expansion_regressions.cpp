#include <memory>
#include <vector>

#include "../include/symbolic_ast.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

using namespace LMCAS;

int main() {
    auto a = LMCAS::detail::make_node<VariableNode>("a");
    auto b = LMCAS::detail::make_node<VariableNode>("b");
    auto a_plus_b = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{a, b});

    auto c = LMCAS::detail::make_node<VariableNode>("c");
    auto d = LMCAS::detail::make_node<VariableNode>("d");
    auto c_plus_d = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{c, d});

    NormalizationVisitor visitor;

    TEST_CASE("Expansion distributes one sum across a factor");
    auto one_sum = visitor.expand_product(a_plus_b, c);
    auto one_sum_add =
        std::dynamic_pointer_cast<const AddNode>(one_sum);
    EXPECT_TRUE(one_sum_add && one_sum_add->operands().size() == 2,
                "(a+b)*c expands to two terms");

    TEST_CASE("Expansion distributes two sums");
    auto two_sums = visitor.expand_product(a_plus_b, c_plus_d);
    auto two_sums_add =
        std::dynamic_pointer_cast<const AddNode>(two_sums);
    EXPECT_TRUE(two_sums_add && two_sums_add->operands().size() == 4,
                "(a+b)*(c+d) expands to four terms");

    return TEST_REPORT();
}
