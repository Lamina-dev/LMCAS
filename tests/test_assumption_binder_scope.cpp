#include "test_common.hpp"
#include "assumption_context.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> normalize_with_ctx(
    const std::shared_ptr<const SymbolicNode> &node,
    const AssumptionContext &ctx) {
    auto expression = LMCAS::detail::expression_from_node(node);
    return LMCAS::detail::node(ctx.simplify(expression));
}

static std::shared_ptr<const SymbolicNode> normalize_no_ctx(
    const std::shared_ptr<const SymbolicNode> &node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

static std::shared_ptr<const SymbolicNode> var(const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> num(int v) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(v));
}

static std::shared_ptr<const SymbolicNode> make_abs(const std::shared_ptr<const SymbolicNode> &arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

TEST(AssumptionBinderScope, FiniteBinderAssumptionShadowing) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("i", Sign::Positive).has_value()) << "outer i is positive";

    auto sum = detail::make_node<SummationNode>(make_abs(var("i")), "i", num(-1), num(1));
    auto renamed_sum = detail::make_node<SummationNode>(make_abs(var("k")), "k", num(-1), num(1));
    EXPECT_TRUE((normalize_no_ctx(sum)->equals(*num(2)))) << "sum(abs(i),i,-1,1) is two";
    EXPECT_TRUE((normalize_with_ctx(sum, ctx)->equals(*num(2)))) << "outer positivity cannot change a finite sum's bound values";
    EXPECT_TRUE((normalize_with_ctx(sum, ctx)->equals(*normalize_with_ctx(renamed_sum, ctx)))) << "alpha-renaming a summation index preserves its value";

    auto product = detail::make_node<ProductNode>(make_abs(var("i")), "i", num(-3), num(-1));
    auto renamed_product = detail::make_node<ProductNode>(make_abs(var("k")), "k", num(-3), num(-1));
    EXPECT_TRUE((normalize_no_ctx(product)->equals(*num(6)))) << "product of three absolute values is six";
    EXPECT_TRUE((normalize_with_ctx(product, ctx)->equals(*num(6)))) << "outer positivity cannot change a finite product's sign";
    EXPECT_TRUE((normalize_with_ctx(product, ctx)->equals(*normalize_with_ctx(renamed_product, ctx)))) << "alpha-renaming a product index preserves its value";

    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "unrelated x is positive";
    auto weighted = detail::make_node<SummationNode>(
        SymbolicFactory::create_multiply({make_abs(var("i")), make_abs(var("x"))}),
        "i", num(-1), num(1));
    auto expected = normalize_no_ctx(SymbolicFactory::create_multiply({num(2), var("x")}));
    EXPECT_TRUE((normalize_with_ctx(weighted, ctx)->equals(*expected))) << "unrelated free-symbol facts remain applicable inside a finite binder";
}

TEST(AssumptionBinderScope, NestedBinderAssumptionShadowing) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("i", Sign::Positive).has_value()) << "outer i is positive";
    ASSERT_TRUE(ctx.assume_sign("j", Sign::Negative).has_value()) << "outer j is negative";

    auto inner = detail::make_node<SummationNode>(
        make_abs(var("i")), "i", num(-1), make_abs(var("i")));
    auto nested = detail::make_node<SummationNode>(inner, "i", num(-1), num(1));
    EXPECT_TRUE((normalize_with_ctx(nested, ctx)->equals(*num(5)))) << "nested same-name bodies and outer-index bounds retain distinct scopes";

    auto mixed_inner = detail::make_node<ProductNode>(
        SymbolicFactory::create_add({make_abs(var("i")), make_abs(var("j"))}),
        "j", num(-1), num(1));
    auto mixed = detail::make_node<SummationNode>(mixed_inner, "i", num(-1), num(1));
    EXPECT_TRUE((normalize_with_ctx(mixed, ctx)->equals(*num(8)))) << "different nested binders mask both enclosing declarations";
}

TEST(AssumptionBinderScope, AllBinderBodyFactScopes) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("i", Sign::Positive).has_value()) << "outer i is positive";
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "free x is positive";
    auto body = SymbolicFactory::create_add({make_abs(var("i")), make_abs(var("x"))});
    auto expected_body = SymbolicFactory::create_add({make_abs(var("i")), var("x")});
    auto outside = make_abs(var("i"));
    auto expected_outside = var("i");
    auto domain = detail::make_node<FiniteSetNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{outside});
    auto expected_domain = detail::make_node<FiniteSetNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{expected_outside});
    auto predicate = detail::make_node<RelationalNode>(body, num(0), RelationOp::GT);
    auto expected_predicate = detail::make_node<RelationalNode>(expected_body, num(0), RelationOp::GT);
    using Node = std::shared_ptr<const SymbolicNode>;
    const std::vector<std::pair<Node, Node>> cases{
        {detail::make_node<SummationNode>(body, "i", num(0), outside),
         detail::make_node<SummationNode>(expected_body, "i", num(0), expected_outside)},
        {detail::make_node<ProductNode>(body, "i", num(0), outside),
         detail::make_node<ProductNode>(expected_body, "i", num(0), expected_outside)},
        {detail::make_node<IntegralNode>(body, "i", num(0), outside),
         detail::make_node<IntegralNode>(expected_body, "i", num(0), expected_outside)},
        {detail::make_node<IntegralNode>(body, "i"),
         detail::make_node<IntegralNode>(expected_body, "i")},
        {detail::make_node<TransformNode>(TransformNode::TransformType::Laplace, body, "i", outside),
         detail::make_node<TransformNode>(TransformNode::TransformType::Laplace, expected_body, "i", expected_outside)},
        {detail::make_node<QuantifierNode>(QuantifierNode::Type::ForAll, "i", domain, predicate),
         detail::make_node<QuantifierNode>(QuantifierNode::Type::ForAll, "i", expected_domain, expected_predicate)},
        {detail::make_node<SetBuilderNode>("i", domain, predicate),
         detail::make_node<SetBuilderNode>("i", expected_domain, expected_predicate)},
        {detail::make_node<LimitNode>(body, "i", outside, LimitDirection::Both),
         detail::make_node<LimitNode>(expected_body, "i", expected_outside, LimitDirection::Both)}};
    for (const auto &[source, expected] : cases) {
        EXPECT_TRUE((normalize_with_ctx(source, ctx)->equals(*normalize_no_ctx(expected)))) << "bound absolute values stay unresolved while free facts simplify body and outside children";
    }
}

TEST(AssumptionBinderScope, BinderRelationDerivedFacts) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("i", Sign::Positive).has_value()) << "logarithm has a positive outer argument";
    auto greater_than_one = detail::expression_from_node(
        detail::make_node<RelationalNode>(var("i"), num(1), RelationOp::GT));
    EXPECT_TRUE((ctx.assume(greater_than_one).has_value())) << "outer i is greater than one";
    auto logarithm = detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<const SymbolicNode>>{var("i")});
    auto absolute_logarithm = make_abs(logarithm);
    EXPECT_TRUE((normalize_with_ctx(absolute_logarithm, ctx)->equals(*logarithm))) << "the outer relation observably proves log(i) positive";
    auto integral = detail::make_node<IntegralNode>(absolute_logarithm, "i", num(1), num(2));
    EXPECT_TRUE((normalize_with_ctx(integral, ctx)->equals(*normalize_no_ctx(integral)))) << "a composite query cannot use an outer relation about its bound argument";

    auto free_integral = detail::make_node<IntegralNode>(absolute_logarithm, "x", num(1), num(2));
    auto expected_free = detail::make_node<IntegralNode>(logarithm, "x", num(1), num(2));
    EXPECT_TRUE((normalize_with_ctx(free_integral, ctx)->equals(*expected_free))) << "relation-derived facts about unrelated free names remain available";
}
