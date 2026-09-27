#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"

using namespace LMCAS;

static SymbolicExpr nested_exponential(const std::string &variable) {
    std::shared_ptr<const SymbolicNode> node = detail::make_node<VariableNode>(variable);
    for (int i = 0; i < 8; ++i) {
        node = detail::make_node<FunctionNode>(FunctionNode::FuncType::Exp,
                                               std::vector<std::shared_ptr<const SymbolicNode>>{node});
    }
    return detail::expression_from_node(node);
}

TEST(AssumptionCycleDepth, SharedDagProofs) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_sign("x", Sign::Positive).has_value())) << "positive variable declared";
    InferenceEngine engine(context);
    auto x = detail::make_node<VariableNode>("x");
    auto shared = detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, x});
    auto expression = detail::expression_from_node(detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{shared, shared}));
    auto positive = engine.query_positive_checked(expression);
    auto negative = engine.query_negative_checked(expression);
    auto real = engine.query_real_checked(expression);
    EXPECT_TRUE((positive && positive.value() == Tribool::True)) << "shared positive factors prove a positive product";
    EXPECT_TRUE((negative && negative.value() == Tribool::False)) << "the product is not negative";
    EXPECT_TRUE((real && real.value() == Tribool::True)) << "independent property queries retain realness";
    auto repeated = engine.query_positive_checked(expression);
    EXPECT_TRUE((repeated && repeated.value() == Tribool::True)) << "completed traversal does not poison a later proof";

    auto distinct = detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::make_node<VariableNode>("x"), detail::make_node<VariableNode>("x")});
    auto equal_operands = detail::expression_from_node(detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{shared, distinct}));
    auto equal_product = engine.query_positive_checked(equal_operands);
    EXPECT_TRUE((equal_product && equal_product.value() == Tribool::True)) << "structurally equal operands with different identities remain provable";
}

TEST(AssumptionCycleDepth, ProofDepthRecovery) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_domain("x", Domain::Real).has_value())) << "real variable declared";
    InferenceEngine engine(context);
    const auto expression = nested_exponential("x");
    engine.set_max_depth(1);
    auto limited = engine.query_positive_checked(expression);
    EXPECT_TRUE((limited && limited.value() == Tribool::Unknown)) << "insufficient proof recursion does not assert a value";
    auto shallow = engine.query_real_checked(*SymbolicExpr::variable("x"));
    EXPECT_TRUE((shallow && shallow.value() == Tribool::True)) << "an exhausted proof does not block an independent shallow proof";
    engine.set_max_depth(32);
    auto proved = engine.query_positive_checked(expression);
    EXPECT_TRUE((proved && proved.value() == Tribool::True)) << "sufficient proof recursion establishes positivity";
    engine.set_max_depth(1);
    auto limited_again = engine.query_positive_checked(expression);
    EXPECT_TRUE((limited_again && limited_again.value() == Tribool::Unknown)) << "a previous deeper proof does not bypass the current recursion policy";
    engine.set_max_depth(0);
    auto invalid_policy = engine.query_positive_checked(expression);
    EXPECT_TRUE((invalid_policy && invalid_policy.value() == Tribool::Unknown)) << "an invalid zero depth does not silently disable the existing bound";
}

TEST(AssumptionCycleDepth, DepthExhaustionPreservesScopedAssumptions) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_sign("x", Sign::Positive).has_value())) << "parent positive assumption declared";
    context.push();
    EXPECT_TRUE((context.assume_sign("x", Sign::Negative).has_value())) << "child negative assumption shadows parent";
    InferenceEngine engine(context);
    engine.set_max_depth(1);
    auto limited = engine.query_positive_checked(nested_exponential("x"));
    EXPECT_TRUE((limited && limited.value() == Tribool::Unknown)) << "child proof reaches the recursion bound";
    auto child = engine.query_positive_checked(*SymbolicExpr::variable("x"));
    EXPECT_TRUE((child && child.value() == Tribool::False)) << "exhaustion does not erase the child assumption";
    EXPECT_TRUE((context.pop().has_value())) << "the original child scope remains available to pop";
    auto parent = engine.query_positive_checked(*SymbolicExpr::variable("x"));
    EXPECT_TRUE((parent && parent.value() == Tribool::True)) << "parent semantics are restored after the exhausted query";
}
