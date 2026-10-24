#include <memory>
#include <vector>
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "test_common.hpp"
#include "internal/rewrite_budget.hpp"

namespace LMCAS::detail {
struct RewriteBudgetTestAccess {
    static std::size_t cache_size(const RewriteBudget& budget) {
        return budget.shapes_.size();
    }
    static void expired_entry(RewriteBudget& budget, const SymbolicNode* address) {
        budget.shapes_[address] = RewriteBudget::Entry{{}, {7, 3}};
    }
};
}

using namespace LMCAS;

std::shared_ptr<const SymbolicNode> normalize(std::shared_ptr<const SymbolicNode> node) {
    ComputationContext context;
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                 context.limits().max_ast_nodes, "test.normalize");
    NormalizationVisitor v(budget, detail::no_facts(), Domain::Real);
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

TEST(LmcasNormalization, RepeatedShapeMeasurementWork) {
    const auto work = [](std::size_t length) {
        ComputationContext context;
        detail::RewriteBudget budget(context, 192, 1024, "test.measure");
        detail::SymbolicNodePtr node = detail::make_node<VariableNode>("x");
        EXPECT_EQ(budget.measure(node), 1u);
        for (std::size_t i = 0; i < length; ++i) {
            node = detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Sin,
                std::vector<detail::SymbolicNodePtr>{node});
            EXPECT_EQ(budget.measure(node), i + 2);
        }
        return context.steps_used();
    };
    const auto small = work(40);
    const auto large = work(80);
    EXPECT_LE(large, 3 * small + 32);
}

TEST(LmcasNormalization, ShapeMeasurementCountsSharedEdges) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 7, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    auto sine = detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin, std::vector<detail::SymbolicNodePtr>{x});
    auto repeated = detail::make_node<AddNode>(
        std::vector<detail::SymbolicNodePtr>{sine, sine, sine});
    EXPECT_EQ(budget.measure(repeated), 7u);
    auto too_large = detail::make_node<AddNode>(
        std::vector<detail::SymbolicNodePtr>{sine, sine, sine, x});
    EXPECT_THROW(budget.measure(too_large), CasError);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(LmcasNormalization, MeasuredShapeObservesActiveDepth) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 3, 32, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    auto sine = detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin, std::vector<detail::SymbolicNodePtr>{x});
    EXPECT_EQ(budget.measure(sine), 2u);
    {
        detail::RewriteScope outer(&budget);
        detail::RewriteScope inner(&budget);
        EXPECT_THROW(budget.measure(sine), CasError);
    }
    EXPECT_EQ(context.recursion_depth(), 0u);
    EXPECT_EQ(budget.measure(sine), 2u);
}

TEST(LmcasNormalization, CandidateMeasurementReleasesNodes) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 8, "test.measure");
    std::weak_ptr<const SymbolicNode> discarded;
    {
        auto candidate = detail::make_node<VariableNode>("discarded");
        discarded = candidate;
        EXPECT_EQ(budget.measure(candidate), 1u);
    }
    EXPECT_TRUE(discarded.expired());
    std::vector<detail::SymbolicNodePtr> retained;
    for (std::size_t i = 0; i < 40; ++i) {
        auto x = detail::make_node<VariableNode>(std::to_string(i));
        auto sine = detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin, std::vector<detail::SymbolicNodePtr>{x});
        EXPECT_EQ(budget.measure(sine), 2u);
        retained.push_back(sine);
        EXPECT_LE(detail::RewriteBudgetTestAccess::cache_size(budget), 8u);
    }
    for (const auto& candidate : retained) {
        EXPECT_EQ(budget.measure(candidate), 2u);
        EXPECT_LE(detail::RewriteBudgetTestAccess::cache_size(budget), 8u);
    }
}

TEST(LmcasNormalization, ShapeChecksObserveCancellation) {
    CancellationToken cancellation;
    ComputationContext context({}, cancellation);
    detail::RewriteBudget budget(context, 64, 32, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    EXPECT_EQ(budget.measure(x), 1u);
    cancellation.cancel();
    EXPECT_THROW(budget.measure(x), CasError);
    EXPECT_THROW(budget.require_nodes(33), CasError);
    try {
        budget.append_size(32, 1);
        FAIL() << "cancelled forest check must fail";
    } catch (const CasError& error) {
        EXPECT_EQ(error.code, CasErrc::Cancelled);
    }
}

TEST(LmcasNormalization, ExpiredCacheAddressIsMeasuredAgain) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 8, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    detail::RewriteBudgetTestAccess::expired_entry(budget, x.get());
    EXPECT_EQ(budget.measure(x), 1u);
}

TEST(LmcasNormalization, ForestMeasurementCountsEveryReference) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 64, 5, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    auto sine = detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin, std::vector<detail::SymbolicNodePtr>{x});
    std::size_t nodes = budget.append_size(0, budget.measure(sine));
    nodes = budget.append_size(nodes, budget.measure(sine));
    EXPECT_EQ(nodes, 4u);
    EXPECT_THROW(budget.append_size(nodes, budget.measure(sine)), CasError);
    EXPECT_EQ(nodes, 4u);
}

TEST(LmcasNormalization, MeasurementPreservesRewriteAllowance) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 1, 64, 8, "test.measure");
    auto x = detail::make_node<VariableNode>("x");
    EXPECT_EQ(budget.measure(x), 1u);
    EXPECT_EQ(budget.measure(x), 1u);
    EXPECT_GT(context.steps_used(), 0u);
    budget.consume();
    EXPECT_THROW(budget.measure(x), CasError);
}

TEST(LmcasNormalization, MeasurementUsesRuntimeDepth) {
    ComputationContext context;
    detail::RewriteBudget budget(context, 256, 256, "test.measure");
    detail::SymbolicNodePtr node = detail::make_node<VariableNode>("x");
    for (std::size_t i = 0; i < 220; ++i) {
        node = detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin, std::vector<detail::SymbolicNodePtr>{node});
    }
    EXPECT_EQ(budget.measure(node), 221u);
    EXPECT_EQ(budget.measure(node), 221u);
    EXPECT_EQ(context.recursion_depth(), 0u);
}
