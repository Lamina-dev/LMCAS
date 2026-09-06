#include <gtest/gtest.h>

#include "computation_context.hpp"

#include <limits>
#include <thread>
#include <type_traits>

namespace LMCAS {
namespace {

TEST(ComputationContextTest, HasExclusiveValueSemantics) {
    EXPECT_FALSE(std::is_copy_constructible_v<ComputationContext>);
    EXPECT_FALSE(std::is_copy_assignable_v<ComputationContext>);
    EXPECT_FALSE(std::is_move_constructible_v<ComputationContext>);
}

TEST(ComputationContextTest, CountersRejectOverflow) {
    ResourceLimits limits;
    limits.max_steps = std::numeric_limits<std::size_t>::max();
    limits.max_ast_nodes = std::numeric_limits<std::size_t>::max();
    ComputationContext context(limits);

    ASSERT_TRUE(context.consume_steps(limits.max_steps, "test.steps"));
    auto step_overflow = context.consume_steps(1, "test.steps");
    ASSERT_FALSE(step_overflow);
    EXPECT_EQ(step_overflow.error().code, CasErrc::ResourceLimit);

    ASSERT_TRUE(context.reserve_nodes(limits.max_ast_nodes, "test.nodes"));
    auto node_overflow = context.reserve_nodes(1, "test.nodes");
    ASSERT_FALSE(node_overflow);
    EXPECT_EQ(node_overflow.error().code, CasErrc::ResourceLimit);
}

TEST(ComputationContextTest, CancellationAppliesToBudgetOperations) {
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext context({}, cancellation);

    auto nodes = context.reserve_nodes(0, "test.cancel");
    ASSERT_FALSE(nodes);
    EXPECT_EQ(nodes.error().code, CasErrc::Cancelled);

    auto diagnostic = context.add_diagnostic(
        {DiagnosticSeverity::Info, "test.cancel", "unused"});
    ASSERT_FALSE(diagnostic);
    EXPECT_EQ(diagnostic.error().code, CasErrc::Cancelled);

    auto assumptions = context.set_assumptions(nullptr, "test.cancel");
    ASSERT_FALSE(assumptions);
    EXPECT_EQ(assumptions.error().code, CasErrc::Cancelled);
}

TEST(ComputationContextTest, DiagnosticEngineEnforcesBudgetAndDispatches) {
    ResourceLimits limits;
    limits.max_diagnostics = 1;
    ComputationContext context(limits);
    std::size_t consumed = 0;
    context.set_diagnostic_consumer([&](const Diagnostic&) { ++consumed; });

    EXPECT_TRUE(context.add_diagnostic(
        {DiagnosticSeverity::Warning, "test.diagnostics", "first"}));
    auto overflow = context.add_diagnostic(
        {DiagnosticSeverity::Warning, "test.diagnostics", "second"});
    ASSERT_FALSE(overflow);
    EXPECT_EQ(overflow.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(context.diagnostics().size(), 1U);
    EXPECT_EQ(consumed, 1U);
}

TEST(ComputationContextTest, RejectsCrossThreadUse) {
    ComputationContext context;
    CasErrc error = CasErrc::Cancelled;
    std::thread worker([&] {
        auto result = context.consume_steps(1, "test.thread");
        ASSERT_FALSE(result);
        error = result.error().code;
    });
    worker.join();

    EXPECT_EQ(error, CasErrc::InternalInvariant);
    EXPECT_EQ(context.steps_used(), 0U);
}

} // namespace
} // namespace LMCAS
