#include "computation_context.hpp"
#include "test_common.hpp"

#include <limits>
#include <thread>
#include <type_traits>

int main() {
    using namespace lamina;

    TEST_CASE("ComputationContext has exclusive value semantics");
    EXPECT_TRUE(!std::is_copy_constructible_v<ComputationContext>,
                "contexts cannot be copied");
    EXPECT_TRUE(!std::is_copy_assignable_v<ComputationContext>,
                "contexts cannot be copy-assigned");
    EXPECT_TRUE(!std::is_move_constructible_v<ComputationContext>,
                "contexts cannot be moved between owners");

    TEST_CASE("Step and node counters reject overflow");
    ResourceLimits maximum_limits;
    maximum_limits.max_steps = std::numeric_limits<std::size_t>::max();
    maximum_limits.max_ast_nodes = std::numeric_limits<std::size_t>::max();
    ComputationContext maximum_context(maximum_limits);
    auto all_steps = maximum_context.consume_steps(maximum_limits.max_steps, "test_steps");
    EXPECT_TRUE(all_steps.has_value(), "maximum step reservation is representable");
    auto step_overflow = maximum_context.consume_steps(1, "test_steps");
    EXPECT_TRUE(!step_overflow && step_overflow.error().code == CasErrc::ResourceLimit,
                "step counter cannot wrap past size_t maximum");
    auto all_nodes = maximum_context.reserve_nodes(maximum_limits.max_ast_nodes, "test_nodes");
    EXPECT_TRUE(all_nodes.has_value(), "maximum node reservation is representable");
    auto node_overflow = maximum_context.reserve_nodes(1, "test_nodes");
    EXPECT_TRUE(!node_overflow && node_overflow.error().code == CasErrc::ResourceLimit,
                "node counter cannot wrap past size_t maximum");

    TEST_CASE("Cancellation applies to every budget operation");
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled_nodes = cancelled_context.reserve_nodes(0, "test_cancel");
    EXPECT_TRUE(!cancelled_nodes && cancelled_nodes.error().code == CasErrc::Cancelled,
                "node reservation observes cancellation");
    auto cancelled_diagnostic = cancelled_context.add_diagnostic(
        {DiagnosticSeverity::Info, "test_cancel", "unused"});
    EXPECT_TRUE(!cancelled_diagnostic &&
                    cancelled_diagnostic.error().code == CasErrc::Cancelled,
                "diagnostic collection observes cancellation");
    auto cancelled_assumptions = cancelled_context.set_assumptions(
        nullptr, "test_cancel");
    EXPECT_TRUE(!cancelled_assumptions &&
                    cancelled_assumptions.error().code == CasErrc::Cancelled,
                "assumption snapshot updates observe cancellation");

    TEST_CASE("Diagnostic collection has an explicit budget");
    ResourceLimits diagnostic_limits;
    diagnostic_limits.max_diagnostics = 1;
    ComputationContext diagnostic_context(diagnostic_limits);
    auto first = diagnostic_context.add_diagnostic(
        {DiagnosticSeverity::Warning, "test_diagnostics", "first"});
    auto second = diagnostic_context.add_diagnostic(
        {DiagnosticSeverity::Warning, "test_diagnostics", "second"});
    EXPECT_TRUE(first.has_value(), "first diagnostic fits the budget");
    EXPECT_TRUE(!second && second.error().code == CasErrc::ResourceLimit,
                "diagnostic overflow returns ResourceLimit");
    EXPECT_TRUE(diagnostic_context.diagnostics().size() == 1,
                "failed diagnostic insertion does not mutate collection");

    TEST_CASE("DiagnosticEngine is the only diagnostic dispatch point");
    std::size_t consumed = 0;
    diagnostic_context.set_diagnostic_consumer(
        [&](const Diagnostic& diagnostic) {
            if (diagnostic.operation == "test_consumer") ++consumed;
        });
    ResourceLimits consumer_limits;
    consumer_limits.max_diagnostics = 2;
    ComputationContext consumer_context(consumer_limits);
    consumer_context.set_diagnostic_consumer(
        [&](const Diagnostic&) { ++consumed; });
    auto emitted = consumer_context.add_diagnostic(
        {DiagnosticSeverity::Info, "test_consumer", "event"});
    EXPECT_TRUE(emitted.has_value() && consumed == 1,
                "context diagnostics are dispatched by DiagnosticEngine");

    TEST_CASE("Diagnostic consumer failures preserve engine state");
    ComputationContext failing_consumer_context;
    failing_consumer_context.set_diagnostic_consumer(
        [](const Diagnostic&) { throw std::runtime_error("consumer failure"); });
    auto consumer_failure = failing_consumer_context.add_diagnostic(
        {DiagnosticSeverity::Error, "test_consumer_failure", "event"});
    EXPECT_TRUE(!consumer_failure &&
                    consumer_failure.error().code == CasErrc::InternalInvariant,
                "consumer exceptions become CasError values");
    EXPECT_TRUE(failing_consumer_context.diagnostics().empty(),
                "failed dispatch does not retain a partial diagnostic");

    TEST_CASE("Cross-thread context use is rejected");
    ComputationContext thread_context;
    CasErrc thread_error = CasErrc::Cancelled;
    bool failed = false;
    std::thread worker([&] {
        auto result = thread_context.consume_steps(1, "test_thread");
        failed = !result;
        if (!result) thread_error = result.error().code;
    });
    worker.join();
    EXPECT_TRUE(failed && thread_error == CasErrc::InternalInvariant,
                "a context cannot be shared with another thread");
    EXPECT_TRUE(thread_context.steps_used() == 0,
                "rejected cross-thread access does not consume budget");

    return TEST_REPORT();
}
