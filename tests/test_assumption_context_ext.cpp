/**
 * @file test_assumption_context_ext.cpp
 * @brief Unit tests for AssumptionContext extensions (task 8.7).
 *
 * Covers:
 * - Conditional assumptions: active when condition satisfied, discarded on pop
 * - with_assumptions: preserves scope depth on success and exception
 * - Serialization round-trip: empty, single-scope, multi-scope contexts
 * - Malformed deserialization: missing END, unknown keyword, DOMAIN before SCOPE
 *
 * Validates: Requirements 5.2, 5.4, 17.2, 17.3, 20.3, 20.4
 */

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <stdexcept>
#include <string>

using namespace lamina;

// ============================================================
// Helpers
// ============================================================

static Interval make_closed_interval(double lo, double hi) {
    auto lower_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lower_val);
    iv.upper = Endpoint::closed(upper_val);
    return iv;
}

// ============================================================
// 1. Conditional assumptions
// ============================================================

static void test_conditional_active_when_condition_satisfied() {
    TEST_CASE("Conditional active when condition satisfied (Req 5.2)");

    AssumptionContext ctx;

    // Declare x as Positive
    ctx.assume_sign("x", Sign::Positive);

    auto x = lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>("x"));
    auto zero = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(0)));
    auto y = lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>("y"));

    // Condition: x > 0 (which is satisfied since x is Positive)
    auto condition = lamina::detail::expression_from_node(lamina::detail::make_node<RelationalNode>(
        lamina::detail::node(x), lamina::detail::node(zero), RelationalNode::Op::GT));
    // Conclusion: y > 0
    auto conclusion = lamina::detail::expression_from_node(lamina::detail::make_node<RelationalNode>(
        lamina::detail::node(y), lamina::detail::node(zero), RelationalNode::Op::GT));

    ctx.assume_conditional(condition, conclusion);

    // The condition x > 0 should evaluate to True
    Tribool cond_result = ctx.evaluate_condition(condition);
    EXPECT_TRUE(cond_result == Tribool::True,
                "Condition x > 0 evaluates to True when x is Positive");

    // Verify the conditional is stored
    auto conditionals = ctx.get_active_conditionals();
    EXPECT_TRUE(conditionals.size() == 1, "One conditional stored");
}

static void test_conditional_discarded_on_pop() {
    TEST_CASE("Conditional discarded on scope pop (Req 5.4)");

    AssumptionContext ctx;

    auto x = lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>("x"));
    auto zero = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(0)));

    // Push a new scope and add a conditional there
    ctx.push();

    auto condition = lamina::detail::expression_from_node(lamina::detail::make_node<RelationalNode>(
        lamina::detail::node(x), lamina::detail::node(zero), RelationalNode::Op::GT));
    auto conclusion = lamina::detail::expression_from_node(lamina::detail::make_node<RelationalNode>(
        lamina::detail::node(x), lamina::detail::node(zero), RelationalNode::Op::GEQ));

    ctx.assume_conditional(condition, conclusion);
    EXPECT_TRUE(ctx.get_active_conditionals().size() == 1,
                "Conditional present in pushed scope");

    // Pop the scope — conditional should be gone
    ctx.pop();
    EXPECT_TRUE(ctx.get_active_conditionals().size() == 0,
                "Conditional discarded after pop");
}

// ============================================================
// 2. with_assumptions
// ============================================================

static void test_with_assumptions_callable_sees_assumptions() {
    TEST_CASE("with_assumptions: callable sees the assumptions (Req 17.2)");

    AssumptionContext ctx;
    int initial_depth = ctx.depth();

    bool saw_positive = with_assumptions(ctx,
        {
            AssumptionDecl::make_domain("x", Domain::Real),
            AssumptionDecl::make_sign("x", Sign::Positive)
        },
        [&]() -> bool {
            return ctx.has_sign("x", Sign::Positive) &&
                   ctx.has_domain("x", Domain::Real);
        });

    EXPECT_TRUE(saw_positive, "Callable sees domain and sign assumptions");
    EXPECT_TRUE(ctx.depth() == initial_depth,
                "Depth restored after with_assumptions");
}

static void test_with_assumptions_preserves_depth_on_success() {
    TEST_CASE("with_assumptions: depth same before and after on success (Req 17.2)");

    AssumptionContext ctx;
    int depth_before = ctx.depth();

    with_assumptions(ctx,
        { AssumptionDecl::make_sign("y", Sign::Negative) },
        [&]() {
            // Just verify we're one level deeper inside
            EXPECT_TRUE(ctx.depth() == depth_before + 1,
                        "Inside with_assumptions, depth is +1");
        });

    EXPECT_TRUE(ctx.depth() == depth_before,
                "After with_assumptions, depth is restored");
}

static void test_with_assumptions_preserves_depth_on_exception() {
    TEST_CASE("with_assumptions: depth restored on exception (Req 17.3)");

    AssumptionContext ctx;
    int depth_before = ctx.depth();

    bool caught = false;
    try {
        with_assumptions(ctx,
            { AssumptionDecl::make_sign("z", Sign::Positive) },
            [&]() {
                throw std::runtime_error("test exception");
            });
    } catch (const std::runtime_error& e) {
        caught = true;
        EXPECT_EQ_STR(std::string(e.what()), "test exception",
                      "Exception propagated correctly");
    }

    EXPECT_TRUE(caught, "Exception was caught");
    EXPECT_TRUE(ctx.depth() == depth_before,
                "Depth restored after exception in with_assumptions");
}

static void test_with_assumptions_checked_success_and_rollback() {
    TEST_CASE("with_assumptions_checked: success returns value and restores scope");

    AssumptionContext ctx;
    int depth_before = ctx.depth();

    auto result = with_assumptions_checked(ctx,
        {
            AssumptionDecl::make_domain("x", Domain::Real),
            AssumptionDecl::make_sign("x", Sign::Positive)
        },
        [&]() -> bool {
            return ctx.has_domain("x", Domain::Real) &&
                   ctx.has_sign("x", Sign::Positive);
        });

    EXPECT_TRUE(result.has_value(), "checked with_assumptions succeeds");
    if (result) {
        EXPECT_TRUE(result.value(), "checked callable sees temporary assumptions");
    }
    EXPECT_TRUE(ctx.depth() == depth_before,
                "checked with_assumptions restores depth on success");
    EXPECT_FALSE(ctx.has_sign("x", Sign::Positive),
                 "checked with_assumptions removes temporary sign after success");
}

static void test_with_assumptions_checked_decl_failure_rolls_back() {
    TEST_CASE("with_assumptions_checked: declaration failure restores scope");

    AssumptionContext ctx;
    int depth_before = ctx.depth();
    bool called = false;

    auto result = with_assumptions_checked(ctx,
        { AssumptionDecl::make_sign("", Sign::Positive) },
        [&]() -> int {
            called = true;
            return 1;
        });

    EXPECT_TRUE(!result.has_value(),
                "checked with_assumptions reports declaration failure");
    EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument,
                "checked with_assumptions preserves declaration error code");
    EXPECT_FALSE(called, "checked with_assumptions does not call body after declaration failure");
    EXPECT_TRUE(ctx.depth() == depth_before,
                "checked with_assumptions restores depth after declaration failure");
}

static void test_with_assumptions_checked_callable_exception() {
    TEST_CASE("with_assumptions_checked: callable exception becomes CasError");

    AssumptionContext ctx;
    int depth_before = ctx.depth();

    auto result = with_assumptions_checked(ctx,
        { AssumptionDecl::make_sign("z", Sign::Positive) },
        [&]() {
            throw std::runtime_error("checked body failure");
        });

    EXPECT_TRUE(!result.has_value(),
                "checked void with_assumptions reports callable exception");
    EXPECT_TRUE(result.error().code == CasErrc::InternalInvariant,
                "checked void with_assumptions maps callable exception to InternalInvariant");
    EXPECT_TRUE(ctx.depth() == depth_before,
                "checked void with_assumptions restores depth after callable exception");
    EXPECT_FALSE(ctx.has_sign("z", Sign::Positive),
                 "checked void with_assumptions removes temporary sign after callable exception");
}

static void test_checked_interval_and_definiteness_queries() {
    TEST_CASE("AssumptionContext checked interval and definiteness queries");

    AssumptionContext ctx;
    Interval unit = make_closed_interval(0.0, 1.0);
    Interval larger = make_closed_interval(0.0, 2.0);

    auto declared = ctx.current_properties().declare_continuous_checked("f", larger);
    EXPECT_TRUE(declared.has_value(), "checked continuous declaration succeeds");

    auto continuous = ctx.is_continuous_checked("f", unit);
    EXPECT_TRUE(continuous.has_value(), "checked continuity query succeeds");
    if (continuous) {
        EXPECT_TRUE(continuous.value() == Tribool::True,
                    "checked continuity query returns True when covered");
    }

    auto differentiable = ctx.is_differentiable_checked("f", unit);
    EXPECT_TRUE(differentiable.has_value(), "checked differentiability query succeeds");
    if (differentiable) {
        EXPECT_TRUE(differentiable.value() == Tribool::Unknown,
                    "checked differentiability query returns Unknown without declaration");
    }

    auto empty_symbol = ctx.is_continuous_checked("", unit);
    EXPECT_TRUE(!empty_symbol.has_value(),
                "checked continuity query rejects empty symbols");
    if (!empty_symbol) {
        EXPECT_TRUE(empty_symbol.error().code == CasErrc::InvalidArgument,
                    "checked continuity empty symbol reports InvalidArgument");
    }
    EXPECT_TRUE(ctx.is_continuous("", unit) == Tribool::Unknown,
                "legacy continuity query unwraps checked failure to Unknown");

    auto positive_def_decl =
        ctx.current_properties().declare_definiteness_checked(
            "M", Definiteness::PositiveDefinite);
    EXPECT_TRUE(positive_def_decl.has_value(),
                "checked positive-definite declaration succeeds");

    auto positive_def = ctx.is_positive_definite_checked("M");
    EXPECT_TRUE(positive_def.has_value(), "checked positive-definite query succeeds");
    if (positive_def) {
        EXPECT_TRUE(positive_def.value() == Tribool::True,
                    "checked positive-definite query returns True");
    }

    auto positive_semidef = ctx.is_positive_semidefinite_checked("M");
    EXPECT_TRUE(positive_semidef.has_value(),
                "checked positive-semidefinite query succeeds");
    if (positive_semidef) {
        EXPECT_TRUE(positive_semidef.value() == Tribool::True,
                    "positive definite implies positive semidefinite");
    }

    auto empty_matrix_symbol = ctx.is_positive_definite_checked("");
    EXPECT_TRUE(!empty_matrix_symbol.has_value(),
                "checked positive-definite query rejects empty symbols");
    if (!empty_matrix_symbol) {
        EXPECT_TRUE(empty_matrix_symbol.error().code == CasErrc::InvalidArgument,
                    "checked positive-definite empty symbol reports InvalidArgument");
    }
    EXPECT_TRUE(ctx.is_positive_definite("") == Tribool::Unknown,
                "legacy positive-definite query unwraps checked failure to Unknown");
}

// ============================================================
// 3. Serialization round-trip
// ============================================================

static void test_serialize_empty_context() {
    TEST_CASE("Serialization round-trip: empty context (Req 20.3)");

    AssumptionContext ctx;
    std::string serialized = ctx.serialize();

    // Deserialize
    AssumptionContext restored = AssumptionContext::deserialize(serialized);

    // Both should have depth 1 (root scope)
    EXPECT_TRUE(restored.depth() == 1, "Restored empty context has depth 1");
}

static void test_serialize_single_scope_with_domain_and_sign() {
    TEST_CASE("Serialization round-trip: single scope with domain+sign (Req 20.3)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    std::string serialized = ctx.serialize();
    AssumptionContext restored = AssumptionContext::deserialize(serialized);

    // Verify same queries
    EXPECT_TRUE(restored.has_domain("x", Domain::Real),
                "Restored context has x as Real");
    EXPECT_TRUE(restored.has_sign("x", Sign::Positive),
                "Restored context has x as Positive");
}

static void test_serialize_multi_scope() {
    TEST_CASE("Serialization round-trip: multi-scope (Req 20.3)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_sign("x", Sign::Positive);

    ctx.push();
    ctx.assume_domain("y", Domain::Real);
    ctx.assume_sign("y", Sign::Negative);

    std::string serialized = ctx.serialize();
    AssumptionContext restored = AssumptionContext::deserialize(serialized);

    // Verify depth
    EXPECT_TRUE(restored.depth() == 2, "Restored multi-scope has depth 2");

    // Verify properties in child scope
    EXPECT_TRUE(restored.has_domain("y", Domain::Real),
                "Restored context has y as Real in child scope");
    EXPECT_TRUE(restored.has_sign("y", Sign::Negative),
                "Restored context has y as Negative in child scope");

    // Verify properties from parent scope (read-through)
    EXPECT_TRUE(restored.has_domain("x", Domain::Integer),
                "Restored context has x as Integer from parent scope");
    EXPECT_TRUE(restored.has_sign("x", Sign::Positive),
                "Restored context has x as Positive from parent scope");
}

// ============================================================
// 4. Malformed deserialization
// ============================================================

static void test_deserialize_missing_end_throws() {
    TEST_CASE("Malformed deserialization: missing END throws with line number (Req 20.4)");

    std::string malformed = "SCOPE 0\nDOMAIN x Real\n";
    // No END terminator

    bool threw = false;
    std::string msg;
    try {
        AssumptionContext::deserialize(malformed);
    } catch (const std::invalid_argument& e) {
        threw = true;
        msg = e.what();
    }

    EXPECT_TRUE(threw, "Missing END throws std::invalid_argument");
    // Should contain line number info
    EXPECT_TRUE(msg.find("Line") != std::string::npos || msg.find("line") != std::string::npos,
                "Error message contains line number info");
}

static void test_deserialize_unknown_keyword_throws() {
    TEST_CASE("Malformed deserialization: unknown keyword throws (Req 20.4)");

    std::string malformed = "SCOPE 0\nFOOBAR x Real\nEND\n";

    bool threw = false;
    std::string msg;
    try {
        AssumptionContext::deserialize(malformed);
    } catch (const std::invalid_argument& e) {
        threw = true;
        msg = e.what();
    }

    EXPECT_TRUE(threw, "Unknown keyword throws std::invalid_argument");
    EXPECT_TRUE(msg.find("FOOBAR") != std::string::npos || msg.find("unknown") != std::string::npos,
                "Error message mentions the unknown keyword");
}

static void test_deserialize_domain_before_scope_throws() {
    TEST_CASE("Malformed deserialization: DOMAIN before SCOPE throws (Req 20.4)");

    std::string malformed = "DOMAIN x Real\nSCOPE 0\nEND\n";

    bool threw = false;
    std::string msg;
    try {
        AssumptionContext::deserialize(malformed);
    } catch (const std::invalid_argument& e) {
        threw = true;
        msg = e.what();
    }

    EXPECT_TRUE(threw, "DOMAIN before SCOPE throws std::invalid_argument");
    EXPECT_TRUE(msg.find("before SCOPE") != std::string::npos,
                "Error message mentions 'before SCOPE'");
}

// ============================================================
// main
// ============================================================

int main() {
    // Conditional assumptions
    test_conditional_active_when_condition_satisfied();
    test_conditional_discarded_on_pop();

    // with_assumptions
    test_with_assumptions_callable_sees_assumptions();
    test_with_assumptions_preserves_depth_on_success();
    test_with_assumptions_preserves_depth_on_exception();
    test_with_assumptions_checked_success_and_rollback();
    test_with_assumptions_checked_decl_failure_rolls_back();
    test_with_assumptions_checked_callable_exception();
    test_checked_interval_and_definiteness_queries();

    // Serialization round-trip
    test_serialize_empty_context();
    test_serialize_single_scope_with_domain_and_sign();
    test_serialize_multi_scope();

    // Malformed deserialization
    test_deserialize_missing_end_throws();
    test_deserialize_unknown_keyword_throws();
    test_deserialize_domain_before_scope_throws();

    return TEST_REPORT();
}
