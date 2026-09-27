
#include "test_interval_support.hpp"
#include "assumption_context.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <string>

using namespace LMCAS;

static AssumptionContext deserialize_success(const std::string &data) {
    auto result = AssumptionContext::deserialize(data);
    EXPECT_TRUE((result.has_value())) << "assumption context deserialization succeeds";
    return result ? std::move(result.value()) : AssumptionContext();
}

TEST(AssumptionContextExt, ConditionalActiveWhenConditionSatisfied) {
    AssumptionContext ctx;

    // Declare x as Positive
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
    auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto y = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("y"));

    // Condition: x > 0 (which is satisfied since x is Positive)
    auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
        LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
    // Conclusion: y > 0
    auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
        LMCAS::detail::node(y), LMCAS::detail::node(zero), RelationalNode::Op::GT));

    EXPECT_TRUE((ctx.assume_conditional(condition, conclusion).has_value())) << "conditional assumption succeeds";

    // The condition x > 0 should evaluate to True
    Tribool cond_result = ctx.evaluate_condition(condition);
    EXPECT_TRUE((cond_result == Tribool::True)) << "Condition x > 0 evaluates to True when x is Positive";

    // Verify the conditional is stored
    auto conditionals = ctx.get_active_conditionals();
    EXPECT_TRUE((conditionals.size() == 1)) << "One conditional stored";
}

TEST(AssumptionContextExt, ConditionalDiscardedOnPop) {
    AssumptionContext ctx;

    auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
    auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

    // Push a new scope and add a conditional there
    ctx.push();

    auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
        LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
    auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
        LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GEQ));

    EXPECT_TRUE((ctx.assume_conditional(condition, conclusion).has_value())) << "conditional assumption succeeds";
    EXPECT_TRUE((ctx.get_active_conditionals().size() == 1)) << "Conditional present in pushed scope";

    // Pop the scope — conditional should be gone
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.get_active_conditionals().size() == 0)) << "Conditional discarded after pop";
}

TEST(AssumptionContextExt, WithAssumptionsCheckedDeclFailureRollsBack) {
    AssumptionContext ctx;
    int depth_before = ctx.depth();
    bool called = false;

    auto result = with_assumptions(ctx,
                                   {AssumptionDecl::make_sign("", Sign::Positive)},
                                   [&]() -> int {
                                       called = true;
                                       return 1;
                                   });

    EXPECT_TRUE((!result.has_value())) << "checked with_assumptions reports declaration failure";
    EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "checked with_assumptions preserves declaration error code";
    EXPECT_FALSE((called)) << "checked with_assumptions does not call body after declaration failure";
    EXPECT_TRUE((ctx.depth() == depth_before)) << "checked with_assumptions restores depth after declaration failure";
}

TEST(AssumptionContextExt, CheckedIntervalAndDefinitenessQueries) {
    AssumptionContext ctx;
    Interval unit = make_closed_interval(0.0, 1.0);
    Interval larger = make_closed_interval(0.0, 2.0);

    auto declared = ctx.current_properties().declare_continuous_checked("f", larger);
    EXPECT_TRUE((declared.has_value())) << "checked continuous declaration succeeds";

    auto continuous = ctx.is_continuous_checked("f", unit);
    ASSERT_TRUE((continuous.has_value())) << "checked continuity query succeeds";
    if (continuous) {
        EXPECT_TRUE((continuous.value() == Tribool::True)) << "checked continuity query returns True when covered";
    }

    auto differentiable = ctx.is_differentiable_checked("f", unit);
    ASSERT_TRUE((differentiable.has_value())) << "checked differentiability query succeeds";
    if (differentiable) {
        EXPECT_TRUE((differentiable.value() == Tribool::Unknown)) << "checked differentiability query returns Unknown without declaration";
    }

    auto empty_symbol = ctx.is_continuous_checked("", unit);
    EXPECT_TRUE((!empty_symbol.has_value())) << "checked continuity query rejects empty symbols";
    if (!empty_symbol) {
        EXPECT_TRUE((empty_symbol.error().code == CasErrc::InvalidArgument)) << "checked continuity empty symbol reports InvalidArgument";
    }
    auto canonical_empty_symbol = ctx.is_continuous("", unit);
    EXPECT_TRUE((!canonical_empty_symbol.has_value())) << "canonical continuity query preserves checked failure";

    auto positive_def_decl =
        ctx.current_properties().declare_definiteness_checked(
            "M", Definiteness::PositiveDefinite);
    EXPECT_TRUE((positive_def_decl.has_value())) << "checked positive-definite declaration succeeds";

    auto positive_def = ctx.is_positive_definite_checked("M");
    ASSERT_TRUE((positive_def.has_value())) << "checked positive-definite query succeeds";
    if (positive_def) {
        EXPECT_TRUE((positive_def.value() == Tribool::True)) << "checked positive-definite query returns True";
    }

    auto positive_semidef = ctx.is_positive_semidefinite_checked("M");
    ASSERT_TRUE((positive_semidef.has_value())) << "checked positive-semidefinite query succeeds";
    if (positive_semidef) {
        EXPECT_TRUE((positive_semidef.value() == Tribool::True)) << "positive definite implies positive semidefinite";
    }

    auto empty_matrix_symbol = ctx.is_positive_definite_checked("");
    EXPECT_TRUE((!empty_matrix_symbol.has_value())) << "checked positive-definite query rejects empty symbols";
    if (!empty_matrix_symbol) {
        EXPECT_TRUE((empty_matrix_symbol.error().code == CasErrc::InvalidArgument)) << "checked positive-definite empty symbol reports InvalidArgument";
    }
    auto canonical_empty_matrix = ctx.is_positive_definite("");
    EXPECT_TRUE((!canonical_empty_matrix.has_value())) << "canonical positive-definite query preserves checked failure";
}

TEST(AssumptionContextExt, SerializeEmptyContext) {
    AssumptionContext ctx;
    std::string serialized = ctx.serialize();

    // Deserialize
    AssumptionContext restored = deserialize_success(serialized);

    // Both should have depth 1 (root scope)
    EXPECT_TRUE((restored.depth() == 1)) << "Restored empty context has depth 1";
}

TEST(AssumptionContextExt, SerializeSingleScopeWithDomainAndSign) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    std::string serialized = ctx.serialize();
    AssumptionContext restored = deserialize_success(serialized);

    // Verify same queries
    EXPECT_TRUE((restored.has_domain("x", Domain::Real))) << "Restored context has x as Real";
    EXPECT_TRUE((restored.has_sign("x", Sign::Positive))) << "Restored context has x as Positive";
}

TEST(AssumptionContextExt, SerializeMultiScope) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());

    std::string serialized = ctx.serialize();
    AssumptionContext restored = deserialize_success(serialized);

    // Verify depth
    EXPECT_TRUE((restored.depth() == 2)) << "Restored multi-scope has depth 2";

    // Verify properties in child scope
    EXPECT_TRUE((restored.has_domain("y", Domain::Real))) << "Restored context has y as Real in child scope";
    EXPECT_TRUE((restored.has_sign("y", Sign::Negative))) << "Restored context has y as Negative in child scope";

    // Verify properties from parent scope (read-through)
    EXPECT_TRUE((restored.has_domain("x", Domain::Integer))) << "Restored context has x as Integer from parent scope";
    EXPECT_TRUE((restored.has_sign("x", Sign::Positive))) << "Restored context has x as Positive from parent scope";
}

TEST(AssumptionContextExt, DeserializeMissingEndThrows) {
    auto result = AssumptionContext::deserialize("SCOPE 0\nDOMAIN x Real\n");
    EXPECT_TRUE((!result.has_value())) << "Missing END returns failure";
    EXPECT_TRUE((result.error().code == CasErrc::ParseError)) << "Missing END reports ParseError";
    EXPECT_TRUE((result.error().message.find("Line") != std::string::npos ||
                 result.error().message.find("line") != std::string::npos))
        << "Error message contains line number info";
}

TEST(AssumptionContextExt, DeserializeUnknownKeywordThrows) {
    auto result = AssumptionContext::deserialize(
        "SCOPE 0\nFOOBAR x Real\nEND\n");
    EXPECT_TRUE((!result.has_value())) << "Unknown keyword returns failure";
    EXPECT_TRUE((result.error().code == CasErrc::ParseError)) << "Unknown keyword reports ParseError";
    EXPECT_TRUE((result.error().message.find("FOOBAR") != std::string::npos ||
                 result.error().message.find("unknown") != std::string::npos))
        << "Error message mentions the unknown keyword";
}

TEST(AssumptionContextExt, DeserializeDomainBeforeScopeThrows) {
    auto result = AssumptionContext::deserialize(
        "DOMAIN x Real\nSCOPE 0\nEND\n");
    EXPECT_TRUE((!result.has_value())) << "DOMAIN before SCOPE returns failure";
    EXPECT_TRUE((result.error().code == CasErrc::ParseError)) << "DOMAIN before SCOPE reports ParseError";
    EXPECT_TRUE((result.error().message.find("before SCOPE") != std::string::npos)) << "Error message mentions 'before SCOPE'";
}
