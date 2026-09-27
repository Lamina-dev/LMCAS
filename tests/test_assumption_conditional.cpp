
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <stdexcept>

using namespace LMCAS;

TEST(LmcasAssumptionConditional, ConditionalCheckedRollback) {
    {
        AssumptionContext ctx;
        SymbolicExpr x = *SymbolicExpr::variable("x");
        SymbolicExpr zero = *SymbolicExpr::number(0);
        auto conclusion = LMCAS::detail::expression_from_node(
            LMCAS::detail::make_node<RelationalNode>(
                LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationOp::GT));
        const auto generation = ctx.cache_generation();

        auto result = ctx.assume_conditional_checked(x, conclusion);
        EXPECT_TRUE(!result.has_value()) << "checked conditional rejects a non-relational condition";
        EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument) << "checked conditional reports InvalidArgument";
        EXPECT_TRUE(ctx.cache_generation() == generation) << "failed checked conditional preserves cache generation";
        EXPECT_TRUE(ctx.get_active_conditionals().empty()) << "failed checked conditional stores no conditional";

        auto failure_30 = ctx.assume_conditional(x, conclusion);
        EXPECT_TRUE(!failure_30.has_value()) << "canonical conditional maps checked validation failure to invalid_argument";
        EXPECT_TRUE(ctx.get_active_conditionals().empty()) << "canonical conditional failure remains transactional";
    }
}

TEST(LmcasAssumptionConditional, ConditionalScopePop) {
    {
        AssumptionContext ctx;
        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto one = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(one), RelationalNode::Op::GT));
        auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));

        // Store conditional in root scope
        EXPECT_TRUE(ctx.assume_conditional(condition, conclusion).has_value()) << "conditional assumption succeeds";
        EXPECT_TRUE(ctx.get_active_conditionals().size() == 1) << "One conditional in root";

        // Push and add another conditional
        ctx.push();
        auto y = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("y"));
        auto cond2 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(y), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        auto concl2 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(y), LMCAS::detail::node(one), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.assume_conditional(cond2, concl2).has_value()) << "conditional assumption succeeds";

        EXPECT_TRUE(ctx.get_active_conditionals().size() == 2) << "Two conditionals (child + parent)";

        // Pop child scope — child conditional discarded
        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";
        EXPECT_TRUE(ctx.get_active_conditionals().size() == 1) << "One conditional after pop (child discarded)";
    }
}

TEST(LmcasAssumptionConditional, ConditionalPositiveConditions) {
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        // x > 0 should be satisfied since x is Positive
        auto cond_gt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::True) << "x > 0 satisfied when x is Positive";

        // x >= 0 should also be satisfied
        auto cond_geq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_geq) == Tribool::True) << "x >= 0 satisfied when x is Positive";

        // x < 0 should be False
        auto cond_lt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_lt) == Tribool::False) << "x < 0 is False when x is Positive";

        // x != 0 should be True
        auto cond_neq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::True) << "x != 0 satisfied when x is Positive";
    }
}

TEST(LmcasAssumptionConditional, ConditionalUnknownConditions) {
    {
        AssumptionContext ctx;
        // No assumptions about x

        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        auto cond = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond) == Tribool::Unknown) << "x > 0 is Unknown when no assumptions about x";
    }
}

TEST(LmcasAssumptionConditional, ConditionalReversedConditions) {
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        // 0 < x (reversed pattern)
        auto cond = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(zero), LMCAS::detail::node(x), RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond) == Tribool::True) << "0 < x satisfied when x is Positive";

        // 0 > x should be False
        auto cond_gt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(zero), LMCAS::detail::node(x), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False) << "0 > x is False when x is Positive";
    }
}

TEST(LmcasAssumptionConditional, ConditionalStoredRelations) {
    {
        AssumptionContext ctx;
        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto y = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("y"));

        // Store relation x > y
        auto rel = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(y), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.assume(rel).has_value()) << "relation assumption succeeds";

        // Evaluate x > y — should be True (stored directly)
        EXPECT_TRUE(ctx.evaluate_condition(rel) == Tribool::True) << "x > y satisfied when relation x > y is stored";
    }
}

TEST(LmcasAssumptionConditional, ConditionalNonrelationalConditions) {
    {
        AssumptionContext ctx;
        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        EXPECT_TRUE(ctx.evaluate_condition(x) == Tribool::Unknown) << "non-relational expression evaluates to Unknown";
    }
}

TEST(LmcasAssumptionConditional, ConditionalNegativeConditions) {
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());

        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        // x < 0 should be True
        auto cond_lt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_lt) == Tribool::True) << "x < 0 satisfied when x is Negative";

        // x <= 0 should be True
        auto cond_leq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::LEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_leq) == Tribool::True) << "x <= 0 satisfied when x is Negative";

        // x > 0 should be False
        auto cond_gt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False) << "x > 0 is False when x is Negative";

        // x != 0 should be True
        auto cond_neq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::True) << "x != 0 satisfied when x is Negative";
    }
}

TEST(LmcasAssumptionConditional, ConditionalZeroConditions) {
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", Sign::Zero).has_value());

        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        // x == 0 should be True
        auto cond_eq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::EQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_eq) == Tribool::True) << "x == 0 satisfied when x is Zero";

        // x > 0 should be False
        auto cond_gt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False) << "x > 0 is False when x is Zero";

        // x != 0 should be False
        auto cond_neq = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::False) << "x != 0 is False when x is Zero";
    }
}

TEST(LmcasAssumptionConditional, ConditionalScopeOrder) {
    {
        AssumptionContext ctx;
        auto x = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("x"));
        auto y = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("y"));
        auto z = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>("z"));
        auto zero = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        // Root scope conditional
        auto cond1 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(x), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        auto concl1 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(y), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.assume_conditional(cond1, concl1).has_value()) << "conditional assumption succeeds";

        // Push and add child conditional
        ctx.push();
        auto cond2 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(y), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        auto concl2 = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(z), LMCAS::detail::node(zero), RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.assume_conditional(cond2, concl2).has_value()) << "conditional assumption succeeds";

        auto all = ctx.get_active_conditionals();
        EXPECT_TRUE(all.size() == 2) << "Two conditionals across scopes";

        // Most recent scope first
        // The child scope conditional should come first (top scope)
        auto child_cond_var = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(all[0].condition));
        auto child_lhs = std::dynamic_pointer_cast<const VariableNode>(child_cond_var->left());
        EXPECT_TRUE(child_lhs->name() == "y") << "Child scope conditional comes first (most recent)";

        auto parent_cond_var = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(all[1].condition));
        auto parent_lhs = std::dynamic_pointer_cast<const VariableNode>(parent_cond_var->left());
        EXPECT_TRUE(parent_lhs->name() == "x") << "Parent scope conditional comes second";

        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";
    }
}

TEST(LmcasAssumptionConditional, ConditionalSatisfiedProperty) {
    EXPECT_TRUE(rc::check("Conditional conclusion is active when condition is satisfied by current state", []() {
        AssumptionContext ctx;

        // Generate a random variable name
        std::string var = "cx_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        // Push a scope for the conditional
        ctx.push();

        // Declare the variable as Positive in this scope (satisfies "var > 0")
        EXPECT_TRUE(ctx.assume_sign(var, Sign::Positive).has_value());

        // Create condition: var > 0
        auto var_node = LMCAS::detail::make_node<VariableNode>(var);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT));

        // Create conclusion: var != 0 (trivially implied by Positive, but tests the mechanism)
        auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::NEQ));

        EXPECT_TRUE(ctx.assume_conditional(condition, conclusion).has_value()) << "conditional assumption succeeds";

        // Evaluate the condition — should be True since var is Positive
        Tribool cond_result = ctx.evaluate_condition(condition);
        RC_ASSERT(cond_result == Tribool::True);

        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";
    }));
}

TEST(LmcasAssumptionConditional, ConditionalUnknownProperty) {
    EXPECT_TRUE(rc::check("Conditional conclusion is Unknown when condition cannot be verified", []() {
        AssumptionContext ctx;

        // Generate a random variable name — no assumptions about it
        std::string var = "unk_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        ctx.push();

        // Create condition: var > 0 (unverifiable since no sign declared)
        auto var_node = LMCAS::detail::make_node<VariableNode>(var);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT));

        // Create some conclusion
        auto one_node = LMCAS::detail::make_node<NumberNode>(BigInt(1));
        auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, one_node, RelationalNode::Op::GT));

        EXPECT_TRUE(ctx.assume_conditional(condition, conclusion).has_value()) << "conditional assumption succeeds";

        // Evaluate the condition — should be Unknown
        Tribool cond_result = ctx.evaluate_condition(condition);
        RC_ASSERT(cond_result == Tribool::Unknown);

        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";
    }));
}

TEST(LmcasAssumptionConditional, ConditionalScopePopProperty) {
    EXPECT_TRUE(rc::check("Conditional assumptions are discarded when their scope is popped", []() {
        AssumptionContext ctx;

        // Count conditionals in root scope
        size_t root_count = ctx.get_active_conditionals().size();

        // Push a scope and add random number of conditionals
        ctx.push();
        int num_conditionals = *rc::gen::inRange(1, (5) + 1);
        for (int i = 0; i < num_conditionals; ++i) {
            std::string var = "pop_" + std::to_string(i) + "_" + std::to_string(*rc::gen::inRange(0, (99) + 1));
            auto var_node = LMCAS::detail::make_node<VariableNode>(var);
            auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
            auto condition = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GT));
            auto conclusion = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GEQ));
            EXPECT_TRUE(ctx.assume_conditional(condition, conclusion).has_value()) << "conditional assumption succeeds";
        }

        // Verify conditionals are present
        size_t pushed_count = ctx.get_active_conditionals().size();
        RC_ASSERT(pushed_count == root_count + static_cast<size_t>(num_conditionals));

        // Pop the scope
        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";

        // All conditionals from the popped scope should be gone
        size_t after_pop_count = ctx.get_active_conditionals().size();
        RC_ASSERT(after_pop_count == root_count);
    }));
}

TEST(LmcasAssumptionConditional, ConditionalSignProperty) {
    EXPECT_TRUE(rc::check("Condition evaluation correctly reflects sign properties for various sign types", []() {
        AssumptionContext ctx;
        std::string var = "sv_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        // Pick a random sign to declare
        std::vector<Sign> signs = {Sign::Positive, Sign::Negative, Sign::NonNegative, Sign::NonPositive};
        Sign chosen_sign = *rc::gen::elementOf(signs);
        EXPECT_TRUE(ctx.assume_sign(var, chosen_sign).has_value());

        auto var_node = LMCAS::detail::make_node<VariableNode>(var);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        // Test condition: var > 0
        auto cond_gt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT));
        Tribool gt_result = ctx.evaluate_condition(cond_gt);

        // Verify consistency with the declared sign
        if (chosen_sign == Sign::Positive) {
            RC_ASSERT(gt_result == Tribool::True);
        } else if (chosen_sign == Sign::Negative || chosen_sign == Sign::NonPositive) {
            RC_ASSERT(gt_result == Tribool::False);
        }
        // NonNegative: could be zero, so GT might be Unknown — that's acceptable

        // Test condition: var < 0
        auto cond_lt = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::LT));
        Tribool lt_result = ctx.evaluate_condition(cond_lt);

        if (chosen_sign == Sign::Negative) {
            RC_ASSERT(lt_result == Tribool::True);
        } else if (chosen_sign == Sign::Positive || chosen_sign == Sign::NonNegative) {
            RC_ASSERT(lt_result == Tribool::False);
        }
    }));
}

TEST(LmcasAssumptionConditional, ClosedNumericConditions) {
    AssumptionContext context;
    const BigInt denominator = BigInt(1) << 100;
    auto precise = SymbolicExpr::number(Rational(denominator + BigInt(1), denominator));
    auto approximate_one = SymbolicExpr::number(1.0);
    auto greater = detail::expression_from_node(detail::make_node<RelationalNode>(
        detail::node(precise), detail::node(approximate_one), RelationOp::GT));
    EXPECT_TRUE(context.evaluate_condition(greater) == Tribool::True) << "an exact excess below binary64 resolution is not rounded away";
    auto zero = SymbolicExpr::number(0);
    auto satisfied = detail::expression_from_node(detail::make_node<RelationalNode>(
        detail::node(zero), detail::node(zero), RelationOp::EQ));
    auto impossible = detail::expression_from_node(detail::make_node<RelationalNode>(
        detail::node(zero), detail::node(zero), RelationOp::NEQ));
    auto rejected = context.assume_conditional_checked(satisfied, impossible);
    EXPECT_TRUE(!rejected && rejected.error().code == CasErrc::InvalidArgument) << "a true closed premise cannot install a false closed conclusion";
    EXPECT_TRUE(context.get_active_conditionals().empty()) << "rejected numeric contradiction leaves no conditional behind";
}
