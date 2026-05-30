/**
 * @file test_assumption_conditional.cpp
 * @brief Unit tests and property tests for conditional assumptions (Requirement 5).
 *
 * Tests:
 * - assume_conditional stores conditionals in the current scope
 * - Condition evaluation against current state (sign properties, relations)
 * - Conditionals are discarded on scope pop
 * - Null condition/conclusion throws std::invalid_argument
 * - get_active_conditionals retrieves from all scopes
 *
 * Property tests (Task 8.6):
 * - Property 5: Conditional assumption scope semantics — For any conditional
 *   assumption in a pushed scope, the conclusion SHALL be active when condition
 *   is satisfied, Unknown when unverifiable, and discarded after scope pop.
 *
 * **Validates: Requirements 5.2, 5.3, 5.4**
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <stdexcept>

using namespace lamina;

int main() {
    // =========================================================
    TEST_CASE("assume_conditional: null condition throws");
    // =========================================================
    {
        AssumptionContext ctx;
        SymbolicExpr null_expr;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));
        SymbolicExpr conclusion(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));

        bool threw = false;
        try {
            ctx.assume_conditional(null_expr, conclusion);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "assume_conditional throws on null condition");
    }

    // =========================================================
    TEST_CASE("assume_conditional: null conclusion throws");
    // =========================================================
    {
        AssumptionContext ctx;
        SymbolicExpr null_expr;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));
        SymbolicExpr condition(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));

        bool threw = false;
        try {
            ctx.assume_conditional(condition, null_expr);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw, "assume_conditional throws on null conclusion");
    }

    // =========================================================
    TEST_CASE("assume_conditional: stores conditional in current scope");
    // =========================================================
    {
        AssumptionContext ctx;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto one = SymbolicExpr(std::make_shared<NumberNode>(BigInt(1)));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // condition: x > 1
        SymbolicExpr condition(std::make_shared<RelationalNode>(
            x.root, one.root, RelationalNode::Op::GT));
        // conclusion: x > 0 (ln(x) > 0 would be more realistic but harder to construct)
        SymbolicExpr conclusion(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));

        ctx.assume_conditional(condition, conclusion);

        auto conditionals = ctx.get_active_conditionals();
        EXPECT_TRUE(conditionals.size() == 1, "One conditional stored");
        EXPECT_TRUE(conditionals[0].condition.root != nullptr, "Condition is not null");
        EXPECT_TRUE(conditionals[0].conclusion.root != nullptr, "Conclusion is not null");
    }

    // =========================================================
    TEST_CASE("Conditionals discarded on scope pop (Req 5.4)");
    // =========================================================
    {
        AssumptionContext ctx;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto one = SymbolicExpr(std::make_shared<NumberNode>(BigInt(1)));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        SymbolicExpr condition(std::make_shared<RelationalNode>(
            x.root, one.root, RelationalNode::Op::GT));
        SymbolicExpr conclusion(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));

        // Store conditional in root scope
        ctx.assume_conditional(condition, conclusion);
        EXPECT_TRUE(ctx.get_active_conditionals().size() == 1, "One conditional in root");

        // Push and add another conditional
        ctx.push();
        auto y = SymbolicExpr(std::make_shared<VariableNode>("y"));
        SymbolicExpr cond2(std::make_shared<RelationalNode>(
            y.root, zero.root, RelationalNode::Op::GT));
        SymbolicExpr concl2(std::make_shared<RelationalNode>(
            y.root, one.root, RelationalNode::Op::GT));
        ctx.assume_conditional(cond2, concl2);

        EXPECT_TRUE(ctx.get_active_conditionals().size() == 2, "Two conditionals (child + parent)");

        // Pop child scope — child conditional discarded
        ctx.pop();
        EXPECT_TRUE(ctx.get_active_conditionals().size() == 1, "One conditional after pop (child discarded)");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: condition satisfied by sign property (Req 5.2)");
    // =========================================================
    {
        AssumptionContext ctx;
        ctx.assume_sign("x", Sign::Positive);

        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // x > 0 should be satisfied since x is Positive
        SymbolicExpr cond_gt(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::True,
                    "x > 0 satisfied when x is Positive");

        // x >= 0 should also be satisfied
        SymbolicExpr cond_geq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_geq) == Tribool::True,
                    "x >= 0 satisfied when x is Positive");

        // x < 0 should be False
        SymbolicExpr cond_lt(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_lt) == Tribool::False,
                    "x < 0 is False when x is Positive");

        // x != 0 should be True
        SymbolicExpr cond_neq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::True,
                    "x != 0 satisfied when x is Positive");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: condition unverifiable returns Unknown (Req 5.3)");
    // =========================================================
    {
        AssumptionContext ctx;
        // No assumptions about x

        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        SymbolicExpr cond(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond) == Tribool::Unknown,
                    "x > 0 is Unknown when no assumptions about x");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: reversed pattern (0 < x) satisfied");
    // =========================================================
    {
        AssumptionContext ctx;
        ctx.assume_sign("x", Sign::Positive);

        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // 0 < x (reversed pattern)
        SymbolicExpr cond(std::make_shared<RelationalNode>(
            zero.root, x.root, RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond) == Tribool::True,
                    "0 < x satisfied when x is Positive");

        // 0 > x should be False
        SymbolicExpr cond_gt(std::make_shared<RelationalNode>(
            zero.root, x.root, RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False,
                    "0 > x is False when x is Positive");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: condition satisfied by stored relation");
    // =========================================================
    {
        AssumptionContext ctx;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto y = SymbolicExpr(std::make_shared<VariableNode>("y"));

        // Store relation x > y
        SymbolicExpr rel(std::make_shared<RelationalNode>(
            x.root, y.root, RelationalNode::Op::GT));
        ctx.assume(rel);

        // Evaluate x > y — should be True (stored directly)
        EXPECT_TRUE(ctx.evaluate_condition(rel) == Tribool::True,
                    "x > y satisfied when relation x > y is stored");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: null expression returns Unknown");
    // =========================================================
    {
        AssumptionContext ctx;
        SymbolicExpr null_expr;
        EXPECT_TRUE(ctx.evaluate_condition(null_expr) == Tribool::Unknown,
                    "null expression evaluates to Unknown");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: non-relational expression returns Unknown");
    // =========================================================
    {
        AssumptionContext ctx;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        EXPECT_TRUE(ctx.evaluate_condition(x) == Tribool::Unknown,
                    "non-relational expression evaluates to Unknown");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: Negative sign checks");
    // =========================================================
    {
        AssumptionContext ctx;
        ctx.assume_sign("x", Sign::Negative);

        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // x < 0 should be True
        SymbolicExpr cond_lt(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::LT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_lt) == Tribool::True,
                    "x < 0 satisfied when x is Negative");

        // x <= 0 should be True
        SymbolicExpr cond_leq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::LEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_leq) == Tribool::True,
                    "x <= 0 satisfied when x is Negative");

        // x > 0 should be False
        SymbolicExpr cond_gt(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False,
                    "x > 0 is False when x is Negative");

        // x != 0 should be True
        SymbolicExpr cond_neq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::True,
                    "x != 0 satisfied when x is Negative");
    }

    // =========================================================
    TEST_CASE("evaluate_condition: Zero sign checks");
    // =========================================================
    {
        AssumptionContext ctx;
        ctx.assume_sign("x", Sign::Zero);

        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // x == 0 should be True
        SymbolicExpr cond_eq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::EQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_eq) == Tribool::True,
                    "x == 0 satisfied when x is Zero");

        // x > 0 should be False
        SymbolicExpr cond_gt(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));
        EXPECT_TRUE(ctx.evaluate_condition(cond_gt) == Tribool::False,
                    "x > 0 is False when x is Zero");

        // x != 0 should be False
        SymbolicExpr cond_neq(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::NEQ));
        EXPECT_TRUE(ctx.evaluate_condition(cond_neq) == Tribool::False,
                    "x != 0 is False when x is Zero");
    }

    // =========================================================
    TEST_CASE("get_active_conditionals: multi-scope ordering");
    // =========================================================
    {
        AssumptionContext ctx;
        auto x = SymbolicExpr(std::make_shared<VariableNode>("x"));
        auto y = SymbolicExpr(std::make_shared<VariableNode>("y"));
        auto z = SymbolicExpr(std::make_shared<VariableNode>("z"));
        auto zero = SymbolicExpr(std::make_shared<NumberNode>(BigInt(0)));

        // Root scope conditional
        SymbolicExpr cond1(std::make_shared<RelationalNode>(
            x.root, zero.root, RelationalNode::Op::GT));
        SymbolicExpr concl1(std::make_shared<RelationalNode>(
            y.root, zero.root, RelationalNode::Op::GT));
        ctx.assume_conditional(cond1, concl1);

        // Push and add child conditional
        ctx.push();
        SymbolicExpr cond2(std::make_shared<RelationalNode>(
            y.root, zero.root, RelationalNode::Op::GT));
        SymbolicExpr concl2(std::make_shared<RelationalNode>(
            z.root, zero.root, RelationalNode::Op::GT));
        ctx.assume_conditional(cond2, concl2);

        auto all = ctx.get_active_conditionals();
        EXPECT_TRUE(all.size() == 2, "Two conditionals across scopes");

        // Most recent scope first
        // The child scope conditional should come first (top scope)
        auto child_cond_var = std::dynamic_pointer_cast<RelationalNode>(all[0].condition.root);
        auto child_lhs = std::dynamic_pointer_cast<VariableNode>(child_cond_var->left);
        EXPECT_TRUE(child_lhs->name == "y", "Child scope conditional comes first (most recent)");

        auto parent_cond_var = std::dynamic_pointer_cast<RelationalNode>(all[1].condition.root);
        auto parent_lhs = std::dynamic_pointer_cast<VariableNode>(parent_cond_var->left);
        EXPECT_TRUE(parent_lhs->name == "x", "Parent scope conditional comes second");

        ctx.pop();
    }

    // =========================================================
    // Property-Based Tests: Property 5 — Conditional assumption scope semantics
    // **Validates: Requirements 5.2, 5.3, 5.4**
    // =========================================================

    // --- Property 5a: Conclusion active when condition satisfied (Req 5.2) ---
    TEST_CASE("Feature: assumption-system-enhancements, Property 5: Conclusion active when condition satisfied");
    rc::check("Conditional conclusion is active when condition is satisfied by current state", []() {
        AssumptionContext ctx;

        // Generate a random variable name
        std::string var = "cx_" + std::to_string(rc::gen::inRange(0, 999));

        // Push a scope for the conditional
        ctx.push();

        // Declare the variable as Positive in this scope (satisfies "var > 0")
        ctx.assume_sign(var, Sign::Positive);

        // Create condition: var > 0
        auto var_node = std::make_shared<VariableNode>(var);
        auto zero_node = std::make_shared<NumberNode>(BigInt(0));
        SymbolicExpr condition(std::make_shared<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT));

        // Create conclusion: var != 0 (trivially implied by Positive, but tests the mechanism)
        SymbolicExpr conclusion(std::make_shared<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::NEQ));

        ctx.assume_conditional(condition, conclusion);

        // Evaluate the condition — should be True since var is Positive
        Tribool cond_result = ctx.evaluate_condition(condition);
        RC_ASSERT(cond_result == Tribool::True);

        // The conditional should be in the active list
        auto conditionals = ctx.get_active_conditionals();
        RC_ASSERT(conditionals.size() >= 1);

        ctx.pop();
    });

    // --- Property 5b: Conclusion Unknown when condition unverifiable (Req 5.3) ---
    TEST_CASE("Feature: assumption-system-enhancements, Property 5: Conclusion Unknown when condition unverifiable");
    rc::check("Conditional conclusion is Unknown when condition cannot be verified", []() {
        AssumptionContext ctx;

        // Generate a random variable name — no assumptions about it
        std::string var = "unk_" + std::to_string(rc::gen::inRange(0, 999));

        ctx.push();

        // Create condition: var > 0 (unverifiable since no sign declared)
        auto var_node = std::make_shared<VariableNode>(var);
        auto zero_node = std::make_shared<NumberNode>(BigInt(0));
        SymbolicExpr condition(std::make_shared<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT));

        // Create some conclusion
        auto one_node = std::make_shared<NumberNode>(BigInt(1));
        SymbolicExpr conclusion(std::make_shared<RelationalNode>(
            var_node, one_node, RelationalNode::Op::GT));

        ctx.assume_conditional(condition, conclusion);

        // Evaluate the condition — should be Unknown
        Tribool cond_result = ctx.evaluate_condition(condition);
        RC_ASSERT(cond_result == Tribool::Unknown);

        ctx.pop();
    });

    // --- Property 5c: Conditionals discarded after scope pop (Req 5.4) ---
    TEST_CASE("Feature: assumption-system-enhancements, Property 5: Conditionals discarded after scope pop");
    rc::check("Conditional assumptions are discarded when their scope is popped", []() {
        AssumptionContext ctx;

        // Count conditionals in root scope
        size_t root_count = ctx.get_active_conditionals().size();

        // Push a scope and add random number of conditionals
        ctx.push();
        int num_conditionals = rc::gen::inRange(1, 5);
        for (int i = 0; i < num_conditionals; ++i) {
            std::string var = "pop_" + std::to_string(i) + "_" + std::to_string(rc::gen::inRange(0, 99));
            auto var_node = std::make_shared<VariableNode>(var);
            auto zero_node = std::make_shared<NumberNode>(BigInt(0));
            SymbolicExpr condition(std::make_shared<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GT));
            SymbolicExpr conclusion(std::make_shared<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GEQ));
            ctx.assume_conditional(condition, conclusion);
        }

        // Verify conditionals are present
        size_t pushed_count = ctx.get_active_conditionals().size();
        RC_ASSERT(pushed_count == root_count + static_cast<size_t>(num_conditionals));

        // Pop the scope
        ctx.pop();

        // All conditionals from the popped scope should be gone
        size_t after_pop_count = ctx.get_active_conditionals().size();
        RC_ASSERT(after_pop_count == root_count);
    });

    // --- Property 5d: Condition satisfied by various sign types ---
    TEST_CASE("Feature: assumption-system-enhancements, Property 5: Condition evaluation with various signs");
    rc::check("Condition evaluation correctly reflects sign properties for various sign types", []() {
        AssumptionContext ctx;
        std::string var = "sv_" + std::to_string(rc::gen::inRange(0, 999));

        // Pick a random sign to declare
        std::vector<Sign> signs = {Sign::Positive, Sign::Negative, Sign::NonNegative, Sign::NonPositive};
        Sign chosen_sign = rc::gen::elementOf(signs);
        ctx.assume_sign(var, chosen_sign);

        auto var_node = std::make_shared<VariableNode>(var);
        auto zero_node = std::make_shared<NumberNode>(BigInt(0));

        // Test condition: var > 0
        SymbolicExpr cond_gt(std::make_shared<RelationalNode>(
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
        SymbolicExpr cond_lt(std::make_shared<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::LT));
        Tribool lt_result = ctx.evaluate_condition(cond_lt);

        if (chosen_sign == Sign::Negative) {
            RC_ASSERT(lt_result == Tribool::True);
        } else if (chosen_sign == Sign::Positive || chosen_sign == Sign::NonNegative) {
            RC_ASSERT(lt_result == Tribool::False);
        }
    });

    return TEST_REPORT();
}
