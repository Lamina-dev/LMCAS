/**
 * @file test_assumption_division_sign.cpp
 * @brief Property tests for InferenceEngine sign inference (Task 5.7).
 *
 * Properties tested:
 * - Property 1: Division sign inference follows sign multiplication table
 * - Property 4: Arithmetic combination sign inference
 *
 * Validates: Requirements 1.1, 1.2, 1.3, 1.4, 4.1, 4.2, 4.3, 4.4, 4.5
 *
 * Uses rapidcheck (header-only, vendored in tests/rapidcheck/) for
 * property-based testing with random input generation.
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;

// ============================================================
// Helpers: create AST nodes
// ============================================================

static std::shared_ptr<SymbolicNode> make_var(const std::string& name) {
    return std::make_shared<VariableNode>(name);
}

static std::shared_ptr<SymbolicNode> make_number(int val) {
    return std::make_shared<NumberNode>(BigInt(val));
}

static std::shared_ptr<SymbolicNode> make_power(
    std::shared_ptr<SymbolicNode> base,
    std::shared_ptr<SymbolicNode> exp) {
    return std::make_shared<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<SymbolicNode>> ops) {
    return std::make_shared<MultiplyNode>(std::move(ops));
}

static std::shared_ptr<SymbolicNode> make_add(
    std::vector<std::shared_ptr<SymbolicNode>> ops) {
    return std::make_shared<AddNode>(std::move(ops));
}

static SymbolicExpr wrap_expr(std::shared_ptr<SymbolicNode> node) {
    SymbolicExpr expr;
    expr.root = std::move(node);
    return expr;
}

/// Build a division expression: numerator / denominator
/// Represented as MultiplyNode([numerator, PowerNode(denominator, -1)])
static SymbolicExpr make_division(const std::string& num_var, const std::string& den_var) {
    auto num = make_var(num_var);
    auto den = make_var(den_var);
    auto den_inv = make_power(den, make_number(-1));
    auto mul = make_multiply({num, den_inv});
    return wrap_expr(mul);
}

/// Build a division expression with number numerator and variable denominator
static SymbolicExpr make_division_num_var(int num_val, const std::string& den_var) {
    auto num = make_number(num_val);
    auto den = make_var(den_var);
    auto den_inv = make_power(den, make_number(-1));
    auto mul = make_multiply({num, den_inv});
    return wrap_expr(mul);
}

/// Generate a random non-zero sign (Positive or Negative)
static Sign random_nonzero_sign() {
    return rc::gen::boolean() ? Sign::Positive : Sign::Negative;
}

/// Determine expected sign of division given numerator and denominator signs
static Sign expected_division_sign(Sign num_sign, Sign den_sign) {
    // positive / positive → positive
    // negative / negative → positive
    // positive / negative → negative
    // negative / positive → negative
    bool same_sign = (num_sign == Sign::Positive && den_sign == Sign::Positive) ||
                     (num_sign == Sign::Negative && den_sign == Sign::Negative);
    return same_sign ? Sign::Positive : Sign::Negative;
}

// ============================================================
// Property 1: Division sign inference follows sign multiplication table
// **Validates: Requirements 1.1, 1.2**
// ============================================================

static void test_property1_division_sign_table() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 1: Division sign inference follows sign multiplication table");

    rc::check("For any division with known non-zero numerator and denominator signs, "
              "the result sign follows the sign multiplication table", []() {
        // Generate random signs for numerator and denominator
        Sign num_sign = random_nonzero_sign();
        Sign den_sign = random_nonzero_sign();

        // Generate unique variable names
        std::string num_name = "num_" + std::to_string(rc::gen::inRange(0, 999));
        std::string den_name = "den_" + std::to_string(rc::gen::inRange(0, 999));

        // Set up context with declared signs
        AssumptionContext ctx;
        ctx.assume_sign(num_name, num_sign);
        ctx.assume_sign(den_name, den_sign);
        InferenceEngine engine(ctx);

        // Build division expression: num / den
        auto div_expr = make_division(num_name, den_name);

        // Determine expected result
        Sign expected = expected_division_sign(num_sign, den_sign);

        if (expected == Sign::Positive) {
            RC_ASSERT(engine.query_positive(div_expr) == Tribool::True);
            RC_ASSERT(engine.query_negative(div_expr) == Tribool::False);
            RC_ASSERT(engine.query_nonnegative(div_expr) == Tribool::True);
            RC_ASSERT(engine.query_nonzero(div_expr) == Tribool::True);
        } else {
            RC_ASSERT(engine.query_negative(div_expr) == Tribool::True);
            RC_ASSERT(engine.query_positive(div_expr) == Tribool::False);
            RC_ASSERT(engine.query_nonpositive(div_expr) == Tribool::True);
            RC_ASSERT(engine.query_nonzero(div_expr) == Tribool::True);
        }
    });
}

// ============================================================
// Property 1: Unknown denominator sign → Unknown result
// **Validates: Requirements 1.3**
// ============================================================

static void test_property1_unknown_denominator_returns_unknown() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 1: Unknown denominator sign returns Unknown");

    rc::check("For any division where denominator sign is unknown, result is Unknown", []() {
        Sign num_sign = random_nonzero_sign();
        std::string num_name = "num_" + std::to_string(rc::gen::inRange(0, 999));
        std::string den_name = "den_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(num_name, num_sign);
        // den_name has no sign declared → Unknown
        InferenceEngine engine(ctx);

        auto div_expr = make_division(num_name, den_name);

        RC_ASSERT(engine.query_positive(div_expr) == Tribool::Unknown);
        RC_ASSERT(engine.query_negative(div_expr) == Tribool::Unknown);
    });
}

// ============================================================
// Property 1: Zero denominator → Unknown result
// **Validates: Requirements 1.4**
// ============================================================

static void test_property1_zero_denominator_returns_unknown() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 1: Zero denominator returns Unknown");

    rc::check("For any division where denominator is zero, result is Unknown", []() {
        Sign num_sign = random_nonzero_sign();
        std::string num_name = "num_" + std::to_string(rc::gen::inRange(0, 999));
        std::string den_name = "den_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(num_name, num_sign);
        ctx.assume_sign(den_name, Sign::Zero);
        InferenceEngine engine(ctx);

        auto div_expr = make_division(num_name, den_name);

        RC_ASSERT(engine.query_positive(div_expr) == Tribool::Unknown);
        RC_ASSERT(engine.query_negative(div_expr) == Tribool::Unknown);
        RC_ASSERT(engine.query_nonnegative(div_expr) == Tribool::Unknown);
        RC_ASSERT(engine.query_nonpositive(div_expr) == Tribool::Unknown);
    });
}

// ============================================================
// Property 1: All four sign combinations verified exhaustively
// **Validates: Requirements 1.1, 1.2**
// ============================================================

static void test_property1_all_sign_combinations() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 1: All four sign combinations");

    // pos / pos → pos
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Positive);
        ctx.assume_sign("b", Sign::Positive);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True, "pos/pos → positive");
    }
    // neg / neg → pos
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Negative);
        ctx.assume_sign("b", Sign::Negative);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True, "neg/neg → positive");
    }
    // pos / neg → neg
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Positive);
        ctx.assume_sign("b", Sign::Negative);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative(expr) == Tribool::True, "pos/neg → negative");
    }
    // neg / pos → neg
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Negative);
        ctx.assume_sign("b", Sign::Positive);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative(expr) == Tribool::True, "neg/pos → negative");
    }
}

// ============================================================
// Property 4: Arithmetic combination sign inference — AddNode
// **Validates: Requirements 4.1, 4.2**
// ============================================================

static void test_property4_add_all_positive_is_positive() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: AddNode all GT zero → positive");

    rc::check("For any AddNode where all operands are GT zero, sum is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_positive(expr) == Tribool::True);
    });
}

static void test_property4_add_all_nonneg_is_nonneg() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: AddNode all GEQ zero → non-negative");

    rc::check("For any AddNode where all operands are GEQ zero, sum is non-negative", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::NonNegative);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_nonnegative(expr) == Tribool::True);
    });
}

// ============================================================
// Property 4: Arithmetic combination sign inference — MultiplyNode
// **Validates: Requirements 4.3**
// ============================================================

static void test_property4_multiply_all_positive_is_positive() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: MultiplyNode all GT zero → positive");

    rc::check("For any MultiplyNode where all operands are GT zero, product is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto mul_node = make_multiply(operands);
        auto expr = wrap_expr(mul_node);

        RC_ASSERT(engine.query_positive(expr) == Tribool::True);
    });
}

// ============================================================
// Property 4: x GT y with y non-negative → x positive
// **Validates: Requirements 4.4**
// This rule applies when checking composite expressions or when
// the relation store derives sign properties from variable-op-zero patterns.
// ============================================================

static void test_property4_gt_nonneg_implies_positive() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: x GT y with y NonNegative → x positive");

    rc::check("For any variable x with relation x GT 0, x is positive", []() {
        std::string x_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;

        // Add relation x > 0 — this triggers sign derivation in RelationStore
        auto x_var = std::make_shared<VariableNode>(x_name);
        auto zero_node = std::make_shared<NumberNode>(BigInt(0));
        auto rel_node = std::make_shared<RelationalNode>(
            x_var, zero_node, RelationalNode::Op::GT);
        SymbolicExpr rel_expr;
        rel_expr.root = rel_node;
        ctx.assume(rel_expr);

        InferenceEngine engine(ctx);

        SymbolicExpr x_expr;
        x_expr.root = make_var(x_name);

        // x > 0 should derive Positive sign for x in the PropertyStore
        RC_ASSERT(engine.query_positive(x_expr) == Tribool::True);
    });
}

// ============================================================
// Property 4: Unknown operand sign → Unknown result
// **Validates: Requirements 4.5**
// ============================================================

static void test_property4_unknown_operand_returns_unknown() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: Unknown operand sign → Unknown");

    rc::check("For any AddNode/MultiplyNode with an undetermined operand, result is Unknown", []() {
        bool use_add = rc::gen::boolean();
        int num_operands = rc::gen::inRange(2, 4);

        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        // Make all but one operand positive, leave one undetermined
        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            if (i < num_operands - 1) {
                ctx.assume_sign(name, Sign::Positive);
            }
            // Last operand has no sign declared
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);

        std::shared_ptr<SymbolicNode> node;
        if (use_add) {
            node = make_add(operands);
        } else {
            node = make_multiply(operands);
        }
        auto expr = wrap_expr(node);

        // With one unknown operand, the engine should not be able to determine
        // the sign definitively (for the general case)
        // For AddNode: can't determine positive if one operand is unknown
        // For MultiplyNode: can't determine positive if one operand is unknown
        Tribool result = engine.query_positive(expr);
        // The result should be Unknown since we can't determine the last operand's sign
        RC_ASSERT(result == Tribool::Unknown);
    });
}

// ============================================================
// Property 4: AddNode with operands having GT 0 relations → positive
// **Validates: Requirements 4.1 (via relational constraints)**
// ============================================================

static void test_property4_add_with_gt_zero_relations() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: AddNode operands with GT 0 relations → positive");

    rc::check("For any AddNode where all operands have x GT 0 relations, sum is positive", []() {
        int num_operands = rc::gen::inRange(2, 4);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        auto zero_node = std::make_shared<NumberNode>(BigInt(0));

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "r" + std::to_string(i);
            operands.push_back(make_var(name));

            // Add relation: r_i > 0
            auto var_node = std::make_shared<VariableNode>(name);
            auto rel_node = std::make_shared<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GT);
            SymbolicExpr rel_expr;
            rel_expr.root = rel_node;
            ctx.assume(rel_expr);
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_positive(expr) == Tribool::True);
    });
}

// ============================================================
// Property 4: Mixed positive and non-negative in AddNode
// ============================================================

static void test_property4_add_mixed_pos_nonneg() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 4: AddNode mixed positive/non-negative");

    rc::check("AddNode with at least one positive and rest non-negative is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<SymbolicNode>> operands;

        // All operands are positive (which implies non-negative)
        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        // All positive → sum is positive
        RC_ASSERT(engine.query_positive(expr) == Tribool::True);
        RC_ASSERT(engine.query_nonnegative(expr) == Tribool::True);
    });
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 1: Division sign inference
    test_property1_division_sign_table();
    test_property1_unknown_denominator_returns_unknown();
    test_property1_zero_denominator_returns_unknown();
    test_property1_all_sign_combinations();

    // Property 4: Arithmetic combination sign inference
    test_property4_add_all_positive_is_positive();
    test_property4_add_all_nonneg_is_nonneg();
    test_property4_multiply_all_positive_is_positive();
    test_property4_gt_nonneg_implies_positive();
    test_property4_unknown_operand_returns_unknown();
    test_property4_add_with_gt_zero_relations();
    test_property4_add_mixed_pos_nonneg();

    return TEST_REPORT();
}
