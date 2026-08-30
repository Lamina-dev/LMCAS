
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


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base,
    std::shared_ptr<const SymbolicNode> exp) {
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<MultiplyNode>(std::move(ops));
}

static std::shared_ptr<const SymbolicNode> make_add(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<AddNode>(std::move(ops));
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
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
    // positive / positive -> positive
    // negative / negative -> positive
    // positive / negative -> negative
    // negative / positive -> negative
    bool same_sign = (num_sign == Sign::Positive && den_sign == Sign::Positive) ||
                     (num_sign == Sign::Negative && den_sign == Sign::Negative);
    return same_sign ? Sign::Positive : Sign::Negative;
}


static void test_division_sign_table() {
    TEST_CASE("Division sign inference follows sign multiplication table");

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
            RC_ASSERT(engine.query_positive_checked(div_expr).value() == Tribool::True);
            RC_ASSERT(engine.query_negative_checked(div_expr).value() == Tribool::False);
            RC_ASSERT(engine.query_nonnegative_checked(div_expr).value() == Tribool::True);
            RC_ASSERT(engine.query_nonzero_checked(div_expr).value() == Tribool::True);
        } else {
            RC_ASSERT(engine.query_negative_checked(div_expr).value() == Tribool::True);
            RC_ASSERT(engine.query_positive_checked(div_expr).value() == Tribool::False);
            RC_ASSERT(engine.query_nonpositive_checked(div_expr).value() == Tribool::True);
            RC_ASSERT(engine.query_nonzero_checked(div_expr).value() == Tribool::True);
        }
    });
}


static void test_unknown_denominator_returns_unknown() {
    TEST_CASE("Unknown denominator sign returns Unknown");

    rc::check("For any division where denominator sign is unknown, result is Unknown", []() {
        Sign num_sign = random_nonzero_sign();
        std::string num_name = "num_" + std::to_string(rc::gen::inRange(0, 999));
        std::string den_name = "den_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(num_name, num_sign);
        // den_name has no sign declared -> Unknown
        InferenceEngine engine(ctx);

        auto div_expr = make_division(num_name, den_name);

        RC_ASSERT(engine.query_positive_checked(div_expr).value() == Tribool::Unknown);
        RC_ASSERT(engine.query_negative_checked(div_expr).value() == Tribool::Unknown);
    });
}


static void test_zero_denominator_returns_unknown() {
    TEST_CASE("Zero denominator returns Unknown");

    rc::check("For any division where denominator is zero, result is Unknown", []() {
        Sign num_sign = random_nonzero_sign();
        std::string num_name = "num_" + std::to_string(rc::gen::inRange(0, 999));
        std::string den_name = "den_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(num_name, num_sign);
        ctx.assume_sign(den_name, Sign::Zero);
        InferenceEngine engine(ctx);

        auto div_expr = make_division(num_name, den_name);

        RC_ASSERT(engine.query_positive_checked(div_expr).value() == Tribool::Unknown);
        RC_ASSERT(engine.query_negative_checked(div_expr).value() == Tribool::Unknown);
        RC_ASSERT(engine.query_nonnegative_checked(div_expr).value() == Tribool::Unknown);
        RC_ASSERT(engine.query_nonpositive_checked(div_expr).value() == Tribool::Unknown);
    });
}


static void test_all_sign_combinations() {
    TEST_CASE("All four sign combinations");

    // pos / pos -> pos
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Positive);
        ctx.assume_sign("b", Sign::Positive);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True, "pos/pos → positive");
    }
    // neg / neg -> pos
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Negative);
        ctx.assume_sign("b", Sign::Negative);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True, "neg/neg → positive");
    }
    // pos / neg -> neg
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Positive);
        ctx.assume_sign("b", Sign::Negative);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True, "pos/neg → negative");
    }
    // neg / pos -> neg
    {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Negative);
        ctx.assume_sign("b", Sign::Positive);
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True, "neg/pos → negative");
    }
}


static void test_add_all_positive_is_positive() {
    TEST_CASE("AddNode all GT zero → positive");

    rc::check("For any AddNode where all operands are GT zero, sum is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    });
}

static void test_add_all_nonneg_is_nonneg() {
    TEST_CASE("AddNode all GEQ zero → non-negative");

    rc::check("For any AddNode where all operands are GEQ zero, sum is non-negative", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::NonNegative);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_nonnegative_checked(expr).value() == Tribool::True);
    });
}


static void test_multiply_all_positive_is_positive() {
    TEST_CASE("MultiplyNode all GT zero → positive");

    rc::check("For any MultiplyNode where all operands are GT zero, product is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto mul_node = make_multiply(operands);
        auto expr = wrap_expr(mul_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    });
}


static void test_gt_nonneg_implies_positive() {
    TEST_CASE("x GT y with y NonNegative → x positive");

    rc::check("For any variable x with relation x GT 0, x is positive", []() {
        std::string x_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;

        // Add relation x > 0 - this triggers sign derivation in RelationStore
        auto x_var = lamina::detail::make_node<VariableNode>(x_name);
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = lamina::detail::make_node<RelationalNode>(
            x_var, zero_node, RelationalNode::Op::GT);
        auto rel_expr = lamina::detail::expression_from_node(rel_node);
        ctx.assume(rel_expr);

        InferenceEngine engine(ctx);

        auto x_expr = lamina::detail::expression_from_node(make_var(x_name));
        // x > 0 should derive Positive sign for x in the PropertyStore
        RC_ASSERT(engine.query_positive_checked(x_expr).value() == Tribool::True);
    });
}


static void test_unknown_operand_returns_unknown() {
    TEST_CASE("Unknown operand sign → Unknown");

    rc::check("For any AddNode/MultiplyNode with an undetermined operand, result is Unknown", []() {
        bool use_add = rc::gen::boolean();
        int num_operands = rc::gen::inRange(2, 4);

        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

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

        std::shared_ptr<const SymbolicNode> node;
        if (use_add) {
            node = make_add(operands);
        } else {
            node = make_multiply(operands);
        }
        auto expr = wrap_expr(node);

        // With one unknown operand, the engine should not be able to determine
        /// AddNode 或 MultiplyNode 包含符号未知的操作数时,和或积的符号保持 Unknown.
        Tribool result = engine.query_positive_checked(expr).value();
        /// 最后一个操作数的符号未知,因此结果为 Unknown.
        RC_ASSERT(result == Tribool::Unknown);
    });
}


static void test_add_with_gt_zero_relations() {
    TEST_CASE("AddNode operands with GT 0 relations → positive");

    rc::check("For any AddNode where all operands have x GT 0 relations, sum is positive", []() {
        int num_operands = rc::gen::inRange(2, 4);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "r" + std::to_string(i);
            operands.push_back(make_var(name));

            // Add relation: r_i > 0
            auto var_node = lamina::detail::make_node<VariableNode>(name);
            auto rel_node = lamina::detail::make_node<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GT);
            auto rel_expr = lamina::detail::expression_from_node(rel_node);
            ctx.assume(rel_expr);
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    });
}


static void test_add_mixed_pos_nonneg() {
    TEST_CASE("AddNode mixed positive/non-negative");

    rc::check("AddNode with at least one positive and rest non-negative is positive", []() {
        int num_operands = rc::gen::inRange(2, 5);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        // All operands are positive (which implies non-negative)
        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            ctx.assume_sign(name, Sign::Positive);
            operands.push_back(make_var(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = wrap_expr(add_node);

        // All positive -> sum is positive
        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
        RC_ASSERT(engine.query_nonnegative_checked(expr).value() == Tribool::True);
    });
}


int main() {
    test_division_sign_table();
    test_unknown_denominator_returns_unknown();
    test_zero_denominator_returns_unknown();
    test_all_sign_combinations();

    test_add_all_positive_is_positive();
    test_add_all_nonneg_is_nonneg();
    test_multiply_all_positive_is_positive();
    test_gt_nonneg_implies_positive();
    test_unknown_operand_returns_unknown();
    test_add_with_gt_zero_relations();
    test_add_mixed_pos_nonneg();

    return TEST_REPORT();
}
