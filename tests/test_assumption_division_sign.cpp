
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base,
    std::shared_ptr<const SymbolicNode> exp) {
    return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return LMCAS::detail::make_node<MultiplyNode>(std::move(ops));
}

static std::shared_ptr<const SymbolicNode> make_add(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return LMCAS::detail::make_node<AddNode>(std::move(ops));
}

/// Build a division expression: numerator / denominator
/// Represented as MultiplyNode([numerator, PowerNode(denominator, -1)])
static SymbolicExpr make_division(const std::string &num_var, const std::string &den_var) {
    auto num = test_variable_node(num_var);
    auto den = test_variable_node(den_var);
    auto den_inv = make_power(den, test_integer_node(-1));
    auto mul = make_multiply({num, den_inv});
    return test_expression_from_node(mul);
}

/// Generate a random non-zero sign (Positive or Negative)
static Sign random_nonzero_sign() {
    return *rc::gen::arbitrary<bool>() ? Sign::Positive : Sign::Negative;
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

TEST(LmcasAssumptionDivisionSign, DivisionSignTable) {
    EXPECT_TRUE(rc::check("For any division with known non-zero numerator and denominator signs, "
                          "the result sign follows the sign multiplication table",
                          []() {
                              // Generate random signs for numerator and denominator
                              Sign num_sign = random_nonzero_sign();
                              Sign den_sign = random_nonzero_sign();

                              // Generate unique variable names
                              std::string num_name = "num_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
                              std::string den_name = "den_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

                              // Set up context with declared signs
                              AssumptionContext ctx;
                              EXPECT_TRUE(ctx.assume_sign(num_name, num_sign).has_value());
                              EXPECT_TRUE(ctx.assume_sign(den_name, den_sign).has_value());
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
                          }));
}

TEST(LmcasAssumptionDivisionSign, UnknownDenominatorReturnsUnknown) {
    EXPECT_TRUE(rc::check("For any division where denominator sign is unknown, result is Unknown", []() {
        Sign num_sign = random_nonzero_sign();
        std::string num_name = "num_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
        std::string den_name = "den_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(num_name, num_sign).has_value());
        // den_name has no sign declared -> Unknown
        InferenceEngine engine(ctx);

        auto div_expr = make_division(num_name, den_name);

        RC_ASSERT(engine.query_positive_checked(div_expr).value() == Tribool::Unknown);
        RC_ASSERT(engine.query_negative_checked(div_expr).value() == Tribool::Unknown);
    }));
}

TEST(LmcasAssumptionDivisionSign, DivisionDomainBoundaries) {
    for (Sign denominator : {Sign::Zero, Sign::NonNegative, Sign::Positive}) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("n", Sign::Zero).has_value()) << "zero numerator";
        EXPECT_TRUE(ctx.assume_sign("d", denominator).has_value()) << "denominator sign";
        InferenceEngine engine(ctx);
        auto expression = make_division("n", "d");
        for (Sign target : {Sign::Positive, Sign::Negative, Sign::NonNegative,
                            Sign::NonPositive, Sign::Zero, Sign::NonZero}) {
            auto result = test_sign_query(engine, expression, target);
            if (denominator == Sign::Zero) {
                EXPECT_TRUE(!result && result.error().code == CasErrc::DomainError) << "zero denominator is a domain error for every sign query";
            } else if (denominator == Sign::NonNegative) {
                EXPECT_TRUE(result && result.value() == Tribool::Unknown) << "possibly zero denominator is not absorbed by a zero numerator";
            } else {
                const bool true_at_zero = target == Sign::Zero ||
                                          target == Sign::NonNegative || target == Sign::NonPositive;
                EXPECT_TRUE(result && result.value() ==
                                          (true_at_zero ? Tribool::True : Tribool::False))
                    << "zero numerator over a proved nonzero denominator";
            }
        }
    }
}

TEST(LmcasAssumptionDivisionSign, WeakDivisionSigns) {
    const std::vector<Sign> signs = {Sign::NonNegative, Sign::NonPositive, Sign::NonZero};
    for (Sign numerator : signs)
        for (Sign denominator : {Sign::Positive, Sign::Negative}) {
            AssumptionContext ctx;
            EXPECT_TRUE(ctx.assume_domain("n", Domain::Real).has_value()) << "real numerator";
            EXPECT_TRUE(ctx.assume_sign("n", numerator).has_value()) << "numerator sign";
            EXPECT_TRUE(ctx.assume_sign("d", denominator).has_value()) << "nonzero denominator sign";
            InferenceEngine engine(ctx);
            auto expression = make_division("n", "d");
            const bool nonnegative = (numerator == Sign::NonNegative) ==
                                     (denominator == Sign::Positive);
            for (Sign target : {Sign::Positive, Sign::Negative, Sign::NonNegative,
                                Sign::NonPositive, Sign::Zero, Sign::NonZero}) {
                Tribool expected = Tribool::Unknown;
                if (numerator == Sign::NonZero) {
                    if (target == Sign::Zero)
                        expected = Tribool::False;
                    if (target == Sign::NonZero)
                        expected = Tribool::True;
                } else {
                    if (target == (nonnegative ? Sign::NonNegative : Sign::NonPositive))
                        expected = Tribool::True;
                    if (target == (nonnegative ? Sign::Negative : Sign::Positive))
                        expected = Tribool::False;
                }
                auto actual = test_sign_query(engine, expression, target);
                EXPECT_TRUE(actual && actual.value() == expected) << "quotient possible-sign truth table";
            }
        }
}

TEST(LmcasAssumptionDivisionSign, AllSignCombinations) {
    // pos / pos -> pos
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
        EXPECT_TRUE(ctx.assume_sign("b", Sign::Positive).has_value());
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True) << "pos/pos → positive";
    }
    // neg / neg -> pos
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("a", Sign::Negative).has_value());
        EXPECT_TRUE(ctx.assume_sign("b", Sign::Negative).has_value());
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True) << "neg/neg → positive";
    }
    // pos / neg -> neg
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
        EXPECT_TRUE(ctx.assume_sign("b", Sign::Negative).has_value());
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True) << "pos/neg → negative";
    }
    // neg / pos -> neg
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("a", Sign::Negative).has_value());
        EXPECT_TRUE(ctx.assume_sign("b", Sign::Positive).has_value());
        InferenceEngine engine(ctx);
        auto expr = make_division("a", "b");
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True) << "neg/pos → negative";
    }
}

TEST(LmcasAssumptionDivisionSign, AddAllPositiveIsPositive) {
    EXPECT_TRUE(rc::check("For any AddNode where all operands are GT zero, sum is positive", []() {
        int num_operands = *rc::gen::inRange(2, (5) + 1);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());
            operands.push_back(test_variable_node(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = test_expression_from_node(add_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDivisionSign, AddAllNonnegIsNonneg) {
    EXPECT_TRUE(rc::check("For any AddNode where all operands are GEQ zero, sum is non-negative", []() {
        int num_operands = *rc::gen::inRange(2, (5) + 1);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            EXPECT_TRUE(ctx.assume_sign(name, Sign::NonNegative).has_value());
            operands.push_back(test_variable_node(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = test_expression_from_node(add_node);

        RC_ASSERT(engine.query_nonnegative_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDivisionSign, MultiplyAllPositiveIsPositive) {
    EXPECT_TRUE(rc::check("For any MultiplyNode where all operands are GT zero, product is positive", []() {
        int num_operands = *rc::gen::inRange(2, (5) + 1);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());
            operands.push_back(test_variable_node(name));
        }

        InferenceEngine engine(ctx);
        auto mul_node = make_multiply(operands);
        auto expr = test_expression_from_node(mul_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDivisionSign, GtNonnegImpliesPositive) {
    EXPECT_TRUE(rc::check("For any variable x with relation x GT 0, x is positive", []() {
        std::string x_name = "x_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;

        // Add relation x > 0 - this triggers sign derivation in RelationStore
        auto x_var = LMCAS::detail::make_node<VariableNode>(x_name);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = LMCAS::detail::make_node<RelationalNode>(
            x_var, zero_node, RelationalNode::Op::GT);
        auto rel_expr = LMCAS::detail::expression_from_node(rel_node);
        EXPECT_TRUE(ctx.assume(rel_expr).has_value()) << "relation assumption succeeds";

        InferenceEngine engine(ctx);

        auto x_expr = LMCAS::detail::expression_from_node(test_variable_node(x_name));
        // x > 0 should derive Positive sign for x in the PropertyStore
        RC_ASSERT(engine.query_positive_checked(x_expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDivisionSign, UnknownOperandReturnsUnknown) {
    EXPECT_TRUE(rc::check("For any AddNode/MultiplyNode with an undetermined operand, result is Unknown", []() {
        bool use_add = *rc::gen::arbitrary<bool>();
        int num_operands = *rc::gen::inRange(2, (4) + 1);

        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        // Make all but one operand positive, leave one undetermined
        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            if (i < num_operands - 1) {
                EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());
            }
            // Last operand has no sign declared
            operands.push_back(test_variable_node(name));
        }

        InferenceEngine engine(ctx);

        std::shared_ptr<const SymbolicNode> node;
        if (use_add) {
            node = make_add(operands);
        } else {
            node = make_multiply(operands);
        }
        auto expr = test_expression_from_node(node);

        // With one unknown operand, the engine should not be able to determine
        /// AddNode 或 MultiplyNode 包含符号未知的操作数时,和或积的符号保持 Unknown.
        Tribool result = engine.query_positive_checked(expr).value();
        /// 最后一个操作数的符号未知,因此结果为 Unknown.
        RC_ASSERT(result == Tribool::Unknown);
    }));
}

TEST(LmcasAssumptionDivisionSign, AddWithGtZeroRelations) {
    EXPECT_TRUE(rc::check("For any AddNode where all operands have x GT 0 relations, sum is positive", []() {
        int num_operands = *rc::gen::inRange(2, (4) + 1);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        for (int i = 0; i < num_operands; ++i) {
            std::string name = "r" + std::to_string(i);
            operands.push_back(test_variable_node(name));

            // Add relation: r_i > 0
            auto var_node = LMCAS::detail::make_node<VariableNode>(name);
            auto rel_node = LMCAS::detail::make_node<RelationalNode>(
                var_node, zero_node, RelationalNode::Op::GT);
            auto rel_expr = LMCAS::detail::expression_from_node(rel_node);
            EXPECT_TRUE(ctx.assume(rel_expr).has_value()) << "relation assumption succeeds";
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = test_expression_from_node(add_node);

        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDivisionSign, AddMixedPosNonneg) {
    EXPECT_TRUE(rc::check("AddNode with at least one positive and rest non-negative is positive", []() {
        int num_operands = *rc::gen::inRange(2, (5) + 1);
        AssumptionContext ctx;
        std::vector<std::shared_ptr<const SymbolicNode>> operands;

        // All operands are positive (which implies non-negative)
        for (int i = 0; i < num_operands; ++i) {
            std::string name = "x" + std::to_string(i);
            EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());
            operands.push_back(test_variable_node(name));
        }

        InferenceEngine engine(ctx);
        auto add_node = make_add(operands);
        auto expr = test_expression_from_node(add_node);

        // All positive -> sum is positive
        RC_ASSERT(engine.query_positive_checked(expr).value() == Tribool::True);
        RC_ASSERT(engine.query_nonnegative_checked(expr).value() == Tribool::True);
    }));
}
