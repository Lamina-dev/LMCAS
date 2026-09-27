
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>

using namespace LMCAS;

/// Create a VariableNode wrapped in a shared_ptr<SymbolicNode>

/// Create a NumberNode from a BigInt value
static std::shared_ptr<const SymbolicNode> make_num(int val) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(val));
}

/// Create an AddNode from a vector of operands (bypasses factory simplification)
static std::shared_ptr<const AddNode> make_add(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return LMCAS::detail::make_node<AddNode>(std::move(ops));
}

/// Wrap an AddNode into a SymbolicExpr for querying

TEST(AssumptionInferenceAdd, AllPositiveOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("c", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    // Two positive operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "a + b is Positive when a, b are Positive";
    }

    // Three positive operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b"), test_variable_node("c")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "a + b + c is Positive when a, b, c are Positive";
    }
}

TEST(AssumptionInferenceAdd, AllNegativeOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("z", Sign::Negative).has_value());

    InferenceEngine engine(ctx);

    // Two negative operands
    {
        auto add = make_add({test_variable_node("x"), test_variable_node("y")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "x + y is Negative when x, y are Negative";
    }

    // Three negative operands
    {
        auto add = make_add({test_variable_node("x"), test_variable_node("y"), test_variable_node("z")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "x + y + z is Negative when x, y, z are Negative";
    }
}

TEST(AssumptionInferenceAdd, AllNonnegativeOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::NonNegative).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::NonNegative).has_value());
    ASSERT_TRUE(ctx.assume_sign("c", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    // Two nonnegative operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "a + b is NonNegative when a, b are NonNegative";
    }

    // Three nonnegative operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b"), test_variable_node("c")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "a + b + c is NonNegative when a, b, c are NonNegative";
    }
}

TEST(AssumptionInferenceAdd, AllNonpositiveOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::NonPositive).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::NonPositive).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "a + b is NonPositive when a, b are NonPositive";
}

TEST(AssumptionInferenceAdd, PositiveImpliesNonnegativeForSum) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);

    // Positive implies NonNegative, so the sum should also be NonNegative
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "a + b is NonNegative when a, b are Positive (Positive implies NonNegative)";
}

TEST(AssumptionInferenceAdd, UnknownOperandYieldsUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
    // "b" has no sign declared -> Unknown

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for Positive when b has Unknown sign";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for Negative when b has Unknown sign";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for NonNegative when b has Unknown sign";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for NonPositive when b has Unknown sign";
}

TEST(AssumptionInferenceAdd, MixedSignsYieldUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("pos", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("neg", Sign::Negative).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("pos"), test_variable_node("neg")});
    auto expr = test_expression_from_node(add);

    /// 正数与负数之和的符号由幅值决定,因此保持 Unknown.
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "pos + neg is Unknown for Positive (mixed signs)";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "pos + neg is Unknown for Negative (mixed signs)";
}

TEST(AssumptionInferenceAdd, SingleOperand) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("x")});
    auto expr = test_expression_from_node(add);

    // Note: AddNode with single operand may be simplified by the factory,
    // but we construct it directly here
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "Single-operand add(x) is Positive when x is Positive";
}

TEST(AssumptionInferenceAdd, SignWithNumberOperands) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Add two positive numbers: 3 + 5
    {
        auto add = make_add({make_num(3), make_num(5)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "3 + 5 is Positive";
    }

    // Add two negative numbers: (-3) + (-5)
    {
        auto add = make_add({make_num(-3), make_num(-5)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "(-3) + (-5) is Negative";
    }

    // Mixed: 3 + (-5) -> Unknown
    {
        auto add = make_add({make_num(3), make_num(-5)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "3 + (-5) is Unknown for Positive (mixed signs)";
    }
}

TEST(AssumptionInferenceAdd, WithVariablesAndNumbersMixed) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    // x + 5 where x is Positive -> sum is Positive
    {
        auto add = make_add({test_variable_node("x"), make_num(5)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "x + 5 is Positive when x is Positive";
    }

    // x + (-3) where x is Positive -> Unknown (mixed)
    {
        auto add = make_add({test_variable_node("x"), make_num(-3)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "x + (-3) is Unknown for Positive (mixed signs)";
    }
}

TEST(AssumptionInferenceAdd, ManyOperandsUniformSign) {
    AssumptionContext ctx;
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    for (int i = 0; i < 10; ++i) {
        std::string name = "v" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());
        ops.push_back(test_variable_node(name));
    }

    InferenceEngine engine(ctx);

    auto add = make_add(ops);
    auto expr = test_expression_from_node(add);
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "Sum of 10 Positive variables is Positive";
}

TEST(AssumptionInferenceAdd, EmptyAddReturnsUnknown) {
    bool rejected = false;
    try {
        (void)LMCAS::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "Empty AddNode violates the AST invariant";
}

TEST(AssumptionInferenceAdd, AllIntegerOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("b", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("c", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    // Two integer operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "a + b is Integer when a, b are Integer";
    }

    // Three integer operands
    {
        auto add = make_add({test_variable_node("a"), test_variable_node("b"), test_variable_node("c")});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "a + b + c is Integer when a, b, c are Integer";
    }
}

TEST(AssumptionInferenceAdd, AllRealOperands) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(add);
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "x + y is Real when x, y are Real";
}

TEST(AssumptionInferenceAdd, IntegerImpliesRealForSum) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("b", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);

    // Integer implies Real, so sum of integers should also be Real
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "a + b is Real when a, b are Integer (Integer implies Real)";
}

TEST(AssumptionInferenceAdd, MixedIntegerAndReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("b", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);

    // Integer + Real -> Real (Integer is subset of Real)
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "a + b is Real when a is Integer and b is Real";

    // But not necessarily Integer (b might not be Integer)
    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for Integer when b is only Real";
}

TEST(AssumptionInferenceAdd, UnknownDomainYieldsUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());
    // "b" has no domain declared (defaults to Complex)

    InferenceEngine engine(ctx);

    auto add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto expr = test_expression_from_node(add);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "a + b is Unknown for Integer when b has no Integer domain";
}

TEST(AssumptionInferenceAdd, DomainWithNumberOperands) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // 3 + 5 (both BigInt -> Integer)
    {
        auto add = make_add({make_num(3), make_num(5)});
        auto expr = test_expression_from_node(add);
        EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "3 + 5 is Integer";
        EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "3 + 5 is Real";
    }
}

TEST(AssumptionInferenceAdd, ManyIntegerOperands) {
    AssumptionContext ctx;
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    for (int i = 0; i < 8; ++i) {
        std::string name = "n" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_domain(name, Domain::Integer).has_value());
        ops.push_back(test_variable_node(name));
    }

    InferenceEngine engine(ctx);

    auto add = make_add(ops);
    auto expr = test_expression_from_node(add);
    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "Sum of 8 Integer variables is Integer";
}

TEST(AssumptionInferenceAdd, RealWithNumbers) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // x + 3 where x is Real -> sum is Real (3 is Integer which implies Real)
    auto add = make_add({test_variable_node("x"), make_num(3)});
    auto expr = test_expression_from_node(add);
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "x + 3 is Real when x is Real";
}

TEST(AssumptionInferenceAdd, NestedAdditionDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("b", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("c", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    // (a + b) + c - the inner add should be Integer, so the outer should too
    auto inner_add = make_add({test_variable_node("a"), test_variable_node("b")});
    auto outer_add = make_add({inner_add, test_variable_node("c")});
    auto expr = test_expression_from_node(outer_add);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "(a + b) + c is Integer when a, b, c are Integer";
}
