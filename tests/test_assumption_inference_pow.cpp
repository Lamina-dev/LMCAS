
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>
#include <vector>
#include <string>
#include <random>

using namespace LMCAS;

/// Create a VariableNode

/// Create a NumberNode from an integer
static std::shared_ptr<const SymbolicNode> make_num(int v) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(v));
}

/// Create a NumberNode from a double
static std::shared_ptr<const SymbolicNode> make_num_d(double v) {
    return LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(v));
}

/// Create a PowerNode expression
static SymbolicExpr make_power_expr(std::shared_ptr<const SymbolicNode> base,
                                    std::shared_ptr<const SymbolicNode> exponent) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exponent)));
    return expr;
}

TEST(AssumptionInferencePow, PositiveBaseRealExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // x^y where x is Positive, y is Real
    auto expr = make_power_expr(test_variable_node("x"), test_variable_node("y"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "x^y is Positive when x>0 and y is Real";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "x^y is NonNegative when x>0 and y is Real";
    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "x^y is NonZero when x>0 and y is Real";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "x^y is not Negative when x>0 and y is Real";
}

TEST(AssumptionInferencePow, PositiveBaseIntegerExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    // Test with multiple integer exponents (integers are Real)
    std::vector<int> exponents = {1, 2, 3, 5, 10, -1, -2, -3};
    for (int exp : exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Positive when x>0";
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << msg;
    }
}

TEST(AssumptionInferencePow, PositiveBaseRealNumberExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    // Real number exponents (non-integer)
    std::vector<double> exponents = {0.5, 1.5, 2.7, -0.5, -1.5, 3.14};
    for (double exp : exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num_d(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Positive when x>0";
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << msg;
    }
}

TEST(AssumptionInferencePow, NumericPositiveBase) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // 2^y, 5^y, 100^y where y is Real
    std::vector<int> bases = {1, 2, 5, 10, 100};
    for (int b : bases) {
        auto expr = make_power_expr(make_num(b), test_variable_node("y"));
        std::string msg = std::to_string(b) + "^y is Positive when y is Real";
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << msg;
    }
}

TEST(AssumptionInferencePow, RealBaseEvenExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // Test with multiple even exponents
    std::vector<int> even_exponents = {2, 4, 6, 8, 10, 20, 100};
    for (int exp : even_exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonNegative when x is Real (even exponent)";
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << msg;
        std::string msg2 = "x^" + std::to_string(exp) +
                           " is not Negative when x is Real (even exponent)";
        EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << msg2;
    }
}

TEST(AssumptionInferencePow, RealBaseOddExponentNotNonneg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // Odd exponents: x^3, x^5 - Real base with odd exponent can be negative
    // so NonNegative should be Unknown (not True)
    std::vector<int> odd_exponents = {1, 3, 5, 7};
    for (int exp : odd_exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Unknown for NonNegative when x is Real (odd exponent)";
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << msg;
    }
}

TEST(AssumptionInferencePow, RandomizedEvenExponents) {
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(1, 50);
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int half = dist(rng);
        int even_exp = half * 2; // Always even
        SCOPED_TRACE(::testing::Message() << "iteration=" << i << " exponent=" << even_exp);

        AssumptionContext ctx;
        std::string var_name = "x" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_domain(var_name, Domain::Real).has_value());

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(test_variable_node(var_name), make_num(even_exp));

        EXPECT_EQ(engine.query_nonnegative_checked(expr).value(), Tribool::True)
            << "random even exponents yield NonNegative";
    }
}

TEST(AssumptionInferencePow, NonnegBasePositiveIntExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    // Test with multiple positive integer exponents
    std::vector<int> pos_exponents = {1, 2, 3, 4, 5, 10, 50};
    for (int exp : pos_exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonNegative when x>=0 (positive int exponent)";
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << msg;
        std::string msg2 = "x^" + std::to_string(exp) +
                           " is not Negative when x>=0";
        EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << msg2;
    }
}

TEST(AssumptionInferencePow, RandomizedPositiveExponents) {
    std::mt19937 rng(123);
    std::uniform_int_distribution<int> dist(1, 100);
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int pos_exp = dist(rng); // Always positive
        SCOPED_TRACE(::testing::Message() << "iteration=" << i << " exponent=" << pos_exp);

        AssumptionContext ctx;
        std::string var_name = "v" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::NonNegative).has_value());

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(test_variable_node(var_name), make_num(pos_exp));

        EXPECT_EQ(engine.query_nonnegative_checked(expr).value(), Tribool::True)
            << "random positive exponents yield NonNegative";
    }
}

TEST(AssumptionInferencePow, NonzeroBaseIntegerExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonZero).has_value());

    InferenceEngine engine(ctx);

    // Test with various integer exponents (positive, negative, even, odd)
    std::vector<int> exponents = {1, 2, 3, -1, -2, -3, 5, 10, -10};
    for (int exp : exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonZero when x!=0 (integer exponent)";
        EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << msg;
    }
}

TEST(AssumptionInferencePow, NonzeroRandomizedIntegerExponents) {
    std::mt19937 rng(456);
    std::uniform_int_distribution<int> dist(-50, 50);
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int exp = dist(rng);
        SCOPED_TRACE(::testing::Message() << "iteration=" << i << " exponent=" << exp);

        AssumptionContext ctx;
        std::string var_name = "z" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::NonZero).has_value());

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(test_variable_node(var_name), make_num(exp));

        EXPECT_EQ(engine.query_nonzero_checked(expr).value(), Tribool::True)
            << "random integer exponents yield NonZero";
    }
}

TEST(AssumptionInferencePow, PositiveBaseEvenExponentIsPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // x^2 where x is Positive and Real -> Positive (rule 16a fires)
    auto expr = make_power_expr(test_variable_node("x"), make_num(2));
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "x^2 is Positive when x>0 (rule 16a)";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "x^2 is NonNegative when x>0";
}

TEST(AssumptionInferencePow, NumericBaseAndExponent) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // 2^3 - NumberNode base (positive), NumberNode exponent (integer)
    auto expr = make_power_expr(make_num(2), make_num(3));
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "2^3 is Positive";
    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "2^3 is NonZero";

    // 3^(-2) - positive base, integer exponent
    auto expr2 = make_power_expr(make_num(3), make_num(-2));
    EXPECT_TRUE((engine.query_positive_checked(expr2).value() == Tribool::True)) << "3^(-2) is Positive";
}

TEST(AssumptionInferencePow, NoRuleMatches) {
    AssumptionContext ctx;
    // No assumptions about x or y
    InferenceEngine engine(ctx);

    auto expr = make_power_expr(test_variable_node("x"), test_variable_node("y"));
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "x^y is Unknown when no assumptions";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "x^y is Unknown for Negative when no assumptions";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "x^y is Unknown for NonNegative when no assumptions";
    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::Unknown)) << "x^y is Unknown for NonZero when no assumptions";
}

TEST(AssumptionInferencePow, NegativeBaseNonIntegerExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_power_expr(test_variable_node("x"), make_num_d(1.5));
    auto positive = engine.query_positive_checked(expr);
    auto real = engine.query_real_checked(expr);
    EXPECT_TRUE((!positive && positive.error().code == CasErrc::DomainError)) << "invalid real power sign";
    EXPECT_TRUE((!real && real.error().code == CasErrc::DomainError)) << "invalid real power domain";
}

TEST(AssumptionInferencePow, RealBaseIntegerExponentDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // Test with multiple integer exponents
    std::vector<int> exponents = {0, 1, 2, 3, -1, -2, -3, 5, 10, -10};
    for (int exp : exponents) {
        auto expr = make_power_expr(test_variable_node("x"), make_num(exp));
        auto result = engine.query_real_checked(expr);
        EXPECT_TRUE((result && result.value() == (exp > 0 ? Tribool::True : Tribool::Unknown))) << "zero and negative exponents need a nonzero base";
    }
}

TEST(AssumptionInferencePow, DomainRandomizedIntegerExponents) {
    std::mt19937 rng(789);
    std::uniform_int_distribution<int> dist(-50, 50);
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int exp = dist(rng);
        SCOPED_TRACE(::testing::Message() << "iteration=" << i << " exponent=" << exp);

        AssumptionContext ctx;
        std::string var_name = "r" + std::to_string(i);
        EXPECT_TRUE(ctx.assume_domain(var_name, Domain::Real).has_value());

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(test_variable_node(var_name), make_num(exp));

        auto result = engine.query_real_checked(expr);
        EXPECT_TRUE(result && result.value() == (exp > 0 ? Tribool::True : Tribool::Unknown))
            << "integer exponent domains respect possible zero bases";
    }
}

TEST(AssumptionInferencePow, IntegerBaseIntegerExponent) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    // Integer is a subset of Real, so this should also yield Real
    std::vector<int> exponents = {1, 2, 3, -1, 5};
    for (int exp : exponents) {
        auto expr = make_power_expr(test_variable_node("n"), make_num(exp));
        auto result = engine.query_real_checked(expr);
        EXPECT_TRUE((result && result.value() == (exp > 0 ? Tribool::True : Tribool::Unknown))) << "Integer alone does not exclude a zero base";
    }
}

TEST(AssumptionInferencePow, NumericBaseIntegerExponent) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Numeric bases are Real, integer exponents are integers
    auto expr1 = make_power_expr(make_num(2), make_num(3));
    EXPECT_TRUE((engine.query_real_checked(expr1).value() == Tribool::True)) << "2^3 is Real";

    auto expr2 = make_power_expr(make_num(-3), make_num(2));
    EXPECT_TRUE((engine.query_real_checked(expr2).value() == Tribool::True)) << "(-3)^2 is Real";

    auto expr3 = make_power_expr(make_num(5), make_num(-1));
    EXPECT_TRUE((engine.query_real_checked(expr3).value() == Tribool::True)) << "5^(-1) is Real";
}

TEST(AssumptionInferencePow, RealBaseNonIntegerExponentUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // x^(1.5) where x is Real - non-integer exponent, rule doesn't apply
    // (could be complex if x < 0)
    auto expr = make_power_expr(test_variable_node("x"), make_num_d(1.5));
    /// 该域推导规则要求整数指数;仅有 Real 基底时结果域保持 Unknown.
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "x^1.5 is Unknown for Real when x is only Real (non-integer exponent)";
}

TEST(AssumptionInferencePow, NoDomainYieldsUnknown) {
    AssumptionContext ctx;
    // No assumptions about x
    InferenceEngine engine(ctx);

    auto expr = make_power_expr(test_variable_node("x"), make_num(2));
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "x^2 is Unknown for Real when x has no domain assumption";
}

TEST(AssumptionInferencePow, PowerDomainBoundaries) {
    for (Sign sign : {Sign::Zero, Sign::NonNegative, Sign::Negative}) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", sign).has_value()) << "base sign";
        InferenceEngine engine(ctx);
        for (int exponent : {0, -2, 2}) {
            auto expression = make_power_expr(test_variable_node("x"), make_num(exponent));
            auto result = engine.query_nonnegative_checked(expression);
            if (sign == Sign::Zero && exponent <= 0) {
                EXPECT_TRUE((!result && result.error().code == CasErrc::DomainError)) << "zero base with zero/negative exponent is undefined";
            } else {
                const auto expected = sign == Sign::NonNegative && exponent <= 0 ? Tribool::Unknown : Tribool::True;
                EXPECT_TRUE((result && result.value() == expected)) << "nonzero obligation precedes even-power nonnegativity";
            }
        }
    }
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "positive logarithm argument";
    InferenceEngine engine(ctx);
    auto logarithm = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                                     std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x")});
    auto square = make_power_expr(logarithm, make_num(2));
    auto good = engine.query_nonnegative_checked(square);
    EXPECT_TRUE((good && good.value() == Tribool::True)) << "defined real logarithm squared is nonnegative";
    auto invalid_log = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                                       std::vector<std::shared_ptr<const SymbolicNode>>{make_num(-1)});
    for (auto expression : {make_power_expr(invalid_log, make_num(2)),
                            make_power_expr(make_num(2), invalid_log),
                            make_power_expr(invalid_log, make_num(0))}) {
        auto result = engine.query_positive_checked(expression);
        EXPECT_TRUE((!result && result.error().code == CasErrc::DomainError)) << "undefined base or exponent is never a positive-power proof";
    }
}
