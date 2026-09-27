
#include "test_assumption_inference_mul_support.hpp"

static std::shared_ptr<const SymbolicNode> make_number_real(double val) {
    return LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val));
}

TEST(AssumptionInferenceMul, SingleZeroOperand) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real operand accepted";
    InferenceEngine engine(ctx);

    // multiply(0, x) - zero operand detected via is_zero()
    auto mul_node = make_multiply({test_integer_node(0), test_variable_node("x")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "0 * x is NonNegative";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "0 * x is NonPositive";
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "0 * x is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "0 * x is not Negative";
    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::False)) << "0 * x is not NonZero";
}

TEST(AssumptionInferenceMul, ZeroAmongMultipleOperands) {
    AssumptionContext ctx;
    for (const auto *name : {"x", "y", "z"})
        EXPECT_TRUE(ctx.assume_domain(name, Domain::Real).has_value()) << "real operand accepted";
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {test_variable_node("x"), test_integer_node(0), test_variable_node("y"), test_variable_node("z")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "x * 0 * y * z is NonNegative";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "x * 0 * y * z is NonPositive";
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "x * 0 * y * z is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "x * 0 * y * z is not Negative";
}

TEST(AssumptionInferenceMul, ZeroWithPositiveNumbers) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(5, 0, 3) - zero among positive numbers
    auto mul_node = make_multiply(
        {test_integer_node(5), test_integer_node(0), test_integer_node(3)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "5 * 0 * 3: not Positive";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "5 * 0 * 3: NonNegative";
}

TEST(AssumptionInferenceMul, ZeroRationalAndFloat) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real operand accepted";
    InferenceEngine engine(ctx);

    // Rational(0) * x
    {
        auto zero_rat = LMCAS::detail::make_node<NumberNode>(Rational(0));
        auto mul_node = make_multiply({zero_rat, test_variable_node("x")});
        auto expr = test_expression_from_node(mul_node);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "Rational(0) * x: not Positive";
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "Rational(0) * x: NonNegative";
    }
    // 0.0 * x
    {
        auto zero_float = LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(0.0));
        auto mul_node = make_multiply({zero_float, test_variable_node("x")});
        auto expr = test_expression_from_node(mul_node);
        EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "0.0 * x: not Positive";
        EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "0.0 * x: NonNegative";
    }
}

TEST(AssumptionInferenceMul, TwoPositivesProductPositive) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, 5) - both positive, 0 negatives (even)
    auto mul_node = make_multiply({test_integer_node(3), test_integer_node(5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "3 * 5 is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "3 * 5 is not Negative";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "3 * 5 is NonNegative";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::False)) << "3 * 5 is not NonPositive";
}

TEST(AssumptionInferenceMul, OneNegativeProductNegative) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, -5) - 1 negative (odd)
    auto mul_node = make_multiply({test_integer_node(3), test_integer_node(-5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "3 * (-5) is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "3 * (-5) is Negative";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::False)) << "3 * (-5) is not NonNegative";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "3 * (-5) is NonPositive";
}

TEST(AssumptionInferenceMul, TwoNegativesProductPositive) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-3, -5) - 2 negatives (even)
    auto mul_node = make_multiply({test_integer_node(-3), test_integer_node(-5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "(-3) * (-5) is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "(-3) * (-5) is not Negative";
}

TEST(AssumptionInferenceMul, ThreeNegativesProductNegative) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-2, -3, -4) - 3 negatives (odd)
    auto mul_node = make_multiply(
        {test_integer_node(-2), test_integer_node(-3), test_integer_node(-4)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "(-2)*(-3)*(-4) is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "(-2)*(-3)*(-4) is Negative";
}

TEST(AssumptionInferenceMul, FourNegativesProductPositive) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-1, -2, -3, -4) - 4 negatives (even)
    auto mul_node = make_multiply(
        {test_integer_node(-1), test_integer_node(-2), test_integer_node(-3), test_integer_node(-4)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "(-1)*(-2)*(-3)*(-4) is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "(-1)*(-2)*(-3)*(-4) is not Negative";
}

TEST(AssumptionInferenceMul, MixedPositiveNegativeEven) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4, -5) - 2 negatives (even)
    auto mul_node = make_multiply(
        {test_integer_node(2), test_integer_node(-3), test_integer_node(4), test_integer_node(-5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "2*(-3)*4*(-5) is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "2*(-3)*4*(-5) is not Negative";
}

TEST(AssumptionInferenceMul, MixedPositiveNegativeOdd) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4) - 1 negative (odd)
    auto mul_node = make_multiply(
        {test_integer_node(2), test_integer_node(-3), test_integer_node(4)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "2*(-3)*4 is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "2*(-3)*4 is Negative";
}

// --- Sign parity with variables that have declared signs ---

TEST(AssumptionInferenceMul, PositiveVariablesProduct) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "pos_x * pos_y is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "pos_x * pos_y is not Negative";
}

TEST(AssumptionInferenceMul, NegativeVariablesProduct) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "neg_x * neg_y is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "neg_x * neg_y is not Negative";
}

TEST(AssumptionInferenceMul, PosNegVariableProduct) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "pos_x * neg_y is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "pos_x * neg_y is Negative";
}

TEST(AssumptionInferenceMul, ThreeNegativeVariables) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("c", Sign::Negative).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {test_variable_node("a"), test_variable_node("b"), test_variable_node("c")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "neg_a * neg_b * neg_c is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "neg_a * neg_b * neg_c is Negative";
}

TEST(AssumptionInferenceMul, NonzeroVariablesProduct) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()); /**< 正数蕴含非零。 */
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value()); /**< 负数蕴含非零。 */
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "pos_x * neg_y is NonZero";
}

TEST(AssumptionInferenceMul, NonzeroNumbersProduct) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_integer_node(3), test_integer_node(-7)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "3 * (-7) is NonZero";
}

TEST(AssumptionInferenceMul, NonnegEvenNegatives) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_sign("z", Sign::Negative).has_value());
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) * z(neg) - 2 negatives (even)
    auto mul_node = make_multiply(
        {test_variable_node("x"), test_variable_node("y"), test_variable_node("z")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "nonneg * neg * neg: NonNegative (even negatives)";
}

TEST(AssumptionInferenceMul, NonposOddNegatives) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) - 1 negative (odd)
    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "nonneg * neg: NonPositive (odd negatives)";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "nonneg * neg may be zero or negative";
}

TEST(AssumptionInferenceMul, UnknownSignNoZero) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(x, y) - no properties declared
    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "x * y: Positive Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "x * y: Negative Unknown";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "x * y: NonNegative Unknown";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::Unknown)) << "x * y: NonPositive Unknown";
}

TEST(AssumptionInferenceMul, OneUnknownAmongKnown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    // y has no sign declared -> Unknown
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "pos_x * unknown_y: Positive Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "pos_x * unknown_y: Negative Unknown";
}

TEST(AssumptionInferenceMul, ZeroOverridesUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real operand accepted";
    InferenceEngine engine(ctx);

    // multiply(x, 0) - x unknown but zero present
    auto mul_node = make_multiply({test_variable_node("x"), test_integer_node(0)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "x * 0: not Positive (zero overrides)";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "x * 0: not Negative (zero overrides)";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "x * 0: NonNegative (zero)";
    EXPECT_TRUE((engine.query_nonpositive_checked(expr).value() == Tribool::True)) << "x * 0: NonPositive (zero)";
}

TEST(AssumptionInferenceMul, EmptyOperands) {
    bool rejected = false;
    try {
        (void)LMCAS::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "Empty MultiplyNode violates the AST invariant";
}

TEST(AssumptionInferenceMul, SinglePositiveNumber) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(7) - single positive number
    auto mul_node = make_multiply({test_integer_node(7)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "multiply(7): Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "multiply(7): not Negative";
}

TEST(AssumptionInferenceMul, SingleNegativeNumber) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-7)
    auto mul_node = make_multiply({test_integer_node(-7)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "multiply(-7): not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "multiply(-7): Negative";
}

TEST(AssumptionInferenceMul, RealNumbersSign) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2.5, -1.5) - 1 negative (odd)
    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "2.5 * (-1.5): not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "2.5 * (-1.5): Negative";
}

TEST(AssumptionInferenceMul, ManyOperandsSignParity) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {test_integer_node(1), test_integer_node(2), test_integer_node(3),
         test_integer_node(4), test_integer_node(5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "1*2*3*4*5: Positive";
    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "1*2*3*4*5: NonZero";
}
