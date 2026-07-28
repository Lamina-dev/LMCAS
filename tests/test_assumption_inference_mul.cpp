/**
 * @file test_assumption_inference_mul.cpp
 * @brief Property tests for InferenceEngine multiplication inference (Properties 14-15).
 *
 * Feature: assumption-system
 * Validates: Requirements 6.1-6.9
 *
 * Property 14: Multiplication sign inference
 *   For any MultiplyNode: (a) if any operand is Zero, the product is Zero;
 *   (b) if all operands have definite sign, the product's sign is determined by
 *   the parity of the count of Negative operands (even -> Positive/NonNegative,
 *   odd -> Negative/NonPositive); (c) if any operand has Unknown sign and none
 *   is Zero, the result is Unknown.
 *
 * Property 15: Multiplication domain closure
 *   For any MultiplyNode where all operands are Integer, the product should be
 *   Integer; where all operands are Real (or Integer), the product should be Real.
 */

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;

// ============================================================
// Helpers: create nodes and expressions
// ============================================================

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_number_real(double val) {
    return lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val));
}

static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<MultiplyNode>(std::move(ops));
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

// ============================================================
// Property 14a: Zero operand detection
// **Validates: Requirements 6.3**
// ============================================================

void test_property14a_single_zero_operand() {
    TEST_CASE("Property 14a: Single zero operand makes product Zero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(0, x) — zero operand detected via is_zero()
    auto mul_node = make_multiply({make_number(0), make_var("x")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "0 * x is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "0 * x is NonPositive");
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "0 * x is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "0 * x is not Negative");
    EXPECT_TRUE(engine.query_nonzero(expr) == Tribool::False,
        "0 * x is not NonZero");
}

void test_property14a_zero_among_multiple_operands() {
    TEST_CASE("Property 14a: Zero among multiple operands makes product Zero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("x"), make_number(0), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "x * 0 * y * z is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "x * 0 * y * z is NonPositive");
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "x * 0 * y * z is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "x * 0 * y * z is not Negative");
}

void test_property14a_zero_with_positive_numbers() {
    TEST_CASE("Property 14a: Zero with positive number operands");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(5, 0, 3) — zero among positive numbers
    auto mul_node = make_multiply(
        {make_number(5), make_number(0), make_number(3)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "5 * 0 * 3: not Positive");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "5 * 0 * 3: NonNegative");
}

void test_property14a_zero_rational_and_float() {
    TEST_CASE("Property 14a: Zero as Rational(0) and 0.0 detected");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Rational(0) * x
    {
        auto zero_rat = lamina::detail::make_node<NumberNode>(Rational(0));
        auto mul_node = make_multiply({zero_rat, make_var("x")});
        auto expr = wrap_expr(mul_node);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
            "Rational(0) * x: not Positive");
        EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
            "Rational(0) * x: NonNegative");
    }
    // 0.0 * x
    {
        auto zero_float = lamina::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(0.0));
        auto mul_node = make_multiply({zero_float, make_var("x")});
        auto expr = wrap_expr(mul_node);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
            "0.0 * x: not Positive");
        EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
            "0.0 * x: NonNegative");
    }
}

// ============================================================
// Property 14b: Sign parity of negatives (using NumberNodes)
// **Validates: Requirements 6.1, 6.2**
// The engine can determine sign of NumberNodes directly.
// ============================================================

void test_property14b_two_positives_product_positive() {
    TEST_CASE("Property 14b: positive * positive = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, 5) — both positive, 0 negatives (even)
    auto mul_node = make_multiply({make_number(3), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "3 * 5 is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "3 * 5 is not Negative");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "3 * 5 is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::False,
        "3 * 5 is not NonPositive");
}

void test_property14b_one_negative_product_negative() {
    TEST_CASE("Property 14b: positive * negative = Negative (odd negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, -5) — 1 negative (odd)
    auto mul_node = make_multiply({make_number(3), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "3 * (-5) is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "3 * (-5) is Negative");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::False,
        "3 * (-5) is not NonNegative");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "3 * (-5) is NonPositive");
}

void test_property14b_two_negatives_product_positive() {
    TEST_CASE("Property 14b: negative * negative = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-3, -5) — 2 negatives (even)
    auto mul_node = make_multiply({make_number(-3), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "(-3) * (-5) is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "(-3) * (-5) is not Negative");
}

void test_property14b_three_negatives_product_negative() {
    TEST_CASE("Property 14b: neg * neg * neg = Negative (odd negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-2, -3, -4) — 3 negatives (odd)
    auto mul_node = make_multiply(
        {make_number(-2), make_number(-3), make_number(-4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "(-2)*(-3)*(-4) is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "(-2)*(-3)*(-4) is Negative");
}

void test_property14b_four_negatives_product_positive() {
    TEST_CASE("Property 14b: 4 negatives = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-1, -2, -3, -4) — 4 negatives (even)
    auto mul_node = make_multiply(
        {make_number(-1), make_number(-2), make_number(-3), make_number(-4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "(-1)*(-2)*(-3)*(-4) is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "(-1)*(-2)*(-3)*(-4) is not Negative");
}

void test_property14b_mixed_positive_negative_even() {
    TEST_CASE("Property 14b: pos * neg * pos * neg = Positive (2 negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4, -5) — 2 negatives (even)
    auto mul_node = make_multiply(
        {make_number(2), make_number(-3), make_number(4), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "2*(-3)*4*(-5) is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "2*(-3)*4*(-5) is not Negative");
}

void test_property14b_mixed_positive_negative_odd() {
    TEST_CASE("Property 14b: pos * neg * pos = Negative (1 negative)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4) — 1 negative (odd)
    auto mul_node = make_multiply(
        {make_number(2), make_number(-3), make_number(4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "2*(-3)*4 is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "2*(-3)*4 is Negative");
}

// --- Sign parity with variables that have declared signs ---

void test_property14b_positive_variables_product() {
    TEST_CASE("Property 14b: Positive variables product is Positive");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "pos_x * pos_y is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "pos_x * pos_y is not Negative");
}

void test_property14b_negative_variables_product() {
    TEST_CASE("Property 14b: Negative variables product is Positive (even neg)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "neg_x * neg_y is Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "neg_x * neg_y is not Negative");
}

void test_property14b_pos_neg_variable_product() {
    TEST_CASE("Property 14b: Positive * Negative variable = Negative (odd neg)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "pos_x * neg_y is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "pos_x * neg_y is Negative");
}

void test_property14b_three_negative_variables() {
    TEST_CASE("Property 14b: 3 Negative variables = Negative (odd neg)");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Negative);
    ctx.assume_sign("c", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("a"), make_var("b"), make_var("c")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "neg_a * neg_b * neg_c is not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "neg_a * neg_b * neg_c is Negative");
}

// --- Req 6.7: All NonZero -> product is NonZero ---

void test_property14b_nonzero_variables_product() {
    TEST_CASE("Property 14b: All NonZero variables -> product NonZero (Req 6.7)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);  // Positive implies NonZero
    ctx.assume_sign("y", Sign::Negative);  // Negative implies NonZero
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonzero(expr) == Tribool::True,
        "pos_x * neg_y is NonZero");
}

void test_property14b_nonzero_numbers_product() {
    TEST_CASE("Property 14b: All nonzero numbers -> product NonZero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_number(3), make_number(-7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonzero(expr) == Tribool::True,
        "3 * (-7) is NonZero");
}

// --- Req 6.8, 6.9: NonNegative/NonPositive with parity ---

void test_property14b_nonneg_even_negatives() {
    TEST_CASE("Property 14b: NonNeg operands + even negatives -> NonNeg (Req 6.8)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::Negative);
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) * z(neg) — 2 negatives (even)
    auto mul_node = make_multiply(
        {make_var("x"), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "nonneg * neg * neg: NonNegative (even negatives)");
}

void test_property14b_nonpos_odd_negatives() {
    TEST_CASE("Property 14b: NonNeg operands + odd negatives -> NonPos (Req 6.9)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) — 1 negative (odd)
    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "nonneg * neg: NonPositive (odd negatives)");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::False,
        "nonneg * neg: not NonNegative");
}

// ============================================================
// Property 14c: Unknown sign propagation
// **Validates: Requirements 6.4**
// ============================================================

void test_property14c_unknown_sign_no_zero() {
    TEST_CASE("Property 14c: Unknown sign operands (no zero) -> Unknown");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(x, y) — no properties declared
    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "x * y: Positive Unknown");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "x * y: Negative Unknown");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::Unknown,
        "x * y: NonNegative Unknown");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::Unknown,
        "x * y: NonPositive Unknown");
}

void test_property14c_one_unknown_among_known() {
    TEST_CASE("Property 14c: One Unknown among known-sign operands -> Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    // y has no sign declared -> Unknown
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "pos_x * unknown_y: Positive Unknown");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "pos_x * unknown_y: Negative Unknown");
}

void test_property14c_zero_overrides_unknown() {
    TEST_CASE("Property 14c: Zero overrides Unknown sign operands");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(x, 0) — x unknown but zero present
    auto mul_node = make_multiply({make_var("x"), make_number(0)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "x * 0: not Positive (zero overrides)");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "x * 0: not Negative (zero overrides)");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "x * 0: NonNegative (zero)");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "x * 0: NonPositive (zero)");
}

// ============================================================
// Property 14: Edge cases
// ============================================================

void test_property14_empty_operands() {
    TEST_CASE("Property 14: Empty MultiplyNode is rejected");
    bool rejected = false;
    try {
        (void)lamina::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "Empty MultiplyNode violates the AST invariant");
}

void test_property14_single_positive_number() {
    TEST_CASE("Property 14: Single positive number operand");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(7) — single positive number
    auto mul_node = make_multiply({make_number(7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "multiply(7): Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::False,
        "multiply(7): not Negative");
}

void test_property14_single_negative_number() {
    TEST_CASE("Property 14: Single negative number operand");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-7)
    auto mul_node = make_multiply({make_number(-7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "multiply(-7): not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "multiply(-7): Negative");
}

void test_property14_real_numbers_sign() {
    TEST_CASE("Property 14: Real number operands sign inference");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2.5, -1.5) — 1 negative (odd)
    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::False,
        "2.5 * (-1.5): not Positive");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
        "2.5 * (-1.5): Negative");
}

void test_property14_many_operands_sign_parity() {
    TEST_CASE("Property 14: 5 positive numbers = Positive");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_number(1), make_number(2), make_number(3),
         make_number(4), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "1*2*3*4*5: Positive");
    EXPECT_TRUE(engine.query_nonzero(expr) == Tribool::True,
        "1*2*3*4*5: NonZero");
}

// ============================================================
// Property 15: Multiplication domain closure
// **Validates: Requirements 6.5, 6.6**
// ============================================================

void test_property15_all_integer_numbers() {
    TEST_CASE("Property 15: All integer NumberNodes -> product is Integer");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, -5, 7) — all BigInt (Integer)
    auto mul_node = make_multiply(
        {make_number(3), make_number(-5), make_number(7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "3 * (-5) * 7: Integer");
}

void test_property15_all_integer_variables() {
    TEST_CASE("Property 15: All Integer-domain variables -> product is Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Integer);
    ctx.assume_domain("z", Domain::Integer);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("x"), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "int_x * int_y * int_z: Integer");
}

void test_property15_mixed_integer_and_number() {
    TEST_CASE("Property 15: Integer variable * integer number -> Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    // multiply(x, 5) — x is Integer, 5 is BigInt (Integer)
    auto mul_node = make_multiply({make_var("x"), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "int_x * 5: Integer");
}

void test_property15_all_real_numbers() {
    TEST_CASE("Property 15: All real NumberNodes -> product is Real");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2.5, -1.5, 3.0) — all finite reals
    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5),
         make_number_real(3.0)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "2.5 * (-1.5) * 3.0: Real");
}

void test_property15_all_real_variables() {
    TEST_CASE("Property 15: All Real-domain variables -> product is Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "real_x * real_y: Real");
}

void test_property15_integer_implies_real() {
    TEST_CASE("Property 15: Integer operands also satisfy Real (Req 6.6)");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Integer);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    // Integer implies Real, so product of Integers is also Real
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "int_x * int_y: also Real (Integer subset of Real)");
}

void test_property15_mixed_integer_real_is_real() {
    TEST_CASE("Property 15: Integer * Real = Real (Req 6.6)");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Real);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "int_x * real_y: Real");
    // But not necessarily Integer
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "int_x * real_y: Integer is Unknown");
}

void test_property15_unknown_domain_propagation() {
    TEST_CASE("Property 15: Unknown domain operand -> Unknown domain result");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    // y has no domain declared -> Unknown (defaults to Complex)
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "int_x * unknown_y: Integer Unknown");
    EXPECT_TRUE(engine.query_real(expr) == Tribool::Unknown,
        "int_x * unknown_y: Real Unknown");
}

void test_property15_number_and_real_variable() {
    TEST_CASE("Property 15: Integer number * Real variable = Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // multiply(5, x) — 5 is BigInt (Integer->Real), x is Real
    auto mul_node = make_multiply({make_number(5), make_var("x")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "5 * real_x: Real");
}

void test_property15_natural_domain_implies_integer() {
    TEST_CASE("Property 15: Natural domain variables -> Integer product");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Natural);
    ctx.assume_domain("y", Domain::Natural);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    // Natural implies Integer, so product should be Integer
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "nat_x * nat_y: Integer (Natural implies Integer)");
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "nat_x * nat_y: Real (Natural implies Real)");
}

void test_property15_empty_multiply_domain() {
    TEST_CASE("Property 15: Empty MultiplyNode has no domain query state");
    bool rejected = false;
    try {
        (void)lamina::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "Invalid empty products are rejected before inference");
}

void test_property15_rational_not_integer() {
    TEST_CASE("Property 15: Rational(1/2) operand -> not Integer");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(Rational(1,2), 3)
    auto rat_node = lamina::detail::make_node<NumberNode>(Rational(1, 2));
    auto mul_node = make_multiply({rat_node, make_number(3)});
    auto expr = wrap_expr(mul_node);

    // Rational(1/2) is not Integer, so product can't be proven Integer
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "Rational(1/2) * 3: Integer Unknown (non-integer operand)");
    // But both are Real
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "Rational(1/2) * 3: Real");
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 14a: Zero operand detection (Req 6.3)
    test_property14a_single_zero_operand();
    test_property14a_zero_among_multiple_operands();
    test_property14a_zero_with_positive_numbers();
    test_property14a_zero_rational_and_float();

    // Property 14b: Sign parity with numbers (Req 6.1, 6.2)
    test_property14b_two_positives_product_positive();
    test_property14b_one_negative_product_negative();
    test_property14b_two_negatives_product_positive();
    test_property14b_three_negatives_product_negative();
    test_property14b_four_negatives_product_positive();
    test_property14b_mixed_positive_negative_even();
    test_property14b_mixed_positive_negative_odd();

    // Property 14b: Sign parity with variables (Req 6.1, 6.2)
    test_property14b_positive_variables_product();
    test_property14b_negative_variables_product();
    test_property14b_pos_neg_variable_product();
    test_property14b_three_negative_variables();

    // Property 14b: NonZero inference (Req 6.7)
    test_property14b_nonzero_variables_product();
    test_property14b_nonzero_numbers_product();

    // Property 14b: NonNeg/NonPos with parity (Req 6.8, 6.9)
    test_property14b_nonneg_even_negatives();
    test_property14b_nonpos_odd_negatives();

    // Property 14c: Unknown sign propagation (Req 6.4)
    test_property14c_unknown_sign_no_zero();
    test_property14c_one_unknown_among_known();
    test_property14c_zero_overrides_unknown();

    // Property 14: Edge cases
    test_property14_empty_operands();
    test_property14_single_positive_number();
    test_property14_single_negative_number();
    test_property14_real_numbers_sign();
    test_property14_many_operands_sign_parity();

    // Property 15: Domain closure (Req 6.5, 6.6)
    test_property15_all_integer_numbers();
    test_property15_all_integer_variables();
    test_property15_mixed_integer_and_number();
    test_property15_all_real_numbers();
    test_property15_all_real_variables();
    test_property15_integer_implies_real();
    test_property15_mixed_integer_real_is_real();
    test_property15_unknown_domain_propagation();
    test_property15_number_and_real_variable();
    test_property15_natural_domain_implies_integer();
    test_property15_empty_multiply_domain();
    test_property15_rational_not_integer();

    return TEST_REPORT();
}
