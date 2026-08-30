
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;


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


void test_single_zero_operand() {
    TEST_CASE("Single zero operand makes product Zero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(0, x) - zero operand detected via is_zero()
    auto mul_node = make_multiply({make_number(0), make_var("x")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "0 * x is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::True,
        "0 * x is NonPositive");
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "0 * x is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "0 * x is not Negative");
    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::False,
        "0 * x is not NonZero");
}

void test_zero_among_multiple_operands() {
    TEST_CASE("Zero among multiple operands makes product Zero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("x"), make_number(0), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "x * 0 * y * z is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::True,
        "x * 0 * y * z is NonPositive");
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "x * 0 * y * z is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "x * 0 * y * z is not Negative");
}

void test_zero_with_positive_numbers() {
    TEST_CASE("Zero with positive number operands");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(5, 0, 3) - zero among positive numbers
    auto mul_node = make_multiply(
        {make_number(5), make_number(0), make_number(3)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "5 * 0 * 3: not Positive");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "5 * 0 * 3: NonNegative");
}

void test_zero_rational_and_float() {
    TEST_CASE("Zero as Rational(0) and 0.0 detected");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Rational(0) * x
    {
        auto zero_rat = lamina::detail::make_node<NumberNode>(Rational(0));
        auto mul_node = make_multiply({zero_rat, make_var("x")});
        auto expr = wrap_expr(mul_node);
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
            "Rational(0) * x: not Positive");
        EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
            "Rational(0) * x: NonNegative");
    }
    // 0.0 * x
    {
        auto zero_float = lamina::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(0.0));
        auto mul_node = make_multiply({zero_float, make_var("x")});
        auto expr = wrap_expr(mul_node);
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
            "0.0 * x: not Positive");
        EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
            "0.0 * x: NonNegative");
    }
}


void test_two_positives_product_positive() {
    TEST_CASE("positive * positive = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, 5) - both positive, 0 negatives (even)
    auto mul_node = make_multiply({make_number(3), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "3 * 5 is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "3 * 5 is not Negative");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "3 * 5 is NonNegative");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::False,
        "3 * 5 is not NonPositive");
}

void test_one_negative_product_negative() {
    TEST_CASE("positive * negative = Negative (odd negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, -5) - 1 negative (odd)
    auto mul_node = make_multiply({make_number(3), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "3 * (-5) is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "3 * (-5) is Negative");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::False,
        "3 * (-5) is not NonNegative");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::True,
        "3 * (-5) is NonPositive");
}

void test_two_negatives_product_positive() {
    TEST_CASE("negative * negative = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-3, -5) - 2 negatives (even)
    auto mul_node = make_multiply({make_number(-3), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "(-3) * (-5) is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "(-3) * (-5) is not Negative");
}

void test_three_negatives_product_negative() {
    TEST_CASE("neg * neg * neg = Negative (odd negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-2, -3, -4) - 3 negatives (odd)
    auto mul_node = make_multiply(
        {make_number(-2), make_number(-3), make_number(-4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "(-2)*(-3)*(-4) is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "(-2)*(-3)*(-4) is Negative");
}

void test_four_negatives_product_positive() {
    TEST_CASE("4 negatives = Positive (even negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-1, -2, -3, -4) - 4 negatives (even)
    auto mul_node = make_multiply(
        {make_number(-1), make_number(-2), make_number(-3), make_number(-4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "(-1)*(-2)*(-3)*(-4) is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "(-1)*(-2)*(-3)*(-4) is not Negative");
}

void test_mixed_positive_negative_even() {
    TEST_CASE("pos * neg * pos * neg = Positive (2 negatives)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4, -5) - 2 negatives (even)
    auto mul_node = make_multiply(
        {make_number(2), make_number(-3), make_number(4), make_number(-5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "2*(-3)*4*(-5) is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "2*(-3)*4*(-5) is not Negative");
}

void test_mixed_positive_negative_odd() {
    TEST_CASE("pos * neg * pos = Negative (1 negative)");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2, -3, 4) - 1 negative (odd)
    auto mul_node = make_multiply(
        {make_number(2), make_number(-3), make_number(4)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "2*(-3)*4 is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "2*(-3)*4 is Negative");
}

// --- Sign parity with variables that have declared signs ---

void test_positive_variables_product() {
    TEST_CASE("Positive variables product is Positive");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Positive);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "pos_x * pos_y is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "pos_x * pos_y is not Negative");
}

void test_negative_variables_product() {
    TEST_CASE("Negative variables product is Positive (even neg)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "neg_x * neg_y is Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "neg_x * neg_y is not Negative");
}

void test_pos_neg_variable_product() {
    TEST_CASE("Positive * Negative variable = Negative (odd neg)");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "pos_x * neg_y is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "pos_x * neg_y is Negative");
}

void test_three_negative_variables() {
    TEST_CASE("3 Negative variables = Negative (odd neg)");
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Negative);
    ctx.assume_sign("b", Sign::Negative);
    ctx.assume_sign("c", Sign::Negative);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("a"), make_var("b"), make_var("c")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "neg_a * neg_b * neg_c is not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "neg_a * neg_b * neg_c is Negative");
}


void test_nonzero_variables_product() {
    TEST_CASE("All NonZero variables -> product NonZero");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);  // Positive implies NonZero
    ctx.assume_sign("y", Sign::Negative);  // Negative implies NonZero
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
        "pos_x * neg_y is NonZero");
}

void test_nonzero_numbers_product() {
    TEST_CASE("All nonzero numbers -> product NonZero");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_number(3), make_number(-7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
        "3 * (-7) is NonZero");
}


void test_nonneg_even_negatives() {
    TEST_CASE("NonNeg operands + even negatives -> NonNeg");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::Negative);
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) * z(neg) - 2 negatives (even)
    auto mul_node = make_multiply(
        {make_var("x"), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "nonneg * neg * neg: NonNegative (even negatives)");
}

void test_nonpos_odd_negatives() {
    TEST_CASE("NonNeg operands + odd negatives -> NonPos");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);
    ctx.assume_sign("y", Sign::Negative);
    InferenceEngine engine(ctx);

    // x(nonneg) * y(neg) - 1 negative (odd)
    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::True,
        "nonneg * neg: NonPositive (odd negatives)");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::False,
        "nonneg * neg: not NonNegative");
}


void test_unknown_sign_no_zero() {
    TEST_CASE("Unknown sign operands (no zero) -> Unknown");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(x, y) - no properties declared
    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
        "x * y: Positive Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
        "x * y: Negative Unknown");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
        "x * y: NonNegative Unknown");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::Unknown,
        "x * y: NonPositive Unknown");
}

void test_one_unknown_among_known() {
    TEST_CASE("One Unknown among known-sign operands -> Unknown");
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    // y has no sign declared -> Unknown
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
        "pos_x * unknown_y: Positive Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
        "pos_x * unknown_y: Negative Unknown");
}

void test_zero_overrides_unknown() {
    TEST_CASE("Zero overrides Unknown sign operands");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(x, 0) - x unknown but zero present
    auto mul_node = make_multiply({make_var("x"), make_number(0)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "x * 0: not Positive (zero overrides)");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "x * 0: not Negative (zero overrides)");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
        "x * 0: NonNegative (zero)");
    EXPECT_TRUE(engine.query_nonpositive_checked(expr).value() == Tribool::True,
        "x * 0: NonPositive (zero)");
}


void test_empty_operands() {
    TEST_CASE("Empty MultiplyNode is rejected");
    bool rejected = false;
    try {
        (void)lamina::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "Empty MultiplyNode violates the AST invariant");
}

void test_single_positive_number() {
    TEST_CASE("Single positive number operand");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(7) - single positive number
    auto mul_node = make_multiply({make_number(7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "multiply(7): Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
        "multiply(7): not Negative");
}

void test_single_negative_number() {
    TEST_CASE("Single negative number operand");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(-7)
    auto mul_node = make_multiply({make_number(-7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "multiply(-7): not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "multiply(-7): Negative");
}

void test_real_numbers_sign() {
    TEST_CASE("Real number operands sign inference");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2.5, -1.5) - 1 negative (odd)
    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::False,
        "2.5 * (-1.5): not Positive");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::True,
        "2.5 * (-1.5): Negative");
}

void test_many_operands_sign_parity() {
    TEST_CASE("5 positive numbers = Positive");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_number(1), make_number(2), make_number(3),
         make_number(4), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
        "1*2*3*4*5: Positive");
    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
        "1*2*3*4*5: NonZero");
}


void test_all_integer_numbers() {
    TEST_CASE("All integer NumberNodes -> product is Integer");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(3, -5, 7) - all BigInt (Integer)
    auto mul_node = make_multiply(
        {make_number(3), make_number(-5), make_number(7)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::True,
        "3 * (-5) * 7: Integer");
}

void test_all_integer_variables() {
    TEST_CASE("All Integer-domain variables -> product is Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Integer);
    ctx.assume_domain("z", Domain::Integer);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_var("x"), make_var("y"), make_var("z")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::True,
        "int_x * int_y * int_z: Integer");
}

void test_mixed_integer_and_number() {
    TEST_CASE("Integer variable * integer number -> Integer");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    // multiply(x, 5) - x is Integer, 5 is BigInt (Integer)
    auto mul_node = make_multiply({make_var("x"), make_number(5)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::True,
        "int_x * 5: Integer");
}

void test_all_real_numbers() {
    TEST_CASE("All real NumberNodes -> product is Real");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(2.5, -1.5, 3.0) - all finite reals
    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5),
         make_number_real(3.0)});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "2.5 * (-1.5) * 3.0: Real");
}

void test_all_real_variables() {
    TEST_CASE("All Real-domain variables -> product is Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "real_x * real_y: Real");
}

void test_integer_implies_real() {
    TEST_CASE("Integer operands also satisfy Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Integer);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    // Integer implies Real, so product of Integers is also Real
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "int_x * int_y: also Real (Integer subset of Real)");
}

void test_mixed_integer_real_is_real() {
    TEST_CASE("Integer * Real = Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    ctx.assume_domain("y", Domain::Real);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "int_x * real_y: Real");
    // But not necessarily Integer
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown,
        "int_x * real_y: Integer is Unknown");
}

void test_unknown_domain_propagation() {
    TEST_CASE("Unknown domain operand -> Unknown domain result");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    // y has no domain declared -> Unknown (defaults to Complex)
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown,
        "int_x * unknown_y: Integer Unknown");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
        "int_x * unknown_y: Real Unknown");
}

void test_number_and_real_variable() {
    TEST_CASE("Integer number * Real variable = Real");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // multiply(5, x) - 5 is BigInt (Integer->Real), x is Real
    auto mul_node = make_multiply({make_number(5), make_var("x")});
    auto expr = wrap_expr(mul_node);

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "5 * real_x: Real");
}

void test_natural_domain_implies_integer() {
    TEST_CASE("Natural domain variables -> Integer product");
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Natural);
    ctx.assume_domain("y", Domain::Natural);
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({make_var("x"), make_var("y")});
    auto expr = wrap_expr(mul_node);

    // Natural implies Integer, so product should be Integer
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::True,
        "nat_x * nat_y: Integer (Natural implies Integer)");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "nat_x * nat_y: Real (Natural implies Real)");
}

void test_empty_multiply_domain() {
    TEST_CASE("Empty MultiplyNode has no domain query state");
    bool rejected = false;
    try {
        (void)lamina::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "Invalid empty products are rejected before inference");
}

void test_rational_not_integer() {
    TEST_CASE("Rational(1/2) operand -> not Integer");
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // multiply(Rational(1,2), 3)
    auto rat_node = lamina::detail::make_node<NumberNode>(Rational(1, 2));
    auto mul_node = make_multiply({rat_node, make_number(3)});
    auto expr = wrap_expr(mul_node);

    /// Rational(1/2) 与整数 3 的积缺少整数性证明,两项仍具有 Real 域.
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown,
        "Rational(1/2) * 3: Integer Unknown (non-integer operand)");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
        "Rational(1/2) * 3: Real");
}


int main() {
    test_single_zero_operand();
    test_zero_among_multiple_operands();
    test_zero_with_positive_numbers();
    test_zero_rational_and_float();

    test_two_positives_product_positive();
    test_one_negative_product_negative();
    test_two_negatives_product_positive();
    test_three_negatives_product_negative();
    test_four_negatives_product_positive();
    test_mixed_positive_negative_even();
    test_mixed_positive_negative_odd();

    test_positive_variables_product();
    test_negative_variables_product();
    test_pos_neg_variable_product();
    test_three_negative_variables();

    test_nonzero_variables_product();
    test_nonzero_numbers_product();

    test_nonneg_even_negatives();
    test_nonpos_odd_negatives();

    test_unknown_sign_no_zero();
    test_one_unknown_among_known();
    test_zero_overrides_unknown();

    test_empty_operands();
    test_single_positive_number();
    test_single_negative_number();
    test_real_numbers_sign();
    test_many_operands_sign_parity();

    test_all_integer_numbers();
    test_all_integer_variables();
    test_mixed_integer_and_number();
    test_all_real_numbers();
    test_all_real_variables();
    test_integer_implies_real();
    test_mixed_integer_real_is_real();
    test_unknown_domain_propagation();
    test_number_and_real_variable();
    test_natural_domain_implies_integer();
    test_empty_multiply_domain();
    test_rational_not_integer();

    return TEST_REPORT();
}
