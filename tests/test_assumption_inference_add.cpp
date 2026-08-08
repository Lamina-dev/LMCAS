
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>

using namespace lamina;


/// Create a VariableNode wrapped in a shared_ptr<SymbolicNode>
static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

/// Create a NumberNode from a BigInt value
static std::shared_ptr<const SymbolicNode> make_num(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

/// Create an AddNode from a vector of operands (bypasses factory simplification)
static std::shared_ptr<const AddNode> make_add(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return lamina::detail::make_node<AddNode>(std::move(ops));
}

/// Wrap an AddNode into a SymbolicExpr for querying
static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}


void test_property12_all_positive_operands() {
    TEST_CASE("Property 12: All Positive operands → sum is Positive");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);
    ctx.assume_sign("c", Sign::Positive);

    InferenceEngine engine(ctx);

    // Two positive operands
    {
        auto add = make_add({make_var("a"), make_var("b")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
            "a + b is Positive when a, b are Positive");
    }

    // Three positive operands
    {
        auto add = make_add({make_var("a"), make_var("b"), make_var("c")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
            "a + b + c is Positive when a, b, c are Positive");
    }
}

void test_property12_all_negative_operands() {
    TEST_CASE("Property 12: All Negative operands → sum is Negative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::Negative);

    InferenceEngine engine(ctx);

    // Two negative operands
    {
        auto add = make_add({make_var("x"), make_var("y")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
            "x + y is Negative when x, y are Negative");
    }

    // Three negative operands
    {
        auto add = make_add({make_var("x"), make_var("y"), make_var("z")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
            "x + y + z is Negative when x, y, z are Negative");
    }
}

void test_property12_all_nonnegative_operands() {
    TEST_CASE("Property 12: All NonNegative operands → sum is NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::NonNegative);
    ctx.assume_sign("b", Sign::NonNegative);
    ctx.assume_sign("c", Sign::NonNegative);

    InferenceEngine engine(ctx);

    // Two nonnegative operands
    {
        auto add = make_add({make_var("a"), make_var("b")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
            "a + b is NonNegative when a, b are NonNegative");
    }

    // Three nonnegative operands
    {
        auto add = make_add({make_var("a"), make_var("b"), make_var("c")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
            "a + b + c is NonNegative when a, b, c are NonNegative");
    }
}

void test_property12_all_nonpositive_operands() {
    TEST_CASE("Property 12: All NonPositive operands → sum is NonPositive");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::NonPositive);
    ctx.assume_sign("b", Sign::NonPositive);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::True,
        "a + b is NonPositive when a, b are NonPositive");
}

void test_property12_positive_implies_nonnegative_for_sum() {
    TEST_CASE("Property 12: All Positive operands → sum is also NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);

    // Positive implies NonNegative, so the sum should also be NonNegative
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::True,
        "a + b is NonNegative when a, b are Positive (Positive implies NonNegative)");
}

void test_property12_unknown_operand_yields_unknown() {
    TEST_CASE("Property 12: Any Unknown operand → result is Unknown");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    // "b" has no sign declared → Unknown

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);

    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "a + b is Unknown for Positive when b has Unknown sign");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "a + b is Unknown for Negative when b has Unknown sign");
    EXPECT_TRUE(engine.query_nonnegative(expr) == Tribool::Unknown,
        "a + b is Unknown for NonNegative when b has Unknown sign");
    EXPECT_TRUE(engine.query_nonpositive(expr) == Tribool::Unknown,
        "a + b is Unknown for NonPositive when b has Unknown sign");
}

void test_property12_mixed_signs_yield_unknown() {
    TEST_CASE("Property 12: Mixed definite signs → result is Unknown");

    AssumptionContext ctx;
    ctx.assume_sign("pos", Sign::Positive);
    ctx.assume_sign("neg", Sign::Negative);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("pos"), make_var("neg")});
    auto expr = wrap_expr(add);

    // Positive + Negative → can't determine sign of sum
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
        "pos + neg is Unknown for Positive (mixed signs)");
    EXPECT_TRUE(engine.query_negative(expr) == Tribool::Unknown,
        "pos + neg is Unknown for Negative (mixed signs)");
}

void test_property12_single_operand() {
    TEST_CASE("Property 12: Single operand AddNode preserves sign");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("x")});
    auto expr = wrap_expr(add);

    // Note: AddNode with single operand may be simplified by the factory,
    // but we construct it directly here
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "Single-operand add(x) is Positive when x is Positive");
}

void test_property12_with_number_operands() {
    TEST_CASE("Property 12: Addition with positive number operands");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Add two positive numbers: 3 + 5
    {
        auto add = make_add({make_num(3), make_num(5)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
            "3 + 5 is Positive");
    }

    // Add two negative numbers: (-3) + (-5)
    {
        auto add = make_add({make_num(-3), make_num(-5)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_negative(expr) == Tribool::True,
            "(-3) + (-5) is Negative");
    }

    // Mixed: 3 + (-5) → Unknown
    {
        auto add = make_add({make_num(3), make_num(-5)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
            "3 + (-5) is Unknown for Positive (mixed signs)");
    }
}

void test_property12_with_variables_and_numbers_mixed() {
    TEST_CASE("Property 12: Addition with variables and numbers");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // x + 5 where x is Positive → sum is Positive
    {
        auto add = make_add({make_var("x"), make_num(5)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
            "x + 5 is Positive when x is Positive");
    }

    // x + (-3) where x is Positive → Unknown (mixed)
    {
        auto add = make_add({make_var("x"), make_num(-3)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_positive(expr) == Tribool::Unknown,
            "x + (-3) is Unknown for Positive (mixed signs)");
    }
}

void test_property12_many_operands_uniform_sign() {
    TEST_CASE("Property 12: Many operands with uniform sign");

    AssumptionContext ctx;
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    for (int i = 0; i < 10; ++i) {
        std::string name = "v" + std::to_string(i);
        ctx.assume_sign(name, Sign::Positive);
        ops.push_back(make_var(name));
    }

    InferenceEngine engine(ctx);

    auto add = make_add(ops);
    auto expr = wrap_expr(add);
    EXPECT_TRUE(engine.query_positive(expr) == Tribool::True,
        "Sum of 10 Positive variables is Positive");
}

void test_property12_empty_add_returns_unknown() {
    TEST_CASE("Property 12: Empty AddNode is rejected");

    bool rejected = false;
    try {
        (void)lamina::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "Empty AddNode violates the AST invariant");
}


void test_property13_all_integer_operands() {
    TEST_CASE("Property 13: All Integer operands → sum is Integer");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Integer);
    ctx.assume_domain("b", Domain::Integer);
    ctx.assume_domain("c", Domain::Integer);

    InferenceEngine engine(ctx);

    // Two integer operands
    {
        auto add = make_add({make_var("a"), make_var("b")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
            "a + b is Integer when a, b are Integer");
    }

    // Three integer operands
    {
        auto add = make_add({make_var("a"), make_var("b"), make_var("c")});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
            "a + b + c is Integer when a, b, c are Integer");
    }
}

void test_property13_all_real_operands() {
    TEST_CASE("Property 13: All Real operands → sum is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("x"), make_var("y")});
    auto expr = wrap_expr(add);
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "x + y is Real when x, y are Real");
}

void test_property13_integer_implies_real_for_sum() {
    TEST_CASE("Property 13: All Integer operands → sum is also Real (Integer ⊂ Real)");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Integer);
    ctx.assume_domain("b", Domain::Integer);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);

    // Integer implies Real, so sum of integers should also be Real
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "a + b is Real when a, b are Integer (Integer implies Real)");
}

void test_property13_mixed_integer_and_real() {
    TEST_CASE("Property 13: Mixed Integer and Real → sum is Real but not necessarily Integer");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Integer);
    ctx.assume_domain("b", Domain::Real);

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);

    // Integer + Real → Real (Integer is subset of Real)
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "a + b is Real when a is Integer and b is Real");

    // But not necessarily Integer (b might not be Integer)
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "a + b is Unknown for Integer when b is only Real");
}

void test_property13_unknown_domain_yields_unknown() {
    TEST_CASE("Property 13: Unknown domain operand → result is Unknown");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Integer);
    // "b" has no domain declared (defaults to Complex)

    InferenceEngine engine(ctx);

    auto add = make_add({make_var("a"), make_var("b")});
    auto expr = wrap_expr(add);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "a + b is Unknown for Integer when b has no Integer domain");
}

void test_property13_with_number_operands() {
    TEST_CASE("Property 13: Addition with integer number operands");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // 3 + 5 (both BigInt → Integer)
    {
        auto add = make_add({make_num(3), make_num(5)});
        auto expr = wrap_expr(add);
        EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
            "3 + 5 is Integer");
        EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
            "3 + 5 is Real");
    }
}

void test_property13_many_integer_operands() {
    TEST_CASE("Property 13: Many Integer operands → sum is Integer");

    AssumptionContext ctx;
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    for (int i = 0; i < 8; ++i) {
        std::string name = "n" + std::to_string(i);
        ctx.assume_domain(name, Domain::Integer);
        ops.push_back(make_var(name));
    }

    InferenceEngine engine(ctx);

    auto add = make_add(ops);
    auto expr = wrap_expr(add);
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "Sum of 8 Integer variables is Integer");
}

void test_property13_real_with_numbers() {
    TEST_CASE("Property 13: Real variable + integer number → Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // x + 3 where x is Real → sum is Real (3 is Integer which implies Real)
    auto add = make_add({make_var("x"), make_num(3)});
    auto expr = wrap_expr(add);
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True,
        "x + 3 is Real when x is Real");
}

void test_property13_nested_addition_domain() {
    TEST_CASE("Property 13: Nested addition preserves Integer domain");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Integer);
    ctx.assume_domain("b", Domain::Integer);
    ctx.assume_domain("c", Domain::Integer);

    InferenceEngine engine(ctx);

    // (a + b) + c — the inner add should be Integer, so the outer should too
    auto inner_add = make_add({make_var("a"), make_var("b")});
    auto outer_add = make_add({inner_add, make_var("c")});
    auto expr = wrap_expr(outer_add);

    EXPECT_TRUE(engine.query_integer(expr) == Tribool::True,
        "(a + b) + c is Integer when a, b, c are Integer");
}


int main() {
    test_property12_all_positive_operands();
    test_property12_all_negative_operands();
    test_property12_all_nonnegative_operands();
    test_property12_all_nonpositive_operands();
    test_property12_positive_implies_nonnegative_for_sum();
    test_property12_unknown_operand_yields_unknown();
    test_property12_mixed_signs_yield_unknown();
    test_property12_single_operand();
    test_property12_with_number_operands();
    test_property12_with_variables_and_numbers_mixed();
    test_property12_many_operands_uniform_sign();
    test_property12_empty_add_returns_unknown();

    test_property13_all_integer_operands();
    test_property13_all_real_operands();
    test_property13_integer_implies_real_for_sum();
    test_property13_mixed_integer_and_real();
    test_property13_unknown_domain_yields_unknown();
    test_property13_with_number_operands();
    test_property13_many_integer_operands();
    test_property13_real_with_numbers();
    test_property13_nested_addition_domain();

    return TEST_REPORT();
}
