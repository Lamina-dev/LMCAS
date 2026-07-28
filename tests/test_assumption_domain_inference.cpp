/**
 * @file test_assumption_domain_inference.cpp
 * @brief Property tests for InferenceEngine composite domain inference (Task 5.7).
 *
 * Properties tested:
 * - Property 2: Composite domain inference follows domain rules
 *
 * Validates: Requirements 2.1, 2.2, 2.3, 2.4, 2.5, 2.6
 *
 * Uses rapidcheck (header-only, vendored in tests/rapidcheck/) for
 * property-based testing with random input generation.
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;

// ============================================================
// Helpers: create AST nodes
// ============================================================

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

static std::shared_ptr<const SymbolicNode> make_function(
    FunctionNode::FuncType type,
    std::shared_ptr<const SymbolicNode> arg) {
    return lamina::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

/// Generate a random domain that is Integer or Real (for trig function arguments)
static Domain random_integer_or_real() {
    return rc::gen::boolean() ? Domain::Integer : Domain::Real;
}

/// Generate a random domain that is Rational or Real (for exp arguments)
static Domain random_rational_or_real() {
    int choice = rc::gen::inRange(0, 2);
    switch (choice) {
        case 0: return Domain::Rational;
        case 1: return Domain::Real;
        default: return Domain::Integer; // Integer implies Rational
    }
}

/// Generate a random Natural or PositiveInt domain (for power exponents)
static Domain random_natural_domain() {
    return rc::gen::boolean() ? Domain::Natural : Domain::PositiveInt;
}

// ============================================================
// Property 2: sin/cos/tan(Integer|Real) → Real
// **Validates: Requirements 2.1**
// ============================================================

static void test_property2_trig_integer_or_real_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: sin/cos/tan(Integer|Real) → Real");

    rc::check("For any trig function with Integer or Real argument, result is Real", []() {
        // Pick a random trig function
        std::vector<FunctionNode::FuncType> trig_funcs = {
            FunctionNode::FuncType::Sin,
            FunctionNode::FuncType::Cos,
            FunctionNode::FuncType::Tan
        };
        auto func_type = rc::gen::elementOf(trig_funcs);

        // Pick a random domain for the argument
        Domain arg_domain = random_integer_or_real();

        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, arg_domain);
        InferenceEngine engine(ctx);

        auto func_node = make_function(func_type, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: exp(Rational|Real) → Real
// **Validates: Requirements 2.2**
// ============================================================

static void test_property2_exp_rational_or_real_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: exp(Rational|Real) → Real");

    rc::check("For exp with Rational or Real argument, result is Real", []() {
        Domain arg_domain = random_rational_or_real();
        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, arg_domain);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Exp, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: ln(Integer) → Real
// **Validates: Requirements 2.3**
// ============================================================

static void test_property2_ln_integer_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: ln(Integer) → Real");

    rc::check("For ln with Integer argument, result is Real", []() {
        std::string var_name = "n_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Ln, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: sqrt(NonNeg Real) → Real
// **Validates: Requirements 2.4**
// ============================================================

static void test_property2_sqrt_nonneg_real_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: sqrt(NonNeg Real) → Real");

    rc::check("For sqrt with non-negative Real argument, result is Real", []() {
        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, Domain::Real);
        ctx.assume_sign(var_name, Sign::NonNegative);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Sqrt, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: Integer^Natural → Integer
// **Validates: Requirements 2.5**
// ============================================================

static void test_property2_integer_power_natural_gives_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: Integer^Natural → Integer");

    rc::check("For Integer base raised to a positive integer exponent, result is Integer", []() {
        std::string base_name = "b_" + std::to_string(rc::gen::inRange(0, 999));
        int exponent = rc::gen::inRange(1, 10); // Positive integer exponent

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(exponent));
        auto expr = wrap_expr(pow_node);

        RC_ASSERT(engine.query_integer(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: Integer^0 → Integer (x^0 = 1)
// **Validates: Requirements 2.5 (edge case)**
// ============================================================

static void test_property2_integer_power_zero_gives_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: Integer^0 → Integer");

    rc::check("For Integer base raised to 0, result is Integer (x^0 = 1)", []() {
        std::string base_name = "b_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(0));
        auto expr = wrap_expr(pow_node);

        RC_ASSERT(engine.query_integer(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: Rational^Integer → Rational (via Real inference)
// **Validates: Requirements 2.6**
// ============================================================

static void test_property2_rational_power_integer_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: Rational^Integer → Real");

    rc::check("For Rational base raised to integer exponent, result is Real", []() {
        std::string base_name = "r_" + std::to_string(rc::gen::inRange(0, 999));
        int exponent = rc::gen::inRange(1, 10);

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Rational);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(exponent));
        auto expr = wrap_expr(pow_node);

        // Rational implies Real, and Real^Integer → Real
        RC_ASSERT(engine.query_real(expr) == Tribool::True);
    });
}

// ============================================================
// Property 2: Unknown domain argument → Unknown result
// ============================================================

static void test_property2_unknown_domain_gives_unknown() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: Unknown domain argument → Unknown");

    rc::check("For functions with unknown domain argument, result domain is Unknown", []() {
        std::vector<FunctionNode::FuncType> funcs = {
            FunctionNode::FuncType::Sin,
            FunctionNode::FuncType::Cos,
            FunctionNode::FuncType::Tan,
            FunctionNode::FuncType::Exp
        };
        auto func_type = rc::gen::elementOf(funcs);
        std::string var_name = "u_" + std::to_string(rc::gen::inRange(0, 999));

        // No domain declared for variable (defaults to Complex)
        AssumptionContext ctx;
        InferenceEngine engine(ctx);

        auto func_node = make_function(func_type, make_var(var_name));
        auto expr = wrap_expr(func_node);

        // Without Real/Integer domain on argument, can't infer Real result
        RC_ASSERT(engine.query_real(expr) == Tribool::Unknown);
    });
}

// ============================================================
// Property 2: Comprehensive domain inference for all trig functions
// ============================================================

static void test_property2_all_trig_with_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: All trig(Integer) → Real");

    // sin(Integer) → Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sin, make_var("n")));
        EXPECT_TRUE(engine.query_real(expr) == Tribool::True, "sin(Integer) → Real");
    }
    // cos(Integer) → Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Cos, make_var("n")));
        EXPECT_TRUE(engine.query_real(expr) == Tribool::True, "cos(Integer) → Real");
    }
    // tan(Integer) → Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Tan, make_var("n")));
        EXPECT_TRUE(engine.query_real(expr) == Tribool::True, "tan(Integer) → Real");
    }
}

static void test_property2_exp_with_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: exp(Integer) → Real");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("n")));
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True, "exp(Integer) → Real");
}

static void test_property2_ln_positive_gives_real() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: ln(Positive) → Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Ln, make_var("x")));
    EXPECT_TRUE(engine.query_real(expr) == Tribool::True, "ln(Positive Real) → Real");
}

static void test_property2_sqrt_without_nonneg_unknown() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: sqrt without NonNeg → Unknown");

    // sqrt of a Real variable without NonNegative sign → Unknown
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    // No sign declared — could be negative
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sqrt, make_var("x")));
    EXPECT_TRUE(engine.query_real(expr) == Tribool::Unknown,
        "sqrt(Real without NonNeg) → Unknown");
}

static void test_property2_power_negative_exponent_not_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 2: Integer^(-1) not necessarily Integer");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);
    InferenceEngine engine(ctx);

    // n^(-1) = 1/n — not necessarily integer
    auto pow_node = make_power(make_var("n"), make_number(-1));
    auto expr = wrap_expr(pow_node);

    // Should NOT be able to infer Integer (e.g., 2^(-1) = 0.5)
    EXPECT_TRUE(engine.query_integer(expr) == Tribool::Unknown,
        "Integer^(-1) → Integer is Unknown (not guaranteed)");
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 2: Composite domain inference
    test_property2_trig_integer_or_real_gives_real();
    test_property2_exp_rational_or_real_gives_real();
    test_property2_ln_integer_gives_real();
    test_property2_sqrt_nonneg_real_gives_real();
    test_property2_integer_power_natural_gives_integer();
    test_property2_integer_power_zero_gives_integer();
    test_property2_rational_power_integer_gives_real();
    test_property2_unknown_domain_gives_unknown();

    // Additional unit-style tests for completeness
    test_property2_all_trig_with_integer();
    test_property2_exp_with_integer();
    test_property2_ln_positive_gives_real();
    test_property2_sqrt_without_nonneg_unknown();
    test_property2_power_negative_exponent_not_integer();

    return TEST_REPORT();
}
