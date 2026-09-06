
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base,
    std::shared_ptr<const SymbolicNode> exp) {
    return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_function(
    FunctionNode::FuncType type,
    std::shared_ptr<const SymbolicNode> arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = LMCAS::detail::expression_from_node(std::move(node));
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


static void test_trig_integer_or_real_gives_real() {
    TEST_CASE("sin/cos/tan(Integer|Real) → Real");

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

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    });
}


static void test_exp_rational_or_real_gives_real() {
    TEST_CASE("exp(Rational|Real) → Real");

    rc::check("For exp with Rational or Real argument, result is Real", []() {
        Domain arg_domain = random_rational_or_real();
        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, arg_domain);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Exp, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    });
}


static void test_ln_integer_gives_real() {
    TEST_CASE("ln(Integer) → Real");

    rc::check("For ln with Integer argument, result is Real", []() {
        std::string var_name = "n_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Ln, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    });
}


static void test_sqrt_nonneg_real_gives_real() {
    TEST_CASE("sqrt(NonNeg Real) → Real");

    rc::check("For sqrt with non-negative Real argument, result is Real", []() {
        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(var_name, Domain::Real);
        ctx.assume_sign(var_name, Sign::NonNegative);
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Sqrt, make_var(var_name));
        auto expr = wrap_expr(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    });
}


static void test_integer_power_natural_gives_integer() {
    TEST_CASE("Integer^Natural → Integer");

    rc::check("For Integer base raised to a positive integer exponent, result is Integer", []() {
        std::string base_name = "b_" + std::to_string(rc::gen::inRange(0, 999));
        int exponent = rc::gen::inRange(1, 10); // Positive integer exponent

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(exponent));
        auto expr = wrap_expr(pow_node);

        RC_ASSERT(engine.query_integer_checked(expr).value() == Tribool::True);
    });
}


static void test_integer_power_zero_gives_integer() {
    TEST_CASE("Integer^0 → Integer");

    rc::check("For Integer base raised to 0, result is Integer (x^0 = 1)", []() {
        std::string base_name = "b_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Integer);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(0));
        auto expr = wrap_expr(pow_node);

        RC_ASSERT(engine.query_integer_checked(expr).value() == Tribool::True);
    });
}


static void test_rational_power_integer_gives_real() {
    TEST_CASE("Rational^Integer → Real");

    rc::check("For Rational base raised to integer exponent, result is Real", []() {
        std::string base_name = "r_" + std::to_string(rc::gen::inRange(0, 999));
        int exponent = rc::gen::inRange(1, 10);

        AssumptionContext ctx;
        ctx.assume_domain(base_name, Domain::Rational);
        InferenceEngine engine(ctx);

        auto pow_node = make_power(make_var(base_name), make_number(exponent));
        auto expr = wrap_expr(pow_node);

        // Rational implies Real, and Real^Integer -> Real
        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    });
}


static void test_unknown_domain_gives_unknown() {
    TEST_CASE("Unknown domain argument → Unknown");

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

        /// 参数域缺少 Real/Integer 证明时,函数值实数性保持 Unknown.
        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::Unknown);
    });
}


static void test_all_trig_with_integer() {
    TEST_CASE("All trig(Integer) → Real");

    // sin(Integer) -> Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sin, make_var("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, "sin(Integer) → Real");
    }
    // cos(Integer) -> Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Cos, make_var("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, "cos(Integer) → Real");
    }
    // tan(Integer) -> Real
    {
        AssumptionContext ctx;
        ctx.assume_domain("n", Domain::Integer);
        InferenceEngine engine(ctx);
        auto expr = wrap_expr(make_function(FunctionNode::FuncType::Tan, make_var("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, "tan(Integer) → Real");
    }
}

static void test_exp_with_integer() {
    TEST_CASE("exp(Integer) → Real");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Exp, make_var("n")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, "exp(Integer) → Real");
}

static void test_ln_positive_gives_real() {
    TEST_CASE("ln(Positive) → Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Ln, make_var("x")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, "ln(Positive Real) → Real");
}

static void test_sqrt_without_nonneg_unknown() {
    TEST_CASE("sqrt without NonNeg → Unknown");

    /// sqrt 的实数性需要 Real 与 NonNegative 共同证明;当前仅声明 Real.
    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    /// 符号属性保持未声明状态.
    InferenceEngine engine(ctx);
    auto expr = wrap_expr(make_function(FunctionNode::FuncType::Sqrt, make_var("x")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
        "sqrt(Real without NonNeg) → Unknown");
}

static void test_power_negative_exponent_not_integer() {
    TEST_CASE("Integer^(-1) not necessarily Integer");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);
    InferenceEngine engine(ctx);

    // n^(-1) = 1/n - not necessarily integer
    auto pow_node = make_power(make_var("n"), make_number(-1));
    auto expr = wrap_expr(pow_node);

    // Should NOT be able to infer Integer (e.g., 2^(-1) = 0.5)
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown,
        "Integer^(-1) → Integer is Unknown (not guaranteed)");
}


int main() {
    test_trig_integer_or_real_gives_real();
    test_exp_rational_or_real_gives_real();
    test_ln_integer_gives_real();
    test_sqrt_nonneg_real_gives_real();
    test_integer_power_natural_gives_integer();
    test_integer_power_zero_gives_integer();
    test_rational_power_integer_gives_real();
    test_unknown_domain_gives_unknown();

    // Additional unit-style tests for completeness
    test_all_trig_with_integer();
    test_exp_with_integer();
    test_ln_positive_gives_real();
    test_sqrt_without_nonneg_unknown();
    test_power_negative_exponent_not_integer();

    return TEST_REPORT();
}
