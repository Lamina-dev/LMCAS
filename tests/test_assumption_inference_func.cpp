
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "symbolic_ast.hpp"
#include <memory>

using namespace lamina;

// Helper: create a SymbolicExpr wrapping a FunctionNode
static SymbolicExpr make_func_expr(FunctionNode::FuncType type,
                                   std::shared_ptr<const SymbolicNode> arg) {
    std::vector<std::shared_ptr<const SymbolicNode>> args = {std::move(arg)};
    return lamina::detail::expression_from_node(lamina::detail::make_node<FunctionNode>(type, std::move(args)));
}

// Helper: create a VariableNode
static std::shared_ptr<const SymbolicNode> var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

// Helper: create a NumberNode from an integer
static std::shared_ptr<const SymbolicNode> num(int v) {
    return lamina::detail::make_node<NumberNode>(BigInt(v));
}

// Helper: create a NumberNode from a double
static std::shared_ptr<const SymbolicNode> num_d(double v) {
    return lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(v));
}


void test_exp_real_arg_positive() {
    TEST_CASE("exp(x) is Positive when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "exp(x) is Positive when x is Real");
}

void test_exp_real_arg_real_domain() {
    TEST_CASE("exp(x) is Real when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "exp(x) is Real when x is Real");
}

void test_exp_real_arg_nonnegative() {
    TEST_CASE("exp(x) is NonNegative when x is Real (implied by Positive)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "exp(x) is NonNegative when x is Real");
}

void test_exp_real_arg_not_negative() {
    TEST_CASE("exp(x) is not Negative when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
                "exp(x) is not Negative when x is Real");
}

void test_exp_real_arg_nonzero() {
    TEST_CASE("exp(x) is NonZero when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
                "exp(x) is NonZero when x is Real");
}

void test_exp_integer_arg() {
    TEST_CASE("exp(n) is Positive and Real when n is Integer (Integer implies Real)");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("n"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "exp(n) is Positive when n is Integer");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "exp(n) is Real when n is Integer");
}

void test_exp_no_assumption() {
    TEST_CASE("exp(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    // No assumptions about x

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "exp(x) is Unknown for Positive when x has no assumptions");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "exp(x) is Unknown for Real when x has no assumptions");
}

void test_exp_numeric_arg() {
    TEST_CASE("exp(2) is Positive and Real (numeric argument is Real)");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // NumberNode 2 is Real (BigInt → Integer → Real)
    auto expr = make_func_expr(FunctionNode::FuncType::Exp, num(2));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "exp(2) is Positive");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "exp(2) is Real");
}


void test_sin_real_arg_real_domain() {
    TEST_CASE("sin(x) is Real when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sin(x) is Real when x is Real");
}

void test_sin_real_arg_sign_unknown() {
    TEST_CASE("sin(x) sign is Unknown when x is Real (can be +/-)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "sin(x) Positive is Unknown (sin can be negative)");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "sin(x) Negative is Unknown (sin can be positive)");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "sin(x) NonNegative is Unknown (sin can be negative)");
}

void test_sin_integer_arg() {
    TEST_CASE("sin(n) is Real when n is Integer (Integer implies Real)");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("n"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sin(n) is Real when n is Integer");
}

void test_sin_no_assumption() {
    TEST_CASE("sin(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "sin(x) Real is Unknown when x has no assumptions");
}

void test_cos_real_arg_real_domain() {
    TEST_CASE("cos(x) is Real when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "cos(x) is Real when x is Real");
}

void test_cos_real_arg_sign_unknown() {
    TEST_CASE("cos(x) sign is Unknown when x is Real (can be +/-)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "cos(x) Positive is Unknown (cos can be negative)");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "cos(x) Negative is Unknown (cos can be positive)");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "cos(x) NonNegative is Unknown (cos can be negative)");
}

void test_cos_integer_arg() {
    TEST_CASE("cos(n) is Real when n is Integer");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("n"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "cos(n) is Real when n is Integer");
}

void test_cos_no_assumption() {
    TEST_CASE("cos(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "cos(x) Real is Unknown when x has no assumptions");
}

void test_sin_numeric_arg() {
    TEST_CASE("sin(1) is Real (numeric argument is Real)");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sin, num(1));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sin(1) is Real");
}

void test_cos_numeric_arg() {
    TEST_CASE("cos(0) is Real (numeric argument is Real)");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Cos, num(0));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "cos(0) is Real");
}


// --- abs() tests ---

void test_abs_real_arg_nonnegative() {
    TEST_CASE("abs(x) is NonNegative when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "abs(x) is NonNegative when x is Real");
}

void test_abs_real_arg_real_domain() {
    TEST_CASE("abs(x) is Real when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "abs(x) is Real when x is Real");
}

void test_abs_real_arg_not_negative() {
    TEST_CASE("abs(x) is not Negative when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
                "abs(x) is not Negative when x is Real");
}

void test_abs_integer_arg() {
    TEST_CASE("abs(n) is NonNegative and Real when n is Integer");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("n"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "abs(n) is NonNegative when n is Integer");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "abs(n) is Real when n is Integer");
}

void test_abs_no_assumption() {
    TEST_CASE("abs(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "abs(x) NonNegative is Unknown when x has no assumptions");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "abs(x) Real is Unknown when x has no assumptions");
}

void test_abs_numeric_arg() {
    TEST_CASE("abs(-3) is NonNegative and Real");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, num(-3));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "abs(-3) is NonNegative");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "abs(-3) is Real");
}

// --- ln() tests ---

void test_ln_positive_arg_real() {
    TEST_CASE("ln(x) is Real when x is Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "ln(x) is Real when x is Positive");
}

void test_ln_positive_arg_sign_unknown() {
    TEST_CASE("ln(x) sign is Unknown when x is Positive (ln can be +/-)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    // ln(x) can be positive (x>1), negative (0<x<1), or zero (x=1)
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "ln(x) Positive is Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "ln(x) Negative is Unknown");
}

void test_ln_nonnegative_arg_not_sufficient() {
    TEST_CASE("ln(x) is Unknown when x is only NonNegative (not Positive)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    /// NonNegative 包含零；ln 的实数性还需要 StrictlyPositive 证明。
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "ln(x) Real is Unknown when x is only NonNegative");
}

void test_ln_no_assumption() {
    TEST_CASE("ln(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "ln(x) Real is Unknown when x has no assumptions");
}

void test_ln_numeric_positive_arg() {
    TEST_CASE("ln(2) is Real (numeric positive argument)");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, num(2));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "ln(2) is Real");
}

// --- sqrt() tests ---

void test_sqrt_nonneg_arg_nonnegative() {
    TEST_CASE("sqrt(x) is NonNegative when x is NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "sqrt(x) is NonNegative when x is NonNegative");
}

void test_sqrt_nonneg_arg_real() {
    TEST_CASE("sqrt(x) is Real when x is NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sqrt(x) is Real when x is NonNegative");
}

void test_sqrt_nonneg_arg_not_negative() {
    TEST_CASE("sqrt(x) is not Negative when x is NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
                "sqrt(x) is not Negative when x is NonNegative");
}

void test_sqrt_positive_arg() {
    TEST_CASE("sqrt(x) is NonNegative and Real when x is Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    // Positive implies NonNegative, so sqrt rule fires
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "sqrt(x) is NonNegative when x is Positive");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sqrt(x) is Real when x is Positive");
}

void test_sqrt_no_assumption() {
    TEST_CASE("sqrt(x) is Unknown when x has no assumptions");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "sqrt(x) NonNegative is Unknown when x has no assumptions");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "sqrt(x) Real is Unknown when x has no assumptions");
}

void test_sqrt_numeric_arg() {
    TEST_CASE("sqrt(4) is NonNegative and Real");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, num(4));

    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "sqrt(4) is NonNegative");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "sqrt(4) is Real");
}


void test_tan_real_arg_real_domain() {
    TEST_CASE("tan(x) is Real when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Tan, var("x"));

    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True,
                "tan(x) is Real when x is Real");
}

void test_tan_real_arg_sign_unknown() {
    TEST_CASE("tan(x) sign is Unknown when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Tan, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "tan(x) Positive is Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "tan(x) Negative is Unknown");
}


void test_unrecognized_function() {
    TEST_CASE("Unrecognized function returns Unknown for all properties");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // LambertW is not in the recognized list for property inference
    auto expr = make_func_expr(FunctionNode::FuncType::LambertW, var("x"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "LambertW(x) Positive is Unknown");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "LambertW(x) Negative is Unknown");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "LambertW(x) NonNegative is Unknown");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "LambertW(x) Real is Unknown");
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown,
                "LambertW(x) Integer is Unknown");
}


void test_insufficient_arg_properties() {
    TEST_CASE("Insufficient argument properties → Unknown");

    AssumptionContext ctx;
    // x has no assumptions — insufficient for any function rule

    InferenceEngine engine(ctx);

    auto exp_expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));
    EXPECT_TRUE(engine.query_positive_checked(exp_expr).value() == Tribool::Unknown,
                "exp(x) Unknown when x has no domain");

    auto sin_expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));
    EXPECT_TRUE(engine.query_real_checked(sin_expr).value() == Tribool::Unknown,
                "sin(x) Unknown when x has no domain");

    auto abs_expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));
    EXPECT_TRUE(engine.query_nonnegative_checked(abs_expr).value() == Tribool::Unknown,
                "abs(x) Unknown when x has no domain");

    auto ln_expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));
    EXPECT_TRUE(engine.query_real_checked(ln_expr).value() == Tribool::Unknown,
                "ln(x) Unknown when x has no sign");

    auto sqrt_expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));
    EXPECT_TRUE(engine.query_nonnegative_checked(sqrt_expr).value() == Tribool::Unknown,
                "sqrt(x) Unknown when x has no sign");
}

int main() {
    test_exp_real_arg_positive();
    test_exp_real_arg_real_domain();
    test_exp_real_arg_nonnegative();
    test_exp_real_arg_not_negative();
    test_exp_real_arg_nonzero();
    test_exp_integer_arg();
    test_exp_no_assumption();
    test_exp_numeric_arg();

    test_sin_real_arg_real_domain();
    test_sin_real_arg_sign_unknown();
    test_sin_integer_arg();
    test_sin_no_assumption();
    test_cos_real_arg_real_domain();
    test_cos_real_arg_sign_unknown();
    test_cos_integer_arg();
    test_cos_no_assumption();
    test_sin_numeric_arg();
    test_cos_numeric_arg();

    test_abs_real_arg_nonnegative();
    test_abs_real_arg_real_domain();
    test_abs_real_arg_not_negative();
    test_abs_integer_arg();
    test_abs_no_assumption();
    test_abs_numeric_arg();
    test_ln_positive_arg_real();
    test_ln_positive_arg_sign_unknown();
    test_ln_nonnegative_arg_not_sufficient();
    test_ln_no_assumption();
    test_ln_numeric_positive_arg();
    test_sqrt_nonneg_arg_nonnegative();
    test_sqrt_nonneg_arg_real();
    test_sqrt_nonneg_arg_not_negative();
    test_sqrt_positive_arg();
    test_sqrt_no_assumption();
    test_sqrt_numeric_arg();

    test_tan_real_arg_real_domain();
    test_tan_real_arg_sign_unknown();

    // Edge cases
    test_unrecognized_function();
    test_insufficient_arg_properties();

    return TEST_REPORT();
}
