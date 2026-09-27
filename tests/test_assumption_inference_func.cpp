
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>

using namespace LMCAS;

// Helper: create a SymbolicExpr wrapping a FunctionNode
static SymbolicExpr make_func_expr(FunctionNode::FuncType type,
                                   std::shared_ptr<const SymbolicNode> arg) {
    std::vector<std::shared_ptr<const SymbolicNode>> args = {std::move(arg)};
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<FunctionNode>(type, std::move(args)));
}

// Helper: create a VariableNode
static std::shared_ptr<const SymbolicNode> var(const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

// Helper: create a NumberNode from an integer
static std::shared_ptr<const SymbolicNode> num(int v) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(v));
}

// Helper: create a NumberNode from a double

TEST(AssumptionInferenceFunc, ExpRealArgPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "exp(x) is Positive when x is Real";
}

TEST(AssumptionInferenceFunc, ExpRealArgRealDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "exp(x) is Real when x is Real";
}

TEST(AssumptionInferenceFunc, ExpRealArgNonnegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "exp(x) is NonNegative when x is Real";
}

TEST(AssumptionInferenceFunc, ExpRealArgNotNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "exp(x) is not Negative when x is Real";
}

TEST(AssumptionInferenceFunc, ExpRealArgNonzero) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_nonzero_checked(expr).value() == Tribool::True)) << "exp(x) is NonZero when x is Real";
}

TEST(AssumptionInferenceFunc, ExpIntegerArg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("n"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "exp(n) is Positive when n is Integer";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "exp(n) is Real when n is Integer";
}

TEST(AssumptionInferenceFunc, ExpNoAssumption) {
    AssumptionContext ctx;
    // No assumptions about x

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "exp(x) is Unknown for Positive when x has no assumptions";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "exp(x) is Unknown for Real when x has no assumptions";
}

TEST(AssumptionInferenceFunc, ExpNumericArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // NumberNode 2 is Real (BigInt -> Integer -> Real)
    auto expr = make_func_expr(FunctionNode::FuncType::Exp, num(2));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "exp(2) is Positive";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "exp(2) is Real";
}

// --- abs() tests ---

TEST(AssumptionInferenceFunc, AbsRealArgNonnegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "abs(x) is NonNegative when x is Real";
}

TEST(AssumptionInferenceFunc, AbsRealArgRealDomain) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "abs(x) is Real when x is Real";
}

TEST(AssumptionInferenceFunc, AbsRealArgNotNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "abs(x) is not Negative when x is Real";
}

TEST(AssumptionInferenceFunc, AbsIntegerArg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("n"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "abs(n) is NonNegative when n is Integer";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "abs(n) is Real when n is Integer";
}

TEST(AssumptionInferenceFunc, AbsNoAssumption) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "abs(x) NonNegative is Unknown when x has no assumptions";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "abs(x) Real is Unknown when x has no assumptions";
}

TEST(AssumptionInferenceFunc, AbsNumericArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Abs, num(-3));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "abs(-3) is NonNegative";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "abs(-3) is Real";
}

// --- ln() tests ---

TEST(AssumptionInferenceFunc, LnPositiveArgReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "ln(x) is Real when x is Positive";
}

TEST(AssumptionInferenceFunc, LnPositiveArgSignUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    // ln(x) can be positive (x>1), negative (0<x<1), or zero (x=1)
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "ln(x) Positive is Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "ln(x) Negative is Unknown";
}

TEST(AssumptionInferenceFunc, LnNonnegativeArgNotSufficient) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    /// NonNegative 包含零;ln 的实数性还需要 StrictlyPositive 证明.
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "ln(x) Real is Unknown when x is only NonNegative";
}

TEST(AssumptionInferenceFunc, LnNoAssumption) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "ln(x) Real is Unknown when x has no assumptions";
}

TEST(AssumptionInferenceFunc, LnNumericPositiveArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Ln, num(2));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "ln(2) is Real";
}

// --- sqrt() tests ---

TEST(AssumptionInferenceFunc, SqrtNonnegArgNonnegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "sqrt(x) is NonNegative when x is NonNegative";
}

TEST(AssumptionInferenceFunc, SqrtNonnegArgReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sqrt(x) is Real when x is NonNegative";
}

TEST(AssumptionInferenceFunc, SqrtNonnegArgNotNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "sqrt(x) is not Negative when x is NonNegative";
}

TEST(AssumptionInferenceFunc, SqrtPositiveArg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    // Positive implies NonNegative, so sqrt rule fires
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "sqrt(x) is NonNegative when x is Positive";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sqrt(x) is Real when x is Positive";
}

TEST(AssumptionInferenceFunc, SqrtNoAssumption) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "sqrt(x) NonNegative is Unknown when x has no assumptions";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "sqrt(x) Real is Unknown when x has no assumptions";
}

TEST(AssumptionInferenceFunc, SqrtNumericArg) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto expr = make_func_expr(FunctionNode::FuncType::Sqrt, num(4));

    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::True)) << "sqrt(4) is NonNegative";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sqrt(4) is Real";
}

TEST(AssumptionInferenceFunc, UnrecognizedFunction) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    // LambertW is not in the recognized list for property inference
    auto expr = make_func_expr(FunctionNode::FuncType::LambertW, var("x"));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "LambertW(x) Positive is Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "LambertW(x) Negative is Unknown";
    EXPECT_TRUE((engine.query_nonnegative_checked(expr).value() == Tribool::Unknown)) << "LambertW(x) NonNegative is Unknown";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "LambertW(x) Real is Unknown";
    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "LambertW(x) Integer is Unknown";
}

TEST(AssumptionInferenceFunc, InsufficientArgProperties) {
    AssumptionContext ctx;
    // x has no assumptions - insufficient for any function rule

    InferenceEngine engine(ctx);

    auto exp_expr = make_func_expr(FunctionNode::FuncType::Exp, var("x"));
    EXPECT_TRUE((engine.query_positive_checked(exp_expr).value() == Tribool::Unknown)) << "exp(x) Unknown when x has no domain";

    auto sin_expr = make_func_expr(FunctionNode::FuncType::Sin, var("x"));
    EXPECT_TRUE((engine.query_real_checked(sin_expr).value() == Tribool::Unknown)) << "sin(x) Unknown when x has no domain";

    auto abs_expr = make_func_expr(FunctionNode::FuncType::Abs, var("x"));
    EXPECT_TRUE((engine.query_nonnegative_checked(abs_expr).value() == Tribool::Unknown)) << "abs(x) Unknown when x has no domain";

    auto ln_expr = make_func_expr(FunctionNode::FuncType::Ln, var("x"));
    EXPECT_TRUE((engine.query_real_checked(ln_expr).value() == Tribool::Unknown)) << "ln(x) Unknown when x has no sign";

    auto sqrt_expr = make_func_expr(FunctionNode::FuncType::Sqrt, var("x"));
    EXPECT_TRUE((engine.query_nonnegative_checked(sqrt_expr).value() == Tribool::Unknown)) << "sqrt(x) Unknown when x has no sign";
}

static void expect_exact_logarithm_signs(InferenceEngine &engine) {
    const Sign targets[] = {Sign::Positive, Sign::Negative, Sign::Zero,
                            Sign::NonNegative, Sign::NonPositive, Sign::NonZero};
    const struct {
        Rational value;
        bool signs[6];
    } cases[] = {
        {Rational(1), {false, false, true, true, true, false}},
        {Rational(1, 10), {false, true, false, false, true, true}},
        {Rational(2), {true, false, false, true, false, true}}};
    for (const auto &item : cases) {
        auto expression = make_func_expr(FunctionNode::FuncType::Ln,
                                         detail::make_node<NumberNode>(item.value));
        for (std::size_t index = 0; index < 6; ++index) {
            const auto target = targets[index];
            const bool expected = item.signs[index];
            auto result = test_sign_query(engine, expression, target);
            EXPECT_TRUE((result && result.value() == (expected ? Tribool::True : Tribool::False))) << "exact logarithm sign";
        }
    }
}

static void expect_logarithm_domain_constraints(InferenceEngine &engine) {
    for (int value : {0, -1}) {
        auto expression = make_func_expr(FunctionNode::FuncType::Ln, num(value));
        auto real = engine.query_real_checked(expression);
        auto positive = engine.query_positive_checked(expression);
        EXPECT_TRUE((!real && real.error().code == CasErrc::DomainError)) << "invalid real logarithm";
        EXPECT_TRUE((!positive && positive.error().code == CasErrc::DomainError)) << "invalid sign";
        auto nested = make_func_expr(FunctionNode::FuncType::Exp, detail::node(expression));
        auto result = engine.query_positive_checked(nested);
        EXPECT_TRUE((!result && result.error().code == CasErrc::DomainError)) << "exp cannot erase invalid ln";
    }
    for (int value : {0, -2}) {
        auto expression = make_func_expr(FunctionNode::FuncType::Exp, num(value));
        auto result = engine.query_positive_checked(expression);
        EXPECT_TRUE((result && result.value() == Tribool::True)) << "exp of every real constant is positive";
    }
}

static void expect_relational_logarithm_signs() {
    for (RelationOp op : {RelationOp::LT, RelationOp::EQ, RelationOp::GT}) {
        AssumptionContext related;
        EXPECT_TRUE((related.assume_sign("x", Sign::Positive).has_value())) << "positive argument";
        auto condition = detail::expression_from_node(detail::make_node<RelationalNode>(var("x"), num(1), op));
        EXPECT_TRUE((related.assume(condition).has_value())) << "relation to one";
        InferenceEngine relational_engine(related);
        auto expression = make_func_expr(FunctionNode::FuncType::Ln, var("x"));
        auto result = test_sign_query(relational_engine, expression, op == RelationOp::LT ? Sign::Negative : op == RelationOp::GT ? Sign::Positive
                                                                                                                                  : Sign::Zero);
        EXPECT_TRUE((result && result.value() == Tribool::True)) << "exact relation determines logarithm sign";
    }
}

TEST(AssumptionInferenceFunc, ExactLogarithmSignAndDomain) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);
    expect_exact_logarithm_signs(engine);
    expect_logarithm_domain_constraints(engine);
    expect_relational_logarithm_signs();
}
