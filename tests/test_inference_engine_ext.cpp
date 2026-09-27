
#include "test_assumption_inference_mul_support.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <type_traits>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exp) {
    return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_function(
    FunctionNode::FuncType type, std::shared_ptr<const SymbolicNode> arg) {
    return LMCAS::detail::make_node<FunctionNode>(type,
                                                  std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static std::shared_ptr<const SymbolicNode> make_add(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return LMCAS::detail::make_node<AddNode>(std::move(ops));
}

/// Create a division expression: num / den represented as MultiplyNode([num, PowerNode(den, -1)])
static std::shared_ptr<const SymbolicNode> make_division(
    std::shared_ptr<const SymbolicNode> num, std::shared_ptr<const SymbolicNode> den) {
    auto den_inv = make_power(std::move(den), test_integer_node(-1));
    return make_multiply({std::move(num), std::move(den_inv)});
}

/// Create a subtraction expression: a - b represented as AddNode([a, MultiplyNode([-1, b])])
static std::shared_ptr<const SymbolicNode> make_subtraction(
    std::shared_ptr<const SymbolicNode> a, std::shared_ptr<const SymbolicNode> b) {
    auto neg_b = make_multiply({test_integer_node(-1), std::move(b)});
    return make_add({std::move(a), std::move(neg_b)});
}

TEST(InferenceEngineExt, DivisionPositiveOverPositive) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "pos / pos is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "pos / pos is not Negative";
}

TEST(InferenceEngineExt, DivisionNegativeOverNegative) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "neg / neg is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "neg / neg is not Negative";
}

TEST(InferenceEngineExt, DivisionPositiveOverNegative) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "pos / neg is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "pos / neg is Negative";
}

TEST(InferenceEngineExt, DivisionNegativeOverPositive) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "neg / pos is not Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "neg / pos is Negative";
}

TEST(InferenceEngineExt, DivisionUnknownDenominator) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    // b has no sign declared
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::Unknown)) << "pos / unknown: Positive is Unknown";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::Unknown)) << "pos / unknown: Negative is Unknown";
}

TEST(InferenceEngineExt, DivisionZeroDenominator) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Zero).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_division(test_variable_node("a"), test_variable_node("b")));

    auto positive = engine.query_positive_checked(expr);
    auto negative = engine.query_negative_checked(expr);
    EXPECT_TRUE((!positive && positive.error().code == CasErrc::DomainError)) << "undefined division has no positive value";
    EXPECT_TRUE((!negative && negative.error().code == CasErrc::DomainError)) << "undefined division has no negative value";
}

TEST(InferenceEngineExt, SubtractionPositiveMinusNegative) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_subtraction(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::True)) << "pos - neg is Positive";
    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::False)) << "pos - neg is not Negative";
}

TEST(InferenceEngineExt, SubtractionNegativeMinusPositive) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_sign("a", Sign::Negative).has_value())) << "sign assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("b", Sign::Positive).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_subtraction(test_variable_node("a"), test_variable_node("b")));

    EXPECT_TRUE((engine.query_negative_checked(expr).value() == Tribool::True)) << "neg - pos is Negative";
    EXPECT_TRUE((engine.query_positive_checked(expr).value() == Tribool::False)) << "neg - pos is not Positive";
}

TEST(InferenceEngineExt, DomainSinInteger) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Integer).has_value())) << "domain assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Sin, test_variable_node("x")));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sin(integer) is Real";
}

TEST(InferenceEngineExt, DomainExpRational) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Rational).has_value())) << "domain assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Exp, test_variable_node("x")));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "exp(rational) is Real";
}

TEST(InferenceEngineExt, DomainLnInteger) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Integer).has_value())) << "domain assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Ln, test_variable_node("x")));

    auto unrestricted = engine.query_real_checked(expr);
    EXPECT_TRUE((unrestricted && unrestricted.value() == Tribool::Unknown)) << "ln(integer) is not proved real without a positive argument";
    EXPECT_TRUE((ctx.assume_sign("x", Sign::Positive).has_value())) << "positive integer declared";
    auto positive = engine.query_real_checked(expr);
    EXPECT_TRUE((positive && positive.value() == Tribool::True)) << "a positive integer argument proves the logarithm real";
}

TEST(InferenceEngineExt, DomainSqrtNonnegReal) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Real).has_value())) << "domain assumption setup succeeds";
    EXPECT_TRUE((ctx.assume_sign("x", Sign::NonNegative).has_value())) << "sign assumption setup succeeds";
    InferenceEngine engine(ctx);

    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Sqrt, test_variable_node("x")));

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "sqrt(nonneg real) is Real";
}

TEST(InferenceEngineExt, DomainIntegerPowNatural) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Integer).has_value())) << "domain assumption setup succeeds";
    InferenceEngine engine(ctx);

    // x^3 where x is Integer and 3 is a positive integer
    auto expr = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(3)));

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "integer^3 is Integer";
}

TEST(InferenceEngineExt, DomainRationalPowInteger) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.assume_domain("x", Domain::Rational).has_value())) << "domain assumption setup succeeds";
    InferenceEngine engine(ctx);

    // x^2 where x is Rational and 2 is an integer exponent
    auto expr = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(2)));

    // Rational  subset  Real, so Rational base + integer exponent -> Real
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "rational^2 is Real";
}
