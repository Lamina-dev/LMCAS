
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>

using namespace lamina;


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exp) {
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_function(
    FunctionNode::FuncType type, std::shared_ptr<const SymbolicNode> arg) {
    return lamina::detail::make_node<FunctionNode>(type,
        std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}


static void test_x_squared_nonnegative_when_real() {
    TEST_CASE("Propagation: x² non-negative when x is Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // Build x² = PowerNode(x, 2)
    auto x_squared = wrap_expr(make_power(make_var("x"), make_number(2)));

    EXPECT_TRUE(engine.query_nonnegative_checked(x_squared).value() == Tribool::True,
                "x² is NonNegative when x is Real");
}

static void test_x_squared_integer_when_integer() {
    TEST_CASE("Propagation: x² Integer when x is Integer");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);
    InferenceEngine engine(ctx);

    // Build x² = PowerNode(x, 2)
    auto x_squared = wrap_expr(make_power(make_var("x"), make_number(2)));

    EXPECT_TRUE(engine.query_integer_checked(x_squared).value() == Tribool::True,
                "x² is Integer when x is Integer");
}

static void test_abs_positive_when_x_positive() {
    TEST_CASE("Propagation: |x| positive when x is Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = wrap_expr(make_function(FunctionNode::FuncType::Abs, make_var("x")));

    EXPECT_TRUE(engine.query_positive_checked(abs_x).value() == Tribool::True,
                "|x| is Positive when x is Positive");
}

static void test_abs_positive_when_x_negative() {
    TEST_CASE("Propagation: |x| positive when x is Negative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = wrap_expr(make_function(FunctionNode::FuncType::Abs, make_var("x")));

    EXPECT_TRUE(engine.query_positive_checked(abs_x).value() == Tribool::True,
                "|x| is Positive when x is Negative");
}

static void test_abs_positive_when_x_nonzero() {
    TEST_CASE("Propagation: |x| positive when x is NonZero");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonZero);
    ctx.assume_domain("x", Domain::Real);
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = wrap_expr(make_function(FunctionNode::FuncType::Abs, make_var("x")));

    EXPECT_TRUE(engine.query_positive_checked(abs_x).value() == Tribool::True,
                "|x| is Positive when x is NonZero");
}


static void test_diagnostic_transcendental_then_integer() {
    TEST_CASE("Diagnostics: Transcendental + Integer contradiction");

    AssumptionContext ctx;
    ctx.current_properties().declare_transcendental("x");

    auto failure_119 = ctx.assume_domain("x", Domain::Integer);
    EXPECT_TRUE(!failure_119.has_value(), "Transcendental + Integer returns InvalidArgument");
    EXPECT_TRUE(failure_119.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    const std::string& msg = failure_119.error().message;
    // Message should contain the symbol name and relevant domain info
    EXPECT_CONTAINS(msg, {"x", "Integer"},
                    "Result error message contains 'x' and 'Integer'");
    // Should also mention Transcendental or Real (the existing constraint)
    bool has_transcendental_or_real =
        (msg.find("Transcendental") != std::string::npos) ||
        (msg.find("Real") != std::string::npos) ||
        (msg.find("transcendental") != std::string::npos);
    EXPECT_TRUE(has_transcendental_or_real,
                "Result error message mentions Transcendental or Real");
}

static void test_diagnostic_positive_then_negative() {
    TEST_CASE("Diagnostics: Positive + Negative contradiction");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    auto failure_147 = ctx.assume_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_147.has_value(), "Positive + Negative returns InvalidArgument");
    EXPECT_TRUE(failure_147.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    const std::string& msg = failure_147.error().message;
    // The message should contain the symbol name and mention Positive.
    // The system may report the implied sign (NonPositive) rather than the
    // literal "Negative" since Negative implies NonPositive which contradicts Positive.
    EXPECT_CONTAINS(msg, {"x", "Positive"},
                    "Result error message contains 'x' and 'Positive'");
    bool has_negative_or_nonpositive =
        (msg.find("Negative") != std::string::npos) ||
        (msg.find("NonPositive") != std::string::npos);
    EXPECT_TRUE(has_negative_or_nonpositive,
                "Result error message mentions Negative or NonPositive");
}

static void test_diagnostic_natural_then_negative() {
    TEST_CASE("Diagnostics: Natural + Negative contradiction");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Natural);

    auto failure_175 = ctx.assume_sign("x", Sign::Negative);
    EXPECT_TRUE(!failure_175.has_value(), "Natural + Negative returns InvalidArgument");
    EXPECT_TRUE(failure_175.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    const std::string& msg = failure_175.error().message;
    EXPECT_CONTAINS(msg, {"Natural", "Negative"},
                    "Result error message contains 'Natural' and 'Negative'");
}

static void test_diagnostic_positiveint_then_zero() {
    TEST_CASE("Diagnostics: PositiveInt + Zero contradiction");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::PositiveInt);

    auto failure_195 = ctx.assume_sign("x", Sign::Zero);
    EXPECT_TRUE(!failure_195.has_value(), "PositiveInt + Zero returns InvalidArgument");
    EXPECT_TRUE(failure_195.error().code == CasErrc::InvalidArgument, "failure reports InvalidArgument");
    const std::string& msg = failure_195.error().message;
    // The message should mention PositiveInt (the domain) and the conflicting sign.
    // The system may report the implied sign (NonPositive) rather than "Zero"
    // since Zero implies NonPositive which contradicts PositiveInt's implied Positive.
    EXPECT_CONTAINS(msg, {"PositiveInt"},
                    "Result error message contains 'PositiveInt'");
    bool has_zero_or_nonpositive =
        (msg.find("Zero") != std::string::npos) ||
        (msg.find("NonPositive") != std::string::npos);
    EXPECT_TRUE(has_zero_or_nonpositive,
                "Result error message mentions Zero or NonPositive");
}


int main() {
    test_x_squared_nonnegative_when_real();
    test_x_squared_integer_when_integer();
    test_abs_positive_when_x_positive();
    test_abs_positive_when_x_negative();
    test_abs_positive_when_x_nonzero();

    test_diagnostic_transcendental_then_integer();
    test_diagnostic_positive_then_negative();
    test_diagnostic_natural_then_negative();
    test_diagnostic_positiveint_then_zero();

    return TEST_REPORT();
}
