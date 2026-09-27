
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>

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

TEST(AssumptionPropagationUnit, XSquaredNonnegativeWhenReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    // Build x^2 = PowerNode(x, 2)
    auto x_squared = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(2)));

    EXPECT_TRUE((engine.query_nonnegative_checked(x_squared).value() == Tribool::True)) << "x² is NonNegative when x is Real";
}

TEST(AssumptionPropagationUnit, XSquaredIntegerWhenInteger) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    // Build x^2 = PowerNode(x, 2)
    auto x_squared = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(2)));

    EXPECT_TRUE((engine.query_integer_checked(x_squared).value() == Tribool::True)) << "x² is Integer when x is Integer";
}

TEST(AssumptionPropagationUnit, AbsPositiveWhenXPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = test_expression_from_node(make_function(FunctionNode::FuncType::Abs, test_variable_node("x")));

    EXPECT_TRUE((engine.query_positive_checked(abs_x).value() == Tribool::True)) << "|x| is Positive when x is Positive";
}

TEST(AssumptionPropagationUnit, AbsPositiveWhenXNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = test_expression_from_node(make_function(FunctionNode::FuncType::Abs, test_variable_node("x")));

    EXPECT_TRUE((engine.query_positive_checked(abs_x).value() == Tribool::True)) << "|x| is Positive when x is Negative";
}

TEST(AssumptionPropagationUnit, AbsPositiveWhenXNonzero) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonZero).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    // Build |x| = FunctionNode::Abs(x)
    auto abs_x = test_expression_from_node(make_function(FunctionNode::FuncType::Abs, test_variable_node("x")));

    EXPECT_TRUE((engine.query_positive_checked(abs_x).value() == Tribool::True)) << "|x| is Positive when x is NonZero";
}

TEST(AssumptionPropagationUnit, DiagnosticTranscendentalThenInteger) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.current_properties().declare_transcendental("x").has_value())) << "transcendental declaration succeeds";

    auto failure_119 = ctx.assume_domain("x", Domain::Integer);
    EXPECT_TRUE((!failure_119.has_value())) << "Transcendental + Integer returns InvalidArgument";
    EXPECT_TRUE((failure_119.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
    const std::string &msg = failure_119.error().message;
    // Message should contain the symbol name and relevant domain info
    {
        const std::string actual_text = (msg);
        for (const auto &token : std::vector<std::string>{"x", "Integer"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Result error message contains 'x' and 'Integer'" << ": missing " << token << " in " << actual_text;
        }
    }
    // Should also mention Transcendental or Real (the existing constraint)
    bool has_transcendental_or_real =
        (msg.find("Transcendental") != std::string::npos) ||
        (msg.find("Real") != std::string::npos) ||
        (msg.find("transcendental") != std::string::npos);
    EXPECT_TRUE((has_transcendental_or_real)) << "Result error message mentions Transcendental or Real";
}

TEST(AssumptionPropagationUnit, DiagnosticPositiveThenNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    auto failure_147 = ctx.assume_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_147.has_value())) << "Positive + Negative returns InvalidArgument";
    EXPECT_TRUE((failure_147.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
    const std::string &msg = failure_147.error().message;
    /// 诊断包含符号名与 Positive;实现可报告由 Negative 推导出的
    /// NonPositive 冲突.
    {
        const std::string actual_text = (msg);
        for (const auto &token : std::vector<std::string>{"x", "Positive"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Result error message contains 'x' and 'Positive'" << ": missing " << token << " in " << actual_text;
        }
    }
    bool has_negative_or_nonpositive =
        (msg.find("Negative") != std::string::npos) ||
        (msg.find("NonPositive") != std::string::npos);
    EXPECT_TRUE((has_negative_or_nonpositive)) << "Result error message mentions Negative or NonPositive";
}

TEST(AssumptionPropagationUnit, DiagnosticNaturalThenNegative) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Natural).has_value());

    auto failure_175 = ctx.assume_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_175.has_value())) << "Natural + Negative returns InvalidArgument";
    EXPECT_TRUE((failure_175.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
    const std::string &msg = failure_175.error().message;
    {
        const std::string actual_text = (msg);
        for (const auto &token : std::vector<std::string>{"Natural", "Negative"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Result error message contains 'Natural' and 'Negative'" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(AssumptionPropagationUnit, DiagnosticPositiveintThenZero) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());

    auto failure_195 = ctx.assume_sign("x", Sign::Zero);
    EXPECT_TRUE((!failure_195.has_value())) << "PositiveInt + Zero returns InvalidArgument";
    EXPECT_TRUE((failure_195.error().code == CasErrc::InvalidArgument)) << "failure reports InvalidArgument";
    const std::string &msg = failure_195.error().message;
    /// 诊断包含 PositiveInt,并可使用 Zero 推导出的 NonPositive
    /// 描述与隐含 Positive 的冲突.
    {
        const std::string actual_text = (msg);
        for (const auto &token : std::vector<std::string>{"PositiveInt"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Result error message contains 'PositiveInt'" << ": missing " << token << " in " << actual_text;
        }
    }
    bool has_zero_or_nonpositive =
        (msg.find("Zero") != std::string::npos) ||
        (msg.find("NonPositive") != std::string::npos);
    EXPECT_TRUE((has_zero_or_nonpositive)) << "Result error message mentions Zero or NonPositive";
}
