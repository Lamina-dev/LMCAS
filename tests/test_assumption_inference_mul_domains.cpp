
#include "test_assumption_inference_mul_support.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_number_real(double val) {
    return LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val));
}

TEST(AssumptionInferenceMulDomains, AllIntegerNumbers) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {test_integer_node(3), test_integer_node(-5), test_integer_node(7)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "3 * (-5) * 7: Integer";
}

TEST(AssumptionInferenceMulDomains, AllIntegerVariables) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("z", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {test_variable_node("x"), test_variable_node("y"), test_variable_node("z")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "int_x * int_y * int_z: Integer";
}

TEST(AssumptionInferenceMulDomains, MixedIntegerAndNumber) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_integer_node(5)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "int_x * 5: Integer";
}

TEST(AssumptionInferenceMulDomains, AllRealNumbers) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply(
        {make_number_real(2.5), make_number_real(-1.5),
         make_number_real(3.0)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "2.5 * (-1.5) * 3.0: Real";
}

TEST(AssumptionInferenceMulDomains, AllRealVariables) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "real_x * real_y: Real";
}

TEST(AssumptionInferenceMulDomains, IntegerImpliesReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "int_x * int_y: also Real (Integer subset of Real)";
}

TEST(AssumptionInferenceMulDomains, MixedIntegerRealIsReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "int_x * real_y: Real";
    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "int_x * real_y: Integer is Unknown";
}

TEST(AssumptionInferenceMulDomains, UnknownDomainPropagation) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "int_x * unknown_y: Integer Unknown";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::Unknown)) << "int_x * unknown_y: Real Unknown";
}

TEST(AssumptionInferenceMulDomains, NumberAndRealVariable) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_integer_node(5), test_variable_node("x")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "5 * real_x: Real";
}

TEST(AssumptionInferenceMulDomains, NaturalDomainImpliesInteger) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Natural).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Natural).has_value());
    InferenceEngine engine(ctx);

    auto mul_node = make_multiply({test_variable_node("x"), test_variable_node("y")});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::True)) << "nat_x * nat_y: Integer (Natural implies Integer)";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "nat_x * nat_y: Real (Natural implies Real)";
}

TEST(AssumptionInferenceMulDomains, EmptyMultiplyDomain) {
    bool rejected = false;
    try {
        (void)LMCAS::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "Invalid empty products are rejected before inference";
}

TEST(AssumptionInferenceMulDomains, RationalNotInteger) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    auto rat_node = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
    auto mul_node = make_multiply({rat_node, test_integer_node(3)});
    auto expr = test_expression_from_node(mul_node);

    EXPECT_TRUE((engine.query_integer_checked(expr).value() == Tribool::Unknown)) << "Rational(1/2) * 3: Integer Unknown (non-integer operand)";
    EXPECT_TRUE((engine.query_real_checked(expr).value() == Tribool::True)) << "Rational(1/2) * 3: Real";
}

TEST(AssumptionInferenceMulDomains, ProductDomainObligations) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);
    auto logarithm = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                                     std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(-1)});
    auto expression = test_expression_from_node(make_multiply({test_variable_node("unknown"), test_integer_node(0), logarithm}));
    for (Domain domain : {Domain::Real, Domain::Integer}) {
        auto result = domain == Domain::Real ? engine.query_real_checked(expression) : engine.query_integer_checked(expression);
        EXPECT_TRUE((!result && result.error().code == CasErrc::DomainError)) << "unknown or zero earlier factors cannot erase later undefinedness";
    }
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonZero).has_value()) << "nonzero assumption";
    auto variable = test_expression_from_node(test_variable_node("x"));
    auto real = engine.query_real_checked(variable);
    EXPECT_TRUE((real && real.value() == Tribool::Unknown)) << "NonZero is not Real";
    for (Sign sign : {Sign::NonNegative, Sign::NonPositive}) {
        AssumptionContext ordered;
        EXPECT_TRUE((ordered.assume_sign("x", sign).has_value())) << "ordered sign";
        InferenceEngine ordered_engine(ordered);
        auto product = test_expression_from_node(make_multiply({test_variable_node("x"), test_integer_node(2)}));
        auto result = ordered_engine.query_real_checked(product);
        EXPECT_TRUE((result && result.value() == Tribool::True)) << "weak ordered signs imply realness";
    }
}
