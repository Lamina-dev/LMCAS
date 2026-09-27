
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;

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

/// Generate a random domain that is Integer or Real (for trig function arguments)
static Domain random_integer_or_real() {
    return *rc::gen::arbitrary<bool>() ? Domain::Integer : Domain::Real;
}

/// Generate a random domain that is Rational or Real (for exp arguments)
static Domain random_rational_or_real() {
    int choice = *rc::gen::inRange(0, (2) + 1);
    switch (choice) {
    case 0:
        return Domain::Rational;
    case 1:
        return Domain::Real;
    default:
        return Domain::Integer; // Integer implies Rational
    }
}

TEST(LmcasAssumptionDomainInference, TrigIntegerOrRealGivesReal) {
    EXPECT_TRUE(rc::check("For any trig function with Integer or Real argument, result is Real", []() {
        // Pick a random trig function
        std::vector<FunctionNode::FuncType> trig_funcs = {
            FunctionNode::FuncType::Sin,
            FunctionNode::FuncType::Cos,
            FunctionNode::FuncType::Tan};
        auto func_type = *rc::gen::elementOf(trig_funcs);

        // Pick a random domain for the argument
        Domain arg_domain = random_integer_or_real();

        std::string var_name = "x_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(var_name, arg_domain).has_value());
        InferenceEngine engine(ctx);

        auto func_node = make_function(func_type, test_variable_node(var_name));
        auto expr = test_expression_from_node(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDomainInference, ExpRationalOrRealGivesReal) {
    EXPECT_TRUE(rc::check("For exp with Rational or Real argument, result is Real", []() {
        Domain arg_domain = random_rational_or_real();
        std::string var_name = "x_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(var_name, arg_domain).has_value());
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Exp, test_variable_node(var_name));
        auto expr = test_expression_from_node(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDomainInference, LnIntegerRequiresPositive) {
    for (Domain domain : {Domain::Integer, Domain::PositiveInt}) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain("n", domain).has_value()) << "domain accepted";
        InferenceEngine engine(ctx);
        auto expression = test_expression_from_node(make_function(FunctionNode::FuncType::Ln, test_variable_node("n")));
        auto result = engine.query_real_checked(expression);
        EXPECT_TRUE(result && result.value() == (domain == Domain::PositiveInt ? Tribool::True : Tribool::Unknown)) << "ln requires strict positivity";
    }
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value()) << "integer accepted";
    EXPECT_TRUE(ctx.assume_sign("n", Sign::Negative).has_value()) << "negative accepted";
    InferenceEngine engine(ctx);
    auto logarithm = make_function(FunctionNode::FuncType::Ln, test_variable_node("n"));
    for (auto node : {logarithm,
                      make_function(FunctionNode::FuncType::Exp, logarithm),
                      make_power(logarithm, test_integer_node(2))}) {
        auto expression = test_expression_from_node(node);
        auto real = engine.query_real_checked(expression);
        auto nonnegative = engine.query_nonnegative_checked(expression);
        EXPECT_TRUE(!real && real.error().code == CasErrc::DomainError) << "undefined nested logarithm cannot establish realness";
        EXPECT_TRUE(!nonnegative && nonnegative.error().code == CasErrc::DomainError) << "undefined nested logarithm cannot establish a sign";
    }
}

TEST(LmcasAssumptionDomainInference, SqrtNonnegRealGivesReal) {
    EXPECT_TRUE(rc::check("For sqrt with non-negative Real argument, result is Real", []() {
        std::string var_name = "x_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(var_name, Domain::Real).has_value());
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::NonNegative).has_value());
        InferenceEngine engine(ctx);

        auto func_node = make_function(FunctionNode::FuncType::Sqrt, test_variable_node(var_name));
        auto expr = test_expression_from_node(func_node);

        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDomainInference, IntegerPowerNaturalGivesInteger) {
    EXPECT_TRUE(rc::check("For Integer base raised to a positive integer exponent, result is Integer", []() {
        std::string base_name = "b_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
        int exponent = *rc::gen::inRange(1, (10) + 1); // Positive integer exponent

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(base_name, Domain::Integer).has_value());
        InferenceEngine engine(ctx);

        auto pow_node = make_power(test_variable_node(base_name), test_integer_node(exponent));
        auto expr = test_expression_from_node(pow_node);

        RC_ASSERT(engine.query_integer_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDomainInference, IntegerPowerZeroRequiresNonzero) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value()) << "integer accepted";
    InferenceEngine engine(ctx);
    const auto expression = test_expression_from_node(make_power(test_variable_node("n"), test_integer_node(0)));
    auto unknown = engine.query_integer_checked(expression);
    EXPECT_TRUE(unknown && unknown.value() == Tribool::Unknown) << "possibly zero remains unknown";
    EXPECT_TRUE(ctx.assume_sign("n", Sign::NonZero).has_value()) << "nonzero accepted";
    auto known = engine.query_integer_checked(expression);
    EXPECT_TRUE(known && known.value() == Tribool::True) << "nonzero integer zero power is integer";
}

TEST(LmcasAssumptionDomainInference, RationalPowerIntegerGivesReal) {
    EXPECT_TRUE(rc::check("For Rational base raised to integer exponent, result is Real", []() {
        std::string base_name = "r_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
        int exponent = *rc::gen::inRange(1, (10) + 1);

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(base_name, Domain::Rational).has_value());
        InferenceEngine engine(ctx);

        auto pow_node = make_power(test_variable_node(base_name), test_integer_node(exponent));
        auto expr = test_expression_from_node(pow_node);

        // Rational implies Real, and Real^Integer -> Real
        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionDomainInference, UnknownDomainGivesUnknown) {
    EXPECT_TRUE(rc::check("For functions with unknown domain argument, result domain is Unknown", []() {
        std::vector<FunctionNode::FuncType> funcs = {
            FunctionNode::FuncType::Sin,
            FunctionNode::FuncType::Cos,
            FunctionNode::FuncType::Tan,
            FunctionNode::FuncType::Exp};
        auto func_type = *rc::gen::elementOf(funcs);
        std::string var_name = "u_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        // No domain declared for variable (defaults to Complex)
        AssumptionContext ctx;
        InferenceEngine engine(ctx);

        auto func_node = make_function(func_type, test_variable_node(var_name));
        auto expr = test_expression_from_node(func_node);

        /// 参数域缺少 Real/Integer 证明时,函数值实数性保持 Unknown.
        RC_ASSERT(engine.query_real_checked(expr).value() == Tribool::Unknown);
    }));
}

TEST(LmcasAssumptionDomainInference, AllTrigWithInteger) {
    // sin(Integer) -> Real
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());
        InferenceEngine engine(ctx);
        auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Sin, test_variable_node("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True) << "sin(Integer) → Real";
    }
    // cos(Integer) -> Real
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());
        InferenceEngine engine(ctx);
        auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Cos, test_variable_node("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True) << "cos(Integer) → Real";
    }
    // tan(Integer) -> Real
    {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());
        InferenceEngine engine(ctx);
        auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Tan, test_variable_node("n")));
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True) << "tan(Integer) → Real";
    }
}

TEST(LmcasAssumptionDomainInference, ExpWithInteger) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());
    InferenceEngine engine(ctx);
    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Exp, test_variable_node("n")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True) << "exp(Integer) → Real";
}

TEST(LmcasAssumptionDomainInference, LnPositiveGivesReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    InferenceEngine engine(ctx);
    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Ln, test_variable_node("x")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True) << "ln(Positive Real) → Real";
}

TEST(LmcasAssumptionDomainInference, SqrtWithoutNonnegUnknown) {
    /// sqrt 的实数性需要 Real 与 NonNegative 共同证明;当前仅声明 Real.
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    /// 符号属性保持未声明状态.
    InferenceEngine engine(ctx);
    auto expr = test_expression_from_node(make_function(FunctionNode::FuncType::Sqrt, test_variable_node("x")));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown) << "sqrt(Real without NonNeg) → Unknown";
}

TEST(LmcasAssumptionDomainInference, PowerNegativeExponentNotInteger) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());
    InferenceEngine engine(ctx);

    // n^(-1) = 1/n - not necessarily integer
    auto pow_node = make_power(test_variable_node("n"), test_integer_node(-1));
    auto expr = test_expression_from_node(pow_node);

    // Should NOT be able to infer Integer (e.g., 2^(-1) = 0.5)
    EXPECT_TRUE(engine.query_integer_checked(expr).value() == Tribool::Unknown) << "Integer^(-1) → Integer is Unknown (not guaranteed)";
}
