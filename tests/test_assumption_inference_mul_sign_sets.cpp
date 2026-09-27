#include "test_assumption_inference_mul_support.hpp"

static bool satisfies_product_sign(int value, Sign target) {
    switch (target) {
    case Sign::Positive:
        return value > 0;
    case Sign::Negative:
        return value < 0;
    case Sign::NonNegative:
        return value >= 0;
    case Sign::NonPositive:
        return value <= 0;
    case Sign::Zero:
        return value == 0;
    case Sign::NonZero:
        return value != 0;
    }
    return false;
}

static Tribool expected_product_sign(
    const std::vector<int> &left, const std::vector<int> &right, Sign target) {
    bool yes = false, no = false;
    for (int x : left) {
        for (int y : right) {
            if (satisfies_product_sign(x * y, target)) {
                yes = true;
            } else {
                no = true;
            }
        }
    }
    return yes && no ? Tribool::Unknown : yes ? Tribool::True
                                              : Tribool::False;
}

TEST(AssumptionInferenceMulSignSets, PossibleSignSets) {
    const std::vector<Sign> targets = {Sign::Positive, Sign::Negative,
                                       Sign::NonNegative, Sign::NonPositive, Sign::Zero, Sign::NonZero};
    const std::vector<std::vector<int>> values = {
        {1}, {-1}, {0}, {0, 1}, {-1, 0}, {-1, 1}, {-1, 0, 1}};
    const std::vector<Sign> declarations = {Sign::Positive, Sign::Negative,
                                            Sign::Zero, Sign::NonNegative, Sign::NonPositive, Sign::NonZero};
    for (std::size_t i = 0; i < values.size(); ++i) {
        for (std::size_t j = 0; j < values.size(); ++j) {
            AssumptionContext ctx;
            EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "x real";
            EXPECT_TRUE(ctx.assume_domain("y", Domain::Real).has_value()) << "y real";
            if (i < declarations.size()) {
                EXPECT_TRUE(ctx.assume_sign("x", declarations[i]).has_value()) << "x sign";
            }
            if (j < declarations.size()) {
                EXPECT_TRUE(ctx.assume_sign("y", declarations[j]).has_value()) << "y sign";
            }
            InferenceEngine engine(ctx);
            auto expression = test_expression_from_node(make_multiply({test_variable_node("x"), test_variable_node("y")}));
            for (Sign target : targets) {
                const auto expected = expected_product_sign(values[i], values[j], target);
                auto actual = test_sign_query(engine, expression, target);
                EXPECT_TRUE((actual && actual.value() == expected)) << "complete product truth table";
            }
        }
    }
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value()) << "weak sign";
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonZero).has_value()) << "nonzero intersection";
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value()) << "negative sign";
    InferenceEngine engine(ctx);
    auto expression = test_expression_from_node(make_multiply({test_variable_node("x"), test_variable_node("y"), test_integer_node(2)}));
    auto negative = engine.query_negative_checked(expression);
    EXPECT_TRUE((negative && negative.value() == Tribool::True)) << "intersection proves negative product";
}

TEST(AssumptionInferenceMulSignSets, ZeroRequiresDefinedRealFactors) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);
    auto unknown = test_expression_from_node(make_multiply({test_integer_node(0), test_variable_node("x")}));
    auto result = engine.query_nonnegative_checked(unknown);
    EXPECT_TRUE((result && result.value() == Tribool::Unknown)) << "unknown realness is not proved";
    auto inverse = detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(-1));
    auto possible_pole = test_expression_from_node(make_multiply({test_integer_node(0), inverse}));
    auto pole = engine.query_nonnegative_checked(possible_pole);
    EXPECT_TRUE((pole && pole.value() == Tribool::Unknown)) << "possibly zero denominator remains unknown";
    auto bad = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                               std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(-1)});
    for (auto operands : {
             std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(0), bad},
             std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x"), test_integer_node(0), bad}}) {
        auto expression = test_expression_from_node(make_multiply(std::move(operands)));
        auto invalid = engine.query_nonnegative_checked(expression);
        EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::DomainError)) << "later undefined factor is not hidden by zero or an earlier unknown";
    }
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonZero).has_value()) << "nonzero accepted";
    auto nonreal = engine.query_nonnegative_checked(unknown);
    EXPECT_TRUE((nonreal && nonreal.value() == Tribool::Unknown)) << "NonZero alone does not imply Real";
}
