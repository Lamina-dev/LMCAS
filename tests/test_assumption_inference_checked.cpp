
#include "test_assumption_inference_mul_support.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "computation_context.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"
#include "expr.hpp"
#include <vector>
#include <string>
#include <memory>
#include <type_traits>

using namespace LMCAS;

using Node = std::shared_ptr<const SymbolicNode>;

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

static std::shared_ptr<const SymbolicNode> make_division(
    std::shared_ptr<const SymbolicNode> num, std::shared_ptr<const SymbolicNode> den) {
    auto den_inv = make_power(std::move(den), test_integer_node(-1));
    return make_multiply({std::move(num), std::move(den_inv)});
}

static std::shared_ptr<const SymbolicNode> make_subtraction(
    std::shared_ptr<const SymbolicNode> a, std::shared_ptr<const SymbolicNode> b) {
    auto neg_b = make_multiply({test_integer_node(-1), std::move(b)});
    return make_add({std::move(a), std::move(neg_b)});
}

static void test_checked_inference_query_contracts_symbol_signs(InferenceEngine &engine) {
    auto x = test_expression_from_node(test_variable_node("x"));
    auto positive = engine.query_positive_checked(x);
    ASSERT_TRUE((positive.has_value())) << "checked query_positive succeeds";
    if (positive) {
        EXPECT_TRUE((positive.value() == Tribool::True)) << "checked query_positive returns True for positive symbol";
    }

    auto nonpositive = engine.query_nonpositive_checked(x);
    ASSERT_TRUE((nonpositive.has_value())) << "checked query_nonpositive succeeds";
    if (nonpositive) {
        EXPECT_TRUE((nonpositive.value() == Tribool::False)) << "checked query_nonpositive returns False for positive symbol";
    }

    auto integer = engine.query_integer_checked(x);
    ASSERT_TRUE((integer.has_value())) << "checked query_integer succeeds";
    if (integer) {
        EXPECT_TRUE((integer.value() == Tribool::True)) << "checked query_integer returns True for integer symbol";
    }

    auto nonzero = engine.query_nonzero_checked(x);
    ASSERT_TRUE((nonzero.has_value())) << "checked query_nonzero succeeds";
    if (nonzero) {
        EXPECT_TRUE((nonzero.value() == Tribool::True)) << "checked query_nonzero returns True for positive symbol";
    }
}

static void test_checked_inference_query_contracts_numbers_and_exponential(InferenceEngine &engine) {
    auto three = test_expression_from_node(test_integer_node(3));
    auto numeric_negative = engine.query_negative_checked(three);
    ASSERT_TRUE((numeric_negative.has_value())) << "checked query_negative succeeds";
    if (numeric_negative) {
        EXPECT_TRUE((numeric_negative.value() == Tribool::False)) << "checked query_negative returns False for positive number";
    }

    auto numeric_nonnegative = engine.query_nonnegative_checked(three);
    ASSERT_TRUE((numeric_nonnegative.has_value())) << "checked query_nonnegative succeeds";
    if (numeric_nonnegative) {
        EXPECT_TRUE((numeric_nonnegative.value() == Tribool::True)) << "checked query_nonnegative returns True for positive number";
    }

    auto numeric_real = engine.query_real_checked(three);
    ASSERT_TRUE((numeric_real.has_value())) << "checked query_real succeeds";
    if (numeric_real) {
        EXPECT_TRUE((numeric_real.value() == Tribool::True)) << "checked query_real returns True for exact number";
    }

    auto exp_rational = test_expression_from_node(make_function(FunctionNode::FuncType::Exp, test_variable_node("q")));
    auto exp_rational_real = engine.query_real_checked(exp_rational);
    ASSERT_TRUE((exp_rational_real.has_value())) << "checked query_real succeeds for exp(rational)";
    if (exp_rational_real) {
        EXPECT_TRUE((exp_rational_real.value() == Tribool::True)) << "checked query_real uses checked internal rational-domain dispatch";
    }
}

static void test_checked_inference_query_contracts_arithmetic_domains(InferenceEngine &engine) {
    auto integer_sum = test_expression_from_node(make_add({test_variable_node("x"), test_variable_node("x")}));
    auto integer_sum_checked = engine.query_integer_checked(integer_sum);
    ASSERT_TRUE((integer_sum_checked.has_value())) << "checked query_integer succeeds for integer addition";
    if (integer_sum_checked) {
        EXPECT_TRUE((integer_sum_checked.value() == Tribool::True)) << "checked query_integer uses checked addition-domain inference";
    }

    auto real_product_domain = test_expression_from_node(make_multiply({test_variable_node("x"), test_variable_node("real_symbol")}));
    auto real_product_checked = engine.query_real_checked(real_product_domain);
    ASSERT_TRUE((real_product_checked.has_value())) << "checked query_real succeeds for mixed integer-real product";
    if (real_product_checked) {
        EXPECT_TRUE((real_product_checked.value() == Tribool::True)) << "checked query_real uses checked multiplication-domain inference";
    }

    auto integer_power_domain = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(2)));
    auto integer_power_domain_checked = engine.query_integer_checked(integer_power_domain);
    ASSERT_TRUE((integer_power_domain_checked.has_value())) << "checked query_integer succeeds for integer nonnegative power";
    if (integer_power_domain_checked) {
        EXPECT_TRUE((integer_power_domain_checked.value() == Tribool::True)) << "checked query_integer uses checked power-domain inference";
    }

    auto abs_integer_domain = test_expression_from_node(make_function(FunctionNode::FuncType::Abs, test_variable_node("x")));
    auto abs_integer_domain_checked = engine.query_integer_checked(abs_integer_domain);
    ASSERT_TRUE((abs_integer_domain_checked.has_value())) << "checked query_integer succeeds for abs(integer)";
    if (abs_integer_domain_checked) {
        EXPECT_TRUE((abs_integer_domain_checked.value() == Tribool::True)) << "checked query_integer uses checked function-domain inference";
    }
}

static void test_checked_inference_query_contracts_addition_signs(AssumptionContext &ctx, InferenceEngine &engine) {
    auto positive_sum = test_expression_from_node(make_add({test_variable_node("x"), test_variable_node("x")}));
    auto positive_sum_checked = engine.query_positive_checked(positive_sum);
    ASSERT_TRUE((positive_sum_checked.has_value())) << "checked query_positive succeeds for addition";
    if (positive_sum_checked) {
        EXPECT_TRUE((positive_sum_checked.value() == Tribool::True)) << "checked query_positive uses checked addition-sign inference";
    }

    EXPECT_TRUE(ctx.assume_sign("neg_symbol", Sign::Negative).has_value());
    auto positive_difference = test_expression_from_node(make_subtraction(test_variable_node("x"), test_variable_node("neg_symbol")));
    auto positive_difference_checked = engine.query_positive_checked(positive_difference);
    ASSERT_TRUE((positive_difference_checked.has_value())) << "checked query_positive succeeds for subtraction-shaped addition";
    if (positive_difference_checked) {
        EXPECT_TRUE((positive_difference_checked.value() == Tribool::True)) << "checked query_positive uses checked subtraction-sign inference";
    }

    auto relation_positive_sum = test_expression_from_node(make_add({test_variable_node("rel_a"), test_variable_node("rel_b")}));
    auto zero_expr = test_expression_from_node(test_integer_node(0));
    auto relation_inserted = ctx.current_relations().add_relation_checked(
        relation_positive_sum, zero_expr, RelationalNode::Op::GT, ctx.current_properties());
    EXPECT_TRUE((relation_inserted.has_value())) << "checked relation insertion succeeds for composite positive sum";
    auto relation_positive_sum_checked = engine.query_positive_checked(relation_positive_sum);
    ASSERT_TRUE((relation_positive_sum_checked.has_value())) << "checked query_positive succeeds for relation-backed composite sum";
    if (relation_positive_sum_checked) {
        EXPECT_TRUE((relation_positive_sum_checked.value() == Tribool::True)) << "checked query_positive uses checked relation sign inference";
    }
}

static void test_checked_inference_query_contracts_product_signs(InferenceEngine &engine) {
    auto positive_product = test_expression_from_node(make_multiply({test_variable_node("x"), test_variable_node("x")}));
    auto positive_product_checked = engine.query_positive_checked(positive_product);
    ASSERT_TRUE((positive_product_checked.has_value())) << "checked query_positive succeeds for multiplication";
    if (positive_product_checked) {
        EXPECT_TRUE((positive_product_checked.value() == Tribool::True)) << "checked query_positive uses checked multiplication-sign inference";
    }

    auto negative_product = test_expression_from_node(make_multiply({test_variable_node("x"), test_variable_node("neg_symbol")}));
    auto negative_product_checked = engine.query_negative_checked(negative_product);
    ASSERT_TRUE((negative_product_checked.has_value())) << "checked query_negative succeeds for multiplication";
    if (negative_product_checked) {
        EXPECT_TRUE((negative_product_checked.value() == Tribool::True)) << "checked query_negative uses checked multiplication-sign inference";
    }

    auto division_expr = test_expression_from_node(make_division(test_variable_node("x"), test_variable_node("x")));
    auto division_positive = engine.query_positive_checked(division_expr);
    ASSERT_TRUE((division_positive.has_value())) << "checked query_positive succeeds for division pattern";
    if (division_positive) {
        EXPECT_TRUE((division_positive.value() == Tribool::True)) << "checked query_positive uses checked division-sign inference";
    }
}

static void test_checked_inference_query_contracts_power_signs(InferenceEngine &engine) {
    auto positive_power = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(2)));
    auto positive_power_checked = engine.query_positive_checked(positive_power);
    ASSERT_TRUE((positive_power_checked.has_value())) << "checked query_positive succeeds for power";
    if (positive_power_checked) {
        EXPECT_TRUE((positive_power_checked.value() == Tribool::True)) << "checked query_positive uses checked positive-base power inference";
    }

    auto even_power = test_expression_from_node(make_power(test_variable_node("unknown_symbol"), test_integer_node(2)));
    auto even_power_nonnegative = engine.query_nonnegative_checked(even_power);
    ASSERT_TRUE((even_power_nonnegative.has_value())) << "checked query_nonnegative succeeds for even power";
    if (even_power_nonnegative) {
        EXPECT_TRUE((even_power_nonnegative.value() == Tribool::Unknown)) << "checked even-power inference preserves Unknown without real-domain proof";
    }
    auto real_even_power = test_expression_from_node(make_power(test_variable_node("real_symbol"), test_integer_node(2)));
    auto real_even_power_nonnegative = engine.query_nonnegative_checked(real_even_power);
    ASSERT_TRUE((real_even_power_nonnegative.has_value())) << "checked query_nonnegative succeeds for real even power";
    if (real_even_power_nonnegative) {
        EXPECT_TRUE((real_even_power_nonnegative.value() == Tribool::True)) << "checked query_nonnegative uses checked even-power inference";
    }

    auto nonzero_power = test_expression_from_node(make_power(test_variable_node("x"), test_integer_node(-1)));
    auto nonzero_power_checked = engine.query_nonzero_checked(nonzero_power);
    ASSERT_TRUE((nonzero_power_checked.has_value())) << "checked query_nonzero succeeds for integer power";
    if (nonzero_power_checked) {
        EXPECT_TRUE((nonzero_power_checked.value() == Tribool::True)) << "checked query_nonzero uses checked nonzero-base power inference";
    }
}

static void test_checked_inference_query_contracts_function_signs(InferenceEngine &engine) {
    auto exp_real = test_expression_from_node(make_function(FunctionNode::FuncType::Exp, test_variable_node("real_symbol")));
    auto exp_real_positive = engine.query_positive_checked(exp_real);
    ASSERT_TRUE((exp_real_positive.has_value())) << "checked query_positive succeeds for exp(real)";
    if (exp_real_positive) {
        EXPECT_TRUE((exp_real_positive.value() == Tribool::True)) << "checked query_positive uses checked exp function-sign inference";
    }

    auto abs_real_pos = test_expression_from_node(make_function(FunctionNode::FuncType::Abs, test_variable_node("real_pos")));
    auto abs_real_pos_positive = engine.query_positive_checked(abs_real_pos);
    ASSERT_TRUE((abs_real_pos_positive.has_value())) << "checked query_positive succeeds for abs(positive real)";
    if (abs_real_pos_positive) {
        EXPECT_TRUE((abs_real_pos_positive.value() == Tribool::True)) << "checked query_positive uses checked abs function-sign inference";
    }

    auto sqrt_nonnegative = test_expression_from_node(make_function(FunctionNode::FuncType::Sqrt, test_variable_node("nn_symbol")));
    auto sqrt_nonnegative_checked = engine.query_nonnegative_checked(sqrt_nonnegative);
    ASSERT_TRUE((sqrt_nonnegative_checked.has_value())) << "checked query_nonnegative succeeds for sqrt(nonnegative)";
    if (sqrt_nonnegative_checked) {
        EXPECT_TRUE((sqrt_nonnegative_checked.value() == Tribool::True)) << "checked query_nonnegative uses checked sqrt function-sign inference";
    }
}

static void test_checked_inference_query_contracts_extended_properties(InferenceEngine &engine) {
    auto finite = engine.query_finite_checked(test_expression_from_node(test_variable_node("finite_symbol")));
    ASSERT_TRUE((finite.has_value())) << "checked query_finite succeeds";
    if (finite) {
        EXPECT_TRUE((finite.value() == Tribool::True)) << "checked query_finite returns True for finite symbol";
    }

    auto periodic_expr = test_expression_from_node(test_variable_node("periodic_symbol"));
    auto periodic = engine.query_periodic_checked(periodic_expr, "x");
    ASSERT_TRUE((periodic.has_value())) << "checked query_periodic succeeds";
    if (periodic) {
        EXPECT_TRUE((periodic.value() == Tribool::True)) << "checked query_periodic returns True for periodic symbol";
    }

    auto period = engine.infer_period_checked(periodic_expr, "x");
    ASSERT_TRUE((period.has_value())) << "checked infer_period succeeds";
    if (period) {
        EXPECT_TRUE((period.value().has_value())) << "checked infer_period returns declared period";
    }

    auto unknown = engine.query_algebraic_checked(test_expression_from_node(test_variable_node("unknown_symbol")));
    ASSERT_TRUE((unknown.has_value())) << "checked query_algebraic accepts valid unknown symbol";
    if (unknown) {
        EXPECT_TRUE((unknown.value() == Tribool::Unknown)) << "checked query_algebraic preserves Unknown for valid unsupported facts";
    }

    auto tau_expr = test_expression_from_node(test_variable_node("tau_symbol"));
    auto transcendental = engine.query_transcendental_checked(tau_expr);
    ASSERT_TRUE((transcendental.has_value())) << "checked query_transcendental succeeds";
    if (transcendental) {
        EXPECT_TRUE((transcendental.value() == Tribool::True)) << "checked query_transcendental returns True for transcendental symbol";
    }

    auto divergent_expr = test_expression_from_node(test_variable_node("divergent_symbol"));
    auto divergent = engine.query_divergent_checked(divergent_expr);
    ASSERT_TRUE((divergent.has_value())) << "checked query_divergent succeeds";
    if (divergent) {
        EXPECT_TRUE((divergent.value() == Tribool::True)) << "checked query_divergent returns True for divergent symbol";
    }
    auto divergent_finite = engine.query_finite_checked(divergent_expr);
    ASSERT_TRUE((divergent_finite.has_value())) << "checked query_finite succeeds for divergent symbol";
    if (divergent_finite) {
        EXPECT_TRUE((divergent_finite.value() == Tribool::False)) << "checked query_finite returns False for divergent symbol";
    }
}

TEST(AssumptionInferenceChecked, CheckedInferenceQueryContracts) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    ASSERT_TRUE(ctx.assume_domain("q", Domain::Rational).has_value());
    ASSERT_TRUE(ctx.assume_domain("real_symbol", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("real_pos", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("real_pos", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("nn_symbol", Sign::NonNegative).has_value());
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("finite_symbol", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("divergent_symbol", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_transcendental("tau_symbol").has_value())) << "transcendental declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_periodic("periodic_symbol", "x", test_expression_from_node(test_integer_node(6))).has_value())) << "period declaration succeeds";
    InferenceEngine engine(ctx);

    test_checked_inference_query_contracts_symbol_signs(engine);
    test_checked_inference_query_contracts_numbers_and_exponential(engine);
    test_checked_inference_query_contracts_arithmetic_domains(engine);
    test_checked_inference_query_contracts_addition_signs(ctx, engine);
    test_checked_inference_query_contracts_product_signs(engine);
    test_checked_inference_query_contracts_power_signs(engine);
    test_checked_inference_query_contracts_function_signs(engine);
    test_checked_inference_query_contracts_extended_properties(engine);
}

TEST(AssumptionInferenceChecked, CancelledContextRejectsLiteralAndDeclaredFacts) {
    AssumptionContext assumptions;
    ASSERT_TRUE(assumptions.assume_sign("x", Sign::Positive));
    InferenceEngine engine(assumptions);
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext context({}, cancellation);
    const auto literal = test_expression_from_node(test_integer_node(1));
    const auto declared = test_expression_from_node(test_variable_node("x"));
    const auto infinity = detail::expression_from_node(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{}));
    using Query = InferenceTriboolResult (InferenceEngine::*)(
        const SymbolicExpr&, ComputationContext&) const;
    const Query queries[] = {
        &InferenceEngine::query_positive_checked,
        &InferenceEngine::query_negative_checked,
        &InferenceEngine::query_nonnegative_checked,
        &InferenceEngine::query_nonpositive_checked,
        &InferenceEngine::query_real_checked,
        &InferenceEngine::query_integer_checked,
        &InferenceEngine::query_nonzero_checked,
        &InferenceEngine::query_algebraic_checked,
        &InferenceEngine::query_transcendental_checked,
        &InferenceEngine::query_finite_checked,
        &InferenceEngine::query_divergent_checked,
    };
    for (const auto query : queries) {
        for (const auto* expression : {&literal, &declared, &infinity}) {
            auto result = (engine.*query)(*expression, context);
            ASSERT_FALSE(result);
            EXPECT_EQ(result.error().code, CasErrc::Cancelled);
        }
    }
    auto adapted = assumptions.is_positive_checked(declared, context);
    ASSERT_FALSE(adapted);
    EXPECT_EQ(adapted.error().code, CasErrc::Cancelled);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(AssumptionInferenceChecked, QueriesShareCumulativeStepBudget) {
    AssumptionContext assumptions;
    InferenceEngine engine(assumptions);
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);
    const auto one = test_expression_from_node(test_integer_node(1));
    auto first = engine.query_positive_checked(one, context);
    ASSERT_TRUE(first);
    EXPECT_EQ(first.value(), Tribool::True);
    auto second = engine.query_positive_checked(one, context);
    ASSERT_FALSE(second);
    EXPECT_EQ(second.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(AssumptionInferenceChecked, ResourceDepthFailureRestoresBothGuards) {
    AssumptionContext assumptions;
    InferenceEngine engine(assumptions);
    engine.set_max_depth(1);
    auto nested = test_integer_node(1);
    for (int i = 0; i < 8; ++i) {
        nested = make_function(FunctionNode::FuncType::Exp, nested);
    }
    ResourceLimits limits;
    limits.max_recursion_depth = 1;
    ComputationContext context(limits);
    auto deep = engine.query_positive_checked(test_expression_from_node(nested), context);
    ASSERT_FALSE(deep);
    EXPECT_EQ(deep.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(context.recursion_depth(), 0u);
    auto shallow = engine.query_positive_checked(test_expression_from_node(test_integer_node(1)), context);
    ASSERT_TRUE(shallow);
    EXPECT_EQ(shallow.value(), Tribool::True);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(AssumptionInferenceChecked, PeriodQueriesRespectCallerLimitsAndCancellation) {
    AssumptionContext assumptions;
    InferenceEngine engine(assumptions);
    const auto sine = test_expression_from_node(make_function(FunctionNode::FuncType::Sin, test_variable_node("x")));
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    auto periodic = engine.query_periodic_checked(sine, "x", exhausted);
    ASSERT_FALSE(periodic);
    EXPECT_EQ(periodic.error().code, CasErrc::ResourceLimit);
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto period = engine.infer_period_checked(sine, "x", cancelled);
    ASSERT_FALSE(period);
    EXPECT_EQ(period.error().code, CasErrc::Cancelled);
    EXPECT_EQ(cancelled.recursion_depth(), 0u);
}

static void expect_integral_exponent_projection(
    const ExprPtr &reciprocal, AssumptionContext &assumptions,
    const FactsQuery &facts, ComputationContext &context) {
    auto domain = detail::domain_constraints(detail::node(reciprocal), facts, Domain::Real, context);
    EXPECT_TRUE((domain && domain.value().has_value())) << "integer arithmetic exponent has an exact projected domain";
    if (domain && domain.value()) {
        for (int value : {-1, 0, 1}) {
            Tribool accepted = Tribool::True;
            for (const auto &condition : *domain.value()) {
                const auto truth = assumptions.evaluate_condition(
                    *condition->substitute("x", SymbolicExpr::number(value)));
                if (truth != Tribool::True) {
                    accepted = truth;
                    break;
                }
            }
            EXPECT_TRUE((accepted == (value ? Tribool::True : Tribool::False))) << "negative bases remain allowed and zero remains excluded";
        }
    }
}

TEST(LmcasDomainPreservingSimplify, RawArithmeticExponentDomains) {
    AssumptionContext assumptions;
    EXPECT_TRUE((assumptions.assume_domain("x", Domain::Real).has_value())) << "real base declared";
    EXPECT_TRUE((assumptions.assume_domain("n", Domain::Real).has_value())) << "real exponent declared";
    detail::AssumptionFacts facts(assumptions);
    ComputationContext context;
    const std::pair<const char *, Tribool> cases[] = {
        {"(-2)^(-1*2)", Tribool::True},
        {"(-2)^(6/3)", Tribool::True},
        {"(-2)^(3/2)", Tribool::False},
        {"0^(-1*2)", Tribool::False},
        {"(-2)^(0/x)", Tribool::Unknown},
        {"(-2)^n", Tribool::Unknown}};
    for (const auto &[source, expected] : cases) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << source;
        if (!expression) {
            continue;
        }
        auto domain = detail::query_definedness(detail::node(expression.value()), facts, Domain::Real, context);
        EXPECT_TRUE((domain && domain.value() == expected)) << source;
    }
    auto reciprocal = parse_expr("x^(-1*2)").value();
    expect_integral_exponent_projection(reciprocal, assumptions, facts, context);
    auto unknown_exponent = parse_expr("(-2)^n").value();
    auto unresolved = detail::domain_constraints(
        detail::node(unknown_exponent), facts, Domain::Real, context);
    EXPECT_TRUE((unresolved && !unresolved.value())) << "unknown integrality cannot be replaced by a positive-base restriction";
    ResourceLimits limits;
    limits.max_integer_bits = 32;
    ComputationContext bounded(limits);
    auto huge_exponent = parse_expr("(-2)^(2^200)").value();
    auto exhausted = detail::domain_constraints(
        detail::node(huge_exponent), facts, Domain::Real, bounded);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "closed exponent arithmetic respects the caller's integer budget";
}

TEST(LmcasDomainPreservingSimplify, SymbolicZeroExponentDomains) {
    AssumptionContext assumptions;
    ASSERT_TRUE(assumptions.assume_sign_checked("z", Sign::Zero));
    detail::AssumptionFacts facts(assumptions);
    ComputationContext context;
    const Node zero = detail::make_node<NumberNode>(BigInt(0));
    const Node minus_one = detail::make_node<NumberNode>(BigInt(-1));
    const Node z = detail::make_node<VariableNode>("z");
    const Node sum = detail::make_node<AddNode>(std::vector<Node>{z, z});
    const Node product = detail::make_node<MultiplyNode>(
        std::vector<Node>{detail::make_node<NumberNode>(BigInt(2)), z});
    for (const auto& exponent : {z, sum, product}) {
        const auto power = make_power(zero, exponent);
        auto complex = detail::query_definedness(power, facts, Domain::Complex, context);
        ASSERT_TRUE(complex);
        EXPECT_EQ(complex.value(), Tribool::True);
        auto real = detail::query_definedness(power, facts, Domain::Real, context);
        ASSERT_TRUE(real);
        EXPECT_EQ(real.value(), Tribool::False);
        auto nonzero = detail::query_nonzero_value(power, facts, Domain::Complex, context);
        ASSERT_TRUE(nonzero);
        EXPECT_EQ(nonzero.value(), Tribool::True);
        for (const auto& dependent : {
                 make_function(FunctionNode::FuncType::Ln, power),
                 make_power(power, minus_one)}) {
            auto defined = detail::query_definedness(dependent, facts, Domain::Complex, context);
            ASSERT_TRUE(defined);
            EXPECT_EQ(defined.value(), Tribool::True);
        }
        auto invalid = detail::query_definedness(
            make_power(make_power(zero, minus_one), exponent),
            facts, Domain::Complex, context);
        ASSERT_TRUE(invalid);
        EXPECT_EQ(invalid.value(), Tribool::False);
    }
}

TEST(LmcasDomainPreservingSimplify, ZeroBaseExponentSignDomains) {
    struct Case {
        std::optional<Sign> sign;
        bool real;
        Tribool complex_domain;
        Tribool real_domain;
    };
    const Case cases[] = {
        {Sign::Zero, true, Tribool::True, Tribool::False},
        {Sign::Positive, true, Tribool::True, Tribool::True},
        {Sign::Negative, true, Tribool::False, Tribool::False},
        {Sign::NonNegative, true, Tribool::True, Tribool::Unknown},
        {Sign::NonPositive, true, Tribool::Unknown, Tribool::False},
        {std::nullopt, true, Tribool::Unknown, Tribool::Unknown},
        {std::nullopt, false, Tribool::Unknown, Tribool::Unknown}};
    const auto power = make_power(detail::make_node<NumberNode>(BigInt(0)),
                                  detail::make_node<VariableNode>("z"));
    for (const auto& entry : cases) {
        AssumptionContext assumptions;
        if (entry.real) { ASSERT_TRUE(assumptions.assume_domain_checked("z", Domain::Real)); }
        if (entry.sign) { ASSERT_TRUE(assumptions.assume_sign_checked("z", *entry.sign)); }
        detail::AssumptionFacts facts(assumptions);
        ComputationContext context;
        auto complex = detail::query_definedness(power, facts, Domain::Complex, context);
        auto real = detail::query_definedness(power, facts, Domain::Real, context);
        ASSERT_TRUE(complex);
        ASSERT_TRUE(real);
        EXPECT_EQ(complex.value(), entry.complex_domain);
        EXPECT_EQ(real.value(), entry.real_domain);
        if (!entry.real) {
            auto projection = detail::domain_constraints(power, facts, Domain::Complex, context);
            ASSERT_TRUE(projection);
            EXPECT_FALSE(projection.value());
        }
    }
}
