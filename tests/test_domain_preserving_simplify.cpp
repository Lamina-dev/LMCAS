#include "test_common.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/facts_query.hpp"
#include "internal/assumption_facts.hpp"
#include "computation_context.hpp"
#include "expr.hpp"
#include "assumption_context.hpp"
#include "numeric_evaluation.hpp"
#include "internal/transcendental_solver_support.hpp"

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> normalize_node(const std::shared_ptr<const SymbolicNode> &node) {
    NormalizationVisitor visitor;
    node->accept(visitor);
    return visitor.get_result();
}

static bool is_bigint_value(const std::shared_ptr<const SymbolicNode> &node, const BigInt &expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    return number &&
           std::holds_alternative<BigInt>(number->value()) &&
           std::get<BigInt>(number->value()) == expected;
}

TEST(LmcasDomainPreservingSimplify, FractionalExpansionAndExactIdentities) {
    auto x = SymbolicExpr::variable("x");
    for (const auto &exponent : {SymbolicExpr::number(Rational(3, 2)),
                                 SymbolicExpr::number(1.5)}) {
        auto expanded = SymbolicExpr::power(x, exponent)->expand();
        auto evaluated = expanded->substitute("x", SymbolicExpr::number(4));
        {
            const double actual_value = evaluated->to_numeric();
            const double expected_value = 8.0;
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, 0.0);
        }
    }
    const BigInt denominator = BigInt(1) << 100;
    auto near_one = SymbolicExpr::number(Rational(denominator + BigInt(1), denominator));
    EXPECT_FALSE((near_one->is_one())) << "an exact rational differing from one below binary64 precision is not one";
}

using Node = std::shared_ptr<const SymbolicNode>;

static Node raw_domain_power(const Node &a, const Node &b) {
    return detail::make_node<PowerNode>(a, b);
}

static Node raw_domain_function(FunctionNode::FuncType type, const Node &a) {
    return detail::make_node<FunctionNode>(type, std::vector<Node>{a});
}

static void expect_node_domain(const Node &expression, Domain domain,
                               Tribool expected, const char *message) {
    ComputationContext context;
    const auto result = detail::query_definedness(expression, detail::no_facts(), domain, context);
    EXPECT_TRUE((result && result.value() == expected)) << message;
}

static void expect_total_function_domains(const Node &x, const Node &minus_one, const Node &zero) {
    expect_node_domain(raw_domain_function(FunctionNode::FuncType::ArcTan, x), Domain::Real, Tribool::True,
                       "arctangent is defined on the whole real independent-variable domain");
    expect_node_domain(raw_domain_function(FunctionNode::FuncType::ArcTan, minus_one), Domain::Complex, Tribool::True,
                       "a certified real arctangent argument avoids complex poles");
    ComputationContext context;
    auto imaginary = LMCAS::imaginary_unit();
    ASSERT_TRUE(imaginary);
    const auto complex_pole = detail::query_definedness(
        raw_domain_function(FunctionNode::FuncType::ArcTan, detail::node(imaginary.value())),
        detail::no_facts(), Domain::Complex, context);
    EXPECT_TRUE((complex_pole && complex_pole.value() != Tribool::True)) << "real totality must not certify the complex arctangent pole";
    expect_node_domain(raw_domain_function(FunctionNode::FuncType::Sgn, x), Domain::Real, Tribool::True,
                       "real sign is defined for every finite real independent-variable value");
    expect_node_domain(raw_domain_function(FunctionNode::FuncType::Sgn, zero), Domain::Real, Tribool::True,
                       "real sign is defined at zero despite its directional discontinuity");
}

static Tribool domain_comparison_value(const RelationalNode &comparison) {
    const auto left = test_numeric_eval(detail::make_expression_ptr(comparison.left()));
    const auto right = test_numeric_eval(detail::make_expression_ptr(comparison.right()));
    if (!left || !right || !std::isfinite(*left) || !std::isfinite(*right)) {
        return Tribool::Unknown;
    }
    bool truth;
    switch (comparison.op()) {
    case RelationOp::EQ: {
        truth = *left == *right;
        break;
    }
    case RelationOp::NEQ: {
        truth = *left != *right;
        break;
    }
    case RelationOp::LT: {
        truth = *left < *right;
        break;
    }
    case RelationOp::LEQ: {
        truth = *left <= *right;
        break;
    }
    case RelationOp::GT: {
        truth = *left > *right;
        break;
    }
    case RelationOp::GEQ: {
        truth = *left >= *right;
        break;
    }
    default: {
        return Tribool::Unknown;
    }
    }
    return truth ? Tribool::True : Tribool::False;
}

static Tribool domain_predicate_value(const Node &node);

static Tribool domain_logical_value(const LogicalNode &logical) {
    const auto left = domain_predicate_value(logical.left());
    const auto right = domain_predicate_value(logical.right());
    if (logical.op() == LogicalNode::Op::And) {
        if (left == Tribool::False || right == Tribool::False) {
            return Tribool::False;
        }
        return left == Tribool::True && right == Tribool::True ? Tribool::True : Tribool::Unknown;
    }
    if (logical.op() == LogicalNode::Op::Or) {
        if (left == Tribool::True || right == Tribool::True) {
            return Tribool::True;
        }
        return left == Tribool::False && right == Tribool::False ? Tribool::False : Tribool::Unknown;
    }
    return Tribool::Unknown;
}

static Tribool domain_predicate_value(const Node &node) {
    if (const auto *logical = dynamic_cast<const LogicalNode *>(node.get())) {
        return domain_logical_value(*logical);
    }
    const auto *comparison = dynamic_cast<const RelationalNode *>(node.get());
    if (!comparison) {
        return Tribool::Unknown;
    }
    return domain_comparison_value(*comparison);
}

static void expect_domain_projection(const Node &expression, int value, Tribool expected) {
    ComputationContext context;
    const auto result = detail::domain_constraints(expression, detail::no_facts(), Domain::Real, context);
    EXPECT_TRUE((result && result.value().has_value())) << "domain is completely expressible";
    if (!result || !result.value()) {
        return;
    }
    Tribool actual = Tribool::True;
    for (const auto &predicate : *result.value()) {
        auto substituted = predicate->substitute("x", SymbolicExpr::number(value));
        const auto truth = domain_predicate_value(detail::node(substituted));
        if (truth == Tribool::False) {
            actual = Tribool::False;
        } else if (truth == Tribool::Unknown && actual != Tribool::False) {
            actual = Tribool::Unknown;
        }
    }
    EXPECT_TRUE((actual == expected)) << "projected domain decides the original boundary correctly";
}

TEST(LmcasDomainPreservingSimplify, SharedDomainProjections) {
    const Node x = detail::make_node<VariableNode>("x");
    const Node zero = detail::make_node<NumberNode>(BigInt(0));
    const Node minus_one = detail::make_node<NumberNode>(BigInt(-1));
    const Node half = detail::make_node<NumberNode>(Rational(1, 2));
    const auto reciprocal = raw_domain_power(x, minus_one);
    const auto undefined = raw_domain_power(zero, minus_one);
    const auto &facts = detail::no_facts();
    expect_domain_projection(reciprocal, 0, Tribool::False);
    expect_domain_projection(reciprocal, -2, Tribool::True);
    expect_domain_projection(raw_domain_function(FunctionNode::FuncType::Ln, x), 0, Tribool::False);
    expect_domain_projection(raw_domain_power(x, half), 0, Tribool::True);
    expect_domain_projection(raw_domain_power(x, half), -1, Tribool::False);
    expect_domain_projection(undefined, 1, Tribool::False);
    ComputationContext context;
    auto unknown = detail::domain_constraints(
        raw_domain_power(x, detail::make_node<VariableNode>("y")), facts, Domain::Real, context);
    EXPECT_TRUE((unknown && !unknown.value())) << "unknown exponent reality cannot be replaced by a narrower domain";
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    auto stopped = detail::domain_constraints(reciprocal, facts, Domain::Real, exhausted);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::ResourceLimit)) << "domain projection propagates the caller's exhausted budget";
}

TEST(LmcasDomainPreservingSimplify, SharedDomainQueries) {
    const auto &facts = detail::no_facts();
    ComputationContext context;
    const Node x = detail::make_node<VariableNode>("x");
    const Node zero = detail::make_node<NumberNode>(BigInt(0));
    const Node one = detail::make_node<NumberNode>(BigInt(1));
    const Node minus_one = detail::make_node<NumberNode>(BigInt(-1));
    const Node half = detail::make_node<NumberNode>(Rational(1, 2));
    expect_node_domain(x, Domain::Real, Tribool::True, "finite placeholders need no bindings");
    expect_node_domain(detail::make_node<VariableNode>("i"), Domain::Real, Tribool::True,
                       "lowercase i remains an ordinary finite variable");
    expect_node_domain(detail::make_node<VariableNode>("I"), Domain::Real, Tribool::True,
                       "uppercase I remains an ordinary finite variable");
    auto imaginary_unit_value = LMCAS::imaginary_unit();
    ASSERT_TRUE(imaginary_unit_value);
    expect_node_domain(detail::node(imaginary_unit_value.value()), Domain::Real, Tribool::False,
                       "the explicit imaginary unit is nonreal");
    expect_node_domain(raw_domain_power(zero, zero), Domain::Real, Tribool::False, "real zero to zero is undefined");
    expect_node_domain(raw_domain_power(zero, zero), Domain::Complex, Tribool::True, "complex integer zero power retains its convention");
    expect_total_function_domains(x, minus_one, zero);
    const auto reciprocal = raw_domain_power(x, minus_one);
    expect_node_domain(reciprocal, Domain::Real, Tribool::Unknown, "unknown reciprocal retains its pole");
    const auto undefined = raw_domain_power(zero, minus_one);
    expect_node_domain(undefined, Domain::Complex, Tribool::False, "zero reciprocal has no finite complex value");
    expect_node_domain(detail::make_node<MultiplyNode>(std::vector<Node>{zero, undefined}),
                       Domain::Real, Tribool::False, "zero multiplication does not skip undefined operands");
    expect_node_domain(detail::make_node<MultiplyNode>(std::vector<Node>{zero, reciprocal}),
                       Domain::Complex, Tribool::Unknown, "zero multiplication does not skip uncertain operands");
    const auto negative_log = raw_domain_function(FunctionNode::FuncType::Ln, minus_one);
    expect_node_domain(negative_log, Domain::Real, Tribool::False, "negative real logarithm is undefined");
    expect_node_domain(negative_log, Domain::Complex, Tribool::True, "negative logarithm has a finite complex principal value");
    const auto invalid_nonzero = detail::query_nonzero_value(undefined, facts, Domain::Real, context);
    EXPECT_TRUE((invalid_nonzero && invalid_nonzero.value() == Tribool::Unknown)) << "undefinedness is neither a zero nor a nonzero value";
    const auto variable_real = detail::query_real_value(x, facts, context);
    EXPECT_TRUE((variable_real && variable_real.value() == Tribool::Unknown)) << "placeholder definedness does not create a global real-valued fact";
    const auto imaginary = detail::make_node<ComplexNode>(zero, one);
    expect_node_domain(imaginary, Domain::Real, Tribool::False, "nonreal complex value is not in the real domain");
    expect_node_domain(imaginary, Domain::Complex, Tribool::True, "explicit finite complex components are defined");
    const auto nonreal_product = detail::make_node<MultiplyNode>(std::vector<Node>{
        detail::make_node<NumberNode>(BigInt(2)), raw_domain_power(minus_one, half)});
    const auto product_real = detail::query_real_value(nonreal_product, facts, context);
    EXPECT_TRUE((product_real && product_real.value() == Tribool::False)) << "nonzero real scaling cannot make a nonreal principal root real";
    auto wrong_domain = detail::query_definedness(x, facts, Domain::Integer, context);
    EXPECT_TRUE((!wrong_domain && wrong_domain.error().code == CasErrc::InvalidArgument)) << "value-domain queries reject variable-assumption domains";
}

static void expect_real_value(const ExprPtr &expression, const NumericBindings &bindings,
                              std::optional<double> expected) {
    auto actual = evaluate_numeric(*expression, bindings);
    if (expected) {
        ASSERT_TRUE((actual.has_value())) << "legal expression remains evaluable";
        if (actual) {
            const double actual_value = actual.value().value;
            const double expected_value = *expected;
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, 1e-12);
        }
    } else {
        EXPECT_TRUE((!actual && actual.error().code == CasErrc::DomainError)) << "undefined expression remains a domain error";
    }
}

static void expect_partial_arithmetic_transformations(
    const ExprPtr &x, const ExprPtr &zero, const ExprPtr &one) {
    auto inverse = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto cancellation = SymbolicExpr::add(inverse,
                                          SymbolicExpr::multiply(SymbolicExpr::number(-1), inverse));
    const std::vector<ExprPtr> expressions{
        SymbolicExpr::multiply(zero, inverse), cancellation,
        SymbolicExpr::power(x, zero), SymbolicExpr::power(one, inverse)};
    for (std::size_t index = 0; index < expressions.size(); ++index) {
        const auto &expression = expressions[index];
        for (const auto &transformed : {expression, expression->simplify(), expression->expand()}) {
            expect_real_value(transformed, {{"x", 0}}, std::nullopt);
            expect_real_value(transformed, {{"x", 2}}, index < 2 ? 0.0 : 1.0);
            expect_real_value(transformed->substitute("x", zero)->simplify(), {}, std::nullopt);
            expect_real_value(transformed->substitute("x", SymbolicExpr::number(2))->simplify(),
                              {}, index < 2 ? 0.0 : 1.0);
        }
    }
}

static void expect_constant_domain_transformations(const ExprPtr &zero) {
    for (const auto &expression : {SymbolicExpr::ln(zero), SymbolicExpr::power(zero, zero)}) {
        for (const auto &transformed : {expression, expression->simplify(), expression->expand()}) {
            expect_real_value(transformed, {}, std::nullopt);
        }
    }
    auto zero_power = SymbolicExpr::power(zero, zero);
    for (const auto &transformed : {zero_power, zero_power->simplify(), zero_power->expand()}) {
        auto value = eval_complex(*transformed);
        EXPECT_TRUE((value && value.value().real.value == 1 && value.value().imag.value == 0)) << "complex integer zero power keeps its separate convention";
    }
}

static void expect_factory_operand_domains(const ExprPtr &zero, const ExprPtr &one) {
    const auto raw_zero = detail::node(zero);
    const auto undefined_exponent = detail::make_node<PowerNode>(raw_zero, raw_zero);
    const auto nested_zero = detail::make_expression_ptr(
        detail::make_node<PowerNode>(raw_zero, undefined_exponent));
    for (const auto &transformed : {nested_zero, nested_zero->simplify(), nested_zero->expand(),
                                    detail::make_expression_ptr(SymbolicFactory::create_power(raw_zero, undefined_exponent))}) {
        expect_real_value(transformed, {}, std::nullopt);
    }
    bool rejected = false;
    try {
        (void)SymbolicFactory::create_multiply({raw_zero, nullptr});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "zero does not skip validation of later factory operands";
    auto factory = detail::make_expression_ptr(SymbolicFactory::create_multiply({raw_zero, SymbolicFactory::create_power(raw_zero, detail::node(SymbolicExpr::number(-1)))}));
    expect_real_value(factory, {}, std::nullopt);
    expect_real_value(factory->simplify(), {}, std::nullopt);
    auto undefined_component = SymbolicFactory::create_complex(detail::node(one),
                                                               SymbolicFactory::create_power(raw_zero, detail::node(SymbolicExpr::number(-1))));
    auto complex_zero_power = detail::make_expression_ptr(SymbolicFactory::create_power(
        undefined_component, raw_zero));
    EXPECT_FALSE((complex_zero_power->is_one())) << "a nonzero real component cannot hide an undefined imaginary component";
}

TEST(LmcasDomainPreservingSimplify, DomainPreservingTransformations) {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    expect_partial_arithmetic_transformations(x, zero, one);
    expect_constant_domain_transformations(zero);
    auto cancelled_polynomial = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    EXPECT_TRUE((cancelled_polynomial->simplify()->is_zero())) << "polynomial cancellation remains valid";
    EXPECT_TRUE((SymbolicExpr::multiply(zero, x)->simplify()->is_zero())) << "zero times a finite variable is zero";
    auto polynomial = SymbolicExpr::power(SymbolicExpr::add(x, one), SymbolicExpr::number(3));
    expect_real_value(polynomial->expand(), {{"x", 2}}, 27.0);
    expect_factory_operand_domains(zero, one);
}

TEST(LmcasDomainPreservingSimplify, PowerAndLogarithmBranches) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto product = SymbolicExpr::multiply(x, y);
    auto root = SymbolicExpr::power(product, SymbolicExpr::number(Rational(1, 2)));
    for (const auto &transformed : {root, root->simplify(), root->expand()}) {
        expect_real_value(transformed, {{"x", -1}, {"y", -4}}, 2.0);
        expect_real_value(transformed, {{"x", 1}, {"y", 4}}, 2.0);
        expect_real_value(transformed, {{"x", 0}, {"y", 4}}, 0.0);
        expect_real_value(transformed, {{"x", -1}, {"y", 4}}, std::nullopt);
    }
    auto negative_power = SymbolicExpr::power(product, SymbolicExpr::number(-2));
    expect_real_value(negative_power->simplify(), {{"x", 0}, {"y", 4}}, std::nullopt);
    expect_real_value(negative_power->simplify(), {{"x", -1}, {"y", 4}}, 1.0 / 16);
    auto nested = SymbolicExpr::power(
        SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 2))), SymbolicExpr::number(2));
    expect_real_value(nested->simplify(), {{"x", -1}}, std::nullopt);
    expect_real_value(nested->simplify(), {{"x", 4}}, 4.0);
    auto logarithm = SymbolicExpr::ln(SymbolicExpr::power(x, SymbolicExpr::number(2)));
    for (const auto &transformed : {logarithm, logarithm->simplify(), logarithm->expand()}) {
        expect_real_value(transformed, {{"x", -2}}, std::log(4.0));
        expect_real_value(transformed, {{"x", 0}}, std::nullopt);
    }
    auto exponential_log = SymbolicExpr::ln(SymbolicExpr::exp(x));
    EXPECT_TRUE((detail::node(exponential_log->simplify())->equals(*detail::node(exponential_log)))) << "unknown complex exponent cannot escape the logarithm principal branch";
    auto imaginary = detail::make_expression_ptr(detail::make_node<ComplexNode>(
        detail::node(SymbolicExpr::number(0)), detail::node(SymbolicExpr::number(4))));
    auto complex_log = SymbolicExpr::ln(SymbolicExpr::exp(imaginary));
    EXPECT_TRUE((dynamic_cast<const FunctionNode *>(detail::node(complex_log->simplify()).get()) != nullptr)) << "explicit nonreal exponent is not returned from ln(exp(z))";

    auto matrix = SymbolicExpr::matrix({{SymbolicExpr::number(1), x}, {SymbolicExpr::number(0), SymbolicExpr::number(1)}});
    auto matrix_product = detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(matrix), detail::node(matrix)});
    auto matrix_power = detail::make_node<PowerNode>(matrix_product, detail::node(SymbolicExpr::number(-1)));
    EXPECT_FALSE((dynamic_cast<const MultiplyNode *>(normalize_node(matrix_power).get()) != nullptr)) << "noncommuting matrix powers are not distributed as scalar powers";
}

TEST(LmcasDomainPreservingSimplify, AssumedTransformationDomains) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_sign("x", Sign::NonZero).has_value())) << "nonzero declaration succeeds";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "context accepts assumptions";
    auto power = SymbolicExpr::power(x, SymbolicExpr::number(0));
    auto simplified = simplify(power, context);
    auto expanded = expand(power, context);
    EXPECT_TRUE((simplified && simplified.value()->is_one())) << "nonzero assumption permits x^0 simplification";
    EXPECT_TRUE((expanded && expanded.value()->is_one())) << "expand applies assumptions after conservative expansion";
    auto invalid = SymbolicExpr::power(SymbolicExpr::ln(SymbolicExpr::number(-1)), SymbolicExpr::number(0));
    auto preserved = simplify(invalid, context);
    ASSERT_TRUE((preserved.has_value())) << "simplification is not a domain rejection API";
    if (preserved)
        expect_real_value(preserved.value(), {}, std::nullopt);
    EXPECT_TRUE((assumptions->assume_sign("y", Sign::Positive).has_value())) << "positive declaration succeeds";
    EXPECT_TRUE((assumptions->assume_sign("x", Sign::Positive).has_value())) << "positive strengthens nonzero";
    auto root = SymbolicExpr::power(SymbolicExpr::multiply(x, y), SymbolicExpr::number(Rational(1, 2)));
    auto distributed = simplify(root, context);
    EXPECT_TRUE((distributed && dynamic_cast<const MultiplyNode *>(detail::node(distributed.value()).get()))) << "positive factors permit real fractional power distribution";
    if (distributed)
        expect_real_value(distributed.value(), {{"x", 4}, {"y", 9}}, 6.0);
    auto logarithm = simplify(SymbolicExpr::ln(SymbolicExpr::exp(y)), context);
    EXPECT_TRUE((logarithm && detail::node(logarithm.value())->equals(*detail::node(y)))) << "positive real exponent permits ln(exp(y)) reduction";
    auto base_two = SymbolicExpr::log(SymbolicExpr::power(SymbolicExpr::number(2), y), SymbolicExpr::number(2));
    auto reduced = simplify(base_two, context);
    EXPECT_TRUE((reduced && detail::node(reduced.value())->equals(*detail::node(y)))) << "positive nonunit base and real exponent permit inverse logarithm";
}

static void expect_conditional_projection_boundaries(
    const std::vector<std::shared_ptr<SymbolicExpr>> &conditions,
    AssumptionContext &assumptions, bool negative_allowed) {
    for (int value : {-2, 0, 2}) {
        Tribool accepted = Tribool::True;
        for (const auto &condition : conditions) {
            auto bound = condition->substitute("x", SymbolicExpr::number(value));
            const auto result = assumptions.evaluate_condition(*bound);
            if (result == Tribool::False) {
                accepted = Tribool::False;
                break;
            }
            if (result == Tribool::Unknown) {
                accepted = Tribool::Unknown;
            }
        }
        const bool expected = value > 0 || (value < 0 && negative_allowed);
        EXPECT_TRUE((accepted == (expected ? Tribool::True : Tribool::False))) << "projected predicates preserve the zero exclusion and logarithm branch";
    }
}

TEST(LmcasDomainPreservingSimplify, ConditionalDomainProjection) {
    AssumptionContext assumptions;
    EXPECT_TRUE((assumptions.assume_domain("x", Domain::Real).has_value())) << "real variable declared";
    detail::AssumptionFacts facts(assumptions);
    ComputationContext context;
    auto x = SymbolicExpr::variable("x");
    auto reciprocal = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>, bool>> cases{
        {SymbolicExpr::power(SymbolicExpr::number(1), reciprocal), true},
        {SymbolicExpr::power(SymbolicExpr::number(1), SymbolicExpr::ln(x)), false},
        {SymbolicExpr::ln(SymbolicExpr::power(x, SymbolicExpr::number(-2))), true}};
    for (const auto &[expression, negative_allowed] : cases) {
        auto domain = detail::query_definedness(detail::node(expression), facts, Domain::Real, context);
        auto real = detail::query_real_value(detail::node(expression), facts, context);
        EXPECT_TRUE((domain && domain.value() == Tribool::Unknown)) << "the original expression is not defined everywhere";
        EXPECT_TRUE((real && real.value() == Tribool::Unknown)) << "conditional realness is not an ordinary value proof";
        auto projected = detail::domain_constraints(detail::node(expression), facts, Domain::Real, context);
        EXPECT_TRUE((projected && projected.value().has_value())) << "the original domain has an expressible conjunction";
        if (!projected || !projected.value()) {
            continue;
        }
        expect_conditional_projection_boundaries(*projected.value(), assumptions, negative_allowed);
    }
}

static void expect_zero_power_projection(
    const Node& power, const FactsQuery& facts, Domain domain,
    const std::vector<std::pair<int, Tribool>>& bindings) {
    ComputationContext context;
    auto projection = detail::domain_constraints(power, facts, domain, context);
    ASSERT_TRUE(projection);
    ASSERT_TRUE(projection.value());
    for (const auto& [value, expected] : bindings) {
        Tribool accepted = Tribool::True;
        for (const auto& condition : *projection.value()) {
            auto bound = condition->substitute("z", SymbolicExpr::number(value));
            const auto truth = domain_predicate_value(detail::node(bound));
            if (truth != Tribool::True) { accepted = truth; break; }
        }
        EXPECT_EQ(accepted, expected) << value;
    }
}

TEST(LmcasDomainPreservingSimplify, ZeroBaseConditionalValueFacts) {
    const auto power = raw_domain_power(detail::make_node<NumberNode>(BigInt(0)),
                                        detail::make_node<VariableNode>("z"));
    AssumptionContext nonnegative;
    ASSERT_TRUE(nonnegative.assume_sign_checked("z", Sign::NonNegative));
    detail::AssumptionFacts nonnegative_facts(nonnegative);
    ComputationContext context;
    auto nonzero = detail::query_nonzero_value(power, nonnegative_facts, Domain::Complex, context);
    auto positive = detail::query_positive_value(power, nonnegative_facts, context);
    auto real = detail::query_real_value(power, nonnegative_facts, context);
    auto nonnegative_value = detail::query_nonnegative_value(power, nonnegative_facts, context);
    ASSERT_TRUE(nonzero);
    ASSERT_TRUE(positive);
    ASSERT_TRUE(real);
    ASSERT_TRUE(nonnegative_value);
    EXPECT_EQ(nonzero.value(), Tribool::Unknown);
    EXPECT_EQ(positive.value(), Tribool::Unknown);
    EXPECT_EQ(real.value(), Tribool::True);
    EXPECT_EQ(nonnegative_value.value(), Tribool::True);
    AssumptionContext nonpositive;
    ASSERT_TRUE(nonpositive.assume_sign_checked("z", Sign::NonPositive));
    detail::AssumptionFacts nonpositive_facts(nonpositive);
    expect_zero_power_projection(power, nonpositive_facts, Domain::Complex,
                                 {{-1, Tribool::False}, {0, Tribool::True}});
    AssumptionContext real_exponent;
    ASSERT_TRUE(real_exponent.assume_domain_checked("z", Domain::Real));
    detail::AssumptionFacts real_facts(real_exponent);
    expect_zero_power_projection(power, real_facts, Domain::Complex,
                                 {{-1, Tribool::False}, {0, Tribool::True}, {1, Tribool::True}});
    expect_zero_power_projection(power, real_facts, Domain::Real,
                                 {{-1, Tribool::False}, {0, Tribool::False}, {1, Tribool::True}});
}

TEST(LmcasDomainPreservingSimplify, SymbolicZeroPowerSolutionFinalization) {
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_sign_checked("z", Sign::Zero));
    const auto power = raw_domain_power(detail::make_node<NumberNode>(BigInt(0)),
                                        detail::make_node<VariableNode>("z"));
    auto equation = detail::make_expression_ptr(detail::make_node<RelationalNode>(
        power, detail::make_node<NumberNode>(BigInt(1)), RelationOp::EQ));
    for (const auto domain : {Domain::Complex, Domain::Real}) {
        ComputationContext context;
        ASSERT_TRUE(context.set_assumptions(assumptions));
        FiniteSolutions candidates{{FiniteSolution{SymbolicExpr::number(0), 1, {}}}};
        auto result = detail::finalize_solution_set(equation, "z", std::move(candidates),
                                                     context, domain, SolveOptions{});
        ASSERT_TRUE(result) << result.error().operation << ": " << result.error().message;
        if (domain == Domain::Complex) {
            ASSERT_TRUE(std::holds_alternative<FiniteSolutions>(result.value()));
            ASSERT_EQ(std::get<FiniteSolutions>(result.value()).values.size(), 1u);
            EXPECT_TRUE(test_proved_equivalent(
                std::get<FiniteSolutions>(result.value()).values[0].value, SymbolicExpr::number(0)));
        } else {
            EXPECT_TRUE(std::holds_alternative<EmptySolutions>(result.value()));
        }
    }
}

TEST(LmcasDomainPreservingSimplify, SimplifyRetainsProductsOverSums) {
    auto x = LMCAS::detail::make_node<VariableNode>("x");
    auto y = LMCAS::detail::make_node<VariableNode>("y");
    auto y_plus_one = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{y, LMCAS::detail::make_node<NumberNode>(BigInt(1))});
    auto product = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, y_plus_one});
    auto simplified_product = normalize_node(product);
    EXPECT_TRUE((std::dynamic_pointer_cast<const MultiplyNode>(simplified_product) != nullptr)) << "simplify keeps x*(y+1) as a product";
}

TEST(LmcasDomainPreservingSimplify, ZeroExponentRequiresNonzeroBase) {
    auto unknown_zero_power = LMCAS::detail::make_node<PowerNode>(
        LMCAS::detail::make_node<VariableNode>("x"), LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto simplified_unknown_zero_power = normalize_node(unknown_zero_power);
    EXPECT_TRUE((std::dynamic_pointer_cast<const PowerNode>(simplified_unknown_zero_power) != nullptr)) << "simplify keeps x^0 without a nonzero assumption";
    auto numeric_zero_power = LMCAS::detail::make_node<PowerNode>(
        LMCAS::detail::make_node<NumberNode>(BigInt(2)), LMCAS::detail::make_node<NumberNode>(BigInt(0)));
    auto simplified_numeric_zero_power = normalize_node(numeric_zero_power);
    EXPECT_TRUE((is_bigint_value(simplified_numeric_zero_power, 1))) << "simplify folds 2^0 because the base is proved nonzero";
}

TEST(LmcasDomainPreservingSimplify, InverseCancellationIsConditional) {
    auto inv_x = LMCAS::detail::make_node<PowerNode>(
        LMCAS::detail::make_node<VariableNode>("x"), LMCAS::detail::make_node<NumberNode>(BigInt(-1)));
    auto x_times_inv_x = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::make_node<VariableNode>("x"), inv_x});
    auto simplified_inverse = normalize_node(x_times_inv_x);
    EXPECT_FALSE((simplified_inverse->is_one())) << "simplify does not turn x*x^-1 into 1 without x != 0";
}

TEST(LmcasDomainPreservingSimplify, PositiveIntegerExponentMerging) {
    auto x_squared = LMCAS::detail::make_node<PowerNode>(
        LMCAS::detail::make_node<VariableNode>("x"), LMCAS::detail::make_node<NumberNode>(BigInt(2)));
    auto x_times_x_squared = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::make_node<VariableNode>("x"), x_squared});
    auto simplified_positive_merge = normalize_node(x_times_x_squared);
    auto merged_power = std::dynamic_pointer_cast<const PowerNode>(simplified_positive_merge);
    EXPECT_TRUE((merged_power != nullptr)) << "x*x^2 merges to a power";
    EXPECT_TRUE((merged_power && is_bigint_value(merged_power->exponent(), 3))) << "x*x^2 merges to x^3";
}
