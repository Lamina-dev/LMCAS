#include "test_assumption_interval_support.hpp"
#include "internal/interval_endpoint.hpp"
#include <limits>
#include <optional>
#include <string>
#include <vector>

static std::shared_ptr<SymbolicExpr> expr_from_node(
    std::shared_ptr<const SymbolicNode> node) {
    return detail::make_expression_ptr(std::move(node));
}

static void declare_expr_bounds(
    AssumptionContext &context, const std::string &variable,
    std::shared_ptr<SymbolicExpr> lower,
    std::shared_ptr<SymbolicExpr> upper) {
    Interval bounds{Endpoint::closed(std::move(lower)),
                    Endpoint::closed(std::move(upper))};
    ASSERT_TRUE(context.current_properties().declare_bounded(
        variable, Boundedness::Bounded, bounds));
}

static Rational exact_endpoint_value(const Endpoint &endpoint) {
    ComputationContext context;
    auto value = detail::comparable_endpoint(endpoint, context, "test.bounds");
    EXPECT_TRUE((value.has_value())) << "propagated endpoint is exactly comparable";
    if (!value)
        return Rational(0);
    EXPECT_TRUE((!value.value().algebraic && value.value().infinity == 0)) << "arithmetic propagation returns finite rational endpoints";
    return value.value().rational;
}

static void expect_exact_bounds(const std::optional<Interval> &bounds,
                                const Rational &lower, const Rational &upper) {
    ASSERT_TRUE((bounds.has_value())) << "exact interval propagation succeeds";
    if (!bounds)
        return;
    EXPECT_TRUE((exact_endpoint_value(bounds->lower) == lower)) << "exact lower endpoint";
    EXPECT_TRUE((exact_endpoint_value(bounds->upper) == upper)) << "exact upper endpoint";
}

static void expect_extreme_exact_enclosures(InferenceEngine &engine) {
    const Rational huge(BigInt("1" + std::string(400, '0')));
    auto huge_node = detail::node(SymbolicExpr::number(huge));
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(huge_node)), huge, huge);
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(detail::make_node<AddNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{huge_node, huge_node}))),
                        huge + huge, huge + huge);
    const double tiny = std::numeric_limits<double>::denorm_min();
    const Rational exact_tiny = Rational::from_double(tiny);
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(detail::make_node<MultiplyNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{make_num(tiny), make_num(tiny)}))),
                        exact_tiny * exact_tiny, exact_tiny * exact_tiny);
    const double maximum = std::numeric_limits<double>::max();
    const Rational exact_maximum = Rational::from_double(maximum);
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(detail::make_node<AddNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{make_num(maximum), make_num(maximum)}))),
                        exact_maximum + exact_maximum, exact_maximum + exact_maximum);
    std::shared_ptr<const SymbolicNode> deep = test_integer_node(1);
    for (std::size_t i = 0; i <= ResourceLimits{}.max_recursion_depth; ++i) {
        deep = detail::make_node<PowerNode>(deep, test_integer_node(2));
    }
    EXPECT_FALSE((engine.propagate_bounds(test_expression_from_node(deep)).has_value())) << "one shared recursion budget rejects the whole result, not a partial sum";
}

TEST(AssumptionIntervalCertified, ExactEnclosures) {
    AssumptionContext context;
    InferenceEngine engine(context);
    auto tenth = SymbolicExpr::number(Rational(1, 10));
    auto fifth = SymbolicExpr::number(Rational(1, 5));
    auto exact = engine.propagate_bounds(*tenth);
    auto approximate = engine.propagate_bounds(test_expression_from_node(make_num(0.1)));
    expect_exact_bounds(exact, Rational(1, 10), Rational(1, 10));
    EXPECT_TRUE((approximate.has_value())) << "binary64 point has exact bounds";
    if (exact && approximate) {
        auto same = IntervalUnion::from_single(*exact).intersect_checked(
            IntervalUnion::from_single(Interval::point(tenth)));
        EXPECT_TRUE((same && !same.value().is_empty())) << "exact tenth intersects itself";
        auto distinct = IntervalUnion::from_single(*exact).intersect_checked(
            IntervalUnion::from_single(*approximate));
        EXPECT_TRUE((distinct && distinct.value().is_empty())) << "binary64 .1 and exact 1/10 are different points";
    }
    auto add = detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(tenth), detail::node(fifth)});
    auto product = detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(tenth), detail::node(fifth)});
    auto quotient = detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(tenth),
                                                         detail::make_node<PowerNode>(detail::node(SymbolicExpr::number(Rational(3, 10))),
                                                                                      test_integer_node(-1))});
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(add)), Rational(3, 10), Rational(3, 10));
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(product)), Rational(1, 50), Rational(1, 50));
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(quotient)), Rational(1, 3), Rational(1, 3));
    declare_expr_bounds(context, "x", SymbolicExpr::number(Rational(-1, 10)), fifth);
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2)))),
                        Rational(0), Rational(1, 25));
    declare_expr_bounds(context, "expression_endpoints",
                        expr_from_node(detail::make_node<AddNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(1), test_integer_node(1)})),
                        expr_from_node(detail::make_node<AddNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(3), test_integer_node(2)})));
    expect_exact_bounds(engine.propagate_bounds(test_expression_from_node(detail::make_node<AddNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{
                                test_variable_node("expression_endpoints"), test_integer_node(1)}))),
                        Rational(3), Rational(6));
    expect_extreme_exact_enclosures(engine);
}

static void expect_contains_exact(const std::optional<Interval> &bounds,
                                  const std::shared_ptr<SymbolicExpr> &value) {
    ASSERT_TRUE((bounds.has_value())) << "algebraic enclosure succeeds";
    if (!bounds)
        return;
    auto intersection = IntervalUnion::from_single(*bounds).intersect_checked(
        IntervalUnion::from_single(Interval::point(value)));
    EXPECT_TRUE((intersection && !intersection.value().is_empty())) << "enclosure contains the exact algebraic value";
}

TEST(AssumptionIntervalCertified, AlgebraicEnclosuresAndOpenness) {
    auto sqrt_two = detail::make_node<FunctionNode>(FunctionNode::FuncType::Sqrt,
                                                    std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(2)});
    auto small = detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{sqrt_two,
                                                         detail::node(SymbolicExpr::number(Rational(-1414213, 1000000)))});
    AssumptionContext context;
    declare_expr_bounds(context, "u", expr_from_node(small), expr_from_node(small));
    InferenceEngine engine(context);
    auto reciprocal = detail::make_node<PowerNode>(test_variable_node("u"), test_integer_node(-1));
    expect_contains_exact(engine.propagate_bounds(test_expression_from_node(reciprocal)),
                          expr_from_node(detail::make_node<PowerNode>(small, test_integer_node(-1))));
    declare_expr_bounds(context, "v", expr_from_node(sqrt_two), expr_from_node(sqrt_two));
    for (bool multiply : {false, true}) {
        std::shared_ptr<const SymbolicNode> expression;
        std::shared_ptr<const SymbolicNode> expected;
        if (multiply) {
            expression = detail::make_node<MultiplyNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("v"), test_integer_node(2)});
            expected = detail::make_node<MultiplyNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{sqrt_two, test_integer_node(2)});
        } else {
            expression = detail::make_node<AddNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("v"), test_integer_node(1)});
            expected = detail::make_node<AddNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{sqrt_two, test_integer_node(1)});
        }
        expect_contains_exact(engine.propagate_bounds(test_expression_from_node(expression)), expr_from_node(expected));
    }
    Interval open{Endpoint::open(SymbolicExpr::number(0)), Endpoint::closed(SymbolicExpr::number(1))};
    EXPECT_TRUE((context.current_properties().declare_bounded("z", Boundedness::Bounded, open).has_value())) << "open interval declaration succeeds";
    auto original = engine.propagate_bounds(test_expression_from_node(test_variable_node("z")));
    EXPECT_TRUE((original && original->lower.is_open && !original->upper.is_open)) << "variable propagation preserves endpoint openness";
    EXPECT_FALSE((engine.propagate_bounds(test_expression_from_node(detail::make_node<PowerNode>(
                                              test_variable_node("z"), test_integer_node(-1))))
                      .has_value()))
        << "open zero endpoint has no finite reciprocal enclosure";

    auto logarithm = test_expression_from_node(detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                                               std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("z")}));
    EXPECT_TRUE((engine.infer_monotonicity(logarithm, "z", open) == Monotonicity::Increasing)) << "logarithm is increasing on an open-zero interval";
    const Rational tiny_negative(BigInt(-1), BigInt("1" + std::string(400, '0')));
    open.lower = Endpoint::open(SymbolicExpr::number(tiny_negative));
    EXPECT_TRUE((engine.infer_monotonicity(logarithm, "z", open) == Monotonicity::Unknown)) << "a tiny negative endpoint cannot round to an open zero";
}

TEST(AssumptionIntervalCertified, CertifiedArctangentEnclosures) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_domain("x", Domain::Real).has_value())) << "real declaration succeeds";
    InferenceEngine engine(context);
    auto atan_x = test_expression_from_node(detail::make_node<FunctionNode>(FunctionNode::FuncType::ArcTan,
                                                            std::vector<std::shared_ptr<const SymbolicNode>>{test_variable_node("x")}));
    const BigInt scale("1" + std::string(50, '0'));
    const BigInt digits("314159265358979323846264338327950288419716939937510");
    const Rational pi_upper(digits + BigInt(1), scale);
    auto global = engine.propagate_bounds(atan_x);
    EXPECT_TRUE((global && global->lower.is_open && global->upper.is_open)) << "global atan enclosure has open endpoints";
    if (global) {
        EXPECT_TRUE((exact_endpoint_value(global->lower) < -pi_upper / Rational(2))) << "global lower endpoint lies outside the pi/2 uncertainty interval";
        EXPECT_TRUE((exact_endpoint_value(global->upper) > pi_upper / Rational(2))) << "global upper endpoint lies outside the pi/2 uncertainty interval";
    }
    declare_expr_bounds(context, "x", SymbolicExpr::number(-1), SymbolicExpr::number(1));
    auto finite = engine.propagate_bounds(atan_x);
    ASSERT_TRUE((finite.has_value())) << "finite atan enclosure succeeds";
    if (finite) {
        EXPECT_TRUE((exact_endpoint_value(finite->lower) < -pi_upper / Rational(4))) << "finite lower endpoint encloses -pi/4";
        EXPECT_TRUE((exact_endpoint_value(finite->upper) > pi_upper / Rational(4))) << "finite upper endpoint encloses pi/4";
    }
    auto invalid = detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln,
                                                   std::vector<std::shared_ptr<const SymbolicNode>>{test_integer_node(-1)});
    auto sine = detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin,
                                                std::vector<std::shared_ptr<const SymbolicNode>>{invalid});
    EXPECT_FALSE((engine.propagate_bounds(test_expression_from_node(sine)).has_value())) << "sin bounds do not hide an undefined argument";
}
