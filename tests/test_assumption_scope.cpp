
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "interval.hpp"
#include "query_interface.hpp"
#include "inference_engine.hpp"
#include <stdexcept>
#include <vector>
#include <string>

using namespace LMCAS;

static const std::vector<Domain> ALL_DOMAINS = {
    Domain::Complex, Domain::Real, Domain::Algebraic, Domain::Rational,
    Domain::Integer, Domain::Natural, Domain::PositiveInt};

static const std::vector<Sign> ALL_SIGNS = {
    Sign::Positive, Sign::Negative, Sign::NonNegative,
    Sign::NonPositive, Sign::Zero, Sign::NonZero};

static const std::vector<std::string> TEST_SYMBOLS = {
    "x", "y", "alpha", "longVar123", "a_b"};

static std::string domain_name(Domain d) {
    switch (d) {
    case Domain::Complex: {
        return "Complex";
    }
    case Domain::Real: {
        return "Real";
    }
    case Domain::Algebraic: {
        return "Algebraic";
    }
    case Domain::Rational: {
        return "Rational";
    }
    case Domain::Integer: {
        return "Integer";
    }
    case Domain::Natural: {
        return "Natural";
    }
    case Domain::PositiveInt: {
        return "PositiveInt";
    }
    }
    return "?";
}

static std::string sign_name(Sign s) {
    switch (s) {
    case Sign::Positive:
        return "Positive";
    case Sign::Negative:
        return "Negative";
    case Sign::NonNegative:
        return "NonNegative";
    case Sign::NonPositive:
        return "NonPositive";
    case Sign::Zero:
        return "Zero";
    case Sign::NonZero:
        return "NonZero";
    }
    return "?";
}

TEST(AssumptionScope, DomainRoundtrip) {
    // For each domain, declare in child scope, verify visible in child,
    // then pop and verify not visible (reverts to default Complex).
    for (const auto &sym : TEST_SYMBOLS) {
        for (Domain d : ALL_DOMAINS) {
            if (d == Domain::Complex)
                continue; // Complex is default, skip

            AssumptionContext ctx;
            // Before push: domain should be Complex (default)
            Domain before = ctx.get_domain(sym);
            EXPECT_TRUE((before == Domain::Complex)) << sym + " domain is Complex before push";

            ctx.push();
            EXPECT_TRUE(ctx.assume_domain(sym, d).has_value());

            // In child scope: domain should be at least d
            EXPECT_TRUE((ctx.has_domain(sym, d))) << sym + " has " + domain_name(d) + " in child scope";

            EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

            // After pop: domain should be back to Complex
            Domain after = ctx.get_domain(sym);
            EXPECT_TRUE((after == before)) << sym + " domain restored to Complex after pop (was " +
                                                  domain_name(d) + " in child)";
        }
    }
}

TEST(AssumptionScope, SignRoundtrip) {
    for (const auto &sym : TEST_SYMBOLS) {
        for (Sign s : ALL_SIGNS) {
            AssumptionContext ctx;
            // Before push: no signs
            auto signs_before = ctx.get_signs(sym);
            EXPECT_TRUE((signs_before.empty())) << sym + " has no signs before push";

            ctx.push();
            EXPECT_TRUE(ctx.assume_sign(sym, s).has_value());

            // In child scope: sign should be present
            EXPECT_TRUE((ctx.has_sign(sym, s))) << sym + " has " + sign_name(s) + " in child scope";

            EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

            // After pop: signs should be empty again
            auto signs_after = ctx.get_signs(sym);
            EXPECT_TRUE((signs_after.empty())) << sym + " signs empty after pop (was " +
                                                      sign_name(s) + " in child)";
        }
    }
}

TEST(AssumptionScope, ParityRoundtrip) {
    AssumptionContext ctx;

    // Before push: parity is Unknown
    EXPECT_TRUE((ctx.get_parity("x") == Parity::Unknown)) << "x parity Unknown before push";

    ctx.push();
    EXPECT_TRUE((ctx.current_properties().declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((ctx.get_parity("x") == Parity::Even)) << "x parity Even in child scope";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.get_parity("x") == Parity::Unknown)) << "x parity restored to Unknown after pop";
}

TEST(AssumptionScope, BoundednessRoundtrip) {
    AssumptionContext ctx;

    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Unknown)) << "x boundedness Unknown before push";

    ctx.push();
    EXPECT_TRUE((ctx.current_properties().declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";
    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Bounded)) << "x boundedness Bounded in child scope";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Unknown)) << "x boundedness restored to Unknown after pop";
}

TEST(AssumptionScope, RelationRoundtrip) {
    AssumptionContext ctx;

    // Create a simple relation: x > 0
    auto var_x = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<VariableNode>("x"));
    auto zero = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(0)));

    // Before push: no relations, no sign for x
    EXPECT_TRUE((ctx.current_relations().get_relations().empty())) << "No relations before push";
    EXPECT_FALSE((ctx.has_sign("x", Sign::Positive))) << "x not Positive before push";

    ctx.push();
    EXPECT_TRUE((ctx.current_relations().add_relation(*var_x, *zero,
                                                      RelationalNode::Op::GT, ctx.current_properties())
                     .has_value()))
        << "relation insertion succeeds";

    // In child scope: relation present, x is Positive
    EXPECT_TRUE((ctx.has_sign("x", Sign::Positive))) << "x is Positive in child scope (from relation x > 0)";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    // After pop: sign should be gone
    EXPECT_FALSE((ctx.has_sign("x", Sign::Positive))) << "x not Positive after pop";
    EXPECT_TRUE((ctx.current_relations().get_relations().empty())) << "No relations after pop";
}

TEST(AssumptionScope, MultipleDeclarationsRoundtrip) {
    AssumptionContext ctx;

    // Set up parent state
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    // Record parent state
    Domain x_domain_before = ctx.get_domain("x");
    auto y_signs_before = ctx.get_signs("y");
    auto z_signs_before = ctx.get_signs("z");
    Parity z_parity_before = ctx.get_parity("z");

    ctx.push();

    // Make multiple declarations in child
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE(ctx.assume_sign("z", Sign::Negative).has_value());
    EXPECT_TRUE((ctx.current_properties().declare_parity("z", Parity::Odd).has_value())) << "parity declaration succeeds";

    // Verify child state
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer in child";
    EXPECT_TRUE((ctx.has_sign("z", Sign::Negative))) << "z is Negative in child";
    EXPECT_TRUE((ctx.get_parity("z") == Parity::Odd)) << "z is Odd in child";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    // Verify all restored
    EXPECT_TRUE((ctx.get_domain("x") == x_domain_before)) << "x domain restored after pop";
    EXPECT_TRUE((ctx.get_signs("y") == y_signs_before)) << "y signs unchanged after pop";
    EXPECT_TRUE((ctx.get_signs("z") == z_signs_before)) << "z signs restored (empty) after pop";
    EXPECT_TRUE((ctx.get_parity("z") == z_parity_before)) << "z parity restored (Unknown) after pop";
}

TEST(AssumptionScope, ParentDeclarationsSurvivePop) {
    AssumptionContext ctx;

    // Declare in root scope
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    ctx.push();
    // Child scope doesn't touch x
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    // x should still have its root declarations
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x domain Integer survives child push/pop";
    EXPECT_TRUE((ctx.has_sign("x", Sign::Positive))) << "x sign Positive survives child push/pop";
}

TEST(AssumptionScope, NestedPushPopRoundtrip) {
    AssumptionContext ctx;

    // Root: x is Real
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    ctx.push(); // depth 2
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer at depth 2";

    ctx.push(); // depth 3
    EXPECT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());
    EXPECT_TRUE((ctx.get_domain("x") == Domain::PositiveInt)) << "x is PositiveInt at depth 3";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds"; /**< 返回深度 2。 */
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer after popping depth 3";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds"; /**< 返回深度 1 的根作用域。 */
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Real)) << "x is Real after popping depth 2";
}

TEST(AssumptionScope, DepthChanges) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.depth() == 1)) << "Initial depth is 1";

    ctx.push();
    EXPECT_TRUE((ctx.depth() == 2)) << "Depth is 2 after push";

    ctx.push();
    EXPECT_TRUE((ctx.depth() == 3)) << "Depth is 3 after second push";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.depth() == 2)) << "Depth is 2 after pop";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.depth() == 1)) << "Depth is 1 after second pop";
}

static void expect_scalar_property_readthrough(
    AssumptionContext &ctx, QueryInterface &query, const ExprPtr &n, const ExprPtr &t) {
    EXPECT_TRUE((query.query_positive(*n).value() == Tribool::True)) << "positive survives empty push";
    EXPECT_TRUE((query.query_real(*n).value() == Tribool::True)) << "real survives empty push";
    EXPECT_TRUE((query.query_integer(*n).value() == Tribool::True)) << "integer survives empty push";
    EXPECT_TRUE((ctx.has_domain("n", Domain::Rational))) << "rational survives empty push";
    EXPECT_TRUE((ctx.has_domain("n", Domain::Natural))) << "natural survives empty push";
    EXPECT_TRUE((query.query_algebraic(*n).value() == Tribool::True)) << "algebraic survives empty push";
    EXPECT_TRUE((query.query_transcendental(*t).value() == Tribool::True)) << "transcendence survives empty push";
}

static void expect_extended_property_shadowing(
    AssumptionContext &ctx, QueryInterface &query, const ExprPtr &matrix,
    const ExprPtr &f, const SymbolicExpr &period, const Interval &interval) {
    EXPECT_TRUE((ctx.current_properties().declare_definiteness("A", Definiteness::NegativeDefinite).has_value())) << "child definiteness overrides";
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("f", Finiteness::Divergent).has_value())) << "child finiteness overrides";
    EXPECT_TRUE((query.query_positive_definite(*matrix).value() == Tribool::False)) << "child definiteness is visible";
    EXPECT_TRUE((query.query_finite(*f).value() == Tribool::False)) << "child finiteness is visible";
    auto child_period = detail::expression_from_node(
        detail::make_node<NumberNode>(BigInt(7)));
    EXPECT_TRUE((ctx.current_properties().declare_periodic("f", "x", child_period).has_value())) << "child period overrides";
    EXPECT_TRUE((ctx.current_properties().declare_monotonicity("f", "x", interval, Monotonicity::Decreasing).has_value())) << "child monotonicity overrides";
    auto overridden_period = query.get_period(*f, "x");
    EXPECT_TRUE((overridden_period && overridden_period.value() &&
                 detail::node(*overridden_period.value())->equals(*detail::node(child_period))))
        << "child period visible";
    EXPECT_TRUE((ctx.get_monotonicity_checked("f", "x", interval).value() == Monotonicity::Decreasing)) << "child interval declaration wins";
    EXPECT_TRUE((ctx.pop().has_value())) << "child scope popped";
    EXPECT_TRUE((query.query_positive_definite(*matrix).value() == Tribool::True)) << "parent definiteness restored";
    EXPECT_TRUE((query.query_finite(*f).value() == Tribool::True)) << "parent finiteness restored";
    auto restored_period = query.get_period(*f, "x");
    EXPECT_TRUE((restored_period && restored_period.value() &&
                 detail::node(*restored_period.value())->equals(*detail::node(period))))
        << "parent period restored";
    EXPECT_TRUE((ctx.get_monotonicity_checked("f", "x", interval).value() == Monotonicity::Increasing)) << "parent monotonicity restored";
}

TEST(AssumptionScope, ExtendedPropertyReadthrough) {
    AssumptionContext ctx;
    auto n = detail::make_expression_ptr(detail::make_node<VariableNode>("n"));
    auto t = detail::make_expression_ptr(detail::make_node<VariableNode>("t"));
    auto f = detail::make_expression_ptr(detail::make_node<VariableNode>("f"));
    auto matrix = detail::make_expression_ptr(detail::make_node<VariableNode>("A"));
    auto period = detail::expression_from_node(
        detail::make_node<NumberNode>(BigInt(5)));
    auto lower = detail::make_expression_ptr(detail::make_node<NumberNode>(BigInt(1)));
    auto upper = detail::make_expression_ptr(detail::make_node<NumberNode>(BigInt(3)));
    Interval interval;
    interval.lower = Endpoint::closed(lower);
    interval.upper = Endpoint::closed(upper);
    auto &properties = ctx.current_properties();
    EXPECT_TRUE((properties.declare_domain("n", Domain::Natural).has_value())) << "natural domain declared";
    EXPECT_TRUE((properties.declare_sign("n", Sign::Positive).has_value())) << "positive declared";
    EXPECT_TRUE((properties.declare_transcendental("t").has_value())) << "transcendence declared";
    EXPECT_TRUE((properties.declare_finiteness("f", Finiteness::Finite).has_value())) << "finiteness declared";
    EXPECT_TRUE((properties.declare_definiteness("A", Definiteness::PositiveDefinite).has_value())) << "definiteness declared";
    EXPECT_TRUE((properties.declare_periodic("f", "x", period).has_value())) << "period declared";
    EXPECT_TRUE((properties.declare_bounded("n", Boundedness::Bounded, interval).has_value())) << "bounds declared";
    EXPECT_TRUE((properties.declare_monotonicity("f", "x", interval, Monotonicity::Increasing).has_value())) << "monotonicity declared";
    QueryInterface query(ctx);
    for (int level = 0; level < 3; ++level) {
        if (level != 0) {
            ctx.push();
        }
        expect_scalar_property_readthrough(ctx, query, n, t);
        EXPECT_TRUE((query.query_finite(*f).value() == Tribool::True)) << "finite survives empty push";
        EXPECT_TRUE((query.query_positive_definite(*matrix).value() == Tribool::True)) << "definiteness survives empty push";
        EXPECT_TRUE((query.query_positive_semidefinite(*matrix).value() == Tribool::True)) << "semidefiniteness survives empty push";
        EXPECT_TRUE((query.query_periodic(*f, "x").value() == Tribool::True)) << "periodicity survives empty push";
        auto inherited_period = query.get_period(*f, "x");
        EXPECT_TRUE((inherited_period && inherited_period.value() &&
                     detail::node(*inherited_period.value())->equals(*detail::node(period))))
            << "exact period survives empty push";
        auto bounds = ctx.get_bounds("n");
        EXPECT_TRUE((bounds && detail::node(*bounds->lower.value)->equals(*detail::node(*lower)) &&
                     detail::node(*bounds->upper.value)->equals(*detail::node(*upper))))
            << "bounds survive empty push";
        InferenceEngine engine(ctx);
        EXPECT_TRUE((engine.infer_monotonicity(*f, "x", interval) == Monotonicity::Increasing)) << "monotonicity reader inherits parent interval";
    }
    expect_extended_property_shadowing(ctx, query, matrix, f, period, interval);
    EXPECT_TRUE((ctx.pop().has_value())) << "empty scope popped";
    EXPECT_TRUE((query.query_integer(*n).value() == Tribool::True)) << "root facts remain intact";
}
