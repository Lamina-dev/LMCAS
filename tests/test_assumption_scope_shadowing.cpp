
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "interval.hpp"
#include "query_interface.hpp"
#include "computation_context.hpp"
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <string>

using namespace LMCAS;

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

TEST(AssumptionScopeShadowing, DomainShadowing) {
    for (const auto &sym : TEST_SYMBOLS) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(sym, Domain::Real).has_value());

        Domain parent_domain = ctx.get_domain(sym);
        EXPECT_TRUE((parent_domain == Domain::Real)) << sym + " is Real in parent";

        ctx.push();
        EXPECT_TRUE(ctx.assume_domain(sym, Domain::Integer).has_value());

        EXPECT_TRUE((ctx.get_domain(sym) == Domain::Integer)) << sym + " is Integer in child (shadows Real)";

        EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

        EXPECT_TRUE((ctx.get_domain(sym) == Domain::Real)) << sym + " is Real again after pop";
    }
}

TEST(AssumptionScopeShadowing, SignShadowing) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    EXPECT_TRUE((ctx.has_sign("x", Sign::NonNegative))) << "x is NonNegative in parent";
    EXPECT_FALSE((ctx.has_sign("x", Sign::Positive))) << "x is NOT Positive in parent (only NonNegative)";

    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    EXPECT_TRUE((ctx.has_sign("x", Sign::Positive))) << "x is Positive in child (shadows parent)";
    EXPECT_TRUE((ctx.has_sign("x", Sign::NonNegative))) << "x is NonNegative in child (implied by Positive)";
    EXPECT_TRUE((ctx.has_sign("x", Sign::NonZero))) << "x is NonZero in child (implied by Positive)";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_TRUE((ctx.has_sign("x", Sign::NonNegative))) << "x is NonNegative after pop (parent value)";
    EXPECT_FALSE((ctx.has_sign("x", Sign::Positive))) << "x is NOT Positive after pop";
    EXPECT_FALSE((ctx.has_sign("x", Sign::NonZero))) << "x is NOT NonZero after pop (was only NonNegative in parent)";
}

TEST(AssumptionScopeShadowing, SignShadowsInheritedDomainPositivity) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value()) << "parent positive integer domain accepted";
    const auto x = SymbolicExpr::variable("x");
    QueryInterface query(ctx);
    auto positive = query.query_positive(*x);
    EXPECT_TRUE((positive && positive.value() == Tribool::True)) << "parent domain proves positivity";

    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value()) << "child negative sign accepted";
    positive = query.query_positive(*x);
    auto negative = query.query_negative(*x);
    EXPECT_TRUE((positive && positive.value() == Tribool::False)) << "child explicit negative sign disproves positivity";
    EXPECT_TRUE((negative && negative.value() == Tribool::True)) << "child explicit negative sign proves negativity";

    EXPECT_TRUE((ctx.pop().has_value())) << "child scope popped";
    positive = query.query_positive(*x);
    EXPECT_TRUE((positive && positive.value() == Tribool::True)) << "popping child restores parent domain positivity";
}

TEST(AssumptionScopeShadowing, ParityShadowing) {
    AssumptionContext ctx;

    EXPECT_TRUE((ctx.current_properties().declare_parity("x", Parity::Even).has_value())) << "parity declaration succeeds";
    EXPECT_TRUE((ctx.get_parity("x") == Parity::Even)) << "x is Even in parent";

    ctx.push();
    EXPECT_TRUE((ctx.current_properties().declare_parity("x", Parity::Odd).has_value())) << "parity declaration succeeds";

    EXPECT_TRUE((ctx.get_parity("x") == Parity::Odd)) << "x is Odd in child (shadows Even)";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_TRUE((ctx.get_parity("x") == Parity::Even)) << "x is Even after pop (parent value restored)";
}

TEST(AssumptionScopeShadowing, BoundednessShadowing) {
    AssumptionContext ctx;

    EXPECT_TRUE((ctx.current_properties().declare_bounded("x", Boundedness::Unbounded).has_value())) << "boundedness declaration succeeds";
    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Unbounded)) << "x is Unbounded in parent";

    ctx.push();
    EXPECT_TRUE((ctx.current_properties().declare_bounded("x", Boundedness::Bounded).has_value())) << "boundedness declaration succeeds";

    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Bounded)) << "x is Bounded in child (shadows Unbounded)";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_TRUE((ctx.get_boundedness("x") == Boundedness::Unbounded)) << "x is Unbounded after pop (parent value restored)";
}

TEST(AssumptionScopeShadowing, DomainShadowingAllPairs) {
    struct DomainPair {
        Domain parent;
        Domain child;
    };
    std::vector<DomainPair> pairs = {
        {Domain::Real, Domain::Integer},
        {Domain::Real, Domain::Natural},
        {Domain::Rational, Domain::PositiveInt},
        {Domain::Integer, Domain::PositiveInt},
        {Domain::Real, Domain::PositiveInt},
    };

    for (const auto &p : pairs) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain("x", p.parent).has_value());

        ctx.push();
        EXPECT_TRUE(ctx.assume_domain("x", p.child).has_value());

        EXPECT_TRUE((ctx.get_domain("x") == p.child)) << "x is " + domain_name(p.child) + " in child (parent was " +
                                                             domain_name(p.parent) + ")";

        EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

        EXPECT_TRUE((ctx.get_domain("x") == p.parent)) << "x is " + domain_name(p.parent) + " after pop";
    }
}

TEST(AssumptionScopeShadowing, SignShadowingVarious) {
    struct SignPair {
        Sign parent;
        Sign child;
    };
    std::vector<SignPair> pairs = {
        {Sign::NonNegative, Sign::Positive},
        {Sign::NonPositive, Sign::Negative},
        {Sign::NonZero, Sign::Positive},
        {Sign::NonZero, Sign::Negative},
        {Sign::NonNegative, Sign::Zero},
    };

    for (const auto &p : pairs) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("x", p.parent).has_value());

        ctx.push();
        EXPECT_TRUE(ctx.assume_sign("x", p.child).has_value());

        EXPECT_TRUE((ctx.has_sign("x", p.child))) << "x has " + sign_name(p.child) + " in child (parent was " +
                                                         sign_name(p.parent) + ")";

        EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

        EXPECT_TRUE((ctx.has_sign("x", p.parent))) << "x has " + sign_name(p.parent) + " after pop";
    }
}

TEST(AssumptionScopeShadowing, ChildDoesNotModifyParent) {
    AssumptionContext ctx;

    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    ctx.push();

    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value());

    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer in child";
    EXPECT_TRUE((ctx.has_sign("y", Sign::Negative))) << "y is Negative in child";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_TRUE((ctx.get_domain("x") == Domain::Real)) << "x is still Real in parent (not modified by child)";
    EXPECT_TRUE((ctx.has_sign("y", Sign::Positive))) << "y is still Positive in parent (not modified by child)";
    EXPECT_FALSE((ctx.has_sign("y", Sign::Negative))) << "y is NOT Negative in parent";
}

TEST(AssumptionScopeShadowing, ReadThroughUndeclaredInChild) {
    AssumptionContext ctx;

    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    ctx.push();

    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x reads through to Integer from parent";
    EXPECT_TRUE((ctx.has_sign("y", Sign::Positive))) << "y reads through to Positive from parent";

    EXPECT_TRUE(ctx.assume_sign("z", Sign::NonNegative).has_value());
    EXPECT_TRUE((ctx.has_sign("z", Sign::NonNegative))) << "z is NonNegative in child";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_FALSE((ctx.has_sign("z", Sign::NonNegative))) << "z not visible after pop";
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x still Integer in root";
    EXPECT_TRUE((ctx.has_sign("y", Sign::Positive))) << "y still Positive in root";
}

TEST(AssumptionScopeShadowing, MultiLevelShadowing) {
    AssumptionContext ctx;

    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer at level 2";

    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());
    EXPECT_TRUE((ctx.get_domain("x") == Domain::PositiveInt)) << "x is PositiveInt at level 3";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Integer)) << "x is Integer at level 2 after popping level 3";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.get_domain("x") == Domain::Real)) << "x is Real at root after popping level 2";
}

TEST(AssumptionScopeShadowing, DifferentSymbolsIndependent) {
    AssumptionContext ctx;

    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    EXPECT_TRUE(ctx.assume_domain("y", Domain::Integer).has_value());

    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Natural).has_value());

    EXPECT_TRUE((ctx.get_domain("x") == Domain::Natural)) << "x is Natural in child (shadowed)";
    EXPECT_TRUE((ctx.get_domain("y") == Domain::Integer)) << "y is Integer in child (read-through, not shadowed)";

    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    EXPECT_TRUE((ctx.get_domain("x") == Domain::Real)) << "x is Real after pop";
    EXPECT_TRUE((ctx.get_domain("y") == Domain::Integer)) << "y is Integer after pop (unchanged)";
}

TEST(AssumptionScopeShadowing, RelationVisibilityAndShadowing) {
    AssumptionContext ctx;
    const auto x = detail::make_node<VariableNode>("x");
    const auto y = detail::make_node<VariableNode>("y");
    const auto sum = detail::expression_from_node(detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, y}));
    const auto one = detail::expression_from_node(detail::make_node<NumberNode>(BigInt(1)));
    EXPECT_TRUE((ctx.current_relations().add_relation(sum, one, RelationOp::GT,
                                                      ctx.current_properties())
                     .has_value()))
        << "parent relation declared";
    QueryInterface query(ctx);
    EXPECT_TRUE((query.query_positive(sum).value() == Tribool::True)) << "parent relation proves positive";
    ctx.push();
    EXPECT_TRUE((ctx.has_relation(sum, one, RelationOp::GT))) << "empty child inherits relation";
    EXPECT_TRUE((query.query_positive(sum).value() == Tribool::True)) << "empty child keeps cached proof";
    EXPECT_TRUE((ctx.current_properties().declare_sign("x", Sign::Negative).has_value())) << "child sign declared";
    EXPECT_FALSE((ctx.has_relation(sum, one, RelationOp::GT))) << "changed dependency shadows parent relation";
    EXPECT_TRUE((query.query_positive(sum).value() == Tribool::Unknown)) << "shadowed relation cannot reprove positivity";
    EXPECT_TRUE((ctx.pop().has_value())) << "child popped";
    EXPECT_TRUE((query.query_positive(sum).value() == Tribool::True)) << "parent relation restored";
    ctx.push();
    EXPECT_TRUE((ctx.current_relations().add_relation(sum, one, RelationOp::LT,
                                                      ctx.current_properties())
                     .has_value()))
        << "child relation overrides same lhs";
    EXPECT_FALSE((ctx.has_relation(sum, one, RelationOp::GT))) << "parent same-lhs relation shadowed";
    EXPECT_TRUE((ctx.has_relation(sum, one, RelationOp::LT))) << "child relation visible";
    EXPECT_TRUE((query.query_positive(sum).value() == Tribool::Unknown)) << "old positive proof not reused";
    const auto condition = detail::expression_from_node(detail::make_node<RelationalNode>(
        detail::node(sum), detail::node(one), RelationOp::GT));
    EXPECT_TRUE((ctx.evaluate_condition(condition) == Tribool::Unknown)) << "condition evaluation respects relation visibility";
    EXPECT_TRUE((ctx.pop().has_value())) << "relation child popped";
    EXPECT_TRUE((ctx.evaluate_condition(condition) == Tribool::True)) << "parent condition restored";
}

TEST(AssumptionScopeShadowing, TransitiveRelationShadowing) {
    AssumptionContext ctx;
    auto sum = SymbolicExpr::add(SymbolicExpr::variable("x"), SymbolicExpr::variable("z"));
    auto intermediate = SymbolicExpr::variable("y");
    auto zero = SymbolicExpr::number(0);
    EXPECT_TRUE((ctx.current_relations().add_relation(*sum, *intermediate, RelationOp::GT,
                                                      ctx.current_properties())
                     .has_value()))
        << "first transitive premise declared";
    EXPECT_TRUE((ctx.current_relations().add_relation(*intermediate, *zero, RelationOp::GT,
                                                      ctx.current_properties())
                     .has_value()))
        << "second transitive premise declared";
    QueryInterface query(ctx);
    auto parent = query.query_positive(*sum);
    EXPECT_TRUE((parent && parent.value() == Tribool::True)) << "parent deduction proves positive";
    ctx.push();
    auto inherited = query.query_positive(*sum);
    EXPECT_TRUE((inherited && inherited.value() == Tribool::True)) << "empty scope inherits deduction";
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Negative).has_value()) << "intermediate variable shadowed";
    EXPECT_FALSE((ctx.has_relation(*sum, *zero, RelationOp::GT))) << "deduction cannot outlive its hidden premises";
    auto shadowed = query.query_positive(*sum);
    EXPECT_TRUE((shadowed && shadowed.value() == Tribool::Unknown)) << "same cached query cannot reuse a shadowed transitive proof";
    EXPECT_TRUE((ctx.pop().has_value())) << "child scope popped";
    auto restored = query.query_positive(*sum);
    EXPECT_TRUE((restored && restored.value() == Tribool::True)) << "parent transitive proof restored";
}

TEST(AssumptionScopeShadowing, ConstantLeftHandSideUsesVariableSubject) {
    AssumptionContext ctx;
    const auto one = SymbolicExpr::number(1);
    const auto two = SymbolicExpr::number(2);
    const auto x = SymbolicExpr::variable("x");
    const auto y = SymbolicExpr::variable("y");
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *one, *x, RelationOp::LT, ctx.current_properties()));
    ctx.push();
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *one, *y, RelationOp::LT, ctx.current_properties()));
    EXPECT_TRUE(ctx.has_relation(*one, *x, RelationOp::LT));
    EXPECT_TRUE(ctx.has_relation(*one, *y, RelationOp::LT));
    ASSERT_TRUE(ctx.pop());

    ctx.push();
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *two, *x, RelationOp::LT, ctx.current_properties()));
    EXPECT_FALSE(ctx.has_relation(*one, *x, RelationOp::LT));
    EXPECT_TRUE(ctx.has_relation(*two, *x, RelationOp::LT));
    ASSERT_TRUE(ctx.pop());
    EXPECT_TRUE(ctx.has_relation(*one, *x, RelationOp::LT));
}

TEST(AssumptionScopeShadowing, ReverseSpellingRebindsEitherEndpoint) {
    for (bool reversed : {false, true}) {
        AssumptionContext ctx;
        const auto one = SymbolicExpr::number(1);
        const auto two = SymbolicExpr::number(2);
        const auto x = SymbolicExpr::variable("x");
        const auto y = SymbolicExpr::variable("y");
        const auto sum = SymbolicExpr::add(x, y);
        const auto* lhs = x.get();
        const auto* rhs = one.get();
        const auto* child_lhs = two.get();
        const auto* child_rhs = x.get();
        auto op = RelationOp::GT;
        auto child_op = RelationOp::LT;
        if (reversed) {
            lhs = one.get();
            rhs = x.get();
            child_lhs = x.get();
            child_rhs = two.get();
            op = RelationOp::LT;
            child_op = RelationOp::GT;
        }
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *lhs, *rhs, op, ctx.current_properties()));
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *one, *sum, RelationOp::LT, ctx.current_properties()));
        ctx.push();
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *child_lhs, *child_rhs, child_op, ctx.current_properties()));
        EXPECT_FALSE(ctx.has_relation(*lhs, *rhs, op));
        EXPECT_FALSE(ctx.has_relation(*one, *sum, RelationOp::LT));
        ASSERT_TRUE(ctx.pop());
        EXPECT_TRUE(ctx.has_relation(*lhs, *rhs, op));
        EXPECT_TRUE(ctx.has_relation(*one, *sum, RelationOp::LT));
    }
}

TEST(AssumptionScopeShadowing, CompositeSubjectDoesNotRebindItsParameters) {
    AssumptionContext ctx;
    const auto a = SymbolicExpr::variable("a");
    const auto b = SymbolicExpr::variable("b");
    const auto c = SymbolicExpr::variable("c");
    const auto sum = SymbolicExpr::add(a, b);
    const auto equal_sum = SymbolicExpr::add(
        SymbolicExpr::variable("a"), SymbolicExpr::variable("b"));
    const auto one = SymbolicExpr::number(1);
    const auto two = SymbolicExpr::number(2);
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *b, *c, RelationOp::GT, ctx.current_properties()));
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *one, *sum, RelationOp::LT, ctx.current_properties()));
    ctx.push();
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *equal_sum, *two, RelationOp::GT, ctx.current_properties()));
    EXPECT_FALSE(ctx.has_relation(*one, *sum, RelationOp::LT));
    EXPECT_TRUE(ctx.has_relation(*b, *c, RelationOp::GT));
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *a, *b, RelationOp::GT, ctx.current_properties()));
    EXPECT_TRUE(ctx.has_relation(*b, *c, RelationOp::GT));
    ASSERT_TRUE(ctx.pop());
    EXPECT_TRUE(ctx.has_relation(*one, *sum, RelationOp::LT));
    ctx.push();
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *one, *a, RelationOp::GT, ctx.current_properties()));
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *a, *b, RelationOp::GT, ctx.current_properties()));
    ASSERT_TRUE(ctx.has_relation(*one, *b, RelationOp::GT));
    EXPECT_TRUE(ctx.has_relation(*b, *c, RelationOp::GT));
}
