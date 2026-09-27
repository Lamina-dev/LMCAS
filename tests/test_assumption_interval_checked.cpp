
#include "test_interval_support.hpp"
#include "assumption.hpp"
#include "property_store.hpp"
#include "interval.hpp"
#include "internal/symbolic_ast.hpp"
#include "assumption_context.hpp"
#include "expr.hpp"
#include <string>
#include <stdexcept>
#include <vector>

using namespace LMCAS;

static void test_checked_interval_property_contracts_continuity(PropertyStore &store, const Interval &interval) {
    auto empty_symbol = store.declare_continuous_checked("", interval);
    EXPECT_TRUE((!empty_symbol && empty_symbol.error().code == CasErrc::InvalidArgument)) << "checked continuity rejects an empty symbol";
    EXPECT_TRUE((store.get_all_symbols().empty())) << "failed checked continuity creates no property record";

    auto empty_interval = store.declare_continuous_checked("f", Interval::empty());
    EXPECT_TRUE((!empty_interval && empty_interval.error().code == CasErrc::InvalidArgument)) << "checked continuity rejects an empty interval";
    EXPECT_TRUE((store.get_all_symbols().empty())) << "empty-interval failure is transactional";

    auto differentiable = store.declare_differentiable_checked("f", interval);
    EXPECT_TRUE((differentiable.has_value())) << "checked differentiability accepts a valid interval";
    const auto declaration_count = store.get_continuity_decls("f").size();
    auto duplicate = store.declare_differentiable_checked("f", interval);
    EXPECT_TRUE((duplicate.has_value())) << "exact duplicate differentiability is idempotent";
    EXPECT_TRUE((store.get_continuity_decls("f").size() == declaration_count)) << "idempotent checked declaration does not duplicate state";

    auto downgrade = store.declare_continuous_checked("f", interval);
    EXPECT_TRUE((!downgrade && downgrade.error().code == CasErrc::InvalidArgument)) << "checked continuity rejects a differentiability downgrade";
    EXPECT_TRUE((store.get_continuity_decls("f").size() == declaration_count)) << "failed checked downgrade preserves declaration state";
}

static void test_checked_interval_property_contracts_boundaries() {
    PropertyStore boundary_store;
    Interval left = make_closed_interval(0.0, 5.0);
    Interval open_touch{
        Endpoint::open(SymbolicExpr::number(5.0)),
        Endpoint::closed(SymbolicExpr::number(10.0))};
    Interval closed_touch = make_closed_interval(5.0, 10.0);
    EXPECT_TRUE((boundary_store.declare_differentiable_checked("edge", left).has_value())) << "boundary test stores differentiability";
    auto disjoint_touch = boundary_store.declare_continuous_checked("edge", open_touch);
    EXPECT_TRUE((disjoint_touch.has_value())) << "open touching boundary is correctly treated as non-overlapping";
    auto overlapping_touch = boundary_store.declare_continuous_checked("edge", closed_touch);
    EXPECT_TRUE((!overlapping_touch &&
                 overlapping_touch.error().code == CasErrc::InvalidArgument))
        << "closed touching boundary is correctly treated as overlapping";
}

static void test_checked_interval_property_contracts_errors_and_budgets(PropertyStore &store, const Interval &interval) {
    Interval symbolic{
        Endpoint::closed(SymbolicExpr::variable("a")),
        Endpoint::closed(SymbolicExpr::variable("b"))};
    auto symbolic_declaration = store.declare_continuous_checked("g", symbolic);
    EXPECT_TRUE((!symbolic_declaration &&
                 symbolic_declaration.error().code == CasErrc::UnboundSymbol))
        << "checked continuity propagates unbound symbolic endpoints";
    auto symbolic_query = store.is_continuous_checked("f", symbolic);
    EXPECT_TRUE((!symbolic_query && symbolic_query.error().code == CasErrc::UnboundSymbol)) << "checked continuity query propagates endpoint errors";

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = store.declare_monotonicity_checked(
        "g", "x", interval, Monotonicity::Increasing, cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << "checked monotonicity observes cancellation";
    EXPECT_TRUE((store.get_monotonicity_decls("g").empty())) << "cancelled monotonicity declaration stores no state";

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited = store.declare_differentiable_checked(
        "h", interval, limited_context);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "checked differentiability observes the step budget";
}

static void test_checked_interval_property_contracts_monotonicity(PropertyStore &store, const Interval &interval) {
    auto empty_variable = store.declare_monotonicity_checked(
        "f", "", interval, Monotonicity::Increasing);
    EXPECT_TRUE((!empty_variable && empty_variable.error().code == CasErrc::InvalidArgument)) << "checked monotonicity rejects an empty variable";
    auto monotonic = store.declare_monotonicity_checked(
        "f", "x", interval, Monotonicity::Increasing);
    EXPECT_TRUE((monotonic.has_value())) << "checked monotonicity stores a valid declaration";
    auto queried = store.get_monotonicity_checked("f", "x", interval);
    EXPECT_TRUE((queried && queried.value() == Monotonicity::Increasing)) << "checked monotonicity query retrieves the proven declaration";
}

TEST(AssumptionIntervalChecked, CheckedIntervalPropertyContracts) {
    PropertyStore store;
    Interval interval = make_closed_interval(0.0, 10.0);

    test_checked_interval_property_contracts_continuity(store, interval);
    test_checked_interval_property_contracts_boundaries();
    test_checked_interval_property_contracts_errors_and_budgets(store, interval);
    test_checked_interval_property_contracts_monotonicity(store, interval);
}

static void expect_monotonicity_roundtrip(const Interval &interval) {
    AssumptionContext ctx;
    auto declared = ctx.current_properties().declare_monotonicity_checked(
        "f", "x", interval, Monotonicity::Increasing);
    ASSERT_TRUE((declared.has_value())) << "monotonicity interval accepted";
    if (!declared) {
        return;
    }
    auto restored = AssumptionContext::deserialize_checked(ctx.serialize());
    ASSERT_TRUE((restored.has_value())) << "serialized approximate interval reparses";
    if (!restored) {
        return;
    }
    auto monotonicity = restored.value().get_monotonicity_checked("f", "x", interval);
    EXPECT_TRUE((monotonicity && monotonicity.value() == Monotonicity::Increasing)) << "monotonicity still covers the declared interval";
    for (const auto *endpoint : {&interval.lower, &interval.upper}) {
        if (!endpoint->value) {
            continue;
        }
        auto boundary = restored.value().get_monotonicity_checked(
            "f", "x", Interval::point(endpoint->value));
        EXPECT_TRUE((boundary && boundary.value() == (endpoint->is_open
                                                          ? Monotonicity::Unknown
                                                          : Monotonicity::Increasing)))
            << "roundtrip preserves open and closed endpoint coverage";
    }
}

TEST(AssumptionIntervalChecked, MonotonicityApproximateIntervalRoundtrip) {
    auto positive = approx_real(0.5);
    auto negative = approx_real(-0.5);
    auto upper = approx_real(2.0);
    EXPECT_TRUE((positive && negative && upper)) << "approximate endpoints constructed";
    if (!positive || !negative || !upper) {
        return;
    }
    for (const auto &interval : std::vector<Interval>{
             {Endpoint::closed(positive.value()), Endpoint::closed(upper.value())},
             {Endpoint::open(negative.value()), Endpoint::open(upper.value())},
             {Endpoint::closed(negative.value()), Endpoint::open(upper.value())},
             {Endpoint::open(negative.value()), Endpoint::closed(upper.value())},
             {Endpoint::neg_inf(), Endpoint::closed(upper.value())},
             {Endpoint::closed(negative.value()), Endpoint::pos_inf()}}) {
        expect_monotonicity_roundtrip(interval);
    }
    for (const char *interval : {"[approx(0.5, approx(2)]", "[approx(0.5), +inf]",
                                 "[-inf, approx(2)]"}) {
        auto restored = AssumptionContext::deserialize_checked(
            std::string("SCOPE 0\nMONOTONICITY f x ") + interval + " Increasing\nEND\n");
        EXPECT_TRUE((!restored && restored.error().code == CasErrc::ParseError)) << "malformed nested endpoints and closed infinities are rejected";
    }
}
