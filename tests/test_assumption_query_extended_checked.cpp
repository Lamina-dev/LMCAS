
#include "test_common.hpp"
#include "query_interface.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>

using namespace LMCAS;

static SymbolicExpr make_var(const std::string &name) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
    return expr;
}

static SymbolicExpr make_int(int v) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(v)));
    return expr;
}

static void test_checked_extended_query_contracts_algebraicity(QueryInterface &qi) {
    auto algebraic = qi.query_algebraic_checked(make_var("a"));
    ASSERT_TRUE((algebraic.has_value())) << "checked query_algebraic succeeds";
    if (algebraic) {
        EXPECT_TRUE((algebraic.value() == Tribool::True)) << "checked query_algebraic returns True for algebraic symbol";
    }

    auto transcendental = qi.query_transcendental_checked(make_var("tau"));
    ASSERT_TRUE((transcendental.has_value())) << "checked query_transcendental succeeds";
    if (transcendental) {
        EXPECT_TRUE((transcendental.value() == Tribool::True)) << "checked query_transcendental returns True for transcendental symbol";
    }
}

static void test_checked_extended_query_contracts_finiteness(QueryInterface &qi) {
    auto finite = qi.query_finite_checked(make_var("finite_symbol"));
    ASSERT_TRUE((finite.has_value())) << "checked query_finite succeeds";
    if (finite) {
        EXPECT_TRUE((finite.value() == Tribool::True)) << "checked query_finite returns True for finite symbol";
    }

    auto divergent = qi.query_divergent_checked(make_var("divergent_symbol"));
    ASSERT_TRUE((divergent.has_value())) << "checked query_divergent succeeds";
    if (divergent) {
        EXPECT_TRUE((divergent.value() == Tribool::True)) << "checked query_divergent returns True for divergent symbol";
    }
}

static void test_checked_extended_query_contracts_periodicity(QueryInterface &qi) {
    auto periodic = qi.query_periodic_checked(make_var("periodic_symbol"), "x");
    ASSERT_TRUE((periodic.has_value())) << "checked query_periodic succeeds";
    if (periodic) {
        EXPECT_TRUE((periodic.value() == Tribool::True)) << "checked query_periodic returns True for periodic symbol";
    }

    auto period = qi.get_period_checked(make_var("periodic_symbol"), "x");
    ASSERT_TRUE((period.has_value())) << "checked get_period succeeds";
    if (period) {
        EXPECT_TRUE((period.value().has_value())) << "checked get_period returns a declared period";
    }
}

static void test_checked_extended_query_contracts_definiteness(QueryInterface &qi) {
    auto positive_definite = qi.query_positive_definite_checked(make_var("M"));
    ASSERT_TRUE((positive_definite.has_value())) << "checked query_positive_definite succeeds";
    if (positive_definite) {
        EXPECT_TRUE((positive_definite.value() == Tribool::True)) << "checked query_positive_definite returns True for PD symbol";
    }

    auto positive_semidefinite = qi.query_positive_semidefinite_checked(make_var("M"));
    ASSERT_TRUE((positive_semidefinite.has_value())) << "checked query_positive_semidefinite succeeds";
    if (positive_semidefinite) {
        EXPECT_TRUE((positive_semidefinite.value() == Tribool::True)) << "checked query_positive_semidefinite returns True for PD symbol";
    }
}

static void test_checked_extended_query_contracts_conditions(QueryInterface &qi) {
    auto conditions = qi.query_conditions_checked(make_var("x"), Sign::Positive);
    ASSERT_TRUE((conditions.has_value())) << "checked query_conditions succeeds";
    if (conditions) {
        EXPECT_TRUE((conditions.value().size() == 1)) << "checked query_conditions preserves simple-variable conditions";
    }

    auto unsupported_conditions = qi.query_conditions_checked(make_int(1), Sign::Positive);
    ASSERT_TRUE((unsupported_conditions.has_value())) << "checked query_conditions accepts valid unsupported expressions";
    if (unsupported_conditions) {
        EXPECT_TRUE((unsupported_conditions.value().empty())) << "checked query_conditions returns empty set for unsupported expressions";
    }
}

TEST(AssumptionQueryExtendedChecked, CheckedExtendedQueryContracts) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("a", Domain::Algebraic).has_value());
    EXPECT_TRUE((ctx.current_properties().declare_transcendental("tau").has_value())) << "transcendental declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("finite_symbol", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("divergent_symbol", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_periodic("periodic_symbol", "x", make_int(6)).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_definiteness("M", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    QueryInterface qi(ctx);

    test_checked_extended_query_contracts_algebraicity(qi);
    test_checked_extended_query_contracts_finiteness(qi);
    test_checked_extended_query_contracts_periodicity(qi);
    test_checked_extended_query_contracts_definiteness(qi);
    test_checked_extended_query_contracts_conditions(qi);
}
