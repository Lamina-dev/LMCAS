
#include "test_common.hpp"
#include "query_interface.hpp"
#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "expr.hpp"
#include "residual_verification.hpp"
#include <memory>

using namespace LMCAS;

static SymbolicExpr make_var(const std::string &name) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
    return expr;
}

/// Build x - y as AddNode([x, MultiplyNode([-1, y])])
static SymbolicExpr make_subtraction(const std::string &lhs, const std::string &rhs) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(lhs);
    auto y_node = LMCAS::detail::make_node<VariableNode>(rhs);
    auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
    auto neg_y = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, y_node});
    auto add = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, neg_y});
    auto expr = LMCAS::detail::expression_from_node(add);
    return expr;
}

/// Build sin(x) as FunctionNode(Sin, [VariableNode(x)])
static SymbolicExpr make_sin(const std::string &var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto sin_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(sin_node);
    return expr;
}

/// Build cos(x) as FunctionNode(Cos, [VariableNode(x)])
static SymbolicExpr make_cos(const std::string &var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto cos_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Cos,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(cos_node);
    return expr;
}

/// Build tan(x) as FunctionNode(Tan, [VariableNode(x)])
static SymbolicExpr make_tan(const std::string &var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto tan_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Tan,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(tan_node);
    return expr;
}

TEST(AssumptionQueryExt, CacheHitReturnsSameResult) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // First query - computes and caches
    Tribool result1 = qi.query_positive(x_expr).value();
    EXPECT_EQ((result1), (Tribool::True)) << "First query: x is Positive";

    // Second query - should return cached result (same value)
    Tribool result2 = qi.query_positive(x_expr).value();
    EXPECT_EQ((result2), (Tribool::True)) << "Second query (cached): x is Positive";

    // Verify consistency across different property queries on same expression
    Tribool neg1 = qi.query_negative(x_expr).value();
    Tribool neg2 = qi.query_negative(x_expr).value();
    EXPECT_EQ((neg1), (Tribool::False)) << "First query: x is not Negative";
    EXPECT_EQ((neg2), (Tribool::False)) << "Second query (cached): x is not Negative";
}

TEST(AssumptionQueryExt, CacheInvalidationOnPush) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Populate cache
    Tribool before = qi.query_positive(x_expr).value();
    EXPECT_EQ((before), (Tribool::True)) << "Before push: x is Positive";

    // Manually invalidate cache (simulating push_scope hook)
    qi.invalidate_cache();

    // After invalidation, query should still return correct result (recomputed)
    Tribool after = qi.query_positive(x_expr).value();
    EXPECT_EQ((after), (Tribool::True)) << "After invalidate (push): x still Positive";
}

TEST(AssumptionQueryExt, CacheInvalidationOnPop) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // x is undeclared -> Unknown
    Tribool before_push = qi.query_positive(x_expr).value();
    EXPECT_EQ((before_push), (Tribool::Unknown)) << "Before push: x is Unknown";

    // Push scope and declare x Positive
    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    qi.invalidate_cache(); // Simulate hook

    Tribool in_scope = qi.query_positive(x_expr).value();
    EXPECT_EQ((in_scope), (Tribool::True)) << "In pushed scope: x is Positive";

    // Pop scope
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    qi.invalidate_cache(); // Simulate hook

    // After pop, x should be Unknown again
    Tribool after_pop = qi.query_positive(x_expr).value();
    EXPECT_EQ((after_pop), (Tribool::Unknown)) << "After pop: x is Unknown again";
}

TEST(AssumptionQueryExt, CacheInvalidationOnAssume) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Initially unknown
    Tribool before = qi.query_positive(x_expr).value();
    EXPECT_EQ((before), (Tribool::Unknown)) << "Before assume: x is Unknown";

    // Declare x Positive
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    qi.invalidate_cache(); // Simulate hook

    // Now should be True
    Tribool after = qi.query_positive(x_expr).value();
    EXPECT_EQ((after), (Tribool::True)) << "After assume + invalidate: x is Positive";
}

TEST(AssumptionQueryExt, QueryConditionsSimpleVariable) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Query: under what conditions is x > 0?
    auto conditions_result = qi.query_conditions(x_expr, Sign::Positive);
    ASSERT_TRUE((conditions_result.has_value())) << "Condition query succeeds";
    if (!conditions_result) {
        return;
    }
    const auto &conditions = conditions_result.value();

    EXPECT_TRUE((conditions.size() == 1)) << "Single variable should return exactly 1 condition set";

    if (!conditions.empty()) {
        const auto &cs = conditions[0];
        EXPECT_TRUE((cs.sign_conditions.size() == 1)) << "Condition set should have 1 sign condition";
        if (!cs.sign_conditions.empty()) {
            EXPECT_TRUE((cs.sign_conditions[0].first == "x")) << "Sign condition variable should be 'x'";
            EXPECT_TRUE((cs.sign_conditions[0].second == Sign::Positive)) << "Sign condition should be Positive";
        }
        EXPECT_TRUE((cs.domain_conditions.empty())) << "No domain conditions for simple sign query";
        EXPECT_TRUE((cs.relational_conditions.empty())) << "No relational conditions for simple variable";
    }

    // Also test with Negative target
    auto neg_conditions_result = qi.query_conditions(x_expr, Sign::Negative);
    ASSERT_TRUE((neg_conditions_result.has_value())) << "Negative condition query succeeds";
    if (!neg_conditions_result) {
        return;
    }
    const auto &neg_conditions = neg_conditions_result.value();
    EXPECT_TRUE((neg_conditions.size() == 1)) << "Negative target: single variable returns 1 condition set";
    if (!neg_conditions.empty()) {
        EXPECT_TRUE((neg_conditions[0].sign_conditions[0].second == Sign::Negative)) << "Negative target: condition should be Negative";
    }
}

TEST(AssumptionQueryExt, QueryConditionsSubtraction) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);

    // Build x - y
    auto expr = make_subtraction("x", "y");

    // Query: under what conditions is (x - y) > 0?
    auto conditions_result = qi.query_conditions(expr, Sign::Positive);
    ASSERT_TRUE((conditions_result.has_value())) << "Composite condition query succeeds";
    if (!conditions_result) {
        return;
    }
    const auto &conditions = conditions_result.value();

    // Should return at least 2 condition sets:
    // 1. {x: Positive, y: Negative}
    // 2. {x GT y, y: NonNegative}
    EXPECT_TRUE((conditions.size() >= 2)) << "x - y Positive should return at least 2 condition sets";

    if (conditions.size() < 2) {
        return;
    }
    // First condition set: x Positive AND y Negative
    const auto &cs1 = conditions[0];
    EXPECT_TRUE((cs1.sign_conditions.size() == 2)) << "First condition set has 2 sign conditions";

    bool has_x_positive = false;
    bool has_y_negative = false;
    for (const auto &sc : cs1.sign_conditions) {
        if (sc.first == "x" && sc.second == Sign::Positive) {
            has_x_positive = true;
        }
        if (sc.first == "y" && sc.second == Sign::Negative) {
            has_y_negative = true;
        }
    }
    EXPECT_TRUE((has_x_positive)) << "CS1: x should be Positive";
    EXPECT_TRUE((has_y_negative)) << "CS1: y should be Negative";
}

TEST(AssumptionQueryExt, QueryPositiveDefinite) {
    AssumptionContext ctx;
    // Declare matrix symbol M as PositiveDefinite
    EXPECT_TRUE((ctx.current_properties().declare_definiteness("M", Definiteness::PositiveDefinite).has_value())) << "definiteness declaration succeeds";

    QueryInterface qi(ctx);
    auto m_expr = make_var("M");

    EXPECT_EQ((qi.query_positive_definite(m_expr).value()), (Tribool::True)) << "M declared PositiveDefinite: query_positive_definite = True";
    EXPECT_EQ((qi.query_positive_semidefinite(m_expr).value()), (Tribool::True)) << "M declared PositiveDefinite: query_positive_semidefinite = True (implied)";
}

TEST(AssumptionQueryExt, QueryPositiveSemidefinite) {
    AssumptionContext ctx;
    // Declare matrix symbol A as PositiveSemiDefinite only
    EXPECT_TRUE((ctx.current_properties().declare_definiteness("A", Definiteness::PositiveSemiDefinite).has_value())) << "definiteness declaration succeeds";

    QueryInterface qi(ctx);
    auto a_expr = make_var("A");

    EXPECT_EQ((qi.query_positive_semidefinite(a_expr).value()), (Tribool::True)) << "A declared PSD: query_positive_semidefinite = True";
    /// PositiveSemiDefinite 对 PositiveDefinite 查询保持 Unknown.
    EXPECT_EQ((qi.query_positive_definite(a_expr).value()), (Tribool::Unknown)) << "A declared PSD: query_positive_definite = Unknown";
}

TEST(AssumptionQueryExt, QueryDefinitenessNegative) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.current_properties().declare_definiteness("N", Definiteness::NegativeDefinite).has_value())) << "definiteness declaration succeeds";

    QueryInterface qi(ctx);
    auto n_expr = make_var("N");

    EXPECT_EQ((qi.query_positive_definite(n_expr).value()), (Tribool::False)) << "N declared NegDef: query_positive_definite = False";
    EXPECT_EQ((qi.query_positive_semidefinite(n_expr).value()), (Tribool::False)) << "N declared NegDef: query_positive_semidefinite = False";
}

TEST(AssumptionQueryExt, QueryDefinitenessUndeclared) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto u_expr = make_var("U");

    EXPECT_EQ((qi.query_positive_definite(u_expr).value()), (Tribool::Unknown)) << "Undeclared: query_positive_definite = Unknown";
    EXPECT_EQ((qi.query_positive_semidefinite(u_expr).value()), (Tribool::Unknown)) << "Undeclared: query_positive_semidefinite = Unknown";
}

TEST(AssumptionQueryExt, QueryAlgebraic) {
    AssumptionContext ctx;
    // Declare x as Algebraic domain
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Algebraic).has_value());

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    EXPECT_EQ((qi.query_algebraic(x_expr).value()), (Tribool::True)) << "x declared Algebraic: query_algebraic = True";
    EXPECT_EQ((qi.query_transcendental(x_expr).value()), (Tribool::False)) << "x declared Algebraic: query_transcendental = False";

    // Undeclared variable
    auto y_expr = make_var("y");
    EXPECT_EQ((qi.query_algebraic(y_expr).value()), (Tribool::Unknown)) << "y undeclared: query_algebraic = Unknown";
}

TEST(AssumptionQueryExt, QueryTranscendental) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.current_properties().declare_transcendental("pi_sym").has_value())) << "transcendental declaration succeeds";

    QueryInterface qi(ctx);
    auto pi_expr = make_var("pi_sym");

    EXPECT_EQ((qi.query_transcendental(pi_expr).value()), (Tribool::True)) << "pi_sym declared Transcendental: query_transcendental = True";
    EXPECT_EQ((qi.query_algebraic(pi_expr).value()), (Tribool::False)) << "pi_sym declared Transcendental: query_algebraic = False";
}

TEST(AssumptionQueryExt, QueryFinite) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("a", Finiteness::Finite).has_value())) << "finiteness declaration succeeds";

    QueryInterface qi(ctx);
    auto a_expr = make_var("a");

    EXPECT_EQ((qi.query_finite(a_expr).value()), (Tribool::True)) << "a declared Finite: query_finite = True";
    EXPECT_EQ((qi.query_divergent(a_expr).value()), (Tribool::False)) << "a declared Finite: query_divergent = False";

    // Undeclared
    auto b_expr = make_var("b");
    EXPECT_EQ((qi.query_finite(b_expr).value()), (Tribool::Unknown)) << "b undeclared: query_finite = Unknown";
    EXPECT_EQ((qi.query_divergent(b_expr).value()), (Tribool::Unknown)) << "b undeclared: query_divergent = Unknown";
}

TEST(AssumptionQueryExt, QueryDivergent) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.current_properties().declare_finiteness("d", Finiteness::Divergent).has_value())) << "finiteness declaration succeeds";

    QueryInterface qi(ctx);
    auto d_expr = make_var("d");

    EXPECT_EQ((qi.query_divergent(d_expr).value()), (Tribool::True)) << "d declared Divergent: query_divergent = True";
    EXPECT_EQ((qi.query_finite(d_expr).value()), (Tribool::False)) << "d declared Divergent: query_finite = False";
}

TEST(AssumptionQueryExt, QueryPeriodicDeclared) {
    AssumptionContext ctx;
    auto period_expr = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(BigInt(6)));
    EXPECT_TRUE((ctx.current_properties().declare_periodic("f", "x", period_expr).has_value())) << "period declaration succeeds";

    QueryInterface qi(ctx);
    auto f_expr = make_var("f");

    EXPECT_EQ((qi.query_periodic(f_expr, "x").value()), (Tribool::True)) << "f declared periodic: query_periodic = True";
}

TEST(AssumptionQueryExt, QueryPeriodicTrig) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto sin_expr = make_sin("x");
    auto cos_expr = make_cos("x");
    auto tan_expr = make_tan("x");

    EXPECT_EQ((qi.query_periodic(sin_expr, "x").value()), (Tribool::True)) << "sin(x): query_periodic = True";
    EXPECT_EQ((qi.query_periodic(cos_expr, "x").value()), (Tribool::True)) << "cos(x): query_periodic = True";
    EXPECT_EQ((qi.query_periodic(tan_expr, "x").value()), (Tribool::True)) << "tan(x): query_periodic = True";
}

TEST(AssumptionQueryExt, GetPeriodTrig) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);
    const auto check = [&](const SymbolicExpr &expression, const char *reference) {
        auto period = qi.get_period(expression, "x");
        auto expected = parse_expr(reference);
        EXPECT_TRUE((period && period.value() && expected)) << "period query and exact reference succeed";
        if (!period || !period.value() || !expected)
            return;
        EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(*period.value()),
                                          expected.value())) << "period is mathematically equal to the exact pi multiple";
    };
    check(make_sin("x"), "2*pi");
    check(make_cos("x"), "2*pi");
    check(make_tan("x"), "pi");
}

TEST(AssumptionQueryExt, GetPeriodNonPeriodic) {
    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto x_expr = make_var("x");
    auto period = qi.get_period(x_expr, "x");
    EXPECT_TRUE((period.has_value())) << "Non-periodic query succeeds";
    EXPECT_TRUE((period.has_value() && !period.value().has_value())) << "Variable x (undeclared periodic) should have no period";
}

TEST(AssumptionQueryExt, CacheMultipleProperties) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Query multiple properties - each should be cached independently
    Tribool pos = qi.query_positive(x_expr).value();
    Tribool intg = qi.query_integer(x_expr).value();
    Tribool neg = qi.query_negative(x_expr).value();

    EXPECT_EQ((pos), (Tribool::True)) << "x Positive cached correctly";
    EXPECT_EQ((intg), (Tribool::True)) << "x Integer cached correctly";
    EXPECT_EQ((neg), (Tribool::False)) << "x not Negative cached correctly";

    // Query again - should hit cache
    EXPECT_EQ((qi.query_positive(x_expr).value()), (Tribool::True)) << "x Positive (cache hit)";
    EXPECT_EQ((qi.query_integer(x_expr).value()), (Tribool::True)) << "x Integer (cache hit)";
    EXPECT_EQ((qi.query_negative(x_expr).value()), (Tribool::False)) << "x not Negative (cache hit)";

    // Invalidate and re-query
    qi.invalidate_cache();
    EXPECT_EQ((qi.query_positive(x_expr).value()), (Tribool::True)) << "x Positive after invalidate";
    EXPECT_EQ((qi.query_integer(x_expr).value()), (Tribool::True)) << "x Integer after invalidate";
}
