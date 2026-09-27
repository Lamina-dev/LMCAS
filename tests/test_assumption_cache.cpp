
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "query_interface.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_var_node(const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

static SymbolicExpr make_var_expr(const std::string &name) {
    auto expr = LMCAS::detail::expression_from_node(make_var_node(name));
    return expr;
}

/// Generate a random Sign from a subset of useful signs
static Sign random_sign() {
    std::vector<Sign> signs = {Sign::Positive, Sign::Negative, Sign::NonNegative, Sign::NonPositive};
    return *rc::gen::elementOf(signs);
}

// --- Test: invalidate_cache() clears the cache ---

TEST(LmcasAssumptionCache, InvalidateCacheClearsCache) {
    EXPECT_TRUE(rc::check("After invalidate_cache(), queries recompute and return correct results", []() {
        std::string var_name = "x_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // First query — populates cache
        Tribool result1 = qi.query_positive(expr).value();
        RC_ASSERT(result1 == Tribool::True);

        // Invalidate cache
        qi.invalidate_cache();

        // Second query — should recompute and still return True
        Tribool result2 = qi.query_positive(expr).value();
        RC_ASSERT(result2 == Tribool::True);

        // Results must be consistent
        RC_ASSERT(result1 == result2);
    }));
}

// --- Test: Cache stores results (same query returns same result) ---

TEST(LmcasAssumptionCache, CacheStoresResults) {
    EXPECT_TRUE(rc::check("Repeated queries on the same expression return the same cached result", []() {
        std::string var_name = "v_" + std::to_string(*rc::gen::inRange(0, (999) + 1));
        Sign sign = random_sign();

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(var_name, sign).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query multiple times — all should return the same result
        Tribool r1 = qi.query_positive(expr).value();
        Tribool r2 = qi.query_positive(expr).value();
        Tribool r3 = qi.query_positive(expr).value();

        RC_ASSERT(r1 == r2);
        RC_ASSERT(r2 == r3);

        // Also test other query types for consistency
        Tribool n1 = qi.query_negative(expr).value();
        Tribool n2 = qi.query_negative(expr).value();
        RC_ASSERT(n1 == n2);

        Tribool nn1 = qi.query_nonnegative(expr).value();
        Tribool nn2 = qi.query_nonnegative(expr).value();
        RC_ASSERT(nn1 == nn2);
    }));
}

TEST(LmcasAssumptionCache, InvalidationAllowsNewResults) {
    EXPECT_TRUE(rc::check("New assumptions are observed without manual invalidation", []() {
        std::string var_name = "z_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        // Create a context where the variable initially has no sign
        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // First query — no sign declared, should be Unknown
        Tribool result_before = qi.query_positive(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Now declare the variable as Positive in the context
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());

        // After invalidation, the query should recompute with new assumptions
        Tribool result_after = qi.query_positive(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    }));
}

TEST(LmcasAssumptionCache, InvalidationOnScopePush) {
    EXPECT_TRUE(rc::check("Scope push and shadowing are observed automatically", []() {
        std::string var_name = "p_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query in initial scope
        Tribool result_initial = qi.query_positive(expr).value();
        RC_ASSERT(result_initial == Tribool::True);

        // Push a new scope and declare the variable as Negative (shadows parent)
        ctx.push();
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Negative).has_value());

        // Query should now reflect the new scope's declaration
        Tribool result_pushed = qi.query_positive(expr).value();
        RC_ASSERT(result_pushed == Tribool::False);

        Tribool result_neg = qi.query_negative(expr).value();
        RC_ASSERT(result_neg == Tribool::True);
    }));
}

TEST(LmcasAssumptionCache, InvalidationOnScopePop) {
    EXPECT_TRUE(rc::check("Scope pop restores cached parent facts automatically", []() {
        std::string var_name = "q_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query in root scope
        Tribool result_root = qi.query_positive(expr).value();
        RC_ASSERT(result_root == Tribool::True);

        // Push scope, declare Negative
        ctx.push();
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Negative).has_value());

        Tribool result_child = qi.query_positive(expr).value();
        RC_ASSERT(result_child == Tribool::False);

        // Pop scope — should revert to parent's Positive
        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";

        Tribool result_popped = qi.query_positive(expr).value();
        RC_ASSERT(result_popped == Tribool::True);
    }));
}

TEST(LmcasAssumptionCache, InvalidationOnAssumeDomain) {
    EXPECT_TRUE(rc::check("Domain transactions update cached results automatically", []() {
        std::string var_name = "d_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Initially no domain declared — query_integer should be Unknown
        Tribool result_before = qi.query_integer(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Declare Integer domain
        EXPECT_TRUE(ctx.assume_domain(var_name, Domain::Integer).has_value());

        // After invalidation, query should reflect new domain
        Tribool result_after = qi.query_integer(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    }));
}

// --- Test: invalidate_cache on assume (add_relation) ---

TEST(LmcasAssumptionCache, InvalidationOnAssumeRelation) {
    EXPECT_TRUE(rc::check("Relation transactions update cached sign results automatically", []() {
        std::string var_name = "r_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Initially no sign — query_positive should be Unknown
        Tribool result_before = qi.query_positive(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Add relation: var > 0 (which derives Positive sign)
        auto var_node = LMCAS::detail::make_node<VariableNode>(var_name);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT);
        auto rel_expr = LMCAS::detail::expression_from_node(rel_node);
        EXPECT_TRUE(ctx.assume(rel_expr).has_value()) << "relation assumption succeeds";

        // After invalidation, query should reflect the new relation
        Tribool result_after = qi.query_positive(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    }));
}

// --- Test: Multiple invalidations maintain correctness ---

TEST(LmcasAssumptionCache, MultipleInvalidationsCorrect) {
    EXPECT_TRUE(rc::check("Successive transactions and scopes update the same query interface", []() {
        std::string var_name = "m_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Round 1: Unknown
        Tribool r1 = qi.query_positive(expr).value();
        RC_ASSERT(r1 == Tribool::Unknown);

        // Round 2: Declare Positive
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());
        Tribool r2 = qi.query_positive(expr).value();
        RC_ASSERT(r2 == Tribool::True);

        // Round 3: Push scope, declare Negative
        ctx.push();
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Negative).has_value());
        Tribool r3 = qi.query_positive(expr).value();
        RC_ASSERT(r3 == Tribool::False);

        // Round 4: Pop scope, revert to Positive
        EXPECT_TRUE(ctx.pop().has_value()) << "scope pop succeeds";
        Tribool r4 = qi.query_positive(expr).value();
        RC_ASSERT(r4 == Tribool::True);
    }));
}

// --- Test: Cache works across different property types independently ---

TEST(LmcasAssumptionCache, CacheDifferentPropertyTypes) {
    EXPECT_TRUE(rc::check("Different property queries on the same expression are cached independently", []() {
        std::string var_name = "t_" + std::to_string(*rc::gen::inRange(0, (999) + 1));

        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(var_name, Sign::Positive).has_value());
        EXPECT_TRUE(ctx.assume_domain(var_name, Domain::Integer).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query different properties
        Tribool pos = qi.query_positive(expr).value();
        Tribool neg = qi.query_negative(expr).value();
        Tribool nonneg = qi.query_nonnegative(expr).value();
        Tribool integer = qi.query_integer(expr).value();
        Tribool real = qi.query_real(expr).value();

        // Verify correctness
        RC_ASSERT(pos == Tribool::True);
        RC_ASSERT(neg == Tribool::False);
        RC_ASSERT(nonneg == Tribool::True);
        RC_ASSERT(integer == Tribool::True);
        RC_ASSERT(real == Tribool::True);

        // Repeat — should return same cached results
        RC_ASSERT(qi.query_positive(expr).value() == pos);
        RC_ASSERT(qi.query_negative(expr).value() == neg);
        RC_ASSERT(qi.query_nonnegative(expr).value() == nonneg);
        RC_ASSERT(qi.query_integer(expr).value() == integer);
        RC_ASSERT(qi.query_real(expr).value() == real);

        // Invalidate and verify all still correct
        qi.invalidate_cache();
        RC_ASSERT(qi.query_positive(expr).value() == Tribool::True);
        RC_ASSERT(qi.query_negative(expr).value() == Tribool::False);
        RC_ASSERT(qi.query_nonnegative(expr).value() == Tribool::True);
        RC_ASSERT(qi.query_integer(expr).value() == Tribool::True);
        RC_ASSERT(qi.query_real(expr).value() == Tribool::True);
    }));
}

// --- Test: invalidate_cache clears ALL entries (not just one property type) ---

TEST(LmcasAssumptionCache, InvalidateClearsAllEntries) {
    EXPECT_TRUE(rc::check("invalidate_cache() clears cache for all expressions and property types", []() {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());
        EXPECT_TRUE(ctx.assume_sign("b", Sign::Negative).has_value());
        EXPECT_TRUE(ctx.assume_domain("a", Domain::Integer).has_value());

        QueryInterface qi(ctx);
        SymbolicExpr expr_a = make_var_expr("a");
        SymbolicExpr expr_b = make_var_expr("b");

        // Populate cache with multiple entries
        Tribool a_pos = qi.query_positive(expr_a).value();
        Tribool b_neg = qi.query_negative(expr_b).value();
        Tribool a_int = qi.query_integer(expr_a).value();

        RC_ASSERT(a_pos == Tribool::True);
        RC_ASSERT(b_neg == Tribool::True);
        RC_ASSERT(a_int == Tribool::True);

        // Invalidate — all entries cleared
        qi.invalidate_cache();

        // Queries still return correct results (recomputed)
        RC_ASSERT(qi.query_positive(expr_a).value() == Tribool::True);
        RC_ASSERT(qi.query_negative(expr_b).value() == Tribool::True);
        RC_ASSERT(qi.query_integer(expr_a).value() == Tribool::True);
    }));
}

TEST(LmcasAssumptionCache, RetainedStoreTransactions) {
    AssumptionContext ctx;
    QueryInterface query(ctx);
    const auto x = make_var_expr("x");
    auto &store = ctx.current_properties();
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "initial sign is unknown";
    EXPECT_TRUE(store.declare_sign("x", Sign::NonNegative).has_value()) << "weak sign commits";
    EXPECT_TRUE(query.query_nonnegative(x).value() == Tribool::True) << "retained reference update visible";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "weak sign is not strict";
    EXPECT_TRUE(store.declare_sign("x", Sign::Positive).has_value()) << "second commit succeeds";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "second retained-reference update visible";
    auto rejected = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE(!rejected && rejected.error().code == CasErrc::InvalidArgument) << "contradiction rejected";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "failed transaction retains facts";
    PropertyStore replacement;
    EXPECT_TRUE(replacement.declare_sign("x", Sign::Negative).has_value()) << "replacement is valid";
    store = replacement;
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::False) << "copy assignment invalidates target cache";
    store = PropertyStore{};
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "move assignment invalidates target cache";
}

TEST(LmcasAssumptionCache, RetainedRelationTransactions) {
    AssumptionContext ctx;
    QueryInterface query(ctx);
    const auto x = detail::expression_from_node(detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{make_var_node("x"), make_var_node("y")}));
    const auto one = detail::expression_from_node(detail::make_node<NumberNode>(BigInt(1)));
    auto &relations = ctx.current_relations();
    auto &properties = ctx.current_properties();
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "initial relation absent";
    EXPECT_TRUE(relations.add_relation(x, one, RelationOp::GT, properties).has_value()) << "relation commits";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "nonzero-right relation proves positivity";
    relations.clear();
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "clear removes cached proof";
    RelationStore replacement;
    EXPECT_TRUE(replacement.add_relation(x, one, RelationOp::GT, properties).has_value()) << "replacement relation commits";
    relations = replacement;
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "relation copy assignment observed";
    relations = RelationStore{};
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "relation move assignment observed";
    const auto p = make_var_expr("p");
    const auto zero = detail::expression_from_node(detail::make_node<NumberNode>(BigInt(0)));
    EXPECT_TRUE(relations.add_relation(p, zero, RelationOp::GT, properties).has_value()) << "derived sign commits";
    EXPECT_TRUE(query.query_positive(p).value() == Tribool::True) << "derived sign cached";
    auto failed = relations.add_relation(p, zero, RelationOp::LT, properties);
    EXPECT_TRUE(!failed && failed.error().code == CasErrc::InvalidArgument) << "contradictory relation rejected";
    EXPECT_TRUE(relations.has_relation(p, zero, RelationOp::GT) &&
                !relations.has_relation(p, zero, RelationOp::LT))
        << "failed relation leaves relation store intact";
    EXPECT_TRUE(query.query_positive(p).value() == Tribool::True) << "failed relation leaves property proof intact";
}

TEST(LmcasAssumptionCache, ContextAssignmentAndCopyIsolation) {
    const auto x = make_var_expr("x");
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "positive context";
    QueryInterface query(ctx);
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "cached positive";
    AssumptionContext other;
    EXPECT_TRUE(other.assume_sign("x", Sign::Negative).has_value()) << "negative context";
    ctx = other;
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::False) << "copy assignment cannot reuse generation";
    AssumptionContext copy(ctx);
    QueryInterface copy_query(copy);
    copy.current_properties() = PropertyStore{};
    EXPECT_TRUE(copy_query.query_positive(x).value() == Tribool::Unknown) << "copy changes independently";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::False) << "source unchanged by copy mutation";
    ctx = std::move(copy);
    EXPECT_TRUE(copy_query.query_positive(x).value() == copy.is_positive_checked(x).value()) << "move source cached query agrees with its current facts";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::Unknown) << "move assignment invalidates destination";
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "source fact before move construction";
    EXPECT_TRUE(query.query_positive(x).value() == Tribool::True) << "source proof cached";
    AssumptionContext moved(std::move(ctx));
    EXPECT_TRUE(moved.is_positive_checked(x).value() == Tribool::True) << "move constructor retains destination facts";
    EXPECT_TRUE(query.query_positive(x).value() == ctx.is_positive_checked(x).value()) << "move constructor invalidates source cached proof";
}

TEST(LmcasAssumptionCache, DepthPolicyInvalidatesUnknown) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real variable declared";
    auto node = make_var_node("x");
    for (int i = 0; i < 8; ++i) {
        node = detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp, std::vector<std::shared_ptr<const SymbolicNode>>{node});
    }
    const auto expression = detail::expression_from_node(node);
    ctx.set_max_query_depth(1);
    QueryInterface query(ctx);
    EXPECT_TRUE(query.query_positive(expression).value() == Tribool::Unknown) << "shallow query remains unknown";
    ctx.set_max_query_depth(32);
    EXPECT_TRUE(query.query_positive(expression).value() == Tribool::True) << "deeper policy recomputes proof";
}
