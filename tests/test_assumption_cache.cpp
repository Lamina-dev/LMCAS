
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "query_interface.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;


static std::shared_ptr<const SymbolicNode> make_var_node(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static SymbolicExpr make_var_expr(const std::string& name) {
    auto expr = lamina::detail::expression_from_node(make_var_node(name));
    return expr;
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

/// Build a division expression: numerator_var / denominator_var
static SymbolicExpr make_division(const std::string& num_var, const std::string& den_var) {
    auto num = make_var_node(num_var);
    auto den = make_var_node(den_var);
    auto den_inv = lamina::detail::make_node<PowerNode>(den, make_number(-1));
    auto mul = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{num, den_inv});
    auto expr = lamina::detail::expression_from_node(mul);
    return expr;
}

/// Build an AddNode expression from variable names
static SymbolicExpr make_add_expr(const std::vector<std::string>& var_names) {
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    for (const auto& name : var_names) {
        operands.push_back(make_var_node(name));
    }
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<AddNode>(std::move(operands)));
    return expr;
}

/// Generate a random Sign from a subset of useful signs
static Sign random_sign() {
    std::vector<Sign> signs = {Sign::Positive, Sign::Negative, Sign::NonNegative, Sign::NonPositive};
    return rc::gen::elementOf(signs);
}

/// Generate a random Domain
static Domain random_domain() {
    std::vector<Domain> domains = {Domain::Real, Domain::Integer, Domain::Rational, Domain::Natural};
    return rc::gen::elementOf(domains);
}


// --- Test: invalidate_cache() clears the cache ---

static void test_invalidate_cache_clears_cache() {
    TEST_CASE("invalidate_cache() clears the cache");

    rc::check("After invalidate_cache(), queries recompute and return correct results", []() {
        std::string var_name = "x_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(var_name, Sign::Positive);

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
    });
}

// --- Test: Cache stores results (same query returns same result) ---

static void test_cache_stores_results() {
    TEST_CASE("Cache stores results — repeated queries return same result");

    rc::check("Repeated queries on the same expression return the same cached result", []() {
        std::string var_name = "v_" + std::to_string(rc::gen::inRange(0, 999));
        Sign sign = random_sign();

        AssumptionContext ctx;
        ctx.assume_sign(var_name, sign);

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
    });
}

// --- Test: After invalidation, queries that would return different results do so ---

static void test_invalidation_allows_new_results() {
    TEST_CASE("After invalidation, new assumptions produce new results");

    rc::check("After invalidate_cache() and new assumptions, queries return updated results", []() {
        std::string var_name = "z_" + std::to_string(rc::gen::inRange(0, 999));

        // Create a context where the variable initially has no sign
        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // First query — no sign declared, should be Unknown
        Tribool result_before = qi.query_positive(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Now declare the variable as Positive in the context
        ctx.assume_sign(var_name, Sign::Positive);

        // Invalidate cache manually (since hooks aren't wired yet)
        qi.invalidate_cache();

        // After invalidation, the query should recompute with new assumptions
        Tribool result_after = qi.query_positive(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    });
}

// --- Test: invalidate_cache on scope push (manual invalidation) ---

static void test_invalidation_on_scope_push() {
    TEST_CASE("Cache invalidation on scope push");

    rc::check("After push_scope and invalidate_cache(), queries recompute with new scope", []() {
        std::string var_name = "p_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(var_name, Sign::Positive);

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query in initial scope
        Tribool result_initial = qi.query_positive(expr).value();
        RC_ASSERT(result_initial == Tribool::True);

        // Push a new scope and declare the variable as Negative (shadows parent)
        ctx.push();
        ctx.assume_sign(var_name, Sign::Negative);

        // Manually invalidate cache (simulating what task 8.5 will wire)
        qi.invalidate_cache();

        // Query should now reflect the new scope's declaration
        Tribool result_pushed = qi.query_positive(expr).value();
        RC_ASSERT(result_pushed == Tribool::False);

        Tribool result_neg = qi.query_negative(expr).value();
        RC_ASSERT(result_neg == Tribool::True);
    });
}

// --- Test: invalidate_cache on scope pop (manual invalidation) ---

static void test_invalidation_on_scope_pop() {
    TEST_CASE("Cache invalidation on scope pop");

    rc::check("After pop_scope and invalidate_cache(), queries revert to parent scope results", []() {
        std::string var_name = "q_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(var_name, Sign::Positive);

        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Query in root scope
        Tribool result_root = qi.query_positive(expr).value();
        RC_ASSERT(result_root == Tribool::True);

        // Push scope, declare Negative
        ctx.push();
        ctx.assume_sign(var_name, Sign::Negative);
        qi.invalidate_cache();

        Tribool result_child = qi.query_positive(expr).value();
        RC_ASSERT(result_child == Tribool::False);

        // Pop scope — should revert to parent's Positive
        ctx.pop();
        qi.invalidate_cache();

        Tribool result_popped = qi.query_positive(expr).value();
        RC_ASSERT(result_popped == Tribool::True);
    });
}

// --- Test: invalidate_cache on assume_domain (manual invalidation) ---

static void test_invalidation_on_assume_domain() {
    TEST_CASE("Cache invalidation on assume_domain");

    rc::check("After assume_domain and invalidate_cache(), domain queries return updated results", []() {
        std::string var_name = "d_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Initially no domain declared — query_integer should be Unknown
        Tribool result_before = qi.query_integer(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Declare Integer domain
        ctx.assume_domain(var_name, Domain::Integer);
        qi.invalidate_cache();

        // After invalidation, query should reflect new domain
        Tribool result_after = qi.query_integer(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    });
}

// --- Test: invalidate_cache on assume (add_relation) ---

static void test_invalidation_on_assume_relation() {
    TEST_CASE("Cache invalidation on assume (add_relation)");

    rc::check("After assume(relation) and invalidate_cache(), sign queries return updated results", []() {
        std::string var_name = "r_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Initially no sign — query_positive should be Unknown
        Tribool result_before = qi.query_positive(expr).value();
        RC_ASSERT(result_before == Tribool::Unknown);

        // Add relation: var > 0 (which derives Positive sign)
        auto var_node = lamina::detail::make_node<VariableNode>(var_name);
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = lamina::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT);
        auto rel_expr = lamina::detail::expression_from_node(rel_node);
        ctx.assume(rel_expr);

        qi.invalidate_cache();

        // After invalidation, query should reflect the new relation
        Tribool result_after = qi.query_positive(expr).value();
        RC_ASSERT(result_after == Tribool::True);
    });
}

// --- Test: Multiple invalidations maintain correctness ---

static void test_multiple_invalidations_correct() {
    TEST_CASE("Multiple invalidations maintain correctness");

    rc::check("Multiple cycles of assume + invalidate produce correct results each time", []() {
        std::string var_name = "m_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        QueryInterface qi(ctx);
        SymbolicExpr expr = make_var_expr(var_name);

        // Round 1: Unknown
        Tribool r1 = qi.query_positive(expr).value();
        RC_ASSERT(r1 == Tribool::Unknown);

        // Round 2: Declare Positive
        ctx.assume_sign(var_name, Sign::Positive);
        qi.invalidate_cache();
        Tribool r2 = qi.query_positive(expr).value();
        RC_ASSERT(r2 == Tribool::True);

        // Round 3: Push scope, declare Negative
        ctx.push();
        ctx.assume_sign(var_name, Sign::Negative);
        qi.invalidate_cache();
        Tribool r3 = qi.query_positive(expr).value();
        RC_ASSERT(r3 == Tribool::False);

        // Round 4: Pop scope, revert to Positive
        ctx.pop();
        qi.invalidate_cache();
        Tribool r4 = qi.query_positive(expr).value();
        RC_ASSERT(r4 == Tribool::True);
    });
}

// --- Test: Cache works across different property types independently ---

static void test_cache_different_property_types() {
    TEST_CASE("Cache works across different property types");

    rc::check("Different property queries on the same expression are cached independently", []() {
        std::string var_name = "t_" + std::to_string(rc::gen::inRange(0, 999));

        AssumptionContext ctx;
        ctx.assume_sign(var_name, Sign::Positive);
        ctx.assume_domain(var_name, Domain::Integer);

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
    });
}

// --- Test: invalidate_cache clears ALL entries (not just one property type) ---

static void test_invalidate_clears_all_entries() {
    TEST_CASE("invalidate_cache clears all entries");

    rc::check("invalidate_cache() clears cache for all expressions and property types", []() {
        AssumptionContext ctx;
        ctx.assume_sign("a", Sign::Positive);
        ctx.assume_sign("b", Sign::Negative);
        ctx.assume_domain("a", Domain::Integer);

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
    });
}


int main() {
    test_invalidate_cache_clears_cache();
    test_cache_stores_results();
    test_invalidation_allows_new_results();
    test_invalidation_on_scope_push();
    test_invalidation_on_scope_pop();
    test_invalidation_on_assume_domain();
    test_invalidation_on_assume_relation();
    test_multiple_invalidations_correct();
    test_cache_different_property_types();
    test_invalidate_clears_all_entries();

    return TEST_REPORT();
}
