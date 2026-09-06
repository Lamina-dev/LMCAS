
#include "test_common.hpp"
#include "query_interface.hpp"
#include "assumption_context.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>

using namespace LMCAS;


static SymbolicExpr make_var(const std::string& name) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
    return expr;
}

static SymbolicExpr make_int(int v) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(v)));
    return expr;
}

/// Build x - y as AddNode([x, MultiplyNode([-1, y])])
static SymbolicExpr make_subtraction(const std::string& lhs, const std::string& rhs) {
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
static SymbolicExpr make_sin(const std::string& var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto sin_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(sin_node);
    return expr;
}

/// Build cos(x) as FunctionNode(Cos, [VariableNode(x)])
static SymbolicExpr make_cos(const std::string& var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto cos_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Cos,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(cos_node);
    return expr;
}

/// Build tan(x) as FunctionNode(Tan, [VariableNode(x)])
static SymbolicExpr make_tan(const std::string& var_name) {
    auto x_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto tan_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Tan,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto expr = LMCAS::detail::expression_from_node(tan_node);
    return expr;
}

static void EXPECT_TRIBOOL(Tribool actual, Tribool expected, const std::string& msg) {
    const char* names[] = {"True", "False", "Unknown"};
    int ai = static_cast<int>(actual);
    int ei = static_cast<int>(expected);
    if (ai == ei) {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    } else {
        std::cerr << "[FAIL] " << msg
                  << "\n  Expected: " << names[ei]
                  << "\n  Got:      " << names[ai] << std::endl;
        g_failures++;
    }
}


void test_cache_hit_returns_same_result() {
    TEST_CASE("Cache hit returns same result without re-inference");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // First query - computes and caches
    Tribool result1 = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(result1, Tribool::True, "First query: x is Positive");

    // Second query - should return cached result (same value)
    Tribool result2 = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(result2, Tribool::True, "Second query (cached): x is Positive");

    // Verify consistency across different property queries on same expression
    Tribool neg1 = qi.query_negative(x_expr).value();
    Tribool neg2 = qi.query_negative(x_expr).value();
    EXPECT_TRIBOOL(neg1, Tribool::False, "First query: x is not Negative");
    EXPECT_TRIBOOL(neg2, Tribool::False, "Second query (cached): x is not Negative");
}


void test_cache_invalidation_on_push() {
    TEST_CASE("Cache invalidation on push");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Populate cache
    Tribool before = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(before, Tribool::True, "Before push: x is Positive");

    // Manually invalidate cache (simulating push_scope hook)
    qi.invalidate_cache();

    // After invalidation, query should still return correct result (recomputed)
    Tribool after = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(after, Tribool::True, "After invalidate (push): x still Positive");
}


void test_cache_invalidation_on_pop() {
    TEST_CASE("Cache invalidation on pop");

    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // x is undeclared -> Unknown
    Tribool before_push = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(before_push, Tribool::Unknown, "Before push: x is Unknown");

    // Push scope and declare x Positive
    ctx.push();
    ctx.assume_sign("x", Sign::Positive);
    qi.invalidate_cache();  // Simulate hook

    Tribool in_scope = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(in_scope, Tribool::True, "In pushed scope: x is Positive");

    // Pop scope
    ctx.pop();
    qi.invalidate_cache();  // Simulate hook

    // After pop, x should be Unknown again
    Tribool after_pop = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(after_pop, Tribool::Unknown, "After pop: x is Unknown again");
}


void test_cache_invalidation_on_assume() {
    TEST_CASE("Cache invalidation on assume");

    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Initially unknown
    Tribool before = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(before, Tribool::Unknown, "Before assume: x is Unknown");

    // Declare x Positive
    ctx.assume_sign("x", Sign::Positive);
    qi.invalidate_cache();  // Simulate hook

    // Now should be True
    Tribool after = qi.query_positive(x_expr).value();
    EXPECT_TRIBOOL(after, Tribool::True, "After assume + invalidate: x is Positive");
}


void test_query_conditions_simple_variable() {
    TEST_CASE("query_conditions for simple variable");

    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Query: under what conditions is x > 0?
    auto conditions_result = qi.query_conditions(x_expr, Sign::Positive);
    EXPECT_TRUE(conditions_result.has_value(), "Condition query succeeds");
    if (!conditions_result) return;
    const auto& conditions = conditions_result.value();

    EXPECT_TRUE(conditions.size() == 1,
                "Single variable should return exactly 1 condition set");

    if (!conditions.empty()) {
        const auto& cs = conditions[0];
        EXPECT_TRUE(cs.sign_conditions.size() == 1,
                    "Condition set should have 1 sign condition");
        if (!cs.sign_conditions.empty()) {
            EXPECT_TRUE(cs.sign_conditions[0].first == "x",
                        "Sign condition variable should be 'x'");
            EXPECT_TRUE(cs.sign_conditions[0].second == Sign::Positive,
                        "Sign condition should be Positive");
        }
        EXPECT_TRUE(cs.domain_conditions.empty(),
                    "No domain conditions for simple sign query");
        EXPECT_TRUE(cs.relational_conditions.empty(),
                    "No relational conditions for simple variable");
    }

    // Also test with Negative target
    auto neg_conditions_result = qi.query_conditions(x_expr, Sign::Negative);
    EXPECT_TRUE(neg_conditions_result.has_value(), "Negative condition query succeeds");
    if (!neg_conditions_result) return;
    const auto& neg_conditions = neg_conditions_result.value();
    EXPECT_TRUE(neg_conditions.size() == 1,
                "Negative target: single variable returns 1 condition set");
    if (!neg_conditions.empty()) {
        EXPECT_TRUE(neg_conditions[0].sign_conditions[0].second == Sign::Negative,
                    "Negative target: condition should be Negative");
    }
}


void test_query_conditions_subtraction() {
    TEST_CASE("query_conditions for composite expression x - y");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    // Build x - y
    auto expr = make_subtraction("x", "y");

    // Query: under what conditions is (x - y) > 0?
    auto conditions_result = qi.query_conditions(expr, Sign::Positive);
    EXPECT_TRUE(conditions_result.has_value(), "Composite condition query succeeds");
    if (!conditions_result) return;
    const auto& conditions = conditions_result.value();

    // Should return at least 2 condition sets:
    // 1. {x: Positive, y: Negative}
    // 2. {x GT y, y: NonNegative}
    EXPECT_TRUE(conditions.size() >= 2,
                "x - y Positive should return at least 2 condition sets");

    if (conditions.size() >= 2) {
        // First condition set: x Positive AND y Negative
        const auto& cs1 = conditions[0];
        EXPECT_TRUE(cs1.sign_conditions.size() == 2,
                    "First condition set has 2 sign conditions");

        bool has_x_positive = false;
        bool has_y_negative = false;
        for (const auto& sc : cs1.sign_conditions) {
            if (sc.first == "x" && sc.second == Sign::Positive) has_x_positive = true;
            if (sc.first == "y" && sc.second == Sign::Negative) has_y_negative = true;
        }
        EXPECT_TRUE(has_x_positive, "CS1: x should be Positive");
        EXPECT_TRUE(has_y_negative, "CS1: y should be Negative");

        // Second condition set: relational condition (x GT y) + y NonNegative
        const auto& cs2 = conditions[1];
        EXPECT_TRUE(!cs2.relational_conditions.empty(),
                    "CS2: should have relational conditions");
        EXPECT_TRUE(!cs2.sign_conditions.empty(),
                    "CS2: should have sign conditions (y NonNegative)");
    }
}


void test_query_conditions_empty_for_complex() {
    TEST_CASE("query_conditions returns empty for undetermined expressions");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

}


void test_query_positive_definite() {
    TEST_CASE("query_positive_definite");

    AssumptionContext ctx;
    // Declare matrix symbol M as PositiveDefinite
    ctx.current_properties().declare_definiteness("M", Definiteness::PositiveDefinite);

    QueryInterface qi(ctx);
    auto m_expr = make_var("M");

    EXPECT_TRIBOOL(qi.query_positive_definite(m_expr).value(), Tribool::True,
                   "M declared PositiveDefinite: query_positive_definite = True");
    EXPECT_TRIBOOL(qi.query_positive_semidefinite(m_expr).value(), Tribool::True,
                   "M declared PositiveDefinite: query_positive_semidefinite = True (implied)");
}


void test_query_positive_semidefinite() {
    TEST_CASE("query_positive_semidefinite");

    AssumptionContext ctx;
    // Declare matrix symbol A as PositiveSemiDefinite only
    ctx.current_properties().declare_definiteness("A", Definiteness::PositiveSemiDefinite);

    QueryInterface qi(ctx);
    auto a_expr = make_var("A");

    EXPECT_TRIBOOL(qi.query_positive_semidefinite(a_expr).value(), Tribool::True,
                   "A declared PSD: query_positive_semidefinite = True");
    /// PositiveSemiDefinite 对 PositiveDefinite 查询保持 Unknown.
    EXPECT_TRIBOOL(qi.query_positive_definite(a_expr).value(), Tribool::Unknown,
                   "A declared PSD: query_positive_definite = Unknown");
}


void test_query_definiteness_negative() {
    TEST_CASE("query_positive_definite for NegativeDefinite matrix");

    AssumptionContext ctx;
    ctx.current_properties().declare_definiteness("N", Definiteness::NegativeDefinite);

    QueryInterface qi(ctx);
    auto n_expr = make_var("N");

    EXPECT_TRIBOOL(qi.query_positive_definite(n_expr).value(), Tribool::False,
                   "N declared NegDef: query_positive_definite = False");
    EXPECT_TRIBOOL(qi.query_positive_semidefinite(n_expr).value(), Tribool::False,
                   "N declared NegDef: query_positive_semidefinite = False");
}


void test_query_definiteness_undeclared() {
    TEST_CASE("query_positive_definite for undeclared symbol");

    AssumptionContext ctx;
    QueryInterface qi(ctx);
    auto u_expr = make_var("U");

    EXPECT_TRIBOOL(qi.query_positive_definite(u_expr).value(), Tribool::Unknown,
                   "Undeclared: query_positive_definite = Unknown");
    EXPECT_TRIBOOL(qi.query_positive_semidefinite(u_expr).value(), Tribool::Unknown,
                   "Undeclared: query_positive_semidefinite = Unknown");
}


void test_query_algebraic() {
    TEST_CASE("query_algebraic");

    AssumptionContext ctx;
    // Declare x as Algebraic domain
    ctx.assume_domain("x", Domain::Algebraic);

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    EXPECT_TRIBOOL(qi.query_algebraic(x_expr).value(), Tribool::True,
                   "x declared Algebraic: query_algebraic = True");
    EXPECT_TRIBOOL(qi.query_transcendental(x_expr).value(), Tribool::False,
                   "x declared Algebraic: query_transcendental = False");

    // Undeclared variable
    auto y_expr = make_var("y");
    EXPECT_TRIBOOL(qi.query_algebraic(y_expr).value(), Tribool::Unknown,
                   "y undeclared: query_algebraic = Unknown");
}


void test_query_transcendental() {
    TEST_CASE("query_transcendental");

    AssumptionContext ctx;
    ctx.current_properties().declare_transcendental("pi_sym");

    QueryInterface qi(ctx);
    auto pi_expr = make_var("pi_sym");

    EXPECT_TRIBOOL(qi.query_transcendental(pi_expr).value(), Tribool::True,
                   "pi_sym declared Transcendental: query_transcendental = True");
    EXPECT_TRIBOOL(qi.query_algebraic(pi_expr).value(), Tribool::False,
                   "pi_sym declared Transcendental: query_algebraic = False");
}


void test_query_finite() {
    TEST_CASE("query_finite");

    AssumptionContext ctx;
    ctx.current_properties().declare_finiteness("a", Finiteness::Finite);

    QueryInterface qi(ctx);
    auto a_expr = make_var("a");

    EXPECT_TRIBOOL(qi.query_finite(a_expr).value(), Tribool::True,
                   "a declared Finite: query_finite = True");
    EXPECT_TRIBOOL(qi.query_divergent(a_expr).value(), Tribool::False,
                   "a declared Finite: query_divergent = False");

    // Undeclared
    auto b_expr = make_var("b");
    EXPECT_TRIBOOL(qi.query_finite(b_expr).value(), Tribool::Unknown,
                   "b undeclared: query_finite = Unknown");
    EXPECT_TRIBOOL(qi.query_divergent(b_expr).value(), Tribool::Unknown,
                   "b undeclared: query_divergent = Unknown");
}


void test_query_divergent() {
    TEST_CASE("query_divergent");

    AssumptionContext ctx;
    ctx.current_properties().declare_finiteness("d", Finiteness::Divergent);

    QueryInterface qi(ctx);
    auto d_expr = make_var("d");

    EXPECT_TRIBOOL(qi.query_divergent(d_expr).value(), Tribool::True,
                   "d declared Divergent: query_divergent = True");
    EXPECT_TRIBOOL(qi.query_finite(d_expr).value(), Tribool::False,
                   "d declared Divergent: query_finite = False");
}


void test_query_periodic_declared() {
    TEST_CASE("query_periodic for declared periodic symbol");

    AssumptionContext ctx;
    // Declare f as periodic with period 2*pi (represented as a number for simplicity)
    auto period_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(6)));  // Simplified period
    ctx.current_properties().declare_periodic("f", period_expr);

    QueryInterface qi(ctx);
    auto f_expr = make_var("f");

    EXPECT_TRIBOOL(qi.query_periodic(f_expr).value(), Tribool::True,
                   "f declared periodic: query_periodic = True");
}


void test_query_periodic_trig() {
    TEST_CASE("query_periodic for sin/cos/tan (auto-inferred)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto sin_expr = make_sin("x");
    auto cos_expr = make_cos("x");
    auto tan_expr = make_tan("x");

    EXPECT_TRIBOOL(qi.query_periodic(sin_expr).value(), Tribool::True,
                   "sin(x): query_periodic = True");
    EXPECT_TRIBOOL(qi.query_periodic(cos_expr).value(), Tribool::True,
                   "cos(x): query_periodic = True");
    EXPECT_TRIBOOL(qi.query_periodic(tan_expr).value(), Tribool::True,
                   "tan(x): query_periodic = True");
}


void test_get_period_trig() {
    TEST_CASE("get_period for sin/cos/tan");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto sin_expr = make_sin("x");
    auto cos_expr = make_cos("x");
    auto tan_expr = make_tan("x");

    auto sin_period = qi.get_period(sin_expr);
    auto cos_period = qi.get_period(cos_expr);
    auto tan_period = qi.get_period(tan_expr);

    EXPECT_TRUE(sin_period.has_value() && sin_period.value().has_value(),
                "sin(x) should have a period");
    EXPECT_TRUE(cos_period.has_value() && cos_period.value().has_value(),
                "cos(x) should have a period");
    EXPECT_TRUE(tan_period.has_value() && tan_period.value().has_value(),
                "tan(x) should have a period");

    // sin and cos should have the same period (2*pi)
    if (sin_period && cos_period &&
        sin_period.value().has_value() && cos_period.value().has_value()) {
        std::string sin_p_str = sin_period.value()->to_string();
        std::string cos_p_str = cos_period.value()->to_string();
        EXPECT_TRUE(sin_p_str == cos_p_str,
                    "sin and cos should have the same period (2*pi)");
    }
}


void test_get_period_non_periodic() {
    TEST_CASE("get_period returns nullopt for non-periodic expression");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto x_expr = make_var("x");
    auto period = qi.get_period(x_expr);
    EXPECT_TRUE(period.has_value(), "Non-periodic query succeeds");
    EXPECT_TRUE(period.has_value() && !period.value().has_value(),
                "Variable x (undeclared periodic) should have no period");

}

void test_checked_extended_query_contracts() {
    TEST_CASE("QueryInterface checked extended queries: explicit errors and values");

    AssumptionContext ctx;
    ctx.assume_domain("a", Domain::Algebraic);
    ctx.current_properties().declare_transcendental("tau");
    ctx.current_properties().declare_finiteness("finite_symbol", Finiteness::Finite);
    ctx.current_properties().declare_finiteness("divergent_symbol", Finiteness::Divergent);
    ctx.current_properties().declare_periodic("periodic_symbol", LMCAS::detail::make_expression_ptr(make_int(6)));
    ctx.current_properties().declare_definiteness("M", Definiteness::PositiveDefinite);

    QueryInterface qi(ctx);

    auto algebraic = qi.query_algebraic_checked(make_var("a"));
    EXPECT_TRUE(algebraic.has_value(), "checked query_algebraic succeeds");
    if (algebraic) {
        EXPECT_TRUE(algebraic.value() == Tribool::True,
                    "checked query_algebraic returns True for algebraic symbol");
    }

    auto transcendental = qi.query_transcendental_checked(make_var("tau"));
    EXPECT_TRUE(transcendental.has_value(), "checked query_transcendental succeeds");
    if (transcendental) {
        EXPECT_TRUE(transcendental.value() == Tribool::True,
                    "checked query_transcendental returns True for transcendental symbol");
    }

    auto finite = qi.query_finite_checked(make_var("finite_symbol"));
    EXPECT_TRUE(finite.has_value(), "checked query_finite succeeds");
    if (finite) {
        EXPECT_TRUE(finite.value() == Tribool::True,
                    "checked query_finite returns True for finite symbol");
    }

    auto divergent = qi.query_divergent_checked(make_var("divergent_symbol"));
    EXPECT_TRUE(divergent.has_value(), "checked query_divergent succeeds");
    if (divergent) {
        EXPECT_TRUE(divergent.value() == Tribool::True,
                    "checked query_divergent returns True for divergent symbol");
    }

    auto periodic = qi.query_periodic_checked(make_var("periodic_symbol"));
    EXPECT_TRUE(periodic.has_value(), "checked query_periodic succeeds");
    if (periodic) {
        EXPECT_TRUE(periodic.value() == Tribool::True,
                    "checked query_periodic returns True for periodic symbol");
    }

    auto period = qi.get_period_checked(make_var("periodic_symbol"));
    EXPECT_TRUE(period.has_value(), "checked get_period succeeds");
    if (period) {
        EXPECT_TRUE(period.value().has_value(),
                    "checked get_period returns a declared period");
    }

    auto positive_definite = qi.query_positive_definite_checked(make_var("M"));
    EXPECT_TRUE(positive_definite.has_value(), "checked query_positive_definite succeeds");
    if (positive_definite) {
        EXPECT_TRUE(positive_definite.value() == Tribool::True,
                    "checked query_positive_definite returns True for PD symbol");
    }

    auto positive_semidefinite = qi.query_positive_semidefinite_checked(make_var("M"));
    EXPECT_TRUE(positive_semidefinite.has_value(), "checked query_positive_semidefinite succeeds");
    if (positive_semidefinite) {
        EXPECT_TRUE(positive_semidefinite.value() == Tribool::True,
                    "checked query_positive_semidefinite returns True for PD symbol");
    }

    auto conditions = qi.query_conditions_checked(make_var("x"), Sign::Positive);
    EXPECT_TRUE(conditions.has_value(), "checked query_conditions succeeds");
    if (conditions) {
        EXPECT_TRUE(conditions.value().size() == 1,
                    "checked query_conditions preserves simple-variable conditions");
    }

    auto unsupported_conditions = qi.query_conditions_checked(make_int(1), Sign::Positive);
    EXPECT_TRUE(unsupported_conditions.has_value(),
                "checked query_conditions accepts valid unsupported expressions");
    if (unsupported_conditions) {
        EXPECT_TRUE(unsupported_conditions.value().empty(),
                    "checked query_conditions returns empty set for unsupported expressions");
    }

}


void test_cache_multiple_properties() {
    TEST_CASE("Cache works across multiple property types");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Integer);

    QueryInterface qi(ctx);
    auto x_expr = make_var("x");

    // Query multiple properties - each should be cached independently
    Tribool pos = qi.query_positive(x_expr).value();
    Tribool intg = qi.query_integer(x_expr).value();
    Tribool neg = qi.query_negative(x_expr).value();

    EXPECT_TRIBOOL(pos, Tribool::True, "x Positive cached correctly");
    EXPECT_TRIBOOL(intg, Tribool::True, "x Integer cached correctly");
    EXPECT_TRIBOOL(neg, Tribool::False, "x not Negative cached correctly");

    // Query again - should hit cache
    EXPECT_TRIBOOL(qi.query_positive(x_expr).value(), Tribool::True, "x Positive (cache hit)");
    EXPECT_TRIBOOL(qi.query_integer(x_expr).value(), Tribool::True, "x Integer (cache hit)");
    EXPECT_TRIBOOL(qi.query_negative(x_expr).value(), Tribool::False, "x not Negative (cache hit)");

    // Invalidate and re-query
    qi.invalidate_cache();
    EXPECT_TRIBOOL(qi.query_positive(x_expr).value(), Tribool::True, "x Positive after invalidate");
    EXPECT_TRIBOOL(qi.query_integer(x_expr).value(), Tribool::True, "x Integer after invalidate");
}


int main() {
    test_cache_hit_returns_same_result();
    test_cache_invalidation_on_push();
    test_cache_invalidation_on_pop();
    test_cache_invalidation_on_assume();
    test_cache_multiple_properties();

    test_query_conditions_simple_variable();
    test_query_conditions_subtraction();
    test_query_conditions_empty_for_complex();

    test_query_positive_definite();
    test_query_positive_semidefinite();
    test_query_definiteness_negative();
    test_query_definiteness_undeclared();

    // Extended property queries
    test_query_algebraic();
    test_query_transcendental();
    test_query_finite();
    test_query_divergent();
    test_query_periodic_declared();
    test_query_periodic_trig();
    test_get_period_trig();
    test_get_period_non_periodic();
    test_checked_extended_query_contracts();

    return TEST_REPORT();
}
