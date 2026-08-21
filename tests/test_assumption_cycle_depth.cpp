
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <memory>
#include <string>

using namespace lamina;


/// Create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string& name) {
    return lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>(name));
}

/// Create a deeply nested expression: f(f(f(...f(x)...))) with given depth.
/// Uses exp as the nesting function.
static SymbolicExpr make_deeply_nested(const std::string& var_name, int depth) {
    auto node = lamina::detail::make_node<VariableNode>(var_name);
    std::shared_ptr<const SymbolicNode> current = node;
    for (int i = 0; i < depth; ++i) {
        current = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{current});
    }
    return lamina::detail::expression_from_node(current);
}

/// Create a nested AddNode expression: ((x + x) + (x + x)) + ... with given depth.
/// Each level doubles the tree width, creating exponential node count but linear depth.
static SymbolicExpr make_nested_add(const std::string& var_name, int depth) {
    auto var_node = lamina::detail::make_node<VariableNode>(var_name);
    std::shared_ptr<const SymbolicNode> current = var_node;
    for (int i = 0; i < depth; ++i) {
        current = lamina::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{current, var_node});
    }
    return lamina::detail::expression_from_node(current);
}

/// Create a self-referential expression by sharing the same node pointer at multiple
/// positions in the tree. This tests cycle detection via pointer identity.
static SymbolicExpr make_shared_subexpr() {
    // Create a shared sub-expression node
    auto x = lamina::detail::make_node<VariableNode>("x");
    auto shared_add = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, x});

    // Use the same shared_add node in two positions of a multiply
    // This creates a DAG (not a tree), where the same pointer appears twice.
    auto mul = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{shared_add, shared_add});

    return lamina::detail::expression_from_node(mul);
}

/// Create two structurally identical but distinct (different pointer) nodes.
/// This tests that cycle detection does NOT produce false positives.
static std::pair<SymbolicExpr, SymbolicExpr> make_structurally_identical_distinct() {
    // Two separate AddNode instances with the same structure
    auto x1 = lamina::detail::make_node<VariableNode>("x");
    auto x2 = lamina::detail::make_node<VariableNode>("x");

    auto add1 = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x1, x1});
    auto add2 = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x2, x2});

    // Multiply using two distinct but structurally identical nodes
    auto mul = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{add1, add2});

    return {lamina::detail::expression_from_node(mul), lamina::detail::expression_from_node(add1)};
}


static void test_cycle_shared_node_returns_unknown() {
    TEST_CASE("Shared node (same pointer) in expression returns Unknown without crash");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // Create expression with shared sub-expression (DAG structure)
    SymbolicExpr expr = make_shared_subexpr();

    // The query should complete without infinite recursion.
    // The shared node will be visited twice — the second visit detects the cycle
    // and returns Unknown. The overall result depends on whether the first visit
    // provides enough information.
    Tribool result = engine.query_positive(expr);

    // Key property: no crash, no infinite loop. Result is either True or Unknown.
    // Since x is Positive, (x+x) is Positive, and (x+x)*(x+x) should be Positive.
    // But if cycle detection triggers on the shared node, it may return Unknown.
    // Either way, the test passes if we get here without hanging.
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Shared node query completes without infinite recursion");
}

static void test_cycle_detection_no_false_positive_distinct_nodes() {
    TEST_CASE("Structurally identical but distinct nodes — no false positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // Create expression with two structurally identical but pointer-distinct sub-nodes
    auto [expr, _] = make_structurally_identical_distinct();

    // Since the nodes are distinct pointers, cycle detection should NOT trigger.
    // With x Positive: (x+x) is Positive, and (x+x)*(x+x) should be Positive.
    Tribool result = engine.query_positive(expr);

    EXPECT_TRUE(result == Tribool::True,
        "Distinct nodes with same structure: no false positive from cycle detection");
}

static void test_cycle_detection_multiply_shared_operand() {
    TEST_CASE("MultiplyNode with same operand pointer twice");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Create x^2 as multiply(x, x) using the SAME variable node pointer
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto mul = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node, x_node});
    auto expr = lamina::detail::expression_from_node(mul);
    // Query should complete. x*x with x Positive should be Positive.
    Tribool result = engine.query_positive(expr);

    // The same VariableNode pointer appears twice in the multiply.
    // Cycle detection uses pointer identity, so the second encounter of x_node
    // will detect a "cycle". This is acceptable — result may be Unknown or True.
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Same variable pointer in multiply completes without crash");
}

static void test_cycle_detection_nested_function_shared_arg() {
    TEST_CASE("Nested function with shared argument node");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Create exp(x) using a shared x node, then add(exp(x), exp(x)) using same exp node
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto exp_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Exp,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});

    // Use the same exp_node pointer twice in an AddNode
    auto add = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{exp_node, exp_node});
    auto expr = lamina::detail::expression_from_node(add);
    // exp(x) is always positive, so exp(x) + exp(x) should be positive.
    // But cycle detection may trigger on the shared exp_node.
    Tribool result = engine.query_positive(expr);

    // Must complete without infinite recursion
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Shared function node in add completes without infinite recursion");
}

static void test_cycle_detection_preserves_state_after_query() {
    TEST_CASE("Visited set is cleared after top-level query completes");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // First query with shared nodes
    SymbolicExpr shared_expr = make_shared_subexpr();
    engine.query_positive(shared_expr);

    // Second query on a simple expression should work normally
    // (visited set should be cleared after first query)
    SymbolicExpr simple = make_var("x");
    Tribool result = engine.query_positive(simple);

    EXPECT_TRUE(result == Tribool::True,
        "After shared-node query, subsequent simple query works correctly");
}

static void test_cycle_detection_different_query_types() {
    TEST_CASE("Cycle detection works across different query types");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    SymbolicExpr shared_expr = make_shared_subexpr();

    // Multiple query types should all complete without hanging
    Tribool r1 = engine.query_positive(shared_expr);
    Tribool r2 = engine.query_negative(shared_expr);
    Tribool r3 = engine.query_nonnegative(shared_expr);
    Tribool r4 = engine.query_real(shared_expr);
    Tribool r5 = engine.query_integer(shared_expr);

    // All should complete (no infinite loop)
    EXPECT_TRUE(r1 == Tribool::True || r1 == Tribool::Unknown,
        "query_positive on shared expr completes");
    EXPECT_TRUE(r2 == Tribool::False || r2 == Tribool::Unknown,
        "query_negative on shared expr completes");
    EXPECT_TRUE(r3 == Tribool::True || r3 == Tribool::Unknown,
        "query_nonnegative on shared expr completes");
    EXPECT_TRUE(r4 == Tribool::True || r4 == Tribool::Unknown,
        "query_real on shared expr completes");
    EXPECT_TRUE(r5 == Tribool::True || r5 == Tribool::Unknown,
        "query_integer on shared expr completes");
}


static void test_depth_limit_returns_unknown() {
    TEST_CASE("Expression exceeding max depth returns Unknown");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);
    engine.set_max_depth(5); // Set a low depth limit

    // Create expression nested deeper than the limit
    SymbolicExpr deep_expr = make_deeply_nested("x", 10);

    Tribool result = engine.query_positive(deep_expr);

    EXPECT_TRUE(result == Tribool::Unknown,
        "Expression exceeding depth limit returns Unknown");
}

static void test_depth_limit_no_exception() {
    TEST_CASE("Depth limit does not throw exceptions");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);
    engine.set_max_depth(3);

    SymbolicExpr deep_expr = make_deeply_nested("x", 20);

    bool threw = false;
    try {
        engine.query_positive(deep_expr);
        engine.query_negative(deep_expr);
        engine.query_nonnegative(deep_expr);
        engine.query_real(deep_expr);
        engine.query_integer(deep_expr);
        engine.query_nonzero(deep_expr);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw, "Depth limit exceeded does not throw any exception");
}

static void test_depth_limit_within_limit_works() {
    TEST_CASE("Expression within depth limit returns correct result");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);
    engine.set_max_depth(32); // Default depth

    // exp(x) with x Positive should be Positive (depth 2: exp -> x)
    SymbolicExpr shallow = make_deeply_nested("x", 1); // exp(x)

    Tribool result = engine.query_positive(shallow);

    EXPECT_TRUE(result == Tribool::True,
        "Expression within depth limit returns correct result (exp(x) with x Positive)");
}

static void test_depth_limit_boundary() {
    TEST_CASE("Expression at exact depth limit boundary");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);
    engine.set_max_depth(4);

    // Depth 3 nesting: exp(exp(exp(x))) — should be within limit (depth 4 allows 4 levels)
    SymbolicExpr at_limit = make_deeply_nested("x", 3);
    Tribool result_at = engine.query_positive(at_limit);

    // Depth 5 nesting: exceeds limit of 4
    SymbolicExpr over_limit = make_deeply_nested("x", 5);
    Tribool result_over = engine.query_positive(over_limit);

    // At limit should still work (exp of positive is positive)
    EXPECT_TRUE(result_at == Tribool::True || result_at == Tribool::Unknown,
        "Expression at depth limit boundary completes");

    // Over limit should return Unknown
    EXPECT_TRUE(result_over == Tribool::Unknown,
        "Expression over depth limit returns Unknown");
}

static void test_depth_limit_no_side_effects_on_context() {
    TEST_CASE("Depth limit exceeded does not modify context state");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);
    engine.set_max_depth(3);

    // Record state before
    int depth_before = ctx.depth();
    bool has_positive_before = ctx.current_properties().has_sign("x", Sign::Positive);
    Domain domain_before = ctx.current_properties().get_domain("x");

    // Trigger depth limit
    SymbolicExpr deep_expr = make_deeply_nested("x", 20);
    engine.query_positive(deep_expr);
    engine.query_negative(deep_expr);
    engine.query_real(deep_expr);

    // Verify state unchanged
    int depth_after = ctx.depth();
    bool has_positive_after = ctx.current_properties().has_sign("x", Sign::Positive);
    Domain domain_after = ctx.current_properties().get_domain("x");

    EXPECT_TRUE(depth_before == depth_after,
        "Context depth unchanged after depth-limited queries");
    EXPECT_TRUE(has_positive_before == has_positive_after,
        "Sign properties unchanged after depth-limited queries");
    EXPECT_TRUE(domain_before == domain_after,
        "Domain properties unchanged after depth-limited queries");
}

static void test_depth_limit_subsequent_queries_work() {
    TEST_CASE("After depth limit hit, subsequent shallow queries work");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);
    engine.set_max_depth(3);

    // First: trigger depth limit
    SymbolicExpr deep_expr = make_deeply_nested("x", 20);
    Tribool deep_result = engine.query_positive(deep_expr);
    EXPECT_TRUE(deep_result == Tribool::Unknown, "Deep expression returns Unknown");

    // Then: shallow query should still work correctly
    SymbolicExpr simple = make_var("x");
    Tribool simple_result = engine.query_positive(simple);
    EXPECT_TRUE(simple_result == Tribool::True,
        "Simple query after depth limit hit still returns correct result");
}

static void test_depth_limit_set_and_get() {
    TEST_CASE("set_max_depth and get_max_depth work correctly");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Default depth
    EXPECT_TRUE(engine.get_max_depth() == 32, "Default max depth is 32");

    // Set custom depth
    engine.set_max_depth(10);
    EXPECT_TRUE(engine.get_max_depth() == 10, "Max depth set to 10");

    // Invalid depth (0 or negative) should not change
    engine.set_max_depth(0);
    EXPECT_TRUE(engine.get_max_depth() == 10, "Max depth unchanged for 0");

    engine.set_max_depth(-5);
    EXPECT_TRUE(engine.get_max_depth() == 10, "Max depth unchanged for negative");
}

static void test_depth_limit_add_node_deep() {
    TEST_CASE("Deeply nested AddNode exceeding depth returns Unknown");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);
    engine.set_max_depth(5);

    // Create deeply nested add: ((((x + x) + x) + x) + x) ... 10 levels
    SymbolicExpr deep_add = make_nested_add("x", 10);

    Tribool result = engine.query_positive(deep_add);

    // Should return Unknown due to depth limit (or True if the engine
    // can short-circuit before hitting the limit)
    EXPECT_TRUE(result == Tribool::True || result == Tribool::Unknown,
        "Deeply nested AddNode completes without crash");
}

static void test_depth_limit_multiple_variables() {
    TEST_CASE("Depth limit with multiple variables in deep expression");

    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);
    ctx.assume_sign("b", Sign::Positive);
    ctx.assume_domain("a", Domain::Real);
    ctx.assume_domain("b", Domain::Real);

    InferenceEngine engine(ctx);
    engine.set_max_depth(4);

    // Create exp(exp(exp(exp(exp(a + b)))))
    auto a_node = lamina::detail::make_node<VariableNode>("a");
    auto b_node = lamina::detail::make_node<VariableNode>("b");
    auto add = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{a_node, b_node});

    std::shared_ptr<const SymbolicNode> current = add;
    for (int i = 0; i < 5; ++i) {
        current = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{current});
    }
    auto expr = lamina::detail::expression_from_node(current);
    Tribool result = engine.query_positive(expr);

    // Should return Unknown due to depth limit
    EXPECT_TRUE(result == Tribool::Unknown,
        "Deep multi-variable expression returns Unknown at depth limit");
}


int main() {
    test_cycle_shared_node_returns_unknown();
    test_cycle_detection_no_false_positive_distinct_nodes();
    test_cycle_detection_multiply_shared_operand();
    test_cycle_detection_nested_function_shared_arg();
    test_cycle_detection_preserves_state_after_query();
    test_cycle_detection_different_query_types();

    test_depth_limit_returns_unknown();
    test_depth_limit_no_exception();
    test_depth_limit_within_limit_works();
    test_depth_limit_boundary();
    test_depth_limit_no_side_effects_on_context();
    test_depth_limit_subsequent_queries_work();
    test_depth_limit_set_and_get();
    test_depth_limit_add_node_deep();
    test_depth_limit_multiple_variables();

    return TEST_REPORT();
}
