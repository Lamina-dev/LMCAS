/**
 * @file test_assumption_propagation.cpp
 * @brief Property tests for automatic propagation at query time (Task 9.3).
 *
 * Property tested:
 * - Property 14: Automatic propagation at query time
 *
 * For any Real variable, is_nonnegative(x²) SHALL be True.
 * For any Integer variable, is_integer(x²) SHALL be True.
 * For any Positive variable, query_positive(|x|) SHALL be True.
 * Propagation is lazy: no properties stored on x² at declaration time.
 *
 * Validates: Requirements 18.1, 18.3, 18.4
 *
 * Uses rapidcheck (header-only, vendored in tests/rapidcheck/) for
 * property-based testing with random input generation.
 */

#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;

// ============================================================
// Helpers: create AST nodes
// ============================================================

static std::shared_ptr<SymbolicNode> make_var(const std::string& name) {
    return std::make_shared<VariableNode>(name);
}

static std::shared_ptr<SymbolicNode> make_number(int val) {
    return std::make_shared<NumberNode>(BigInt(val));
}

static std::shared_ptr<SymbolicNode> make_power(
    std::shared_ptr<SymbolicNode> base,
    std::shared_ptr<SymbolicNode> exp) {
    return std::make_shared<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<SymbolicNode> make_abs(std::shared_ptr<SymbolicNode> arg) {
    return std::make_shared<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<SymbolicNode>>{std::move(arg)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<SymbolicNode> node) {
    SymbolicExpr expr;
    expr.root = std::move(node);
    return expr;
}

/// Build x² expression: PowerNode(VariableNode(name), NumberNode(2))
static SymbolicExpr make_x_squared(const std::string& var_name) {
    auto var = make_var(var_name);
    auto two = make_number(2);
    auto pow = make_power(var, two);
    return wrap_expr(pow);
}

/// Build |x| expression: FunctionNode(Abs, [VariableNode(name)])
static SymbolicExpr make_abs_x(const std::string& var_name) {
    auto var = make_var(var_name);
    auto abs_node = make_abs(var);
    return wrap_expr(abs_node);
}

/// Generate a random variable name for property tests
static std::string random_var_name() {
    static const std::vector<std::string> prefixes = {"x", "y", "z", "a", "b", "t", "u", "v", "w"};
    std::string prefix = rc::gen::elementOf(prefixes);
    return prefix + "_" + std::to_string(rc::gen::inRange(0, 999));
}

/// Generate a random Real-or-more-specific domain (Real, Algebraic, Rational, Integer, Natural, PositiveInt)
static Domain random_real_domain() {
    static const std::vector<Domain> real_domains = {
        Domain::Real, Domain::Algebraic, Domain::Rational,
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(real_domains);
}

/// Generate a random Integer-or-more-specific domain (Integer, Natural, PositiveInt)
static Domain random_integer_domain() {
    static const std::vector<Domain> int_domains = {
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(int_domains);
}

// ============================================================
// Property 14: For any variable declared Real, query_nonnegative on x² returns True
// **Validates: Requirements 18.1**
// ============================================================

static void test_property14_real_var_x_squared_nonnegative() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 14: Real variable x² is non-negative");

    rc::check("For any variable declared Real (or more specific), "
              "query_nonnegative on x² returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        InferenceEngine engine(ctx);

        auto x_squared = make_x_squared(var_name);

        // x² should be non-negative for any Real variable
        RC_ASSERT(engine.query_nonnegative(x_squared) == Tribool::True);
    });
}

// ============================================================
// Property 14: For any Integer variable, query_integer on x² returns True
// **Validates: Requirements 18.3**
// ============================================================

static void test_property14_integer_var_x_squared_integer() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 14: Integer variable x² has Integer domain");

    rc::check("For any variable declared Integer (or more specific), "
              "query_integer on x² returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_integer_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        InferenceEngine engine(ctx);

        auto x_squared = make_x_squared(var_name);

        // x² should have Integer domain when x is Integer
        RC_ASSERT(engine.query_integer(x_squared) == Tribool::True);
    });
}

// ============================================================
// Property 14: For any Positive variable, query_positive on |x| returns True
// **Validates: Requirements 18.1 (propagation of sign through abs)**
// ============================================================

static void test_property14_positive_var_abs_positive() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 14: Positive variable |x| is positive");

    rc::check("For any variable declared Positive and Real, "
              "query_positive on |x| returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        ctx.assume_sign(var_name, Sign::Positive);
        InferenceEngine engine(ctx);

        auto abs_x = make_abs_x(var_name);

        // |x| should be positive when x is positive and real
        RC_ASSERT(engine.query_positive(abs_x) == Tribool::True);
        RC_ASSERT(engine.query_nonnegative(abs_x) == Tribool::True);
    });
}

// ============================================================
// Property 14: Propagation is lazy — no properties stored on x² at declaration time
// **Validates: Requirements 18.4**
// ============================================================

static void test_property14_propagation_is_lazy() {
    TEST_CASE("Feature: assumption-system-enhancements, Property 14: Propagation is lazy (query-time only)");

    rc::check("After declaring a variable Real, the PropertyStore does not "
              "contain any entry for x² — propagation happens at query time only", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);

        // The PropertyStore should only have the declared variable, not x²
        // Check that the property store does not have a symbol named "x²" or
        // any composite expression stored eagerly
        const PropertyStore& store = ctx.current_properties();

        // The store should have the variable we declared
        RC_ASSERT(store.get_domain(var_name) == domain);

        // A synthetic name like "var_name^2" should NOT be in the store
        // (since propagation is lazy, no derived expressions are stored)
        std::string squared_name = var_name + "^2";
        RC_ASSERT(store.get_domain(squared_name) == Domain::Complex); // default = not stored

        // But querying through the InferenceEngine should still work
        InferenceEngine engine(ctx);
        auto x_squared = make_x_squared(var_name);
        RC_ASSERT(engine.query_nonnegative(x_squared) == Tribool::True);
    });
}

// ============================================================
// main
// ============================================================

int main() {
    // Property 14: Automatic propagation at query time
    test_property14_real_var_x_squared_nonnegative();
    test_property14_integer_var_x_squared_integer();
    test_property14_positive_var_abs_positive();
    test_property14_propagation_is_lazy();

    return TEST_REPORT();
}
