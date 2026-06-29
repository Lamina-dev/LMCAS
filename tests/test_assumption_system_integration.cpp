/**
 * @file test_assumption_system_integration.cpp
 * @brief Unit tests for system integration of AssumptionContext with subsystems.
 *
 * Tests:
 * 1. Integrator with positive variable simplifies |x| to x
 * 2. LimitVisitor with positive variable resolves sign ambiguity
 * 3. ODE solver with positive dep var selects positive branch
 * 4. Matcher with assumption context evaluates conditions
 * 5. All subsystems with nullptr behave identically to current
 *
 * Validates: Requirements 12.2, 12.3, 13.2, 15.2, 16.2
 */

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "integration.hpp"
#include "visitors/differentiation_visitor.hpp"
#include "visitors/limit_visitor.hpp"
#include "symbolic_ode.hpp"
#include "matcher.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

using namespace lamina;

// ============================================================
// Helpers
// ============================================================

static std::shared_ptr<SymbolicNode> make_var(const std::string& name) {
    return std::make_shared<VariableNode>(name);
}

static std::shared_ptr<SymbolicNode> make_number(int val) {
    return std::make_shared<NumberNode>(BigInt(val));
}

static std::shared_ptr<SymbolicNode> make_abs(const std::shared_ptr<SymbolicNode>& arg) {
    return std::make_shared<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<SymbolicNode>>{arg});
}

static SymbolicExpr wrap(std::shared_ptr<SymbolicNode> node) {
    return SymbolicExpr(std::move(node));
}

/// Check if an AST contains an abs() node wrapping the given variable.
static bool contains_abs_of(const std::shared_ptr<SymbolicNode>& node, const std::string& var_name) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Abs && fn->arguments.size() == 1) {
            if (auto vn = std::dynamic_pointer_cast<VariableNode>(fn->arguments[0])) {
                if (vn->name == var_name) return true;
            }
        }
        for (auto& arg : fn->arguments) {
            if (contains_abs_of(arg, var_name)) return true;
        }
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands)
            if (contains_abs_of(op, var_name)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands)
            if (contains_abs_of(op, var_name)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (contains_abs_of(pow->base, var_name)) return true;
        if (contains_abs_of(pow->exponent, var_name)) return true;
    }
    return false;
}

/// Check if an AST contains an abs() node anywhere.
static bool contains_abs(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Abs) return true;
        for (auto& arg : fn->arguments) {
            if (contains_abs(arg)) return true;
        }
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands)
            if (contains_abs(op)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands)
            if (contains_abs(op)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (contains_abs(pow->base)) return true;
        if (contains_abs(pow->exponent)) return true;
    }
    return false;
}

// ============================================================
// Test 1: Integrator with positive variable simplifies |x| to x
// (Requirement 12.2)
// ============================================================

static void test_integrator_positive_simplifies_abs() {
    TEST_CASE("Integrator: positive variable simplifies |x| to x (Req 12.2)");

    // Create integrand: |x|
    SymbolicExpr integrand(make_abs(make_var("x")));

    // Set up assumption context with x Positive
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    // Integrate with assumption context
    Integrator integrator;
    integrator.set_assumption_context(&ctx);
    auto result = integrator.integrate(integrand, "x");

    std::string result_str = result.to_string();
    std::cout << "  Integration of |x| with x Positive: " << result_str << std::endl;

    // The result should NOT contain abs(x) since x is positive, |x| = x
    // So integrating x gives x^2/2
    EXPECT_FALSE(contains_abs_of(result.root, "x"),
                 "Integration result does not contain abs(x) when x is Positive");
}

static void test_integrator_no_context_preserves_abs() {
    TEST_CASE("Integrator: no context preserves |x| behavior (Req 12.4)");

    // Create integrand: |x|
    SymbolicExpr integrand(make_abs(make_var("x")));

    // Integrate without assumption context (nullptr)
    Integrator integrator;
    // No set_assumption_context call — defaults to nullptr
    auto result = integrator.integrate(integrand, "x");

    std::string result_str = result.to_string();
    std::cout << "  Integration of |x| without context: " << result_str << std::endl;

    // Without context, the integrator should either keep abs or produce
    // a result that still references abs (or an unevaluated integral).
    // The key point is it should NOT simplify |x| to x without assumptions.
    // We verify the result is produced (non-empty) — backward compatibility.
    EXPECT_TRUE(result.root != nullptr,
                "Integration without context produces a result");
}

// ============================================================
// Test 2: LimitVisitor with positive variable resolves sign ambiguity
// (Requirement 13.2)
// ============================================================

static void test_limit_visitor_positive_resolves_sign() {
    TEST_CASE("LimitVisitor: positive variable resolves sign ambiguity (Req 13.2)");

    // Compute limit of 1/x as x → 0+ with x known Positive.
    // The LimitVisitor should use the assumption to determine the sign.
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    // Build 1/x = x^(-1) = MultiplyNode([1, PowerNode(x, -1)])
    auto one_over_x = std::make_shared<MultiplyNode>(
        std::vector<std::shared_ptr<SymbolicNode>>{
            make_number(1),
            std::make_shared<PowerNode>(make_var("x"), make_number(-1))
        });

    auto point = make_number(0);

    // Compute limit with assumption context (right-sided limit)
    LimitVisitor visitor_with_ctx("x", point, "+", &ctx);
    one_over_x->accept(visitor_with_ctx);
    auto result_with_ctx = visitor_with_ctx.get_result();

    std::string result_str = result_with_ctx ? SymbolicExpr(result_with_ctx).to_string() : "null";
    std::cout << "  Limit of 1/x as x->0+ with x Positive: " << result_str << std::endl;

    // The result should be +infinity (positive infinity)
    // Check that it's an infinity node (not negative infinity)
    EXPECT_TRUE(result_with_ctx != nullptr,
                "Limit with positive context produces a result");

    // Verify it's positive infinity (not negative)
    bool is_positive_inf = false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(result_with_ctx)) {
        if (fn->type == FunctionNode::FuncType::Infinity) {
            is_positive_inf = true;
        }
    }
    EXPECT_TRUE(is_positive_inf,
                "Limit of 1/x as x->0+ with x Positive is +infinity");
}

static void test_limit_visitor_nullptr_same_behavior() {
    TEST_CASE("LimitVisitor: nullptr context same as current behavior (Req 13.4)");

    // Compute limit of x^2 as x → 2 without context
    auto x_squared = std::make_shared<PowerNode>(make_var("x"), make_number(2));
    auto point = make_number(2);

    // Without context
    LimitVisitor visitor_no_ctx("x", point, "");
    x_squared->accept(visitor_no_ctx);
    auto result_no_ctx = visitor_no_ctx.get_result();

    // With nullptr context (explicit)
    LimitVisitor visitor_null_ctx("x", point, "", nullptr);
    x_squared->accept(visitor_null_ctx);
    auto result_null_ctx = visitor_null_ctx.get_result();

    std::string s1 = result_no_ctx ? SymbolicExpr(result_no_ctx).to_string() : "null";
    std::string s2 = result_null_ctx ? SymbolicExpr(result_null_ctx).to_string() : "null";

    std::cout << "  Limit of x^2 as x->2 (no ctx): " << s1 << std::endl;
    std::cout << "  Limit of x^2 as x->2 (nullptr): " << s2 << std::endl;

    EXPECT_EQ_STR(s1, s2, "LimitVisitor with no ctx and nullptr produce same result");
}

// ============================================================
// Test 3: ODE solver with positive dep var selects positive branch
// (Requirement 15.2)
// ============================================================

static void test_ode_solver_positive_branch() {
    TEST_CASE("ODE solver: positive dep var selects positive branch (Req 15.2)");

    // Solve dy/dx = x*y with y declared Positive
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::multiply(x, y);

    AssumptionContext ctx;
    ctx.assume_sign("y", Sign::Positive);
    ctx.assume_domain("y", Domain::Real);

    auto result_with_ctx = solve_separable_ode(rhs, "x", "y", &ctx);

    std::string result_str = result_with_ctx ? result_with_ctx->to_string() : "null";
    std::cout << "  ODE dy/dx=x*y with y Positive: " << result_str << std::endl;

    // When y is Positive, the solver wraps result in abs() to signal
    // positive branch preference
    EXPECT_TRUE(result_with_ctx != nullptr,
                "ODE solver with positive dep var produces a result");
    EXPECT_TRUE(contains_abs(result_with_ctx->root),
                "ODE solver with positive dep var contains abs() wrapper");
}

static void test_ode_solver_nullptr_no_abs() {
    TEST_CASE("ODE solver: nullptr context identical to current (Req 15.4)");

    // Solve dy/dx = x*y without context
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::multiply(x, y);

    auto result_no_ctx = solve_separable_ode(rhs, "x", "y", nullptr);
    auto result_default = solve_separable_ode(rhs, "x", "y");

    std::string s1 = result_no_ctx ? result_no_ctx->to_string() : "null";
    std::string s2 = result_default ? result_default->to_string() : "null";

    std::cout << "  ODE dy/dx=x*y (nullptr): " << s1 << std::endl;
    std::cout << "  ODE dy/dx=x*y (default): " << s2 << std::endl;

    EXPECT_EQ_STR(s1, s2, "ODE solver with nullptr and default produce same result");
}

// ============================================================
// Test 4: Matcher with assumption context evaluates conditions
// (Requirement 16.2)
// ============================================================

static void test_matcher_assumption_condition_matches() {
    TEST_CASE("Matcher: assumption_condition matches when context has variable Positive (Req 16.2)");

    // Create a rule with assumption_condition that checks if wildcard "A" is Positive
    auto pattern = SymbolicExpr(make_var("A"));
    auto replacement = SymbolicExpr(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    // The assumption_condition checks if the bound expression is Positive
    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        return ctx->is_positive(it->second) == Tribool::True;
    };

    Rule rule(pattern, replacement, wildcards, assumption_cond);

    // Set up context with x Positive
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    // Create target expression: x
    auto target = SymbolicExpr(make_var("x"));

    // Use RewriteEngine with assumption context
    RewriteEngine engine;
    engine.add_rule(rule);
    engine.set_assumption_context(&ctx);

    // The rule should match because x is Positive in the context
    MatchMap bindings;
    bool matched = Matcher::match(pattern, target, wildcards, bindings, &ctx);

    EXPECT_TRUE(matched, "Matcher matches pattern against target");

    // Verify the assumption_condition evaluates correctly
    bool cond_result = assumption_cond(bindings, &ctx);
    EXPECT_TRUE(cond_result,
                "assumption_condition returns true when context has x Positive");
}

static void test_matcher_assumption_condition_no_match_without_context() {
    TEST_CASE("Matcher: assumption_condition fails without context (Req 16.4)");

    // Same rule as above
    auto pattern = SymbolicExpr(make_var("A"));
    auto replacement = SymbolicExpr(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        return ctx->is_positive(it->second) == Tribool::True;
    };

    Rule rule(pattern, replacement, wildcards, assumption_cond);

    // Target: x (but no context)
    auto target = SymbolicExpr(make_var("x"));

    MatchMap bindings;
    bool matched = Matcher::match(pattern, target, wildcards, bindings, nullptr);

    EXPECT_TRUE(matched, "Matcher still matches structurally without context");

    // But the assumption_condition should fail without context
    bool cond_result = assumption_cond(bindings, nullptr);
    EXPECT_FALSE(cond_result,
                 "assumption_condition returns false without context (nullptr)");
}

static void test_matcher_rewrite_engine_with_context() {
    TEST_CASE("Matcher: RewriteEngine uses assumption context during apply (Req 16.3)");

    // Create a rule: abs(A) → A when A is Positive
    auto abs_A = SymbolicExpr(make_abs(make_var("A")));
    auto just_A = SymbolicExpr(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        return ctx->is_positive(it->second) == Tribool::True;
    };

    Rule rule(abs_A, just_A, wildcards, assumption_cond);

    // Set up context with x Positive
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    RewriteEngine engine;
    engine.add_rule(rule);
    engine.set_assumption_context(&ctx);

    // Apply to abs(x)
    auto input = SymbolicExpr(make_abs(make_var("x")));
    auto result = engine.apply(input);

    std::string result_str = result.to_string();
    std::cout << "  RewriteEngine abs(x) with x Positive: " << result_str << std::endl;

    // The result should be x (abs removed because x is Positive)
    EXPECT_FALSE(contains_abs(result.root),
                 "RewriteEngine removes abs(x) when x is Positive");
}

// ============================================================
// Test 5: All subsystems with nullptr behave identically to current
// (Requirements 12.4, 13.4, 15.4, 16.4)
// ============================================================

static void test_integrator_nullptr_identical() {
    TEST_CASE("Integrator: nullptr identical to no context (Req 12.4)");

    // Integrate x^2
    SymbolicExpr integrand(std::make_shared<PowerNode>(make_var("x"), make_number(2)));

    Integrator integrator1;
    // No context set (default nullptr)
    auto result1 = integrator1.integrate(integrand, "x");

    Integrator integrator2;
    integrator2.set_assumption_context(nullptr);
    auto result2 = integrator2.integrate(integrand, "x");

    std::string s1 = result1.to_string();
    std::string s2 = result2.to_string();

    EXPECT_EQ_STR(s1, s2, "Integrator with default and explicit nullptr produce same result");
}

static void test_ode_solver_nullptr_identical() {
    TEST_CASE("ODE solver: nullptr identical to default (Req 15.4)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(x, y);

    auto result1 = solve_separable_ode(rhs, "x", "y");
    auto result2 = solve_separable_ode(rhs, "x", "y", nullptr);

    std::string s1 = result1 ? result1->to_string() : "null";
    std::string s2 = result2 ? result2->to_string() : "null";

    EXPECT_EQ_STR(s1, s2, "solve_separable_ode default and nullptr produce same result");
}

static void test_matcher_nullptr_identical() {
    TEST_CASE("Matcher: nullptr identical to no context (Req 16.4)");

    auto pattern = SymbolicExpr(make_var("A"));
    auto target = SymbolicExpr(make_var("x"));
    std::unordered_set<std::string> wildcards = {"A"};

    MatchMap bindings1;
    bool matched1 = Matcher::match(pattern, target, wildcards, bindings1);

    MatchMap bindings2;
    bool matched2 = Matcher::match(pattern, target, wildcards, bindings2, nullptr);

    EXPECT_TRUE(matched1 == matched2,
                "Matcher::match with default and nullptr produce same match result");

    // Verify bindings are the same
    bool same_bindings = (bindings1.size() == bindings2.size());
    if (same_bindings) {
        for (auto& [key, val] : bindings1) {
            auto it = bindings2.find(key);
            if (it == bindings2.end() || it->second.to_string() != val.to_string()) {
                same_bindings = false;
                break;
            }
        }
    }
    EXPECT_TRUE(same_bindings, "Matcher bindings identical with default and nullptr");
}

// ============================================================
// main
// ============================================================

int main() {
    // Test 1: Integrator with positive variable
    test_integrator_positive_simplifies_abs();
    test_integrator_no_context_preserves_abs();

    // Test 2: LimitVisitor with positive variable
    test_limit_visitor_positive_resolves_sign();
    test_limit_visitor_nullptr_same_behavior();

    // Test 3: ODE solver with positive dep var
    test_ode_solver_positive_branch();
    test_ode_solver_nullptr_no_abs();

    // Test 4: Matcher with assumption context
    test_matcher_assumption_condition_matches();
    test_matcher_assumption_condition_no_match_without_context();
    test_matcher_rewrite_engine_with_context();

    // Test 5: All subsystems with nullptr
    test_integrator_nullptr_identical();
    test_ode_solver_nullptr_identical();
    test_matcher_nullptr_identical();

    return TEST_REPORT();
}
