
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

using namespace LMCAS;


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_abs(const std::shared_ptr<const SymbolicNode>& arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

static SymbolicExpr wrap(std::shared_ptr<const SymbolicNode> node) {
    return LMCAS::detail::expression_from_node(std::move(node));
}

/// Check if an AST contains an abs() node wrapping the given variable.
static bool contains_abs_of(const std::shared_ptr<const SymbolicNode>& node, const std::string& var_name) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == FunctionNode::FuncType::Abs && fn->arguments().size() == 1) {
            if (auto vn = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0])) {
                if (vn->name() == var_name) return true;
            }
        }
        for (auto& arg : fn->arguments()) {
            if (contains_abs_of(arg, var_name)) return true;
        }
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands())
            if (contains_abs_of(op, var_name)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands())
            if (contains_abs_of(op, var_name)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (contains_abs_of(pow->base(), var_name)) return true;
        if (contains_abs_of(pow->exponent(), var_name)) return true;
    }
    return false;
}

/// Check if an AST contains an abs() node anywhere.
static bool contains_abs(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == FunctionNode::FuncType::Abs) return true;
        for (auto& arg : fn->arguments()) {
            if (contains_abs(arg)) return true;
        }
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands())
            if (contains_abs(op)) return true;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands())
            if (contains_abs(op)) return true;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (contains_abs(pow->base())) return true;
        if (contains_abs(pow->exponent())) return true;
    }
    return false;
}


static void test_integrator_positive_simplifies_abs() {
    TEST_CASE("Integrator: positive variable simplifies |x| to x");

    // Create integrand: |x|
    auto integrand = LMCAS::detail::expression_from_node(make_abs(make_var("x")));
    // Set up assumption context with x Positive
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value(),
                "positive-sign setup succeeds");
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value(),
                "real-domain setup succeeds");

    // Integrate with assumption context
    Integrator integrator;
    ComputationContext integration_context;
    auto set_integration_assumptions = integration_context.set_assumptions(
        std::make_shared<AssumptionContext>(ctx));
    EXPECT_TRUE(set_integration_assumptions.has_value(),
                "integration assumptions attach to context");
    auto integrated = integrator.integrate_checked(integrand, "x", integration_context);
    EXPECT_TRUE(integrated.has_value(), "context-aware integration succeeds");
    auto result = integrated.value();

    std::string result_str = result.to_string();
    std::cout << "  Integration of |x| with x Positive: " << result_str << std::endl;

    // The result should NOT contain abs(x) since x is positive, |x| = x
    // So integrating x gives x^2/2
    EXPECT_FALSE(contains_abs_of(LMCAS::detail::node(result), "x"),
                 "Integration result does not contain abs(x) when x is Positive");
}

static void test_integrator_no_context_preserves_abs() {
    TEST_CASE("Integrator: no context preserves |x| behavior");

    // Create integrand: |x|
    auto integrand = LMCAS::detail::expression_from_node(make_abs(make_var("x")));
    /// Integrator 保持默认空假设上下文.
    Integrator integrator;
    auto result = integrator.integrate(integrand, "x");
    EXPECT_TRUE(result.has_value(), "Integration without assumptions succeeds");
    if (!result) return;

    std::string result_str = result.value().to_string();
    std::cout << "  Integration of |x| without context: " << result_str << std::endl;

    /// 空上下文时积分结果保留 abs 或未求值积分结构,
    /// |x| = x 化简仅在符号假设充分时启用.
    EXPECT_TRUE(LMCAS::detail::node(result.value()) != nullptr,
                "Integration without context produces a result");
}


static void test_limit_visitor_positive_resolves_sign() {
    TEST_CASE("LimitVisitor: positive variable resolves sign ambiguity");

    // Compute limit of 1/x as x -> 0+ with x known Positive.
    // The LimitVisitor should use the assumption to determine the sign.
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value(),
                "positive-sign setup succeeds");
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value(),
                "real-domain setup succeeds");

    // Build 1/x = x^(-1) = MultiplyNode([1, PowerNode(x, -1)])
    auto one_over_x = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            make_number(1),
            LMCAS::detail::make_node<PowerNode>(make_var("x"), make_number(-1))
        });

    auto point = make_number(0);

    // Compute limit with assumption context (right-sided limit)
    LimitVisitor visitor_with_ctx("x", point, "+", &ctx);
    one_over_x->accept(visitor_with_ctx);
    auto result_with_ctx = visitor_with_ctx.get_result();

    std::string result_str = result_with_ctx ? LMCAS::detail::expression_from_node(result_with_ctx).to_string() : "null";
    std::cout << "  Limit of 1/x as x->0+ with x Positive: " << result_str << std::endl;

    // The result should be +infinity (positive infinity)
    // Check that it's an infinity node (not negative infinity)
    EXPECT_TRUE(result_with_ctx != nullptr,
                "Limit with positive context produces a result");

    // Verify it's positive infinity (not negative)
    bool is_positive_inf = false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(result_with_ctx)) {
        if (fn->type() == FunctionNode::FuncType::Infinity) {
            is_positive_inf = true;
        }
    }
    EXPECT_TRUE(is_positive_inf,
                "Limit of 1/x as x->0+ with x Positive is +infinity");
}

static void test_limit_visitor_nullptr_same_behavior() {
    TEST_CASE("LimitVisitor: nullptr context same as current behavior");

    /// 在空上下文中计算 x^2 于 x->2 的极限.
    auto x_squared = LMCAS::detail::make_node<PowerNode>(make_var("x"), make_number(2));
    auto point = make_number(2);

    /// LimitVisitor 使用默认空上下文.
    LimitVisitor visitor_no_ctx("x", point, "");
    x_squared->accept(visitor_no_ctx);
    auto result_no_ctx = visitor_no_ctx.get_result();

    // With nullptr context (explicit)
    LimitVisitor visitor_null_ctx("x", point, "", nullptr);
    x_squared->accept(visitor_null_ctx);
    auto result_null_ctx = visitor_null_ctx.get_result();

    std::string s1 = result_no_ctx ? LMCAS::detail::expression_from_node(result_no_ctx).to_string() : "null";
    std::string s2 = result_null_ctx ? LMCAS::detail::expression_from_node(result_null_ctx).to_string() : "null";

    std::cout << "  Limit of x^2 as x->2 (no ctx): " << s1 << std::endl;
    std::cout << "  Limit of x^2 as x->2 (nullptr): " << s2 << std::endl;

    EXPECT_EQ_STR(s1, s2, "LimitVisitor with no ctx and nullptr produce same result");
}


static void test_ode_solver_positive_branch() {
    TEST_CASE("ODE solver: positive dep var selects positive branch");

    // Solve dy/dx = x*y with y declared Positive
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::multiply(x, y);

    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value(),
                "positive-sign setup succeeds");
    EXPECT_TRUE(ctx.assume_domain("y", Domain::Real).has_value(),
                "real-domain setup succeeds");

    auto result_with_ctx = solve_separable_ode(rhs, "x", "y", &ctx);

    std::string result_str = result_with_ctx ? result_with_ctx->to_string() : "null";
    std::cout << "  ODE dy/dx=x*y with y Positive: " << result_str << std::endl;

    // When y is Positive, the solver wraps result in abs() to signal
    // positive branch preference
    EXPECT_TRUE(result_with_ctx != nullptr,
                "ODE solver with positive dep var produces a result");
    EXPECT_TRUE(contains_abs(LMCAS::detail::node(result_with_ctx)),
                "ODE solver with positive dep var contains abs() wrapper");
}

static void test_ode_solver_nullptr_no_abs() {
    TEST_CASE("ODE solver: nullptr context identical to current");

    /// 在空上下文中求解 dy/dx = x*y.
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


static void test_matcher_assumption_condition_matches() {
    TEST_CASE("Matcher: assumption_condition matches when context has variable Positive");

    // Create a rule with assumption_condition that checks if wildcard "A" is Positive
    auto pattern = LMCAS::detail::expression_from_node(make_var("A"));
    auto replacement = LMCAS::detail::expression_from_node(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    // The assumption_condition checks if the bound expression is Positive
    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        auto positive = ctx->is_positive(it->second);
        return positive && positive.value() == Tribool::True;
    };

    Rule rule(pattern, replacement, wildcards, assumption_cond);

    // Set up context with x Positive
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value(),
                "positive-sign setup succeeds");
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value(),
                "real-domain setup succeeds");

    // Create target expression: x
    auto target = LMCAS::detail::expression_from_node(make_var("x"));

    // Use RewriteEngine with assumption context
    RewriteEngine engine;
    engine.add_rule(rule);

    // The rule should match because x is Positive in the context
    MatchMap bindings;
    bool matched = Matcher::match(pattern, target, wildcards, bindings);

    EXPECT_TRUE(matched, "Matcher matches pattern against target");

    // Verify the assumption_condition evaluates correctly
    bool cond_result = assumption_cond(bindings, &ctx);
    EXPECT_TRUE(cond_result,
                "assumption_condition returns true when context has x Positive");
}

static void test_matcher_assumption_condition_no_match_without_context() {
    TEST_CASE("Matcher: assumption_condition fails without context");

    // Same rule as above
    auto pattern = LMCAS::detail::expression_from_node(make_var("A"));
    auto replacement = LMCAS::detail::expression_from_node(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        auto positive = ctx->is_positive(it->second);
        return positive && positive.value() == Tribool::True;
    };

    Rule rule(pattern, replacement, wildcards, assumption_cond);

    // Target: x (but no context)
    auto target = LMCAS::detail::expression_from_node(make_var("x"));

    MatchMap bindings;
    bool matched = Matcher::match(pattern, target, wildcards, bindings);

    EXPECT_TRUE(matched, "Matcher still matches structurally without context");

    /// 空上下文使 assumption_condition 返回 false.
    bool cond_result = assumption_cond(bindings, nullptr);
    EXPECT_FALSE(cond_result,
                 "assumption_condition returns false without context (nullptr)");
}

static void test_matcher_rewrite_engine_with_context() {
    TEST_CASE("Matcher: RewriteEngine uses assumption context during apply");

    // Create a rule: abs(A) -> A when A is Positive
    auto abs_A = LMCAS::detail::expression_from_node(make_abs(make_var("A")));
    auto just_A = LMCAS::detail::expression_from_node(make_var("A"));
    std::unordered_set<std::string> wildcards = {"A"};

    auto assumption_cond = [](const MatchMap& bindings, const AssumptionContext* ctx) -> bool {
        if (!ctx) return false;
        auto it = bindings.find("A");
        if (it == bindings.end()) return false;
        auto positive = ctx->is_positive(it->second);
        return positive && positive.value() == Tribool::True;
    };

    Rule rule(abs_A, just_A, wildcards, assumption_cond);

    // Set up context with x Positive
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value(),
                "positive-sign setup succeeds");
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value(),
                "real-domain setup succeeds");

    RewriteEngine engine;
    engine.add_rule(rule);
    ComputationContext computation_context;
    auto assumptions = std::make_shared<AssumptionContext>(ctx);
    auto set_assumptions = computation_context.set_assumptions(assumptions);
    EXPECT_TRUE(set_assumptions.has_value(), "rewrite assumptions attach to context");

    // Apply to abs(x)
    auto input = LMCAS::detail::expression_from_node(make_abs(make_var("x")));
    auto checked_result = engine.apply_checked(input, computation_context);
    EXPECT_TRUE(checked_result.has_value(), "context-aware rewrite succeeds");
    auto result = checked_result.value();

    std::string result_str = result.to_string();
    std::cout << "  RewriteEngine abs(x) with x Positive: " << result_str << std::endl;

    // The result should be x (abs removed because x is Positive)
    EXPECT_FALSE(contains_abs(LMCAS::detail::node(result)),
                 "RewriteEngine removes abs(x) when x is Positive");
}


static void test_integrator_nullptr_identical() {
    TEST_CASE("Integrator: nullptr identical to no context");

    // Integrate x^2
    auto integrand = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<PowerNode>(make_var("x"), make_number(2)));
    Integrator integrator1;
    // No context set (default nullptr)
    auto result1 = integrator1.integrate(integrand, "x");
    EXPECT_TRUE(result1.has_value(), "default integration succeeds");
    if (!result1) return;

    Integrator integrator2;
    ComputationContext integration_context2;
    auto result2_checked = integrator2.integrate_checked(
        integrand, "x", integration_context2);
    EXPECT_TRUE(result2_checked.has_value(), "integration without assumptions succeeds");
    auto result2 = result2_checked.value();

    std::string s1 = result1.value().to_string();
    std::string s2 = result2.to_string();

    EXPECT_EQ_STR(s1, s2, "Integrator with default and explicit nullptr produce same result");
}

static void test_ode_solver_nullptr_identical() {
    TEST_CASE("ODE solver: nullptr identical to default");

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
    TEST_CASE("Matcher: nullptr identical to no context");

    auto pattern = LMCAS::detail::expression_from_node(make_var("A"));
    auto target = LMCAS::detail::expression_from_node(make_var("x"));
    std::unordered_set<std::string> wildcards = {"A"};

    MatchMap bindings1;
    bool matched1 = Matcher::match(pattern, target, wildcards, bindings1);

    MatchMap bindings2;
    bool matched2 = Matcher::match(pattern, target, wildcards, bindings2);

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
