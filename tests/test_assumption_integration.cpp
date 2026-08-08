
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "integration.hpp"
#include "visitors/differentiation_visitor.hpp"
#include "visitors/limit_visitor.hpp"
#include "symbolic_ode.hpp"
#include "matcher.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "visitors/print_visitor.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

/// Convert a SymbolicNode to string via PrintVisitor
static std::string node_to_string(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return "null";
    PrintVisitor pv;
    node->accept(pv);
    return pv.get_result();
}


static void test_integrator_nullptr_backward_compat() {
    TEST_CASE("Integrator with nullptr context produces same results as before");

    rc::check("Integrator with nullptr context is identical to default Integrator for polynomial integrands", []() {
        // Generate a random polynomial integrand: a*x^n where a in [1,5], n in [0,4]
        int a = rc::gen::inRange(1, 6);
        int n = rc::gen::inRange(0, 5);

        // Build expression: a * x^n
        auto x_node = make_var("x");
        auto coeff_node = make_number(a);
        auto exp_node = make_number(n);
        auto power_node = lamina::detail::make_node<PowerNode>(x_node, exp_node);
        std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {coeff_node, power_node};
        auto mul_node = lamina::detail::make_node<MultiplyNode>(mul_ops);
        auto integrand = lamina::detail::expression_from_node(mul_node);
        // Integrator without context (default)
        Integrator integrator_default;
        auto result_default = integrator_default.integrate(integrand, "x");

        // Integrator with explicit nullptr context
        Integrator integrator_nullptr;
        integrator_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = integrator_nullptr.integrate(integrand, "x");

        // Results must be identical
        RC_ASSERT(result_default.to_string() == result_nullptr.to_string());
    });
}

static void test_integrator_trig_nullptr() {
    TEST_CASE("Integrator nullptr for trig functions");

    // Test sin(x), cos(x), exp(x) — known integrals
    {
        auto x = SymbolicExpr::variable("x");
        auto sin_x = SymbolicExpr::sin(x);

        Integrator integrator_default;
        auto result_default = integrator_default.integrate(*sin_x, "x");

        Integrator integrator_nullptr;
        integrator_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = integrator_nullptr.integrate(*sin_x, "x");

        EXPECT_EQ_STR(result_default.to_string(), result_nullptr.to_string(),
                      "sin(x) integral: default == nullptr context");
    }
    {
        auto x = SymbolicExpr::variable("x");
        auto cos_x = SymbolicExpr::cos(x);

        Integrator integrator_default;
        auto result_default = integrator_default.integrate(*cos_x, "x");

        Integrator integrator_nullptr;
        integrator_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = integrator_nullptr.integrate(*cos_x, "x");

        EXPECT_EQ_STR(result_default.to_string(), result_nullptr.to_string(),
                      "cos(x) integral: default == nullptr context");
    }
    {
        auto x = SymbolicExpr::variable("x");
        auto exp_x = SymbolicExpr::exp(x);

        Integrator integrator_default;
        auto result_default = integrator_default.integrate(*exp_x, "x");

        Integrator integrator_nullptr;
        integrator_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = integrator_nullptr.integrate(*exp_x, "x");

        EXPECT_EQ_STR(result_default.to_string(), result_nullptr.to_string(),
                      "exp(x) integral: default == nullptr context");
    }
}


static void test_limit_nullptr_backward_compat() {
    TEST_CASE("LimitVisitor with nullptr context produces same results as before");

    rc::check("LimitVisitor with nullptr context is identical to default LimitVisitor", []() {
        // Generate a random polynomial: a*x + b where a in [1,5], b in [-3,3]
        int a = rc::gen::inRange(1, 6);
        int b = rc::gen::inRange(-3, 4);

        // Build expression: a*x + b
        auto x_node = make_var("x");
        auto a_node = make_number(a);
        auto b_node = make_number(b);
        std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {a_node, x_node};
        auto ax = lamina::detail::make_node<MultiplyNode>(mul_ops);
        std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {ax, b_node};
        auto expr_node = lamina::detail::make_node<AddNode>(add_ops);

        // Limit point: x → 1
        auto point = lamina::detail::make_node<NumberNode>(BigInt(1));

        // LimitVisitor without context (default)
        LimitVisitor visitor_default("x", point, "");
        expr_node->accept(visitor_default);
        auto result_default = visitor_default.get_result();

        // LimitVisitor with explicit nullptr context
        LimitVisitor visitor_nullptr("x", point, "", nullptr);
        expr_node->accept(visitor_nullptr);
        auto result_nullptr = visitor_nullptr.get_result();

        std::string str_default = node_to_string(result_default);
        std::string str_nullptr = node_to_string(result_nullptr);

        RC_ASSERT(str_default == str_nullptr);
    });
}

static void test_limit_trig_nullptr() {
    TEST_CASE("LimitVisitor nullptr for trig at 0");

    // lim x→0 sin(x) = 0
    {
        auto x_node = make_var("x");
        std::vector<std::shared_ptr<const SymbolicNode>> args = {x_node};
        auto sin_x = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, args);
        auto point = lamina::detail::make_node<NumberNode>(BigInt(0));

        LimitVisitor visitor_default("x", point, "");
        sin_x->accept(visitor_default);
        auto result_default = visitor_default.get_result();

        LimitVisitor visitor_nullptr("x", point, "", nullptr);
        sin_x->accept(visitor_nullptr);
        auto result_nullptr = visitor_nullptr.get_result();

        std::string str_default = node_to_string(result_default);
        std::string str_nullptr = node_to_string(result_nullptr);

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "lim x->0 sin(x): default == nullptr context");
    }
    // lim x→0 cos(x) = 1
    {
        auto x_node = make_var("x");
        std::vector<std::shared_ptr<const SymbolicNode>> args = {x_node};
        auto cos_x = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, args);
        auto point = lamina::detail::make_node<NumberNode>(BigInt(0));

        LimitVisitor visitor_default("x", point, "");
        cos_x->accept(visitor_default);
        auto result_default = visitor_default.get_result();

        LimitVisitor visitor_nullptr("x", point, "", nullptr);
        cos_x->accept(visitor_nullptr);
        auto result_nullptr = visitor_nullptr.get_result();

        std::string str_default = node_to_string(result_default);
        std::string str_nullptr = node_to_string(result_nullptr);

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "lim x->0 cos(x): default == nullptr context");
    }
}


static void test_series_nullptr_backward_compat() {
    TEST_CASE("Series expansion with nullptr context produces same results as before");

    rc::check("Series expansion with nullptr context is identical to default for exp(x)", []() {
        int order = rc::gen::inRange(2, 6);

        auto x = SymbolicExpr::variable("x");
        auto exp_x = SymbolicExpr::exp(x);
        auto zero = SymbolicExpr::number(0);

        // Series without context (default)
        auto result_default = exp_x->series("x", zero, order);

        // Series with explicit nullptr context
        auto result_nullptr = exp_x->series("x", zero, order, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        RC_ASSERT(str_default == str_nullptr);
    });
}

static void test_series_sin_cos_nullptr() {
    TEST_CASE("Series nullptr for sin/cos");

    // sin(x) series at 0, order 5
    {
        auto x = SymbolicExpr::variable("x");
        auto sin_x = SymbolicExpr::sin(x);
        auto zero = SymbolicExpr::number(0);

        auto result_default = sin_x->series("x", zero, 5);
        auto result_nullptr = sin_x->series("x", zero, 5, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "sin(x) series order 5: default == nullptr context");
    }
    // cos(x) series at 0, order 4
    {
        auto x = SymbolicExpr::variable("x");
        auto cos_x = SymbolicExpr::cos(x);
        auto zero = SymbolicExpr::number(0);

        auto result_default = cos_x->series("x", zero, 4);
        auto result_nullptr = cos_x->series("x", zero, 4, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "cos(x) series order 4: default == nullptr context");
    }
}


static void test_ode_nullptr_backward_compat() {
    TEST_CASE("ODE solver with nullptr context produces same results as before");

    // Test solve_separable_ode: dy/dx = x/y
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(x, y);

        auto result_default = solve_separable_ode(rhs, "x", "y");
        auto result_nullptr = solve_separable_ode(rhs, "x", "y", nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "separable ODE dy/dx=x/y: default == nullptr context");
    }

    // Test solve_separable_ode: dy/dx = x*y
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::multiply(x, y);

        auto result_default = solve_separable_ode(rhs, "x", "y");
        auto result_nullptr = solve_separable_ode(rhs, "x", "y", nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "separable ODE dy/dx=x*y: default == nullptr context");
    }
}

static void test_ode_linear_nullptr() {
    TEST_CASE("ODE linear solvers with nullptr");

    // Test solve_linear1_ode: dy/dx + 2*y = 0
    {
        auto Px = SymbolicExpr::number(2);
        auto Qx = SymbolicExpr::number(0);

        auto result_default = solve_linear1_ode(Px, Qx, "x", "y");
        auto result_nullptr = solve_linear1_ode(Px, Qx, "x", "y", nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "linear1 ODE P=2 Q=0: default == nullptr context");
    }

    // Test solve_linear2_ode: y'' + y = 0 (a=1, b=0, c=1, f=0)
    {
        auto fx = SymbolicExpr::number(0);

        auto result_default = solve_linear2_ode(1.0, 0.0, 1.0, fx, "x", "y");
        auto result_nullptr = solve_linear2_ode(1.0, 0.0, 1.0, fx, "x", "y", nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ_STR(str_default, str_nullptr,
                      "linear2 ODE y''+y=0: default == nullptr context");
    }
}


static void test_matcher_nullptr_backward_compat() {
    TEST_CASE("Matcher with nullptr context produces same results as before");

    rc::check("Matcher::match with nullptr context is identical to default match", []() {
        // Generate a random pattern: wildcard + number
        int num_val = rc::gen::inRange(1, 10);

        // Pattern: _a + num_val
        auto wc = lamina::detail::make_node<VariableNode>("_a");
        auto num = make_number(num_val);
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {wc, num};
        auto pat_node = lamina::detail::make_node<AddNode>(pat_ops);
        auto pattern = lamina::detail::expression_from_node(pat_node);
        // Target: x + num_val
        auto x_node = make_var("x");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {x_node, num};
        auto tgt_node = lamina::detail::make_node<AddNode>(tgt_ops);
        auto target = lamina::detail::expression_from_node(tgt_node);
        std::unordered_set<std::string> wildcards = {"_a"};

        // Match without context (default)
        MatchMap results_default;
        bool matched_default = Matcher::match(pattern, target, wildcards, results_default);

        // Match with explicit nullptr context
        MatchMap results_nullptr;
        bool matched_nullptr = Matcher::match(pattern, target, wildcards, results_nullptr, nullptr);

        RC_ASSERT(matched_default == matched_nullptr);

        if (matched_default && matched_nullptr) {
            // Verify bindings are identical
            RC_ASSERT(results_default.size() == results_nullptr.size());
            for (auto& [key, val] : results_default) {
                RC_ASSERT(results_nullptr.count(key) > 0);
                RC_ASSERT(val.to_string() == results_nullptr.at(key).to_string());
            }
        }
    });
}

static void test_rewrite_engine_nullptr() {
    TEST_CASE("RewriteEngine with nullptr context");

    // Test that RewriteEngine with nullptr context behaves identically to default
    // Rule: x + 0 -> x
    {
        auto wc_a = wildcard("_a");
        auto zero_expr = lamina::detail::expression_from_node(make_number(0));
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {lamina::detail::node(wc_a), lamina::detail::node(zero_expr)};
        auto pattern = lamina::detail::expression_from_node(lamina::detail::make_node<AddNode>(pat_ops));
        SymbolicExpr replacement = wc_a;
        std::unordered_set<std::string> wcs = {"_a"};

        Rule rule(pattern, replacement, wcs);

        // Target: y + 0
        auto y_node = make_var("y");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {y_node, make_number(0)};
        auto target = lamina::detail::expression_from_node(lamina::detail::make_node<AddNode>(tgt_ops));
        // RewriteEngine without context (default)
        RewriteEngine engine_default;
        engine_default.add_rule(rule);
        auto result_default = engine_default.apply(target, 10);

        // RewriteEngine with explicit nullptr context
        RewriteEngine engine_nullptr;
        engine_nullptr.add_rule(rule);
        engine_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = engine_nullptr.apply(target, 10);

        EXPECT_EQ_STR(result_default.to_string(), result_nullptr.to_string(),
                      "RewriteEngine y+0: default == nullptr context");
    }

    // Rule: _a * 1 -> _a
    {
        auto wc_a = wildcard("_a");
        auto one_expr = lamina::detail::expression_from_node(make_number(1));
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {lamina::detail::node(wc_a), lamina::detail::node(one_expr)};
        auto pattern = lamina::detail::expression_from_node(lamina::detail::make_node<MultiplyNode>(pat_ops));
        SymbolicExpr replacement = wc_a;
        std::unordered_set<std::string> wcs = {"_a"};

        Rule rule(pattern, replacement, wcs);

        // Target: z * 1
        auto z_node = make_var("z");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {z_node, make_number(1)};
        auto target = lamina::detail::expression_from_node(lamina::detail::make_node<MultiplyNode>(tgt_ops));
        RewriteEngine engine_default;
        engine_default.add_rule(rule);
        auto result_default = engine_default.apply(target, 10);

        RewriteEngine engine_nullptr;
        engine_nullptr.add_rule(rule);
        engine_nullptr.set_assumption_context(nullptr);
        auto result_nullptr = engine_nullptr.apply(target, 10);

        EXPECT_EQ_STR(result_default.to_string(), result_nullptr.to_string(),
                      "RewriteEngine z*1: default == nullptr context");
    }
}


int main() {
    test_integrator_nullptr_backward_compat();
    test_integrator_trig_nullptr();

    test_limit_nullptr_backward_compat();
    test_limit_trig_nullptr();

    test_series_nullptr_backward_compat();
    test_series_sin_cos_nullptr();

    test_ode_nullptr_backward_compat();
    test_ode_linear_nullptr();

    test_matcher_nullptr_backward_compat();
    test_rewrite_engine_nullptr();

    return TEST_REPORT();
}
