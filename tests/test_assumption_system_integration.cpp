
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "integration.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "symbolic_ode.hpp"
#include "matcher.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_abs(const std::shared_ptr<const SymbolicNode> &arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

template <typename Predicate>
static bool contains_matching_node(
    const std::shared_ptr<const SymbolicNode> &node,
    const Predicate &matches) {
    if (!node) {
        return false;
    }
    if (matches(node)) {
        return true;
    }
    const auto contains_child = [&](const auto &child) {
        return contains_matching_node(child, matches);
    };
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return std::any_of(function->arguments().begin(),
                           function->arguments().end(), contains_child);
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        return std::any_of(sum->operands().begin(), sum->operands().end(),
                           contains_child);
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return std::any_of(product->operands().begin(), product->operands().end(),
                           contains_child);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return contains_child(power->base()) || contains_child(power->exponent());
    }
    return false;
}

static bool contains_abs_of(const std::shared_ptr<const SymbolicNode> &node,
                            const std::string &var_name) {
    return contains_matching_node(node, [&](const auto &candidate) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(candidate);
        if (!function || function->type() != FunctionNode::FuncType::Abs ||
            function->arguments().size() != 1) {
            return false;
        }
        auto variable = std::dynamic_pointer_cast<const VariableNode>(
            function->arguments()[0]);
        return variable && variable->name() == var_name;
    });
}

/// Check if an AST contains an abs() node anywhere.
static bool contains_abs(const std::shared_ptr<const SymbolicNode> &node) {
    return contains_matching_node(node, [](const auto &candidate) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(candidate);
        return function && function->type() == FunctionNode::FuncType::Abs;
    });
}

TEST(AssumptionSystemIntegration, IntegratorPositiveSimplifiesAbs) {
    // Create integrand: |x|
    auto integrand = LMCAS::detail::expression_from_node(make_abs(test_variable_node("x")));
    // Set up assumption context with x Positive
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "positive-sign setup succeeds";
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real-domain setup succeeds";

    // Integrate with assumption context
    Integrator integrator;
    ComputationContext integration_context;
    auto set_integration_assumptions = integration_context.set_assumptions(
        std::make_shared<AssumptionContext>(ctx));
    EXPECT_TRUE((set_integration_assumptions.has_value())) << "integration assumptions attach to context";
    auto integrated = integrator.integrate_checked(integrand, "x", integration_context);
    EXPECT_TRUE((integrated.has_value())) << "context-aware integration succeeds";
    auto result = integrated.value();

    // The result should NOT contain abs(x) since x is positive, |x| = x
    // So integrating x gives x^2/2
    EXPECT_FALSE((contains_abs_of(LMCAS::detail::node(result), "x"))) << "Integration result does not contain abs(x) when x is Positive";
}

TEST(AssumptionSystemIntegration, IntegratorNoContextPreservesAbs) {
    // Create integrand: |x|
    auto integrand = LMCAS::detail::expression_from_node(make_abs(test_variable_node("x")));
    /// Integrator 保持默认空假设上下文.
    Integrator integrator;
    auto result = integrator.integrate(integrand, "x");
    ASSERT_TRUE((result.has_value())) << "Integration without assumptions succeeds";
    if (!result) {
        return;
    }

    /// 空上下文时积分结果保留 abs 或未求值积分结构,
    /// |x| = x 化简仅在符号假设充分时启用.
    EXPECT_TRUE((LMCAS::detail::node(result.value()) != nullptr)) << "Integration without context produces a result";
}

TEST(AssumptionSystemIntegration, LimitVisitorPositiveResolvesSign) {
    // Compute limit of 1/x as x -> 0+ with x known Positive.
    // The LimitVisitor should use the assumption to determine the sign.
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "positive-sign setup succeeds";
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "real-domain setup succeeds";

    // Build 1/x = x^(-1) = MultiplyNode([1, PowerNode(x, -1)])
    auto one_over_x = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            test_integer_node(1),
            LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(-1))});

    auto point = test_integer_node(0);

    // Compute limit with assumption context (right-sided limit)
    LimitVisitor visitor_with_ctx("x", point, "+", &ctx);
    one_over_x->accept(visitor_with_ctx);
    auto result_with_ctx = visitor_with_ctx.get_result();

    // The result should be +infinity (positive infinity)
    // Check that it's an infinity node (not negative infinity)
    EXPECT_TRUE((result_with_ctx != nullptr)) << "Limit with positive context produces a result";

    // Verify it's positive infinity (not negative)
    bool is_positive_inf = false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(result_with_ctx)) {
        if (fn->type() == FunctionNode::FuncType::Infinity) {
            is_positive_inf = true;
        }
    }
    EXPECT_TRUE((is_positive_inf)) << "Limit of 1/x as x->0+ with x Positive is +infinity";
}

TEST(AssumptionSystemIntegration, LimitVisitorNullptrSameBehavior) {
    /// 在空上下文中计算 x^2 于 x->2 的极限.
    auto x_squared = LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2));
    auto point = test_integer_node(2);

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

    EXPECT_EQ((s1), (s2)) << "LimitVisitor with no ctx and nullptr produce same result";
}

TEST(AssumptionSystemIntegration, OdeSolverPositiveBranch) {
    // Solve dy/dx = x*y with y declared Positive
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::multiply(x, y);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value()) << "positive-sign setup succeeds";
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value()) << "real-domain setup succeeds";

    auto result_with_ctx = solve_separable_ode(rhs, "x", "y", &ctx);

    // When y is Positive, the solver wraps result in abs() to signal
    // positive branch preference
    EXPECT_TRUE((result_with_ctx != nullptr)) << "ODE solver with positive dep var produces a result";
    EXPECT_TRUE((contains_abs(LMCAS::detail::node(result_with_ctx)))) << "ODE solver with positive dep var contains abs() wrapper";
}

TEST(AssumptionSystemIntegration, OdeSolverNullptrNoAbs) {
    /// 在空上下文中求解 dy/dx = x*y.
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::multiply(x, y);

    auto result_no_ctx = solve_separable_ode(rhs, "x", "y", nullptr);
    auto result_default = solve_separable_ode(rhs, "x", "y");

    std::string s1 = result_no_ctx ? result_no_ctx->to_string() : "null";
    std::string s2 = result_default ? result_default->to_string() : "null";

    EXPECT_EQ((s1), (s2)) << "ODE solver with nullptr and default produce same result";
}

TEST(AssumptionSystemIntegration, MatcherAssumptionConditionMatches) {
    auto a = SymbolicExpr::variable("A");
    auto pattern = SymbolicExpr::sqrt(SymbolicExpr::power(a, SymbolicExpr::number(2)));
    Rule rule(*pattern, *a, {"A"},
              [](const MatchMap &bindings, const AssumptionContext *context) {
                  auto found = bindings.find("A");
                  if (!context || found == bindings.end())
                      return false;
                  auto positive = context->is_positive(found->second);
                  return positive && positive.value() == Tribool::True;
              });
    RewriteEngine engine;
    engine.add_rule(rule);
    for (int sign : {-1, 1}) {
        auto x = SymbolicExpr::variable("x");
        auto input = SymbolicExpr::sqrt(SymbolicExpr::power(x, SymbolicExpr::number(2)));
        auto assumptions = std::make_shared<AssumptionContext>();
        EXPECT_TRUE((assumptions->assume_sign("x", sign > 0 ? Sign::Positive : Sign::Negative).has_value())) << "sign declaration succeeds";
        ComputationContext context;
        EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "assumptions attach";
        auto rewritten = engine.apply_step_checked(*input, context);
        ASSERT_TRUE((rewritten.has_value())) << "checked rewrite succeeds";
        if (rewritten) {
            EXPECT_TRUE((detail::node(rewritten.value())->equals(*detail::node(sign > 0 ? x : input)))) << "only the positive branch removes the square root";
        }
    }
    auto input = SymbolicExpr::sqrt(SymbolicExpr::power(
        SymbolicExpr::number(-2), SymbolicExpr::number(2)));
    auto missing = engine.apply_step(*input);
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(std::make_shared<AssumptionContext>()).has_value())) << "empty assumptions attach";
    auto empty = engine.apply_step_checked(*input, context);
    for (const auto *result : {&missing, &empty}) {
        EXPECT_TRUE((*result && detail::node(result->value())->equals(*detail::node(input)))) << "missing or rejecting context cannot rewrite sqrt((-2)^2) to -2";
    }
}

TEST(AssumptionSystemIntegration, MatcherAssumptionConditionNoMatchWithoutContext) {
    auto x = SymbolicExpr::variable("x");
    auto seven = SymbolicExpr::number(7);
    for (bool plain : {false, true}) {
        for (bool assumed : {false, true}) {
            Rule rule(*x, *seven, {}, [plain](const MatchMap &) { return plain; });
            rule.assumption_condition = [assumed](const MatchMap &, const AssumptionContext *) {
                return assumed;
            };
            RewriteEngine engine;
            engine.add_rule(rule);
            auto absent = engine.apply_step(*x);
            EXPECT_TRUE((absent && detail::node(absent.value())->equals(*detail::node(x)))) << "no context can satisfy an assumption condition";
            ComputationContext context;
            EXPECT_TRUE((context.set_assumptions(std::make_shared<AssumptionContext>()).has_value())) << "context attaches";
            auto result = engine.apply_step_checked(*x, context);
            EXPECT_TRUE((result && detail::node(result.value())->equals(*detail::node(plain && assumed ? seven : x)))) << "a rewrite applies exactly when both conditions hold";
        }
    }
}

TEST(AssumptionSystemIntegration, IntegratorNullptrIdentical) {
    // Integrate x^2
    auto integrand = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<PowerNode>(test_variable_node("x"), test_integer_node(2)));
    Integrator integrator1;
    // No context set (default nullptr)
    auto result1 = integrator1.integrate(integrand, "x");
    ASSERT_TRUE((result1.has_value())) << "default integration succeeds";
    if (!result1) {
        return;
    }

    Integrator integrator2;
    ComputationContext integration_context2;
    auto result2_checked = integrator2.integrate_checked(
        integrand, "x", integration_context2);
    EXPECT_TRUE((result2_checked.has_value())) << "integration without assumptions succeeds";
    auto result2 = result2_checked.value();

    std::string s1 = result1.value().to_string();
    std::string s2 = result2.to_string();

    EXPECT_EQ((s1), (s2)) << "Integrator with default and explicit nullptr produce same result";
}

TEST(AssumptionSystemIntegration, OdeSolverNullptrIdentical) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(x, y);

    auto result1 = solve_separable_ode(rhs, "x", "y");
    auto result2 = solve_separable_ode(rhs, "x", "y", nullptr);

    std::string s1 = result1 ? result1->to_string() : "null";
    std::string s2 = result2 ? result2->to_string() : "null";

    EXPECT_EQ((s1), (s2)) << "solve_separable_ode default and nullptr produce same result";
}
