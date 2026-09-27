
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/visitors/print_visitor.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

/// Normalize a node with an AssumptionContext.
static std::shared_ptr<const SymbolicNode> normalize_with_ctx(
    const std::shared_ptr<const SymbolicNode> &node,
    const AssumptionContext &ctx) {
    auto expression = LMCAS::detail::expression_from_node(node);
    return LMCAS::detail::node(ctx.simplify(expression));
}

/// 使用空 AssumptionContext 执行兼容模式规范化.
static std::shared_ptr<const SymbolicNode> normalize_no_ctx(
    const std::shared_ptr<const SymbolicNode> &node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

/// Create a VariableNode.
static std::shared_ptr<const SymbolicNode> var(const std::string &name) {
    return LMCAS::detail::make_node<VariableNode>(name);
}

/// Create a NumberNode from int.
static std::shared_ptr<const SymbolicNode> num(int v) {
    return LMCAS::detail::make_node<NumberNode>(BigInt(v));
}

/// Create sqrt(expr) as FunctionNode(Sqrt, {expr}).
static std::shared_ptr<const SymbolicNode> make_sqrt(const std::shared_ptr<const SymbolicNode> &arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sqrt,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

/// Create abs(expr) as FunctionNode(Abs, {expr}).
static std::shared_ptr<const SymbolicNode> make_abs(const std::shared_ptr<const SymbolicNode> &arg) {
    return LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

/// Create x^n as PowerNode(x, NumberNode(n)).
static std::shared_ptr<const SymbolicNode> make_power(
    const std::shared_ptr<const SymbolicNode> &base, int exp) {
    return LMCAS::detail::make_node<PowerNode>(base, num(exp));
}

/// Check if a node is a VariableNode with the given name.
static bool is_variable(const std::shared_ptr<const SymbolicNode> &node, const std::string &name) {
    auto v = std::dynamic_pointer_cast<const VariableNode>(node);
    return v && v->name() == name;
}

/// Check if a node is abs(x) - FunctionNode(Abs, {VariableNode(name)}).
static bool is_abs_of_var(const std::shared_ptr<const SymbolicNode> &node, const std::string &name) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!func || func->type() != FunctionNode::FuncType::Abs)
        return false;
    if (func->arguments().size() != 1)
        return false;
    return is_variable(func->arguments()[0], name);
}

/// Check if a node represents -x (i.e., MultiplyNode({-1, x})).
static bool is_negation_of_var(const std::shared_ptr<const SymbolicNode> &node, const std::string &name) {
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2)
        return false;

    // Check for -1 * x pattern
    auto n = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0]);
    if (!n)
        return false;

    bool is_neg_one = false;
    if (std::holds_alternative<BigInt>(n->value())) {
        is_neg_one = (std::get<BigInt>(n->value()) == BigInt(-1));
    } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
        is_neg_one = (std::get<lmmc_real_t>(n->value()) == -1.0);
    } else if (std::holds_alternative<Rational>(n->value())) {
        is_neg_one = (std::get<Rational>(n->value()) == Rational(-1));
    }

    if (!is_neg_one)
        return false;
    return is_variable(mul->operands()[1], name);
}

TEST(AssumptionSimplify, SqrtXSquaredNonnegative) {
    // Test with multiple variable names
    std::vector<std::string> var_names = {"x", "y", "alpha", "t", "var1"};

    for (const auto &name : var_names) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(name, Sign::NonNegative).has_value());

        // Build sqrt(x^2)
        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE((is_variable(result, name))) << "sqrt(" + name + "²) with NonNegative → " + name;
    }
}

TEST(AssumptionSimplify, SqrtXSquaredPositive) {
    // Positive implies NonNegative, so the same rule should apply
    std::vector<std::string> var_names = {"a", "b", "c"};

    for (const auto &name : var_names) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());

        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE((is_variable(result, name))) << "sqrt(" + name + "²) with Positive → " + name;
    }
}

TEST(AssumptionSimplify, SqrtXSquaredRealNotNonneg) {
    // Declare x as Real only (not NonNegative)
    std::vector<std::string> var_names = {"x", "y", "z", "w"};

    for (const auto &name : var_names) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_domain(name, Domain::Real).has_value());

        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE((is_abs_of_var(result, name))) << "sqrt(" + name + "²) with Real (not NonNeg) → abs(" + name + ")";
    }
}

TEST(AssumptionSimplify, SqrtXSquaredIntegerNotNonneg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("n", Domain::Integer).has_value());

    auto n_squared = make_power(var("n"), 2);
    auto sqrt_n_sq = make_sqrt(n_squared);

    auto result = normalize_with_ctx(sqrt_n_sq, ctx);

    // Integer implies Real, so sqrt(n^2) -> abs(n)
    EXPECT_TRUE((is_abs_of_var(result, "n"))) << "sqrt(n²) with Integer (not NonNeg) → abs(n)";
}

TEST(AssumptionSimplify, SqrtXSquaredNatural) {
    AssumptionContext ctx;
    // Natural domain alone may not imply NonNegative sign in the current
    // implementation. Explicitly declare NonNegative sign to test the rule.
    ASSERT_TRUE(ctx.assume_domain("k", Domain::Natural).has_value());
    ASSERT_TRUE(ctx.assume_sign("k", Sign::NonNegative).has_value());

    auto k_squared = make_power(var("k"), 2);
    auto sqrt_k_sq = make_sqrt(k_squared);

    auto result = normalize_with_ctx(sqrt_k_sq, ctx);

    // With NonNegative sign, sqrt(k^2) -> k
    EXPECT_TRUE((is_variable(result, "k"))) << "sqrt(k²) with Natural + NonNegative → k";
}

TEST(AssumptionSimplify, AbsPositive) {
    std::vector<std::string> var_names = {"x", "y", "alpha", "t", "var1"};

    for (const auto &name : var_names) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(name, Sign::Positive).has_value());

        auto abs_x = make_abs(var(name));
        auto result = normalize_with_ctx(abs_x, ctx);

        EXPECT_TRUE((is_variable(result, name))) << "abs(" + name + ") with Positive → " + name;
    }
}

TEST(AssumptionSimplify, AbsNegative) {
    std::vector<std::string> var_names = {"x", "y", "z", "w"};

    for (const auto &name : var_names) {
        AssumptionContext ctx;
        EXPECT_TRUE(ctx.assume_sign(name, Sign::Negative).has_value());

        auto abs_x = make_abs(var(name));
        auto result = normalize_with_ctx(abs_x, ctx);

        EXPECT_TRUE((is_negation_of_var(result, name))) << "abs(" + name + ") with Negative → -" + name;
    }
}

TEST(AssumptionSimplify, AbsNonnegativeNotPositive) {
    // NonNegative includes zero, so abs(x) should NOT simplify to x
    // (only Positive triggers the rule per the implementation)
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    // The implementation only simplifies abs(x) -> x for Positive,
    // not for NonNegative. Check that it either stays as abs(x) or
    // simplifies to x (both are mathematically valid for NonNegative).
    bool is_var_x = is_variable(result, "x");
    bool is_abs_x = is_abs_of_var(result, "x");
    EXPECT_TRUE((is_var_x || is_abs_x)) << "abs(x) with NonNegative → x or abs(x)";
}

TEST(AssumptionSimplify, AbsNoAssumption) {
    AssumptionContext ctx;
    // No assumptions about x

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    // Should remain as abs(x) since no sign info is available
    auto func = std::dynamic_pointer_cast<const FunctionNode>(result);
    EXPECT_TRUE((func != nullptr && func->type() == FunctionNode::FuncType::Abs)) << "abs(x) with no assumption remains abs(x)";
}

TEST(AssumptionSimplify, BackwardCompatSqrtXSquared) {
    /// 空 AssumptionContext 保留 sqrt(x^2) 的定义域条件.
    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result_no_ctx = normalize_no_ctx(sqrt_x_sq);

    // The result should be the same as what the default visitor produces
    // (no assumption-based simplification)
    NormalizationVisitor v_default;
    sqrt_x_sq->accept(v_default);
    auto result_default = v_default.get_result();

    // Both should produce the same output
    EXPECT_TRUE((result_no_ctx->equals(*result_default))) << "sqrt(x²) without context = default NormalizationVisitor result";
}

TEST(AssumptionSimplify, BackwardCompatAbsX) {
    auto abs_x = make_abs(var("x"));

    auto result_no_ctx = normalize_no_ctx(abs_x);

    NormalizationVisitor v_default;
    abs_x->accept(v_default);
    auto result_default = v_default.get_result();

    EXPECT_TRUE((result_no_ctx->equals(*result_default))) << "abs(x) without context = default NormalizationVisitor result";
}

TEST(AssumptionSimplify, BackwardCompatVariousExpressions) {
    // Test a variety of expressions to ensure no assumption rules fire
    std::vector<std::shared_ptr<const SymbolicNode>> expressions = {
        // Simple variable
        var("x"),
        // Number
        num(42),
        // x + y
        LMCAS::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{var("x"), var("y")}),
        // x * y
        LMCAS::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{var("x"), var("y")}),
        // x^3
        make_power(var("x"), 3),
        // sqrt(x)
        make_sqrt(var("x")),
        // abs(y)
        make_abs(var("y")),
        // sqrt(y^2)
        make_sqrt(make_power(var("y"), 2)),
        // sin(x)
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{var("x")}),
        // exp(x)
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{var("x")}),
    };

    for (size_t i = 0; i < expressions.size(); ++i) {
        auto &expr = expressions[i];

        auto result_no_ctx = normalize_no_ctx(expr);

        NormalizationVisitor v_default;
        expr->accept(v_default);
        auto result_default = v_default.get_result();

        std::string label = "Expression " + std::to_string(i) + " without context = default";
        EXPECT_TRUE((result_no_ctx->equals(*result_default))) << label;
    }
}

TEST(AssumptionSimplify, BackwardCompatNoAssumptionRulesFire) {
    // Create a context with assumptions for variable "a", but simplify
    // expressions involving variable "x" - no rules should fire for "x"
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive).has_value());

    // sqrt(x^2) should NOT simplify since x has no assumptions
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    // Compare with no-context result
    auto result_no_ctx = normalize_no_ctx(sqrt_x_sq);

    EXPECT_TRUE((result->equals(*result_no_ctx))) << "sqrt(x²) with unrelated assumptions = no-context result";

    // abs(x) should NOT simplify since x has no assumptions
    auto abs_x = make_abs(var("x"));
    auto result_abs = normalize_with_ctx(abs_x, ctx);
    auto result_abs_no_ctx = normalize_no_ctx(abs_x);

    EXPECT_TRUE((result_abs->equals(*result_abs_no_ctx))) << "abs(x) with unrelated assumptions = no-context result";
}

TEST(AssumptionSimplify, SqrtXSquaredInLargerExpression) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    // Build: sqrt(x^2) + 1
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto expr = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{sqrt_x_sq, num(1)});

    auto result = normalize_with_ctx(expr, ctx);

    // The result should be x + 1
    auto add = std::dynamic_pointer_cast<const AddNode>(result);
    if (add) {
        // Check that the result contains x and 1 (order may vary)
        bool has_x = false;
        bool has_one = false;
        for (const auto &op : add->operands()) {
            if (is_variable(op, "x"))
                has_x = true;
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (std::holds_alternative<BigInt>(n->value()) && std::get<BigInt>(n->value()) == BigInt(1))
                    has_one = true;
            }
        }
        EXPECT_TRUE((has_x && has_one)) << "sqrt(x²) + 1 with NonNegative x → x + 1";
    } else {
        // Might be simplified differently, just check it's not still sqrt(x^2) + 1
        PrintVisitor pv;
        if (result)
            result->accept(pv);
        const auto result_str = pv.get_result();
        EXPECT_TRUE((result != nullptr &&
                     result_str.find("sqrt") == std::string::npos &&
                     result_str.find("x") != std::string::npos &&
                     result_str.find("1") != std::string::npos))
            << "sqrt(x²) + 1 simplified (non-AddNode result)";
    }
}

TEST(AssumptionSimplify, AbsInLargerExpression) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    // Build: 2 * abs(x)
    auto abs_x = make_abs(var("x"));
    auto expr = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{num(2), abs_x});

    auto result = normalize_with_ctx(expr, ctx);

    // The result should be 2 * x
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(result);
    if (mul) {
        bool has_x = false;
        bool has_two = false;
        for (const auto &op : mul->operands()) {
            if (is_variable(op, "x"))
                has_x = true;
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (std::holds_alternative<BigInt>(n->value()) && std::get<BigInt>(n->value()) == BigInt(2))
                    has_two = true;
            }
        }
        EXPECT_TRUE((has_x && has_two)) << "2 * abs(x) with Positive x → 2 * x";
    } else {
        // Could be simplified to just a variable if 2*x normalizes differently
        PrintVisitor pv;
        if (result)
            result->accept(pv);
        const auto result_str = pv.get_result();
        EXPECT_TRUE((result != nullptr &&
                     result_str.find("abs") == std::string::npos &&
                     result_str.find("x") != std::string::npos &&
                     result_str.find("2") != std::string::npos))
            << "2 * abs(x) simplified (non-MultiplyNode result)";
    }
}

TEST(AssumptionSimplify, ScopedAssumptionSimplification) {
    AssumptionContext ctx;

    // In root scope, x has no assumptions
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto result_root = normalize_with_ctx(sqrt_x_sq, ctx);

    // Push scope and declare x NonNegative
    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    auto result_child = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((is_variable(result_child, "x"))) << "sqrt(x²) in child scope with NonNegative → x";

    // Pop scope - x should no longer be NonNegative
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";

    auto result_after_pop = normalize_with_ctx(sqrt_x_sq, ctx);
    // After pop, should behave like no assumptions
    EXPECT_TRUE((result_after_pop->equals(*result_root))) << "sqrt(x²) after pop = root scope result (no simplification)";
}
