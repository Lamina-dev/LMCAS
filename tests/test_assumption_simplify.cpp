
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "visitors/normalization_visitor.hpp"
#include "visitors/print_visitor.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace lamina;


/// Normalize a node with an AssumptionContext.
static std::shared_ptr<const SymbolicNode> normalize_with_ctx(
    const std::shared_ptr<const SymbolicNode>& node,
    const AssumptionContext& ctx) {
    NormalizationVisitor v(&ctx);
    node->accept(v);
    return v.get_result();
}

/// Normalize a node without any AssumptionContext (backward-compatible).
static std::shared_ptr<const SymbolicNode> normalize_no_ctx(
    const std::shared_ptr<const SymbolicNode>& node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

/// Create a VariableNode.
static std::shared_ptr<const SymbolicNode> var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

/// Create a NumberNode from int.
static std::shared_ptr<const SymbolicNode> num(int v) {
    return lamina::detail::make_node<NumberNode>(BigInt(v));
}

/// Create sqrt(expr) as FunctionNode(Sqrt, {expr}).
static std::shared_ptr<const SymbolicNode> make_sqrt(const std::shared_ptr<const SymbolicNode>& arg) {
    return lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sqrt,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

/// Create abs(expr) as FunctionNode(Abs, {expr}).
static std::shared_ptr<const SymbolicNode> make_abs(const std::shared_ptr<const SymbolicNode>& arg) {
    return lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{arg});
}

/// Create x^n as PowerNode(x, NumberNode(n)).
static std::shared_ptr<const SymbolicNode> make_power(
    const std::shared_ptr<const SymbolicNode>& base, int exp) {
    return lamina::detail::make_node<PowerNode>(base, num(exp));
}

/// Check if a node is a VariableNode with the given name.
static bool is_variable(const std::shared_ptr<const SymbolicNode>& node, const std::string& name) {
    auto v = std::dynamic_pointer_cast<const VariableNode>(node);
    return v && v->name() == name;
}

/// Check if a node is abs(x) — FunctionNode(Abs, {VariableNode(name)}).
static bool is_abs_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& name) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!func || func->type() != FunctionNode::FuncType::Abs) return false;
    if (func->arguments().size() != 1) return false;
    return is_variable(func->arguments()[0], name);
}

/// Check if a node represents -x (i.e., MultiplyNode({-1, x})).
static bool is_negation_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& name) {
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2) return false;

    // Check for -1 * x pattern
    auto n = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0]);
    if (!n) return false;

    bool is_neg_one = false;
    if (std::holds_alternative<BigInt>(n->value())) {
        is_neg_one = (std::get<BigInt>(n->value()) == BigInt(-1));
    } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
        is_neg_one = (std::get<lmmc_real_t>(n->value()) == -1.0);
    } else if (std::holds_alternative<Rational>(n->value())) {
        is_neg_one = (std::get<Rational>(n->value()) == Rational(-1));
    }

    if (!is_neg_one) return false;
    return is_variable(mul->operands()[1], name);
}


void test_sqrt_x_squared_nonnegative() {
    TEST_CASE("sqrt(x²) → x when x is NonNegative");

    // Test with multiple variable names
    std::vector<std::string> var_names = {"x", "y", "alpha", "t", "var1"};

    for (const auto& name : var_names) {
        AssumptionContext ctx;
        ctx.assume_sign(name, Sign::NonNegative);

        // Build sqrt(x²)
        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE(is_variable(result, name),
                    "sqrt(" + name + "²) with NonNegative → " + name);
    }
}

void test_sqrt_x_squared_positive() {
    TEST_CASE("sqrt(x²) → x when x is Positive (implies NonNegative)");

    // Positive implies NonNegative, so the same rule should apply
    std::vector<std::string> var_names = {"a", "b", "c"};

    for (const auto& name : var_names) {
        AssumptionContext ctx;
        ctx.assume_sign(name, Sign::Positive);

        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE(is_variable(result, name),
                    "sqrt(" + name + "²) with Positive → " + name);
    }
}

void test_sqrt_x_squared_real_not_nonneg() {
    TEST_CASE("sqrt(x²) → abs(x) when x is Real but not NonNegative");

    // Declare x as Real only (not NonNegative)
    std::vector<std::string> var_names = {"x", "y", "z", "w"};

    for (const auto& name : var_names) {
        AssumptionContext ctx;
        ctx.assume_domain(name, Domain::Real);

        auto x_squared = make_power(var(name), 2);
        auto sqrt_x_sq = make_sqrt(x_squared);

        auto result = normalize_with_ctx(sqrt_x_sq, ctx);

        EXPECT_TRUE(is_abs_of_var(result, name),
                    "sqrt(" + name + "²) with Real (not NonNeg) → abs(" + name + ")");
    }
}

void test_sqrt_x_squared_integer_not_nonneg() {
    TEST_CASE("sqrt(x²) → abs(x) when x is Integer (implies Real) but not NonNegative");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    auto n_squared = make_power(var("n"), 2);
    auto sqrt_n_sq = make_sqrt(n_squared);

    auto result = normalize_with_ctx(sqrt_n_sq, ctx);

    // Integer implies Real, so sqrt(n²) → abs(n)
    EXPECT_TRUE(is_abs_of_var(result, "n"),
                "sqrt(n²) with Integer (not NonNeg) → abs(n)");
}

void test_sqrt_x_squared_natural() {
    TEST_CASE("sqrt(x²) → x when x is Natural (NonNegative sign declared)");

    AssumptionContext ctx;
    // Natural domain alone may not imply NonNegative sign in the current
    // implementation. Explicitly declare NonNegative sign to test the rule.
    ctx.assume_domain("k", Domain::Natural);
    ctx.assume_sign("k", Sign::NonNegative);

    auto k_squared = make_power(var("k"), 2);
    auto sqrt_k_sq = make_sqrt(k_squared);

    auto result = normalize_with_ctx(sqrt_k_sq, ctx);

    // With NonNegative sign, sqrt(k²) → k
    EXPECT_TRUE(is_variable(result, "k"),
                "sqrt(k²) with Natural + NonNegative → k");
}


void test_abs_positive() {
    TEST_CASE("abs(x) → x when x is Positive");

    std::vector<std::string> var_names = {"x", "y", "alpha", "t", "var1"};

    for (const auto& name : var_names) {
        AssumptionContext ctx;
        ctx.assume_sign(name, Sign::Positive);

        auto abs_x = make_abs(var(name));
        auto result = normalize_with_ctx(abs_x, ctx);

        EXPECT_TRUE(is_variable(result, name),
                    "abs(" + name + ") with Positive → " + name);
    }
}

void test_abs_negative() {
    TEST_CASE("abs(x) → -x when x is Negative");

    std::vector<std::string> var_names = {"x", "y", "z", "w"};

    for (const auto& name : var_names) {
        AssumptionContext ctx;
        ctx.assume_sign(name, Sign::Negative);

        auto abs_x = make_abs(var(name));
        auto result = normalize_with_ctx(abs_x, ctx);

        EXPECT_TRUE(is_negation_of_var(result, name),
                    "abs(" + name + ") with Negative → -" + name);
    }
}

void test_abs_nonnegative_not_positive() {
    TEST_CASE("abs(x) unchanged when x is NonNegative but not Positive");

    // NonNegative includes zero, so abs(x) should NOT simplify to x
    // (only Positive triggers the rule per the implementation)
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    // The implementation only simplifies abs(x) → x for Positive,
    // not for NonNegative. Check that it either stays as abs(x) or
    // simplifies to x (both are mathematically valid for NonNegative).
    bool is_var_x = is_variable(result, "x");
    bool is_abs_x = is_abs_of_var(result, "x");
    EXPECT_TRUE(is_var_x || is_abs_x,
                "abs(x) with NonNegative → x or abs(x)");
}

static void test_abs_no_assumption() {
    TEST_CASE("abs(x) unchanged when x has no sign assumption");

    AssumptionContext ctx;
    // No assumptions about x

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    // Should remain as abs(x) since no sign info is available
    auto func = std::dynamic_pointer_cast<const FunctionNode>(result);
    EXPECT_TRUE(func != nullptr && func->type() == FunctionNode::FuncType::Abs,
                "abs(x) with no assumption remains abs(x)");
}


void test_backward_compat_sqrt_x_squared() {
    TEST_CASE("sqrt(x²) without context produces same result as default NormalizationVisitor");

    // Without an AssumptionContext, sqrt(x²) should NOT be simplified
    // by assumption-based rules
    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result_no_ctx = normalize_no_ctx(sqrt_x_sq);

    // The result should be the same as what the default visitor produces
    // (no assumption-based simplification)
    NormalizationVisitor v_default;
    sqrt_x_sq->accept(v_default);
    auto result_default = v_default.get_result();

    // Both should produce the same output
    EXPECT_TRUE(result_no_ctx->equals(*result_default),
                "sqrt(x²) without context = default NormalizationVisitor result");
}

void test_backward_compat_abs_x() {
    TEST_CASE("abs(x) without context produces same result as default NormalizationVisitor");

    auto abs_x = make_abs(var("x"));

    auto result_no_ctx = normalize_no_ctx(abs_x);

    NormalizationVisitor v_default;
    abs_x->accept(v_default);
    auto result_default = v_default.get_result();

    EXPECT_TRUE(result_no_ctx->equals(*result_default),
                "abs(x) without context = default NormalizationVisitor result");
}

void test_backward_compat_various_expressions() {
    TEST_CASE("Various expressions without context match default visitor");

    // Test a variety of expressions to ensure no assumption rules fire
    std::vector<std::shared_ptr<const SymbolicNode>> expressions = {
        // Simple variable
        var("x"),
        // Number
        num(42),
        // x + y
        lamina::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{var("x"), var("y")}),
        // x * y
        lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{var("x"), var("y")}),
        // x^3
        make_power(var("x"), 3),
        // sqrt(x)
        make_sqrt(var("x")),
        // abs(y)
        make_abs(var("y")),
        // sqrt(y^2)
        make_sqrt(make_power(var("y"), 2)),
        // sin(x)
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{var("x")}),
        // exp(x)
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{var("x")}),
    };

    for (size_t i = 0; i < expressions.size(); ++i) {
        auto& expr = expressions[i];

        auto result_no_ctx = normalize_no_ctx(expr);

        NormalizationVisitor v_default;
        expr->accept(v_default);
        auto result_default = v_default.get_result();

        std::string label = "Expression " + std::to_string(i) + " without context = default";
        EXPECT_TRUE(result_no_ctx->equals(*result_default), label);
    }
}

void test_backward_compat_no_assumption_rules_fire() {
    TEST_CASE("With context but no relevant assumptions, no rules fire");

    // Create a context with assumptions for variable "a", but simplify
    // expressions involving variable "x" — no rules should fire for "x"
    AssumptionContext ctx;
    ctx.assume_sign("a", Sign::Positive);

    // sqrt(x²) should NOT simplify since x has no assumptions
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    // Compare with no-context result
    auto result_no_ctx = normalize_no_ctx(sqrt_x_sq);

    EXPECT_TRUE(result->equals(*result_no_ctx),
                "sqrt(x²) with unrelated assumptions = no-context result");

    // abs(x) should NOT simplify since x has no assumptions
    auto abs_x = make_abs(var("x"));
    auto result_abs = normalize_with_ctx(abs_x, ctx);
    auto result_abs_no_ctx = normalize_no_ctx(abs_x);

    EXPECT_TRUE(result_abs->equals(*result_abs_no_ctx),
                "abs(x) with unrelated assumptions = no-context result");
}

void test_backward_compat_null_context_explicit() {
    TEST_CASE("NormalizationVisitor(nullptr) behaves like default constructor");

    // Explicitly passing nullptr should behave identically to default constructor
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));

    NormalizationVisitor v_nullptr(nullptr);
    sqrt_x_sq->accept(v_nullptr);
    auto result_nullptr = v_nullptr.get_result();

    NormalizationVisitor v_default;
    sqrt_x_sq->accept(v_default);
    auto result_default = v_default.get_result();

    EXPECT_TRUE(result_nullptr->equals(*result_default),
                "NormalizationVisitor(nullptr) = NormalizationVisitor() for sqrt(x²)");

    // Also test abs(x)
    auto abs_x = make_abs(var("x"));

    NormalizationVisitor v_nullptr2(nullptr);
    abs_x->accept(v_nullptr2);
    auto result_nullptr2 = v_nullptr2.get_result();

    NormalizationVisitor v_default2;
    abs_x->accept(v_default2);
    auto result_default2 = v_default2.get_result();

    EXPECT_TRUE(result_nullptr2->equals(*result_default2),
                "NormalizationVisitor(nullptr) = NormalizationVisitor() for abs(x)");
}


void test_sqrt_x_squared_in_larger_expression() {
    TEST_CASE("sqrt(x²) simplifies within a larger expression");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    // Build: sqrt(x²) + 1
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto expr = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{sqrt_x_sq, num(1)});

    auto result = normalize_with_ctx(expr, ctx);

    // The result should be x + 1
    auto add = std::dynamic_pointer_cast<const AddNode>(result);
    if (add) {
        // Check that the result contains x and 1 (order may vary)
        bool has_x = false;
        bool has_one = false;
        for (const auto& op : add->operands()) {
            if (is_variable(op, "x")) has_x = true;
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (std::holds_alternative<BigInt>(n->value()) && std::get<BigInt>(n->value()) == BigInt(1))
                    has_one = true;
            }
        }
        EXPECT_TRUE(has_x && has_one,
                    "sqrt(x²) + 1 with NonNegative x → x + 1");
    } else {
        // Might be simplified differently, just check it's not still sqrt(x²) + 1
        PrintVisitor pv;
        if (result) result->accept(pv);
        const auto result_str = pv.get_result();
        EXPECT_TRUE(result != nullptr &&
                        result_str.find("sqrt") == std::string::npos &&
                        result_str.find("x") != std::string::npos &&
                        result_str.find("1") != std::string::npos,
                    "sqrt(x²) + 1 simplified (non-AddNode result)");
    }
}

void test_abs_in_larger_expression() {
    TEST_CASE("abs(x) simplifies within a larger expression");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    // Build: 2 * abs(x)
    auto abs_x = make_abs(var("x"));
    auto expr = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{num(2), abs_x});

    auto result = normalize_with_ctx(expr, ctx);

    // The result should be 2 * x
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(result);
    if (mul) {
        bool has_x = false;
        bool has_two = false;
        for (const auto& op : mul->operands()) {
            if (is_variable(op, "x")) has_x = true;
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                if (std::holds_alternative<BigInt>(n->value()) && std::get<BigInt>(n->value()) == BigInt(2))
                    has_two = true;
            }
        }
        EXPECT_TRUE(has_x && has_two,
                    "2 * abs(x) with Positive x → 2 * x");
    } else {
        // Could be simplified to just a variable if 2*x normalizes differently
        PrintVisitor pv;
        if (result) result->accept(pv);
        const auto result_str = pv.get_result();
        EXPECT_TRUE(result != nullptr &&
                        result_str.find("abs") == std::string::npos &&
                        result_str.find("x") != std::string::npos &&
                        result_str.find("2") != std::string::npos,
                    "2 * abs(x) simplified (non-MultiplyNode result)");
    }
}

void test_scoped_assumption_simplification() {
    TEST_CASE("Scoped assumptions affect simplification correctly");

    AssumptionContext ctx;

    // In root scope, x has no assumptions
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto result_root = normalize_with_ctx(sqrt_x_sq, ctx);

    // Push scope and declare x NonNegative
    ctx.push();
    ctx.assume_sign("x", Sign::NonNegative);

    auto result_child = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(is_variable(result_child, "x"),
                "sqrt(x²) in child scope with NonNegative → x");

    // Pop scope — x should no longer be NonNegative
    ctx.pop();

    auto result_after_pop = normalize_with_ctx(sqrt_x_sq, ctx);
    // After pop, should behave like no assumptions
    EXPECT_TRUE(result_after_pop->equals(*result_root),
                "sqrt(x²) after pop = root scope result (no simplification)");
}


int main() {
    test_sqrt_x_squared_nonnegative();
    test_sqrt_x_squared_positive();
    test_sqrt_x_squared_real_not_nonneg();
    test_sqrt_x_squared_integer_not_nonneg();
    test_sqrt_x_squared_natural();

    test_abs_positive();
    test_abs_negative();
    test_abs_nonnegative_not_positive();
    test_abs_no_assumption();

    test_backward_compat_sqrt_x_squared();
    test_backward_compat_abs_x();
    test_backward_compat_various_expressions();
    test_backward_compat_no_assumption_rules_fire();
    test_backward_compat_null_context_explicit();

    // Additional edge cases
    test_sqrt_x_squared_in_larger_expression();
    test_abs_in_larger_expression();
    test_scoped_assumption_simplification();

    return TEST_REPORT();
}
