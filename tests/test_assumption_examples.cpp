
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "query_interface.hpp"
#include "solver.hpp"
#include "visitors/normalization_visitor.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>

using namespace lamina;


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

/// Normalize a node with an AssumptionContext.
static std::shared_ptr<const SymbolicNode> normalize_with_ctx(
    const std::shared_ptr<const SymbolicNode>& node,
    const AssumptionContext& ctx) {
    NormalizationVisitor v(&ctx);
    node->accept(v);
    return v.get_result();
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

    auto n = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0]);
    if (!n) return false;

    bool is_neg_one = false;
    if (std::holds_alternative<BigInt>(n->value()))
        is_neg_one = (std::get<BigInt>(n->value()) == BigInt(-1));
    else if (std::holds_alternative<lmmc_real_t>(n->value()))
        is_neg_one = (std::get<lmmc_real_t>(n->value()) == -1.0);
    else if (std::holds_alternative<Rational>(n->value()))
        is_neg_one = (std::get<Rational>(n->value()) == Rational(-1));

    if (!is_neg_one) return false;
    return is_variable(mul->operands()[1], name);
}

/// Try to extract a numeric double value from a solution expression.
static bool try_numeric(const std::shared_ptr<SymbolicExpr>& expr, double& out) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto n = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!n) return false;

    if (std::holds_alternative<BigInt>(n->value())) {
        out = std::get<BigInt>(n->value()).to_double();
        return true;
    }
    if (std::holds_alternative<Rational>(n->value())) {
        out = std::get<Rational>(n->value()).to_double();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(n->value())) {
        out = std::get<lmmc_real_t>(n->value());
        return true;
    }
    return false;
}

/// Check if a solution set contains a numeric value (within tolerance).
static bool solutions_contain_value(
    const std::vector<std::shared_ptr<SymbolicExpr>>& solutions,
    double target, double tol = 1e-9)
{
    for (const auto& sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            if (std::abs(v - target) < tol) return true;
        }
    }
    return false;
}


void test_end_to_end_sqrt_x_squared_nonneg() {
    TEST_CASE("Integration: assume x >= 0, simplify sqrt(x^2), get x");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    EXPECT_TRUE(is_variable(result, "x"),
                "End-to-end: sqrt(x^2) with x>=0 simplifies to x");
}

void test_end_to_end_abs_negative_var() {
    TEST_CASE("Integration: assume x < 0, simplify abs(x), get -x");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(result);
    EXPECT_TRUE(mul != nullptr, "End-to-end: abs(x) with x<0 is a MultiplyNode");
    if (mul && mul->operands().size() == 2) {
        auto coeff = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0]);
        bool is_neg_one = false;
        if (coeff) {
            if (std::holds_alternative<BigInt>(coeff->value()))
                is_neg_one = (std::get<BigInt>(coeff->value()) == BigInt(-1));
        }
        EXPECT_TRUE(is_neg_one, "End-to-end: abs(x) with x<0 has coefficient -1");
        EXPECT_TRUE(is_variable(mul->operands()[1], "x"),
                    "End-to-end: abs(x) with x<0 has variable x");
    }
}

void test_end_to_end_sqrt_x_squared_real() {
    TEST_CASE("Integration: assume x is Real, simplify sqrt(x^2), get abs(x)");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    EXPECT_TRUE(is_abs_of_var(result, "x"),
                "End-to-end: sqrt(x^2) with x Real simplifies to abs(x)");
}

void test_end_to_end_query_after_assumption() {
    TEST_CASE("Integration: assume x > 0, query properties");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    auto x_expr = lamina::detail::expression_from_node(var("x"));
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::True,
                "End-to-end: x is Positive after assume_sign(Positive)");
    EXPECT_TRUE(ctx.is_nonnegative(x_expr).value() == Tribool::True,
                "End-to-end: x is NonNegative (implied by Positive)");
    EXPECT_TRUE(ctx.is_nonzero(x_expr).value() == Tribool::True,
                "End-to-end: x is NonZero (implied by Positive)");
    EXPECT_TRUE(ctx.is_negative(x_expr).value() == Tribool::False,
                "End-to-end: x is NOT Negative when Positive");
}


void test_solver_with_positive_int_domain() {
    TEST_CASE("Integration: solve x^2 - 4 = 0 with PositiveInt, only x=2");

    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-4));

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::PositiveInt);

    auto solutions = lamina::detail::propagate_result(solve_with_assumptions_checked(eq, "x", &ctx));

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "Solver+assumptions: x^2-4=0 PositiveInt contains x=2");
    EXPECT_FALSE(has_neg2, "Solver+assumptions: x^2-4=0 PositiveInt excludes x=-2");
}

void test_solver_with_nonnegative_sign() {
    TEST_CASE("Integration: solve x^2 - 9 = 0 with NonNegative, only x=3");

    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-9));

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::NonNegative);

    auto solutions = lamina::detail::propagate_result(solve_with_assumptions_checked(eq, "x", &ctx));

    bool has_3 = solutions_contain_value(solutions, 3.0);
    bool has_neg3 = solutions_contain_value(solutions, -3.0);

    EXPECT_TRUE(has_3, "Solver+assumptions: x^2-9=0 NonNeg contains x=3");
    EXPECT_FALSE(has_neg3, "Solver+assumptions: x^2-9=0 NonNeg excludes x=-3");
}


void test_nested_scopes_query_roundtrip() {
    TEST_CASE("Integration: nested scopes push/assume/query/pop/query");

    AssumptionContext ctx;
    auto x_expr = lamina::detail::expression_from_node(var("x"));
    // Root scope: x has no assumptions
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::Unknown,
                "Nested scopes: x is Unknown in root scope");

    // Push scope 1: assume x > 0
    ctx.push();
    ctx.assume_sign("x", Sign::Positive);

    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::True,
                "Nested scopes: x is Positive in scope 1");

    // Push scope 2: assume y is Integer (x should still be accessible via has_sign)
    ctx.push();
    ctx.assume_domain("y", Domain::Integer);

    auto y_expr = lamina::detail::expression_from_node(var("y"));
    EXPECT_TRUE(ctx.is_integer(y_expr).value() == Tribool::True,
                "Nested scopes: y is Integer in scope 2");
    // x Positive is visible via read-through (has_sign reads all scopes)
    EXPECT_TRUE(ctx.has_sign("x", Sign::Positive),
                "Nested scopes: x still Positive in scope 2 (read-through)");

    // Pop scope 2
    ctx.pop();
    EXPECT_TRUE(ctx.is_integer(y_expr).value() == Tribool::Unknown,
                "Nested scopes: y is Unknown after popping scope 2");
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::True,
                "Nested scopes: x still Positive in scope 1");

    // Pop scope 1
    ctx.pop();
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::Unknown,
                "Nested scopes: x is Unknown after popping scope 1");
}

void test_nested_scopes_simplification_changes() {
    TEST_CASE("Integration: nested scopes affect simplification differently");

    AssumptionContext ctx;
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));

    // Root scope: no assumptions, sqrt(x^2) stays as-is
    auto result_root = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_FALSE(is_variable(result_root, "x"),
                 "Nested scopes: sqrt(x^2) does NOT simplify to x in root");

    // Push and assume x is Real
    ctx.push();
    ctx.assume_domain("x", Domain::Real);
    auto result_real = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(is_abs_of_var(result_real, "x"),
                "Nested scopes: sqrt(x^2) -> abs(x) with Real assumption");

    // Push deeper and assume x >= 0
    ctx.push();
    ctx.assume_sign("x", Sign::NonNegative);
    auto result_nonneg = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(is_variable(result_nonneg, "x"),
                "Nested scopes: sqrt(x^2) -> x with NonNegative assumption");

    // Pop back to Real-only scope
    ctx.pop();
    auto result_after_pop = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(is_abs_of_var(result_after_pop, "x"),
                "Nested scopes: sqrt(x^2) -> abs(x) after popping NonNeg scope");

    // Pop back to root
    ctx.pop();
    auto result_final = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(result_final->equals(*result_root),
                "Nested scopes: sqrt(x^2) back to original after all pops");
}


void test_unrecognized_function_returns_unknown() {
    TEST_CASE("Unrecognized function returns Unknown for all properties");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    // LambertW is not in the recognized built-in list for property inference
    auto lambert_w = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::LambertW,
        std::vector<std::shared_ptr<const SymbolicNode>>{var("x")});
    auto expr = lamina::detail::expression_from_node(lambert_w);
    EXPECT_TRUE(ctx.is_positive(expr).value() == Tribool::Unknown,
                "LambertW(x) is_positive -> Unknown");
    EXPECT_TRUE(ctx.is_negative(expr).value() == Tribool::Unknown,
                "LambertW(x) is_negative -> Unknown");
    EXPECT_TRUE(ctx.is_nonnegative(expr).value() == Tribool::Unknown,
                "LambertW(x) is_nonnegative -> Unknown");
    EXPECT_TRUE(ctx.is_real(expr).value() == Tribool::Unknown,
                "LambertW(x) is_real -> Unknown");
    EXPECT_TRUE(ctx.is_integer(expr).value() == Tribool::Unknown,
                "LambertW(x) is_integer -> Unknown");
    EXPECT_TRUE(ctx.is_nonzero(expr).value() == Tribool::Unknown,
                "LambertW(x) is_nonzero -> Unknown");
}

void test_unrecognized_function_erf() {
    TEST_CASE("Erf (special function) returns Unknown for all properties");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto erf_x = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Erf,
        std::vector<std::shared_ptr<const SymbolicNode>>{var("x")});
    auto expr = lamina::detail::expression_from_node(erf_x);
    EXPECT_TRUE(ctx.is_positive(expr).value() == Tribool::Unknown,
                "Erf(x) is_positive -> Unknown");
    EXPECT_TRUE(ctx.is_negative(expr).value() == Tribool::Unknown,
                "Erf(x) is_negative -> Unknown");
    EXPECT_TRUE(ctx.is_real(expr).value() == Tribool::Unknown,
                "Erf(x) is_real -> Unknown");
    EXPECT_TRUE(ctx.is_integer(expr).value() == Tribool::Unknown,
                "Erf(x) is_integer -> Unknown");
}


void test_domain_filtering_all_excluded_empty_set() {
    TEST_CASE("Domain filtering excludes all solutions -> empty set");

    // x^2 + 1 = 0 has only imaginary solutions (x = i, x = -i)
    // With Real domain, all solutions should be excluded
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(1));

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto solutions = lamina::detail::propagate_result(solve_with_assumptions_checked(eq, "x", &ctx));

    EXPECT_TRUE(solutions.empty(),
                "x^2+1=0 with Real domain -> empty set");
}

void test_domain_filtering_positive_int_excludes_all() {
    TEST_CASE("x^2 - 2 = 0 with PositiveInt -> empty set (sqrt(2) not integer)");

    // x^2 - 2 = 0 -> x = sqrt(2), x = -sqrt(2)
    // Neither is a positive integer
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-2));

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::PositiveInt);

    auto solutions = lamina::detail::propagate_result(solve_with_assumptions_checked(eq, "x", &ctx));

    // All numeric solutions should be excluded (sqrt(2) is not an integer)
    for (const auto& sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            EXPECT_TRUE(false,
                        "No numeric solutions should pass PositiveInt filter");
        }
    }
}


void test_pop_on_root_scope_throws() {
    TEST_CASE("pop() on root scope returns InvalidArgument");

    AssumptionContext ctx;
    EXPECT_TRUE(ctx.depth() == 1, "Initial depth is 1 (root scope)");

    auto failure_416 = ctx.pop();
    EXPECT_TRUE(!failure_416.has_value(), "pop() on root scope returns InvalidArgument");
    EXPECT_TRUE(ctx.depth() == 1, "depth unchanged after failed pop");
}


void test_nesting_depth_128() {
    TEST_CASE("Nesting depth of 128 push() calls supported");

    AssumptionContext ctx;

    // Push 128 times
    for (int i = 0; i < 128; ++i) {
        ctx.push();
    }

    EXPECT_TRUE(ctx.depth() == 129, "Depth is 129 after 128 pushes (root + 128)");

    // Declare something in the deepest scope
    ctx.assume_sign("x", Sign::Positive);
    auto x_expr = lamina::detail::expression_from_node(var("x"));
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::True,
                "Can declare and query at depth 129");

    // Pop all 128 scopes
    for (int i = 0; i < 128; ++i) {
        ctx.pop();
    }

    EXPECT_TRUE(ctx.depth() == 1, "Depth is 1 after popping all 128 scopes");
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::Unknown,
                "x is Unknown after popping all scopes");
}


static void test_nan_handling() {
    TEST_CASE("NaN NumberNode construction is rejected");

    bool rejected = false;
    try {
        (void)lamina::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(
                std::numeric_limits<double>::quiet_NaN()));
    } catch (const std::invalid_argument& error) {
        rejected =
            std::string(error.what()) == "approximate number must be finite";
    }
    EXPECT_TRUE(rejected,
                "NaN cannot enter the symbolic AST");
}

static void test_infinity_handling() {
    TEST_CASE("Infinity -> derive sign, False for integer");

    AssumptionContext ctx;

    // Positive infinity: FunctionNode(Infinity, {})
    auto pos_inf = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    auto pos_inf_expr = lamina::detail::expression_from_node(pos_inf);
    EXPECT_TRUE(ctx.is_integer(pos_inf_expr).value() == Tribool::False,
                "+Infinity is_integer -> False");
    EXPECT_TRUE(ctx.is_positive(pos_inf_expr).value() == Tribool::True,
                "+Infinity is_positive -> True");

    // Negative infinity: MultiplyNode(-1, Infinity)
    auto neg_inf = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{num(-1), pos_inf});
    auto neg_inf_expr = lamina::detail::expression_from_node(neg_inf);
    EXPECT_TRUE(ctx.is_integer(neg_inf_expr).value() == Tribool::False,
                "-Infinity is_integer -> False");
    EXPECT_TRUE(ctx.is_negative(neg_inf_expr).value() == Tribool::True,
                "-Infinity is_negative -> True");
}


void test_combined_pipeline() {
    TEST_CASE("Integration: combined pipeline - assume, simplify, query, solve");

    AssumptionContext ctx;

    // Step 1: Declare x > 0 and x is Real
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    // Step 2: Query properties
    auto x_expr = lamina::detail::expression_from_node(var("x"));
    EXPECT_TRUE(ctx.is_positive(x_expr).value() == Tribool::True,
                "Pipeline: x is Positive");
    EXPECT_TRUE(ctx.is_real(x_expr).value() == Tribool::True,
                "Pipeline: x is Real");

    // Step 3: Simplify abs(x) -> x (since x > 0)
    auto abs_x = make_abs(var("x"));
    auto simplified = normalize_with_ctx(abs_x, ctx);
    EXPECT_TRUE(is_variable(simplified, "x"),
                "Pipeline: abs(x) simplifies to x when x > 0");

    // Step 4: Simplify sqrt(x^2) -> x (since x >= 0 implied by Positive)
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto simplified2 = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE(is_variable(simplified2, "x"),
                "Pipeline: sqrt(x^2) simplifies to x when x > 0");

    // Step 5: Solve x^2 - 1 = 0 with Positive constraint -> only x=1
    auto eq = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    auto solutions = lamina::detail::propagate_result(solve_with_assumptions_checked(eq, "x", &ctx));

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);
    EXPECT_TRUE(has_1, "Pipeline: x^2-1=0 with Positive contains x=1");
    EXPECT_FALSE(has_neg1, "Pipeline: x^2-1=0 with Positive excludes x=-1");
}


int main() {
    // Integration Test 1: End-to-end assumption -> simplification
    test_end_to_end_sqrt_x_squared_nonneg();
    test_end_to_end_abs_negative_var();
    test_end_to_end_sqrt_x_squared_real();
    test_end_to_end_query_after_assumption();

    // Integration Test 2: Solver with assumptions
    test_solver_with_positive_int_domain();
    test_solver_with_nonnegative_sign();

    // Integration Test 3: Nested scopes
    test_nested_scopes_query_roundtrip();
    test_nested_scopes_simplification_changes();

    test_unrecognized_function_returns_unknown();
    test_unrecognized_function_erf();

    test_domain_filtering_all_excluded_empty_set();
    test_domain_filtering_positive_int_excludes_all();

    // Edge Cases
    test_pop_on_root_scope_throws();
    test_nesting_depth_128();
    test_nan_handling();
    test_infinity_handling();

    // Combined pipeline test
    test_combined_pipeline();

    return TEST_REPORT();
}
