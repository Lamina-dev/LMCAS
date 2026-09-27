
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "query_interface.hpp"
#include "solver.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>

using namespace LMCAS;

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

/// Normalize a node with an AssumptionContext.
static std::shared_ptr<const SymbolicNode> normalize_with_ctx(
    const std::shared_ptr<const SymbolicNode> &node,
    const AssumptionContext &ctx) {
    auto expression = LMCAS::detail::expression_from_node(node);
    return LMCAS::detail::node(ctx.simplify(expression));
}

/// Check if a node is a VariableNode with the given name.
static bool is_variable(const std::shared_ptr<const SymbolicNode> &node, const std::string &name) {
    auto v = std::dynamic_pointer_cast<const VariableNode>(node);
    return v && v->name() == name;
}

/// Check if a node is abs(x) — FunctionNode(Abs, {VariableNode(name)}).
static bool is_abs_of_var(const std::shared_ptr<const SymbolicNode> &node, const std::string &name) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!func || func->type() != FunctionNode::FuncType::Abs)
        return false;
    if (func->arguments().size() != 1)
        return false;
    return is_variable(func->arguments()[0], name);
}

/// Try to extract a numeric double value from a solution expression.
static bool try_numeric(const std::shared_ptr<SymbolicExpr> &expr, double &out) {
    if (!expr || !LMCAS::detail::node(expr))
        return false;
    auto n = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
    if (!n)
        return false;

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
    const std::vector<std::shared_ptr<SymbolicExpr>> &solutions,
    double target, double tol = 1e-9) {
    for (const auto &sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            if (std::abs(v - target) < tol)
                return true;
        }
    }
    return false;
}

TEST(AssumptionExamples, EndToEndSqrtXSquaredNonneg) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    EXPECT_TRUE((is_variable(result, "x"))) << "End-to-end: sqrt(x^2) with x>=0 simplifies to x";
}

TEST(AssumptionExamples, EndToEndAbsNegativeVar) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());

    auto abs_x = make_abs(var("x"));
    auto result = normalize_with_ctx(abs_x, ctx);

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(result);
    EXPECT_TRUE((mul != nullptr)) << "End-to-end: abs(x) with x<0 is a MultiplyNode";
    if (mul && mul->operands().size() == 2) {
        auto coeff = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[0]);
        bool is_neg_one = false;
        if (coeff) {
            if (std::holds_alternative<BigInt>(coeff->value()))
                is_neg_one = (std::get<BigInt>(coeff->value()) == BigInt(-1));
        }
        EXPECT_TRUE((is_neg_one)) << "End-to-end: abs(x) with x<0 has coefficient -1";
        EXPECT_TRUE((is_variable(mul->operands()[1], "x"))) << "End-to-end: abs(x) with x<0 has variable x";
    }
}

TEST(AssumptionExamples, EndToEndSqrtXSquaredReal) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    auto x_squared = make_power(var("x"), 2);
    auto sqrt_x_sq = make_sqrt(x_squared);

    auto result = normalize_with_ctx(sqrt_x_sq, ctx);

    EXPECT_TRUE((is_abs_of_var(result, "x"))) << "End-to-end: sqrt(x^2) with x Real simplifies to abs(x)";
}

TEST(AssumptionExamples, EndToEndQueryAfterAssumption) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    auto x_expr = LMCAS::detail::expression_from_node(var("x"));
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::True)) << "End-to-end: x is Positive after assume_sign(Positive)";
    EXPECT_TRUE((ctx.is_nonnegative(x_expr).value() == Tribool::True)) << "End-to-end: x is NonNegative (implied by Positive)";
    EXPECT_TRUE((ctx.is_nonzero(x_expr).value() == Tribool::True)) << "End-to-end: x is NonZero (implied by Positive)";
    EXPECT_TRUE((ctx.is_negative(x_expr).value() == Tribool::False)) << "End-to-end: x is NOT Negative when Positive";
}

TEST(AssumptionExamples, SolverWithPositiveIntDomain) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-4));

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE((has_2)) << "Solver+assumptions: x^2-4=0 PositiveInt contains x=2";
    EXPECT_FALSE((has_neg2)) << "Solver+assumptions: x^2-4=0 PositiveInt excludes x=-2";
}

TEST(AssumptionExamples, SolverWithNonnegativeSign) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-9));

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_3 = solutions_contain_value(solutions, 3.0);
    bool has_neg3 = solutions_contain_value(solutions, -3.0);

    EXPECT_TRUE((has_3)) << "Solver+assumptions: x^2-9=0 NonNeg contains x=3";
    EXPECT_FALSE((has_neg3)) << "Solver+assumptions: x^2-9=0 NonNeg excludes x=-3";
}

TEST(AssumptionExamples, NestedScopesQueryRoundtrip) {
    AssumptionContext ctx;
    auto x_expr = LMCAS::detail::expression_from_node(var("x"));
    // Root scope: x has no assumptions
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::Unknown)) << "Nested scopes: x is Unknown in root scope";

    // Push scope 1: assume x > 0
    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::True)) << "Nested scopes: x is Positive in scope 1";

    // Push scope 2: assume y is Integer (x should still be accessible via has_sign)
    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("y", Domain::Integer).has_value());

    auto y_expr = LMCAS::detail::expression_from_node(var("y"));
    EXPECT_TRUE((ctx.is_integer(y_expr).value() == Tribool::True)) << "Nested scopes: y is Integer in scope 2";
    // x Positive is visible via read-through (has_sign reads all scopes)
    EXPECT_TRUE((ctx.has_sign("x", Sign::Positive))) << "Nested scopes: x still Positive in scope 2 (read-through)";

    // Pop scope 2
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.is_integer(y_expr).value() == Tribool::Unknown)) << "Nested scopes: y is Unknown after popping scope 2";
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::True)) << "Nested scopes: x still Positive in scope 1";

    // Pop scope 1
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::Unknown)) << "Nested scopes: x is Unknown after popping scope 1";
}

TEST(AssumptionExamples, NestedScopesSimplificationChanges) {
    AssumptionContext ctx;
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));

    // Root scope: no assumptions, sqrt(x^2) stays as-is
    auto result_root = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_FALSE((is_variable(result_root, "x"))) << "Nested scopes: sqrt(x^2) does NOT simplify to x in root";

    // Push and assume x is Real
    ctx.push();
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    auto result_real = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((is_abs_of_var(result_real, "x"))) << "Nested scopes: sqrt(x^2) -> abs(x) with Real assumption";

    // Push deeper and assume x >= 0
    ctx.push();
    EXPECT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());
    auto result_nonneg = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((is_variable(result_nonneg, "x"))) << "Nested scopes: sqrt(x^2) -> x with NonNegative assumption";

    // Pop back to Real-only scope
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    auto result_after_pop = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((is_abs_of_var(result_after_pop, "x"))) << "Nested scopes: sqrt(x^2) -> abs(x) after popping NonNeg scope";

    // Pop back to root
    EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    auto result_final = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((result_final->equals(*result_root))) << "Nested scopes: sqrt(x^2) back to original after all pops";
}

TEST(AssumptionExamples, UnrecognizedFunctionReturnsUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    // LambertW is not in the recognized built-in list for property inference
    auto lambert_w = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::LambertW,
        std::vector<std::shared_ptr<const SymbolicNode>>{var("x")});
    auto expr = LMCAS::detail::expression_from_node(lambert_w);
    EXPECT_TRUE((ctx.is_positive(expr).value() == Tribool::Unknown)) << "LambertW(x) is_positive -> Unknown";
    EXPECT_TRUE((ctx.is_negative(expr).value() == Tribool::Unknown)) << "LambertW(x) is_negative -> Unknown";
    EXPECT_TRUE((ctx.is_nonnegative(expr).value() == Tribool::Unknown)) << "LambertW(x) is_nonnegative -> Unknown";
    EXPECT_TRUE((ctx.is_real(expr).value() == Tribool::Unknown)) << "LambertW(x) is_real -> Unknown";
    EXPECT_TRUE((ctx.is_integer(expr).value() == Tribool::Unknown)) << "LambertW(x) is_integer -> Unknown";
    EXPECT_TRUE((ctx.is_nonzero(expr).value() == Tribool::Unknown)) << "LambertW(x) is_nonzero -> Unknown";
}

TEST(AssumptionExamples, UnrecognizedFunctionErf) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    auto erf_x = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Erf,
        std::vector<std::shared_ptr<const SymbolicNode>>{var("x")});
    auto expr = LMCAS::detail::expression_from_node(erf_x);
    EXPECT_TRUE((ctx.is_positive(expr).value() == Tribool::Unknown)) << "Erf(x) is_positive -> Unknown";
    EXPECT_TRUE((ctx.is_negative(expr).value() == Tribool::Unknown)) << "Erf(x) is_negative -> Unknown";
    EXPECT_TRUE((ctx.is_real(expr).value() == Tribool::Unknown)) << "Erf(x) is_real -> Unknown";
    EXPECT_TRUE((ctx.is_integer(expr).value() == Tribool::Unknown)) << "Erf(x) is_integer -> Unknown";
}

TEST(AssumptionExamples, DomainFilteringAllExcludedEmptySet) {
    // x^2 + 1 = 0 has only imaginary solutions (x = i, x = -i)
    // With Real domain, all solutions should be excluded
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(1));

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    EXPECT_TRUE((solutions.empty())) << "x^2+1=0 with Real domain -> empty set";
}

TEST(AssumptionExamples, DomainFilteringPositiveIntExcludesAll) {
    // x^2 - 2 = 0 -> x = sqrt(2), x = -sqrt(2)
    // Neither is a positive integer
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-2));

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    // All numeric solutions should be excluded (sqrt(2) is not an integer)
    for (const auto &sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            ADD_FAILURE() << "No numeric solutions should pass PositiveInt filter";
        }
    }
}

TEST(AssumptionExamples, PopOnRootScopeThrows) {
    AssumptionContext ctx;
    EXPECT_TRUE((ctx.depth() == 1)) << "Initial depth is 1 (root scope)";

    auto failure_416 = ctx.pop();
    EXPECT_TRUE((!failure_416.has_value())) << "pop() on root scope returns InvalidArgument";
    EXPECT_TRUE((ctx.depth() == 1)) << "depth unchanged after failed pop";
}

TEST(AssumptionExamples, NestingDepth128) {
    AssumptionContext ctx;

    // Push 128 times
    for (int i = 0; i < 128; ++i) {
        ctx.push();
    }

    EXPECT_TRUE((ctx.depth() == 129)) << "Depth is 129 after 128 pushes (root + 128)";

    // Declare something in the deepest scope
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    auto x_expr = LMCAS::detail::expression_from_node(var("x"));
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::True)) << "Can declare and query at depth 129";

    // Pop all 128 scopes
    for (int i = 0; i < 128; ++i) {
        EXPECT_TRUE((ctx.pop().has_value())) << "scope pop succeeds";
    }

    EXPECT_TRUE((ctx.depth() == 1)) << "Depth is 1 after popping all 128 scopes";
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::Unknown)) << "x is Unknown after popping all scopes";
}

TEST(AssumptionExamples, NanHandling) {
    bool rejected = false;
    try {
        (void)LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(
                std::numeric_limits<double>::quiet_NaN()));
    } catch (const std::invalid_argument &error) {
        rejected =
            std::string(error.what()) == "approximate number must be finite";
    }
    EXPECT_TRUE((rejected)) << "NaN cannot enter the symbolic AST";
}

TEST(AssumptionExamples, InfinityHandling) {
    AssumptionContext ctx;

    // Positive infinity: FunctionNode(Infinity, {})
    auto pos_inf = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    auto pos_inf_expr = LMCAS::detail::expression_from_node(pos_inf);
    EXPECT_TRUE((ctx.is_integer(pos_inf_expr).value() == Tribool::False)) << "+Infinity is_integer -> False";
    EXPECT_TRUE((ctx.is_positive(pos_inf_expr).value() == Tribool::True)) << "+Infinity is_positive -> True";

    // Negative infinity: MultiplyNode(-1, Infinity)
    auto neg_inf = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{num(-1), pos_inf});
    auto neg_inf_expr = LMCAS::detail::expression_from_node(neg_inf);
    EXPECT_TRUE((ctx.is_integer(neg_inf_expr).value() == Tribool::False)) << "-Infinity is_integer -> False";
    EXPECT_TRUE((ctx.is_negative(neg_inf_expr).value() == Tribool::True)) << "-Infinity is_negative -> True";
}

TEST(AssumptionExamples, CombinedPipeline) {
    AssumptionContext ctx;

    // Step 1: Declare x > 0 and x is Real
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    // Step 2: Query properties
    auto x_expr = LMCAS::detail::expression_from_node(var("x"));
    EXPECT_TRUE((ctx.is_positive(x_expr).value() == Tribool::True)) << "Pipeline: x is Positive";
    EXPECT_TRUE((ctx.is_real(x_expr).value() == Tribool::True)) << "Pipeline: x is Real";

    // Step 3: Simplify abs(x) -> x (since x > 0)
    auto abs_x = make_abs(var("x"));
    auto simplified = normalize_with_ctx(abs_x, ctx);
    EXPECT_TRUE((is_variable(simplified, "x"))) << "Pipeline: abs(x) simplifies to x when x > 0";

    // Step 4: Simplify sqrt(x^2) -> x (since x >= 0 implied by Positive)
    auto sqrt_x_sq = make_sqrt(make_power(var("x"), 2));
    auto simplified2 = normalize_with_ctx(sqrt_x_sq, ctx);
    EXPECT_TRUE((is_variable(simplified2, "x"))) << "Pipeline: sqrt(x^2) simplifies to x when x > 0";

    // Step 5: Solve x^2 - 1 = 0 with Positive constraint -> only x=1
    auto eq = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);
    EXPECT_TRUE((has_1)) << "Pipeline: x^2-1=0 with Positive contains x=1";
    EXPECT_FALSE((has_neg1)) << "Pipeline: x^2-1=0 with Positive excludes x=-1";
}
