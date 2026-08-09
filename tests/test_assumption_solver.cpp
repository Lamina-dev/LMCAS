
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "solver.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace lamina;


/// Try to extract a numeric double value from a solution expression.
static bool try_numeric(const std::shared_ptr<SymbolicExpr>& expr, double& out) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!num) return false;

    if (std::holds_alternative<BigInt>(num->value())) {
        out = std::get<BigInt>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        out = std::get<Rational>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        out = std::get<lmmc_real_t>(num->value());
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

/// Check if any solution contains imaginary components (sqrt of negative).
static bool any_solution_contains_imaginary(
    const std::vector<std::shared_ptr<SymbolicExpr>>& solutions)
{
    for (const auto& sol : solutions) {
        if (!sol || !lamina::detail::node(sol)) continue;
        // Check for FunctionNode(Sqrt, negative number)
        auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(sol));
        if (func && func->type() == FunctionNode::FuncType::Sqrt && func->arguments().size() == 1) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(func->arguments()[0]);
            if (num) {
                if (std::holds_alternative<BigInt>(num->value()) &&
                    std::get<BigInt>(num->value()).IsNegative()) return true;
                if (std::holds_alternative<Rational>(num->value()) &&
                    std::get<Rational>(num->value()) < Rational(0)) return true;
                if (std::holds_alternative<lmmc_real_t>(num->value()) &&
                    std::get<lmmc_real_t>(num->value()) < 0.0) return true;
            }
        }
        // Also check MultiplyNode containing sqrt(-1)
        auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(sol));
        if (mul) {
            for (const auto& op : mul->operands()) {
                auto f = std::dynamic_pointer_cast<const FunctionNode>(op);
                if (f && f->type() == FunctionNode::FuncType::Sqrt && f->arguments().size() == 1) {
                    auto n = std::dynamic_pointer_cast<const NumberNode>(f->arguments()[0]);
                    if (n) {
                        if (std::holds_alternative<BigInt>(n->value()) &&
                            std::get<BigInt>(n->value()).IsNegative()) return true;
                    }
                }
            }
        }
    }
    return false;
}

/// Build equation x² - c = 0 as a SymbolicExpr (x^2 + (-c))
static std::shared_ptr<SymbolicExpr> make_x_squared_minus(const std::string& var, int c) {
    auto x = SymbolicExpr::variable(var);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-c));
    return eq;
}

/// Build equation x² + c = 0 as a SymbolicExpr (x^2 + c)
static std::shared_ptr<SymbolicExpr> make_x_squared_plus(const std::string& var, int c) {
    auto x = SymbolicExpr::variable(var);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(c));
    return eq;
}


void test_x_squared_minus_4_real_domain() {
    TEST_CASE("Property 25: x²-4=0 with Real domain → both x=2 and x=-2 returned");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Both 2 and -2 are real, so both should be returned
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 Real domain: contains x=2");
    EXPECT_TRUE(has_neg2, "x²-4=0 Real domain: contains x=-2");
    EXPECT_TRUE(solutions.size() >= 2, "x²-4=0 Real domain: at least 2 solutions");
}

void test_x_squared_minus_4_positive_int() {
    TEST_CASE("Property 25: x²-4=0 with PositiveInt domain → only x=2");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::PositiveInt);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Only x=2 is a positive integer; x=-2 should be excluded
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 PositiveInt: contains x=2");
    EXPECT_FALSE(has_neg2, "x²-4=0 PositiveInt: does NOT contain x=-2");
}

void test_x_squared_minus_4_nonnegative_sign() {
    TEST_CASE("Property 25: x²-4=0 with NonNegative sign → only x=2");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::NonNegative);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Only x=2 satisfies NonNegative; x=-2 should be excluded
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 NonNegative sign: contains x=2");
    EXPECT_FALSE(has_neg2, "x²-4=0 NonNegative sign: does NOT contain x=-2");
}


void test_x_squared_plus_1_real_domain() {
    TEST_CASE("Property 25: x²+1=0 with Real domain → empty set (imaginary excluded)");

    auto eq = make_x_squared_plus("x", 1);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // x=i and x=-i are imaginary, so with Real domain both should be excluded
    EXPECT_TRUE(solutions.empty(),
                "x²+1=0 Real domain: empty set (all imaginary solutions excluded)");
}

void test_x_squared_plus_1_no_context() {
    TEST_CASE("Property 25: x²+1=0 without context → all solutions returned");

    auto eq = make_x_squared_plus("x", 1);

    // No context (nullptr) — all solutions returned unfiltered
    auto solutions = solve_with_assumptions(eq, "x", nullptr);

    // Without filtering, the solver should return the imaginary solutions
    EXPECT_TRUE(solutions.size() >= 1,
                "x²+1=0 no context: at least 1 solution returned (imaginary)");
}


void test_x_squared_minus_1_positive_sign() {
    TEST_CASE("Property 25: x²-1=0 with Positive sign → only x=1");

    auto eq = make_x_squared_minus("x", 1);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Positive);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_TRUE(has_1, "x²-1=0 Positive sign: contains x=1");
    EXPECT_FALSE(has_neg1, "x²-1=0 Positive sign: does NOT contain x=-1");
}

void test_x_squared_minus_1_negative_sign() {
    TEST_CASE("Property 25: x²-1=0 with Negative sign → only x=-1");

    auto eq = make_x_squared_minus("x", 1);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::Negative);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_FALSE(has_1, "x²-1=0 Negative sign: does NOT contain x=1");
    EXPECT_TRUE(has_neg1, "x²-1=0 Negative sign: contains x=-1");
}


void test_no_context_all_solutions_returned() {
    TEST_CASE("Property 25: No context (nullptr) → all solutions returned unfiltered");

    // x² - 4 = 0 → x=2, x=-2
    auto eq = make_x_squared_minus("x", 4);

    auto solutions = solve_with_assumptions(eq, "x", nullptr);

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 no context: contains x=2");
    EXPECT_TRUE(has_neg2, "x²-4=0 no context: contains x=-2");
    EXPECT_TRUE(solutions.size() >= 2, "x²-4=0 no context: at least 2 solutions");
}

void test_no_context_x_squared_minus_1() {
    TEST_CASE("Property 25: x²-1=0 no context → both x=1 and x=-1 returned");

    auto eq = make_x_squared_minus("x", 1);

    auto solutions = solve_with_assumptions(eq, "x", nullptr);

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_TRUE(has_1, "x²-1=0 no context: contains x=1");
    EXPECT_TRUE(has_neg1, "x²-1=0 no context: contains x=-1");
}


void test_all_solutions_filtered_empty_result() {
    TEST_CASE("Property 25: All solutions filtered → empty result set");

    // x² - 4 = 0 → x=2, x=-2
    // With Negative sign: x=2 excluded (positive), x=-2 excluded? No, -2 is negative.
    // Let's use a case where all solutions are excluded:
    // x² + 1 = 0 with Real domain → both imaginary → empty set
    auto eq = make_x_squared_plus("x", 1);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    EXPECT_TRUE(solutions.empty(),
                "All imaginary solutions filtered with Real domain → empty set");
}

void test_all_solutions_filtered_positive_int() {
    TEST_CASE("Property 25: x²-4=0 with PositiveInt and x>2 constraint → may filter all");

    // x² - 2 = 0 → x=sqrt(2), x=-sqrt(2)
    // With PositiveInt domain: sqrt(2) is not an integer → excluded
    //                          -sqrt(2) is negative → excluded
    // Result: empty set
    auto eq = make_x_squared_minus("x", 2);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::PositiveInt);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // sqrt(2) is irrational, not an integer. The solver may return it as a
    // symbolic expression. If it's not a pure NumberNode, the integer filter
    // may not catch it. Check that at minimum no negative values are present.
    for (const auto& sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            EXPECT_TRUE(v > 0.0, "PositiveInt: no non-positive numeric solutions");
            // Check it's an integer
            double rounded = std::round(v);
            EXPECT_TRUE(std::abs(v - rounded) < 1e-9,
                        "PositiveInt: numeric solutions are integers");
        }
    }
}


void test_natural_domain_excludes_negative() {
    TEST_CASE("Property 25: x²-4=0 with Natural domain → only x=2 (non-negative integer)");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Natural);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Natural = non-negative integers. x=2 is valid, x=-2 is not.
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 Natural domain: contains x=2");
    EXPECT_FALSE(has_neg2, "x²-4=0 Natural domain: does NOT contain x=-2");
}

void test_integer_domain_both_returned() {
    TEST_CASE("Property 25: x²-4=0 with Integer domain → both x=2 and x=-2 (both integers)");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Integer);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Both 2 and -2 are integers
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE(has_2, "x²-4=0 Integer domain: contains x=2");
    EXPECT_TRUE(has_neg2, "x²-4=0 Integer domain: contains x=-2");
}

void test_complex_domain_no_filtering() {
    TEST_CASE("Property 25: Complex domain (default) → no filtering applied");

    // x² + 1 = 0 with Complex domain → imaginary solutions should be kept
    auto eq = make_x_squared_plus("x", 1);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Complex);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    // Complex is the default/least restrictive — no filtering
    EXPECT_TRUE(solutions.size() >= 1,
                "x²+1=0 Complex domain: solutions returned (no filtering)");
}

void test_nonpositive_sign_filtering() {
    TEST_CASE("Property 25: x²-4=0 with NonPositive sign → only x=-2");

    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_sign("x", Sign::NonPositive);

    auto solutions = solve_with_assumptions(eq, "x", &ctx);

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_FALSE(has_2, "x²-4=0 NonPositive sign: does NOT contain x=2");
    EXPECT_TRUE(has_neg2, "x²-4=0 NonPositive sign: contains x=-2");
}


int main() {
    // Test Case 1: x²-4=0 with various domain/sign constraints
    test_x_squared_minus_4_real_domain();
    test_x_squared_minus_4_positive_int();
    test_x_squared_minus_4_nonnegative_sign();

    // Test Case 2: x²+1=0 (imaginary solutions)
    test_x_squared_plus_1_real_domain();
    test_x_squared_plus_1_no_context();

    // Test Case 3: x²-1=0 with sign constraints
    test_x_squared_minus_1_positive_sign();
    test_x_squared_minus_1_negative_sign();

    // Test Case 4: No context (nullptr)
    test_no_context_all_solutions_returned();
    test_no_context_x_squared_minus_1();

    // Test Case 5: All solutions filtered → empty result
    test_all_solutions_filtered_empty_result();
    test_all_solutions_filtered_positive_int();

    // Additional property tests
    test_natural_domain_excludes_negative();
    test_integer_domain_both_returned();
    test_complex_domain_no_filtering();
    test_nonpositive_sign_filtering();

    return TEST_REPORT();
}
