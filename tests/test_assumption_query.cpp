/**
 * @file test_assumption_query.cpp
 * @brief Property tests for QueryInterface (Property 21).
 *
 * Feature: assumption-system, Property 21: NumberNode direct evaluation
 * Validates: Requirements 10.7
 *
 * For any finite numeric value (BigInt, Rational, or floating-point),
 * the QueryInterface should determine sign and domain properties directly
 * from the numeric value without consulting the PropertyStore, and the
 * results should be mathematically correct.
 */

#include "test_common.hpp"
#include "query_interface.hpp"
#include "assumption_context.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <cmath>

using namespace lamina;

// Helper: create a SymbolicExpr wrapping a NumberNode from BigInt
static SymbolicExpr make_bigint_expr(int v) {
    SymbolicExpr expr;
    expr.root = std::make_shared<NumberNode>(BigInt(v));
    return expr;
}

// Helper: create a SymbolicExpr wrapping a NumberNode from Rational
static SymbolicExpr make_rational_expr(int num, int den) {
    SymbolicExpr expr;
    expr.root = std::make_shared<NumberNode>(Rational(BigInt(num), BigInt(den)));
    return expr;
}

// Helper: create a SymbolicExpr wrapping a NumberNode from double
static SymbolicExpr make_double_expr(double v) {
    SymbolicExpr expr;
    expr.root = std::make_shared<NumberNode>(static_cast<lmmc_real_t>(v));
    return expr;
}

// ============================================================
// Test: Positive integers (BigInt)
// query_positive=True, query_negative=False, query_nonnegative=True,
// query_integer=True, query_real=True, query_nonzero=True
// ============================================================

void test_positive_bigint() {
    TEST_CASE("Property 21: Positive BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    int values[] = {1, 2, 7, 42, 100};
    for (int v : values) {
        auto expr = make_bigint_expr(v);
        std::string label = "BigInt(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                    label + " query_positive=True");
        EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                    label + " query_negative=False");
        EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                    label + " query_nonnegative=True");
        EXPECT_TRUE(qi.query_integer(expr) == Tribool::True,
                    label + " query_integer=True");
        EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                    label + " query_nonzero=True");
    }
}

// ============================================================
// Test: Negative integers (BigInt)
// query_positive=False, query_negative=True, query_nonnegative=False,
// query_integer=True, query_real=True, query_nonzero=True
// ============================================================

void test_negative_bigint() {
    TEST_CASE("Property 21: Negative BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    int values[] = {-1, -2, -7, -42, -100};
    for (int v : values) {
        auto expr = make_bigint_expr(v);
        std::string label = "BigInt(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr) == Tribool::False,
                    label + " query_positive=False");
        EXPECT_TRUE(qi.query_negative(expr) == Tribool::True,
                    label + " query_negative=True");
        EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::False,
                    label + " query_nonnegative=False");
        EXPECT_TRUE(qi.query_integer(expr) == Tribool::True,
                    label + " query_integer=True");
        EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                    label + " query_nonzero=True");
    }
}

// ============================================================
// Test: Zero (BigInt(0))
// query_positive=False, query_negative=False, query_nonnegative=True,
// query_integer=True, query_real=True, query_nonzero=False
// ============================================================

void test_zero_bigint() {
    TEST_CASE("Property 21: Zero BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_bigint_expr(0);

    EXPECT_TRUE(qi.query_positive(expr) == Tribool::False,
                "BigInt(0) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                "BigInt(0) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                "BigInt(0) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr) == Tribool::True,
                "BigInt(0) query_integer=True");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "BigInt(0) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::False,
                "BigInt(0) query_nonzero=False");
}

// ============================================================
// Test: Positive Rational (e.g., 3/2)
// query_positive=True, query_integer=False (non-integer rational),
// query_real=True
// ============================================================

void test_positive_rational() {
    TEST_CASE("Property 21: Positive non-integer Rational — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(3, 2);

    EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                "Rational(3/2) query_positive=True");
    EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                "Rational(3/2) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                "Rational(3/2) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr) == Tribool::False,
                "Rational(3/2) query_integer=False");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "Rational(3/2) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                "Rational(3/2) query_nonzero=True");
}

// ============================================================
// Test: Negative Rational
// query_negative=True, query_integer=False, query_real=True
// ============================================================

void test_negative_rational() {
    TEST_CASE("Property 21: Negative non-integer Rational — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(-5, 3);

    EXPECT_TRUE(qi.query_positive(expr) == Tribool::False,
                "Rational(-5/3) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr) == Tribool::True,
                "Rational(-5/3) query_negative=True");
    EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::False,
                "Rational(-5/3) query_nonnegative=False");
    EXPECT_TRUE(qi.query_integer(expr) == Tribool::False,
                "Rational(-5/3) query_integer=False");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "Rational(-5/3) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                "Rational(-5/3) query_nonzero=True");
}

// ============================================================
// Test: Integer Rational (e.g., 4/1)
// query_integer=True
// ============================================================

void test_integer_rational() {
    TEST_CASE("Property 21: Integer Rational (4/1) — query_integer=True");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(4, 1);

    EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                "Rational(4/1) query_positive=True");
    EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                "Rational(4/1) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                "Rational(4/1) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr) == Tribool::True,
                "Rational(4/1) query_integer=True");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "Rational(4/1) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                "Rational(4/1) query_nonzero=True");
}

// ============================================================
// Test: Positive double
// query_positive=True, query_real=True
// ============================================================

void test_positive_double() {
    TEST_CASE("Property 21: Positive double — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    double values[] = {0.5, 1.0, 3.14, 100.0};
    for (double v : values) {
        auto expr = make_double_expr(v);
        std::string label = "double(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                    label + " query_positive=True");
        EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                    label + " query_negative=False");
        EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                    label + " query_nonnegative=True");
        EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                    label + " query_nonzero=True");
    }
}

// ============================================================
// Test: Negative double
// query_negative=True, query_real=True
// ============================================================

void test_negative_double() {
    TEST_CASE("Property 21: Negative double — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    double values[] = {-0.5, -1.0, -3.14, -100.0};
    for (double v : values) {
        auto expr = make_double_expr(v);
        std::string label = "double(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr) == Tribool::False,
                    label + " query_positive=False");
        EXPECT_TRUE(qi.query_negative(expr) == Tribool::True,
                    label + " query_negative=True");
        EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::False,
                    label + " query_nonnegative=False");
        EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::True,
                    label + " query_nonzero=True");
    }
}

// ============================================================
// Test: Zero double (0.0)
// query_nonzero=False
// ============================================================

void test_zero_double() {
    TEST_CASE("Property 21: Zero double (0.0) — query_nonzero=False");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(0.0);

    EXPECT_TRUE(qi.query_positive(expr) == Tribool::False,
                "double(0.0) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr) == Tribool::False,
                "double(0.0) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr) == Tribool::True,
                "double(0.0) query_nonnegative=True");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "double(0.0) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr) == Tribool::False,
                "double(0.0) query_nonzero=False");
}

// ============================================================
// Test: Non-integer double (e.g., 2.5)
// query_integer=False
// ============================================================

void test_non_integer_double() {
    TEST_CASE("Property 21: Non-integer double (2.5) — query_integer=False");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(2.5);

    EXPECT_TRUE(qi.query_integer(expr) == Tribool::False,
                "double(2.5) query_integer=False");
    EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                "double(2.5) query_positive=True");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "double(2.5) query_real=True");

    // Also test a negative non-integer double
    auto expr2 = make_double_expr(-1.7);

    EXPECT_TRUE(qi.query_integer(expr2) == Tribool::False,
                "double(-1.7) query_integer=False");
    EXPECT_TRUE(qi.query_negative(expr2) == Tribool::True,
                "double(-1.7) query_negative=True");
    EXPECT_TRUE(qi.query_real(expr2) == Tribool::True,
                "double(-1.7) query_real=True");
}

// ============================================================
// Additional: Integer double (e.g., 3.0) should report query_integer=True
// ============================================================

void test_integer_double() {
    TEST_CASE("Property 21: Integer double (3.0) — query_integer=True");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(3.0);

    EXPECT_TRUE(qi.query_integer(expr) == Tribool::True,
                "double(3.0) query_integer=True");
    EXPECT_TRUE(qi.query_positive(expr) == Tribool::True,
                "double(3.0) query_positive=True");
    EXPECT_TRUE(qi.query_real(expr) == Tribool::True,
                "double(3.0) query_real=True");
}

// ============================================================
// Edge Case Tests (Task 8.3)
// ============================================================

// Helper: create a SymbolicExpr from a node
static SymbolicExpr make_expr(std::shared_ptr<SymbolicNode> node) {
    SymbolicExpr expr;
    expr.root = std::move(node);
    return expr;
}

// Helper: create a VariableNode
static std::shared_ptr<SymbolicNode> var(const std::string& name) {
    return std::make_shared<VariableNode>(name);
}

// Helper: create a NumberNode from an integer (shared_ptr)
static std::shared_ptr<SymbolicNode> num(int v) {
    return std::make_shared<NumberNode>(BigInt(v));
}

// Helper: check Tribool equality with descriptive output
static void EXPECT_TRIBOOL(Tribool actual, Tribool expected, const std::string& msg) {
    const char* names[] = {"True", "False", "Unknown"};
    int ai = static_cast<int>(actual);
    int ei = static_cast<int>(expected);
    if (ai == ei) {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    } else {
        std::cerr << "[FAIL] " << msg
                  << "\n  Expected: " << names[ei]
                  << "\n  Got:      " << names[ai] << std::endl;
        g_failures++;
    }
}

// ============================================================
// Test: NaN handling (Req 10.11)
// query_integer → False
// query_positive, query_negative, query_nonnegative, query_nonzero → Unknown
// ============================================================

void test_nan_handling() {
    TEST_CASE("NaN handling (Req 10.11)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto nan_node = std::make_shared<NumberNode>(static_cast<lmmc_real_t>(std::nan("")));
    SymbolicExpr nan_expr = make_expr(nan_node);

    EXPECT_TRIBOOL(qi.query_integer(nan_expr), Tribool::False,
                   "NaN: query_integer should return False");
    EXPECT_TRIBOOL(qi.query_positive(nan_expr), Tribool::Unknown,
                   "NaN: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(nan_expr), Tribool::Unknown,
                   "NaN: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(nan_expr), Tribool::Unknown,
                   "NaN: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(nan_expr), Tribool::Unknown,
                   "NaN: query_nonzero should return Unknown");
}

// ============================================================
// Test: Infinity handling (Req 10.11)
// Positive infinity: query_positive=True, query_negative=False,
//   query_integer=False, query_nonzero=True
// Negative infinity: query_negative=True, query_positive=False
// ============================================================

void test_infinity_handling() {
    TEST_CASE("Infinity handling (Req 10.11)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    // Positive infinity: FunctionNode with FuncType::Infinity
    auto inf_node = std::make_shared<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<SymbolicNode>>{});
    SymbolicExpr pos_inf_expr = make_expr(inf_node);

    EXPECT_TRIBOOL(qi.query_positive(pos_inf_expr), Tribool::True,
                   "+Infinity: query_positive should return True");
    EXPECT_TRIBOOL(qi.query_negative(pos_inf_expr), Tribool::False,
                   "+Infinity: query_negative should return False");
    EXPECT_TRIBOOL(qi.query_integer(pos_inf_expr), Tribool::False,
                   "+Infinity: query_integer should return False");
    EXPECT_TRIBOOL(qi.query_nonzero(pos_inf_expr), Tribool::True,
                   "+Infinity: query_nonzero should return True");

    // Negative infinity: MultiplyNode(-1, Infinity)
    auto neg_one = std::make_shared<NumberNode>(BigInt(-1));
    auto inf_node2 = std::make_shared<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<SymbolicNode>>{});
    auto neg_inf_node = std::make_shared<MultiplyNode>(
        std::vector<std::shared_ptr<SymbolicNode>>{neg_one, inf_node2});
    SymbolicExpr neg_inf_expr = make_expr(neg_inf_node);

    EXPECT_TRIBOOL(qi.query_negative(neg_inf_expr), Tribool::True,
                   "-Infinity: query_negative should return True");
    EXPECT_TRIBOOL(qi.query_positive(neg_inf_expr), Tribool::False,
                   "-Infinity: query_positive should return False");
}

// ============================================================
// Test: Null root node (Req 10.10)
// All queries → Unknown
// ============================================================

void test_null_root_node() {
    TEST_CASE("Null root node (Req 10.10)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    SymbolicExpr null_expr;

    EXPECT_TRIBOOL(qi.query_positive(null_expr), Tribool::Unknown,
                   "Null root: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(null_expr), Tribool::Unknown,
                   "Null root: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(null_expr), Tribool::Unknown,
                   "Null root: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(null_expr), Tribool::Unknown,
                   "Null root: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(null_expr), Tribool::Unknown,
                   "Null root: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(null_expr), Tribool::Unknown,
                   "Null root: query_nonzero should return Unknown");
}

// ============================================================
// Test: Undeclared variable (Req 10.8)
// All queries → Unknown
// ============================================================

void test_undeclared_variable() {
    TEST_CASE("Undeclared variable (Req 10.8)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    SymbolicExpr var_expr = make_expr(var("undeclared"));

    EXPECT_TRIBOOL(qi.query_positive(var_expr), Tribool::Unknown,
                   "Undeclared var: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(var_expr), Tribool::Unknown,
                   "Undeclared var: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(var_expr), Tribool::Unknown,
                   "Undeclared var: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(var_expr), Tribool::Unknown,
                   "Undeclared var: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(var_expr), Tribool::Unknown,
                   "Undeclared var: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(var_expr), Tribool::Unknown,
                   "Undeclared var: query_nonzero should return Unknown");
}

// ============================================================
// Test: MatrixNode returns Unknown (Req 10.10)
// ============================================================

void test_matrix_node() {
    TEST_CASE("MatrixNode returns Unknown (Req 10.10)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    std::vector<std::vector<std::shared_ptr<SymbolicNode>>> grid = {
        {num(1)}
    };
    auto mat_node = std::make_shared<MatrixNode>(grid);
    SymbolicExpr mat_expr = make_expr(mat_node);

    EXPECT_TRIBOOL(qi.query_positive(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(mat_expr), Tribool::Unknown,
                   "MatrixNode: query_nonzero should return Unknown");
}

// ============================================================
// Test: RelationalNode returns Unknown (Req 10.10)
// ============================================================

void test_relational_node() {
    TEST_CASE("RelationalNode returns Unknown (Req 10.10)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto rel_node = std::make_shared<RelationalNode>(
        var("x"), num(0), RelationalNode::Op::GT);
    SymbolicExpr rel_expr = make_expr(rel_node);

    EXPECT_TRIBOOL(qi.query_positive(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(rel_expr), Tribool::Unknown,
                   "RelationalNode: query_nonzero should return Unknown");
}

// ============================================================
// Test: LogicalNode returns Unknown (Req 10.10)
// ============================================================

void test_logical_node() {
    TEST_CASE("LogicalNode returns Unknown (Req 10.10)");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto left_rel = std::make_shared<RelationalNode>(
        var("x"), num(0), RelationalNode::Op::GT);
    auto right_rel = std::make_shared<RelationalNode>(
        var("y"), num(0), RelationalNode::Op::GT);
    auto logical_node = std::make_shared<LogicalNode>(
        left_rel, right_rel, LogicalNode::Op::And);
    SymbolicExpr logical_expr = make_expr(logical_node);

    EXPECT_TRIBOOL(qi.query_positive(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(logical_expr), Tribool::Unknown,
                   "LogicalNode: query_nonzero should return Unknown");
}

// ============================================================
// main
// ============================================================

int main() {
    // Task 8.2: Property 21 — NumberNode direct evaluation
    test_positive_bigint();
    test_negative_bigint();
    test_zero_bigint();
    test_positive_rational();
    test_negative_rational();
    test_integer_rational();
    test_positive_double();
    test_negative_double();
    test_zero_double();
    test_non_integer_double();
    test_integer_double();

    // Task 8.3: Edge cases
    test_nan_handling();
    test_infinity_handling();
    test_null_root_node();
    test_undeclared_variable();
    test_matrix_node();
    test_relational_node();
    test_logical_node();

    return TEST_REPORT();
}
