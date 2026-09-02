
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
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(v)));
    return expr;
}

// Helper: create a SymbolicExpr wrapping a NumberNode from Rational
static SymbolicExpr make_rational_expr(int num, int den) {
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(Rational(BigInt(num), BigInt(den))));
    return expr;
}

// Helper: create a SymbolicExpr wrapping a NumberNode from double
static SymbolicExpr make_double_expr(double v) {
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(v)));
    return expr;
}


void test_positive_bigint() {
    TEST_CASE("Positive BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    int values[] = {1, 2, 7, 42, 100};
    for (int v : values) {
        auto expr = make_bigint_expr(v);
        std::string label = "BigInt(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                    label + " query_positive=True");
        EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                    label + " query_negative=False");
        EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                    label + " query_nonnegative=True");
        EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::True,
                    label + " query_integer=True");
        EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                    label + " query_nonzero=True");
    }
}


void test_negative_bigint() {
    TEST_CASE("Negative BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    int values[] = {-1, -2, -7, -42, -100};
    for (int v : values) {
        auto expr = make_bigint_expr(v);
        std::string label = "BigInt(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::False,
                    label + " query_positive=False");
        EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::True,
                    label + " query_negative=True");
        EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::False,
                    label + " query_nonnegative=False");
        EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::True,
                    label + " query_integer=True");
        EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                    label + " query_nonzero=True");
    }
}


void test_zero_bigint() {
    TEST_CASE("Zero BigInt — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_bigint_expr(0);

    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::False,
                "BigInt(0) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                "BigInt(0) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                "BigInt(0) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::True,
                "BigInt(0) query_integer=True");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "BigInt(0) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::False,
                "BigInt(0) query_nonzero=False");
}


void test_positive_rational() {
    TEST_CASE("Positive non-integer Rational — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(3, 2);

    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                "Rational(3/2) query_positive=True");
    EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                "Rational(3/2) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                "Rational(3/2) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::False,
                "Rational(3/2) query_integer=False");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "Rational(3/2) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                "Rational(3/2) query_nonzero=True");
}


void test_negative_rational() {
    TEST_CASE("Negative non-integer Rational — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(-5, 3);

    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::False,
                "Rational(-5/3) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::True,
                "Rational(-5/3) query_negative=True");
    EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::False,
                "Rational(-5/3) query_nonnegative=False");
    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::False,
                "Rational(-5/3) query_integer=False");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "Rational(-5/3) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                "Rational(-5/3) query_nonzero=True");
}


void test_integer_rational() {
    TEST_CASE("Integer Rational (4/1) — query_integer=True");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_rational_expr(4, 1);

    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                "Rational(4/1) query_positive=True");
    EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                "Rational(4/1) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                "Rational(4/1) query_nonnegative=True");
    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::True,
                "Rational(4/1) query_integer=True");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "Rational(4/1) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                "Rational(4/1) query_nonzero=True");
}


void test_positive_double() {
    TEST_CASE("Positive double — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    double values[] = {0.5, 1.0, 3.14, 100.0};
    for (double v : values) {
        auto expr = make_double_expr(v);
        std::string label = "double(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                    label + " query_positive=True");
        EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                    label + " query_negative=False");
        EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                    label + " query_nonnegative=True");
        EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                    label + " query_nonzero=True");
    }
}


void test_negative_double() {
    TEST_CASE("Negative double — sign and domain properties");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    double values[] = {-0.5, -1.0, -3.14, -100.0};
    for (double v : values) {
        auto expr = make_double_expr(v);
        std::string label = "double(" + std::to_string(v) + ")";

        EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::False,
                    label + " query_positive=False");
        EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::True,
                    label + " query_negative=True");
        EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::False,
                    label + " query_nonnegative=False");
        EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                    label + " query_real=True");
        EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::True,
                    label + " query_nonzero=True");
    }
}


void test_zero_double() {
    TEST_CASE("Zero double (0.0) — query_nonzero=False");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(0.0);

    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::False,
                "double(0.0) query_positive=False");
    EXPECT_TRUE(qi.query_negative(expr).value() == Tribool::False,
                "double(0.0) query_negative=False");
    EXPECT_TRUE(qi.query_nonnegative(expr).value() == Tribool::True,
                "double(0.0) query_nonnegative=True");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "double(0.0) query_real=True");
    EXPECT_TRUE(qi.query_nonzero(expr).value() == Tribool::False,
                "double(0.0) query_nonzero=False");
}


void test_non_integer_double() {
    TEST_CASE("Non-integer double (2.5) — query_integer=False");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(2.5);

    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::False,
                "double(2.5) query_integer=False");
    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                "double(2.5) query_positive=True");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "double(2.5) query_real=True");

    // Also test a negative non-integer double
    auto expr2 = make_double_expr(-1.7);

    EXPECT_TRUE(qi.query_integer(expr2).value() == Tribool::False,
                "double(-1.7) query_integer=False");
    EXPECT_TRUE(qi.query_negative(expr2).value() == Tribool::True,
                "double(-1.7) query_negative=True");
    EXPECT_TRUE(qi.query_real(expr2).value() == Tribool::True,
                "double(-1.7) query_real=True");
}


void test_integer_double() {
    TEST_CASE("Integer double (3.0) — query_integer=True");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto expr = make_double_expr(3.0);

    EXPECT_TRUE(qi.query_integer(expr).value() == Tribool::True,
                "double(3.0) query_integer=True");
    EXPECT_TRUE(qi.query_positive(expr).value() == Tribool::True,
                "double(3.0) query_positive=True");
    EXPECT_TRUE(qi.query_real(expr).value() == Tribool::True,
                "double(3.0) query_real=True");
}


// Helper: create a SymbolicExpr from a node
static SymbolicExpr make_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

// Helper: create a VariableNode
static std::shared_ptr<const SymbolicNode> var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

// Helper: create a NumberNode from an integer (shared_ptr)
static std::shared_ptr<const SymbolicNode> num(int v) {
    return lamina::detail::make_node<NumberNode>(BigInt(v));
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


void test_nan_handling() {
    TEST_CASE("NaN rejection");

    bool rejected = false;
    try {
        (void)lamina::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(std::nan("")));
    } catch (const std::invalid_argument& error) {
        rejected =
            std::string(error.what()) == "approximate number must be finite";
    }
    EXPECT_TRUE(rejected,
                "NaN NumberNode construction rejects the invalid AST state");
}


void test_infinity_handling() {
    TEST_CASE("Infinity handling");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    // Positive infinity: FunctionNode with FuncType::Infinity
    auto inf_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    SymbolicExpr pos_inf_expr = make_expr(inf_node);

    EXPECT_TRIBOOL(qi.query_positive(pos_inf_expr).value(), Tribool::True,
                   "+Infinity: query_positive should return True");
    EXPECT_TRIBOOL(qi.query_negative(pos_inf_expr).value(), Tribool::False,
                   "+Infinity: query_negative should return False");
    EXPECT_TRIBOOL(qi.query_integer(pos_inf_expr).value(), Tribool::False,
                   "+Infinity: query_integer should return False");
    EXPECT_TRIBOOL(qi.query_nonzero(pos_inf_expr).value(), Tribool::True,
                   "+Infinity: query_nonzero should return True");

    // Negative infinity: MultiplyNode(-1, Infinity)
    auto neg_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
    auto inf_node2 = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    auto neg_inf_node = lamina::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, inf_node2});
    SymbolicExpr neg_inf_expr = make_expr(neg_inf_node);

    EXPECT_TRIBOOL(qi.query_negative(neg_inf_expr).value(), Tribool::True,
                   "-Infinity: query_negative should return True");
    EXPECT_TRIBOOL(qi.query_positive(neg_inf_expr).value(), Tribool::False,
                   "-Infinity: query_positive should return False");
}


void test_undeclared_variable() {
    TEST_CASE("Undeclared variable");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    SymbolicExpr var_expr = make_expr(var("undeclared"));

    EXPECT_TRIBOOL(qi.query_positive(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(var_expr).value(), Tribool::Unknown,
                   "Undeclared var: query_nonzero should return Unknown");
}


void test_matrix_node() {
    TEST_CASE("MatrixNode returns Unknown");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    std::vector<std::vector<std::shared_ptr<const SymbolicNode>>> grid = {
        {num(1)}
    };
    auto mat_node = lamina::detail::make_node<MatrixNode>(grid);
    SymbolicExpr mat_expr = make_expr(mat_node);

    EXPECT_TRIBOOL(qi.query_positive(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(mat_expr).value(), Tribool::Unknown,
                   "MatrixNode: query_nonzero should return Unknown");
}


void test_relational_node() {
    TEST_CASE("RelationalNode returns Unknown");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto rel_node = lamina::detail::make_node<RelationalNode>(
        var("x"), num(0), RelationalNode::Op::GT);
    SymbolicExpr rel_expr = make_expr(rel_node);

    EXPECT_TRIBOOL(qi.query_positive(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(rel_expr).value(), Tribool::Unknown,
                   "RelationalNode: query_nonzero should return Unknown");
}


void test_logical_node() {
    TEST_CASE("LogicalNode returns Unknown");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto left_rel = lamina::detail::make_node<RelationalNode>(
        var("x"), num(0), RelationalNode::Op::GT);
    auto right_rel = lamina::detail::make_node<RelationalNode>(
        var("y"), num(0), RelationalNode::Op::GT);
    auto logical_node = lamina::detail::make_node<LogicalNode>(
        left_rel, right_rel, LogicalNode::Op::And);
    SymbolicExpr logical_expr = make_expr(logical_node);

    EXPECT_TRIBOOL(qi.query_positive(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_positive should return Unknown");
    EXPECT_TRIBOOL(qi.query_negative(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_negative should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonnegative(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_nonnegative should return Unknown");
    EXPECT_TRIBOOL(qi.query_real(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_real should return Unknown");
    EXPECT_TRIBOOL(qi.query_integer(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_integer should return Unknown");
    EXPECT_TRIBOOL(qi.query_nonzero(logical_expr).value(), Tribool::Unknown,
                   "LogicalNode: query_nonzero should return Unknown");
}

void test_checked_query_interface_contracts() {
    TEST_CASE("QueryInterface checked core queries: explicit errors and values");

    AssumptionContext ctx;
    QueryInterface qi(ctx);

    auto positive = make_bigint_expr(3);
    auto positive_result = qi.query_positive_checked(positive);
    EXPECT_TRUE(positive_result.has_value(), "checked query_positive succeeds");
    if (positive_result) {
        EXPECT_TRUE(positive_result.value() == Tribool::True,
                    "checked query_positive returns True for positive integer");
    }

    auto real_result = qi.query_real_checked(positive);
    EXPECT_TRUE(real_result.has_value(), "checked query_real succeeds");
    if (real_result) {
        EXPECT_TRUE(real_result.value() == Tribool::True,
                    "checked query_real returns True for integer");
    }

    auto integer_result = qi.query_integer_checked(positive);
    EXPECT_TRUE(integer_result.has_value(), "checked query_integer succeeds");
    if (integer_result) {
        EXPECT_TRUE(integer_result.value() == Tribool::True,
                    "checked query_integer returns True for integer");
    }

    auto negative = make_bigint_expr(-2);
    auto negative_result = qi.query_negative_checked(negative);
    EXPECT_TRUE(negative_result.has_value(), "checked query_negative succeeds");
    if (negative_result) {
        EXPECT_TRUE(negative_result.value() == Tribool::True,
                    "checked query_negative returns True for negative integer");
    }

    auto zero = make_bigint_expr(0);
    auto nonnegative_result = qi.query_nonnegative_checked(zero);
    EXPECT_TRUE(nonnegative_result.has_value(), "checked query_nonnegative succeeds");
    if (nonnegative_result) {
        EXPECT_TRUE(nonnegative_result.value() == Tribool::True,
                    "checked query_nonnegative returns True for zero");
    }

    auto nonzero_result = qi.query_nonzero_checked(zero);
    EXPECT_TRUE(nonzero_result.has_value(), "checked query_nonzero succeeds");
    if (nonzero_result) {
        EXPECT_TRUE(nonzero_result.value() == Tribool::False,
                    "checked query_nonzero returns False for zero");
    }

    auto rel_node = lamina::detail::make_node<RelationalNode>(
        var("x"), num(0), RelationalNode::Op::GT);
    SymbolicExpr rel_expr = make_expr(rel_node);
    auto relation_result = qi.query_real_checked(rel_expr);
    EXPECT_TRUE(relation_result.has_value(),
                "checked query_real accepts handled compatibility expression types");
    if (relation_result) {
        EXPECT_TRUE(relation_result.value() == Tribool::Unknown,
                    "checked query_real preserves Unknown for relational expressions");
    }

    auto relation_positive = qi.query_positive_checked(rel_expr);
    EXPECT_TRUE(relation_positive.has_value(),
                "checked query_positive accepts relational compatibility expressions");
    if (relation_positive) {
        EXPECT_TRUE(relation_positive.value() == Tribool::Unknown,
                    "checked query_positive reports Unknown for relational expressions");
    }
}


int main() {
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

    test_nan_handling();
    test_infinity_handling();
    test_undeclared_variable();
    test_matrix_node();
    test_relational_node();
    test_logical_node();
    test_checked_query_interface_contracts();

    return TEST_REPORT();
}
