
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "symbolic_ast.hpp"
#include <memory>
#include <vector>
#include <string>
#include <random>

using namespace lamina;


/// Create a VariableNode
static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

/// Create a NumberNode from an integer
static std::shared_ptr<const SymbolicNode> make_num(int v) {
    return lamina::detail::make_node<NumberNode>(BigInt(v));
}

/// Create a NumberNode from a double
static std::shared_ptr<const SymbolicNode> make_num_d(double v) {
    return lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(v));
}

/// Create a PowerNode expression
static SymbolicExpr make_power_expr(std::shared_ptr<const SymbolicNode> base,
                                    std::shared_ptr<const SymbolicNode> exponent) {
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent)));
    return expr;
}


void test_positive_base_real_exponent() {
    TEST_CASE("Positive base + Real exponent -> Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    // x^y where x is Positive, y is Real
    auto expr = make_power_expr(make_var("x"), make_var("y"));

    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "x^y is Positive when x>0 and y is Real");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "x^y is NonNegative when x>0 and y is Real");
    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
                "x^y is NonZero when x>0 and y is Real");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False,
                "x^y is not Negative when x>0 and y is Real");
}

void test_positive_base_integer_exponent() {
    TEST_CASE("Positive base + integer exponent -> Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // Test with multiple integer exponents (integers are Real)
    std::vector<int> exponents = {1, 2, 3, 5, 10, -1, -2, -3};
    for (int exp : exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Positive when x>0";
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True, msg);
    }
}

void test_positive_base_real_number_exponent() {
    TEST_CASE("Positive base + real number exponent -> Positive");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);

    InferenceEngine engine(ctx);

    // Real number exponents (non-integer)
    std::vector<double> exponents = {0.5, 1.5, 2.7, -0.5, -1.5, 3.14};
    for (double exp : exponents) {
        auto expr = make_power_expr(make_var("x"), make_num_d(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Positive when x>0";
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True, msg);
    }
}

void test_numeric_positive_base() {
    TEST_CASE("Numeric positive base + variable Real exponent");

    AssumptionContext ctx;
    ctx.assume_domain("y", Domain::Real);

    InferenceEngine engine(ctx);

    // 2^y, 5^y, 100^y where y is Real
    std::vector<int> bases = {1, 2, 5, 10, 100};
    for (int b : bases) {
        auto expr = make_power_expr(make_num(b), make_var("y"));
        std::string msg = std::to_string(b) + "^y is Positive when y is Real";
        EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True, msg);
    }
}


void test_real_base_even_exponent() {
    TEST_CASE("Real base + even integer exponent -> NonNegative");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Test with multiple even exponents
    std::vector<int> even_exponents = {2, 4, 6, 8, 10, 20, 100};
    for (int exp : even_exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonNegative when x is Real (even exponent)";
        EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True, msg);
        std::string msg2 = "x^" + std::to_string(exp) +
                           " is not Negative when x is Real (even exponent)";
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False, msg2);
    }
}

void test_real_base_odd_exponent_not_nonneg() {
    TEST_CASE("Real base + odd exponent -> NOT necessarily NonNegative");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Odd exponents: x^3, x^5 — Real base with odd exponent can be negative
    // so NonNegative should be Unknown (not True)
    std::vector<int> odd_exponents = {1, 3, 5, 7};
    for (int exp : odd_exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Unknown for NonNegative when x is Real (odd exponent)";
        EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown, msg);
    }
}

void test_randomized_even_exponents() {
    TEST_CASE("Randomized even exponents (100 iterations)");

    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(1, 50);
    int pass_count = 0;
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int half = dist(rng);
        int even_exp = half * 2;  // Always even

        AssumptionContext ctx;
        std::string var_name = "x" + std::to_string(i);
        ctx.assume_domain(var_name, Domain::Real);

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(make_var(var_name), make_num(even_exp));

        if (engine.query_nonnegative_checked(expr).value() == Tribool::True) {
            pass_count++;
        }
    }

    std::string msg = std::to_string(pass_count) +
                      "/" + std::to_string(NUM_ITERATIONS) +
                      " random even exponents yield NonNegative";
    EXPECT_TRUE(pass_count == NUM_ITERATIONS, msg);
}


void test_nonneg_base_positive_int_exponent() {
    TEST_CASE("NonNegative base + positive integer exponent -> NonNegative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonNegative);

    InferenceEngine engine(ctx);

    // Test with multiple positive integer exponents
    std::vector<int> pos_exponents = {1, 2, 3, 4, 5, 10, 50};
    for (int exp : pos_exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonNegative when x>=0 (positive int exponent)";
        EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True, msg);
        std::string msg2 = "x^" + std::to_string(exp) +
                           " is not Negative when x>=0";
        EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::False, msg2);
    }
}

void test_randomized_positive_exponents() {
    TEST_CASE("Randomized positive exponents (100 iterations)");

    std::mt19937 rng(123);
    std::uniform_int_distribution<int> dist(1, 100);
    int pass_count = 0;
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int pos_exp = dist(rng);  // Always positive

        AssumptionContext ctx;
        std::string var_name = "v" + std::to_string(i);
        ctx.assume_sign(var_name, Sign::NonNegative);

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(make_var(var_name), make_num(pos_exp));

        if (engine.query_nonnegative_checked(expr).value() == Tribool::True) {
            pass_count++;
        }
    }

    std::string msg = std::to_string(pass_count) +
                      "/" + std::to_string(NUM_ITERATIONS) +
                      " random positive exponents yield NonNegative";
    EXPECT_TRUE(pass_count == NUM_ITERATIONS, msg);
}


void test_nonzero_base_integer_exponent() {
    TEST_CASE("NonZero base + integer exponent -> NonZero");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::NonZero);

    InferenceEngine engine(ctx);

    // Test with various integer exponents (positive, negative, even, odd)
    std::vector<int> exponents = {1, 2, 3, -1, -2, -3, 5, 10, -10};
    for (int exp : exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is NonZero when x!=0 (integer exponent)";
        EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True, msg);
    }
}

void test_nonzero_randomized_integer_exponents() {
    TEST_CASE("Randomized integer exponents (100 iterations)");

    std::mt19937 rng(456);
    std::uniform_int_distribution<int> dist(-50, 50);
    int pass_count = 0;
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int exp = dist(rng);

        AssumptionContext ctx;
        std::string var_name = "z" + std::to_string(i);
        ctx.assume_sign(var_name, Sign::NonZero);

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(make_var(var_name), make_num(exp));

        if (engine.query_nonzero_checked(expr).value() == Tribool::True) {
            pass_count++;
        }
    }

    std::string msg = std::to_string(pass_count) +
                      "/" + std::to_string(NUM_ITERATIONS) +
                      " random integer exponents yield NonZero";
    EXPECT_TRUE(pass_count == NUM_ITERATIONS, msg);
}


void test_positive_base_even_exponent_is_positive() {
    TEST_CASE("Positive base + even exponent -> Positive (via 16a)");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // x^2 where x is Positive and Real → Positive (rule 16a fires)
    auto expr = make_power_expr(make_var("x"), make_num(2));
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "x^2 is Positive when x>0 (rule 16a)");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::True,
                "x^2 is NonNegative when x>0");
}

void test_numeric_base_and_exponent() {
    TEST_CASE("Numeric positive base + integer exponent");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // 2^3 — NumberNode base (positive), NumberNode exponent (integer)
    auto expr = make_power_expr(make_num(2), make_num(3));
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::True,
                "2^3 is Positive");
    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::True,
                "2^3 is NonZero");

    // 3^(-2) — positive base, integer exponent
    auto expr2 = make_power_expr(make_num(3), make_num(-2));
    EXPECT_TRUE(engine.query_positive_checked(expr2).value() == Tribool::True,
                "3^(-2) is Positive");
}

void test_no_rule_matches() {
    TEST_CASE("No rule matches -> Unknown");

    AssumptionContext ctx;
    // No assumptions about x or y
    InferenceEngine engine(ctx);

    auto expr = make_power_expr(make_var("x"), make_var("y"));
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "x^y is Unknown when no assumptions");
    EXPECT_TRUE(engine.query_negative_checked(expr).value() == Tribool::Unknown,
                "x^y is Unknown for Negative when no assumptions");
    EXPECT_TRUE(engine.query_nonnegative_checked(expr).value() == Tribool::Unknown,
                "x^y is Unknown for NonNegative when no assumptions");
    EXPECT_TRUE(engine.query_nonzero_checked(expr).value() == Tribool::Unknown,
                "x^y is Unknown for NonZero when no assumptions");
}

void test_negative_base_non_integer_exponent() {
    TEST_CASE("Negative base + non-integer exponent -> Unknown");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Negative);

    InferenceEngine engine(ctx);

    // x^(1.5) where x is Negative — could be complex, no rule matches
    auto expr = make_power_expr(make_var("x"), make_num_d(1.5));
    EXPECT_TRUE(engine.query_positive_checked(expr).value() == Tribool::Unknown,
                "(-x)^1.5 is Unknown for Positive");
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "(-x)^1.5 is Unknown for Real");
}


void test_real_base_integer_exponent_domain() {
    TEST_CASE("Real base + integer exponent -> Real domain");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // Test with multiple integer exponents
    std::vector<int> exponents = {0, 1, 2, 3, -1, -2, -3, 5, 10, -10};
    for (int exp : exponents) {
        auto expr = make_power_expr(make_var("x"), make_num(exp));
        std::string msg = "x^" + std::to_string(exp) +
                          " is Real when x is Real (integer exponent)";
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, msg);
    }
}

void test_domain_randomized_integer_exponents() {
    TEST_CASE("Randomized integer exponents (100 iterations)");

    std::mt19937 rng(789);
    std::uniform_int_distribution<int> dist(-50, 50);
    int pass_count = 0;
    const int NUM_ITERATIONS = 100;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        int exp = dist(rng);

        AssumptionContext ctx;
        std::string var_name = "r" + std::to_string(i);
        ctx.assume_domain(var_name, Domain::Real);

        InferenceEngine engine(ctx);
        auto expr = make_power_expr(make_var(var_name), make_num(exp));

        if (engine.query_real_checked(expr).value() == Tribool::True) {
            pass_count++;
        }
    }

    std::string msg = std::to_string(pass_count) +
                      "/" + std::to_string(NUM_ITERATIONS) +
                      " random integer exponents yield Real domain";
    EXPECT_TRUE(pass_count == NUM_ITERATIONS, msg);
}

void test_integer_base_integer_exponent() {
    TEST_CASE("Integer base + integer exponent -> Real (Integer implies Real)");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);

    InferenceEngine engine(ctx);

    // Integer is a subset of Real, so this should also yield Real
    std::vector<int> exponents = {1, 2, 3, -1, 5};
    for (int exp : exponents) {
        auto expr = make_power_expr(make_var("n"), make_num(exp));
        std::string msg = "n^" + std::to_string(exp) +
                          " is Real when n is Integer (Integer implies Real)";
        EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::True, msg);
    }
}

void test_numeric_base_integer_exponent() {
    TEST_CASE("Numeric base + integer exponent -> Real");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    // Numeric bases are Real, integer exponents are integers
    auto expr1 = make_power_expr(make_num(2), make_num(3));
    EXPECT_TRUE(engine.query_real_checked(expr1).value() == Tribool::True,
                "2^3 is Real");

    auto expr2 = make_power_expr(make_num(-3), make_num(2));
    EXPECT_TRUE(engine.query_real_checked(expr2).value() == Tribool::True,
                "(-3)^2 is Real");

    auto expr3 = make_power_expr(make_num(5), make_num(-1));
    EXPECT_TRUE(engine.query_real_checked(expr3).value() == Tribool::True,
                "5^(-1) is Real");
}

void test_real_base_non_integer_exponent_unknown() {
    TEST_CASE("Real base + non-integer exponent -> Unknown for Real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);

    InferenceEngine engine(ctx);

    // x^(1.5) where x is Real — non-integer exponent, rule doesn't apply
    // (could be complex if x < 0)
    auto expr = make_power_expr(make_var("x"), make_num_d(1.5));
    /// 该域推导规则要求整数指数；仅有 Real 基底时结果域保持 Unknown。
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "x^1.5 is Unknown for Real when x is only Real (non-integer exponent)");
}

void test_no_domain_yields_unknown() {
    TEST_CASE("No domain assumption -> Unknown for Real");

    AssumptionContext ctx;
    // No assumptions about x
    InferenceEngine engine(ctx);

    auto expr = make_power_expr(make_var("x"), make_num(2));
    EXPECT_TRUE(engine.query_real_checked(expr).value() == Tribool::Unknown,
                "x^2 is Unknown for Real when x has no domain assumption");
}


int main() {
    test_positive_base_real_exponent();
    test_positive_base_integer_exponent();
    test_positive_base_real_number_exponent();
    test_numeric_positive_base();
    test_real_base_even_exponent();
    test_real_base_odd_exponent_not_nonneg();
    test_randomized_even_exponents();
    test_nonneg_base_positive_int_exponent();
    test_randomized_positive_exponents();
    test_nonzero_base_integer_exponent();
    test_nonzero_randomized_integer_exponents();
    test_positive_base_even_exponent_is_positive();
    test_numeric_base_and_exponent();
    test_no_rule_matches();
    test_negative_base_non_integer_exponent();

    test_real_base_integer_exponent_domain();
    test_domain_randomized_integer_exponents();
    test_integer_base_integer_exponent();
    test_numeric_base_integer_exponent();
    test_real_base_non_integer_exponent_unknown();
    test_no_domain_yields_unknown();

    return TEST_REPORT();
}
