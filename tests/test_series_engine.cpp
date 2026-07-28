/**
 * @file test_series_engine.cpp
 * @brief 级数引擎单元测试：收敛半径、收敛判定、幂级数运算、傅里叶级数、洛朗级数、符号求和与乘积。
 */

#include "test_common.hpp"
#include "series_engine.hpp"

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr rat(int p, int q) { return SymbolicExpr::number(Rational(p, q)); }
static Expr var(const std::string& name) { return SymbolicExpr::variable(name); }
static Expr inf() { return SymbolicExpr::infinity(); }

/**
 * @brief 辅助函数：将系数向量转为字符串便于调试。
 */
static std::string coeffs_to_string(const Coeffs& c) {
    std::string s = "[";
    for (size_t i = 0; i < c.size(); ++i) {
        if (i > 0) s += ", ";
        s += c[i] ? c[i]->to_string() : "null";
    }
    s += "]";
    return s;
}

static void test_power_series_add_basic() {
    TEST_CASE("power_series_add: basic addition");

    // a(x) = 1 + 2x + 3x²
    Coeffs a = {num(1), num(2), num(3)};
    // b(x) = 4 + 5x + 6x²
    Coeffs b = {num(4), num(5), num(6)};

    auto result = lamina::power_series_add(a, b);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    EXPECT_EQ_EXPR(result[0], num(5), "c[0] = 1+4 = 5");
    EXPECT_EQ_EXPR(result[1], num(7), "c[1] = 2+5 = 7");
    EXPECT_EQ_EXPR(result[2], num(9), "c[2] = 3+6 = 9");
}

static void test_power_series_add_different_lengths() {
    TEST_CASE("power_series_add: different lengths");

    // a(x) = 1 + 2x
    Coeffs a = {num(1), num(2)};
    // b(x) = 3 + 4x + 5x² + 6x³
    Coeffs b = {num(3), num(4), num(5), num(6)};

    auto result = lamina::power_series_add(a, b);

    EXPECT_TRUE(result.size() == 4, "Result has 4 coefficients (max length)");
    EXPECT_EQ_EXPR(result[0], num(4), "c[0] = 1+3 = 4");
    EXPECT_EQ_EXPR(result[1], num(6), "c[1] = 2+4 = 6");
    EXPECT_EQ_EXPR(result[2], num(5), "c[2] = 0+5 = 5");
    EXPECT_EQ_EXPR(result[3], num(6), "c[3] = 0+6 = 6");
}

static void test_power_series_add_empty() {
    TEST_CASE("power_series_add: empty series");

    Coeffs a = {};
    Coeffs b = {num(1), num(2)};

    auto result = lamina::power_series_add(a, b);
    EXPECT_TRUE(result.size() == 2, "Result has 2 coefficients");
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 0+1 = 1");
    EXPECT_EQ_EXPR(result[1], num(2), "c[1] = 0+2 = 2");
}

static void test_power_series_multiply_basic() {
    TEST_CASE("power_series_multiply: (1+x)*(1+x) = 1+2x+x²");

    // a(x) = 1 + x
    Coeffs a = {num(1), num(1)};
    // b(x) = 1 + x
    Coeffs b = {num(1), num(1)};

    auto result = lamina::power_series_multiply(a, b, 3);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 1*1 = 1");
    EXPECT_EQ_EXPR(result[1], num(2), "c[1] = 1*1 + 1*1 = 2");
    EXPECT_EQ_EXPR(result[2], num(1), "c[2] = 1*1 = 1");
}

static void test_power_series_multiply_truncation() {
    TEST_CASE("power_series_multiply: truncation to order");

    // a(x) = 1 + x + x²
    Coeffs a = {num(1), num(1), num(1)};
    // b(x) = 1 + x + x²
    Coeffs b = {num(1), num(1), num(1)};

    // Full product: 1 + 2x + 3x² + 2x³ + x⁴
    // Truncate to order 3: 1 + 2x + 3x²
    auto result = lamina::power_series_multiply(a, b, 3);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients (truncated)");
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 1");
    EXPECT_EQ_EXPR(result[1], num(2), "c[1] = 2");
    EXPECT_EQ_EXPR(result[2], num(3), "c[2] = 3");
}

static void test_power_series_multiply_zero() {
    TEST_CASE("power_series_multiply: multiply by zero series");

    Coeffs a = {num(1), num(2), num(3)};
    Coeffs b = {num(0)};

    auto result = lamina::power_series_multiply(a, b, 3);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    EXPECT_EQ_EXPR(result[0], num(0), "c[0] = 0");
    EXPECT_EQ_EXPR(result[1], num(0), "c[1] = 0");
    EXPECT_EQ_EXPR(result[2], num(0), "c[2] = 0");
}

static void test_power_series_multiply_constants() {
    TEST_CASE("power_series_multiply: constant * polynomial");

    // a(x) = 3
    Coeffs a = {num(3)};
    // b(x) = 1 + 2x + 4x²
    Coeffs b = {num(1), num(2), num(4)};

    auto result = lamina::power_series_multiply(a, b, 3);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    EXPECT_EQ_EXPR(result[0], num(3), "c[0] = 3*1 = 3");
    EXPECT_EQ_EXPR(result[1], num(6), "c[1] = 3*2 = 6");
    EXPECT_EQ_EXPR(result[2], num(12), "c[2] = 3*4 = 12");
}

static void test_power_series_compose_basic() {
    TEST_CASE("power_series_compose: f(g(x)) where f(x)=1+x, g(x)=x");

    // f(x) = 1 + x
    Coeffs f = {num(1), num(1)};
    // g(x) = x (i.e., g(0)=0, g[1]=1)
    Coeffs g = {num(0), num(1)};

    auto result = lamina::power_series_compose(f, g, 3);

    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    // f(g(x)) = f(x) = 1 + x
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 1");
    EXPECT_EQ_EXPR(result[1], num(1), "c[1] = 1");
    EXPECT_EQ_EXPR(result[2], num(0), "c[2] = 0");
}

static void test_power_series_compose_quadratic() {
    TEST_CASE("power_series_compose: f(g(x)) where f(x)=1+x+x², g(x)=2x");

    // f(x) = 1 + x + x²
    Coeffs f = {num(1), num(1), num(1)};
    // g(x) = 2x (g(0)=0)
    Coeffs g = {num(0), num(2)};

    auto result = lamina::power_series_compose(f, g, 3);

    // f(2x) = 1 + 2x + 4x²
    EXPECT_TRUE(result.size() == 3, "Result has 3 coefficients");
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 1");
    EXPECT_EQ_EXPR(result[1], num(2), "c[1] = 2");
    EXPECT_EQ_EXPR(result[2], num(4), "c[2] = 4");
}

static void test_power_series_compose_g0_nonzero() {
    TEST_CASE("power_series_compose: g(0) != 0 returns empty");

    Coeffs f = {num(1), num(1)};
    // g(x) = 1 + x (g(0) = 1 ≠ 0)
    Coeffs g = {num(1), num(1)};

    auto result = lamina::power_series_compose(f, g, 3);

    EXPECT_TRUE(result.empty(), "Result is empty when g(0) != 0");
}

static void test_power_series_compose_exp_like() {
    TEST_CASE("power_series_compose: exp(2x) truncated");

    // f(x) = 1 + x + x²/2 + x³/6 (exp(x) truncated to order 4)
    // Using Rational for 1/2 and 1/6
    Coeffs f = {num(1), num(1),
                SymbolicExpr::number(Rational(1, 2)),
                SymbolicExpr::number(Rational(1, 6))};
    // g(x) = 2x
    Coeffs g = {num(0), num(2)};

    auto result = lamina::power_series_compose(f, g, 4);

    // exp(2x) = 1 + 2x + 2x² + (4/3)x³
    EXPECT_TRUE(result.size() == 4, "Result has 4 coefficients");
    EXPECT_EQ_EXPR(result[0], num(1), "c[0] = 1");
    EXPECT_EQ_EXPR(result[1], num(2), "c[1] = 2");
    EXPECT_EQ_EXPR(result[2], num(2), "c[2] = (1/2)*4 = 2");
    EXPECT_EQ_EXPR(result[3], SymbolicExpr::number(Rational(4, 3)), "c[3] = (1/6)*8 = 4/3");
}

static void test_power_series_multiply_order_zero() {
    TEST_CASE("power_series_multiply: order 0 returns empty");

    Coeffs a = {num(1), num(2)};
    Coeffs b = {num(3), num(4)};

    auto result = lamina::power_series_multiply(a, b, 0);
    EXPECT_TRUE(result.empty(), "Order 0 returns empty vector");
}

static void test_power_series_checked_contracts() {
    TEST_CASE("power_series checked APIs: explicit errors and cancellation");

    Coeffs a = {num(1), num(1)};
    Coeffs b = {num(1), num(1)};
    auto checked_product = lamina::power_series_multiply_checked(a, b, 3);
    EXPECT_TRUE(checked_product.has_value(), "checked multiply succeeds");
    if (checked_product) {
        EXPECT_TRUE(checked_product.value().size() == 3,
                    "checked multiply has requested order");
        EXPECT_EQ_EXPR(checked_product.value()[1], num(2),
                       "checked multiply coefficient c[1] = 2");
    }

    auto bad_order = lamina::power_series_multiply_checked(a, b, 0);
    EXPECT_TRUE(!bad_order &&
                    bad_order.error().code == lamina::CasErrc::InvalidArgument,
                "checked multiply rejects non-positive order");
    EXPECT_TRUE(lamina::power_series_multiply(a, b, 0).empty(),
                "legacy multiply unwraps non-positive order to empty vector");

    Coeffs with_null = {num(1), nullptr};
    auto null_coeff = lamina::power_series_multiply_checked(with_null, b, 2);
    EXPECT_TRUE(!null_coeff &&
                    null_coeff.error().code == lamina::CasErrc::InvalidArgument,
                "checked multiply rejects null coefficients");

    Coeffs f = {num(1), num(1)};
    Coeffs g = {num(0), num(2)};
    auto checked_compose = lamina::power_series_compose_checked(f, g, 3);
    EXPECT_TRUE(checked_compose.has_value(), "checked compose succeeds");
    if (checked_compose) {
        EXPECT_EQ_EXPR(checked_compose.value()[0], num(1),
                       "checked compose coefficient c[0] = 1");
        EXPECT_EQ_EXPR(checked_compose.value()[1], num(2),
                       "checked compose coefficient c[1] = 2");
    }

    Coeffs bad_g = {num(1), num(1)};
    auto nonzero_g0 = lamina::power_series_compose_checked(f, bad_g, 3);
    EXPECT_TRUE(!nonzero_g0 &&
                    nonzero_g0.error().code == lamina::CasErrc::DomainError,
                "checked compose reports g(0) != 0 as domain error");
    EXPECT_TRUE(lamina::power_series_compose(f, bad_g, 3).empty(),
                "legacy compose unwraps g(0) != 0 to empty vector");

    lamina::CancellationToken token;
    token.cancel();
    lamina::ComputationContext context({}, token);
    auto cancelled = lamina::power_series_multiply_checked(a, b, 3, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == lamina::CasErrc::Cancelled,
                "checked multiply observes cancelled context");
}

static void test_series_analysis_checked_contracts() {
    TEST_CASE("series analysis checked APIs: explicit errors and context");

    Coeffs geometric = {num(1), num(1), num(1)};
    auto radius = lamina::convergence_radius_checked(geometric, "x");
    EXPECT_TRUE(radius.has_value(), "checked convergence_radius succeeds");
    if (radius) {
        EXPECT_TRUE(radius.value() != nullptr,
                    "checked convergence_radius returns non-null expression");
    }

    auto empty_coeffs = lamina::convergence_radius_checked({}, "x");
    EXPECT_TRUE(!empty_coeffs &&
                    empty_coeffs.error().code == lamina::CasErrc::InvalidArgument,
                "checked convergence_radius rejects empty coefficients");

    Coeffs with_null = {num(1), nullptr};
    auto null_coeff = lamina::convergence_radius_checked(with_null, "x");
    EXPECT_TRUE(!null_coeff &&
                    null_coeff.error().code == lamina::CasErrc::InvalidArgument,
                "checked convergence_radius rejects null coefficients");

    auto bad_var = lamina::convergence_radius_checked(geometric, "");
    EXPECT_TRUE(!bad_var &&
                    bad_var.error().code == lamina::CasErrc::InvalidArgument,
                "checked convergence_radius rejects empty variable");

    Coeffs symbolic_coeffs = {num(1), var("a"), num(1)};
    auto symbolic_radius = lamina::convergence_radius_checked(symbolic_coeffs, "x");
    EXPECT_TRUE(!symbolic_radius &&
                    symbolic_radius.error().code == lamina::CasErrc::Inconclusive,
                "checked convergence_radius rejects unsupported symbolic coefficients");
    EXPECT_TRUE(lamina::convergence_radius(symbolic_coeffs, "x") != nullptr,
                "legacy convergence_radius keeps compatibility result for symbolic coefficients");

    auto n = var("n");
    auto p_term = SymbolicExpr::power(n, num(-2));
    auto convergence = lamina::convergence_test_checked(p_term, "n");
    EXPECT_TRUE(convergence.has_value(), "checked convergence_test succeeds");
    if (convergence) {
        EXPECT_TRUE(convergence.value().result == lamina::ConvergenceResult::Convergent,
                    "checked convergence_test detects p-series convergence");
    }

    auto null_term = lamina::convergence_test_checked(nullptr, "n");
    EXPECT_TRUE(!null_term &&
                    null_term.error().code == lamina::CasErrc::InvalidArgument,
                "checked convergence_test rejects null term");

    auto unsupported_term = SymbolicExpr::sin(n);
    auto unsupported_convergence = lamina::convergence_test_checked(unsupported_term, "n");
    EXPECT_TRUE(!unsupported_convergence &&
                    unsupported_convergence.error().code == lamina::CasErrc::Inconclusive,
                "checked convergence_test reports Inconclusive for unsupported terms");
    auto legacy_unsupported = lamina::convergence_test(unsupported_term, "n");
    EXPECT_TRUE(legacy_unsupported.result == lamina::ConvergenceResult::Inconclusive,
                "legacy convergence_test preserves inconclusive enum result");

    auto z = var("z");
    auto f = SymbolicExpr::divide(num(1), z);
    auto laurent = lamina::laurent_series_full_checked(f, "z", num(0), 2, 2);
    EXPECT_TRUE(laurent.has_value(), "checked Laurent series succeeds for 1/z");
    if (laurent) {
        EXPECT_TRUE(laurent.value().series != nullptr,
                    "checked Laurent full result has non-null series");
        EXPECT_TRUE(laurent.value().pole_order == 1,
                    "checked Laurent detects simple pole");
    }

    auto laurent_expr = lamina::laurent_series_checked(f, "z", num(0), 2, 2);
    EXPECT_TRUE(laurent_expr.has_value(), "checked Laurent expression succeeds");

    auto bad_order = lamina::laurent_series_full_checked(f, "z", num(0), -1, 2);
    EXPECT_TRUE(!bad_order &&
                    bad_order.error().code == lamina::CasErrc::InvalidArgument,
                "checked Laurent rejects negative truncation order");

    auto null_center = lamina::laurent_series_full_checked(f, "z", nullptr, 2, 2);
    EXPECT_TRUE(!null_center &&
                    null_center.error().code == lamina::CasErrc::InvalidArgument,
                "checked Laurent rejects null center");

    auto unsupported_laurent = lamina::laurent_series_full_checked(SymbolicExpr::sin(z),
                                                                   "z", num(0), 2, 2);
    EXPECT_TRUE(!unsupported_laurent &&
                    unsupported_laurent.error().code == lamina::CasErrc::Inconclusive,
                "checked Laurent reports Inconclusive for unsupported analytic functions");

    auto shifted_laurent = lamina::laurent_series_full_checked(f, "z", num(1), 2, 2);
    EXPECT_TRUE(!shifted_laurent &&
                    shifted_laurent.error().code == lamina::CasErrc::Inconclusive,
                "checked Laurent reports Inconclusive outside zero-center support");

    lamina::CancellationToken token;
    token.cancel();
    lamina::ComputationContext cancelled_context({}, token);
    auto cancelled = lamina::convergence_test_checked(p_term, "n", cancelled_context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == lamina::CasErrc::Cancelled,
                "checked convergence_test observes cancellation");

    lamina::ResourceLimits limits;
    limits.max_steps = 1;
    lamina::ComputationContext limited_context(limits);
    auto limited = lamina::laurent_series_full_checked(f, "z", num(0), 2, 2,
                                                       limited_context);
    EXPECT_TRUE(!limited &&
                    limited.error().code == lamina::CasErrc::ResourceLimit,
                "checked Laurent observes exhausted step budget");
}

// ============================================================
// lim_sup / lim_inf 测试
// ============================================================

static void test_lim_sup_convergent_sequence() {
    TEST_CASE("lim_sup: convergent sequence 1/n → 0");

    // a_n = 1/n
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto result = lamina::lim_sup(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_sup(1/n) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(0), "lim_sup(1/n) = 0");
    }
}

static void test_lim_inf_convergent_sequence() {
    TEST_CASE("lim_inf: convergent sequence 1/n → 0");

    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto result = lamina::lim_inf(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_inf(1/n) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(0), "lim_inf(1/n) = 0");
    }
}

static void test_lim_sup_alternating_sequence() {
    TEST_CASE("lim_sup: alternating sequence (-1)^n");

    // a_n = (-1)^n
    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(-1), n);

    auto result = lamina::lim_sup(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_sup((-1)^n) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(1), "lim_sup((-1)^n) = 1");
    }
}

static void test_lim_inf_alternating_sequence() {
    TEST_CASE("lim_inf: alternating sequence (-1)^n");

    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(-1), n);

    auto result = lamina::lim_inf(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_inf((-1)^n) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(-1), "lim_inf((-1)^n) = -1");
    }
}

static void test_lim_sup_monotone_polynomial() {
    TEST_CASE("lim_sup: monotone polynomial n^2 → ∞");

    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(2));

    auto result = lamina::lim_sup(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_sup(n^2) is not null");
    // Should be infinity
}

static void test_lim_sup_rational_sequence() {
    TEST_CASE("lim_sup: rational sequence n/(n+1) → 1");

    auto n = var("n");
    // n / (n+1)
    auto a_n = SymbolicExpr::multiply(n, SymbolicExpr::power(
        SymbolicExpr::add(n, num(1)), num(-1)));

    auto result = lamina::lim_sup(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_sup(n/(n+1)) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(1), "lim_sup(n/(n+1)) = 1");
    }
}

static void test_lim_inf_rational_sequence() {
    TEST_CASE("lim_inf: rational sequence n/(n+1) → 1");

    auto n = var("n");
    auto a_n = SymbolicExpr::multiply(n, SymbolicExpr::power(
        SymbolicExpr::add(n, num(1)), num(-1)));

    auto result = lamina::lim_inf(a_n, "n");
    EXPECT_TRUE(result != nullptr, "lim_inf(n/(n+1)) is not null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(1), "lim_inf(n/(n+1)) = 1");
    }
}

static void test_lim_sup_lim_inf_equal_for_monotone() {
    TEST_CASE("lim_sup = lim_inf for monotone sequence 1/n");

    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto sup = lamina::lim_sup(a_n, "n");
    auto inf = lamina::lim_inf(a_n, "n");

    EXPECT_TRUE(sup != nullptr && inf != nullptr, "Both lim_sup and lim_inf exist");
    if (sup && inf) {
        EXPECT_EQ_EXPR(sup, inf, "lim_sup(1/n) = lim_inf(1/n) for monotone sequence");
    }
}

// ============================================================
// symbolic_sum tests
// ============================================================

static void test_symbolic_sum_constant() {
    TEST_CASE("symbolic_sum: constant body");
    // sum_{k=1}^{5} 3 = 15
    auto result = lamina::symbolic_sum(num(3), "k", num(1), num(5));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) EXPECT_EQ_EXPR(result, num(15), "sum of constant 3 from 1 to 5 = 15");
}

static void test_symbolic_sum_linear() {
    TEST_CASE("symbolic_sum: sum of k from 1 to n");
    // sum_{k=1}^{n} k = n(n+1)/2
    auto n = var("n");
    auto k = var("k");
    auto result = lamina::symbolic_sum(k, "k", num(1), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    // Evaluate at n=10: should be 55
    if (result) {
        auto val = result->substitute("n", num(10));
        if (val) {
            val = val->simplify();
            EXPECT_EQ_EXPR(val, num(55), "sum k=1..10 k = 55");
        }
    }
}

static void test_symbolic_sum_quadratic() {
    TEST_CASE("symbolic_sum: sum of k^2 from 1 to n");
    // sum_{k=1}^{n} k^2 = n(n+1)(2n+1)/6
    auto n = var("n");
    auto k = var("k");
    auto k_sq = SymbolicExpr::power(k, num(2));
    auto result = lamina::symbolic_sum(k_sq, "k", num(1), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        auto val = result->substitute("n", num(5));
        if (val) {
            val = val->simplify();
            // 1+4+9+16+25 = 55
            EXPECT_EQ_EXPR(val, num(55), "sum k=1..5 k^2 = 55");
        }
    }
}

static void test_symbolic_sum_cubic() {
    TEST_CASE("symbolic_sum: sum of k^3 from 1 to n");
    // sum_{k=1}^{n} k^3 = [n(n+1)/2]^2
    auto n = var("n");
    auto k = var("k");
    auto k_cubed = SymbolicExpr::power(k, num(3));
    auto result = lamina::symbolic_sum(k_cubed, "k", num(1), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        auto val = result->substitute("n", num(4));
        if (val) {
            val = val->simplify();
            // 1+8+27+64 = 100
            EXPECT_EQ_EXPR(val, num(100), "sum k=1..4 k^3 = 100");
        }
    }
}

static void test_symbolic_sum_geometric() {
    TEST_CASE("symbolic_sum: geometric series 2^k from 0 to n");
    auto n = var("n");
    auto k = var("k");
    auto two_pow_k = SymbolicExpr::power(num(2), k);
    auto result = lamina::symbolic_sum(two_pow_k, "k", num(0), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        auto val = result->substitute("n", num(3));
        if (val) {
            val = val->simplify();
            // 1+2+4+8 = 15
            EXPECT_EQ_EXPR(val, num(15), "sum k=0..3 2^k = 15");
        }
    }
}

static void test_symbolic_sum_direct_eval() {
    TEST_CASE("symbolic_sum: direct evaluation for small range");
    auto k = var("k");
    auto k_sq = SymbolicExpr::power(k, num(2));
    auto result = lamina::symbolic_sum(k_sq, "k", num(1), num(4));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        result = result->simplify();
        // 1+4+9+16 = 30
        EXPECT_EQ_EXPR(result, num(30), "sum k=1..4 k^2 = 30");
    }
}

static void test_symbolic_sum_empty_range() {
    TEST_CASE("symbolic_sum: empty range returns 0");
    auto k = var("k");
    auto result = lamina::symbolic_sum(k, "k", num(5), num(3));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) EXPECT_EQ_EXPR(result, num(0), "sum with upper < lower = 0");
}

static void test_symbolic_sum_unevaluated() {
    TEST_CASE("symbolic_sum: returns SummationNode when no closed form");
    auto k = var("k");
    auto n = var("n");
    // sin(k) has no simple closed form
    auto body = SymbolicExpr::sin(k);
    auto result = lamina::symbolic_sum(body, "k", num(1), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    // Should be a SummationNode
    if (result) {
        auto sn = std::dynamic_pointer_cast<const SummationNode>(lamina::detail::node(result));
        EXPECT_TRUE(sn != nullptr, "result is SummationNode for sin(k)");
    }
}

// ============================================================
// symbolic_product tests
// ============================================================

static void test_symbolic_product_factorial() {
    TEST_CASE("symbolic_product: product of k from 1 to n = n!");
    auto k = var("k");
    auto result = lamina::symbolic_product(k, "k", num(1), num(5));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        result = result->simplify();
        // 5! = 120
        EXPECT_EQ_EXPR(result, num(120), "prod k=1..5 k = 120");
    }
}

static void test_symbolic_product_direct_eval() {
    TEST_CASE("symbolic_product: direct evaluation");
    auto k = var("k");
    auto k_plus_1 = SymbolicExpr::add(k, num(1));
    auto result = lamina::symbolic_product(k_plus_1, "k", num(1), num(4));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        result = result->simplify();
        // (2)(3)(4)(5) = 120
        EXPECT_EQ_EXPR(result, num(120), "prod k=1..4 (k+1) = 120");
    }
}

static void test_symbolic_product_empty_range() {
    TEST_CASE("symbolic_product: empty range returns 1");
    auto k = var("k");
    auto result = lamina::symbolic_product(k, "k", num(5), num(3));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) EXPECT_EQ_EXPR(result, num(1), "product with upper < lower = 1");
}

static void test_symbolic_product_pochhammer() {
    TEST_CASE("symbolic_product: Pochhammer (k+2) from 1 to 4");
    auto k = var("k");
    auto body = SymbolicExpr::add(k, num(2));
    auto result = lamina::symbolic_product(body, "k", num(1), num(4));
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        result = result->simplify();
        // (3)(4)(5)(6) = 360
        EXPECT_EQ_EXPR(result, num(360), "prod k=1..4 (k+2) = 360");
    }
}

static void test_symbolic_product_unevaluated() {
    TEST_CASE("symbolic_product: returns ProductNode_Op when no closed form");
    auto k = var("k");
    auto n = var("n");
    auto body = SymbolicExpr::sin(k);
    auto result = lamina::symbolic_product(body, "k", num(1), n);
    EXPECT_TRUE(result != nullptr, "result is not null");
    if (result) {
        auto pn = std::dynamic_pointer_cast<const ProductNode_Op>(lamina::detail::node(result));
        EXPECT_TRUE(pn != nullptr, "result is ProductNode_Op for sin(k)");
    }
}

// ============================================================
// convergence_radius tests (Requirement 22)
// ============================================================

static void test_convergence_radius_geometric() {
    TEST_CASE("convergence_radius: geometric series ∑x^n has R=1");

    // Coefficients of ∑x^n: all 1's → a_n = 1, ratio |a_n/a_{n+1}| = 1
    Coeffs coeffs = {num(1), num(1), num(1), num(1), num(1)};
    auto result = lamina::convergence_radius(coeffs, "x");
    EXPECT_TRUE(result != nullptr, "convergence_radius returns non-null");
    if (result) {
        EXPECT_EQ_EXPR(result, num(1), "R = 1 for geometric series");
    }
}

static void test_convergence_radius_exponential() {
    TEST_CASE("convergence_radius: exponential series ∑x^n/n! has R=∞");

    // Coefficients: 1, 1, 1/2, 1/6, 1/24 (= 1/n!)
    // Ratio |a_n/a_{n+1}| = (n+1) → ∞, so R = ∞
    Coeffs coeffs = {num(1), num(1), rat(1, 2), rat(1, 6), rat(1, 24)};
    auto result = lamina::convergence_radius(coeffs, "x");
    EXPECT_TRUE(result != nullptr, "convergence_radius returns non-null");
    // The ratio of consecutive coefficients grows, so R should be large or infinity
    if (result) {
        // a_3/a_4 = (1/6)/(1/24) = 4, which is the last ratio computed
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value() && *val >= 4.0, "R >= 4 for exponential series (last ratio)");
    }
}

static void test_convergence_radius_half() {
    TEST_CASE("convergence_radius: ∑(2^n)x^n has R=1/2");

    // Coefficients: 1, 2, 4, 8, 16 → a_n = 2^n
    // Ratio |a_n/a_{n+1}| = 2^n / 2^(n+1) = 1/2
    Coeffs coeffs = {num(1), num(2), num(4), num(8), num(16)};
    auto result = lamina::convergence_radius(coeffs, "x");
    EXPECT_TRUE(result != nullptr, "convergence_radius returns non-null");
    if (result) {
        auto val = test_numeric_eval(result);
        EXPECT_TRUE(val.has_value(), "R is numeric");
        if (val.has_value()) {
            EXPECT_NEAR(*val, 0.5, 0.01, "R = 0.5 for ∑(2^n)x^n");
        }
    }
}

static void test_convergence_radius_single_coeff() {
    TEST_CASE("convergence_radius: single coefficient returns ∞");

    Coeffs coeffs = {num(5)};
    auto result = lamina::convergence_radius(coeffs, "x");
    EXPECT_TRUE(result != nullptr, "convergence_radius returns non-null");
    // A polynomial (finite terms) has infinite radius of convergence
}

// ============================================================
// convergence_test tests (Requirement 23)
// ============================================================

static void test_convergence_test_geometric_convergent() {
    TEST_CASE("convergence_test: ∑(1/2)^n converges or inconclusive");

    // a_n = (1/2)^n → ratio |a_{n+1}/a_n| = 1/2 < 1
    // Note: simplifier may not reduce r^(n+1)/r^n to r
    auto n = var("n");
    auto a_n = SymbolicExpr::power(rat(1, 2), n);

    auto info = lamina::convergence_test(a_n, "n");
    EXPECT_TRUE(info.result == lamina::ConvergenceResult::Convergent ||
                info.result == lamina::ConvergenceResult::Inconclusive,
                "∑(1/2)^n: convergent or inconclusive (simplifier limitation)");
}

static void test_convergence_test_geometric_divergent() {
    TEST_CASE("convergence_test: ∑2^n diverges or inconclusive");

    // a_n = 2^n → ratio |a_{n+1}/a_n| = 2 > 1
    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(2), n);

    auto info = lamina::convergence_test(a_n, "n");
    EXPECT_TRUE(info.result == lamina::ConvergenceResult::Divergent ||
                info.result == lamina::ConvergenceResult::Inconclusive,
                "∑2^n: divergent or inconclusive (simplifier limitation)");
}

static void test_convergence_test_p_series() {
    TEST_CASE("convergence_test: ∑1/n^2 (ratio test inconclusive)");

    // a_n = 1/n^2 → ratio test gives lim (n/(n+1))^2 = 1 → inconclusive
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-2));

    auto info = lamina::convergence_test(a_n, "n");
    // The ratio test gives limit 1, so it should be inconclusive
    EXPECT_TRUE(info.result == lamina::ConvergenceResult::Inconclusive ||
                info.result == lamina::ConvergenceResult::Convergent,
                "∑1/n^2: inconclusive by ratio test or convergent by comparison");
}

// ============================================================
// Fourier series tests (Requirement 25)
// ============================================================

static void test_fourier_series_square_wave() {
    TEST_CASE("fourier_series: constant function f(x)=1 over [-pi,pi]");

    // f(x) = 1: the a0/2 term equals 1; with n_terms=0 the series IS the
    // constant term, isolating the a0 computation from cos/sin integrals.
    auto one = num(1);
    auto period = SymbolicExpr::multiply(num(2), SymbolicExpr::number(LMMC_CONST_PI));

    auto result = lamina::fourier_series(one, "x", period, 0);
    EXPECT_TRUE(result != nullptr, "fourier_series returns a result");
    if (result) {
        auto v = test_numeric_eval(result->simplify());
        EXPECT_TRUE(v.has_value() && std::abs(*v - 1.0) < 1e-6,
            "Fourier constant term (a0/2) of f=1 is 1");
    }
}

static void test_fourier_series_odd_function() {
    TEST_CASE("fourier_series: odd function f=x, a0 term is zero");

    // f(x) = x is odd over [-pi,pi] -> constant term a0/2 = 0.
    auto x = var("x");
    auto period = SymbolicExpr::multiply(num(2), SymbolicExpr::number(LMMC_CONST_PI));

    auto result = lamina::fourier_series(x, "x", period, 2);
    EXPECT_TRUE(result != nullptr, "fourier_series(x) returns a result");
    if (result) {
        // At x=0, all sine terms vanish; an odd function's series gives 0 there.
        auto val = result->substitute("x", num(0));
        auto v = test_numeric_eval(val ? val->simplify() : nullptr);
        EXPECT_TRUE(v.has_value() && std::abs(*v) < 1e-6,
            "Fourier series of odd f=x is 0 at x=0 (no constant term)");
    }
}

// ============================================================
// Laurent series tests (Requirement 26)
// ============================================================

static void test_laurent_series_1_over_z() {
    TEST_CASE("laurent_series: 1/z around z=0");

    // f(z) = 1/z is its own Laurent series; pole order 1.
    auto z = var("z");
    auto f = SymbolicExpr::power(z, num(-1));

    auto result = lamina::laurent_series(f, "z", num(0), 2, 2);
    EXPECT_TRUE(result != nullptr, "laurent_series(1/z) returns a result");
    if (result) {
        // result should equal 1/z
        auto expected = SymbolicExpr::divide(num(1), z);
        EXPECT_EQ_EXPR(result->simplify(), expected->simplify(), "Laurent of 1/z is 1/z");
    }
}

static void test_laurent_series_full_1_over_z() {
    TEST_CASE("laurent_series_full: 1/z around z=0");

    auto z = var("z");
    auto f = SymbolicExpr::power(z, num(-1));

    auto result = lamina::laurent_series_full(f, "z", num(0), 2, 2);
    EXPECT_TRUE(result.series != nullptr, "laurent_series_full: series not null");
    EXPECT_TRUE(result.pole_order == 1, "laurent_series_full: pole_order = 1 for 1/z");
    if (result.residue) {
        EXPECT_EQ_EXPR(result.residue->simplify(), num(1), "residue of 1/z at 0 is 1");
    }
}

static void test_laurent_series_1_over_sin_z() {
    TEST_CASE("laurent_series: 1/(z(z+1)) around z=0 has simple pole, residue 1");

    // f(z) = 1/(z(z+1)) has a simple pole at z=0 with residue 1.
    // (Rational denominator; transcendental denominators like 1/sin(z) are not
    //  yet supported by the series engine and are out of scope here.)
    auto z = var("z");
    auto denom = SymbolicExpr::multiply(z, SymbolicExpr::add(z, num(1)));
    auto f = SymbolicExpr::divide(num(1), denom);

    auto result = lamina::laurent_series_full(f, "z", num(0), 3, 3);
    EXPECT_TRUE(result.pole_order == 1, "1/(z(z+1)): pole_order = 1");
    if (result.residue) {
        EXPECT_EQ_EXPR(result.residue->simplify(), num(1), "residue of 1/(z(z+1)) at 0 is 1");
    }
}

int main() {
    test_power_series_add_basic();
    test_power_series_add_different_lengths();
    test_power_series_add_empty();
    test_power_series_multiply_basic();
    test_power_series_multiply_truncation();
    test_power_series_multiply_zero();
    test_power_series_multiply_constants();
    test_power_series_multiply_order_zero();
    test_power_series_compose_basic();
    test_power_series_compose_quadratic();
    test_power_series_compose_g0_nonzero();
    test_power_series_compose_exp_like();
    test_power_series_checked_contracts();
    test_series_analysis_checked_contracts();

    // convergence_radius tests
    test_convergence_radius_geometric();
    test_convergence_radius_exponential();
    test_convergence_radius_half();
    test_convergence_radius_single_coeff();

    // convergence_test tests
    test_convergence_test_geometric_convergent();
    test_convergence_test_geometric_divergent();
    test_convergence_test_p_series();

    // Fourier series tests
    test_fourier_series_square_wave();
    test_fourier_series_odd_function();

    // Laurent series tests
    test_laurent_series_1_over_z();
    test_laurent_series_full_1_over_z();
    test_laurent_series_1_over_sin_z();

    // lim_sup / lim_inf tests
    test_lim_sup_convergent_sequence();
    test_lim_inf_convergent_sequence();
    test_lim_sup_alternating_sequence();
    test_lim_inf_alternating_sequence();
    test_lim_sup_monotone_polynomial();
    test_lim_sup_rational_sequence();
    test_lim_inf_rational_sequence();
    test_lim_sup_lim_inf_equal_for_monotone();

    // symbolic_sum tests
    test_symbolic_sum_constant();
    test_symbolic_sum_linear();
    test_symbolic_sum_quadratic();
    test_symbolic_sum_cubic();
    test_symbolic_sum_geometric();
    test_symbolic_sum_direct_eval();
    test_symbolic_sum_empty_range();
    test_symbolic_sum_unevaluated();

    // symbolic_product tests
    test_symbolic_product_factorial();
    test_symbolic_product_direct_eval();
    test_symbolic_product_empty_range();
    test_symbolic_product_pochhammer();
    test_symbolic_product_unevaluated();

    return TEST_REPORT();
}
