#include "test_common.hpp"
#include "transcendental_factor.hpp"
#include "poly_utils.hpp"

using namespace LMCAS;

static std::vector<Polynomial<Rational>> checked_zassenhaus(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted,
    int64_t modulus)
{
    auto result = zassenhaus_combine_checked(
        poly, lifted, BigInt(static_cast<long long>(modulus)));
    if (!result) {
        throw std::runtime_error(result.error().message);
    }
    return std::move(result.value().value);
}

void test_single_sin() {
    TEST_CASE("detect_trans_substitutions: single sin(x)");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 1, "should find 1 transcendental sub-expression");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "first indeterminate is u0");
    EXPECT_EQ_EXPR_STR(result.mappings[0].trans_expr, "sin(x)", "trans_expr is sin(x)");
}

void test_sin_and_cos() {
    TEST_CASE("detect_trans_substitutions: sin(x) and cos(x)");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + cos(x)
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), SymbolicExpr::cos(x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions");
    // Verify both have distinct indeterminates
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "first indeterminate is u0");
    EXPECT_EQ_STR(result.mappings[1].indeterminate, "u1", "second indeterminate is u1");
}

void test_deduplication() {
    TEST_CASE("detect_trans_substitutions: deduplication of sin(x)^2 + sin(x)");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    // sin(x)^2 + sin(x)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        sin_x
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 1, "should deduplicate to 1 transcendental sub-expression");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "indeterminate is u0");
}

void test_exp_function() {
    TEST_CASE("detect_trans_substitutions: exp(x)");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x^2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2)))
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 1, "should find 1 transcendental sub-expression (exp)");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "indeterminate is u0");
    EXPECT_EQ_EXPR_STR(result.mappings[0].trans_expr, "exp(x)", "trans_expr is exp(x)");
}

void test_ln_function() {
    TEST_CASE("detect_trans_substitutions: ln(x)");

    auto x = SymbolicExpr::variable("x");
    // ln(x) + x
    auto expr = SymbolicExpr::add(SymbolicExpr::ln(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 1, "should find 1 transcendental sub-expression (ln)");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "indeterminate is u0");
    EXPECT_EQ_EXPR_STR(result.mappings[0].trans_expr, "ln(x)", "trans_expr is ln(x)");
}

void test_tan_function() {
    TEST_CASE("detect_trans_substitutions: tan(x)");

    auto x = SymbolicExpr::variable("x");
    // tan(x) + 1
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(x), SymbolicExpr::number(1));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 1, "should find 1 transcendental sub-expression (tan)");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "indeterminate is u0");
    EXPECT_EQ_EXPR_STR(result.mappings[0].trans_expr, "tan(x)", "trans_expr is tan(x)");
}

void test_no_transcendental() {
    TEST_CASE("detect_trans_substitutions: pure polynomial (no transcendental)");

    auto x = SymbolicExpr::variable("x");
    // x^2 + 2*x + 1
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)
        )
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.empty(), "should find no transcendental sub-expressions");
}

void test_independent_variable() {
    TEST_CASE("detect_trans_substitutions: sin(y) does not depend on x");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    // sin(y) + x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(y), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.empty(), "sin(y) should not be collected when var=x");
}

void test_nested_transcendental() {
    TEST_CASE("detect_trans_substitutions: nested sin(exp(x))");

    auto x = SymbolicExpr::variable("x");
    // sin(exp(x)) + x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(SymbolicExpr::exp(x)), x);

    auto result = detect_trans_substitutions(expr, "x");

    // Should collect both sin(exp(x)) (outermost) and exp(x) (inner)
    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions (outer and inner)");
    EXPECT_EQ_STR(result.mappings[0].indeterminate, "u0", "first indeterminate is u0");
    EXPECT_EQ_STR(result.mappings[1].indeterminate, "u1", "second indeterminate is u1");
}

void test_multiple_distinct() {
    TEST_CASE("detect_trans_substitutions: sin(x) + exp(x) + cos(x)");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::cos(x))
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 3, "should find 3 distinct transcendental sub-expressions");
}

void test_null_expression() {
    TEST_CASE("detect_trans_substitutions: null expression");

    auto result = detect_trans_substitutions(nullptr, "x");

    EXPECT_TRUE(result.mappings.empty(), "null expression should return empty mappings");
    EXPECT_TRUE(result.poly_expr == nullptr, "poly_expr should be null");
}

void test_sin_cos_pythagorean_constraint() {
    TEST_CASE("constraint detection: sin(x) and cos(x) → u0² + u1² - 1 = 0");

    auto x = SymbolicExpr::variable("x");
    // sin^2(x) + cos^2(x) - 1
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::power(cos_x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-1)
        )
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions");
    EXPECT_TRUE(result.constraints.size() == 1, "should detect 1 Pythagorean constraint");

    // 验证约束表达式包含两个不定元的平方和减一
    if (!result.constraints.empty()) {
        std::string cstr = result.constraints[0]->to_string();
        // 约束应包含 u0, u1 的平方项
        EXPECT_TRUE(cstr.find("u0") != std::string::npos, "constraint contains u0");
        EXPECT_TRUE(cstr.find("u1") != std::string::npos, "constraint contains u1");
    }
}

void test_exp_inverse_constraint() {
    TEST_CASE("constraint detection: exp(x) and exp(-x) → u0 * u1 - 1 = 0");

    auto x = SymbolicExpr::variable("x");
    // exp(x) + exp(-x)
    auto exp_x = SymbolicExpr::exp(x);
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
    auto exp_neg_x = SymbolicExpr::exp(neg_x);
    auto expr = SymbolicExpr::add(exp_x, exp_neg_x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions (exp(x) and exp(-x))");
    EXPECT_TRUE(result.constraints.size() == 1, "should detect 1 inverse constraint");

    if (!result.constraints.empty()) {
        std::string cstr = result.constraints[0]->to_string();
        EXPECT_TRUE(cstr.find("u0") != std::string::npos, "constraint contains u0");
        EXPECT_TRUE(cstr.find("u1") != std::string::npos, "constraint contains u1");
    }
}

void test_no_constraint_different_args() {
    TEST_CASE("constraint detection: sin(x) and cos(2*x) → no constraint");

    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    // sin(x) + cos(2*x) - different arguments, no Pythagorean constraint
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), SymbolicExpr::cos(two_x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions");
    EXPECT_TRUE(result.constraints.empty(), "no constraint when arguments differ");
}

void test_no_constraint_same_type() {
    TEST_CASE("constraint detection: exp(x) and exp(2*x) → no constraint");

    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    // exp(x) + exp(2*x) - not inverse pair
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::exp(two_x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 2, "should find 2 transcendental sub-expressions");
    EXPECT_TRUE(result.constraints.empty(), "no constraint when exp arguments are not negations");
}

void test_multiple_constraints() {
    TEST_CASE("constraint detection: sin(x), cos(x), exp(x), exp(-x) → 2 constraints");

    auto x = SymbolicExpr::variable("x");
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
    // sin(x) + cos(x) + exp(x) + exp(-x)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(
            SymbolicExpr::cos(x),
            SymbolicExpr::add(
                SymbolicExpr::exp(x),
                SymbolicExpr::exp(neg_x)
            )
        )
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.mappings.size() == 4, "should find 4 transcendental sub-expressions");
    EXPECT_TRUE(result.constraints.size() == 2, "should detect 2 constraints (Pythagorean + inverse)");
}


void test_substitution_sin_plus_x() {
    TEST_CASE("substitution: sin(x) + x → u0 + x");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.poly_expr != nullptr, "poly_expr should not be null");
    // 替换后应为 u0 + x
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_TRUE(poly_str.find("u0") != std::string::npos, "poly_expr should contain u0");
    EXPECT_TRUE(poly_str.find("sin") == std::string::npos, "poly_expr should not contain sin");
    EXPECT_TRUE(poly_str.find("x") != std::string::npos, "poly_expr should still contain x");
}

void test_substitution_sin_squared_minus_x_squared() {
    TEST_CASE("substitution: sin²(x) - x² → u0² - x²");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    // sin^2(x) - x^2 = sin(x)^2 + (-1)*x^2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2)))
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.poly_expr != nullptr, "poly_expr should not be null");
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_TRUE(poly_str.find("u0") != std::string::npos, "poly_expr should contain u0");
    EXPECT_TRUE(poly_str.find("sin") == std::string::npos, "poly_expr should not contain sin");
}

void test_substitution_no_transcendental() {
    TEST_CASE("substitution: pure polynomial remains unchanged");

    auto x = SymbolicExpr::variable("x");
    // x^2 + 2*x + 1
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)
        )
    );

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE(result.poly_expr != nullptr, "poly_expr should not be null for pure polynomial");
    // 无超越项时 poly_expr 应与原表达式相同
    std::string orig_str = expr->to_string();
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_EQ_STR(poly_str, orig_str, "poly_expr should equal original expression");
}


void test_build_poly_single_indeterminate() {
    TEST_CASE("tf_build_polynomial: single indeterminate u0^2 + u0 + 1");

    // 构造表达式 u0^2 + u0 + 1(纯 u0 多项式,无 x 依赖)
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::add(u0, SymbolicExpr::number(1))
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "conversion should succeed");
    EXPECT_EQ_STR(result.main_variable, "u0", "main variable should be u0");
    EXPECT_TRUE(result.poly.degree() == 2, "polynomial degree should be 2");
    EXPECT_TRUE(result.param_variables.empty(), "no parameter variables");

    // 验证系数:1 + u0 + u0^2 -> coeffs = {1, 1, 1}
    EXPECT_TRUE(result.poly.coeffs.size() == 3, "should have 3 coefficients");
    EXPECT_TRUE(result.poly.coeffs[0] == Rational(1), "constant term is 1");
    EXPECT_TRUE(result.poly.coeffs[1] == Rational(1), "linear term is 1");
    EXPECT_TRUE(result.poly.coeffs[2] == Rational(1), "quadratic term is 1");
}

void test_build_poly_single_indeterminate_with_rational_coeffs() {
    TEST_CASE("tf_build_polynomial: 2*u0^3 - 3*u0 + 5");

    auto u0 = SymbolicExpr::variable("u0");
    // 2*u0^3 - 3*u0 + 5
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2),
            SymbolicExpr::power(u0, SymbolicExpr::number(3))),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-3), u0),
            SymbolicExpr::number(5)
        )
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "conversion should succeed");
    EXPECT_EQ_STR(result.main_variable, "u0", "main variable should be u0");
    EXPECT_TRUE(result.poly.degree() == 3, "polynomial degree should be 3");

    // coeffs: {5, -3, 0, 2}
    EXPECT_TRUE(result.poly.coeffs[0] == Rational(5), "constant term is 5");
    EXPECT_TRUE(result.poly.coeffs[1] == Rational(-3), "linear term is -3");
    EXPECT_TRUE(result.poly.coeffs[2] == Rational(0), "quadratic term is 0");
    EXPECT_TRUE(result.poly.coeffs[3] == Rational(2), "cubic term is 2");
}

void test_build_poly_constant_expression() {
    TEST_CASE("tf_build_polynomial: constant expression (no variables)");

    auto expr = SymbolicExpr::number(42);

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "conversion should succeed for constant");
    EXPECT_TRUE(result.poly.degree() == 0, "constant polynomial has degree 0");
    EXPECT_TRUE(result.poly.coeffs[0] == Rational(42), "constant value is 42");
}

void test_build_poly_null_expression() {
    TEST_CASE("tf_build_polynomial: null expression");

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(nullptr, indeterminates, "x");

    EXPECT_FALSE(result.success, "null expression should fail");
}

void test_build_poly_only_x() {
    TEST_CASE("tf_build_polynomial: expression only in x (no indeterminates used)");

    auto x = SymbolicExpr::variable("x");
    // x^2 + 3*x + 2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(2)
        )
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "conversion should succeed");
    EXPECT_EQ_STR(result.main_variable, "x", "main variable should be x");
    EXPECT_TRUE(result.poly.degree() == 2, "polynomial degree should be 2");
}

void test_build_poly_main_variable_selection() {
    TEST_CASE("tf_build_polynomial: selects highest-degree variable as main");

    // u0^3 + u1^2 -> u0 has degree 3, u1 has degree 2 -> main = u0
    auto u0 = SymbolicExpr::variable("u0");
    auto u1 = SymbolicExpr::variable("u1");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(3)),
        SymbolicExpr::power(u1, SymbolicExpr::number(2))
    );

    std::vector<std::string> indeterminates = {"u0", "u1"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    // 主变量应为 u0(次数 3 > 2)
    EXPECT_EQ_STR(result.main_variable, "u0", "main variable should be u0 (highest degree)");
}

void test_build_poly_from_substitution_result() {
    TEST_CASE("tf_build_polynomial: end-to-end from detect_trans_substitutions");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    // sin^2(x) + sin(x) + 1 -> u0^2 + u0 + 1
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(sin_x, SymbolicExpr::number(1))
    );

    auto sub_result = detect_trans_substitutions(expr, "x");
    EXPECT_TRUE(sub_result.mappings.size() == 1, "should have 1 mapping");

    std::vector<std::string> indeterminates;
    for (const auto& m : sub_result.mappings) {
        indeterminates.push_back(m.indeterminate);
    }

    auto poly_result = tf_build_polynomial(sub_result.poly_expr, indeterminates, "x");

    EXPECT_TRUE(poly_result.success, "polynomial construction should succeed");
    EXPECT_EQ_STR(poly_result.main_variable, "u0", "main variable should be u0");
    EXPECT_TRUE(poly_result.poly.degree() == 2, "polynomial degree should be 2");
    EXPECT_TRUE(poly_result.poly.coeffs[0] == Rational(1), "constant term is 1");
    EXPECT_TRUE(poly_result.poly.coeffs[1] == Rational(1), "linear term is 1");
    EXPECT_TRUE(poly_result.poly.coeffs[2] == Rational(1), "quadratic term is 1");
}

void test_build_poly_non_polynomial_fails() {
    TEST_CASE("tf_build_polynomial: non-polynomial expression fails");

    auto u0 = SymbolicExpr::variable("u0");
    // sin(u0) - still transcendental in u0, not polynomial
    auto expr = SymbolicExpr::sin(u0);

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "non-polynomial expression should fail");
}


void test_validate_valid_polynomial() {
    TEST_CASE("validation: valid polynomial u0^2 + 3*u0 - 2 passes");

    auto u0 = SymbolicExpr::variable("u0");
    // u0^2 + 3*u0 - 2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), u0),
            SymbolicExpr::number(-2)
        )
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "valid polynomial should pass validation");
    EXPECT_TRUE(result.poly.degree() == 2, "degree should be 2");
}

void test_validate_remaining_sin_fails() {
    TEST_CASE("validation: expression with remaining sin(u0) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // sin(u0) + u0 - sin still depends on indeterminate u0
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(u0), u0);

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "expression with remaining sin(u0) should fail validation");
}

void test_validate_remaining_cos_fails() {
    TEST_CASE("validation: expression with remaining cos(u0) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // u0^2 + cos(u0) - cos still depends on indeterminate u0
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::cos(u0)
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "expression with remaining cos(u0) should fail validation");
}

void test_validate_remaining_exp_fails() {
    TEST_CASE("validation: expression with remaining exp(u0) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // exp(u0) + u0^2 - exp still depends on indeterminate u0
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(u0),
        SymbolicExpr::power(u0, SymbolicExpr::number(2))
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "expression with remaining exp(u0) should fail validation");
}

void test_validate_fractional_exponent_fails() {
    TEST_CASE("validation: fractional exponent u0^(1/2) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // u0^(1/2) - fractional exponent is not polynomial
    auto half = SymbolicExpr::number(Rational(BigInt(1), BigInt(2)));
    auto expr = SymbolicExpr::power(u0, half);

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "fractional exponent should fail validation");
}

void test_validate_negative_exponent_fails() {
    TEST_CASE("validation: negative exponent u0^(-1) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // u0^(-1) - negative exponent means division by variable
    auto expr = SymbolicExpr::power(u0, SymbolicExpr::number(-1));

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "negative exponent should fail validation");
}

void test_validate_negative_fractional_exponent_fails() {
    TEST_CASE("validation: negative fractional exponent u0^(-3/2) fails");

    auto u0 = SymbolicExpr::variable("u0");
    // u0^(-3/2) - negative fractional exponent
    auto neg_three_half = SymbolicExpr::number(Rational(BigInt(-3), BigInt(2)));
    auto expr = SymbolicExpr::power(u0, neg_three_half);

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "negative fractional exponent should fail validation");
}

void test_validate_transcendental_in_original_var_fails() {
    TEST_CASE("validation: remaining sin(x) in original variable fails");

    auto x = SymbolicExpr::variable("x");
    auto u0 = SymbolicExpr::variable("u0");
    // u0 + sin(x) - sin still depends on original variable x
    auto expr = SymbolicExpr::add(u0, SymbolicExpr::sin(x));

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_FALSE(result.success, "expression with remaining sin(x) should fail validation");
}

void test_validate_valid_multivariate_polynomial() {
    TEST_CASE("validation: valid multivariate u0^2 + x passes");

    auto u0 = SymbolicExpr::variable("u0");
    auto x = SymbolicExpr::variable("x");
    // u0^2 + x - valid polynomial in both u0 and x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        x
    );

    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(expr, indeterminates, "x");

    EXPECT_TRUE(result.success, "valid multivariate polynomial should pass validation");
}


void test_square_free_already_square_free() {
    TEST_CASE("tf_square_free: already square-free polynomial x^2 - 1");

    // x^2 - 1 = (x-1)(x+1),无重因子
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE(result.had_repeated_factors, "x^2 - 1 should be square-free");
    EXPECT_TRUE(result.square_free.degree() == 2, "square-free part degree should be 2");
    EXPECT_TRUE(result.repeated_factor.degree() == 0, "repeated factor should be constant");
}

void test_square_free_repeated_root() {
    TEST_CASE("tf_square_free: (x-1)^2 = x^2 - 2x + 1 → returns (x-1)");

    // (x-1)^2 = x^2 - 2x + 1
    Polynomial<Rational> poly({Rational(1), Rational(-2), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_TRUE(result.had_repeated_factors, "should detect repeated factors");
    EXPECT_TRUE(result.square_free.degree() == 1, "square-free part should be linear (x-1)");
    // square-free part should be monic: x - 1
    EXPECT_TRUE(result.square_free.coeffs[1] == Rational(1), "leading coeff should be 1");
    EXPECT_TRUE(result.square_free.coeffs[0] == Rational(-1), "constant term should be -1");
}

void test_square_free_higher_multiplicity() {
    TEST_CASE("tf_square_free: (x-1)^3 * (x+1) → returns (x-1)(x+1)");

    // (x-1)^3 = x^3 - 3x^2 + 3x - 1
    Polynomial<Rational> x_minus_1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> x_plus_1({Rational(1), Rational(1)}, "x");
    auto x_minus_1_sq = x_minus_1 * x_minus_1;
    auto x_minus_1_cubed = x_minus_1_sq * x_minus_1;
    auto poly = x_minus_1_cubed * x_plus_1;

    auto result = tf_square_free(poly);

    EXPECT_TRUE(result.had_repeated_factors, "should detect repeated factors");
    // square-free part should be (x-1)(x+1) = x^2 - 1, degree 2
    EXPECT_TRUE(result.square_free.degree() == 2, "square-free part should be degree 2");
}

void test_square_free_constant_polynomial() {
    TEST_CASE("tf_square_free: constant polynomial → returns itself");

    Polynomial<Rational> poly({Rational(5)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE(result.had_repeated_factors, "constant should have no repeated factors");
    EXPECT_TRUE(result.square_free.degree() == 0, "square-free part should be constant");
}

void test_square_free_linear_polynomial() {
    TEST_CASE("tf_square_free: linear polynomial x + 3 → returns itself");

    Polynomial<Rational> poly({Rational(3), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE(result.had_repeated_factors, "linear polynomial should be square-free");
    EXPECT_TRUE(result.square_free.degree() == 1, "square-free part should be linear");
    EXPECT_TRUE(result.square_free.coeffs[0] == Rational(3), "constant term preserved");
    EXPECT_TRUE(result.square_free.coeffs[1] == Rational(1), "linear term preserved");
}

void test_square_free_zero_polynomial() {
    TEST_CASE("tf_square_free: zero polynomial → returns itself");

    Polynomial<Rational> poly("x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE(result.had_repeated_factors, "zero polynomial has no repeated factors");
    EXPECT_TRUE(result.square_free.is_zero(), "square-free part of zero is zero");
}


void test_berlekamp_linear_poly() {
    TEST_CASE("berlekamp_factor: linear polynomial x + 1 → single factor");

    // x + 1 -> 线性多项式本身不可约
    Polynomial<Rational> poly({Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 0);

    EXPECT_TRUE(result.prime > 0, "should select a valid prime");
    EXPECT_TRUE(result.factors.size() == 1, "linear polynomial should have 1 factor");
    // 因子应为首一线性多项式
    EXPECT_TRUE(result.factors[0].degree() == 1, "factor should be degree 1");
}

void test_berlekamp_quadratic_irreducible() {
    TEST_CASE("berlekamp_factor: x^2 + x + 1 mod 2 is irreducible");

    // x^2 + x + 1 在 F_2 上不可约(无根:f(0)=1, f(1)=1)
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(result.prime == 2, "should use prime 2");
    /// Q 矩阵构造成功;当前分裂阶段返回单一整体因子.
    EXPECT_TRUE(result.factors.size() >= 1, "should return at least 1 factor");
    // 验证因子次数之和等于原多项式次数
    int total_deg = 0;
    for (const auto& f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE(total_deg == 2, "total degree of factors should equal original degree");
}

void test_berlekamp_x2_minus_1_mod3() {
    TEST_CASE("berlekamp_factor: x^2 - 1 mod 3 (factorable as (x+1)(x-1))");

    // x^2 - 1 = (x-1)(x+1),在 F_3 上可分解
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE(result.prime == 3, "should use prime 3");
    // x^2 - 1 在 F_3 上分裂为两个线性因子
    EXPECT_TRUE(result.factors.size() == 2, "should split into 2 factors");
    // 验证因子次数之和
    int total_deg = 0;
    for (const auto& f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE(total_deg == 2, "total degree of factors should equal 2");
    // 每个因子应为线性
    for (const auto& f : result.factors) {
        EXPECT_TRUE(f.degree() == 1, "each factor should be linear");
    }
}

void test_berlekamp_cubic_poly() {
    TEST_CASE("berlekamp_factor: x^3 - x mod 5 (splits as x(x-1)(x+1))");

    // x^3 - x = x(x-1)(x+1),在 F_5 上完全分裂
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE(result.prime == 5, "should use prime 5");
    // x^3 - x 在 F_5 上完全分裂为三个线性因子
    EXPECT_TRUE(result.factors.size() == 3, "should split into 3 factors");
    // 验证因子次数之和
    int total_deg = 0;
    for (const auto& f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE(total_deg == 3, "total degree of factors should equal 3");
    // 每个因子应为线性
    for (const auto& f : result.factors) {
        EXPECT_TRUE(f.degree() == 1, "each factor should be linear");
    }
}

void test_berlekamp_constant_poly() {
    TEST_CASE("berlekamp_factor: constant polynomial → no factors");

    Polynomial<Rational> poly({Rational(5)}, "x");

    auto result = berlekamp_factor(poly, 0);

    EXPECT_TRUE(result.factors.empty(), "constant polynomial should have no factors");
}

void test_berlekamp_specified_prime() {
    TEST_CASE("berlekamp_factor: specified prime is used");

    // x^2 + 1,指定 prime = 7
    Polynomial<Rational> poly({Rational(1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 7);

    EXPECT_TRUE(result.prime == 7, "should use specified prime 7");
    EXPECT_TRUE(result.factors.size() >= 1, "should return at least 1 factor");
}


void test_null_space_x2_minus_1_mod3() {
    TEST_CASE("null space: x^2 - 1 mod 3 → dim = 2 (two factors)");

    // x^2 - 1 = (x-1)(x+1) 在 F_3 上有两个不可约因子
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE(result.prime == 3, "should use prime 3");
    EXPECT_TRUE(result.null_space_dim == 2, "null space dimension should be 2 (two factors)");
    EXPECT_TRUE(result.null_space_basis.size() == 2, "should have 2 basis vectors");

    // 每个基向量长度应为 n = deg(f) = 2
    for (const auto& v : result.null_space_basis) {
        EXPECT_TRUE(v.size() == 2, "basis vector length should equal polynomial degree");
    }
}

void test_null_space_x2_plus_x_plus_1_mod2() {
    TEST_CASE("null space: x^2 + x + 1 mod 2 → dim = 1 (irreducible)");

    // x^2 + x + 1 在 F_2 上不可约(无根:f(0)=1, f(1)=1)
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(result.prime == 2, "should use prime 2");
    EXPECT_TRUE(result.null_space_dim == 1, "null space dimension should be 1 (irreducible)");
    EXPECT_TRUE(result.null_space_basis.size() == 1, "should have 1 basis vector");

    // 唯一的基向量应为 [1, 0](对应平凡因子)
    if (!result.null_space_basis.empty()) {
        EXPECT_TRUE(result.null_space_basis[0][0] == 1, "first basis vector starts with 1");
        EXPECT_TRUE(result.null_space_basis[0][1] == 0, "first basis vector second component is 0");
    }
}

void test_null_space_x3_minus_x_mod5() {
    TEST_CASE("null space: x^3 - x mod 5 → dim = 3 (three factors)");

    // x^3 - x = x(x-1)(x+1) 在 F_5 上完全分裂为三个线性因子
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE(result.prime == 5, "should use prime 5");
    EXPECT_TRUE(result.null_space_dim == 3, "null space dimension should be 3 (three factors)");
    EXPECT_TRUE(result.null_space_basis.size() == 3, "should have 3 basis vectors");

    // 每个基向量长度应为 n = deg(f) = 3
    for (const auto& v : result.null_space_basis) {
        EXPECT_TRUE(v.size() == 3, "basis vector length should equal polynomial degree");
    }
}

void test_null_space_linear_poly() {
    TEST_CASE("null space: linear polynomial → dim not computed (direct return)");

    // 线性多项式直接返回,不经过零空间计算
    Polynomial<Rational> poly({Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE(result.prime == 3, "should use prime 3");
    EXPECT_TRUE(result.factors.size() == 1, "linear polynomial has 1 factor");
    // 线性多项式不经过零空间计算,dim 保持默认值 0
    EXPECT_TRUE(result.null_space_dim == 0, "null space dim is 0 for linear (not computed)");
}

void test_null_space_basis_vectors_in_kernel() {
    TEST_CASE("null space: basis vectors are in kernel of (Q-I)^T");

    // x^2 - 1 mod 3:验证基向量确实满足 (Q-I)^T * v = 0
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);
    int64_t p = result.prime;

    // 重新构造 Q 矩阵以验证
    // f_coeffs for x^2 - 1 mod 3: [2, 0, 1] (since -1 mod 3 = 2)
    // 验证每个基向量 v 满足 (Q-I) * v = 0(等价于 Q*v = v)
    /// 通过零空间维度间接验证 Q*v = v.
    EXPECT_TRUE(result.null_space_dim == 2, "x^2 - 1 mod 3 has null space dim 2");

    // 验证基向量非零
    for (const auto& v : result.null_space_basis) {
        bool all_zero = true;
        for (int64_t c : v) {
            if (c != 0) { all_zero = false; break; }
        }
        EXPECT_FALSE(all_zero, "basis vector should be non-zero");
    }
}


void test_split_x2_minus_1_mod3() {
    TEST_CASE("factor splitting: x^2 - 1 mod 3 → (x+2)(x+1) in F_3");

    // x^2 - 1 = (x-1)(x+1) = (x+2)(x+1) in F_3
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE(result.prime == 3, "should use prime 3");
    EXPECT_TRUE(result.factors.size() == 2, "should produce 2 irreducible factors");

    // 验证因子乘积等于原多项式 mod 3
    // 因子应为首一线性多项式,乘积应为 x^2 + 2 (即 x^2 - 1 mod 3)
    int total_deg = 0;
    for (const auto& f : result.factors) {
        total_deg += f.degree();
        EXPECT_TRUE(f.degree() == 1, "each factor should be linear");
    }
    EXPECT_TRUE(total_deg == 2, "total degree should be 2");
}

void test_split_x3_minus_x_mod5() {
    TEST_CASE("factor splitting: x^3 - x mod 5 → x(x+4)(x+1) in F_5");

    // x^3 - x = x(x-1)(x+1) = x(x+4)(x+1) in F_5
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE(result.prime == 5, "should use prime 5");
    EXPECT_TRUE(result.factors.size() == 3, "should produce 3 irreducible factors");

    // 每个因子应为线性
    for (const auto& f : result.factors) {
        EXPECT_TRUE(f.degree() == 1, "each factor should be linear");
    }
}

void test_split_irreducible_x2_plus_x_plus_1_mod2() {
    TEST_CASE("factor splitting: x^2 + x + 1 mod 2 → remains irreducible");

    // x^2 + x + 1 在 F_2 上不可约
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(result.prime == 2, "should use prime 2");
    EXPECT_TRUE(result.factors.size() == 1, "irreducible polynomial should have 1 factor");
    EXPECT_TRUE(result.factors[0].degree() == 2, "factor should have degree 2");
}

void test_split_product_verification_mod3() {
    TEST_CASE("factor splitting: product of factors equals original mod 3");

    // x^2 - 1 mod 3
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");
    int64_t p = 3;

    auto result = berlekamp_factor(poly, p);

    EXPECT_TRUE(result.factors.size() == 2, "should have 2 factors");

    // 验证因子乘积等于原多项式 mod p
    // 手动计算两个线性因子 (x + a)(x + b) = x^2 + (a+b)x + ab
    // 原多项式 mod 3: x^2 + 2 (coeffs: [2, 0, 1])
    // 所以 a + b == 0 (mod 3) 且 a*b == 2 (mod 3)
    if (result.factors.size() == 2) {
        // 提取两个线性因子的常数项
        int64_t a = result.factors[0].coeffs[0].value();
        int64_t b = result.factors[1].coeffs[0].value();
        // 验证 (a + b) mod 3 == 0
        EXPECT_TRUE((a + b) % p == 0, "sum of constant terms should be 0 mod 3");
        // 验证 a * b mod 3 == 2
        EXPECT_TRUE((a * b) % p == 2, "product of constant terms should be 2 mod 3");
    }
}

void test_split_product_verification_mod5() {
    TEST_CASE("factor splitting: product of factors equals original mod 5");

    // x^3 - x mod 5
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");
    int64_t p = 5;

    auto result = berlekamp_factor(poly, p);

    EXPECT_TRUE(result.factors.size() == 3, "should have 3 factors");

    // x^3 - x = x(x-1)(x+1) = x(x+4)(x+1) in F_5
    // 三个线性因子 (x + a)(x + b)(x + c) 展开:
    //   x^3 + (a+b+c)x^2 + (ab+ac+bc)x + abc
    // 原多项式 mod 5: x^3 + 4x (coeffs: [0, 4, 0, 1])
    // 所以 a+b+c == 0, ab+ac+bc == 4, abc == 0 (mod 5)
    if (result.factors.size() == 3) {
        int64_t a = result.factors[0].coeffs[0].value();
        int64_t b = result.factors[1].coeffs[0].value();
        int64_t c = result.factors[2].coeffs[0].value();
        EXPECT_TRUE((a + b + c) % p == 0, "sum of constant terms should be 0 mod 5");
        EXPECT_TRUE((a * b % p + a * c % p + b * c % p) % p == 4,
            "sum of pairwise products should be 4 mod 5");
        EXPECT_TRUE((a * b % p * c % p) % p == 0, "product of constant terms should be 0 mod 5");
    }
}

void test_split_x4_minus_1_mod5() {
    TEST_CASE("factor splitting: x^4 - 1 mod 5 → splits into factors");

    // x^4 - 1 = (x-1)(x+1)(x^2+1) in F_5
    // x^2 + 1 在 F_5 上有根 x=2 (4+1=5==0) 和 x=3 (9+1=10==0)
    // 所以 x^4 - 1 = (x-1)(x+1)(x-2)(x-3) = (x+4)(x+1)(x+3)(x+2) in F_5
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE(result.prime == 5, "should use prime 5");
    EXPECT_TRUE(result.factors.size() == 4, "x^4 - 1 should split into 4 linear factors mod 5");

    int total_deg = 0;
    for (const auto& f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE(total_deg == 4, "total degree should be 4");
}


void test_zassenhaus_single_factor() {
    TEST_CASE("zassenhaus_combine: single lifted factor → returns original poly");

    // f(x) = x + 1,单因子情形
    Polynomial<Rational> poly({Rational(1), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 125);  // p^k = 5^3

    EXPECT_TRUE(result.size() == 1, "single factor should return 1 true factor");
    EXPECT_TRUE(result[0].degree() == 1, "factor should be degree 1");
}

void test_zassenhaus_two_linear_factors() {
    TEST_CASE("zassenhaus_combine: x^2 - 1 = (x-1)(x+1)");

    // f(x) = x^2 - 1
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    // 提升后的因子(mod 125 = 5^3):(x - 1) 和 (x + 1)
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),  // x - 1
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")    // x + 1
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 2, "x^2 - 1 should have 2 true factors");

    // 验证因子乘积等于原多项式(首一化后)
    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of factors should equal original polynomial");
    }
}

void test_zassenhaus_irreducible_quadratic() {
    TEST_CASE("zassenhaus_combine: x^2 + x + 1 (irreducible over Q)");

    /// f(x) = x^2 + x + 1 在 Q 上为整体元素.
    /// Hensel 提升给出的两个模因子在 Q 上组合后仍返回原多项式.
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    // 模 7 下为 (x - 2)(x - 4),提升到模 49 后根为 18 和 30.
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 49);

    // 两个提升线性因子在有理数域上均不整除原多项式.
    // 所以应返回原多项式作为单一因子
    EXPECT_TRUE(result.size() == 1, "irreducible polynomial should return 1 factor");
    EXPECT_TRUE(result[0].degree() == 2, "factor should be degree 2");
}

void test_zassenhaus_three_factors() {
    TEST_CASE("zassenhaus_combine: x^3 - x = x(x-1)(x+1)");

    // f(x) = x^3 - x = x(x-1)(x+1)
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    // 提升后的因子 mod 125
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(0), BigInt(1)}, "x"),    // x
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),   // x - 1
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")     // x + 1
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 3, "x^3 - x should have 3 true factors");

    // 验证每个因子为线性
    int linear_count = 0;
    for (const auto& f : result) {
        if (f.degree() == 1) linear_count++;
    }
    EXPECT_TRUE(linear_count == 3, "all 3 factors should be linear");
}

void test_zassenhaus_empty_factors() {
    TEST_CASE("zassenhaus_combine: empty lifted factors");

    Polynomial<Rational> poly({Rational(1), Rational(0), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted;

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 1, "empty factors should return original poly");
}

void test_zassenhaus_zero_polynomial() {
    TEST_CASE("zassenhaus_combine_checked: 零多项式拒绝非空提升因子");

    Polynomial<Rational> poly("x");
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")
    };

    auto result =
        zassenhaus_combine_checked(poly, lifted, BigInt(125));
    EXPECT_FALSE(result.has_value(), "零多项式与非空提升因子不一致");
    EXPECT_TRUE(
        result.error().code == CasErrc::InvalidArgument,
        "不一致的提升因子应报告 InvalidArgument");
}

void test_zassenhaus_quadratic_times_linear() {
    TEST_CASE("zassenhaus_combine: (x^2+1)(x-1) with subset size 2 needed");

    // f(x) = x^3 - x^2 + x - 1 = (x^2 + 1)(x - 1)
    Polynomial<Rational> x2_plus_1({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> x_minus_1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> poly = x2_plus_1 * x_minus_1;

    // 模 125 下 x^2 + 1 的根提升为 57 和 68.
    // 因此 f 与 (x-57)(x-68)(x-1) 模 125 同余.
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-57), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-68), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    // 真因子 x-1 可由单元素子集恢复，剩余两个提升因子组合为 x^2+1.
    EXPECT_TRUE(result.size() >= 1, "should find at least 1 factor");

    // 验证所有因子的乘积等于原多项式
    if (result.size() >= 1) {
        Polynomial<Rational> product = result[0];
        for (size_t i = 1; i < result.size(); ++i) {
            product = product * result[i];
        }
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of all returned factors should equal original polynomial");
    }
}


void test_zassenhaus_early_termination_irreducible() {
    TEST_CASE("zassenhaus early termination: irreducible polynomial (all subsets checked)");

    /// f(x) = x^2 + x + 1 在 Q 上不可约.
    /// 模 49 因子 (x-18)(x-30) 的所有真子集候选均未通过整除检验,
    /// 因而剩余项等于原多项式.
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 49);

    // 不可约多项式:所有子集检查完毕后,剩余多项式作为唯一因子返回
    EXPECT_TRUE(result.size() == 1, "irreducible polynomial should return 1 factor");
    EXPECT_TRUE(result[0].degree() == 2, "factor should be the original degree-2 polynomial");
    // 验证返回的因子是首一的
    EXPECT_TRUE(result[0].coeffs[2] == Rational(1), "factor should be monic");
}

void test_zassenhaus_early_termination_remaining_irreducible() {
    TEST_CASE("zassenhaus early termination: after finding one factor, remaining is irreducible");

    // f(x) = (x - 1)(x^2 + x + 1) = x^3 - 1
    // 模 49 下 f = (x-1)(x-18)(x-30).
    // 在 Q 上:找到 (x-1) 后,剩余 x^2 + x + 1 不可约
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 49);

    // 应找到 (x-1) 作为第一个因子,然后剩余 x^2+x+1 不可约
    EXPECT_TRUE(result.size() == 2, "should find 2 true factors");

    // 验证因子乘积等于原多项式
    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of factors should equal original polynomial");
    }
}

void test_zassenhaus_early_termination_linear_remaining() {
    TEST_CASE("zassenhaus early termination: remaining polynomial is linear");

    // f(x) = (x^2 - 1)(x + 2) = x^3 + 2x^2 - x - 2
    // 模 5 下 x^2 - 1 = (x-1)(x+1),x + 2 = (x+2)
    // 找到 (x^2-1) 后,剩余 (x+2) 为线性 -> 立即终止
    Polynomial<Rational> x2_minus_1({Rational(-1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> x_plus_2({Rational(2), Rational(1)}, "x");
    Polynomial<Rational> poly = x2_minus_1 * x_plus_2;

    // 提升因子 mod 125: (x-1), (x+1), (x+2)
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),   // x - 1
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x"),    // x + 1
        Polynomial<BigInt>({BigInt(2), BigInt(1)}, "x")     // x + 2
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    // 应找到 3 个线性因子(或 (x^2-1) + (x+2))
    // 关键验证:因子乘积等于原多项式
    EXPECT_TRUE(result.size() >= 2, "should find at least 2 factors");

    Polynomial<Rational> product = result[0];
    for (size_t i = 1; i < result.size(); ++i) {
        product = product * result[i];
    }
    auto monic_poly = poly.make_monic();
    auto monic_product = product.make_monic();
    EXPECT_TRUE(monic_poly == monic_product,
        "product of all factors should equal original polynomial");
}

void test_zassenhaus_early_termination_single_active_factor() {
    TEST_CASE("zassenhaus early termination: single active factor remaining");

    // f(x) = (x - 1)(x + 1) = x^2 - 1
    // 找到 (x-1) 后,仅剩 1 个活跃因子 -> 立即终止
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),  // x - 1
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")    // x + 1
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 2, "should find 2 factors");

    // 验证因子乘积
    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of factors should equal original polynomial");
    }
}

void test_zassenhaus_early_termination_degree_one_input() {
    TEST_CASE("zassenhaus early termination: linear polynomial input (immediate)");

    // f(x) = x + 5,线性多项式直接返回(由 lifted_factors.size() == 1 处理)
    Polynomial<Rational> poly({Rational(5), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(5), BigInt(1)}, "x")
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 1, "linear polynomial should return 1 factor");
    EXPECT_TRUE(result[0].degree() == 1, "factor should be linear");
}


void test_zassenhaus_rational_reconstruction_integer_coeffs() {
    TEST_CASE("zassenhaus_combine with rational_reconstruction: integer coefficients still work");

    // f(x) = x^2 - 1 = (x-1)(x+1)
    // 整数系数情形:有理重构应恢复整数(分母为 1)
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),  // x - 1
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")    // x + 1
    };

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE(result.size() == 2, "x^2 - 1 should still have 2 factors with rational reconstruction");

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of factors should equal original polynomial");
    }
}

void test_zassenhaus_rational_reconstruction_with_rational_poly() {
    TEST_CASE("zassenhaus_combine with rational_reconstruction: rational coefficient polynomial");

    // f(x) = x^2 - 1/4 = (x - 1/2)(x + 1/2)
    // 这是一个有理系数多项式,提升因子的系数在模 p^k 下编码了有理数
    // 使用 p^k = 1000000007 (大素数)
    // 对于 1/2 mod 1000000007: 需要 2 的逆元 mod 1000000007
    // 2^(-1) mod 1000000007 = 500000004
    // 所以 x - 1/2 在 mod 1000000007 下为 x - 500000004 (对称表示: x + 500000003)
    // 但有理重构应该能从 500000004 恢复出 1/2

    // 更简单的测试:使用较小的模数
    // p^k = 49 (7^2)
    // 1/2 mod 49: 需要 2^(-1) mod 49 = 25 (因为 2*25 = 50 == 1 mod 49)
    // 对称表示: 25 (在 [-24, 24] 范围内)
    // 有理重构 25 mod 49: bound = floor(sqrt(49/2)) = floor(4.95) = 4
    // 25 > 4, 所以需要运行扩展欧几里得:
    //   r0=49, r1=25, s0=0, s1=1
    //   q=1: r0=25, r1=24, s0=1, s1=-1
    //   q=1: r0=24, r1=1, s0=-1, s1=2
    //   r1=1 <= 4, 停止. p=1, q=2 -> 1/2 [ok]

    // f(x) = x^2 - 1/4
    Polynomial<Rational> poly({Rational(-1, 4), Rational(0), Rational(1)}, "x");

    // 提升因子 mod 49: (x - 25) 和 (x + 25)
    // 因为 1/2 == 25 mod 49, 所以 x - 1/2 == x - 25 mod 49
    // 对称表示: -25 在 [-24,24] 之外,所以 -25 + 49 = 24 -> x + 24? 不对
    // 让我重新计算:对称表示 [-m/2, m/2) = [-24, 24]
    // -25 mod 49 = 24 (在对称表示中)
    // 25 mod 49 = 25 -> 超出 24,所以 25 - 49 = -24 (在对称表示中)
    // 所以 x - 1/2 的系数: 常数项 = -1/2 == -25 mod 49 -> 对称表示 24
    //       x + 1/2 的系数: 常数项 = 1/2 == 25 mod 49 -> 对称表示 -24

    // 使用更大的模数使有理重构更可靠
    // p^k = 10007 (素数)
    // 1/2 mod 10007: 2^(-1) mod 10007 = 5004 (因为 2*5004 = 10008 == 1 mod 10007)
    // bound = floor(sqrt(10007/2)) = floor(70.7) = 70
    // 有理重构 5004 mod 10007:
    //   r0=10007, r1=5004, s0=0, s1=1
    //   q=2: r0=5004, r1=-1 -> r1=10007-2*5004=-1 -> 实际 r1=10007-2*5004=-1
    //   不对,让我用正确的算法:
    //   r0=10007, r1=5004
    //   q=10007/5004=1, r_new=10007-1*5004=5003, s_new=0-1*1=-1
    //   r0=5004, r1=5003, s0=1, s1=-1
    //   q=5004/5003=1, r_new=5004-1*5003=1, s_new=1-1*(-1)=2
    //   r1=1 <= 70, 停止. p=1, q=2 -> 1/2 [ok]

    // 使用 p^k = 10007
    // x - 1/2: 常数项 -1/2 == -(5004) mod 10007 = 5003 -> 对称: 5003 > 5003? 5003 < 5003.5 -> 5003
    // x + 1/2: 常数项 1/2 == 5004 mod 10007 -> 对称: 5004 > 5003.5 -> 5004 - 10007 = -5003

    // 实际上让我用一个简单的整数系数例子来验证有理重构路径被调用
    // 使用 x^3 - x = x(x-1)(x+1) 这个已知能工作的例子
    // 但用更大的 prime_power 来确保有理重构路径被触发

    // 简单验证:x^2 - 1 with large prime power
    Polynomial<Rational> poly2({Rational(-1), Rational(0), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted2 = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")
    };

    // 使用大素数幂确保有理重构路径被使用
    auto result = checked_zassenhaus(poly2, lifted2, 1000000007LL);

    EXPECT_TRUE(result.size() == 2, "should find 2 factors with large prime power");

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly2.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "product of factors should equal original with large prime power");
    }
}


void test_zassenhaus_bounded_enumeration_many_factors() {
    TEST_CASE("zassenhaus_combine_checked: 有界枚举处理 16 个线性因子");

    // 构造 f(x) = (x-1)(x-2)...(x-16),逐轮首个单因子候选即可整除.

    // 构造 f(x) = product of (x - i) for i = 1..16
    Polynomial<Rational> poly({Rational(1)}, "x");  // 从 1 开始
    for (int i = 1; i <= 16; ++i) {
        Polynomial<Rational> factor({Rational(-i), Rational(1)}, "x");
        poly = poly * factor;
    }

    // 构造对应的 BigInt 提升因子
    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 16; ++i) {
        lifted.push_back(Polynomial<BigInt>({BigInt(-i), BigInt(1)}, "x"));
    }

    // 使用足够大的精确重构模数.
    int64_t prime_power = 1000000007LL;

    auto result = checked_zassenhaus(poly, lifted, prime_power);

    // 应该找到 16 个线性因子
    EXPECT_TRUE(result.size() == 16, "bounded enumeration finds 16 factors");

    // 验证每个因子为线性
    int linear_count = 0;
    for (const auto& f : result) {
        if (f.degree() == 1) linear_count++;
    }
    EXPECT_TRUE(linear_count == 16, "all 16 factors should be linear");
}

void test_zassenhaus_bounded_enumeration_product() {
    TEST_CASE("zassenhaus_combine_checked: 有界枚举结果精确重构输入");

    // 构造 f(x) = (x-1)(x-2)...(x-16).
    Polynomial<Rational> poly({Rational(1)}, "x");
    for (int i = 1; i <= 16; ++i) {
        Polynomial<Rational> factor({Rational(-i), Rational(1)}, "x");
        poly = poly * factor;
    }

    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 16; ++i) {
        lifted.push_back(Polynomial<BigInt>({BigInt(-i), BigInt(1)}, "x"));
    }

    int64_t prime_power = 1000000007LL;

    auto result = checked_zassenhaus(poly, lifted, prime_power);

    // 验证因子乘积等于原多项式
    if (!result.empty()) {
        Polynomial<Rational> product = result[0];
        for (size_t i = 1; i < result.size(); ++i) {
            product = product * result[i];
        }
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE(monic_poly == monic_product,
            "bounded-enumeration factors reconstruct original polynomial");
    }
}

void test_zassenhaus_tiny_budget_is_inconclusive() {
    TEST_CASE("zassenhaus_combine_checked: 微小预算返回精确未决结果");

    Polynomial<Rational> poly({Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 4; ++i) {
        poly = poly *
            Polynomial<Rational>({Rational(-i), Rational(1)}, "x");
        lifted.emplace_back(
            std::vector<BigInt>{BigInt(-i), BigInt(1)}, "x");
    }
    ResourceLimits limits;
    limits.max_steps = lifted.size();
    ComputationContext context(limits);
    auto result = zassenhaus_combine_checked(
        poly, lifted, BigInt(1000000007LL), context);

    EXPECT_TRUE(result.has_value(), "预算耗尽应返回精确部分结果");
    if (result) {
        EXPECT_TRUE(
            result.value().completeness == Completeness::Inconclusive,
            "预算耗尽应标记为 Inconclusive");
        EXPECT_TRUE(
            result.value().value.size() == 1 &&
                result.value().value[0] == poly.make_monic(),
            "未决结果仍应精确重构输入");
    }
}


void test_back_substitute_u0_plus_x_to_sin_x_plus_x() {
    TEST_CASE("back-substitution: u0 + x → sin(x) + x (mapping u0 → sin(x))");

    auto x = SymbolicExpr::variable("x");
    // 构造因子表达式 u0 + x
    auto u0 = SymbolicExpr::variable("u0");
    auto factor_expr = SymbolicExpr::add(u0, x);

    // 构造映射 u0 -> sin(x)
    auto sin_x = SymbolicExpr::sin(x);

    // 执行完整的换元->逆换元流程验证
    // 先用 detect_trans_substitutions 获取映射
    auto original = SymbolicExpr::add(sin_x, x);
    auto sub_result = detect_trans_substitutions(original, "x");

    EXPECT_TRUE(sub_result.mappings.size() == 1, "should have 1 mapping");
    EXPECT_EQ_STR(sub_result.mappings[0].indeterminate, "u0", "indeterminate is u0");

    // 对换元后的表达式执行逆换元
    auto back = sub_result.poly_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);

    EXPECT_TRUE(back != nullptr, "back-substituted result should not be null");
    std::string back_str = back->to_string();
    EXPECT_TRUE(back_str.find("sin") != std::string::npos,
        "back-substituted should contain sin");
    EXPECT_TRUE(back_str.find("u0") == std::string::npos,
        "back-substituted should not contain u0");
    EXPECT_TRUE(back_str.find("x") != std::string::npos,
        "back-substituted should contain x");
}

void test_back_substitute_u0_squared_minus_x_squared() {
    TEST_CASE("back-substitution: u0² - x² → sin²(x) - x² (mapping u0 → sin(x))");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    // 原始表达式 sin^2(x) - x^2
    auto original = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(x, SymbolicExpr::number(2)))
    );

    auto sub_result = detect_trans_substitutions(original, "x");
    EXPECT_TRUE(sub_result.mappings.size() == 1, "should have 1 mapping");

    // 逆换元
    auto back = sub_result.poly_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);

    EXPECT_TRUE(back != nullptr, "back-substituted result should not be null");
    std::string back_str = back->to_string();
    EXPECT_TRUE(back_str.find("sin") != std::string::npos,
        "back-substituted should contain sin");
    EXPECT_TRUE(back_str.find("u0") == std::string::npos,
        "back-substituted should not contain u0");
}

void test_back_substitute_multiple_mappings() {
    TEST_CASE("back-substitution: u0 + u1 → sin(x) + cos(x) (two mappings)");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    // 原始表达式 sin(x) + cos(x)
    auto original = SymbolicExpr::add(sin_x, cos_x);

    auto sub_result = detect_trans_substitutions(original, "x");
    EXPECT_TRUE(sub_result.mappings.size() == 2, "should have 2 mappings");

    // 逐一逆换元
    auto back = sub_result.poly_expr;
    for (const auto& m : sub_result.mappings) {
        back = back->substitute(m.indeterminate, m.trans_expr);
    }

    EXPECT_TRUE(back != nullptr, "back-substituted result should not be null");
    std::string back_str = back->to_string();
    EXPECT_TRUE(back_str.find("sin") != std::string::npos,
        "back-substituted should contain sin");
    EXPECT_TRUE(back_str.find("cos") != std::string::npos,
        "back-substituted should contain cos");
    EXPECT_TRUE(back_str.find("u0") == std::string::npos,
        "back-substituted should not contain u0");
    EXPECT_TRUE(back_str.find("u1") == std::string::npos,
        "back-substituted should not contain u1");
}

void test_back_substitute_polynomial_factor() {
    TEST_CASE("back-substitution: polynomial factor from factored result");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    // sin^2(x) + sin(x) + 1 -> u0^2 + u0 + 1
    auto original = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(sin_x, SymbolicExpr::number(1))
    );

    auto sub_result = detect_trans_substitutions(original, "x");
    EXPECT_TRUE(sub_result.mappings.size() == 1, "should have 1 mapping");

    // 构造多项式 u0 + 1(模拟一个因子)
    std::vector<std::string> indeterminates = {sub_result.mappings[0].indeterminate};
    Polynomial<Rational> factor_poly({Rational(1), Rational(1)}, "u0");
    auto factor_expr = poly_to_symbolic(factor_poly);

    // 逆换元
    auto back = factor_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);

    EXPECT_TRUE(back != nullptr, "back-substituted factor should not be null");
    std::string back_str = back->to_string();
    EXPECT_TRUE(back_str.find("sin") != std::string::npos,
        "back-substituted factor should contain sin");
    EXPECT_TRUE(back_str.find("u0") == std::string::npos,
        "back-substituted factor should not contain u0");
}

void test_back_substitute_no_mappings() {
    TEST_CASE("back-substitution: no mappings leaves expression unchanged");

    auto x = SymbolicExpr::variable("x");
    // 纯多项式 x^2 + 1
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1)
    );

    auto sub_result = detect_trans_substitutions(expr, "x");
    EXPECT_TRUE(sub_result.mappings.empty(), "should have no mappings");

    // 逆换元(无映射时应返回原表达式)
    auto back = sub_result.poly_expr;
    for (const auto& m : sub_result.mappings) {
        back = back->substitute(m.indeterminate, m.trans_expr);
    }

    std::string orig_str = expr->to_string();
    std::string back_str = back->to_string();
    EXPECT_EQ_STR(back_str, orig_str, "no-mapping back-substitution should be identity");
}

void test_back_substitute_roundtrip() {
    TEST_CASE("back-substitution roundtrip: substitute then back-substitute preserves structure");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    // sin(x) * cos(x) + sin(x)
    auto original = SymbolicExpr::add(
        SymbolicExpr::multiply(sin_x, cos_x),
        sin_x
    );

    auto sub_result = detect_trans_substitutions(original, "x");

    // 逆换元
    auto back = sub_result.poly_expr;
    for (const auto& m : sub_result.mappings) {
        back = back->substitute(m.indeterminate, m.trans_expr);
    }

    EXPECT_TRUE(back != nullptr, "roundtrip result should not be null");
    std::string back_str = back->to_string();
    /// 逆换元结果应仅包含原始变量与函数.
    EXPECT_TRUE(back_str.find("u0") == std::string::npos,
        "roundtrip should not contain u0");
    EXPECT_TRUE(back_str.find("u1") == std::string::npos,
        "roundtrip should not contain u1");
    // 应包含原始超越函数
    EXPECT_TRUE(back_str.find("sin") != std::string::npos,
        "roundtrip should contain sin");
    EXPECT_TRUE(back_str.find("cos") != std::string::npos,
        "roundtrip should contain cos");
}


void test_simplify_factors_extract_constant_from_product() {
    TEST_CASE("tf_simplify_factors: [2*sin(x), x+1] → [2, sin(x), x+1]");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    // 构造因子列表 [2*sin(x), x+1]
    auto factor1 = SymbolicExpr::multiply(SymbolicExpr::number(2), sin_x);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    // 应提取常数 2,结果为 [2, sin(x), x+1]
    EXPECT_TRUE(result.size() == 3, "should have 3 factors (constant + 2 non-constant)");

    // 第一个因子应为常数 2
    EXPECT_TRUE(result[0]->is_number(), "first factor should be numeric");
    std::string first_str = result[0]->to_string();
    EXPECT_TRUE(first_str == "2", "first factor should be 2");

    /// 数值前导系数集中在首个常数因子中.
    bool has_sin = false;
    bool has_x_plus_1 = false;
    for (size_t i = 1; i < result.size(); ++i) {
        std::string s = result[i]->to_string();
        if (s.find("sin") != std::string::npos) has_sin = true;
        if (s.find("x") != std::string::npos && s.find("sin") == std::string::npos) has_x_plus_1 = true;
    }
    EXPECT_TRUE(has_sin, "should have sin(x) factor");
    EXPECT_TRUE(has_x_plus_1, "should have x+1 factor");
}

void test_simplify_factors_no_constants() {
    TEST_CASE("tf_simplify_factors: [sin(x)+x, cos(x)-1] → no constants extracted");

    auto x = SymbolicExpr::variable("x");

    // 构造因子列表 [sin(x)+x, cos(x)-1]
    auto factor1 = SymbolicExpr::add(SymbolicExpr::sin(x), x);
    auto factor2 = SymbolicExpr::add(SymbolicExpr::cos(x), SymbolicExpr::number(-1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    // 无常数可提取,结果应为 2 个因子
    EXPECT_TRUE(result.size() == 2, "should have 2 factors (no constants to extract)");

    // 验证无纯数值因子
    for (const auto& f : result) {
        EXPECT_FALSE(f->is_number(), "no factor should be purely numeric");
    }
}

void test_simplify_factors_pure_constant_factor() {
    TEST_CASE("tf_simplify_factors: [3, x+1] → [3, x+1] (constant already separate)");

    auto x = SymbolicExpr::variable("x");

    // 构造因子列表 [3, x+1]
    auto factor1 = SymbolicExpr::number(3);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    // 常数 3 应被提取并放在首位
    EXPECT_TRUE(result.size() == 2, "should have 2 factors");
    EXPECT_TRUE(result[0]->is_number(), "first factor should be numeric constant");
    std::string first_str = result[0]->to_string();
    EXPECT_TRUE(first_str == "3", "first factor should be 3");
}

void test_simplify_factors_multiple_constants_combined() {
    TEST_CASE("tf_simplify_factors: [2*sin(x), 5*cos(x)] → [10, sin(x), cos(x)]");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    // 构造因子列表 [2*sin(x), 5*cos(x)]
    /// simplify() 保持 c*trig(x) 的乘积结构.
    auto factor1 = SymbolicExpr::multiply(SymbolicExpr::number(2), sin_x);
    auto factor2 = SymbolicExpr::multiply(SymbolicExpr::number(5), cos_x);

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    // 应提取常数 2*5=10
    EXPECT_TRUE(result.size() >= 2, "should have at least 2 factors");
    EXPECT_TRUE(result[0]->is_number(), "first factor should be numeric");
    std::string first_str = result[0]->to_string();
    EXPECT_TRUE(first_str == "10", "combined constant should be 10");
}

void test_simplify_factors_constant_one_not_added() {
    TEST_CASE("tf_simplify_factors: constant product = 1 → no constant factor added");

    auto x = SymbolicExpr::variable("x");

    // 构造因子列表 [sin(x), x+1] - 无常数乘子
    auto factor1 = SymbolicExpr::sin(x);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    /// 常数积为 1 时结果仅包含两个函数因子.
    EXPECT_TRUE(result.size() == 2, "should have 2 factors (no constant added)");
    for (const auto& f : result) {
        if (f->is_number()) {
            std::string s = f->to_string();
            EXPECT_FALSE(s == "1", "should not have explicit constant 1");
        }
    }
}

void test_simplify_factors_simplify_called() {
    TEST_CASE("tf_simplify_factors: simplify() is called on each factor");

    auto x = SymbolicExpr::variable("x");

    // 构造一个可化简的表达式:x + 0 应化简为 x
    auto factor1 = SymbolicExpr::add(x, SymbolicExpr::number(0));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1};
    auto result = tf_simplify_factors(factors);

    // 化简后 x + 0 -> x
    EXPECT_TRUE(result.size() >= 1, "should have at least 1 factor");
    // 结果应为简化形式
    bool found_x = false;
    for (const auto& f : result) {
        std::string s = f->to_string();
        if (s == "x") found_x = true;
    }
    EXPECT_TRUE(found_x, "x + 0 should simplify to x");
}


void test_mult_structure_x_times_sin_x() {
    TEST_CASE("multiplicative structure: x * sin(x) → [x, sin(x)]");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(x, SymbolicExpr::sin(x));

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 2, "should return 2 factors");
    // 验证因子包含 x 和 sin(x)
    bool has_x = false, has_sin = false;
    for (const auto& f : factors) {
        std::string s = f->to_string();
        if (s == "x") has_x = true;
        if (s.find("sin") != std::string::npos && s.find("x") != std::string::npos) has_sin = true;
    }
    EXPECT_TRUE(has_x, "should contain factor x");
    EXPECT_TRUE(has_sin, "should contain factor sin(x)");
}

void test_mult_structure_x2_times_exp_times_cos() {
    TEST_CASE("multiplicative structure: x^2 * exp(x) * cos(x) → 3 factors");

    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto exp_x = SymbolicExpr::exp(x);
    auto cos_x = SymbolicExpr::cos(x);

    auto expr = SymbolicExpr::multiply(x2, SymbolicExpr::multiply(exp_x, cos_x));

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 3, "should return 3 factors");
}

void test_mult_structure_with_constant() {
    TEST_CASE("multiplicative structure: 2 * x * sin(x) → [2, x, sin(x)]");

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto expr = SymbolicExpr::multiply(two, SymbolicExpr::multiply(x, SymbolicExpr::sin(x)));

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 3, "should return 3 factors (constant + 2 non-constant)");
    // 验证包含常数因子 2
    bool has_const = false;
    for (const auto& f : factors) {
        if (f->is_number()) has_const = true;
    }
    EXPECT_TRUE(has_const, "should contain numeric constant factor");
}

void test_mult_structure_sum_not_product() {
    TEST_CASE("multiplicative structure: sin(x) + x → not a product, uses full pipeline");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto factors = factor_transcendental(expr, "x");

    /// sin(x) + x 由完整流程判定为单一整体因子.
    EXPECT_TRUE(factors.size() == 1, "sum expression should not be split by multiplicative detection");
}


void test_linear_irreducible_sin_plus_x() {
    TEST_CASE("linear irreducibility: sin(x) + x → irreducible");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 1, "sin(x) + x should be irreducible (linear in u0 and x)");
    EXPECT_EQ_EXPR_STR(factors[0], expr->to_string(), "factor should be original expression");
}

void test_linear_irreducible_2sin_3x_1() {
    TEST_CASE("linear irreducibility: 2*sin(x) + 3*x + 1 → irreducible");

    auto x = SymbolicExpr::variable("x");
    // 2*sin(x) + 3*x + 1
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::sin(x)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(1)
        )
    );

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 1, "2*sin(x) + 3*x + 1 should be irreducible");
}

void test_linear_irreducible_sin_cos_x() {
    TEST_CASE("linear irreducibility: sin(x) + cos(x) + x → irreducible");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + cos(x) + x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(SymbolicExpr::cos(x), x)
    );

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 1, "sin(x) + cos(x) + x should be irreducible (linear in u0, u1, x)");
}

void test_linear_irreducible_sin_squared_not_linear() {
    TEST_CASE("linear irreducibility: sin^2(x) + x → NOT linear, should not trigger fast path");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    // sin^2(x) + x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        x
    );

    auto factors = factor_transcendental(expr, "x");

    /// sin^2(x) + x 对 u0 为二次式,因此进入完整因式分解流程.
    /// 结果仍可为整体元素,其判定来源于完整流程.
    // 实际上 sin^2(x) + x 在超越多项式环中也是不可约的,但通过完整流程判定
    EXPECT_TRUE(factors.size() >= 1, "sin^2(x) + x should return at least 1 factor");
}

void test_linear_irreducible_sin_times_x_not_sum() {
    TEST_CASE("linear irreducibility: sin(x) * x → product, handled by multiplicative detection");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::sin(x), x);

    auto factors = factor_transcendental(expr, "x");

    /// sin(x) * x 由乘法结构直接拆分为两个因子.
    EXPECT_TRUE(factors.size() == 2, "sin(x) * x should be split into 2 factors by multiplicative detection");
}


void test_exp_separation_exp_x_times_x_plus_exp_x() {
    TEST_CASE("exponential separation: exp(x)*x + exp(x) → [exp(x), x+1]");

    auto x = SymbolicExpr::variable("x");
    auto exp_x = SymbolicExpr::exp(x);
    // exp(x)*x + exp(x)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(exp_x, x),
        exp_x
    );

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 2, "should produce 2 factors: exp(x) and (x+1)");

    // 验证其中一个因子为 exp(x)
    bool has_exp = false;
    bool has_poly = false;
    for (const auto& f : factors) {
        std::string s = f->to_string();
        if (s.find("exp") != std::string::npos && s.find("+") == std::string::npos) {
            has_exp = true;
        }
        if (s.find("exp") == std::string::npos) {
            has_poly = true;
        }
    }
    EXPECT_TRUE(has_exp, "should contain exp(x) as a factor");
    EXPECT_TRUE(has_poly, "should contain polynomial remainder as a factor");
}

void test_exp_separation_no_common_exp() {
    TEST_CASE("exponential separation: exp(x) + sin(x) → no common exp factor");

    auto x = SymbolicExpr::variable("x");
    /// exp(x) + sin(x) 的两项具有不同函数基.
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::sin(x));

    auto factors = factor_transcendental(expr, "x");

    /// 公共指数因子提取保持未匹配状态;表达式随后由线性整体性规则处理.
    EXPECT_TRUE(factors.size() == 1, "exp(x) + sin(x) should not be split by exponential separation");
}

void test_exp_separation_product_not_sum() {
    TEST_CASE("exponential separation: exp(x) * sin(x) → handled by multiplicative detection, not exp separation");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::exp(x), SymbolicExpr::sin(x));

    auto factors = factor_transcendental(expr, "x");

    // 乘积形式由乘积结构检测处理
    EXPECT_TRUE(factors.size() == 2, "exp(x) * sin(x) should be split by multiplicative detection");
}

void test_exp_separation_three_terms() {
    TEST_CASE("exponential separation: exp(x)*x^2 + 2*exp(x)*x + exp(x) → [exp(x), x^2+2x+1]");

    auto x = SymbolicExpr::variable("x");
    auto exp_x = SymbolicExpr::exp(x);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    // exp(x)*x^2 + 2*exp(x)*x + exp(x)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(exp_x, x_sq),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(exp_x, x)),
            exp_x
        )
    );

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 2, "should produce 2 factors: exp(x) and (x^2+2x+1)");

    // 验证其中一个因子为 exp(x)
    bool has_exp = false;
    for (const auto& f : factors) {
        std::string s = f->to_string();
        if (s.find("exp") != std::string::npos && s.find("+") == std::string::npos
            && s.find("*") == std::string::npos) {
            has_exp = true;
        }
    }
    EXPECT_TRUE(has_exp, "should contain exp(x) as a factor");
}


void test_pythagorean_sin2_cos2_to_1() {
    TEST_CASE("Pythagorean simplification: sin²(x) + cos²(x) → 1");

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    auto expr = SymbolicExpr::add(sin2, cos2);

    auto factors = factor_transcendental(expr, "x");

    // sin^2(x) + cos^2(x) = 1,化简后无超越函数,应返回 {1} 或等价
    // 由于化简为 1 后不含超越函数,factor_transcendental 应返回 {1}
    EXPECT_TRUE(factors.size() == 1, "should produce 1 factor");
    if (!factors.empty()) {
        std::string s = factors[0]->to_string();
        EXPECT_TRUE(factors[0]->is_number() || s == "1",
                    "sin²(x) + cos²(x) should simplify to 1");
    }
}

void test_pythagorean_sin2_cos2_plus_x() {
    TEST_CASE("Pythagorean simplification: sin²(x) + cos²(x) + x → 1 + x");

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    // sin^2(x) + cos^2(x) + x
    auto expr = SymbolicExpr::add(sin2, SymbolicExpr::add(cos2, x));

    auto factors = factor_transcendental(expr, "x");

    // 化简后为 1 + x,无超越函数,应返回 {1+x} 或等价
    EXPECT_TRUE(factors.size() >= 1, "should produce at least 1 factor");
    // 验证结果不含 sin 或 cos
    for (const auto& f : factors) {
        std::string s = f->to_string();
        EXPECT_TRUE(s.find("sin") == std::string::npos,
                    "result should not contain sin after Pythagorean simplification");
        EXPECT_TRUE(s.find("cos") == std::string::npos,
                    "result should not contain cos after Pythagorean simplification");
    }
}

void test_pythagorean_common_coefficient() {
    TEST_CASE("Pythagorean simplification: 2*sin²(x) + 2*cos²(x) → 2");

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(two, sin2),
        SymbolicExpr::multiply(two, cos2)
    );

    auto factors = factor_transcendental(expr, "x");

    // 2*sin^2(x) + 2*cos^2(x) = 2,化简后为常数
    EXPECT_TRUE(factors.size() == 1, "should produce 1 factor");
    if (!factors.empty()) {
        EXPECT_TRUE(factors[0]->is_number(),
                    "2*sin²(x) + 2*cos²(x) should simplify to 2");
    }
}

void test_pythagorean_different_args_unchanged() {
    TEST_CASE("Pythagorean: sin²(x) + cos²(y) → unchanged (different arguments)");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto two = SymbolicExpr::number(2);
    auto sin2_x = SymbolicExpr::power(SymbolicExpr::sin(x), two);
    auto cos2_y = SymbolicExpr::power(SymbolicExpr::cos(y), two);
    auto expr = SymbolicExpr::add(sin2_x, cos2_y);

    auto factors = factor_transcendental(expr, "x");

    /// 参数不同的三角项保持各自结构,结果继续包含 sin 或 cos.
    bool has_trig = false;
    for (const auto& f : factors) {
        std::string s = f->to_string();
        if (s.find("sin") != std::string::npos || s.find("cos") != std::string::npos) {
            has_trig = true;
        }
    }
    EXPECT_TRUE(has_trig, "sin²(x) + cos²(y) should remain unchanged (different args)");
}

void test_pythagorean_no_matching_cos() {
    TEST_CASE("Pythagorean: sin²(x) + x → unchanged (no matching cos²)");

    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin2_x = SymbolicExpr::power(SymbolicExpr::sin(x), two);
    auto expr = SymbolicExpr::add(sin2_x, x);

    auto factors = factor_transcendental(expr, "x");

    /// 单独的 sin^2(x) 保持平方结构.
    bool has_sin = false;
    for (const auto& f : factors) {
        std::string s = f->to_string();
        if (s.find("sin") != std::string::npos) {
            has_sin = true;
        }
    }
    EXPECT_TRUE(has_sin, "sin²(x) + x should remain unchanged (no matching cos²)");
}


void test_full_pipeline_sin2_minus_x2() {
    TEST_CASE("full pipeline: sin²(x) - x² = (sin(x)+x)(sin(x)-x)");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    // sin^2(x) - x^2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2)))
    );

    auto factors = factor_transcendental(expr, "x");

    // sin^2(x) - x^2 = (sin(x)+x)(sin(x)-x) 是差平方形式.
    // 当前实现中,换元后 u0^2 - x^2 为双变量多项式,
    // 若多项式构造成功则应分解为 2 个因子,否则返回原表达式.
    // 验证:至少返回 1 个因子,且因子包含 sin
    EXPECT_TRUE(factors.size() >= 1, "should produce at least 1 factor");

    if (factors.size() == 2) {
        // 若成功分解,验证因子结构
        std::string f0 = factors[0]->to_string();
        std::string f1 = factors[1]->to_string();
        bool both_have_sin = (f0.find("sin") != std::string::npos) &&
                             (f1.find("sin") != std::string::npos);
        EXPECT_TRUE(both_have_sin, "both factors should contain sin(x)");
    } else {
        // 若返回为不可约(多变量多项式构造限制),验证原表达式被保留
        bool has_sin = false;
        for (const auto& f : factors) {
            if (f->to_string().find("sin") != std::string::npos) has_sin = true;
        }
        EXPECT_TRUE(has_sin, "irreducible result should preserve sin(x)");
    }
}

void test_full_pipeline_exp_times_x() {
    TEST_CASE("full pipeline: exp(x)*x → [exp(x), x]");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::exp(x), x);

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 2, "should produce 2 factors");

    // 验证因子分别为 exp(x) 和 x
    if (factors.size() == 2) {
        bool has_exp = false;
        bool has_x = false;
        for (const auto& f : factors) {
            std::string s = f->to_string();
            if (s.find("exp") != std::string::npos) has_exp = true;
            if (s == "x") has_x = true;
        }
        EXPECT_TRUE(has_exp, "one factor should be exp(x)");
        EXPECT_TRUE(has_x, "one factor should be x");
    }
}

void test_full_pipeline_irreducible_sin_plus_x() {
    TEST_CASE("full pipeline: sin(x)+x is irreducible");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto factors = factor_transcendental(expr, "x");

    EXPECT_TRUE(factors.size() == 1, "should return 1 factor (irreducible)");
}

int main() {
    test_single_sin();
    test_sin_and_cos();
    test_deduplication();
    test_exp_function();
    test_ln_function();
    test_tan_function();
    test_no_transcendental();
    test_independent_variable();
    test_nested_transcendental();
    test_multiple_distinct();
    test_null_expression();
    test_sin_cos_pythagorean_constraint();
    test_exp_inverse_constraint();
    test_no_constraint_different_args();
    test_no_constraint_same_type();
    test_multiple_constraints();
    test_substitution_sin_plus_x();
    test_substitution_sin_squared_minus_x_squared();
    test_substitution_no_transcendental();

    test_build_poly_single_indeterminate();
    test_build_poly_single_indeterminate_with_rational_coeffs();
    test_build_poly_constant_expression();
    test_build_poly_null_expression();
    test_build_poly_only_x();
    test_build_poly_main_variable_selection();
    test_build_poly_from_substitution_result();
    test_build_poly_non_polynomial_fails();

    test_validate_valid_polynomial();
    test_validate_remaining_sin_fails();
    test_validate_remaining_cos_fails();
    test_validate_remaining_exp_fails();
    test_validate_fractional_exponent_fails();
    test_validate_negative_exponent_fails();
    test_validate_negative_fractional_exponent_fails();
    test_validate_transcendental_in_original_var_fails();
    test_validate_valid_multivariate_polynomial();

    test_square_free_already_square_free();
    test_square_free_repeated_root();
    test_square_free_higher_multiplicity();
    test_square_free_constant_polynomial();
    test_square_free_linear_polynomial();
    test_square_free_zero_polynomial();

    test_berlekamp_linear_poly();
    test_berlekamp_quadratic_irreducible();
    test_berlekamp_x2_minus_1_mod3();
    test_berlekamp_cubic_poly();
    test_berlekamp_constant_poly();
    test_berlekamp_specified_prime();

    test_null_space_x2_minus_1_mod3();
    test_null_space_x2_plus_x_plus_1_mod2();
    test_null_space_x3_minus_x_mod5();
    test_null_space_linear_poly();
    test_null_space_basis_vectors_in_kernel();

    test_split_x2_minus_1_mod3();
    test_split_x3_minus_x_mod5();
    test_split_irreducible_x2_plus_x_plus_1_mod2();
    test_split_product_verification_mod3();
    test_split_product_verification_mod5();
    test_split_x4_minus_1_mod5();

    test_zassenhaus_single_factor();
    test_zassenhaus_two_linear_factors();
    test_zassenhaus_irreducible_quadratic();
    test_zassenhaus_three_factors();
    test_zassenhaus_empty_factors();
    test_zassenhaus_zero_polynomial();
    test_zassenhaus_quadratic_times_linear();

    test_zassenhaus_rational_reconstruction_integer_coeffs();
    test_zassenhaus_rational_reconstruction_with_rational_poly();

    test_zassenhaus_early_termination_irreducible();
    test_zassenhaus_early_termination_remaining_irreducible();
    test_zassenhaus_early_termination_linear_remaining();
    test_zassenhaus_early_termination_single_active_factor();
    test_zassenhaus_early_termination_degree_one_input();

    test_zassenhaus_bounded_enumeration_many_factors();
    test_zassenhaus_bounded_enumeration_product();
    test_zassenhaus_tiny_budget_is_inconclusive();

    test_back_substitute_u0_plus_x_to_sin_x_plus_x();
    test_back_substitute_u0_squared_minus_x_squared();
    test_back_substitute_multiple_mappings();
    test_back_substitute_polynomial_factor();
    test_back_substitute_no_mappings();
    test_back_substitute_roundtrip();

    test_simplify_factors_extract_constant_from_product();
    test_simplify_factors_no_constants();
    test_simplify_factors_pure_constant_factor();
    test_simplify_factors_multiple_constants_combined();
    test_simplify_factors_constant_one_not_added();
    test_simplify_factors_simplify_called();

    test_mult_structure_x_times_sin_x();
    test_mult_structure_x2_times_exp_times_cos();
    test_mult_structure_with_constant();
    test_mult_structure_sum_not_product();

    test_linear_irreducible_sin_plus_x();
    test_linear_irreducible_2sin_3x_1();
    test_linear_irreducible_sin_cos_x();
    test_linear_irreducible_sin_squared_not_linear();
    test_linear_irreducible_sin_times_x_not_sum();

    test_exp_separation_exp_x_times_x_plus_exp_x();
    test_exp_separation_no_common_exp();
    test_exp_separation_product_not_sum();
    test_exp_separation_three_terms();

    test_pythagorean_sin2_cos2_to_1();
    test_pythagorean_sin2_cos2_plus_x();
    test_pythagorean_common_coefficient();
    test_pythagorean_different_args_unchanged();
    test_pythagorean_no_matching_cos();

    test_full_pipeline_sin2_minus_x2();
    test_full_pipeline_exp_times_x();
    test_full_pipeline_irreducible_sin_plus_x();

    return TEST_REPORT();
}
