#include "../include/lsr_expr.hpp"
#include "../include/assumption_context.hpp"
#include "test_common.hpp"

int main() {
    TEST_CASE("LSR symbols reject reserved mathematical constants");

    auto x = lamina::lsr::sym("x");
    EXPECT_TRUE(x.has_value(), "ordinary symbol can be created");

    auto reserved_i = lamina::lsr::sym("i");
    EXPECT_TRUE(!reserved_i &&
                    reserved_i.error().code == lamina::CasErrc::InvalidArgument,
                "i is the imaginary unit and cannot be shadowed as a symbol");
    EXPECT_TRUE(!reserved_i &&
                    std::string(lamina::lsr::error_name(reserved_i.error())) ==
                        "ImaginaryUnitReserved",
                "LSR diagnostic name preserves reserved imaginary unit errors");

    auto reserved_I = lamina::lsr::sym("I");
    EXPECT_TRUE(!reserved_I &&
                    reserved_I.error().code == lamina::CasErrc::InvalidArgument,
                "I is the imaginary unit alias and cannot be shadowed as a symbol");
    EXPECT_TRUE(!reserved_I &&
                    std::string(lamina::lsr::error_name(reserved_I.error())) ==
                        "ImaginaryUnitReserved",
                "LSR diagnostic name preserves reserved imaginary unit alias errors");

    auto reserved_pi = lamina::lsr::sym("pi");
    EXPECT_TRUE(!reserved_pi &&
                    reserved_pi.error().code == lamina::CasErrc::InvalidArgument,
                "pi is a constant and cannot be shadowed as a symbol");
    EXPECT_TRUE(!reserved_pi &&
                    std::string(lamina::lsr::error_name(reserved_pi.error())) ==
                        "InvalidArgument",
                "non-imaginary reserved constants keep the generic argument diagnostic");

    auto reserved_unicode_pi = lamina::lsr::sym("\xCF\x80");
    EXPECT_TRUE(!reserved_unicode_pi &&
                    reserved_unicode_pi.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "unicode pi is a constant alias and cannot be shadowed");

    auto reserved_e = lamina::lsr::sym("e");
    EXPECT_TRUE(!reserved_e &&
                    reserved_e.error().code == lamina::CasErrc::InvalidArgument,
                "e is a constant and cannot be shadowed as a symbol");

    auto reserved_phi = lamina::lsr::sym("phi");
    EXPECT_TRUE(!reserved_phi &&
                    reserved_phi.error().code == lamina::CasErrc::InvalidArgument,
                "phi is a constant and cannot be shadowed as a symbol");

    TEST_CASE("LSR std.math constants have explicit Expr constructors");

    auto pi_constant = lamina::lsr::pi();
    auto e_constant = lamina::lsr::e();
    auto phi_constant = lamina::lsr::phi();
    EXPECT_TRUE(pi_constant.has_value(), "std.math.pi Expr can be constructed");
    EXPECT_TRUE(e_constant.has_value(), "std.math.e Expr can be constructed");
    EXPECT_TRUE(phi_constant.has_value(), "std.math.phi Expr can be constructed");
    auto pi_value = pi_constant ? lamina::lsr::evalf(*pi_constant.value())
                                : lamina::Result<lamina::ApproxReal>::failure(
                                      lamina::CasErrc::InternalInvariant,
                                      "pi construction failed", "test");
    auto e_value = e_constant ? lamina::lsr::evalf(*e_constant.value())
                              : lamina::Result<lamina::ApproxReal>::failure(
                                    lamina::CasErrc::InternalInvariant,
                                    "e construction failed", "test");
    auto phi_value = phi_constant ? lamina::lsr::evalf(*phi_constant.value())
                                  : lamina::Result<lamina::ApproxReal>::failure(
                                        lamina::CasErrc::InternalInvariant,
                                        "phi construction failed", "test");
    EXPECT_TRUE(pi_value && pi_value.value().is_finite(),
                "std.math.pi explicitly evaluates through evalf");
    EXPECT_TRUE(e_value && e_value.value().is_finite(),
                "std.math.e explicitly evaluates through evalf");
    EXPECT_TRUE(phi_value && phi_value.value().is_finite(),
                "std.math.phi explicitly evaluates through evalf");
    EXPECT_NEAR(pi_value.value().value, LMMC_CONST_PI, 1e-15,
                "std.math.pi numeric value matches LMMC");
    EXPECT_NEAR(e_value.value().value, std::exp(1.0), 1e-15,
                "std.math.e numeric value matches exp(1)");
    EXPECT_NEAR(phi_value.value().value, (1.0 + std::sqrt(5.0)) / 2.0,
                1e-15, "std.math.phi numeric value matches golden ratio");
    auto unicode_pi_value = lamina::lsr::evalf(
        *SymbolicExpr::variable("\xCF\x80"));
    EXPECT_TRUE(unicode_pi_value && unicode_pi_value.value().is_finite(),
                "unicode pi compatibility alias explicitly evaluates through evalf");
    EXPECT_NEAR(unicode_pi_value.value().value, LMMC_CONST_PI, 1e-15,
                "unicode pi compatibility alias matches std.math.pi");

    TEST_CASE("LSR approximate real construction is explicit and finite");

    auto approximate_half = lamina::lsr::approx_real(0.5);
    EXPECT_TRUE(approximate_half.has_value(),
                "approx_real constructs an explicit approximate Expr");
    auto approximate_half_value =
        approximate_half ? lamina::lsr::evalf(*approximate_half.value())
                         : lamina::Result<lamina::ApproxReal>::failure(
                               lamina::CasErrc::InternalInvariant,
                               "approx_real construction failed", "test");
    EXPECT_TRUE(approximate_half_value &&
                    approximate_half_value.value().is_finite(),
                "explicit approximate Expr can be evaluated with evalf");
    EXPECT_NEAR(approximate_half_value.value().value, 0.5, 0.0,
                "approx_real preserves the requested finite value");

    auto nan_approx = lamina::lsr::approx_real(NAN);
    auto inf_approx = lamina::lsr::approx_real(INFINITY);
    EXPECT_TRUE(!nan_approx &&
                    nan_approx.error().code == lamina::CasErrc::InvalidArgument,
                "approx_real rejects NaN");
    EXPECT_TRUE(!inf_approx &&
                    inf_approx.error().code == lamina::CasErrc::InvalidArgument,
                "approx_real rejects infinity");

    TEST_CASE("LSR Expr arithmetic wrappers return Result values");

    auto expr_x = lamina::lsr::sym("expr_x");
    auto expr_two = lamina::lsr::integer(2);
    auto expr_three = lamina::lsr::integer(3);
    auto expr_five = lamina::lsr::integer(5);
    auto expr_sum = lamina::lsr::add(expr_two.value(), expr_three.value());
    auto expr_product =
        lamina::lsr::mul(expr_sum.value(), SymbolicExpr::number(4));
    auto expr_quotient =
        lamina::lsr::div(expr_product.value(), SymbolicExpr::number(2));
    auto expr_difference =
        lamina::lsr::sub(expr_quotient.value(), SymbolicExpr::number(5));
    auto expr_negated = lamina::lsr::neg(expr_difference.value());
    EXPECT_TRUE(expr_sum && expr_product && expr_quotient &&
                    expr_difference && expr_negated,
                "Expr arithmetic wrappers construct symbolic expressions");

    auto expr_negated_value =
        expr_negated ? lamina::lsr::evalf(*expr_negated.value())
                     : lamina::Result<lamina::ApproxReal>::failure(
                           lamina::CasErrc::InternalInvariant,
                           "arithmetic construction failed", "test");
    EXPECT_NEAR(expr_negated_value.value().value, -5.0, 0.0,
                "Expr arithmetic wrappers evaluate explicitly");

    auto expr_polynomial =
        lamina::lsr::add(expr_x.value(), expr_five.value());
    auto expr_equation =
        lamina::lsr::eq(expr_polynomial.value(), SymbolicExpr::number(0));
    auto expr_solved =
        lamina::lsr::solve_expr_set(expr_polynomial.value(), "expr_x");
    EXPECT_TRUE(expr_equation && !expr_equation.value()->to_string().empty(),
                "Expr eq wrapper constructs a relational expression");
    EXPECT_TRUE(expr_solved &&
                    expr_solved.value().contains(*SymbolicExpr::number(-5)),
                "Expr arithmetic wrappers feed the LSR solve set facade");

    auto null_add = lamina::lsr::add(nullptr, expr_two.value());
    auto null_div = lamina::lsr::div(expr_two.value(), nullptr);
    auto null_neg = lamina::lsr::neg(nullptr);
    auto null_eq = lamina::lsr::eq(expr_two.value(), nullptr);
    EXPECT_TRUE(!null_add &&
                    null_add.error().code == lamina::CasErrc::InvalidArgument,
                "Expr add rejects null input");
    EXPECT_TRUE(!null_div &&
                    null_div.error().code == lamina::CasErrc::InvalidArgument,
                "Expr div rejects null input");
    EXPECT_TRUE(!null_neg &&
                    null_neg.error().code == lamina::CasErrc::InvalidArgument,
                "Expr neg rejects null input");
    EXPECT_TRUE(!null_eq &&
                    null_eq.error().code == lamina::CasErrc::InvalidArgument,
                "Expr eq rejects null input");

    lamina::ResourceLimits exhausted_expr_limits;
    exhausted_expr_limits.max_steps = 0;
    lamina::ComputationContext exhausted_expr_context(exhausted_expr_limits);
    auto exhausted_expr_add =
        lamina::lsr::add(expr_two.value(), expr_three.value(),
                         exhausted_expr_context);
    EXPECT_TRUE(!exhausted_expr_add &&
                    exhausted_expr_add.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "Expr arithmetic wrappers observe the computation budget");

    TEST_CASE("LSR Expr transform wrappers return Result values");

    auto transform_x = lamina::lsr::sym("transform_x");
    auto transform_zero = lamina::lsr::integer(0);
    auto transform_one = lamina::lsr::integer(1);
    auto transform_two = lamina::lsr::integer(2);
    auto transform_three = lamina::lsr::integer(3);
    auto transform_x_plus_zero =
        lamina::lsr::add(transform_x.value(), transform_zero.value());
    auto transform_simplified =
        lamina::lsr::simplify(transform_x_plus_zero.value());
    auto transform_simplified_value =
        transform_simplified
            ? lamina::lsr::evalf(*transform_simplified.value(),
                                 lamina::NumericBindings{{"transform_x", 7.0}})
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "simplify construction failed", "test");
    EXPECT_TRUE(transform_simplified &&
                    transform_simplified_value &&
                    transform_simplified_value.value().value == 7.0,
                "simplify lowers through the LSR Result facade");

    auto transform_left =
        lamina::lsr::add(transform_x.value(), transform_one.value());
    auto transform_right =
        lamina::lsr::add(transform_x.value(), transform_two.value());
    auto transform_product =
        lamina::lsr::mul(transform_left.value(), transform_right.value());
    auto transform_expanded =
        lamina::lsr::expand(transform_product.value());
    auto transform_expanded_value =
        transform_expanded
            ? lamina::lsr::evalf(*transform_expanded.value(),
                                 lamina::NumericBindings{{"transform_x", 3.0}})
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "expand construction failed", "test");
    EXPECT_TRUE(transform_expanded && transform_expanded_value &&
                    transform_expanded_value.value().value == 20.0,
                "expand lowers through the LSR Result facade");

    auto transform_x_cubed =
        lamina::lsr::pow(transform_x.value(), transform_three.value());
    auto transform_derivative =
        lamina::lsr::differentiate(transform_x_cubed.value(), "transform_x");
    auto transform_derivative_value =
        transform_derivative
            ? lamina::lsr::evalf(*transform_derivative.value(),
                                 lamina::NumericBindings{{"transform_x", 2.0}})
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "differentiate construction failed", "test");
    EXPECT_TRUE(transform_derivative && transform_derivative_value &&
                    transform_derivative_value.value().value == 12.0,
                "differentiate lowers through the LSR Result facade");

    auto null_simplify = lamina::lsr::simplify(nullptr);
    auto null_expand = lamina::lsr::expand(nullptr);
    auto null_differentiate = lamina::lsr::differentiate(nullptr, "x");
    auto empty_variable =
        lamina::lsr::differentiate(transform_x.value(), "");
    EXPECT_TRUE(!null_simplify &&
                    null_simplify.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "simplify rejects null input");
    EXPECT_TRUE(!null_expand &&
                    null_expand.error().code == lamina::CasErrc::InvalidArgument,
                "expand rejects null input");
    EXPECT_TRUE(!null_differentiate &&
                    null_differentiate.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "differentiate rejects null input");
    EXPECT_TRUE(!empty_variable &&
                    empty_variable.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "differentiate rejects empty variable names");

    lamina::ResourceLimits exhausted_transform_limits;
    exhausted_transform_limits.max_steps = 0;
    lamina::ComputationContext exhausted_transform_context(
        exhausted_transform_limits);
    auto exhausted_transform =
        lamina::lsr::simplify(transform_x.value(), exhausted_transform_context);
    EXPECT_TRUE(!exhausted_transform &&
                    exhausted_transform.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "Expr transform wrappers observe the computation budget");

    TEST_CASE("LSR std.math Expr wrappers return Result values");

    auto quarter_pi = SymbolicExpr::divide(pi_constant.value(),
                                           SymbolicExpr::number(4));
    auto math_sin_quarter_pi = lamina::lsr::sin(quarter_pi);
    auto math_cos_zero = lamina::lsr::cos(SymbolicExpr::number(0));
    auto math_tan_zero = lamina::lsr::tan(SymbolicExpr::number(0));
    auto math_sqrt_four = lamina::lsr::sqrt(SymbolicExpr::number(4));
    auto math_pow_two_three =
        lamina::lsr::pow(SymbolicExpr::number(2), SymbolicExpr::number(3));
    auto math_asin_half = lamina::lsr::asin(SymbolicExpr::number(0.5));
    auto math_acos_one = lamina::lsr::acos(SymbolicExpr::number(1));
    auto math_atan_one = lamina::lsr::atan(SymbolicExpr::number(1));
    auto math_exp_zero = lamina::lsr::exp(SymbolicExpr::number(0));
    auto math_log_e = lamina::lsr::log(e_constant.value());
    auto math_log10_hundred = lamina::lsr::log10(SymbolicExpr::number(100));
    auto math_floor = lamina::lsr::floor(SymbolicExpr::number(2.75));
    auto math_ceil = lamina::lsr::ceil(SymbolicExpr::number(2.25));
    auto math_round = lamina::lsr::round(SymbolicExpr::number(-2.5));
    auto math_clamp = lamina::lsr::clamp(SymbolicExpr::number(7),
                                         SymbolicExpr::number(0),
                                         SymbolicExpr::number(5));
    EXPECT_TRUE(math_sin_quarter_pi && math_cos_zero && math_tan_zero &&
                    math_sqrt_four && math_pow_two_three && math_exp_zero &&
                    math_asin_half && math_acos_one && math_atan_one &&
                    math_log_e && math_log10_hundred && math_floor &&
                    math_ceil && math_round && math_clamp,
                "std.math Expr wrappers construct symbolic expressions");

    auto sin_quarter_pi_value =
        math_sin_quarter_pi
            ? lamina::lsr::evalf(*math_sin_quarter_pi.value())
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "sin construction failed", "test");
    auto cos_zero_value =
        math_cos_zero ? lamina::lsr::evalf(*math_cos_zero.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "cos construction failed", "test");
    auto tan_zero_value =
        math_tan_zero ? lamina::lsr::evalf(*math_tan_zero.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "tan construction failed", "test");
    auto sqrt_four_value =
        math_sqrt_four ? lamina::lsr::evalf(*math_sqrt_four.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "sqrt construction failed", "test");
    auto pow_two_three_value =
        math_pow_two_three
            ? lamina::lsr::evalf(*math_pow_two_three.value())
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "pow construction failed", "test");
    auto asin_half_value =
        math_asin_half ? lamina::lsr::evalf(*math_asin_half.value())
                       : lamina::Result<lamina::ApproxReal>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "asin construction failed", "test");
    auto acos_one_value =
        math_acos_one ? lamina::lsr::evalf(*math_acos_one.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "acos construction failed", "test");
    auto atan_one_value =
        math_atan_one ? lamina::lsr::evalf(*math_atan_one.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "atan construction failed", "test");
    auto exp_zero_value =
        math_exp_zero ? lamina::lsr::evalf(*math_exp_zero.value())
                      : lamina::Result<lamina::ApproxReal>::failure(
                            lamina::CasErrc::InternalInvariant,
                            "exp construction failed", "test");
    auto log_e_value =
        math_log_e ? lamina::lsr::evalf(*math_log_e.value())
                   : lamina::Result<lamina::ApproxReal>::failure(
                         lamina::CasErrc::InternalInvariant,
                         "log construction failed", "test");
    auto log10_hundred_value =
        math_log10_hundred
            ? lamina::lsr::evalf(*math_log10_hundred.value())
            : lamina::Result<lamina::ApproxReal>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "log10 construction failed", "test");
    auto floor_value =
        math_floor ? lamina::lsr::evalf(*math_floor.value())
                   : lamina::Result<lamina::ApproxReal>::failure(
                         lamina::CasErrc::InternalInvariant,
                         "floor construction failed", "test");
    auto ceil_value =
        math_ceil ? lamina::lsr::evalf(*math_ceil.value())
                  : lamina::Result<lamina::ApproxReal>::failure(
                        lamina::CasErrc::InternalInvariant,
                        "ceil construction failed", "test");
    auto round_value =
        math_round ? lamina::lsr::evalf(*math_round.value())
                   : lamina::Result<lamina::ApproxReal>::failure(
                         lamina::CasErrc::InternalInvariant,
                         "round construction failed", "test");
    auto clamp_value =
        math_clamp ? lamina::lsr::evalf(*math_clamp.value())
                   : lamina::Result<lamina::ApproxReal>::failure(
                         lamina::CasErrc::InternalInvariant,
                         "clamp construction failed", "test");

    EXPECT_NEAR(sin_quarter_pi_value.value().value, std::sqrt(0.5), 1e-12,
                "sin(pi / 4) evaluates explicitly");
    EXPECT_NEAR(cos_zero_value.value().value, 1.0, 0.0,
                "cos(0) evaluates explicitly");
    EXPECT_NEAR(tan_zero_value.value().value, 0.0, 0.0,
                "tan(0) evaluates explicitly");
    EXPECT_NEAR(sqrt_four_value.value().value, 2.0, 0.0,
                "sqrt(4) evaluates explicitly");
    EXPECT_NEAR(pow_two_three_value.value().value, 8.0, 0.0,
                "pow(2, 3) evaluates explicitly");
    EXPECT_NEAR(asin_half_value.value().value, std::asin(0.5), 1e-12,
                "asin(0.5) evaluates explicitly");
    EXPECT_NEAR(acos_one_value.value().value, 0.0, 0.0,
                "acos(1) evaluates explicitly");
    EXPECT_NEAR(atan_one_value.value().value, std::atan(1.0), 1e-12,
                "atan(1) evaluates explicitly");
    EXPECT_NEAR(exp_zero_value.value().value, 1.0, 0.0,
                "exp(0) evaluates explicitly");
    EXPECT_NEAR(log_e_value.value().value, 1.0, 1e-12,
                "log(e) evaluates explicitly");
    EXPECT_NEAR(log10_hundred_value.value().value, 2.0, 1e-12,
                "log10(100) evaluates explicitly");
    EXPECT_NEAR(floor_value.value().value, 2.0, 0.0,
                "floor(2.75) evaluates explicitly");
    EXPECT_NEAR(ceil_value.value().value, 3.0, 0.0,
                "ceil(2.25) evaluates explicitly");
    EXPECT_NEAR(round_value.value().value, -3.0, 0.0,
                "round(-2.5) evaluates explicitly");
    EXPECT_NEAR(clamp_value.value().value, 5.0, 0.0,
                "clamp(7, 0, 5) evaluates explicitly");

    auto null_sin = lamina::lsr::sin(nullptr);
    auto null_pow = lamina::lsr::pow(SymbolicExpr::number(2), nullptr);
    auto null_clamp = lamina::lsr::clamp(SymbolicExpr::number(1), nullptr,
                                         SymbolicExpr::number(2));
    EXPECT_TRUE(!null_sin &&
                    null_sin.error().code == lamina::CasErrc::InvalidArgument,
                "std.math Expr unary wrappers reject null input");
    EXPECT_TRUE(!null_pow &&
                    null_pow.error().code == lamina::CasErrc::InvalidArgument,
                "std.math Expr binary wrappers reject null input");
    EXPECT_TRUE(!null_clamp &&
                    null_clamp.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "std.math Expr clamp rejects null input");

    lamina::ResourceLimits exhausted_math_limits;
    exhausted_math_limits.max_steps = 0;
    lamina::ComputationContext exhausted_math_context(exhausted_math_limits);
    auto exhausted_sin =
        lamina::lsr::sin(SymbolicExpr::number(1), exhausted_math_context);
    EXPECT_TRUE(!exhausted_sin &&
                    exhausted_sin.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "std.math Expr wrappers observe the computation budget");

    TEST_CASE("LSR imaginary unit is canonical complex zero plus one i");

    auto i = lamina::lsr::imaginary_unit();
    EXPECT_TRUE(i.has_value(), "imaginary unit can be constructed");
    auto lower_i = lamina::lsr::i();
    auto upper_i = lamina::lsr::I();
    EXPECT_TRUE(lower_i.has_value(), "std.math.i Expr alias can be constructed");
    EXPECT_TRUE(upper_i.has_value(), "std.math.I Expr alias can be constructed");
    EXPECT_TRUE(i && lower_i &&
                    lamina::lsr::structurally_equal(*i.value(),
                                                    *lower_i.value()),
                "std.math.i aliases the canonical imaginary unit");
    EXPECT_TRUE(i && upper_i &&
                    lamina::lsr::structurally_equal(*i.value(),
                                                    *upper_i.value()),
                "std.math.I aliases the canonical imaginary unit");

    auto explicit_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(1));
    EXPECT_TRUE(explicit_i.has_value(), "complex(0, 1) can be constructed");
    EXPECT_TRUE(i && explicit_i &&
                    lamina::lsr::structurally_equal(*i.value(),
                                                    *explicit_i.value()),
                "imaginary unit is structurally complex(0, 1)");
    auto null_real_complex =
        lamina::lsr::complex(nullptr, SymbolicExpr::number(1));
    auto null_imag_complex =
        lamina::lsr::complex(SymbolicExpr::number(0), nullptr);
    EXPECT_TRUE(!null_real_complex &&
                    null_real_complex.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "complex(nullptr, 1) rejects an incompatible real part");
    EXPECT_TRUE(!null_imag_complex &&
                    null_imag_complex.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "complex(0, nullptr) rejects an incompatible imaginary part");
    EXPECT_TRUE(!null_real_complex &&
                    std::string(lamina::lsr::error_name(
                        null_real_complex.error())) == "ComplexTypeMismatch",
                "complex(nullptr, 1) exposes ComplexTypeMismatch");
    EXPECT_TRUE(!null_imag_complex &&
                    std::string(lamina::lsr::error_name(
                        null_imag_complex.error())) == "ComplexTypeMismatch",
                "complex(0, nullptr) exposes ComplexTypeMismatch");

    auto variable_i = SymbolicExpr::variable("i");
    EXPECT_TRUE(i && !lamina::lsr::structurally_equal(*i.value(), *variable_i),
                "legacy variable(\"i\") is not structurally the imaginary unit");

    if (i) {
        auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
        lamina::ComputationContext complex_equivalence_context;
        auto complex_equivalent = lamina::lsr::equivalent_core(
            *i_squared, *SymbolicExpr::number(-1), complex_equivalence_context);
        EXPECT_TRUE(complex_equivalent && complex_equivalent.value(),
                    "i * i is equivalent to -1 in the LSR core profile");

        auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
        lamina::ComputationContext complex_power_equivalence_context;
        auto complex_power_equivalent = lamina::lsr::equivalent_core(
            *i_power_two, *SymbolicExpr::number(-1),
            complex_power_equivalence_context);
        EXPECT_TRUE(complex_power_equivalent &&
                        complex_power_equivalent.value(),
                    "i^2 is equivalent to -1 in the LSR core profile");

        auto legacy_i_squared =
            SymbolicExpr::multiply(variable_i, variable_i);
        lamina::ComputationContext legacy_i_equivalence_context;
        auto legacy_i_equivalent = lamina::lsr::equivalent_core(
            *legacy_i_squared, *SymbolicExpr::number(-1),
            legacy_i_equivalence_context);
        EXPECT_TRUE(legacy_i_equivalent && legacy_i_equivalent.value(),
                    "legacy variable(\"i\") follows the LSR i*i == -1 rule in equivalence");

        auto legacy_i_plus_one = SymbolicExpr::add(variable_i,
                                                   SymbolicExpr::number(1));
        auto one_plus_canonical_i =
            SymbolicExpr::add(SymbolicExpr::number(1), i.value());
        lamina::ComputationContext legacy_i_add_context;
        auto legacy_i_add_equivalent = lamina::lsr::equivalent_core(
            *legacy_i_plus_one, *one_plus_canonical_i,
            legacy_i_add_context);
        EXPECT_TRUE(legacy_i_add_equivalent &&
                        legacy_i_add_equivalent.value(),
                    "legacy variable(\"i\") normalizes inside LSR additive equivalence");
    }

    TEST_CASE("LSR evalf is explicit and propagates missing bindings");

    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    auto unbound = lamina::lsr::evalf(*linear);
    EXPECT_TRUE(!unbound &&
                    unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "evalf without a required binding reports UnboundSymbol");
    EXPECT_TRUE(!unbound &&
                    std::string(lamina::lsr::error_name(unbound.error())) ==
                        "UnboundSymbol",
                "evalf keeps the generic unbound symbol diagnostic");

    lamina::NumericBindings bindings{{"x", 3.0}};
    auto evaluated = lamina::lsr::evalf(*linear, bindings);
    EXPECT_TRUE(evaluated && evaluated.value().is_finite(),
                "evalf with a binding succeeds");
    EXPECT_NEAR(evaluated.value().value, 5.0, 0.0,
                "evalf computes the numeric value");

    auto nonfinite_binding = lamina::lsr::evalf(
        *x.value(), lamina::NumericBindings{{"x", INFINITY}});
    EXPECT_TRUE(!nonfinite_binding &&
                    nonfinite_binding.error().code ==
                        lamina::CasErrc::NumericFailure,
                "evalf rejects non-finite numeric bindings");
    EXPECT_TRUE(!nonfinite_binding &&
                    std::string(lamina::lsr::error_name(
                        nonfinite_binding.error())) == "NumericFailure",
                "evalf reports the LSR numeric failure diagnostic for non-finite bindings");

    auto nonfinite_expression = lamina::lsr::evalf(*SymbolicExpr::infinity());
    EXPECT_TRUE(!nonfinite_expression &&
                    nonfinite_expression.error().code ==
                        lamina::CasErrc::NumericFailure,
                "evalf rejects expressions that evaluate to infinity");

    lamina::ResourceLimits exhausted_evalf_limits;
    exhausted_evalf_limits.max_steps = 0;
    lamina::ComputationContext exhausted_evalf_context(exhausted_evalf_limits);
    auto exhausted_evalf =
        lamina::lsr::evalf(*linear, bindings, exhausted_evalf_context);
    EXPECT_TRUE(!exhausted_evalf &&
                    exhausted_evalf.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "evalf reports ResourceLimit when the computation budget is exhausted");
    EXPECT_TRUE(!exhausted_evalf &&
                    std::string(lamina::lsr::error_name(
                        exhausted_evalf.error())) == "ResourceLimit",
                "evalf exposes the LSR resource-limit diagnostic");

    TEST_CASE("LSR substitution explicitly rewrites Expr before evaluation");

    auto substituted = lamina::lsr::substitute(
        linear, "x", SymbolicExpr::number(3));
    EXPECT_TRUE(substituted.has_value(), "substitute(x + 2, x => 3) succeeds");
    auto substituted_value =
        substituted ? lamina::lsr::evalf(*substituted.value())
                    : lamina::Result<lamina::ApproxReal>::failure(
                          lamina::CasErrc::InternalInvariant,
                          "substitution failed", "test");
    EXPECT_TRUE(substituted_value && substituted_value.value().is_finite(),
                "substituted Expr can be explicitly evaluated");
    EXPECT_NEAR(substituted_value.value().value, 5.0, 0.0,
                "substitute(x + 2, x => 3) evaluates to 5");

    auto unchanged_substitution = lamina::lsr::substitute(
        linear, "y", SymbolicExpr::number(9));
    EXPECT_TRUE(unchanged_substitution &&
                    lamina::lsr::structurally_equal(
                        *unchanged_substitution.value(), *linear),
                "substituting an absent symbol leaves the Expr unchanged");

    auto null_substitution_expr = lamina::lsr::substitute(
        nullptr, "x", SymbolicExpr::number(1));
    EXPECT_TRUE(!null_substitution_expr &&
                    null_substitution_expr.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "substitute rejects a null Expr");

    auto empty_substitution_var = lamina::lsr::substitute(
        linear, "", SymbolicExpr::number(1));
    EXPECT_TRUE(!empty_substitution_var &&
                    empty_substitution_var.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "substitute rejects an empty variable name");

    auto null_substitution_value = lamina::lsr::substitute(linear, "x", nullptr);
    EXPECT_TRUE(!null_substitution_value &&
                    null_substitution_value.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "substitute rejects a null replacement Expr");

    lamina::ResourceLimits exhausted_substitution_limits;
    exhausted_substitution_limits.max_steps = 0;
    lamina::ComputationContext exhausted_substitution_context(
        exhausted_substitution_limits);
    auto exhausted_substitution = lamina::lsr::substitute(
        linear, "x", SymbolicExpr::number(1), exhausted_substitution_context);
    EXPECT_TRUE(!exhausted_substitution &&
                    exhausted_substitution.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "substitute observes the computation budget");

    TEST_CASE("LSR expression matching exposes deterministic bindings");

    auto wildcard_a = SymbolicExpr::variable("A");
    auto match_pattern =
        SymbolicExpr::add(wildcard_a, SymbolicExpr::number(1));
    auto match_target =
        SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::variable("x"));
    auto matched =
        lamina::lsr::expr_match(match_pattern, match_target, {"A"});
    EXPECT_TRUE(matched && matched.value().matched,
                "expr_match succeeds for a commutative additive pattern");
    EXPECT_TRUE(matched && matched.value().bindings.size() == 1,
                "expr_match returns one wildcard binding");
    EXPECT_TRUE(matched && matched.value().bindings[0].name == "A",
                "expr_match sorts bindings by wildcard name");
    EXPECT_TRUE(matched && matched.value().bindings[0].value &&
                    lamina::lsr::structurally_equal(
                        *matched.value().bindings[0].value,
                        *SymbolicExpr::variable("x")),
                "expr_match binds A to x");

    auto multiply_pattern = SymbolicExpr::multiply(SymbolicExpr::variable("B"),
                                                   SymbolicExpr::number(2));
    auto multiply_target = SymbolicExpr::multiply(SymbolicExpr::number(2),
                                                  SymbolicExpr::variable("y"));
    auto multiply_matched =
        lamina::lsr::expr_match(multiply_pattern, multiply_target, {"B"});
    EXPECT_TRUE(multiply_matched && multiply_matched.value().matched &&
                    multiply_matched.value().bindings.size() == 1 &&
                    lamina::lsr::structurally_equal(
                        *multiply_matched.value().bindings[0].value,
                        *SymbolicExpr::variable("y")),
                "expr_match succeeds for a commutative multiplicative pattern");

    auto power_pattern = SymbolicExpr::power(SymbolicExpr::variable("U"),
                                             SymbolicExpr::variable("N"));
    auto power_target = SymbolicExpr::power(
        SymbolicExpr::sin(SymbolicExpr::variable("theta")),
        SymbolicExpr::number(2));
    auto power_matched =
        lamina::lsr::expr_match(power_pattern, power_target, {"N", "U"});
    EXPECT_TRUE(power_matched && power_matched.value().matched &&
                    power_matched.value().bindings.size() == 2 &&
                    power_matched.value().bindings[0].name == "N" &&
                    power_matched.value().bindings[1].name == "U" &&
                    lamina::lsr::structurally_equal(
                        *power_matched.value().bindings[0].value,
                        *SymbolicExpr::number(2)) &&
                    lamina::lsr::structurally_equal(
                        *power_matched.value().bindings[1].value,
                        *SymbolicExpr::sin(SymbolicExpr::variable("theta"))),
                "expr_match exposes LSR power pattern bindings deterministically");

    auto function_matched = lamina::lsr::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("U")),
        SymbolicExpr::sin(SymbolicExpr::add(SymbolicExpr::variable("theta"),
                                            SymbolicExpr::number(1))),
        {"U"});
    EXPECT_TRUE(function_matched && function_matched.value().matched &&
                    function_matched.value().bindings.size() == 1,
                "expr_match supports function-node patterns");

    auto repeated_wildcard_match = lamina::lsr::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(SymbolicExpr::variable("z"),
                          SymbolicExpr::variable("z")),
        {"A"});
    auto inconsistent_wildcard_match = lamina::lsr::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(SymbolicExpr::variable("z"),
                          SymbolicExpr::variable("w")),
        {"A"});
    EXPECT_TRUE(repeated_wildcard_match &&
                    repeated_wildcard_match.value().matched,
                "expr_match accepts repeated wildcards with identical bindings");
    EXPECT_TRUE(inconsistent_wildcard_match &&
                    !inconsistent_wildcard_match.value().matched,
                "expr_match rejects repeated wildcards with inconsistent bindings");

    auto unmatched = lamina::lsr::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("A")),
        SymbolicExpr::cos(SymbolicExpr::variable("x")), {"A"});
    EXPECT_TRUE(unmatched && !unmatched.value().matched &&
                    unmatched.value().bindings.empty(),
                "expr_match reports structural non-matches without error");

    auto null_pattern_match =
        lamina::lsr::expr_match(nullptr, match_target, {"A"});
    auto null_target_match =
        lamina::lsr::expr_match(match_pattern, nullptr, {"A"});
    auto empty_wildcard_match =
        lamina::lsr::expr_match(match_pattern, match_target, {""});
    auto duplicate_wildcard_match =
        lamina::lsr::expr_match(match_pattern, match_target, {"A", "A"});
    EXPECT_TRUE(!null_pattern_match &&
                    null_pattern_match.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "expr_match rejects a null pattern");
    EXPECT_TRUE(!null_target_match &&
                    null_target_match.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "expr_match rejects a null target");
    EXPECT_TRUE(!empty_wildcard_match &&
                    empty_wildcard_match.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "expr_match rejects an empty wildcard name");
    EXPECT_TRUE(!duplicate_wildcard_match &&
                    duplicate_wildcard_match.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "expr_match rejects duplicate wildcard names");

    lamina::ResourceLimits exhausted_match_limits;
    exhausted_match_limits.max_steps = 0;
    lamina::ComputationContext exhausted_match_context(exhausted_match_limits);
    auto exhausted_match = lamina::lsr::expr_match(
        match_pattern, match_target, {"A"}, exhausted_match_context);
    EXPECT_TRUE(!exhausted_match &&
                    exhausted_match.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "expr_match observes the computation budget");

    TEST_CASE("LSR eval_complex explicitly lowers Expr to complex");

    auto real_as_complex = lamina::lsr::eval_complex(*SymbolicExpr::number(5));
    EXPECT_TRUE(real_as_complex && real_as_complex.value().is_finite(),
                "eval_complex accepts real expressions explicitly");
    EXPECT_NEAR(real_as_complex.value().real.value, 5.0, 0.0,
                "eval_complex preserves real component");
    EXPECT_NEAR(real_as_complex.value().imag.value, 0.0, 0.0,
                "eval_complex promotes real expression with zero imaginary component");

    auto legacy_i_complex =
        lamina::lsr::eval_complex(*SymbolicExpr::variable("i"));
    auto legacy_upper_i_complex =
        lamina::lsr::eval_complex(*SymbolicExpr::variable("I"));
    EXPECT_TRUE(legacy_i_complex && legacy_i_complex.value().is_finite() &&
                    legacy_upper_i_complex &&
                    legacy_upper_i_complex.value().is_finite(),
                "eval_complex treats legacy i/I Expr variables as the imaginary unit");
    EXPECT_NEAR(legacy_i_complex.value().real.value, 0.0, 0.0,
                "legacy i has zero real component");
    EXPECT_NEAR(legacy_i_complex.value().imag.value, 1.0, 0.0,
                "legacy i has unit imaginary component");
    EXPECT_NEAR(legacy_upper_i_complex.value().imag.value, 1.0, 0.0,
                "legacy I aliases the imaginary unit during complex evaluation");

    auto four_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = lamina::lsr::eval_complex(*three_plus_four_i);
    EXPECT_TRUE(lowered_complex && lowered_complex.value().is_finite(),
                "eval_complex lowers 3 + 4i");
    EXPECT_NEAR(lowered_complex.value().real.value, 3.0, 0.0,
                "eval_complex computes real part of 3 + 4i");
    EXPECT_NEAR(lowered_complex.value().imag.value, 4.0, 0.0,
                "eval_complex computes imaginary part of 3 + 4i");

    if (i) {
        auto ordinary_multiply_complex = SymbolicExpr::add(
            SymbolicExpr::number(3),
            SymbolicExpr::multiply(SymbolicExpr::number(4), i.value()));
        auto lowered_ordinary_multiply =
            lamina::lsr::eval_complex(*ordinary_multiply_complex);
        EXPECT_TRUE(lowered_ordinary_multiply &&
                        lowered_ordinary_multiply.value().is_finite(),
                    "eval_complex lowers the LSR 3 + 4 * i ordinary multiplication form");
        EXPECT_NEAR(lowered_ordinary_multiply.value().real.value, 3.0, 0.0,
                    "ordinary multiplication complex form preserves real part");
        EXPECT_NEAR(lowered_ordinary_multiply.value().imag.value, 4.0, 0.0,
                    "ordinary multiplication complex form preserves imaginary part");

        auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
        auto lowered_i_squared = lamina::lsr::eval_complex(*i_power_two);
        EXPECT_TRUE(lowered_i_squared && lowered_i_squared.value().is_finite(),
                    "eval_complex supports the LSR i^2 rule");
        EXPECT_NEAR(lowered_i_squared.value().real.value, -1.0, 0.0,
                    "eval_complex computes i^2 real part");
        EXPECT_NEAR(lowered_i_squared.value().imag.value, 0.0, 0.0,
                    "eval_complex computes i^2 imaginary part");
    }

    auto complex_unbound = lamina::lsr::eval_complex(*linear);
    EXPECT_TRUE(!complex_unbound &&
                    complex_unbound.error().code == lamina::CasErrc::UnboundSymbol,
                "eval_complex rejects unbound symbols during Expr to complex lowering");
    EXPECT_TRUE(!complex_unbound &&
                    std::string(lamina::lsr::error_name(complex_unbound.error())) ==
                        "ComplexEvalUnboundSymbol",
                "eval_complex exposes the LSR complex unbound-symbol diagnostic");

    if (i) {
        auto fractional_power = SymbolicExpr::power(i.value(), SymbolicExpr::number(0.5));
        auto unsupported_power = lamina::lsr::eval_complex(*fractional_power);
        EXPECT_TRUE(!unsupported_power &&
                        unsupported_power.error().code ==
                            lamina::CasErrc::UnsupportedExpression,
                    "eval_complex does not silently approximate unsupported complex powers");

        auto nonfinite_power = SymbolicExpr::power(
            i.value(), SymbolicExpr::infinity());
        auto nonfinite_power_result =
            lamina::lsr::eval_complex(*nonfinite_power);
        EXPECT_TRUE(!nonfinite_power_result &&
                        nonfinite_power_result.error().code ==
                            lamina::CasErrc::NumericFailure,
                    "eval_complex rejects non-finite complex power exponents");
        EXPECT_TRUE(!nonfinite_power_result &&
                        std::string(lamina::lsr::error_name(
                            nonfinite_power_result.error())) ==
                            "NumericFailure",
                    "eval_complex reports NumericFailure for non-finite exponents");
    }

    auto zero_inverse = lamina::lsr::eval_complex(
        *SymbolicExpr::power(SymbolicExpr::number(0), SymbolicExpr::number(-1)));
    EXPECT_TRUE(!zero_inverse &&
                    zero_inverse.error().code == lamina::CasErrc::DomainError,
                "eval_complex reports DomainError for complex reciprocal of zero");

    lamina::ResourceLimits exhausted_complex_limits;
    exhausted_complex_limits.max_steps = 0;
    lamina::ComputationContext exhausted_complex_context(exhausted_complex_limits);
    auto exhausted_complex = lamina::lsr::eval_complex(
        *three_plus_four_i, {}, exhausted_complex_context);
    EXPECT_TRUE(!exhausted_complex &&
                    exhausted_complex.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "eval_complex reports ResourceLimit when the computation budget is exhausted");
    EXPECT_TRUE(!exhausted_complex &&
                    std::string(lamina::lsr::error_name(
                        exhausted_complex.error())) == "ResourceLimit",
                "eval_complex exposes the LSR resource-limit diagnostic");

    TEST_CASE("LSR complex part functions expose std.math boundaries");

    auto real_part = lamina::lsr::real(three_plus_four_i);
    EXPECT_TRUE(real_part && lamina::lsr::structurally_equal(
                                 *real_part.value(), *SymbolicExpr::number(3)),
                "real(3 + 4i) returns 3");

    auto imag_part = lamina::lsr::imag(three_plus_four_i);
    EXPECT_TRUE(imag_part && lamina::lsr::structurally_equal(
                                 *imag_part.value(), *SymbolicExpr::number(4)),
                "imag(3 + 4i) returns 4");

    auto conjugated = lamina::lsr::conj(three_plus_four_i);
    EXPECT_TRUE(conjugated.has_value(), "conj(3 + 4i) succeeds");
    auto expected_conj = lamina::lsr::complex(SymbolicExpr::number(3),
                                             SymbolicExpr::number(-4));
    lamina::ComputationContext conj_context;
    auto conj_equiv = lamina::lsr::equivalent_core(
        *conjugated.value(), *expected_conj.value(), conj_context);
    EXPECT_TRUE(conj_equiv && conj_equiv.value(),
                "conj(3 + 4i) returns 3 - 4i");

    auto complex_abs = lamina::lsr::abs(three_plus_four_i);
    EXPECT_TRUE(complex_abs.has_value(), "abs(3 + 4i) succeeds");
    auto abs_value = lamina::lsr::evalf(*complex_abs.value());
    EXPECT_TRUE(abs_value && abs_value.value().is_finite(),
                "abs(3 + 4i) can be explicitly numerically evaluated");
    EXPECT_NEAR(abs_value.value().value, 5.0, 1e-12,
                "abs(3 + 4i) evaluates to 5");

    auto null_real = lamina::lsr::real(nullptr);
    auto null_imag = lamina::lsr::imag(nullptr);
    auto null_conj = lamina::lsr::conj(nullptr);
    auto null_abs = lamina::lsr::abs(nullptr);
    EXPECT_TRUE(!null_real &&
                    null_real.error().code == lamina::CasErrc::InvalidArgument,
                "real(nullptr) rejects null LSR Expr");
    EXPECT_TRUE(!null_imag &&
                    null_imag.error().code == lamina::CasErrc::InvalidArgument,
                "imag(nullptr) rejects null LSR Expr");
    EXPECT_TRUE(!null_conj &&
                    null_conj.error().code == lamina::CasErrc::InvalidArgument,
                "conj(nullptr) rejects null LSR Expr");
    EXPECT_TRUE(!null_abs &&
                    null_abs.error().code == lamina::CasErrc::InvalidArgument,
                "abs(nullptr) rejects null LSR Expr");

    auto real_number = SymbolicExpr::number(-5);
    auto real_number_part = lamina::lsr::real(real_number);
    EXPECT_TRUE(real_number_part &&
                    lamina::lsr::structurally_equal(*real_number_part.value(),
                                                    *real_number),
                "real(-5) preserves a real value under R subset C");

    auto real_number_imag = lamina::lsr::imag(real_number);
    EXPECT_TRUE(real_number_imag &&
                    lamina::lsr::structurally_equal(*real_number_imag.value(),
                                                    *SymbolicExpr::number(0)),
                "imag(-5) returns zero for a real value under R subset C");

    auto real_number_conj = lamina::lsr::conj(real_number);
    EXPECT_TRUE(real_number_conj &&
                    lamina::lsr::structurally_equal(*real_number_conj.value(),
                                                    *real_number),
                "conj(-5) preserves a real value under R subset C");

    auto real_number_abs = lamina::lsr::abs(real_number);
    auto real_number_abs_value =
        real_number_abs ? lamina::lsr::evalf(*real_number_abs.value())
                        : lamina::Result<lamina::ApproxReal>::failure(
                              lamina::CasErrc::InternalInvariant,
                              "abs(-5) construction failed", "test");
    EXPECT_TRUE(real_number_abs_value &&
                    real_number_abs_value.value().is_finite(),
                "abs(-5) can be explicitly evaluated under R subset C");
    EXPECT_NEAR(real_number_abs_value.value().value, 5.0, 1e-12,
                "abs(-5) evaluates to 5 under R subset C");

    auto unsupported_real = lamina::lsr::real(
        SymbolicExpr::sin(three_plus_four_i));
    EXPECT_TRUE(!unsupported_real &&
                    unsupported_real.error().code == lamina::CasErrc::Inconclusive,
                "real(sin(3 + 4i)) reports unsupported complex function split");

    TEST_CASE("LSR solve_set returns mathematical sets");

    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    lamina::ComputationContext context;
    auto solved = lamina::lsr::solve_set(equation, "x", context);
    EXPECT_TRUE(solved &&
                    solved.value().kind() == lamina::SolutionSet::Kind::Finite,
                "solve_set returns a finite solution set for x^2 - 1");
    EXPECT_TRUE(solved && solved.value().finite_solutions().size() == 2,
                "solve_set preserves both roots");

    auto empty = lamina::lsr::solve_set(SymbolicExpr::number(1), "x");
    EXPECT_TRUE(empty &&
                    empty.value().kind() == lamina::SolutionSet::Kind::Empty,
                "solve_set represents mathematical no-solution as Empty");

    TEST_CASE("LSR ExprSet implements set<Expr> finite collection semantics");

    auto one_a = SymbolicExpr::number(1);
    auto one_b = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto base_set = lamina::lsr::expr_set({one_a, one_b, two});
    EXPECT_TRUE(base_set && base_set.value().size() == 2,
                "set<Expr> removes structurally equal duplicates");
    EXPECT_TRUE(base_set && base_set.value().contains(*SymbolicExpr::number(1)),
                "set<Expr> membership uses structural equality");
    auto facade_contains_one = base_set
        ? lamina::lsr::expr_set_contains(base_set.value(),
                                         SymbolicExpr::number(1))
        : lamina::Result<bool>::failure(lamina::CasErrc::InternalInvariant,
                                        "set construction failed", "test");
    EXPECT_TRUE(facade_contains_one && facade_contains_one.value(),
                "set<Expr> in operator facade reports membership");
    auto facade_not_contains_three = base_set
        ? lamina::lsr::expr_set_not_contains(base_set.value(),
                                             SymbolicExpr::number(3))
        : lamina::Result<bool>::failure(lamina::CasErrc::InternalInvariant,
                                        "set construction failed", "test");
    EXPECT_TRUE(facade_not_contains_three && facade_not_contains_three.value(),
                "set<Expr> not in operator facade reports non-membership");
    auto empty_expr_set = lamina::lsr::expr_set({});
    EXPECT_TRUE(empty_expr_set && empty_expr_set.value().empty(),
                "set<Expr> can represent the empty finite set");

    auto rhs_set = lamina::lsr::expr_set({SymbolicExpr::number(2),
                                          SymbolicExpr::number(3)});
    EXPECT_TRUE(rhs_set.has_value(), "second set<Expr> can be created");
    if (base_set && rhs_set && empty_expr_set) {
        auto union_set = base_set.value().set_union(rhs_set.value());
        auto facade_union =
            lamina::lsr::expr_set_union(base_set.value(), rhs_set.value());
        EXPECT_TRUE(facade_union && facade_union.value().size() == 3,
                    "set<Expr> union facade returns deduplicated elements");
        EXPECT_TRUE(union_set.size() == 3,
                    "set<Expr> union returns deduplicated elements");
        auto intersection_set = base_set.value().intersection(rhs_set.value());
        auto facade_intersection = lamina::lsr::expr_set_intersection(
            base_set.value(), rhs_set.value());
        EXPECT_TRUE(facade_intersection &&
                        facade_intersection.value().size() == 1 &&
                        facade_intersection.value().contains(
                            *SymbolicExpr::number(2)),
                    "set<Expr> intersection facade keeps common elements");
        EXPECT_TRUE(intersection_set.size() == 1 &&
                        intersection_set.contains(*SymbolicExpr::number(2)),
                    "set<Expr> intersection keeps common elements");
        auto difference_set = base_set.value().difference(rhs_set.value());
        auto facade_difference = lamina::lsr::expr_set_difference(
            base_set.value(), rhs_set.value());
        EXPECT_TRUE(facade_difference &&
                        facade_difference.value().size() == 1 &&
                        facade_difference.value().contains(
                            *SymbolicExpr::number(1)),
                    "set<Expr> difference facade removes right-hand elements");
        EXPECT_TRUE(difference_set.size() == 1 &&
                        difference_set.contains(*SymbolicExpr::number(1)),
                    "set<Expr> difference removes right-hand elements");
        auto symmetric = base_set.value().symmetric_difference(rhs_set.value());
        auto facade_symmetric = lamina::lsr::expr_set_symmetric_difference(
            base_set.value(), rhs_set.value());
        EXPECT_TRUE(facade_symmetric &&
                        facade_symmetric.value().size() == 2 &&
                        facade_symmetric.value().contains(
                            *SymbolicExpr::number(1)) &&
                        facade_symmetric.value().contains(
                            *SymbolicExpr::number(3)),
                    "set<Expr> xor facade follows symmetric difference semantics");
        EXPECT_TRUE(symmetric.size() == 2 &&
                        symmetric.contains(*SymbolicExpr::number(1)) &&
                        symmetric.contains(*SymbolicExpr::number(3)),
                    "set<Expr> symmetric difference follows xor semantics");
        auto facade_subset =
            lamina::lsr::expr_set_subset(intersection_set, union_set);
        EXPECT_TRUE(facade_subset && facade_subset.value(),
                    "set<Expr> subset facade checks membership of every element");
        EXPECT_TRUE(intersection_set.subset_of(union_set),
                    "set<Expr> subset checks membership of every element");
        EXPECT_TRUE(empty_expr_set.value().subset_of(base_set.value()),
                    "empty set<Expr> is a subset of every set<Expr>");
        EXPECT_TRUE(base_set.value().set_union(empty_expr_set.value()).size() ==
                        base_set.value().size(),
                    "set<Expr> union with empty preserves the left set");
        EXPECT_TRUE(base_set.value().intersection(empty_expr_set.value()).empty(),
                    "set<Expr> intersection with empty is empty");
        EXPECT_TRUE(base_set.value().difference(empty_expr_set.value()).size() ==
                        base_set.value().size(),
                    "set<Expr> difference by empty preserves the left set");
    }

    auto null_set = lamina::lsr::expr_set({nullptr});
    EXPECT_TRUE(!null_set &&
                    null_set.error().code == lamina::CasErrc::InvalidArgument,
                "set<Expr> rejects null elements");
    EXPECT_TRUE(!null_set &&
                    std::string(lamina::lsr::error_name(null_set.error())) ==
                        "SetElementTypeMismatch",
                "set<Expr> construction exposes the LSR element type diagnostic");
    if (base_set) {
        auto null_membership =
            lamina::lsr::expr_set_contains(base_set.value(), nullptr);
        EXPECT_TRUE(!null_membership &&
                        null_membership.error().code ==
                            lamina::CasErrc::InvalidArgument,
                    "set<Expr> membership rejects null elements");
        EXPECT_TRUE(!null_membership &&
                        std::string(lamina::lsr::error_name(
                            null_membership.error())) == "SetElementTypeMismatch",
                    "set<Expr> membership exposes the element type diagnostic");
    }

    TEST_CASE("LSR predefined number domains expose Z Q R C set semantics");

    auto domain_z = lamina::lsr::integers();
    auto domain_q = lamina::lsr::rationals();
    auto domain_r = lamina::lsr::reals();
    auto domain_c = lamina::lsr::complexes();
    auto domain_expr = lamina::lsr::expressions();
    EXPECT_TRUE(std::string(domain_z.name()) == "Z" &&
                    std::string(domain_q.name()) == "Q" &&
                    std::string(domain_r.name()) == "R" &&
                    std::string(domain_c.name()) == "C" &&
                    std::string(domain_expr.name()) == "Expr",
                "predefined number domain sets use LSR names");
    EXPECT_TRUE(domain_z.subset_of(domain_q) &&
                    domain_q.subset_of(domain_r) &&
                    domain_r.subset_of(domain_c) &&
                    domain_c.subset_of(domain_expr),
                "predefined number domains follow Z subset Q subset R subset C subset Expr");
    auto facade_z_q = lamina::lsr::domain_subset(domain_z, domain_q);
    auto facade_c_r = lamina::lsr::domain_subset(domain_c, domain_r);
    auto facade_c_expr = lamina::lsr::domain_subset(domain_c, domain_expr);
    EXPECT_TRUE(facade_z_q && facade_z_q.value(),
                "number domain subset facade accepts Z subset Q");
    EXPECT_TRUE(facade_c_r && !facade_c_r.value(),
                "number domain subset facade rejects C subset R");
    EXPECT_TRUE(facade_c_expr && facade_c_expr.value(),
                "number domain subset facade accepts C subset Expr");

    auto exact_two = lamina::lsr::integer(2);
    auto exact_half = lamina::lsr::rational(
        Rational(BigInt(1), BigInt(2)));
    auto domain_approx_half = lamina::lsr::approx_real(0.5);
    auto explicit_i_for_domain = lamina::lsr::imaginary_unit();
    auto z_contains_two = lamina::lsr::domain_contains(domain_z,
                                                       exact_two.value());
    auto z_contains_half = lamina::lsr::domain_contains(domain_z,
                                                        exact_half.value());
    auto q_contains_half = lamina::lsr::domain_contains(domain_q,
                                                        exact_half.value());
    auto q_contains_approx = lamina::lsr::domain_contains(
        domain_q, domain_approx_half.value());
    auto r_contains_approx = lamina::lsr::domain_contains(
        domain_r, domain_approx_half.value());
    auto r_contains_i = lamina::lsr::domain_contains(
        domain_r, explicit_i_for_domain.value());
    auto c_contains_i = lamina::lsr::domain_contains(
        domain_c, explicit_i_for_domain.value());
    auto r_contains_legacy_i = lamina::lsr::domain_contains(
        domain_r, SymbolicExpr::variable("i"));
    auto c_contains_legacy_i = lamina::lsr::domain_contains(
        domain_c, SymbolicExpr::variable("i"));
    auto c_contains_legacy_upper_i = lamina::lsr::domain_contains(
        domain_c, SymbolicExpr::variable("I"));
    auto legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), legacy_four_i);
    auto c_contains_legacy_four_i =
        lamina::lsr::domain_contains(domain_c, legacy_four_i);
    auto c_contains_legacy_three_plus_four_i =
        lamina::lsr::domain_contains(domain_c, legacy_three_plus_four_i);
    auto expr_contains_symbol = lamina::lsr::domain_contains(
        domain_expr, lamina::lsr::sym("domain_expr_symbol").value());
    EXPECT_TRUE(z_contains_two && z_contains_two.value(),
                "Z contains exact integer literals");
    EXPECT_TRUE(z_contains_half && !z_contains_half.value(),
                "Z rejects non-integer rationals");
    EXPECT_TRUE(q_contains_half && q_contains_half.value(),
                "Q contains exact rational literals");
    EXPECT_TRUE(q_contains_approx && !q_contains_approx.value(),
                "Q does not claim approximate reals as exact rationals");
    EXPECT_TRUE(r_contains_approx && r_contains_approx.value(),
                "R contains finite approximate real literals");
    EXPECT_TRUE(r_contains_i && !r_contains_i.value(),
                "R rejects explicit non-real complex values");
    EXPECT_TRUE(c_contains_i && c_contains_i.value(),
                "C contains explicit complex values");
    EXPECT_TRUE(r_contains_legacy_i && !r_contains_legacy_i.value(),
                "R rejects legacy i variables under the LSR complex boundary");
    EXPECT_TRUE(c_contains_legacy_i && c_contains_legacy_i.value() &&
                    c_contains_legacy_upper_i &&
                    c_contains_legacy_upper_i.value(),
                "C contains legacy i/I variables under the LSR complex boundary");
    EXPECT_TRUE(c_contains_legacy_four_i &&
                    c_contains_legacy_four_i.value() &&
                    c_contains_legacy_three_plus_four_i &&
                    c_contains_legacy_three_plus_four_i.value(),
                "C contains legacy arithmetic expressions built from i/I");
    EXPECT_TRUE(expr_contains_symbol && expr_contains_symbol.value(),
                "Expr contains symbolic expressions");

    auto numeric_domain_set = lamina::lsr::expr_set({
        exact_two.value(), exact_half.value(), domain_approx_half.value()});
    auto complex_domain_set = lamina::lsr::expr_set({
        exact_two.value(), explicit_i_for_domain.value()});
    auto legacy_complex_domain_set = lamina::lsr::expr_set({
        SymbolicExpr::variable("i")});
    auto legacy_complex_arithmetic_domain_set = lamina::lsr::expr_set({
        legacy_three_plus_four_i});
    auto unknown_domain_set = lamina::lsr::expr_set({
        lamina::lsr::sym("domain_set_unknown").value()});
    auto numeric_subset_r =
        numeric_domain_set ? lamina::lsr::expr_set_subset_domain(
                                 numeric_domain_set.value(), domain_r)
                           : lamina::Result<bool>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "numeric domain set construction failed",
                                 "test_lsr_expr");
    auto complex_subset_r =
        complex_domain_set ? lamina::lsr::expr_set_subset_domain(
                                 complex_domain_set.value(), domain_r)
                           : lamina::Result<bool>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "complex domain set construction failed",
                                 "test_lsr_expr");
    auto complex_subset_c =
        complex_domain_set ? lamina::lsr::expr_set_subset_domain(
                                 complex_domain_set.value(), domain_c)
                           : lamina::Result<bool>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "complex domain set construction failed",
                                 "test_lsr_expr");
    auto legacy_complex_subset_c =
        legacy_complex_domain_set ? lamina::lsr::expr_set_subset_domain(
                                        legacy_complex_domain_set.value(),
                                        domain_c)
                                  : lamina::Result<bool>::failure(
                                        lamina::CasErrc::InternalInvariant,
                                        "legacy complex domain set construction failed",
                                        "test_lsr_expr");
    auto legacy_complex_arithmetic_subset_c =
        legacy_complex_arithmetic_domain_set
            ? lamina::lsr::expr_set_subset_domain(
                  legacy_complex_arithmetic_domain_set.value(), domain_c)
            : lamina::Result<bool>::failure(
                  lamina::CasErrc::InternalInvariant,
                  "legacy complex arithmetic domain set construction failed",
                  "test_lsr_expr");
    auto unknown_subset_r =
        unknown_domain_set ? lamina::lsr::expr_set_subset_domain(
                                 unknown_domain_set.value(), domain_r)
                           : lamina::Result<bool>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "unknown domain set construction failed",
                                 "test_lsr_expr");
    auto unknown_subset_expr =
        unknown_domain_set ? lamina::lsr::expr_set_subset_domain(
                                 unknown_domain_set.value(), domain_expr)
                           : lamina::Result<bool>::failure(
                                 lamina::CasErrc::InternalInvariant,
                                 "unknown domain set construction failed",
                                 "test_lsr_expr");
    EXPECT_TRUE(numeric_subset_r && numeric_subset_r.value(),
                "set<Expr> subset facade accepts numeric real sets");
    EXPECT_TRUE(complex_subset_r && !complex_subset_r.value(),
                "set<Expr> subset facade rejects non-real complex members");
    EXPECT_TRUE(complex_subset_c && complex_subset_c.value(),
                "set<Expr> subset facade accepts explicit complex members in C");
    EXPECT_TRUE(legacy_complex_subset_c && legacy_complex_subset_c.value(),
                "set<Expr> subset facade accepts legacy i members in C");
    EXPECT_TRUE(legacy_complex_arithmetic_subset_c &&
                    legacy_complex_arithmetic_subset_c.value(),
                "set<Expr> subset facade accepts legacy complex arithmetic members in C");
    EXPECT_TRUE(!unknown_subset_r &&
                    unknown_subset_r.error().code ==
                        lamina::CasErrc::Inconclusive,
                "set<Expr> subset facade propagates undecidable domain membership");
    EXPECT_TRUE(unknown_subset_expr && unknown_subset_expr.value(),
                "set<Expr> subset facade accepts arbitrary Expr members in Expr");

    auto unknown_domain_member = lamina::lsr::domain_contains(
        domain_r, lamina::lsr::sym("domain_unknown").value());
    auto null_domain_member = lamina::lsr::domain_contains(domain_r, nullptr);
    EXPECT_TRUE(!unknown_domain_member &&
                    unknown_domain_member.error().code ==
                        lamina::CasErrc::Inconclusive,
                "number domain membership does not guess symbolic variables");
    EXPECT_TRUE(!null_domain_member &&
                    null_domain_member.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "number domain membership rejects null Expr values");

    TEST_CASE("LSR solve_expr_set lowers only complete finite CAS results");

    auto finite_solved = lamina::lsr::solve_expr_set(equation, "x");
    EXPECT_TRUE(finite_solved && finite_solved.value().size() == 2,
                "solve_expr_set lowers finite solutions to set<Expr>");

    auto repeated_root_equation =
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2));
    auto repeated_roots = lamina::lsr::roots(repeated_root_equation, "x");
    EXPECT_TRUE(repeated_roots && repeated_roots.value().size() == 1 &&
                    repeated_roots.value().contains(*SymbolicExpr::number(0)),
                "roots lowers repeated roots to one set<Expr> member");

    auto cubic_equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(3)),
        SymbolicExpr::number(-2));
    auto cubic_roots = lamina::lsr::roots(cubic_equation, "x");
    EXPECT_TRUE(cubic_roots && cubic_roots.value().size() == 3,
                "roots lowers exact higher-degree polynomial roots to finite set<Expr>");
    if (cubic_roots) {
        bool all_root_of = true;
        for (const auto& root : cubic_roots.value().elements()) {
            auto function =
                root ? std::dynamic_pointer_cast<const FunctionNode>(
                           lamina::detail::node(root))
                     : nullptr;
            all_root_of = all_root_of && function &&
                          function->type() == FunctionNode::FuncType::RootOf;
        }
        EXPECT_TRUE(all_root_of,
                    "higher-degree finite polynomial roots remain explicit RootOf Expr values");
    }

    auto roots_solved = lamina::lsr::roots(equation, "x");
    EXPECT_TRUE(roots_solved && roots_solved.value().size() == 2 &&
                    roots_solved.value().contains(*SymbolicExpr::number(-1)) &&
                    roots_solved.value().contains(*SymbolicExpr::number(1)),
                "roots returns the LSR set<Expr> finite root collection");

    auto solve_solved = lamina::lsr::solve(equation, "x");
    EXPECT_TRUE(solve_solved && solve_solved.value().size() == 2 &&
                    solve_solved.value().contains(*SymbolicExpr::number(-1)) &&
                    solve_solved.value().contains(*SymbolicExpr::number(1)),
                "solve returns the LSR set<Expr> finite solution collection");

    auto complex_equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto complex_solved = lamina::lsr::solve_expr_set(complex_equation, "x");
    auto negative_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(-1));
    EXPECT_TRUE(complex_solved && complex_solved.value().size() == 2,
                "solve_expr_set returns both complex roots for x^2 + 1");
    auto complex_roots_subset_c =
        complex_solved ? lamina::lsr::expr_set_subset_domain(
                             complex_solved.value(), domain_c)
                       : lamina::Result<bool>::failure(
                             lamina::CasErrc::InternalInvariant,
                             "complex solve set construction failed",
                             "test_lsr_expr");
    EXPECT_TRUE(complex_roots_subset_c && complex_roots_subset_c.value(),
                "solve(x^2 + 1, x) subset C is directly checkable");
    EXPECT_TRUE(complex_solved && i && negative_i &&
                    complex_solved.value().contains(*i.value()) &&
                    complex_solved.value().contains(*negative_i.value()),
                "solve_expr_set lowers x^2 + 1 roots to explicit LSR complex expressions");
    if (complex_solved) {
        bool saw_positive_i = false;
        bool saw_negative_i = false;
        for (const auto& root : complex_solved.value().elements()) {
            auto lowered_root = lamina::lsr::eval_complex(*root);
            EXPECT_TRUE(lowered_root && lowered_root.value().is_finite(),
                        "complex solve roots explicitly lower to complex values");
            if (lowered_root && lowered_root.value().real.value == 0.0 &&
                lowered_root.value().imag.value == 1.0) {
                saw_positive_i = true;
            }
            if (lowered_root && lowered_root.value().real.value == 0.0 &&
                lowered_root.value().imag.value == -1.0) {
                saw_negative_i = true;
            }
        }
        EXPECT_TRUE(saw_positive_i && saw_negative_i,
                    "complex solve roots evaluate to -i and i");
    }

    auto shifted_complex_equation = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
                          SymbolicExpr::multiply(SymbolicExpr::number(2), x.value())),
        SymbolicExpr::number(2));
    auto shifted_complex_solved =
        lamina::lsr::solve_expr_set(shifted_complex_equation, "x");
    EXPECT_TRUE(shifted_complex_solved &&
                    shifted_complex_solved.value().size() == 2,
                "solve_expr_set returns both complex roots for x^2 + 2x + 2");
    if (shifted_complex_solved) {
        bool saw_negative_one_plus_i = false;
        bool saw_negative_one_minus_i = false;
        for (const auto& root : shifted_complex_solved.value().elements()) {
            auto lowered_root = lamina::lsr::eval_complex(*root);
            EXPECT_TRUE(lowered_root && lowered_root.value().is_finite(),
                        "shifted complex roots explicitly lower to complex values");
            if (lowered_root && lowered_root.value().real.value == -1.0 &&
                lowered_root.value().imag.value == 1.0) {
                saw_negative_one_plus_i = true;
            }
            if (lowered_root && lowered_root.value().real.value == -1.0 &&
                lowered_root.value().imag.value == -1.0) {
                saw_negative_one_minus_i = true;
            }
        }
        EXPECT_TRUE(saw_negative_one_plus_i && saw_negative_one_minus_i,
                    "complex solve roots preserve nonzero real components");
    }

    auto empty_solved = lamina::lsr::solve_expr_set(SymbolicExpr::number(1), "x");
    EXPECT_TRUE(empty_solved && empty_solved.value().empty(),
                "solve_expr_set lowers mathematical no-solution to empty set<Expr>");

    auto universal_solved = lamina::lsr::solve_expr_set(SymbolicExpr::number(0), "x");
    EXPECT_TRUE(!universal_solved &&
                    universal_solved.error().code == lamina::CasErrc::Inconclusive,
                "solve_expr_set does not pretend universal solutions are finite set<Expr>");
    EXPECT_TRUE(!universal_solved &&
                    std::string(lamina::lsr::error_name(universal_solved.error())) ==
                        "SetResultInconclusive",
                "non-finite solution lowering exposes the LSR set inconclusive diagnostic");

    auto null_solve_input = lamina::lsr::solve_expr_set(nullptr, "x");
    auto null_solve_set_input = lamina::lsr::solve_set(nullptr, "x");
    auto null_roots_input = lamina::lsr::roots(nullptr, "x");
    auto null_solve_alias_input = lamina::lsr::solve(nullptr, "x");
    EXPECT_TRUE(!null_solve_input &&
                    null_solve_input.error().code == lamina::CasErrc::InvalidArgument,
                "solve_expr_set rejects null equations before lowering");
    EXPECT_TRUE(!null_solve_set_input &&
                    null_solve_set_input.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "solve_set rejects null equations");
    EXPECT_TRUE(!null_roots_input &&
                    null_roots_input.error().code == lamina::CasErrc::InvalidArgument,
                "roots rejects null expressions");
    EXPECT_TRUE(!null_solve_alias_input &&
                    null_solve_alias_input.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "solve rejects null equations");

    auto empty_solve_set_variable = lamina::lsr::solve_set(equation, "");
    auto empty_solve_variable = lamina::lsr::solve_expr_set(equation, "");
    auto empty_roots_variable = lamina::lsr::roots(equation, "");
    auto empty_solve_alias_variable = lamina::lsr::solve(equation, "");
    EXPECT_TRUE(!empty_solve_set_variable &&
                    empty_solve_set_variable.error().code ==
                        lamina::CasErrc::InvalidArgument,
                "solve_set rejects empty variable names");
    EXPECT_TRUE(!empty_solve_variable &&
                    empty_solve_variable.error().code == lamina::CasErrc::InvalidArgument,
                "solve_expr_set rejects empty variable names");
    EXPECT_TRUE(!empty_roots_variable &&
                    empty_roots_variable.error().code == lamina::CasErrc::InvalidArgument,
                "roots rejects empty variable names");
    EXPECT_TRUE(!empty_solve_alias_variable &&
                    empty_solve_alias_variable.error().code == lamina::CasErrc::InvalidArgument,
                "solve rejects empty variable names");

    TEST_CASE("LSR equivalence core handles local exact identities");

    auto one_plus_x = SymbolicExpr::add(SymbolicExpr::number(1), x.value());
    auto x_plus_one = SymbolicExpr::add(x.value(), SymbolicExpr::number(1));
    lamina::ComputationContext equivalent_context;
    auto equivalent = lamina::lsr::equivalent_core(*one_plus_x,
                                                   *x_plus_one,
                                                   equivalent_context);
    EXPECT_TRUE(equivalent && equivalent.value(),
                "equivalent_core proves normalized additive equality");

    auto x_times_one = SymbolicExpr::multiply(x.value(), SymbolicExpr::number(1));
    auto x_times_zero = SymbolicExpr::multiply(x.value(), SymbolicExpr::number(0));
    auto x_minus_x = SymbolicExpr::add(
        x.value(), SymbolicExpr::multiply(SymbolicExpr::number(-1), x.value()));
    lamina::ComputationContext multiply_identity_context;
    auto multiply_identity = lamina::lsr::equivalent_core(
        *x_times_one, *x.value(), multiply_identity_context);
    lamina::ComputationContext multiply_zero_context;
    auto multiply_zero = lamina::lsr::equivalent_core(
        *x_times_zero, *SymbolicExpr::number(0), multiply_zero_context);
    lamina::ComputationContext subtract_self_context;
    auto subtract_self = lamina::lsr::equivalent_core(
        *x_minus_x, *SymbolicExpr::number(0), subtract_self_context);
    EXPECT_TRUE(multiply_identity && multiply_identity.value(),
                "equivalent_core proves Core x * 1 identity");
    EXPECT_TRUE(multiply_zero && multiply_zero.value(),
                "equivalent_core proves Core x * 0 identity");
    EXPECT_TRUE(subtract_self && subtract_self.value(),
                "equivalent_core proves Core x - x identity");

    auto x_plus_one_squared = SymbolicExpr::power(x_plus_one, SymbolicExpr::number(2));
    auto x_squared = SymbolicExpr::power(x.value(), SymbolicExpr::number(2));
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x.value());
    auto expanded_square = SymbolicExpr::add(
        SymbolicExpr::add(x_squared, two_x), SymbolicExpr::number(1));
    lamina::ComputationContext polynomial_eqv_context;
    auto polynomial_equivalent = lamina::lsr::equivalent_core(
        *x_plus_one_squared, *expanded_square, polynomial_eqv_context);
    EXPECT_TRUE(polynomial_equivalent && polynomial_equivalent.value(),
                "equivalent_core proves the LSR polynomial Core example");

    auto y = lamina::lsr::sym("y");
    auto y_plus_one = SymbolicExpr::add(y.value(), SymbolicExpr::number(1));
    auto y_plus_one_squared = SymbolicExpr::power(y_plus_one, SymbolicExpr::number(2));
    auto y_squared = SymbolicExpr::power(y.value(), SymbolicExpr::number(2));
    auto two_y = SymbolicExpr::multiply(SymbolicExpr::number(2), y.value());
    auto expanded_y_square = SymbolicExpr::add(
        SymbolicExpr::add(y_squared, two_y), SymbolicExpr::number(1));
    lamina::ComputationContext y_polynomial_eqv_context;
    auto y_polynomial_equivalent = lamina::lsr::equivalent_core(
        *y_plus_one_squared, *expanded_y_square, y_polynomial_eqv_context);
    EXPECT_TRUE(y_polynomial_equivalent && y_polynomial_equivalent.value(),
                "equivalent_core polynomial proof is not hard-coded to x");

    lamina::lsr::EqvOptions no_budget;
    auto invalid_budget = lamina::lsr::set_eqv_budget(no_budget, 0, 64, 4);
    lamina::lsr::EqvOptions no_depth_budget;
    auto invalid_depth_budget =
        lamina::lsr::set_eqv_budget(no_depth_budget, 256, 0, 4);
    lamina::lsr::EqvOptions no_growth_budget;
    auto invalid_growth_budget =
        lamina::lsr::set_eqv_budget(no_growth_budget, 256, 64, 0);
    EXPECT_TRUE(!invalid_budget &&
                    invalid_budget.error().code == lamina::CasErrc::ResourceLimit,
                "set_eqv_budget rejects zero rewrite steps");
    EXPECT_TRUE(!invalid_depth_budget &&
                    invalid_depth_budget.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "set_eqv_budget rejects zero rewrite depth");
    EXPECT_TRUE(!invalid_growth_budget &&
                    invalid_growth_budget.error().code ==
                        lamina::CasErrc::ResourceLimit,
                "set_eqv_budget rejects zero node growth factor");
    EXPECT_TRUE(!invalid_budget &&
                    std::string(lamina::lsr::error_name(invalid_budget.error())) ==
                        "EqvBudgetExceeded",
                "set_eqv_budget exposes EqvBudgetExceeded for invalid budgets");
    EXPECT_TRUE(!invalid_depth_budget &&
                    std::string(lamina::lsr::error_name(
                        invalid_depth_budget.error())) == "EqvBudgetExceeded",
                "zero rewrite depth exposes EqvBudgetExceeded");
    EXPECT_TRUE(!invalid_growth_budget &&
                    std::string(lamina::lsr::error_name(
                        invalid_growth_budget.error())) == "EqvBudgetExceeded",
                "zero node growth factor exposes EqvBudgetExceeded");
    no_budget.budget.max_rewrite_steps = 0;
    lamina::ComputationContext no_budget_context;
    auto exhausted_eqv = lamina::lsr::equivalent_core(
        *one_plus_x, *x_plus_one, no_budget_context, no_budget);
    EXPECT_TRUE(!exhausted_eqv &&
                    exhausted_eqv.error().code == lamina::CasErrc::ResourceLimit,
                "equivalent_core observes explicit LSR rewrite budgets");
    EXPECT_TRUE(!exhausted_eqv &&
                    std::string(lamina::lsr::error_name(exhausted_eqv.error())) ==
                        "EqvBudgetExceeded",
                "LSR equivalence budget exhaustion exposes EqvBudgetExceeded");
    lamina::ComputationContext lsr_no_budget_context;
    auto lsr_exhausted_eqv = lamina::lsr::equivalent(
        *one_plus_x, *x_plus_one, lsr_no_budget_context, no_budget);
    EXPECT_TRUE(lsr_exhausted_eqv && !lsr_exhausted_eqv.value(),
                "LSR equivalent returns false when the proof budget is exhausted");
    EXPECT_TRUE(!lsr_no_budget_context.diagnostics().empty(),
                "LSR equivalent records a diagnostic for budget exhaustion");

    auto sin_y_squared = SymbolicExpr::power(
        SymbolicExpr::sin(y.value()), SymbolicExpr::number(2));
    auto cos_y_squared = SymbolicExpr::power(
        SymbolicExpr::cos(y.value()), SymbolicExpr::number(2));
    auto trig_identity = SymbolicExpr::add(sin_y_squared, cos_y_squared);
    lamina::lsr::EqvOptions trig_profile;
    auto trig_configured = lamina::lsr::set_eqv_profile(trig_profile, "Trig-Basic");
    EXPECT_TRUE(trig_configured.has_value(),
                "set_eqv_profile accepts the LSR Trig-Basic profile name");
    lamina::ComputationContext trig_profile_context;
    auto trig_equivalent = lamina::lsr::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), trig_profile_context,
        trig_profile);
    EXPECT_TRUE(trig_equivalent && trig_equivalent.value(),
                "equivalent_core proves the LSR Trig-Basic sin^2+cos^2 identity");

    auto negative_y = SymbolicExpr::multiply(SymbolicExpr::number(-1), y.value());
    auto sin_negative_y = SymbolicExpr::sin(negative_y);
    auto negative_sin_y = SymbolicExpr::multiply(
        SymbolicExpr::number(-1), SymbolicExpr::sin(y.value()));
    lamina::ComputationContext trig_sin_odd_context;
    auto sin_odd_equivalent = lamina::lsr::equivalent_core(
        *sin_negative_y, *negative_sin_y, trig_sin_odd_context, trig_profile);
    EXPECT_TRUE(sin_odd_equivalent && sin_odd_equivalent.value(),
                "equivalent_core proves the LSR Trig-Basic sin odd identity");

    auto cos_negative_y = SymbolicExpr::cos(negative_y);
    auto cos_y = SymbolicExpr::cos(y.value());
    lamina::ComputationContext trig_cos_even_context;
    auto cos_even_equivalent = lamina::lsr::equivalent_core(
        *cos_negative_y, *cos_y, trig_cos_even_context, trig_profile);
    EXPECT_TRUE(cos_even_equivalent && cos_even_equivalent.value(),
                "equivalent_core proves the LSR Trig-Basic cos even identity");

    lamina::ComputationContext core_trig_context;
    auto core_trig_equivalent = lamina::lsr::equivalent_core(
        *trig_identity, *SymbolicExpr::number(1), core_trig_context);
    EXPECT_TRUE(core_trig_equivalent && !core_trig_equivalent.value(),
                "Core profile does not silently enable Trig-Basic rules");
    lamina::ComputationContext lsr_trig_profile_context;
    auto lsr_trig_equivalent = lamina::lsr::equivalent(
        *trig_identity, *SymbolicExpr::number(1), lsr_trig_profile_context,
        trig_profile);
    EXPECT_TRUE(lsr_trig_equivalent && lsr_trig_equivalent.value(),
                "LSR equivalent proves enabled Trig-Basic rules");

    lamina::lsr::EqvOptions exp_log_profile;
    auto exp_log_configured =
        lamina::lsr::set_eqv_profile(exp_log_profile, "ExpLog-Basic");
    EXPECT_TRUE(exp_log_configured.has_value(),
                "set_eqv_profile accepts the LSR ExpLog-Basic profile name");
    auto unsupported_profile =
        lamina::lsr::eqv_profile_from_name("Richardson-Complete");
    EXPECT_TRUE(!unsupported_profile &&
                    unsupported_profile.error().code ==
                        lamina::CasErrc::UnsupportedExpression,
                "eqv_profile_from_name rejects unsupported profile names");
    EXPECT_TRUE(!unsupported_profile &&
                    std::string(lamina::lsr::error_name(
                        unsupported_profile.error())) == "EqvRuleDisabled",
                "unsupported equivalence profile exposes EqvRuleDisabled");
    lamina::ComputationContext exp_log_profile_context;
    auto exp_zero = SymbolicExpr::exp(SymbolicExpr::number(0));
    auto exp_log_equivalent = lamina::lsr::equivalent_core(
        *exp_zero, *SymbolicExpr::number(1), exp_log_profile_context,
        exp_log_profile);
    EXPECT_TRUE(exp_log_equivalent && exp_log_equivalent.value(),
                "equivalent_core proves the LSR ExpLog-Basic exp(0) identity");

    auto ln_one = SymbolicExpr::ln(SymbolicExpr::number(1));
    lamina::ComputationContext ln_one_context;
    auto ln_one_equivalent = lamina::lsr::equivalent_core(
        *ln_one, *SymbolicExpr::number(0), ln_one_context, exp_log_profile);
    EXPECT_TRUE(ln_one_equivalent && ln_one_equivalent.value(),
                "equivalent_core proves the LSR ExpLog-Basic ln(1) identity");

    auto exp_ln_y = SymbolicExpr::exp(SymbolicExpr::ln(y.value()));
    lamina::ComputationContext exp_ln_unproven_context;
    auto exp_ln_unproven = lamina::lsr::equivalent_core(
        *exp_ln_y, *y.value(), exp_ln_unproven_context, exp_log_profile);
    EXPECT_TRUE(exp_ln_unproven && !exp_ln_unproven.value(),
                "ExpLog-Basic does not prove exp(ln(y)) without domain evidence");

    auto positive_assumptions = std::make_shared<lamina::AssumptionContext>();
    positive_assumptions->assume_sign("y", lamina::Sign::Positive);
    lamina::ComputationContext exp_ln_positive_context;
    auto set_positive_assumptions =
        exp_ln_positive_context.set_assumptions(positive_assumptions);
    EXPECT_TRUE(set_positive_assumptions.has_value(),
                "equivalence context accepts positive assumptions");
    auto exp_ln_positive = lamina::lsr::equivalent_core(
        *exp_ln_y, *y.value(), exp_ln_positive_context, exp_log_profile);
    EXPECT_TRUE(exp_ln_positive && exp_ln_positive.value(),
                "equivalent_core proves ExpLog-Basic exp(ln(y)) for positive y");

    return TEST_REPORT();
}
