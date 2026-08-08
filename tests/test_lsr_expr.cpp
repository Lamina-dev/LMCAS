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

    TEST_CASE("LSR eval_complex explicitly lowers Expr to complex");

    auto real_as_complex = lamina::lsr::eval_complex(*SymbolicExpr::number(5));
    EXPECT_TRUE(real_as_complex && real_as_complex.value().is_finite(),
                "eval_complex accepts real expressions explicitly");
    EXPECT_NEAR(real_as_complex.value().real.value, 5.0, 0.0,
                "eval_complex preserves real component");
    EXPECT_NEAR(real_as_complex.value().imag.value, 0.0, 0.0,
                "eval_complex promotes real expression with zero imaginary component");

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
    auto empty_expr_set = lamina::lsr::expr_set({});
    EXPECT_TRUE(empty_expr_set && empty_expr_set.value().empty(),
                "set<Expr> can represent the empty finite set");

    auto rhs_set = lamina::lsr::expr_set({SymbolicExpr::number(2),
                                          SymbolicExpr::number(3)});
    EXPECT_TRUE(rhs_set.has_value(), "second set<Expr> can be created");
    if (base_set && rhs_set && empty_expr_set) {
        auto union_set = base_set.value().set_union(rhs_set.value());
        EXPECT_TRUE(union_set.size() == 3,
                    "set<Expr> union returns deduplicated elements");
        auto intersection_set = base_set.value().intersection(rhs_set.value());
        EXPECT_TRUE(intersection_set.size() == 1 &&
                        intersection_set.contains(*SymbolicExpr::number(2)),
                    "set<Expr> intersection keeps common elements");
        auto difference_set = base_set.value().difference(rhs_set.value());
        EXPECT_TRUE(difference_set.size() == 1 &&
                        difference_set.contains(*SymbolicExpr::number(1)),
                    "set<Expr> difference removes right-hand elements");
        auto symmetric = base_set.value().symmetric_difference(rhs_set.value());
        EXPECT_TRUE(symmetric.size() == 2 &&
                        symmetric.contains(*SymbolicExpr::number(1)) &&
                        symmetric.contains(*SymbolicExpr::number(3)),
                    "set<Expr> symmetric difference follows xor semantics");
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
