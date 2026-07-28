#include "../include/lsr_expr.hpp"
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

    auto reserved_pi = lamina::lsr::sym("pi");
    EXPECT_TRUE(!reserved_pi &&
                    reserved_pi.error().code == lamina::CasErrc::InvalidArgument,
                "pi is a constant and cannot be shadowed as a symbol");
    EXPECT_TRUE(!reserved_pi &&
                    std::string(lamina::lsr::error_name(reserved_pi.error())) ==
                        "InvalidArgument",
                "non-imaginary reserved constants keep the generic argument diagnostic");

    TEST_CASE("LSR imaginary unit is canonical complex zero plus one i");

    auto i = lamina::lsr::imaginary_unit();
    EXPECT_TRUE(i.has_value(), "imaginary unit can be constructed");

    auto explicit_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(1));
    EXPECT_TRUE(explicit_i.has_value(), "complex(0, 1) can be constructed");
    EXPECT_TRUE(i && explicit_i &&
                    lamina::lsr::structurally_equal(*i.value(),
                                                    *explicit_i.value()),
                "imaginary unit is structurally complex(0, 1)");

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
    }

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

    auto rhs_set = lamina::lsr::expr_set({SymbolicExpr::number(2),
                                          SymbolicExpr::number(3)});
    EXPECT_TRUE(rhs_set.has_value(), "second set<Expr> can be created");
    if (base_set && rhs_set) {
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

    TEST_CASE("LSR equivalence core handles local exact identities");

    auto one_plus_x = SymbolicExpr::add(SymbolicExpr::number(1), x.value());
    auto x_plus_one = SymbolicExpr::add(x.value(), SymbolicExpr::number(1));
    lamina::ComputationContext equivalent_context;
    auto equivalent = lamina::lsr::equivalent_core(*one_plus_x,
                                                   *x_plus_one,
                                                   equivalent_context);
    EXPECT_TRUE(equivalent && equivalent.value(),
                "equivalent_core proves normalized additive equality");

    return TEST_REPORT();
}
