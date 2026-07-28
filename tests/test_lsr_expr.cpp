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

    auto reserved_pi = lamina::lsr::sym("pi");
    EXPECT_TRUE(!reserved_pi &&
                    reserved_pi.error().code == lamina::CasErrc::InvalidArgument,
                "pi is a constant and cannot be shadowed as a symbol");

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

    lamina::NumericBindings bindings{{"x", 3.0}};
    auto evaluated = lamina::lsr::evalf(*linear, bindings);
    EXPECT_TRUE(evaluated && evaluated.value().is_finite(),
                "evalf with a binding succeeds");
    EXPECT_NEAR(evaluated.value().value, 5.0, 0.0,
                "evalf computes the numeric value");

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
