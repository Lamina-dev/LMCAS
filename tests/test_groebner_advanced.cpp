#include "solver.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

TEST(LmcasGroebnerAdvanced, ReducedBasis) {
    using namespace LMCAS;
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto p1 = *SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::number(-1));
    auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

    auto rgb = Solver::reduced_groebner_basis({p1, p2}, {"x", "y"});
    ASSERT_FALSE(rgb.empty());
    EXPECT_TRUE(Solver::ideal_membership(p1, rgb, {"x", "y"}));
    EXPECT_TRUE(Solver::ideal_membership(p2, rgb, {"x", "y"}));
}

TEST(LmcasGroebnerAdvanced, IdealMembership) {
    using namespace LMCAS;
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto p1 = *SymbolicExpr::add(x, y);
    auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

    auto gb = Solver::groebner_basis({p1, p2}, {"x", "y"});

    ASSERT_FALSE(gb.empty());
    EXPECT_TRUE(Solver::ideal_membership(p1, gb, {"x", "y"}));
    EXPECT_TRUE(Solver::ideal_membership(p2, gb, {"x", "y"}));

    auto test_poly = *SymbolicExpr::add(x, SymbolicExpr::number(1));
    EXPECT_FALSE(Solver::ideal_membership(
        test_poly, gb, {"x", "y"}));
}

TEST(LmcasGroebnerAdvanced, EliminationIdeal) {
    using namespace LMCAS;
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto p1 = *SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::number(-1));
    auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

    auto gb = Solver::groebner_basis({p1, p2}, {"x", "y"});
    auto elim = Solver::elimination_ideal(gb, {"x", "y"}, 1);

    ASSERT_FALSE(elim.empty());
    const auto eliminated_relation = *SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::power(y, SymbolicExpr::number(2))), SymbolicExpr::number(-1));
    EXPECT_TRUE(Solver::ideal_membership(eliminated_relation, elim, {"x", "y"}));
}
