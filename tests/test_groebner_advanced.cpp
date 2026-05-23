#include "solver.hpp"
#include "symbolic.hpp"
#include <iostream>
#include <cassert>

int main() {
    using namespace lamina;

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    std::cout << "Test 1: Reduced Groebner Basis" << std::endl;
    {

        auto p1 = *SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                              SymbolicExpr::power(y, SymbolicExpr::number(2))),
            SymbolicExpr::number(-1));
        auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

        auto rgb = Solver::reduced_groebner_basis({p1, p2}, {"x", "y"});
        std::cout << "  Reduced basis size: " << rgb.size() << std::endl;
        for (auto& g : rgb) {
            std::cout << "  " << g.to_string() << std::endl;
        }
        assert(!rgb.empty());
        std::cout << "  [PASS]" << std::endl;
    }

    std::cout << "Test 2: Ideal Membership" << std::endl;
    {

        auto p1 = *SymbolicExpr::add(x, y);
        auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

        auto gb = Solver::groebner_basis({p1, p2}, {"x", "y"});
        std::cout << "  Basis: ";
        for (auto& g : gb) std::cout << g.to_string() << " ; ";
        std::cout << std::endl;

        bool member1 = Solver::ideal_membership(p1, gb, {"x", "y"});
        std::cout << "  x+y in ideal: " << (member1 ? "true" : "false") << std::endl;
        assert(member1);

        auto test_poly = *SymbolicExpr::add(x, SymbolicExpr::number(1));
        bool member2 = Solver::ideal_membership(test_poly, gb, {"x", "y"});
        std::cout << "  x+1 in ideal: " << (member2 ? "true" : "false") << std::endl;

        std::cout << "  [PASS]" << std::endl;
    }

    std::cout << "Test 3: Elimination Ideal" << std::endl;
    {

        auto p1 = *SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                              SymbolicExpr::power(y, SymbolicExpr::number(2))),
            SymbolicExpr::number(-1));
        auto p2 = *SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));

        auto gb = Solver::groebner_basis({p1, p2}, {"x", "y"});
        auto elim = Solver::elimination_ideal(gb, {"x", "y"}, 1);

        std::cout << "  Full basis size: " << gb.size() << std::endl;
        std::cout << "  Elimination ideal (eliminate x): " << elim.size() << " elements" << std::endl;
        for (auto& e : elim) {
            std::cout << "    " << e.to_string() << std::endl;
        }

        assert(!elim.empty());
        std::cout << "  [PASS]" << std::endl;
    }

    std::cout << "\nAll advanced Groebner tests passed!" << std::endl;
    return 0;
}
