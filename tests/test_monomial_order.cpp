#include "monomial_order.hpp"
#include <iostream>
#include <cassert>

int main() {
    using namespace lamina;

    Monomial a = {2, 1, 0};
    Monomial b = {1, 0, 2};
    Monomial c = {3, 0, 0};
    Monomial d = {1, 1, 0};

    assert(total_degree(a) == 3);
    assert(total_degree(b) == 3);
    assert(total_degree(d) == 2);
    std::cout << "[PASS] total_degree" << std::endl;

    MonomialOrder lex = MonomialOrder::lex();
    assert(lex(a, b) == true);
    assert(lex(b, a) == false);
    assert(lex(c, a) == true);
    std::cout << "[PASS] Lex order" << std::endl;

    MonomialOrder grevlex = MonomialOrder::grevlex();
    assert(grevlex(a, b) == true);
    assert(grevlex(b, a) == false);

    assert(grevlex(a, d) == true);
    std::cout << "[PASS] GrevLex order" << std::endl;

    MonomialOrder deglex = MonomialOrder::deglex();
    assert(deglex(a, b) == true);
    assert(deglex(a, d) == true);
    std::cout << "[PASS] DegLex order" << std::endl;

    Monomial div = {1, 0, 0};
    assert(divides_monomial(div, a) == true);
    assert(divides_monomial(a, div) == false);
    std::cout << "[PASS] divides_monomial" << std::endl;

    Monomial lcm = lcm_monomial(a, b);
    assert(lcm[0] == 2 && lcm[1] == 1 && lcm[2] == 2);
    std::cout << "[PASS] lcm_monomial" << std::endl;

    std::cout << "\nAll monomial order tests passed!" << std::endl;
    return 0;
}
