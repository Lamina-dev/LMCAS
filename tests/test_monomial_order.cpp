#include "monomial_order.hpp"
#include <iostream>
#include <cassert>

int main() {
    using namespace lamina;
    
    // Test total_degree
    Monomial a = {2, 1, 0}; // x^2*y, degree 3
    Monomial b = {1, 0, 2}; // x*z^2, degree 3
    Monomial c = {3, 0, 0}; // x^3, degree 3
    Monomial d = {1, 1, 0}; // x*y, degree 2
    
    assert(total_degree(a) == 3);
    assert(total_degree(b) == 3);
    assert(total_degree(d) == 2);
    std::cout << "[PASS] total_degree" << std::endl;
    
    // Test Lex order: x^2*y > x*z^2 (first component 2 > 1)
    MonomialOrder lex = MonomialOrder::lex();
    assert(lex(a, b) == true);   // {2,1,0} > {1,0,2}
    assert(lex(b, a) == false);
    assert(lex(c, a) == true);   // {3,0,0} > {2,1,0}
    std::cout << "[PASS] Lex order" << std::endl;
    
    // Test GrevLex order: same total degree, compare from right
    // a={2,1,0} vs b={1,0,2}: rightmost diff at index 2: 0-2=-2 < 0, so a > b
    MonomialOrder grevlex = MonomialOrder::grevlex();
    assert(grevlex(a, b) == true);  // a > b under grevlex
    assert(grevlex(b, a) == false);
    // Higher total degree always wins
    assert(grevlex(a, d) == true);  // degree 3 > degree 2
    std::cout << "[PASS] GrevLex order" << std::endl;
    
    // Test DegLex order: total degree first, then lex
    MonomialOrder deglex = MonomialOrder::deglex();
    assert(deglex(a, b) == true);   // same degree, lex: {2,1,0} > {1,0,2}
    assert(deglex(a, d) == true);   // degree 3 > degree 2
    std::cout << "[PASS] DegLex order" << std::endl;
    
    // Test divides_monomial
    Monomial div = {1, 0, 0};
    assert(divides_monomial(div, a) == true);   // x | x^2*y
    assert(divides_monomial(a, div) == false);  // x^2*y does not divide x
    std::cout << "[PASS] divides_monomial" << std::endl;
    
    // Test lcm_monomial
    Monomial lcm = lcm_monomial(a, b); // lcm({2,1,0}, {1,0,2}) = {2,1,2}
    assert(lcm[0] == 2 && lcm[1] == 1 && lcm[2] == 2);
    std::cout << "[PASS] lcm_monomial" << std::endl;
    
    std::cout << "\nAll monomial order tests passed!" << std::endl;
    return 0;
}
