#include "monomial_order.hpp"
#include "test_common.hpp"
#include <iostream>

int main() {
    using namespace lamina;

    Monomial a = {2, 1, 0};
    Monomial b = {1, 0, 2};
    Monomial c = {3, 0, 0};
    Monomial d = {1, 1, 0};

    EXPECT_TRUE(total_degree(a) == 3, "total_degree(a) is 3");
    EXPECT_TRUE(total_degree(b) == 3, "total_degree(b) is 3");
    EXPECT_TRUE(total_degree(d) == 2, "total_degree(d) is 2");
    std::cout << "[PASS] total_degree" << std::endl;

    MonomialOrder lex = MonomialOrder::lex();
    EXPECT_TRUE(lex(a, b), "lex orders a before b");
    EXPECT_FALSE(lex(b, a), "lex does not order b before a");
    EXPECT_TRUE(lex(c, a), "lex orders c before a");
    std::cout << "[PASS] Lex order" << std::endl;

    MonomialOrder grevlex = MonomialOrder::grevlex();
    EXPECT_TRUE(grevlex(a, b), "grevlex orders a before b");
    EXPECT_FALSE(grevlex(b, a), "grevlex does not order b before a");

    EXPECT_TRUE(grevlex(a, d), "grevlex orders higher degree a before d");
    std::cout << "[PASS] GrevLex order" << std::endl;

    MonomialOrder deglex = MonomialOrder::deglex();
    EXPECT_TRUE(deglex(a, b), "deglex orders a before b");
    EXPECT_TRUE(deglex(a, d), "deglex orders higher degree a before d");
    std::cout << "[PASS] DegLex order" << std::endl;

    Monomial div = {1, 0, 0};
    EXPECT_TRUE(divides_monomial(div, a), "div divides a");
    EXPECT_FALSE(divides_monomial(a, div), "a does not divide div");
    std::cout << "[PASS] divides_monomial" << std::endl;

    Monomial lcm = lcm_monomial(a, b);
    EXPECT_TRUE(lcm[0] == 2 && lcm[1] == 1 && lcm[2] == 2,
                "lcm_monomial combines max exponents");
    std::cout << "[PASS] lcm_monomial" << std::endl;

    std::cout << "\nAll monomial order tests passed!" << std::endl;
    return TEST_REPORT();
}
