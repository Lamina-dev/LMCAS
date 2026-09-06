#include "test_common.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace LMCAS;

int main() {

    BigInt pow20("100000000000000000000");
    BigInt pow20_plus_1("100000000000000000001");

    auto node1 = SymbolicExpr::number(pow20);
    auto node2 = SymbolicExpr::number(pow20_plus_1);

    EXPECT_TRUE(node1->compare(node2) < 0, "10^20 should be less than 10^20 + 1");
    EXPECT_TRUE(node2->compare(node1) > 0, "10^20 + 1 should be greater than 10^20");
    EXPECT_TRUE(node1->compare(node1) == 0, "10^20 should equal 10^20");

    Rational one_third(1, 3);
    lmmc_real_t decimal_third = 0.33333333;

    auto node3 = SymbolicExpr::number(one_third);
    auto node4 = SymbolicExpr::number(decimal_third);

    EXPECT_TRUE(node3->compare(node4) > 0, "1/3 should be greater than 0.33333333");
    EXPECT_TRUE(node4->compare(node3) < 0, "0.33333333 should be less than 1/3");

    auto node5 = SymbolicExpr::number(Rational(1, 3));
    auto node6 = SymbolicExpr::number(Rational(33333333, 100000000));

    EXPECT_TRUE(node5->compare(node6) > 0, "1/3 should be greater than 33333333/100000000 exactly");

    auto approximate = SymbolicExpr::number(1.5);
    auto approximate_value = approximate->get_number();
    EXPECT_TRUE(std::holds_alternative<Rational>(approximate_value) &&
                    std::get<Rational>(approximate_value) ==
                        Rational::from_double(1.5),
                "get_number preserves a non-integral approximate value");

    bool large_bigint_threw = false;
    try {
        (void)SymbolicExpr::number(
            BigInt("999999999999999999999999999999"))->get_int();
    } catch (const std::out_of_range&) {
        large_bigint_threw = true;
    }
    EXPECT_TRUE(large_bigint_threw,
                "get_int rejects BigInt values outside the int range");

    bool large_approximate_threw = false;
    try {
        (void)SymbolicExpr::number(
            static_cast<double>(std::numeric_limits<int>::max()) * 2.0)->get_int();
    } catch (const std::out_of_range&) {
        large_approximate_threw = true;
    }
    EXPECT_TRUE(large_approximate_threw,
                "get_int rejects approximate integers outside the int range");

    return TEST_REPORT();
}
