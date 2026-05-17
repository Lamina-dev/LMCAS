#include "test_common.hpp"
#include <iostream>

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

    return TEST_REPORT();
}